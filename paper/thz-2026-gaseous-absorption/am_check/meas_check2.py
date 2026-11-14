"""Further comparisons with measurement (round 3, review D13).

Two kinds of published measured data are used, each independent of the
Pampa la Bola FTS slopes of meas_check.py:

1. Zenith opacity ratios against a 225 GHz tipper at Mauna Kea (4.1 km),
   from the skydip calibrations of the SCUBA and SCUBA-2 cameras on the JCMT
   (DATA["zenith"]). Each reports tau_band = a (tau_225 - b); the slope a is
   a measured ratio of water-vapor opacities. The camera bands are broad, so
   each model gives the slope of the band-effective opacity,
   tau_eff = -ln(mean_band exp(-m tau(f))) / m, at airmass m, between two
   water amounts, divided by the change of tau(225 GHz).

2. Horizontal-path specific attenuation measured at stated temperature,
   pressure and humidity (DATA["horizontal"]); each model gives gamma (dB/km)
   at the same state.

Models: P.676-13 Annex 1 (the paper's reference), am 14.0 and the PyRTlib
models. Output: meas_check2.json and meas2_macros.tex (macros Mtwo*).
"""
import json, os, re, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, SIM)
import thzatm as A  # noqa: E402

AM = os.environ.get("AM_BIN", "am")
os.environ.setdefault("AM_CACHE_PATH", os.path.join(HERE, "am_cache"))
os.makedirs(os.environ["AM_CACHE_PATH"], exist_ok=True)
DB = 10 * np.log10(np.e)
MODELS = [("R98", "R98", "R98"), ("R17", "R17", "R17"), ("R22", "R22SD", "R22"),
          ("R24", "R24", "R24"), ("MWL24", "MWL24", "R24")]
DATA = json.load(open(os.path.join(HERE, "meas_data.json")))
SCALES = (0.5, 1.0)


def atm(pr, s, h0):
    h, d = A.annex1_layers_from(h0)
    hm = h + 0.5 * d
    if pr == "p835-std":
        T, P, rho = A.state(hm, rho0=7.5 * s)
        Tb, Pb, _ = A.state(h, rho0=7.5 * s)
    else:
        T, P, rho = A.p835_seasonal(hm, pr)
        Tb, Pb, _ = A.p835_seasonal(h, pr)
        rho = rho * s
    return h, d, hm, T, P, rho, Tb, Pb


_MEMO = {}


def zenith_tau(model, pr, s, h0, f):
    """Zenith opacity (nepers) on frequency grid f (GHz), memoized (the
    SCUBA and SCUBA-2 datasets share their bands)."""
    key = (model, pr, round(float(s), 10), h0, float(f[0]), float(f[-1]), len(f))
    if key not in _MEMO:
        _MEMO[key] = _zenith_tau(model, pr, s, h0, f)
    return _MEMO[key]


def _zenith_tau(model, pr, s, h0, f):
    h, d, hm, T, P, rho, Tb, Pb = atm(pr, s, h0)
    if model.startswith("P676"):
        return (d @ A.gamma_grid(f, hm, T, P, rho, **shortcut_kw(model, f))) / DB
    if model == "am":
        return tau_am(pr, s, h0, h, T, P, rho, Tb, Pb, f)
    lab = {m[0]: m for m in MODELS}[model]
    # PyRTlib takes one frequency per call; it runs on every third layer
    # (thickness summed), as the slant study runs it on every ninth
    n3 = len(d) // 3 * 3
    d3 = d[:n3].reshape(-1, 3).sum(1)
    k = np.arange(1, n3, 3)
    return np.array([d3 @ np.asarray(A.gamma_pyrtlib(fi, T[k], P[k], rho[k], lab[1], lab[2])).ravel()
                     for fi in np.atleast_1d(f)]) / DB


def shortcut_kw(model, f):
    """Line-set shortcuts of Sec. V-E applied to the P.676-13 kernel:
    P676np drops the 1780 GHz pseudo-line, P676loc keeps only lines within
    +-50 GHz of the band center (the dry continuum is kept in both)."""
    if model == "P676np":
        return {"h2o_mask": A.H2O[:, 0] != 1780.0}
    if model == "P676loc":
        fc = 0.5 * (float(np.min(f)) + float(np.max(f)))
        return {"o2_mask": np.abs(A.O2[:, 0] - fc) <= 50.0, "h2o_mask": np.abs(A.H2O[:, 0] - fc) <= 50.0}
    return {}


def tau_am(pr, s, h0, h, T, P, rho, Tb, Pb, f):
    lines = ["f %1 %2 %3 %4 %5 %6", "output f GHz tau", "za 0 deg", "tol 1e-5", "T0 2.7 K", ""]
    # every third Annex 1 layer, as for PyRTlib: the base of each group of
    # three and the state of its middle layer
    for i in range(0, len(h) // 3 * 3, 3)[::-1]:
        e = max(rho[i + 1], 0.0) * T[i + 1] / 216.7
        lines += ["layer", f"Pbase {Pb[i]:.6g} mbar", f"Tbase {Tb[i]:.5g} K", "column dry_air vmr",
                  f"column h2o vmr {e / P[i + 1]:.6e}", ""]
    pth = os.path.join(HERE, "slant_cfg", f"meas2_{pr}_{s:g}_{h0:g}.amc")
    os.makedirs(os.path.dirname(pth), exist_ok=True)
    open(pth, "w").write("\n".join(lines))
    df = max(float(np.min(np.diff(f))) if len(f) > 1 else 0.1, 0.01)
    r = subprocess.run([AM, pth, f"{f.min() - df:.4f}", "GHz", f"{f.max() + df:.4f}", "GHz",
                        f"{min(df, 0.5) * 1000:.3f}", "MHz"], capture_output=True, text=True)
    a = np.loadtxt(r.stdout.splitlines())
    return np.interp(f, a[:, 0], a[:, 1])


def band_grid(lo, hi, step=0.5):
    return np.arange(lo, hi + 1e-9, step)


def tau_eff(tau, m):
    return -np.log(np.mean(np.exp(-m * tau))) / m


def scale_for(pr, V, h0):
    """Water-vapor scale that gives column V (mm) above the site."""
    h, d, hm, T, P, rho, *_ = atm(pr, 1.0, h0)
    return V / A.column_mm(h, d, rho)


def zenith_slopes(model, ds):
    """Slope of the band-effective opacity against tau(225 GHz)."""
    out = {}
    h0 = ds["site_km"]
    for pr in ds["profiles"]:
        sc = [scale_for(pr, V, h0) for V in ds["V_mm"]]
        res = {}
        t225 = [float(zenith_tau(model, pr, s, h0, np.array([225.0]))[0]) for s in sc]
        for b in ds["bands"]:
            f = band_grid(b["lo_ghz"], b["hi_ghz"])
            tb = [tau_eff(zenith_tau(model, pr, s, h0, f), ds["airmass"]) for s in sc]
            res[b["name"]] = (tb[1] - tb[0]) / (t225[1] - t225[0])
        out[pr] = res
    return out


def pwv_slopes(model, ds):
    """Double-sideband zenith opacity against the column: dry term a (column
    -> 0) and slope b per mm, and the 225 GHz slope."""
    out = {}
    h0 = ds["site_km"]
    f = np.concatenate([band_grid(c - ds["if_bw_ghz"] / 2, c + ds["if_bw_ghz"] / 2, 0.05) for c in ds["dsb_ghz"]])
    for pr in ds["profiles"]:
        V = ds["V_mm"]
        sc = [scale_for(pr, v, h0) for v in V]
        # each sideband carries half the gain: average the transmissions
        t = [float(-np.log(np.mean(np.exp(-zenith_tau(model, pr, s, h0, f))))) for s in sc]
        t225 = [float(zenith_tau(model, pr, s, h0, np.array([225.0]))[0]) for s in sc]
        b = (t[1] - t[0]) / (V[1] - V[0])
        b225 = (t225[1] - t225[0]) / (V[1] - V[0])
        out[pr] = {"a": t[0] - b * V[0], "b": b, "b225": b225}
    return out


def horizontal_gamma(model, st, f, kind):
    """Specific attenuation (dB/km) on grid f; kind 'water' removes the dry-air part."""
    T = np.array([st["T_K"]])
    if "rho_gm3" in st:
        rho = np.array([st["rho_gm3"]])
    else:
        rho = np.array([st["rh_pct"] / 100 * float(A.sat_vapor_pressure(T)[0]) * 216.7 / st["T_K"]])
    P = np.array([st["P_hpa"]])
    e = float(rho[0] * T[0] / 216.7)

    def g(r):
        if model.startswith("P676"):
            return A.gamma_grid(f, np.array([0.0]), T, P, np.array([r]), **shortcut_kw(model, f))[0]
        if model == "am":
            return gamma_am_horiz(T[0], P[0], r * T[0] / 216.7, f, water_only=False)
        lab = {m[0]: m for m in MODELS}[model]
        return np.array([float(np.asarray(A.gamma_pyrtlib(fi, T, P, np.array([r]), lab[1], lab[2])).ravel()[0])
                         for fi in f])
    if model == "am" and kind == "water":
        return gamma_am_horiz(T[0], P[0], e, f, water_only=True)
    out = g(float(rho[0]))
    return out - g(0.0) if kind == "water" else out


def gamma_am_horiz(T, P, e, f, water_only):
    """am specific attenuation (dB/km) from a uniform 1 km horizontal layer."""
    lines = ["f %1 %2 %3 %4 %5 %6", "output f GHz tau", "za 0 deg", "tol 1e-5", "T0 2.7 K", "",
             "layer", f"P {P:.6g} mbar", f"T {T:.5g} K", "h 1000 m"]
    if not water_only:
        lines.append("column dry_air vmr")
    lines += [f"column h2o vmr {e / P:.6e}", ""]
    pth = os.path.join(HERE, "slant_cfg", f"meas2_h_{T:.1f}_{P:.1f}_{e:.3f}_{int(water_only)}.amc")
    os.makedirs(os.path.dirname(pth), exist_ok=True)
    open(pth, "w").write("\n".join(lines))
    r = subprocess.run([AM, pth, f"{f.min() - 0.1:.3f}", "GHz", f"{f.max() + 0.1:.3f}", "GHz", "50", "MHz"],
                       capture_output=True, text=True)
    a = np.loadtxt(r.stdout.splitlines())
    return np.interp(f, a[:, 0], a[:, 1]) * DB


ALL = ["P676", "am", "R98", "R17", "R22", "R24", "MWL24", "P676np", "P676loc"]


def run_model(mo):
    """One model over every dataset. Results are saved after each dataset,
    and a rerun resumes from the saved part (delete slant_cfg/meas2_<model>.json
    to recompute from scratch)."""
    pth = os.path.join(HERE, "slant_cfg", f"meas2_{mo}.json")
    res = json.load(open(pth)) if os.path.exists(pth) else {"zenith": {}, "pwv": {}, "horizontal": {}}

    def save():
        json.dump(res, open(pth, "w"), indent=1)
    for ds in DATA["zenith"]:
        if ds["id"] not in res["zenith"]:
            res["zenith"][ds["id"]] = zenith_slopes(mo, ds)
            save()
        print(mo, ds["id"], json.dumps(res["zenith"][ds["id"]]), flush=True)
    for ds in DATA["pwv"]:
        if ds["id"] not in res["pwv"]:
            res["pwv"][ds["id"]] = pwv_slopes(mo, ds)
            save()
        print(mo, ds["id"], json.dumps(res["pwv"][ds["id"]]), flush=True)
    for ds in DATA["horizontal"]:
        if ds["id"] not in res["horizontal"]:
            vals = []
            for pt in ds["points"]:
                f = band_grid(pt["lo_ghz"], pt["hi_ghz"], 0.1)
                vals.append(float(np.min(horizontal_gamma(mo, ds["state"], f, ds["kind"]))))
            res["horizontal"][ds["id"]] = vals
            save()
        print(mo, ds["id"], json.dumps(res["horizontal"][ds["id"]]), flush=True)


def merge():
    res = {"zenith": {}, "pwv": {}, "horizontal": {}}
    for mo in ALL:
        r = json.load(open(os.path.join(HERE, "slant_cfg", f"meas2_{mo}.json")))
        for g in res:
            for k, v in r[g].items():
                res[g].setdefault(k, {})[mo] = v
    json.dump(res, open(os.path.join(HERE, "meas_check2.json"), "w"), indent=1)


DS = {"archibald": "Arch", "dempsey": "Demp", "radford_mk": "RadMk", "radford_cj": "RadCj",
      "chamberlin": "Cham", "yang2014": "Yfourteen", "yang2011": "Yeleven"}
BAND = {"850um": "Eight", "450um": "Four", "350um": "Three", "350": "ThreeFifty", "410": "FourTen",
        "680": "SixEighty", "850": "EightFifty"}
NM = {"P676": "Pitu", "am": "Am", "R98": "RNinetyEight", "R17": "RSeventeen", "R22": "RTwentyTwo",
      "R24": "RTwentyFour", "MWL24": "Mwl", "P676np": "Nopseudo", "P676loc": "Local"}


def _rng(vals, fmt):
    lo, hi = min(vals), max(vals)
    a, b = (fmt.format(x).replace("-0", "+0") if fmt.format(x).lstrip("+-").strip("0.") == "" else fmt.format(x)
            for x in (lo, hi))
    return a if a == b else f"{a} to {b}"


def summarize():
    res = json.load(open(os.path.join(HERE, "meas_check2.json")))
    m = {}
    worst = {mo: {"low": [], "high": [], "hor": []} for mo in ALL}
    for ds in DATA["zenith"]:
        for b in ds["bands"]:
            key = f"Mtwo{DS[ds['id']]}{BAND[b['name']]}"
            m[key + "Meas"] = f"{b['slope']:g}"
            # stated uncertainty of the source; "n/s" where none is stated (round 3, item 6)
            m[key + "Err"] = f"{b['err']:g}" if b["err"] else "n/s"
            if b["err"]:
                m[key + "ErrPct"] = f"{100 * b['err'] / b['slope']:.0f}"
                for mo in ("P676",):
                    vals = [res["zenith"][ds["id"]][mo][pr][b["name"]] for pr in ds["profiles"]]
                    m[key + NM[mo] + "Sigma"] = _rng([(v - b["slope"]) / b["err"] for v in vals], "{:+.1f}")
            for mo in ALL:
                vals = [res["zenith"][ds["id"]][mo][pr][b["name"]] for pr in ds["profiles"]]
                rel = [100 * (v - b["slope"]) / b["slope"] for v in vals]
                m[key + NM[mo] + "Val"] = _rng(vals, "{:.1f}" if b["slope"] > 10 else "{:.2f}")
                m[key + NM[mo] + "Pct"] = _rng(rel, "{:+.0f}")
                worst[mo]["low" if b["hi_ghz"] < 400 else "high"] += rel
    for ds in DATA["pwv"]:
        key = f"Mtwo{DS[ds['id']]}"
        m[key + "AMeas"], m[key + "BMeas"] = f"{ds['a']:g}", f"{ds['b_per_mm']:g}"
        m[key + "BErr"], m[key + "AErr"] = f"{ds['b_err']:g}", f"{ds['a_err']:g}"
        m[key + "RatioMeas"] = f"{ds['b_per_mm'] / ds['b225_per_mm']:.1f}"
        for mo in ALL:
            r = res["pwv"][ds["id"]][mo]
            bs = [r[pr]["b"] for pr in ds["profiles"]]
            As = [r[pr]["a"] for pr in ds["profiles"]]
            rat = [r[pr]["b"] / r[pr]["b225"] for pr in ds["profiles"]]
            m[key + "B" + NM[mo] + "Val"] = _rng(bs, "{:.2f}")
            m[key + "B" + NM[mo] + "Pct"] = _rng([100 * (x - ds["b_per_mm"]) / ds["b_per_mm"] for x in bs], "{:+.0f}")
            m[key + "A" + NM[mo] + "Val"] = _rng(As, "{:.2f}")
            m[key + "Ratio" + NM[mo] + "Val"] = _rng(rat, "{:.1f}")
            worst[mo]["high"] += [100 * (x - ds["b_per_mm"]) / ds["b_per_mm"] for x in bs]
    for ds in DATA["horizontal"]:
        if "rho_err_gm3" in ds["state"]:
            m[f"Mtwo{DS[ds['id']]}RhoErr"] = f"{ds['state']['rho_err_gm3']:g}"
        for k, pt in enumerate(ds["points"]):
            m[f"Mtwo{DS[ds['id']]}{BAND[pt['name']]}Err"] = f"{pt['err']:g}" if pt["err"] else "n/s"
            key = f"Mtwo{DS[ds['id']]}{BAND[pt['name']]}"
            m[key + "Meas"] = f"{pt['meas']:g}"
            for mo in ALL:
                v = res["horizontal"][ds["id"]][mo][k]
                m[key + NM[mo] + "Val"] = f"{v:.1f}"
                m[key + NM[mo] + "Pct"] = f"{100 * (v - pt['meas']) / pt['meas']:+.0f}"
                if ds["id"] == "yang2014":
                    worst[mo]["hor"].append(100 * (v - pt["meas"]) / pt["meas"])
    for mo in ALL:
        for g in ("low", "high", "hor"):
            if worst[mo][g]:
                m[f"Mtwo{NM[mo]}{g.capitalize()}Range"] = _rng(worst[mo][g], "{:+.0f}")
                m[f"Mtwo{NM[mo]}{g.capitalize()}AbsMax"] = f"{max(abs(x) for x in worst[mo][g]):.0f}"
    # sensitivity of the band-effective slopes to the assumed airmass (P.676)
    import copy
    spread = {}
    for ds0 in DATA["zenith"]:
        if ds0["id"] == "dempsey":
            continue
        base = zenith_slopes("P676", ds0)
        for mass in (1.0, 2.5):
            ds = copy.deepcopy(ds0)
            ds["airmass"] = mass
            alt = zenith_slopes("P676", ds)
            for pr in ds["profiles"]:
                for b in ds["bands"]:
                    k = "Eight" if b["name"] == "850um" else "Other"
                    spread.setdefault(k, []).append(abs(alt[pr][b["name"]] / base[pr][b["name"]] - 1) * 100)
    m["MtwoAirmassEightPct"] = f"{max(spread['Eight']):.0f}"
    m["MtwoAirmassOtherPct"] = f"{max(spread['Other']):.0f}"
    # cross-dataset summaries quoted in the text. The FTS slopes of
    # meas_check.py (Matsushita and Matsuo 2003) are included.
    fts = json.load(open(os.path.join(HERE, "meas_check.json")))
    FMEAS = {"345.0": 3.65, "410.0": 7.53, "675.0": 22.4, "875.0": 24.2}
    FPROF = ["mid-winter", "p835-std", "high-winter"]

    def fts_rel(mo, f):
        key = {"R22": "R22", "MWL24": "MWL24", "P676": "P676", "am": "am"}[mo]
        return [100 * (fts[pr][key][f] - FMEAS[f]) / FMEAS[f] for pr in FPROF]

    def zrel(mo, did, band):
        ds = next(x for x in DATA["zenith"] if x["id"] == did)
        b = next(x for x in ds["bands"] if x["name"] == band)
        return [100 * (res["zenith"][did][mo][pr][band] - b["slope"]) / b["slope"] for pr in ds["profiles"]]

    def hrel(mo, did, k):
        ds = next(x for x in DATA["horizontal"] if x["id"] == did)
        return [100 * (res["horizontal"][did][mo][k] - ds["points"][k]["meas"]) / ds["points"][k]["meas"]]
    for mo in ("P676", "am", "R22", "MWL24"):
        hi = zrel(mo, "radford_mk", "350um") + zrel(mo, "radford_cj", "350um") + fts_rel(mo, "875.0") + hrel(mo, "yang2014", 2)
        m[f"MtwoEightFifty{NM[mo]}MinPct"] = f"{min(hi):.0f}"
        m[f"MtwoEightFifty{NM[mo]}MaxPct"] = f"{max(hi):.0f}"
        m[f"MtwoEightFifty{NM[mo]}AbsMax"] = f"{max(abs(x) for x in hi):.0f}"
        lo = (fts_rel(mo, "345.0") + fts_rel(mo, "410.0") + zrel(mo, "dempsey", "850um")
              + hrel(mo, "yang2014", 0))
        m[f"MtwoLowExArch{NM[mo]}AbsMax"] = f"{max(abs(x) for x in lo):.0f}"
    # scatter of the data around P.676-13 at the 340, 410 and 670 GHz windows
    scat = (fts_rel("P676", "345.0") + fts_rel("P676", "410.0") + fts_rel("P676", "675.0")
            + zrel("P676", "archibald", "850um") + zrel("P676", "dempsey", "850um")
            + zrel("P676", "archibald", "450um") + zrel("P676", "dempsey", "450um")
            + hrel("P676", "yang2014", 0) + hrel("P676", "yang2014", 1))
    m["MtwoWinScatterMaxPct"] = f"{max(abs(x) for x in scat):.0f}"
    sa, sd = (next(b["slope"] for b in next(x for x in DATA["zenith"] if x["id"] == i)["bands"]
                   if b["name"] == "850um") for i in ("archibald", "dempsey"))
    tp = [100 * b["err"] / b["slope"] for x in DATA["zenith"] if x["id"].startswith("radford") for b in x["bands"]]
    m["MtwoTipperErrMinPct"], m["MtwoTipperErrMaxPct"] = f"{min(tp):.0f}", f"{max(tp):.0f}"
    m["MtwoScubaEightDiffPct"] = f"{100 * (sd - sa) / sa:.0f}"
    refs = {re.match(r"^(.*?\d{4})", ds["ref"]).group(1) for g in ("zenith", "pwv", "horizontal") for ds in DATA[g]}
    m["MtwoNumStudies"] = str(len(refs) + 1)  # plus Matsushita and Matsuo (meas_check.py)
    m["MtwoNumDatasets"] = str(len(DATA["zenith"]) + len(DATA["pwv"]) + len(DATA["horizontal"]))
    with open(os.path.join(SIM, "..", "meas2_macros.tex"), "w") as fh:
        fh.write("% Generated by sim/am_check/meas_check2.py. Do not edit by hand.\n")
        for k in sorted(m):
            fh.write(f"\\newcommand{{\\{k}}}{{{A.tex_minus(m[k])}}}\n")
    print(json.dumps(m, indent=1))


if __name__ == "__main__":
    # --model M computes one model (the run_all.sh stage runs them in
    # parallel); without arguments the cached per-model results are merged
    # and summarized; --macros-only only summarizes meas_check2.json
    if "--model" in sys.argv:
        run_model(sys.argv[sys.argv.index("--model") + 1])
    else:
        if "--macros-only" not in sys.argv:
            merge()
        summarize()
