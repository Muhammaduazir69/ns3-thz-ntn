"""am 14.0 as a second reference on slant paths (round 2, review D7).

For each P.835-6 atmosphere at the driest and wettest realizable columns of
the six-atmosphere study (out/atm_meta.json), am computes the slant opacity
in spherical geometry with radio refraction, launched at the apparent zenith
angle of the P.676 ray, and the zenith opacity. The shortcuts are then scored
against am instead of P.676-13 Annex 1:
  cosecant law      am zenith / sin(theta) - am slant
  Annex 2 forms     shortcut - am slant (P.676-13 at 341.5 GHz; P.676-12 forms)
  column-only       am on the Sec. 1.1 profile at the same column - am slant
Output: am_slant.json and macros AmSl* in am_slant_macros.tex.
"""
import json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, SIM)
import thzatm as A  # noqa: E402

AM = os.environ.get("AM_BIN", "am")
OUT = os.path.join(SIM, "out")
WIN = json.load(open(os.path.join(OUT, "windows.json")))["centers"]
ELS = [10, 30, 90]
DB = 10 * np.log10(np.e)
WORK = os.path.join(HERE, "slant_cfg")
AMAX = 40.0
STEP = 6
os.environ.setdefault("AM_CACHE_PATH", os.path.join(HERE, "am_cache"))
os.makedirs(os.environ["AM_CACHE_PATH"], exist_ok=True)


def config(h, d, T, P, rho, za):
    hm = h + 0.5 * d
    lines = ["f %1 %2 %3 %4 %5 %6", "output f GHz tau", f"za {za:.6f} deg", "tol 1e-5",
             "T0 2.7 K", "geometry spherical", "refract radio", ""]
    Tb, Pb = T["base"], P["base"]
    for i in range(len(h))[::-1]:
        e = max(rho[i], 0.0) * T["mid"][i] / 216.7
        lines += ["layer", f"Pbase {Pb[i]:.6g} mbar", f"Tbase {Tb[i]:.5g} K",
                  "column dry_air vmr", f"column h2o vmr {e / P['mid'][i]:.6e}", ""]
    return "\n".join(lines)


def run(path, fc):
    r = subprocess.run([AM, path, f"{fc - 0.05}", "GHz", f"{fc + 0.05}", "GHz", "10", "MHz"],
                       capture_output=True, text=True)
    arr = np.loadtxt(r.stdout.splitlines())
    return DB * float(np.interp(fc, arr[:, 0], arr[:, 1]))


def am_case(tag, pr, V):
    """am slant loss (dB) at ELS x WIN for atmosphere pr scaled to column V."""
    h, d = A.annex1_layers()
    n = A.refractive_index(*A.state(h + 0.5 * d, pwv_mm=V, profile=pr))
    za_of = {e: 90.0 - A.apparent_elevation(e, h, d, n) for e in ELS}
    # am runs on every STEP-th Annex 1 boundary (cost scales with the layer
    # count); the P.676 kernel on the same coarse grid moves by under 1 percent.
    idx = np.arange(0, len(h), STEP)
    top = h[-1] + d[-1]
    h = h[idx]
    d = np.append(h[1:], top) - h
    hm = h + 0.5 * d
    T, P, rho = A.state(hm, pwv_mm=V, profile=pr)
    Tb, Pb, _ = A.state(h, pwv_mm=V, profile=pr)
    out = np.zeros((len(ELS), len(WIN)))
    jobs = []
    for j, e in enumerate(ELS):
        za = za_of[e]
        pth = os.path.join(WORK, f"{tag}_{e}.amc")
        open(pth, "w").write(config(h, d, {"base": Tb, "mid": T}, {"base": Pb, "mid": P}, rho, za))
        for k, fc in enumerate(WIN):
            jobs.append((j, k, pth, fc))
    with ThreadPoolExecutor(max_workers=int(os.environ.get("AM_WORKERS", "2"))) as ex:
        for (j, k, _, _), v in zip(jobs, ex.map(lambda a: run(a[2], a[3]), jobs)):
            out[j, k] = v
    return out


def main():
    os.makedirs(WORK, exist_ok=True)
    meta = json.load(open(os.path.join(OUT, "atm_meta.json")))
    Z = np.load(os.path.join(OUT, "atm.npz"))
    sys.path.insert(0, SIM)
    import analysis as AN
    floor = AN.RH_FLOOR
    cases = []
    for pr in AN.ATM_PROFS:
        cols = sorted(float(k.split("|")[1]) for k, m in meta.items()
                      if not k.startswith("native|") and k.split("|")[0] == pr and k.split("|")[2] == "0"
                      and m["rh_s_sea"] >= floor and m["rh_max"] <= 1.0)
        for V in sorted({cols[0], cols[-1]}):
            cases.append((pr, V))
    res = {}
    std_cache = {}
    # SKIP_AM=1 keeps the am 14.0 outputs of am_slant.json and recomputes only
    # the P.676 side (see run_all.sh)
    cached = json.load(open(os.path.join(HERE, "am_slant.json"))) if os.environ.get("SKIP_AM") == "1" else None
    for pr, V in cases:
        key = f"{pr}|{V:g}"
        if cached is not None:
            am = np.array(cached[key]["am"])
            std_cache[V] = np.array(cached[key]["am_std"])
        else:
            am = am_case(key.replace("|", "_"), pr, V)
            if V not in std_cache:
                std_cache[V] = am if pr == "p835-std" else am_case(f"std_{V:g}", "p835-std", V)
        X = Z[f"{pr}|{V:g}|0"]
        ei = [AN.ELEVS.index(e) for e in ELS]
        res[key] = {"am": am.tolist(), "am_std": std_cache[V].tolist(),
                    "p676": X[0][ei].tolist(), "csc": X[3][ei].tolist(), "a2": X[4][ei].tolist(),
                    "a2v": X[5][ei].tolist(), "a2r": X[6][ei].tolist(), "a2vr": X[7][ei].tolist()}
        # P.676 reference on the Sec. 1.1 profile at the same column (column-only shortcut)
        Xs = Z[f"p835-std|{V:g}|0"] if f"p835-std|{V:g}|0" in Z.files else None
        if Xs is None:
            h, d = A.annex1_layers()
            Xs, _ = AN.atm_case("p835-std", V, 0.0, np.array(WIN))
        res[key]["p676_std"] = Xs[0][ei].tolist()
        print(key, np.round(am[:, 0], 2), np.round(X[0][ei][:, 0], 2), flush=True)
    json.dump(res, open(os.path.join(HERE, "am_slant.json"), "w"), indent=1)


def summarize():
    d = json.load(open(os.path.join(HERE, "am_slant.json")))
    rows = {"Csc": [], "AnnexTwo": [], "AnnexTwoVt": [], "AnnexTwoRtw": [], "AnnexTwoVtRtw": [], "Prof": []}
    rows_p = {k: [] for k in rows}
    refdiff = []
    for key, r in d.items():
        am, am_std, p676 = (np.array(r[k]) for k in ("am", "am_std", "p676"))
        zen_am = am[-1]
        csc_am = zen_am[None, :] / np.sin(np.radians(ELS))[:, None]
        dom_am, dom_p = am <= AMAX, p676 <= AMAX
        sh = {"Csc": (csc_am, np.array(r["csc"])), "AnnexTwo": (np.array(r["a2"]),) * 2,
              "AnnexTwoVt": (np.array(r["a2v"]),) * 2, "AnnexTwoRtw": (np.array(r["a2r"]),) * 2,
              "AnnexTwoVtRtw": (np.array(r["a2vr"]),) * 2, "Prof": (am_std, np.array(r["p676_std"]))}
        for k, (s_am, s_p) in sh.items():
            ok = dom_am & np.isfinite(s_am)
            rows[k] += list((s_am - am)[ok])
            ok2 = dom_p & np.isfinite(s_p)
            rows_p[k] += list((s_p - p676)[ok2])
        refdiff += list(((p676 - am) / am)[dom_am])
    m = {}
    for k in rows:
        m[f"AmSl{k}WorstDb"] = np.max(np.abs(rows[k])) if rows[k] else float("nan")
        m[f"AmSl{k}WorstDbPsix"] = np.max(np.abs(rows_p[k])) if rows_p[k] else float("nan")
    order_am = sorted(rows, key=lambda k: m[f"AmSl{k}WorstDb"])
    order_p = sorted(rows, key=lambda k: m[f"AmSl{k}WorstDbPsix"])
    m["AmSlRankSame"] = "the same" if order_am == order_p else "different"
    m["AmSlCases"] = len(d)
    m["AmSlRefDiffMinPct"] = 100 * min(refdiff)
    m["AmSlRefDiffMaxPct"] = 100 * max(refdiff)
    out = {}
    for k, v in m.items():
        if isinstance(v, str):
            out[k] = v
        elif isinstance(v, int):
            out[k] = f"{v:d}"
        elif k.startswith("AmSlRefDiff"):
            # a range "above am": unsigned, and a rounded -0 prints as 0
            out[k] = f"{v:.0f}".replace("-0", "0") if abs(v) < 0.5 else f"{v:.0f}"
        elif "Pct" in k:
            out[k] = f"{v:+.0f}"
        else:
            out[k] = f"{v:.2f}" if v < 1 else f"{v:.1f}"
    with open(os.path.join(SIM, "..", "am_slant_macros.tex"), "w") as f:
        f.write("% Generated by sim/am_check/am_slant.py. Do not edit by hand.\n")
        for k in sorted(out):
            f.write(f"\\newcommand{{\\{k}}}{{{A.tex_minus(out[k])}}}\n")
    json.dump({"order_am": order_am, "order_p676": order_p, **out},
              open(os.path.join(HERE, "am_slant_summary.json"), "w"), indent=1)
    print(json.dumps({"order_am": order_am, "order_p676": order_p, **out}, indent=1))


if __name__ == "__main__":
    if "--macros-only" not in sys.argv:
        main()
    summarize()
