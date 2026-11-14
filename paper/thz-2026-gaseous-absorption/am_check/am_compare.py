"""Independent line-by-line check of the P.676-13 Annex 1 reference against
am 14.0 (Paine; HITRAN-based line catalog, MT_CKD continuum).

Both codes see the same atmosphere: the reference engine's P.835 profile
(thzatm.state) discretized into the same layers, with dry air and water vapor
only. Compared quantity: zenith gaseous loss (dB) at each window center.
"""
import json, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import thzatm as A  # noqa: E402

AM = os.environ.get("AM_BIN", "am")
WIN = json.load(open(os.path.join(os.path.dirname(HERE), "out", "windows.json")))["centers"]


def am_config(T, P, rho, h_top_first, parts=("dry", "wet")):
    """Layers top-down: am wants Pbase/Tbase per layer, with vmr at mid-layer.
    parts selects dry air (O2, N2/O2 collision-induced absorption and am's
    trace gases), water vapor (lines plus MT_CKD continuum), or both."""
    lines = ["f %1 %2 %3 %4 %5 %6", "output f GHz tau", "za 0 deg", "tol 1e-5", "T0 2.7 K", ""]
    for Pi, Ti, ri in h_top_first:
        e = max(ri, 0.0) * Ti / 216.7
        vmr = e / Pi
        lines += ["layer", f"Pbase {Pi:.6g} mbar", f"Tbase {Ti:.5g} K"]
        if "dry" in parts:
            lines += ["column dry_air vmr"]
        if "wet" in parts:
            lines += [f"column h2o vmr {vmr:.6e}"]
        lines += [""]
    return "\n".join(lines)


def am_tau(cfg_path, fc):
    r = subprocess.run([AM, cfg_path, f"{fc-0.05}", "GHz", f"{fc+0.05}", "GHz", "10", "MHz"],
                       capture_output=True, text=True)
    arr = np.loadtxt(r.stdout.splitlines())
    return float(np.interp(fc, arr[:, 0], arr[:, 1]))


def main():
    # SKIP_AM=1 keeps the am 14.0 outputs of am_compare.json and recomputes
    # only the P.676 side (see run_all.sh)
    cached = json.load(open(os.path.join(HERE, "am_compare.json"))) if os.environ.get("SKIP_AM") == "1" else None
    h, d = A.annex1_layers()  # the full 922-layer P.676 grid, 0-100 km
    hm = h + 0.5 * d
    out = {}
    for pwv in (1.0, 2.0, 5.0, 10.0):
        sc, prof = A.scale_for_pwv(pwv)
        T, P, rho = A.atm_state(hm, prof, sc)
        Tb, Pb, _ = A.atm_state(h, prof, sc)  # base (bottom) of each layer
        # am layer i: base at the bottom pressure, mid-layer humidity
        layers = [(Pb[i], Tb[i], rho[i]) for i in range(len(h))][::-1]
        db = 10 * np.log10(np.e)
        paths = {}
        for part, sel in (("all", ("dry", "wet")), ("dry", ("dry",)), ("wet", ("wet",))):
            pth = os.path.join(HERE, f"ref_{pwv:g}mm_{part}.amc")
            open(pth, "w").write(am_config(T, P, rho, layers, sel))
            paths[part] = pth
        a_z = A.annex1_path_lengths(90.0, h, d)
        ref_zen = a_z @ A.gamma_grid(np.array(WIN), hm, T, P, rho)
        gd = np.stack([A.gamma_p676(np.array(WIN), T[k], P[k], rho[k], split=True)[0] for k in range(len(hm))])
        ref_dry = a_z @ gd
        row = []
        for k_, (fc, rz, rd) in enumerate(zip(WIN, ref_zen, ref_dry)):
            if cached is not None:
                c_ = cached[str(pwv)][k_]
                assert abs(c_["f"] - fc) < 1e-9
                am_db, am_dry, am_wet = c_["am_db"], c_["am_dry"], c_["am_wet"]
            else:
                am_db = db * am_tau(paths["all"], fc)
                am_dry = db * am_tau(paths["dry"], fc)
                am_wet = db * am_tau(paths["wet"], fc)
            row.append({"f": fc, "ref_db": float(rz), "am_db": am_db, "diff_db": float(rz - am_db),
                        "rel": float((rz - am_db) / am_db), "ref_dry": float(rd), "ref_wet": float(rz - rd),
                        "am_dry": am_dry, "am_wet": am_wet})
        out[pwv] = row
        print(pwv, [(x["f"], round(x["ref_db"], 2), round(x["am_db"], 2)) for x in row], flush=True)
    json.dump(out, open(os.path.join(HERE, "am_compare.json"), "w"), indent=1)




def write_macros():
    """am_macros.tex: agreement of the reference with am, per window."""
    d = json.load(open(os.path.join(HERE, "am_compare.json")))
    names = dict(zip(WIN, ["ThreeForty", "FourTen", "FourSixty", "SixSeventy", "EightFifty"]))
    m = {}
    for pwv, rows in d.items():
        tg = {"1.0": "One", "2.0": "Two", "5.0": "Five", "10.0": "Ten"}[pwv]
        for r in rows:
            w = names[r["f"]]
            m[f"AmDb{w}{tg}Mm"] = f"{r['am_db']:.1f}"
            m[f"AmRefDb{w}{tg}Mm"] = f"{r['ref_db']:.1f}"
            m[f"AmDiffPct{w}{tg}Mm"] = f"{100 * r['rel']:+.0f}"
            m[f"AmDiffDb{w}{tg}Mm"] = f"{r['diff_db']:+.1f}"
    lo = [r["rel"] for rows in d.values() for r in rows if r["f"] < 500]
    hi = [r["rel"] for rows in d.values() for r in rows if r["f"] > 500]
    m["AmLowMinPct"] = f"{100 * min(lo):.0f}"
    m["AmLowMaxPct"] = f"{100 * max(lo):.0f}"
    m["AmHighMinPct"] = f"{100 * min(hi):.0f}"
    m["AmHighMaxPct"] = f"{100 * max(hi):.0f}"
    dr = [(r["ref_dry"] - r["am_dry"]) / r["am_dry"] for rows in d.values() for r in rows]
    wr = [(r["ref_wet"] - r["am_wet"]) / r["am_wet"] for rows in d.values() for r in rows]
    m["AmDryMinPct"] = f"{100 * min(dr):.0f}"
    m["AmDryBelowLoPct"] = f"{-100 * max(dr):.0f}"
    m["AmDryBelowHiPct"] = f"{-100 * min(dr):.0f}"
    m["AmDryMaxPct"] = f"{100 * max(dr):.0f}"
    m["AmWetMinPct"] = f"{100 * min(wr):.0f}"
    m["AmWetMaxPct"] = f"{100 * max(wr):.0f}"
    m["AmVersion"] = "14.0"
    with open(os.path.join(os.path.dirname(HERE), "..", "am_macros.tex"), "w") as f:
        f.write("% Generated by sim/am_check/am_compare.py. Do not edit by hand.\n")
        for k in sorted(m):
            f.write(f"\\newcommand{{\\{k}}}{{{A.tex_minus(m[k])}}}\n")
    print("wrote", len(m), "am macros")


if __name__ == "__main__":
    if "--macros-only" not in sys.argv:
        main()
    write_macros()
