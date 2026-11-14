"""Comparison with measurement (round 2, review D13).

Matsushita and Matsuo (PASJ 55, 325, 2003, Table 1) report, from Fourier
transform spectrometer (FTS) zenith opacities at Pampa la Bola (4800 m,
northern Chile), the slope a_WV = d tau(nu) / d tau(220 GHz) of the water-vapor
("normal") correlation at 345, 410, 492, 675 and 875 GHz. The slope is a
measured ratio of water-vapor opacities and needs no column calibration.

Each model gives the same ratio from two water amounts on a zenith path
from 4.8 km through a P.835-6 atmosphere: P.676-13 Annex 1 (the paper's
reference), am 14.0, and the PyRTlib models. Output: meas_check.json and
meas_macros.tex (macros Meas*).
"""
import json, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, SIM)
import thzatm as A  # noqa: E402

AM = os.environ.get("AM_BIN", "am")
os.environ.setdefault("AM_CACHE_PATH", os.path.join(HERE, "am_cache"))
os.makedirs(os.environ["AM_CACHE_PATH"], exist_ok=True)
FREQ = [220.0, 345.0, 410.0, 492.0, 675.0, 875.0]
MEAS = {345.0: (3.65, 0.02), 410.0: (7.53, 0.05), 492.0: (23.6, 0.3), 675.0: (22.4, 0.2),
        875.0: (24.2, 0.2)}
H0 = 4.8
PROFS = ["mid-winter", "p835-std", "high-winter"]
SCALES = (0.5, 1.0)
DB = 10 * np.log10(np.e)
MODELS = [("R98", "R98", "R98"), ("R17", "R17", "R17"), ("R22", "R22SD", "R22"),
          ("R24", "R24", "R24"), ("MWL24", "MWL24", "R24")]


def atm(pr, s):
    h, d = A.annex1_layers_from(H0)
    hm = h + 0.5 * d
    if pr == "p835-std":
        T, P, rho = A.state(hm, rho0=7.5 * s)
        Tb, Pb, _ = A.state(h, rho0=7.5 * s)
    else:
        T, P, rho = A.p835_seasonal(hm, pr)
        Tb, Pb, _ = A.p835_seasonal(h, pr)
        rho = rho * s
    return h, d, hm, T, P, rho, Tb, Pb


def tau_p676(pr, s):
    h, d, hm, T, P, rho, *_ = atm(pr, s)
    return (d @ A.gamma_grid(np.array(FREQ), hm, T, P, rho)) / DB  # nepers (opacity)


def tau_pyrtlib(pr, s, h2o, o2):
    h, d, hm, T, P, rho, *_ = atm(pr, s)
    return np.array([d @ A.gamma_pyrtlib(f, T, P, rho, h2o, o2).ravel() for f in FREQ]) / DB


def tau_am(pr, s):
    h, d, hm, T, P, rho, Tb, Pb = atm(pr, s)
    lines = ["f %1 %2 %3 %4 %5 %6", "output f GHz tau", "za 0 deg", "tol 1e-5", "T0 2.7 K", ""]
    for i in range(len(h))[::-1]:
        e = max(rho[i], 0.0) * T[i] / 216.7
        lines += ["layer", f"Pbase {Pb[i]:.6g} mbar", f"Tbase {Tb[i]:.5g} K", "column dry_air vmr",
                  f"column h2o vmr {e / P[i]:.6e}", ""]
    pth = os.path.join(HERE, "slant_cfg", f"meas_{pr}_{s:g}.amc")
    os.makedirs(os.path.dirname(pth), exist_ok=True)
    open(pth, "w").write("\n".join(lines))
    out = []
    for f in FREQ:
        r = subprocess.run([AM, pth, f"{f - 0.05}", "GHz", f"{f + 0.05}", "GHz", "10", "MHz"],
                           capture_output=True, text=True)
        a = np.loadtxt(r.stdout.splitlines())
        out.append(float(np.interp(f, a[:, 0], a[:, 1])))
    return np.array(out)


def ratio(t_lo, t_hi):
    dt = t_hi - t_lo
    return {f: float(dt[k] / dt[0]) for k, f in enumerate(FREQ) if f in MEAS}


def main():
    res = {}
    # SKIP_AM=1 keeps the am 14.0 ratios of meas_check.json (see run_all.sh)
    cached = json.load(open(os.path.join(HERE, "meas_check.json"))) if os.environ.get("SKIP_AM") == "1" else None
    for pr in PROFS:
        h, d, hm, T, P, rho, *_ = atm(pr, 1.0)
        res[pr] = {"V": A.column_mm(h, d, rho)}
        res[pr]["P676"] = ratio(tau_p676(pr, SCALES[0]), tau_p676(pr, SCALES[1]))
        res[pr]["am"] = cached[pr]["am"] if cached else ratio(tau_am(pr, SCALES[0]), tau_am(pr, SCALES[1]))
        for lab, h2o, o2 in MODELS:
            res[pr][lab] = ratio(tau_pyrtlib(pr, SCALES[0], h2o, o2), tau_pyrtlib(pr, SCALES[1], h2o, o2))
        print(pr, json.dumps(res[pr]), flush=True)
    json.dump(res, open(os.path.join(HERE, "meas_check.json"), "w"), indent=1)


def summarize():
    res = json.load(open(os.path.join(HERE, "meas_check.json")))
    m = {}
    tags = {345.0: "ThreeFortyFive", 410.0: "FourTen", 492.0: "FourNinetyTwo", 675.0: "SixSeventyFive",
            875.0: "EightSeventyFive"}
    models = ["P676", "am", "R98", "R17", "R22", "R24", "MWL24"]
    NM = {"P676": "Pitu", "am": "Am", "R98": "RNinetyEight", "R17": "RSeventeen", "R22": "RTwentyTwo",
          "R24": "RTwentyFour", "MWL24": "Mwl"}
    for f, (a, s) in MEAS.items():
        m[f"MeasA{tags[f]}"] = f"{a:g}"
        m[f"MeasErr{tags[f]}"] = f"{s:g}"  # stated uncertainty, Table 1 of the source (round 3, item 6)
        for mo in models:
            vals = [res[pr][mo][str(f)] for pr in PROFS]
            rel = [100 * (v - a) / a for v in vals]
            lo, hi = round(min(rel)), round(max(rel))
            m[f"Meas{NM[mo]}{tags[f]}Pct"] = f"{lo:+d}" if lo == hi else f"{lo:+d} to {hi:+d}"
            m[f"Meas{NM[mo]}{tags[f]}Lo"] = f"{min(rel):+.0f}"
            m[f"Meas{NM[mo]}{tags[f]}Hi"] = f"{max(rel):+.0f}"
    # worst absolute departure of each model from the measured slopes above 600 GHz
    for mo in models:
        w = max(abs(100 * ((res[pr][mo][str(f)]) - MEAS[f][0]) / MEAS[f][0]) for pr in PROFS for f in (675.0, 875.0))
        wl = max(abs(100 * ((res[pr][mo][str(f)]) - MEAS[f][0]) / MEAS[f][0]) for pr in PROFS for f in (345.0, 410.0))
        m[f"Meas{NM[mo]}HighWorstPct"] = f"{w:.0f}"
        m[f"Meas{NM[mo]}LowWorstPct"] = f"{wl:.0f}"
    m["MeasSiteKm"] = f"{H0:g}"
    with open(os.path.join(SIM, "..", "meas_macros.tex"), "w") as fh:
        fh.write("% Generated by sim/am_check/meas_check.py. Do not edit by hand.\n")
        for k in sorted(m):
            fh.write(f"\\newcommand{{\\{k}}}{{{A.tex_minus(m[k])}}}\n")
    print(json.dumps(m, indent=1))


if __name__ == "__main__":
    if "--macros-only" not in sys.argv:
        main()
    summarize()
