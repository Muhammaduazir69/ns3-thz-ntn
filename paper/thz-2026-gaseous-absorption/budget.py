#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""Per-window deviation budget against three references (T-TST paper, Table VII).

The loss of each shortcut at the 670 and 850 GHz window centers is scored against P.676-13 Annex 1,
against R22 and against am 14.0, because the published data put Annex 1 itself 16 to 36 percent high at
850 GHz. Geometric and profile shortcuts are applied to each reference's own model; the P.676-12 forms
and the no-pseudo-line set are fixed losses, so each is compared directly with that reference's loss.
Domain, every column: Annex 1 at most 40 dB at sea level on the realizable atmospheres; R22 on every
ninth layer; am on the driest and wettest realizable column of each
atmosphere and 10, 30, 90 degrees (am_check/am_slant.json, every sixth layer).
Also: worst Annex 1 deviation per shortcut in all five windows, the cosecant law at 5 degrees, the
total-pressure variant of the ITU-Rpy score, and whether each shortcut clears 1 dB everywhere.
Writes budget_macros.tex and out/budget.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analysis as AN  # noqa: E402
import thzatm as A  # noqa: E402

OUT = os.path.join(HERE, "out")
ELEVS = AN.ELEVS
AMAX = AN.AMAX
WN = ["ThreeForty", "FourTen", "FourSixty", "SixSeventy", "EightFifty"]
ROWS = ["Csc", "Prof", "Ptwh", "Ptwc", "Pseudo"]


def worst(vals):
    """Signed value of largest magnitude."""
    v = np.asarray(vals, float)
    return float(v[int(np.argmax(np.abs(v)))]) if v.size else float("nan")


def main():
    Z = np.load(os.path.join(OUT, "atm.npz"))
    L = np.load(os.path.join(OUT, "lines.npz"))
    M = np.load(os.path.join(OUT, "models.npz"))
    MS = np.load(os.path.join(OUT, "models_slant.npz"))
    SL = json.load(open(os.path.join(HERE, "am_check", "am_slant.json")))
    cs = np.array(json.load(open(os.path.join(OUT, "windows.json")))["centers"])
    ip = int(np.nanargmin(np.abs(L["f0"] - 1780.0)))
    cases = [(c.split("|")[0], float(c.split("|")[1])) for c in M["cases"]]
    h, d = A.annex1_layers()
    hm = h + 0.5 * d

    def pseudo(pr, p, e, k):
        return float(L[f"{AN.WINDOWS[k][0]}_{pr}_{p:g}_{e}"][ip])

    # Section 1.1 profile at the same column, Annex 1 (full grid)
    std_memo = {}

    def std_annex1(V):
        if V not in std_memo:
            T, P, rho = A.state(hm, pwv_mm=V)
            n = A.refractive_index(T, P, rho)
            G = A.gamma_grid(cs, hm, T, P, rho)
            std_memo[V] = {e: A.annex1_paths_geo(e, h, d, n) @ G for e in ELEVS}
        return std_memo[V]

    # Section 1.1 profile at the same column, R22 (same coarse layers as models.npz)
    hl, dl = M["h_lo"], M["delta"]
    hc = hl + 0.5 * dl
    r22_memo = {}

    def std_r22(V):
        if V not in r22_memo:
            import run_models2 as RM
            RM.setup("R22SD", "R22")
            T, P, rho = A.state(hc, pwv_mm=V)
            n = A.refractive_index(T, P, rho)
            G = RM.grid(cs[3:], T, P, rho)
            r22_memo[V] = {e: A.annex1_paths_geo(e, hl, dl, n) @ G for e in ELEVS}
        return r22_memo[V]

    # one domain for every column: Annex 1 at most 40 dB, the operable domain of the paper
    dev = {r: {k: [] for k in range(5)} for r in ROWS}
    devr = {r: {k: [] for k in (3, 4)} for r in ROWS}
    refr = {k: [] for k in range(5)}
    nop = {k: 0 for k in range(5)}
    nr22 = {k: 0 for k in range(5)}
    for i, (pr, p) in enumerate(cases):
        X = Z[f"{pr}|{p:g}|0"]
        for j, e in enumerate(ELEVS):
            s = np.sin(np.radians(e))
            for k in range(5):
                a1 = X[0][j, k]
                sh = {"Csc": X[3][j, k], "Ptwh": X[6][j, k], "Ptwc": X[7][j, k], "Pseudo": a1 - pseudo(pr, p, e, k)}
                if a1 <= AMAX:
                    nop[k] += 1
                    dev["Prof"][k].append(std_annex1(p)[e][k] - a1)
                    for r, v in sh.items():
                        dev[r][k].append(v - a1)
                r22 = MS[f"R22_{i}_{e}"][k]
                if a1 <= AMAX:
                    nr22[k] += 1
                    refr[k].append(MS[f"p676_{i}_{e}"][k] - r22)  # same coarse layers as R22
                    if k >= 3:
                        devr["Csc"][k].append(MS[f"R22_{i}_90"][k] / s - r22)
                        devr["Prof"][k].append(std_r22(p)[e][k - 3] - r22)
                        for r in ("Ptwh", "Ptwc", "Pseudo"):
                            devr[r][k].append(sh[r] - r22)
    # against am
    deva = {r: {k: [] for k in (3, 4)} for r in ROWS}
    refa = {k: [] for k in (3, 4)}
    nam = {k: 0 for k in (3, 4)}
    common = {k: [] for k in (3, 4)}
    ELS = [10, 30, 90]
    for key, r in SL.items():
        pr, p = key.split("|")
        p = float(p)
        am, am_std = np.array(r["am"]), np.array(r["am_std"])
        for j, e in enumerate(ELS):
            for k in (3, 4):
                a1 = r["p676"][j][k]
                if a1 > AMAX:
                    continue
                nam[k] += 1
                a = am[j, k]
                ci = cases.index((pr, p))
                r22c = MS[f"R22_{ci}_{e}"][k]
                common[k].append((a1, a, r22c))
                refa[k].append(a1 - a)
                deva["Csc"][k].append(am[2, k] / np.sin(np.radians(e)) - a)
                deva["Prof"][k].append(am_std[j, k] - a)
                deva["Ptwh"][k].append(r["a2r"][j][k] - a)
                deva["Ptwc"][k].append(r["a2vr"][j][k] - a)
                deva["Pseudo"][k].append(a1 - pseudo(pr, p, e, k) - a)

    m = {}
    for k in (3, 4):
        w = WN[k]
        m[f"BudRefRtwo{w}Db"] = f"{worst(refr[k]):+.1f}"
        m[f"BudRefAm{w}Db"] = f"{worst(refa[k]):+.1f}"
        m[f"BudRefRtwo{w}AbsDb"] = f"{abs(worst(refr[k])):.1f}"
        m[f"BudRefAm{w}AbsDb"] = f"{abs(worst(refa[k])):.1f}"
        m[f"BudN{w}"] = f"{nop[k]:d}"
        c_ = np.array(common[k])
        m[f"BudCommonAm{w}Pct"] = f"{100 * np.max((c_[:, 0] - c_[:, 1]) / c_[:, 1]):.0f}"
        m[f"BudCommonRtwo{w}Pct"] = f"{100 * np.max((c_[:, 0] - c_[:, 2]) / c_[:, 2]):.0f}"
        m[f"BudCommonAm{w}Db"] = f"{np.max(c_[:, 0] - c_[:, 1]):.1f}"
        m[f"BudCommonRtwo{w}Db"] = f"{np.max(c_[:, 0] - c_[:, 2]):.1f}"
        m[f"BudNAm{w}"] = f"{nam[k]:d}"
        for r in ROWS:
            m[f"Bud{r}{w}AnnexDb"] = f"{worst(dev[r][k]):+.1f}"
            m[f"Bud{r}{w}RtwoDb"] = f"{worst(devr[r][k]):+.1f}"
            m[f"Bud{r}{w}AmDb"] = f"{worst(deva[r][k]):+.1f}"
    # every-window check: smallest per-window worst deviation of each shortcut against Annex 1
    table = {r: [worst(dev[r][k]) for k in range(5)] for r in ROWS}
    table["Annex1_minus_R22"] = [worst(refr[k]) for k in range(5)]
    for r, nm in (("Prof", "Prof"), ("Ptwh", "Ptwh"), ("Ptwc", "Ptwc"), ("Pseudo", "Pseudo")):
        m[f"Bud{nm}MinWinDb"] = f"{min(abs(x) for x in table[r]):.1f}"
    m["BudRefMinHighDb"] = f"{min(abs(worst(q[k])) for q in (refr, refa) for k in (3, 4)):.1f}"
    m["BudCscAllRefsMaxDb"] = f"{max(abs(worst(q[r][k])) for q in (dev, devr, deva) for r in ['Csc'] for k in (3, 4)):.1f}"
    m["BudRefMinWinDb"] = f"{min(abs(x) for x in table['Annex1_minus_R22']):.1f}"
    # cosecant law at 5 degrees, inside the same 40 dB cap (the elevation floor of the claim)
    five = []
    for pr, p in cases:
        T, P, rho = A.state(hm, pwv_mm=p, profile=pr)
        n = A.refractive_index(T, P, rho)
        G = A.gamma_grid(cs, hm, T, P, rho)
        a5 = A.annex1_paths_geo(5.0, h, d, n) @ G
        z = A.annex1_paths_geo(90.0, h, d, n) @ G
        c = z / np.sin(np.radians(5.0))
        ok = a5 <= AMAX
        five += list((c - a5)[ok])
    m["BudCscFiveDegDb"] = f"{np.max(np.abs(five)):.1f}"
    # ITU-Rpy score with total instead of dry pressure at the headline point of the P.676-12 heights form
    w_ = AN.atm_worst(Z, json.load(open(os.path.join(OUT, "atm_meta.json"))), 6)
    pr, p, e, f = w_[3]
    ref = w_[2]
    Ts_, Ps_, rs_ = A.state(np.array([0.0]), pwv_mm=p, profile=pr)
    tot = float(A.itur_annex2(np.array([f]), 90.0, float(rs_[0]), float(Ps_[0]), float(Ts_[0]))[0]) / np.sin(np.radians(e))
    m["BudRtwTotalPressureDb"] = f"{abs(tot - ref):.1f}"
    m["BudRtwDryPressureDb"] = f"{abs(w_[0]):.1f}"
    with open(os.path.join(HERE, "..", "budget_macros.tex"), "w") as fh:
        fh.write("% Generated by sim/budget.py. Do not edit by hand.\n")
        for k_ in sorted(m):
            fh.write(f"\\newcommand{{\\{k_}}}{{{A.tex_minus(m[k_])}}}\n")
    json.dump({"macros": m, "annex1_all_windows": table, "n": nop, "worst_point": str(w_)},
              open(os.path.join(OUT, "budget.json"), "w"), indent=1)
    print(json.dumps({"macros": m, "annex1_all_windows": table, "n": nop}, indent=1))


if __name__ == "__main__":
    main()
