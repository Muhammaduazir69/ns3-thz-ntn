#!/usr/bin/env python3
"""Absorbing-height analysis of the Annex 2 heights form (Section III-F of the T-TST manuscript).

Proposition 1: the zenith deviation of the heights form is
    (gamma_o h_o - A_o + gamma_w h_w - A_w) / A,  with  A_w = gamma_w,s H_gamma,
so on the water-vapor part it is (h_w / H_gamma - 1) A_w / A, where H_gamma = A_w / gamma_w,s is the
absorbing height of the profile.
Corollary 1: in a window, H_gamma ~ (1 - sigma) / (1/h_eff + 1/H_PT) + sigma h_eff / 2, with
h_eff = V / rho_s, sigma the self-broadened share of the surface water-vapor absorption and H_PT a
pressure-temperature scale height. H_PT is fitted once on the six P.835-6 atmospheres and then
tested, unchanged, on the Ny-Alesund radiosondes (held out). The heights form is unbiased at
h_eff* = 1 / (1/h_w - 1/H_PT) for sigma -> 0.
Writes out/heff_model.json and heff_macros.tex (run with sim/.venv/bin/python).
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thzatm as A  # noqa: E402

OUT = os.path.join(HERE, "out")


def h_w(fc):
    return 5.6585e-5 * fc + 1.8348 + sum(ai / ((fc - fi) ** 2 + bi) for fi, ai, bi in
                                         ((22.235080, 2.6846, 2.7649), (183.310087, 5.8905, 4.9219),
                                          (325.152888, 2.9810, 3.0748)))


def profile_terms(h, d, T, P, rho, fc):
    """Exact Annex 1 wet and dry zenith losses, surface gamma_w, self share and h_eff."""
    g = [A.gamma_p676(fc, T[k], P[k], rho[k], split=True) for k in range(len(h))]
    go = np.array([x[0][0] for x in g])
    gw = np.array([x[1][0] for x in g])
    ts, ps, rs = float(T[0]), float(P[0]), float(rho[0])
    es = rs * ts / 216.7
    gself = A.gamma_p676(fc, ts, es + 1e-6, rs, split=True)[1][0]
    V = A.column_mm(h, d, rho)
    return dict(Aw=float((gw * d).sum()), Ao=float((go * d).sum()), gws=float(gw[0]),
                sig=float(gself / gw[0]), heff=float(V / rs), V=float(V))


def closed_H(heff, sig, H):
    return (1 - sig) / (1 / heff + 1 / H) + sig * heff / 2


def p835_terms(fc):
    h, d = A.annex1_layers_from(1e-3)
    hm = h + 0.5 * d
    rows = {}
    for kind in ["low", "mid-summer", "mid-winter", "high-summer", "high-winter"]:
        T, P, rho = A.p835_seasonal(hm, kind)
        rows[kind] = profile_terms(h, d, T, P, rho, fc)
    T, P = A.p835_TP(hm)
    rows["p835-std"] = profile_terms(h, d, T, P, A.p835_rho(hm), fc)
    return rows


def sonde_terms(fc):
    S = np.load(os.path.join(HERE, "data", "sondes_nyalesund_year.npz"))
    rows = json.load(open(os.path.join(OUT, "sonde.json")))
    out = []
    for i, date in enumerate(S["date"]):
        z0 = float(S["z0"][i])
        h, d = A.annex1_layers_from(max(z0, 1e-3))
        hm = h + 0.5 * d
        zz = z0 + S["zg"]
        ok = np.isfinite(S["T"][i])
        top = zz[ok][-1]
        season = "high-winter" if str(date)[5:7] in ("10", "11", "12", "01", "02", "03") else "high-summer"
        Tp, Pp, rp = A.p835_seasonal(hm, season)
        _, Ptop, _ = A.p835_seasonal(np.array([top]), season)
        below = hm <= top
        T = np.where(below, np.interp(hm, zz[ok], S["T"][i][ok]), Tp)
        P = np.where(below, np.exp(np.interp(hm, zz[ok], np.log(S["P"][i][ok]))),
                     Pp * S["P"][i][ok][-1] / float(Ptop[0]))
        rho = np.where(below, np.interp(hm, zz[ok], S["rho"][i][ok]), np.maximum(rp, 2e-6 * P * 216.7 / T))
        t = profile_terms(h, d, T, P, rho, fc)
        r = rows[i]
        t["zrel"] = (r["a2"][-1][0] - r["ref"][-1][0]) / r["ref"][-1][0]
        t["a2z"] = r["a2"][-1][0]
        t["season"] = season
        t["group"] = r["group"]
        t["date"] = r["date"]
        out.append(t)
    return out


def main():
    fc = float(json.load(open(os.path.join(OUT, "windows.json")))["centers"][0])
    hw = h_w(fc)
    ref = p835_terms(fc)
    loss = lambda H: sum((closed_H(r["heff"], r["sig"], H) / (r["Aw"] / r["gws"]) - 1) ** 2 for r in ref.values())
    H = float(minimize_scalar(loss, bounds=(2, 200), method="bounded").x)
    ref_err = [abs(closed_H(r["heff"], r["sig"], H) / (r["Aw"] / r["gws"]) - 1) for r in ref.values()]
    so = sonde_terms(fc)
    Hg = np.array([r["Aw"] / r["gws"] for r in so])
    Hc = np.array([closed_H(r["heff"], r["sig"], H) for r in so])
    wet = np.array([r["Aw"] / (r["Aw"] + r["Ao"]) for r in so])
    zrel = np.array([r["zrel"] for r in so])
    ident = (hw / Hg - 1) * wet
    # fully closed prediction: the water share also comes from the closed form, not from Annex 1
    a2z = np.array([r["a2z"] for r in so])
    gws = np.array([r["gws"] for r in so])
    Ac = a2z - gws * hw + gws * Hc
    pred = (a2z - Ac) / Ac
    he = np.array([r["heff"] for r in so])
    res = dict(fc=fc, hw=hw, H_PT=H, ref={k: dict(v, Hg=v["Aw"] / v["gws"]) for k, v in ref.items()},
               sondes=dict(heff=he.tolist(), zrel=zrel.tolist(), pred=pred.tolist(), ident=ident.tolist(),
                           Hg=Hg.tolist(), Hc=Hc.tolist(), sig=[r["sig"] for r in so], wet=wet.tolist()))
    json.dump(res, open(os.path.join(OUT, "heff_model.json"), "w"), indent=1)
    sig_med = float(np.median([r["sig"] for r in so]))
    mac = {
        "HmHw": f"{hw:.2f}",
        "HmHpt": f"{H:.1f}",
        "HmRefMaxErrPct": f"{100 * max(ref_err):.0f}",
        "HmIdentRmsPts": f"{100 * np.sqrt(np.mean((ident - zrel) ** 2)):.2f}",
        "HmHgMeanErrPct": f"{100 * np.mean(np.abs(Hc / Hg - 1)):.1f}",
        "HmHgCorr": f"{np.corrcoef(Hg, Hc)[0, 1]:.3f}",
        "HmDevRmsPts": f"{100 * np.sqrt(np.mean((pred - zrel) ** 2)):.1f}",
        "HmDevCorr": f"{np.corrcoef(pred, zrel)[0, 1]:.3f}",
        "HmSlopePred": f"{np.polyfit(he, 100 * pred, 1)[0]:.0f}",
        "HmHeffStar": f"{1 / (1 / hw - 1 / H):.1f}",
        "HmSigMinPct": f"{100 * min(r['sig'] for r in so):.0f}",
        "HmSigMaxPct": f"{100 * max(r['sig'] for r in so):.0f}",
        "HmSigMedPct": f"{100 * sig_med:.0f}",
        "HmWetMinPct": f"{100 * wet.min():.0f}",
        "HmFc": f"{fc:.1f}",
        "HmHighWinPredPct": f"{100 * (hw / closed_H(ref['high-winter']['heff'], ref['high-winter']['sig'], H) - 1):.0f}",
        "HmHighWinHeff": f"{ref['high-winter']['heff']:.1f}",
    }
    # the P.835-6 high-latitude winter profile against the winter soundings it stands for
    win = np.array([r["group"] == "win" for r in so])
    jja = np.array([r["group"] == "sum" for r in so])
    # out-of-sample test: the closed form was written on 29 Oct 2026, after the 109 soundings of
    # data/sondes_nyalesund.npz had been seen and before the other days of the year were fetched
    seen = {str(x) for x in np.load(os.path.join(HERE, "data", "sondes_nyalesund.npz"))["date"]}
    old = np.array([r["date"] in seen for r in so])
    months = np.array([r["date"][5:7] for r in so])
    mrms = [np.sqrt(np.mean((pred[months == m] - zrel[months == m]) ** 2)) for m in sorted(set(months))]
    mac.update({
        "HmSeenN": f"{int(old.sum()):d}",
        "HmSeenRmsPts": f"{100 * np.sqrt(np.mean((pred[old] - zrel[old]) ** 2)):.1f}",
        "HmOosN": f"{int((~old).sum()):d}",
        "HmOosRmsPts": f"{100 * np.sqrt(np.mean((pred[~old] - zrel[~old]) ** 2)):.1f}",
        "HmOosCorr": f"{np.corrcoef(pred[~old], zrel[~old])[0, 1]:.3f}",
        "HmOosBiasPts": f"{100 * np.mean(pred[~old] - zrel[~old]):+.1f}",
        "HmMonRmsMinPts": f"{100 * min(mrms):.1f}",
        "HmMonRmsMaxPts": f"{100 * max(mrms):.1f}",
        "HmHeffMin": f"{he.min():.2f}",
        "HmHeffMax": f"{he.max():.1f}",
    })
    zw, hw_ = 100 * zrel[win], he[win]
    mac.update({
        "HmWinZenLoPct": f"{np.percentile(zw, 25):+.0f}",
        "HmWinZenHiPct": f"{np.percentile(zw, 75):+.0f}",
        "HmWinNegPct": f"{100 * np.mean(zw < 0):.0f}",
        "HmWinNHeffAboveProf": f"{int(np.sum(hw_ >= ref['high-winter']['heff'] - 0.005)):d}",
        "HmWinNBelowTwentyFive": f"{int(np.sum(zw < -25)):d}",
        "HmZenMinPct": f"{100 * zrel.min():+.0f}",
        "HmZenMaxPct": f"{100 * zrel.max():+.0f}",
    })
    # a straight line is a weak test of a 1/h_eff law: fit quality and slope per season
    lin = np.polyfit(he, 100 * zrel, 1)
    inv = np.polyfit(1 / he, 100 * zrel, 1)
    mac.update({
        "HmLinRmsPts": f"{np.sqrt(np.mean((np.polyval(lin, he) - 100 * zrel) ** 2)):.1f}",
        "HmInvRmsPts": f"{np.sqrt(np.mean((np.polyval(inv, 1 / he) - 100 * zrel) ** 2)):.1f}",
        "HmSlopeWin": f"{np.polyfit(he[win], 100 * zrel[win], 1)[0]:.0f}",
        "HmSlopeSum": f"{np.polyfit(he[jja], 100 * zrel[jja], 1)[0]:.0f}",
        "HmSlopePredWin": f"{np.polyfit(he[win], 100 * pred[win], 1)[0]:.0f}",
        "HmSlopePredSum": f"{np.polyfit(he[jja], 100 * pred[jja], 1)[0]:.0f}",
    })
    # first-order slope of Corollary 1 against the exact derivative of the closed form, at h_eff*
    hs, sg = 1 / (1 / hw - 1 / H), sig_med
    fo = -(hw / hs ** 2) * (1 + sg * (H ** 2 + hs ** 2) / (2 * H ** 2))
    dh = 1e-5
    ex = (hw / closed_H(hs + dh, sg, H) - hw / closed_H(hs - dh, sg, H)) / (2 * dh)
    assert abs(fo - ex) / abs(ex) < 0.01, (fo, ex)
    mac["HmCorSlopeExact"] = f"{ex:.2f}"
    mac["HmCorSlopeFirst"] = f"{fo:.2f}"
    with open(os.path.join(HERE, "..", "heff_macros.tex"), "w") as f:
        f.write("% written by sim/heff_model.py; do not edit\n")
        for k in sorted(mac):
            f.write(f"\\newcommand{{\\{k}}}{{{A.tex_minus(mac[k])}}}\n")
    print(json.dumps(mac, indent=1))


if __name__ == "__main__":
    main()
