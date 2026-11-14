#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""Sounding-driven profile rule of the availability study at Svalbard (340 GHz window).

The three rules of analysis.py all use P.835-6 profiles scaled to the P.836 column. Here each
Ny-Alesund sounding of one year (1 Dec 2023 to 30 Nov 2024, data/sondes_nyalesund_year.npz) is itself
the atmosphere: full 10 GHz channel in 20 sub-bands, refracted ray on the Annex 1 layers, T_mr from
radiative transfer through that sounding, the Svalbard elevation histogram, the same terminal. The
heights form and the column form of P.676-13 Annex 2 take the surface state of the sounding (and, for
the column form, its column). Each month carries the same weight, whatever its number of soundings.
The December to February plus July subset gives the seasonal-weight sensitivity.
Writes sonde_rate_macros.tex, out/sonde_rate.npz and out/sonde_rate.json.
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
EG = np.arange(10.0, 90.01, 2.0)
SUB = AN.SUB
DATA = os.path.join(HERE, "data", "sondes_nyalesund_year.npz")


def sounding(S, i):
    """Layers and state of sounding i, spliced to P.835-6 above the burst as in analysis.stage_sonde.
    The splice takes the winter profile from October to March and the summer one otherwise; December
    to February and July therefore match stage_sonde."""
    z0 = float(S["z0"][i])
    h, d = A.annex1_layers_from(max(z0, 1e-3))
    hm = h + 0.5 * d
    zz = z0 + S["zg"]
    ok = np.isfinite(S["T"][i])
    top = zz[ok][-1]
    mon = int(str(S["date"][i])[5:7])
    season = "high-winter" if mon in (10, 11, 12, 1, 2, 3) else "high-summer"
    Tp, Pp, rp = A.p835_seasonal(hm, season)
    _, Ptop, _ = A.p835_seasonal(np.array([top]), season)
    below = hm <= top
    T = np.where(below, np.interp(hm, zz[ok], S["T"][i][ok]), Tp)
    P = np.where(below, np.exp(np.interp(hm, zz[ok], np.log(S["P"][i][ok]))), Pp * S["P"][i][ok][-1] / float(Ptop[0]))
    rho = np.where(below, np.interp(hm, zz[ok], S["rho"][i][ok]), np.maximum(rp, 2e-6 * P * 216.7 / T))
    return h, d, hm, T, P, rho, float(S["T"][i][0]), float(S["P"][i][0]), float(S["rho"][i][0])


def fill(S, fc):
    n = len(S["date"])
    shp = (n, len(EG), len(SUB))
    out = {w: np.zeros(shp) for w in ("Aref", "Tmr", "Aa2", "Aa2v", "Aa2r12", "Aa2s", "Aa2l")}
    V_ = np.zeros(n)
    he = np.zeros(n)
    inv = 1.0 / np.sin(np.radians(EG))[:, None]
    for i in range(n):
        h, d, hm, T, P, rho, Ts, Ps, rs = sounding(S, i)
        f = fc + SUB
        G = A.gamma_grid(f, hm, T, P, rho)
        nr = A.refractive_index(T, P, rho)
        paths = np.array([A.annex1_paths_geo(e, h, d, nr) for e in EG])
        out["Aref"][i] = paths @ G
        out["Tmr"][i] = A.sky_temperature_rt(paths, G, T, f)[1]
        V = A.column_mm(h, d, rho)
        out["Aa2"][i] = A.p676_13_annex2(f, 90.0, rs, Ps, Ts)[None, :] * inv
        out["Aa2v"][i] = A.p676_13_annex2(f, 90.0, rs, Ps, Ts, V_t=V)[None, :] * inv
        # the heights form with rho_s from the second record and from the lowest 100 m; reference fixed
        out["Aa2s"][i] = A.p676_13_annex2(f, 90.0, float(S["rs2"][i]), Ps, Ts)[None, :] * inv
        out["Aa2l"][i] = A.p676_13_annex2(f, 90.0, float(S["rs100"][i]), Ps, Ts)[None, :] * inv
        out["Aa2r12"][i] = A.itur_annex2(f, 90.0, rs, Ps - rs * Ts / 216.7, Ts)[None, :] * inv
        V_[i], he[i] = V, V / rs
        if i % 25 == 0:
            print("sonde_rate", S["date"][i], round(V, 2), flush=True)
    out.update(V=V_, heff=he, date=np.array([str(x) for x in S["date"]]))
    return out


def elev_weights(res):
    eh = np.histogram(res["Svalbard_el"], bins=np.arange(10, 91, 1.0))
    return 0.5 * (eh[1][1:] + eh[1][:-1]), eh[0] / eh[0].sum()


def rates(Att, Tmr, ec, fc):
    """Rate (Gbit/s) per (sounding, elevation bin): the same link as analysis.rate_samples, sky='rt'."""
    re_, hs = 6371.0, AN.TERM["h_sat"]
    dkm = np.sqrt((re_ + hs) ** 2 - (re_ * np.cos(np.radians(ec))) ** 2) - re_ * np.sin(np.radians(ec))
    df = (SUB[1] - SUB[0]) * 1e9
    T_rx = 290.0 * (10 ** (AN.TERM["nf_db"] / 10) - 1)
    rate = np.zeros((Att.shape[0], len(ec)))
    for q, off in enumerate(SUB):
        att = np.stack([np.interp(ec, EG, Att[i, :, q]) for i in range(Att.shape[0])])
        tm = np.stack([np.interp(ec, EG, Tmr[i, :, q]) for i in range(Att.shape[0])])
        t = 10 ** (-att / 10)
        T_sky = tm * (1 - t) + A.planck_temperature(2.725, fc + off) * t
        noise = 10 * np.log10(1.380649e-23 * (T_rx + T_sky) * df) + 30.0
        f = fc + off
        lam = 0.299792458e-3 / f
        fspl = 20 * np.log10(4 * np.pi * dkm / lam)
        snr = (AN.TERM["pt_dbm"] - 10 * np.log10(len(SUB)) + AN.gain_db(AN.TERM["d_sat"], f)
               + AN.gain_db(AN.TERM["d_gnd"], f) - fspl[None, :] - att - noise)
        rate += df * np.log2(1 + 10 ** (snr / 10)) / 1e9
    return rate


def calibrate(res, fc, ec, ew):
    """rates() must reproduce analysis.rate_samples on the P.835 grid at one tabulated column."""
    i = 8
    r2 = dict(res)
    r2["Svalbard_V"] = np.array([res["Svalbard_pg"][i]])
    r2["Svalbard_p"] = np.array([50.0])
    r0, _ = AN.rate_samples(r2, "Svalbard", "Aref", 0, fc)
    att = res["Svalbard_Aref"][i:i + 1, :, 0, :]
    tmr = res["Svalbard_Tmr"][i:i + 1, :, 0, :]
    r1 = rates(att, tmr, ec, fc).ravel()
    err = float(np.max(np.abs(r1 - r0) / np.maximum(r0, 1e-6)))
    assert err < 1e-3, err
    return err


def stat(rate, wi, ew):
    W = wi[:, None] * ew[None, :]
    med = AN.quantile_rate(rate.ravel(), W.ravel(), 0.5)
    av = 100 * float(W[rate >= AN.RATE_TARGET].sum())
    return med, av


def main():
    res = dict(np.load(os.path.join(OUT, "avail.npz")))
    fc = float(json.load(open(os.path.join(OUT, "windows.json")))["centers"][0])
    ec, ew = elev_weights(res)
    cal = calibrate(res, fc, ec, ew)
    S = np.load(DATA)
    F = fill(S, fc)
    np.savez_compressed(os.path.join(OUT, "sonde_rate.npz"), **F)
    Rr = rates(F["Aref"], F["Tmr"], ec, fc)
    Rh = rates(F["Aa2"], F["Tmr"], ec, fc)
    Rc = rates(F["Aa2v"], F["Tmr"], ec, fc)
    Rp = rates(F["Aa2r12"], F["Tmr"], ec, fc)
    mon = np.array([int(d[5:7]) for d in F["date"]])
    nm = {m: int((mon == m).sum()) for m in range(1, 13)}
    assert all(v > 0 for v in nm.values()), nm
    wy = np.array([1 / (12 * nm[m]) for m in mon])

    def group(ms):
        w = np.where(np.isin(mon, ms), 1.0, 0.0) * wy
        return w / w.sum()

    groups = {"Year": wy, "Win": group([12, 1, 2]), "Spr": group([3, 4, 5]), "Sum": group([6, 7, 8]),
              "Aut": group([9, 10, 11]), "Jul": group([7])}
    # the seasonal-weight sensitivity on the subset the first estimate used
    for lab, a in (("Half", 0.5), ("Quarter", 0.75)):
        groups[f"Sub{lab}"] = a * groups["Win"] + (1 - a) * groups["Jul"]
    m = {"SrN": f"{len(mon):d}", "SrNMonMin": f"{min(nm.values()):d}", "SrNMonMax": f"{max(nm.values()):d}",
         "SrCalErr": AN.texsci(cal)}
    tab = {}
    for g, w in groups.items():
        mr, ar = stat(Rr, w, ew)
        mh, ah = stat(Rh, w, ew)
        mc, ac = stat(Rc, w, ew)
        tab[g] = dict(ref=mr, heights=mh, column=mc, av_ref=ar, av_heights=ah, av_column=ac)
        m[f"Sr{g}RefMed"] = f"{mr:.1f}"
        m[f"Sr{g}HeightsMed"] = f"{mh:.1f}"
        m[f"Sr{g}ColumnMed"] = f"{mc:.1f}"
        m[f"Sr{g}HeightsBiasPct"] = f"{100 * (mh - mr) / mr:+.0f}"
        m[f"Sr{g}ColumnBiasPct"] = f"{100 * (mc - mr) / mr:+.0f}"
        m[f"Sr{g}RefAv"] = f"{ar:.0f}"
        m[f"Sr{g}HeightsAv"] = f"{ah:.0f}"
        m[f"Sr{g}ColumnAv"] = f"{ac:.0f}"
    m["SrYearPtwAv"] = f"{stat(Rp, wy, ew)[1]:.0f}"
    # sensitivity to the surface record: year median and mean bias, winter median bias, summer availability
    wmean = lambda R, w: float((w[:, None] * ew[None, :] * R).sum())
    Rs_, Rl_ = rates(F["Aa2s"], F["Tmr"], ec, fc), rates(F["Aa2l"], F["Tmr"], ec, fc)
    for tag_, R_ in (("First", Rh), ("Second", Rs_), ("Low", Rl_)):
        m[f"SrRule{tag_}YearMed"] = f"{100 * (stat(R_, wy, ew)[0] - stat(Rr, wy, ew)[0]) / stat(Rr, wy, ew)[0]:+.0f}"
        m[f"SrRule{tag_}YearMean"] = f"{100 * (wmean(R_, wy) - wmean(Rr, wy)) / wmean(Rr, wy):+.1f}"
        m[f"SrRule{tag_}WinMed"] = f"{100 * (stat(R_, groups['Win'], ew)[0] - stat(Rr, groups['Win'], ew)[0]) / stat(Rr, groups['Win'], ew)[0]:+.0f}"
        m[f"SrRule{tag_}SumAv"] = f"{stat(R_, groups['Sum'], ew)[1]:.0f}"
        m[f"SrRule{tag_}YearAv"] = f"{stat(R_, wy, ew)[1]:.0f}"
    _av = [stat(R_, wy, ew)[1] for R_ in (Rh, Rs_, Rl_)]
    m["SrYearHtAvLo"], m["SrYearHtAvHi"] = f"{min(_av):.0f}", f"{max(_av):.0f}"
    # 95 percent bootstrap over soundings, stratified by month
    rg = np.random.default_rng(20261103)
    idx = {k: np.flatnonzero(mon == k) for k in range(1, 13)}
    bh, bc, bd = [], [], []
    bw, bs, bl = [], [], []
    for _ in range(2000):
        ii = np.concatenate([rg.choice(idx[k], len(idx[k])) for k in range(1, 13)])
        mr, ar = stat(Rr[ii], wy[ii], ew)
        mh, ah = stat(Rh[ii], wy[ii], ew)
        mc, _ = stat(Rc[ii], wy[ii], ew)
        bh.append(100 * (mh - mr) / mr)
        bc.append(100 * (mc - mr) / mr)
        bd.append(ah - ar)
        mo = mon[ii]
        for ms_, lst, kind in (([12, 1, 2], bw, "med"), ([6, 7, 8], bs, "av")):
            sel = np.isin(mo, ms_)
            wi = wy[ii][sel] / wy[ii][sel].sum()
            a_, b_ = stat(Rr[ii][sel], wi, ew), stat(Rh[ii][sel], wi, ew)
            lst.append(100 * (b_[0] - a_[0]) / a_[0] if kind == "med" else b_[1] - a_[1])
        sel = np.isin(mo, [12, 1, 2])
        wi = wy[ii][sel] / wy[ii][sel].sum()
        a_ = stat(Rr[ii][sel], wi, ew)[0]
        bl.append(100 * (stat(Rl_[ii][sel], wi, ew)[0] - a_) / a_)
    for nm_, arr, fmt in (("Heights", bh, "{:+.0f}"), ("Column", bc, "{:+.0f}"), ("AvDiff", bd, "{:+.0f}"), ("WinBias", bw, "{:+.0f}"), ("SumAvDiff", bs, "{:+.0f}"), ("LowWinBias", bl, "{:+.0f}")):
        m[f"SrYear{nm_}Lo"] = fmt.format(np.percentile(arr, 2.5))
        m[f"SrYear{nm_}Hi"] = fmt.format(np.percentile(arr, 97.5))
    # the column statistics of the soundings against P.836
    W1 = wy / wy.sum()
    o = np.argsort(F["V"])
    m["SrVMed"] = f"{F['V'][o][np.searchsorted(np.cumsum(W1[o]), 0.5)]:.1f}"
    o = np.argsort(F["heff"])
    m["SrHeffMed"] = f"{F['heff'][o][np.searchsorted(np.cumsum(W1[o]), 0.5)]:.1f}"
    with open(os.path.join(HERE, "..", "sonde_rate_macros.tex"), "w") as f:
        f.write("% Generated by sim/sonde_rate.py. Do not edit by hand.\n")
        for k in sorted(m):
            f.write(f"\\newcommand{{\\{k}}}{{{A.tex_minus(m[k])}}}\n")
    json.dump({"macros": m, "groups": tab, "n_month": nm}, open(os.path.join(OUT, "sonde_rate.json"), "w"), indent=1)
    print(json.dumps({"macros": m, "groups": tab, "n_month": nm}, indent=1))


if __name__ == "__main__":
    main()
