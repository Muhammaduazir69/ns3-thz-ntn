#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""Reach of the heights-form analysis beyond Ny-Alesund: six other Arctic stations.

Same processing as analysis.stage_sonde and heff_model.py, H_PT unchanged (fitted on the six P.835-6
atmospheres). Per station: h_eff, the Annex 2 zenith deviation at 341.5 GHz in December to February,
the rms of Corollary 1 over the year's sample, and the sounding-driven rate bias of sonde_rate.py with
that site's own elevation statistics. Sensitivity of rho_s: first record, second record, mean over the lowest 100 m, with h_eff, the self share, gamma_w,s and the limit recomputed per rule. Season scope:
median deviation, loss and bootstrap interval for December to February and June to August. Pooled with the
Ny-Alesund year: rms per h_eff bin, a validity limit from a stated rule, the fit-free rival
(heights form / column form - 1), inversion strength against ln h_eff, and a test of assumption (A1), the ratio Q = (integral of rho^2 dh / rho_s^2) / (h_eff / 2), which is 1 for an exponential.
Usage: arctic.py --prep <dir> writes data/sondes_arctic.npz from the raw files; arctic.py analyzes it.
Writes arctic_macros.tex, out/arctic.json.
"""
import csv
import glob
import json
import os
import sys
import tarfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analysis as AN  # noqa: E402
import heff_model as HM  # noqa: E402
import prep_sondes as PS  # noqa: E402
import sonde_rate as SR  # noqa: E402
import thzatm as A  # noqa: E402

OUT = os.path.join(HERE, "out")
DATA = os.path.join(HERE, "data", "sondes_arctic.npz")
# id, tag, name, lat, lon
ST = [("01028", "Bj", "Bj{\\o}rn{\\o}ya", 74.50, 19.00), ("04320", "Da", "Danmarkshavn", 76.77, -18.67),
      ("70026", "Ut", "Utqiagvik", 71.28, -156.79), ("71082", "Al", "Alert", 82.50, -62.33),
      ("71917", "Eu", "Eureka", 79.98, -85.93), ("02836", "So", "Sodankyl{\\\"a}", 67.37, 26.63)]


def prep(src):
    out = {"zg": PS.ZG}
    for sid, tag, *_ in ST:
        dates, Z0, Ps, Ts, Rs, R2, R100 = [], [], [], [], [], [], []
        for pth in sorted(glob.glob(os.path.join(src, f"s{sid}_*_12.csv"))):
            s = PS.read(pth)
            if s is None:
                continue
            z, P, T, rho = s
            zz = z[0] + PS.ZG
            m = zz <= z[-1]
            Ps.append(np.where(m, np.exp(np.interp(zz, z, np.log(P))), np.nan))
            Ts.append(np.where(m, np.interp(zz, z, T), np.nan))
            Rs.append(np.where(m, np.interp(zz, z, rho), np.nan))
            Z0.append(z[0])
            a_, b_ = PS.surface_rules(z, rho)
            R2.append(a_)
            R100.append(b_)
            dates.append(os.path.basename(pth)[7:17])
        out.update({f"{tag}_date": np.array(dates), f"{tag}_z0": np.array(Z0), f"{tag}_P": np.array(Ps),
                    f"{tag}_T": np.array(Ts), f"{tag}_rho": np.array(Rs), f"{tag}_rs2": np.array(R2), f"{tag}_rs100": np.array(R100)})
        print(sid, len(dates))
    np.savez_compressed(DATA, **out)


def terms(S, i, fc, hw, H):
    h, d, hm, T, P, rho, Ts, Ps_, rs = SR.sounding(S, i)
    t = HM.profile_terms(h, d, T, P, rho, fc)
    G = A.gamma_grid(np.array([fc]), hm, T, P, rho)
    n = A.refractive_index(T, P, rho)
    ref = float((A.annex1_paths_geo(90.0, h, d, n) @ G)[0])
    a2f = lambda r_: float(A.p676_13_annex2(np.array([fc]), 90.0, r_, Ps_, Ts)[0])
    a2 = a2f(rs)
    a2v = float(A.p676_13_annex2(np.array([fc]), 90.0, rs, Ps_, Ts, V_t=t["V"])[0])
    Hc = HM.closed_H(t["heff"], t["sig"], H)
    Ac = a2 - t["gws"] * hw + t["gws"] * Hc  # fully closed: no Annex 1 share
    Q = float((rho ** 2 * d).sum() / rho[0] ** 2 / (t["heff"] / 2))
    zg = np.asarray(S["zg"])
    Tp = np.asarray(S["T"][i])
    inv = float(np.nanmax(Tp[zg <= 1.0]) - Tp[0])
    return dict(heff=t["heff"], zrel=(a2 - ref) / ref, zrelv=(a2v - ref) / ref, pred=(a2 - Ac) / Ac, Q=Q, V=t["V"],
                rival=a2 / a2v - 1, inv=inv, ref=ref, a2=a2, gws=t["gws"], sig=t["sig"], r1=float(rho[0]), **rule_terms(S, i, fc, hw, H, a2f, ref, t["V"], Ts, Ps_))


def rule_terms(S, i, fc, hw, H, a2f, ref, V, Ts, Ps_):
    """Deviation, closed form and rival with rho_s from the second record and from the lowest 100 m."""
    out = {}
    for tag, key in (("2", "rs2"), ("l", "rs100")):
        rr = max(float(S[key][i]), 1e-6)
        a2 = a2f(rr)
        a2v = float(A.p676_13_annex2(np.array([fc]), 90.0, rr, Ps_, Ts, V_t=V)[0])
        gw = float(A.gamma_p676(fc, Ts, Ps_, rr, split=True)[1][0])
        gself = float(A.gamma_p676(fc, Ts, rr * Ts / 216.7 + 1e-6, rr, split=True)[1][0])
        Ac = a2 - gw * hw + gw * HM.closed_H(V / rr, gself / gw, H)
        out.update({f"z{tag}": (a2 - ref) / ref, f"he{tag}": V / rr, f"p{tag}": (a2 - Ac) / Ac, f"rv{tag}": a2 / a2v - 1,
                    f"r{tag}": rr})
    return out


def raw_levels():
    """Levels with humidity in the lowest 3 km of each raw winter or summer file, read from the raw tar."""
    lev = {}
    with tarfile.open(os.path.join(HERE, "data", "raw_soundings.tar.xz")) as tf:
        for mem in tf:
            base = os.path.basename(mem.name)
            if not (mem.isfile() and base.endswith("_12.csv") and mem.name.startswith("arctic_raw/s")):
                continue
            rows = list(csv.DictReader(line.decode("utf-8", "replace") for line in tf.extractfile(mem)))
            zz = []
            for r_ in rows:
                try:
                    if r_["mixing ratio_g/kg"] not in ("", None):
                        zz.append(float(r_["geopotential height_m"]))
                except (KeyError, TypeError, ValueError):
                    continue
            lev[(base[1:6], base[7:17])] = int(np.sum(np.array(zz) < min(zz) + 3000.0)) if zz else 0
    return lev


def view(D, tag):
    """A sonde_rate.sounding-compatible view of one station."""
    return {"z0": D[f"{tag}_z0"], "zg": D["zg"], "T": D[f"{tag}_T"], "P": D[f"{tag}_P"], "rho": D[f"{tag}_rho"], "date": D[f"{tag}_date"],
            "rs2": D[f"{tag}_rs2"], "rs100": D[f"{tag}_rs100"]}


def station(S, fc, hw, H, lat=None, lon=None):
    R = [terms(S, i, fc, hw, H) for i in range(len(S["date"]))]
    a = {k: np.array([r[k] for r in R]) for k in R[0]}
    a["date"] = np.array([str(x) for x in S["date"]])
    mon = np.array([int(str(x)[5:7]) for x in S["date"]])
    w = np.isin(mon, [12, 1, 2])
    t = dict(n=len(R), nwin=int(w.sum()), heffwin=float(np.median(a["heff"][w])), zwin=float(np.median(a["zrel"][w])),
             zwin2=float(np.median(a["z2"][w])), zwinl=float(np.median(a["zl"][w])), neg=float(np.mean(a["zrel"][w] < 0)),
             neg2=float(np.mean(a["z2"][w] < 0)), negl=float(np.mean(a["zl"][w] < 0)),
             step=float(np.mean(np.abs(S["rho"][:, 0] / S["rs2"] - 1)[w] > 0.1)),
             rms=float(np.sqrt(np.mean((a["pred"] - a["zrel"]) ** 2))), rival=float(np.sqrt(np.mean((a["rival"] - a["zrel"]) ** 2))),
             hemin=float(a["heff"].min()), hemax=float(a["heff"].max()), zv=[float(a["zrelv"].min()), float(a["zrelv"].max())],
             z=[float(a["zrel"].min()), float(a["zrel"].max())], invcorr=float(np.corrcoef(a["inv"], np.log(a["heff"]))[0, 1]))
    if lat is not None:
        el = AN.elevation_samples(lat, lon)
        eh = np.histogram(el, bins=np.arange(10, 91, 1.0))
        ec, ew = 0.5 * (eh[1][1:] + eh[1][:-1]), eh[0] / eh[0].sum()
        F = SR.fill(S, fc)
        Rr = SR.rates(F["Aref"], F["Tmr"], ec, fc)
        nm = {k: max(int((mon == k).sum()), 1) for k in range(1, 13)}
        wy = np.array([1 / (len(set(mon)) * nm[k]) for k in mon])
        t["mr"], t["ar"] = SR.stat(Rr, wy, ew)
        for key, lab in (("Aa2", "h"), ("Aa2s", "s"), ("Aa2l", "l")):
            t[f"m{lab}"], t[f"a{lab}"] = SR.stat(SR.rates(F[key], F["Tmr"], ec, fc), wy, ew)
    return t, a


def main():
    fc = float(json.load(open(os.path.join(OUT, "windows.json")))["centers"][0])
    hmj = json.load(open(os.path.join(OUT, "heff_model.json")))
    hw, H = hmj["hw"], hmj["H_PT"]
    D = dict(np.load(DATA))
    keys = ("heff", "zrel", "pred", "Q", "rival", "inv", "ref", "a2", "gws", "sig", "r1", "r2", "rl", "z2", "zl", "p2", "pl", "he2", "hel",
            "rv2", "rvl", "date")
    table, pts = {}, {k: [] for k in keys}
    Y = dict(np.load(os.path.join(HERE, "data", "sondes_nyalesund_year.npz")))
    t, a = station(Y, fc, hw, H)
    t["name"] = "Ny-{\\AA}lesund"
    table["01004"] = t
    ny = {k: a[k].tolist() for k in pts}
    lev = raw_levels()
    for sid, tag, name, lat, lon in ST:
        t, a = station(view(D, tag), fc, hw, H, lat, lon)
        t["name"] = name
        t["lev"] = [lev[(sid, x)] for x in a["date"].tolist()]
        table[sid] = t
        for k in pts:
            pts[k] += a[k].tolist()
        print(sid, json.dumps(t), flush=True)
    json.dump({"macros": {}, "table": table, "points": pts, "ny": ny}, open(os.path.join(OUT, "arctic.json"), "w"), indent=1)
    summarize()


def summarize():
    d = json.load(open(os.path.join(OUT, "arctic.json")))
    table, pts, ny = d["table"], d["points"], d["ny"]
    m = {}
    order = [s[0] for s in ST]
    sid_all = np.array(sum([[s] * table[s]["n"] for s in order], []) + ["01004"] * table["01004"]["n"])
    All = {k: np.array(pts[k] + ny[k]) for k in pts}
    mon_all = np.array([int(x[5:7]) for x in All["date"]])
    rg = np.random.default_rng(20261103)
    med_ci = lambda x: np.percentile([np.median(rg.choice(x, len(x))) for _ in range(2000)], [2.5, 97.5])
    season = {}
    for sid, tag, name in [("01004", "Ny", "Ny-{\\AA}lesund")] + [(s[0], s[1], s[2]) for s in ST]:
        t = table[sid]
        m[f"Arc{tag}N"] = f"{t['n']:d}"
        m[f"Arc{tag}NWin"] = f"{t['nwin']:d}"
        m[f"Arc{tag}Name"] = name
        m[f"Arc{tag}Heff"] = f"{t['heffwin']:.1f}"
        m[f"Arc{tag}Zwin"] = f"{100 * t['zwin']:+.0f}"
        m[f"Arc{tag}ZwinTwo"] = f"{100 * t['zwin2']:+.0f}"
        m[f"Arc{tag}ZwinLow"] = f"{100 * t['zwinl']:+.0f}"
        m[f"Arc{tag}Rms"] = f"{100 * t['rms']:.1f}"
        m[f"Arc{tag}Rival"] = f"{100 * t['rival']:.1f}"
        m[f"Arc{tag}Step"] = f"{100 * t['step']:.0f}"
        m[f"Arc{tag}InvCorr"] = f"{t['invcorr']:.2f}"
        q = sid_all == sid
        for lab, mm in (("Win", (12, 1, 2)), ("Sum", (6, 7, 8))):
            s_ = q & np.isin(mon_all, mm)
            dA = All["a2"][s_] - All["ref"][s_]
            lo, hi = med_ci(100 * All["zrel"][s_])
            season[(tag, lab)] = dict(z=[float(np.median(All[k][s_])) for k in ("zrel", "z2", "zl")], db=float(np.median(dA)),
                                      ref=float(np.median(All["ref"][s_])), n=int(s_.sum()), lo=float(lo), hi=float(hi),
                                      heff=float(np.median(All["heff"][s_])))
            m[f"Arc{tag}Db{lab}"] = f"{np.median(dA):+.1f}"
            m[f"Arc{tag}Loss{lab}"] = f"{np.median(All['ref'][s_]):.1f}"
            m[f"Arc{tag}Ci{lab}Lo"] = f"{lo:+.0f}"
            m[f"Arc{tag}Ci{lab}Hi"] = f"{hi:+.0f}"
        m[f"Arc{tag}NSum"] = f"{season[(tag, 'Sum')]['n']:d}"
        m[f"Arc{tag}Zsum"] = f"{100 * season[(tag, 'Sum')]['z'][0]:+.0f}"
        m[f"Arc{tag}HeffSum"] = f"{season[(tag, 'Sum')]['heff']:.1f}"
        if "ar" in t:
            m[f"Arc{tag}RefAv"] = f"{t['ar']:.0f}"
            m[f"Arc{tag}HtAv"] = f"{t['ah']:.0f}"
            m[f"Arc{tag}HtAvLo"] = f"{min(t['ah'], t['as'], t['al']):.0f}"
            m[f"Arc{tag}HtAvHi"] = f"{max(t['ah'], t['as'], t['al']):.0f}"
    six = [table[s[0]] for s in ST]
    m["ArcRmsMin"] = f"{100 * min(t['rms'] for t in six):.1f}"
    m["ArcRmsMax"] = f"{100 * max(t['rms'] for t in six):.1f}"
    m["ArcZwinMoveMax"] = f"{100 * max(max(t['zwin'], t['zwin2'], t['zwinl']) - min(t['zwin'], t['zwin2'], t['zwinl']) for t in table.values()):.0f}"
    m["ArcColMin"] = f"{100 * min(t['zv'][0] for t in six):+.0f}"
    m["ArcColMax"] = f"{100 * max(t['zv'][1] for t in six):+.0f}"
    # season scope of the sign: every June to August median, every record rule
    sm = [z for (tag, lab), v in season.items() if lab == "Sum" for z in v["z"]]
    assert max(sm) < 0, "a June to August median is not negative"
    m["ArcZsumMin"] = f"{100 * min(sm):+.0f}"
    m["ArcZsumMax"] = f"{100 * max(sm):+.0f}"
    six_tags = [s[1] for s in ST]
    ns_ = [season[(t_, "Sum")]["n"] for t_ in six_tags]
    m["ArcNSumMin"], m["ArcNSumMax"] = f"{min(ns_):d}", f"{max(ns_):d}"
    hs_ = [v["heff"] for (tag, lab), v in season.items() if lab == "Sum"]
    m["ArcHeffSumMin"], m["ArcHeffSumMax"] = f"{min(hs_):.1f}", f"{max(hs_):.1f}"
    # decibels beside the winter percent range, at the two ends, and the largest summer deviation
    wtags = [k for k in season if k[1] == "Win"]
    end_hi, end_lo = max(wtags, key=lambda k: season[k]["z"][0]), min(wtags, key=lambda k: season[k]["z"][0])
    m["ArcDbWinAtMax"], m["ArcDbWinAtMin"] = f"{season[end_hi]['db']:+.1f}", f"{season[end_lo]['db']:+.1f}"
    m["ArcLossWinAtMin"] = f"{season[end_lo]['ref']:.1f}"
    big_s = max([k for k in season if k[1] == "Sum"], key=lambda k: abs(season[k]["db"]))
    m["ArcDbSumBig"], m["ArcLossSumBig"] = f"{season[big_s]['db']:+.1f}", f"{season[big_s]['ref']:.1f}"
    m["ArcDbSumBigName"] = dict([(s[1], s[2]) for s in ST] + [("Ny", "Ny-{\\AA}lesund")])[big_s[0]]
    # Eureka winter files with few levels in the lowest 3 km
    eu = sid_all == "71917"
    ew_ = np.isin(mon_all[eu], (12, 1, 2))
    lev_ = np.array(table["71917"]["lev"])[ew_]
    z_ = All["zrel"][eu][ew_]
    m["ArcEuCoarseN"], m["ArcEuFineN"] = f"{int((lev_ < 30).sum()):d}", f"{int((lev_ >= 30).sum()):d}"
    m["ArcEuZwinFine"] = f"{100 * np.median(z_[lev_ >= 30]):+.0f}"
    H_all = np.array(pts["heff"] + ny["heff"])
    Z_all = np.array(pts["zrel"] + ny["zrel"])
    P_all = np.array(pts["pred"] + ny["pred"])
    Rv_all = np.array(pts["rival"] + ny["rival"])
    rms = lambda s: 100 * float(np.sqrt(np.mean((P_all[s] - Z_all[s]) ** 2)))
    # validity limit from a stated rule: the lowest edge above which every bin (0.2 km up to 2 km, then
    # coarser) has at least 10 soundings or none, and an rms under 6 points
    edges = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 1e3]

    def limit(H_, P_, Z_):
        lim_ = edges[-2]
        for k in range(len(edges) - 2, -1, -1):
            s_ = (H_ >= edges[k]) & (H_ < edges[k + 1])
            if s_.sum() >= 10 and 100 * np.sqrt(np.mean((P_[s_] - Z_[s_]) ** 2)) < 6.0:
                lim_ = edges[k]
            else:
                break
        return lim_

    lim = limit(H_all, P_all, Z_all)
    # the same rule and the closed-form rms with rho_s from the other two records
    srms = lambda P_, Z_, q_: 100 * float(np.sqrt(np.mean((P_[q_] - Z_[q_]) ** 2)))
    all_six = []
    for lab, hk, pk, zk in (("Two", "he2", "p2", "z2"), ("Low", "hel", "pl", "zl")):
        m[f"ArcHeffValid{lab}"] = f"{limit(All[hk], All[pk], All[zk]):.1f}"
        m[f"ArcPoolRms{lab}"] = f"{srms(All[pk], All[zk], np.ones(len(H_all), bool)):.1f}"
        all_six += [srms(All[pk], All[zk], sid_all == s_) for s_ in order]
    m["ArcHeffValidHi"] = f"{max(lim, limit(All['he2'], All['p2'], All['z2']), limit(All['hel'], All['pl'], All['zl'])):.1f}"
    m["ArcRmsAllMin"] = f"{min(all_six + [100 * t['rms'] for t in six]):.1f}"
    m["ArcRmsAllMax"] = f"{max(all_six + [100 * t['rms'] for t in six]):.1f}"
    abv = [srms(P_all, Z_all, (sid_all == s_) & (H_all >= lim)) for s_ in order]
    m["ArcRmsAboveMin"], m["ArcRmsAboveMax"] = f"{min(abv):.1f}", f"{max(abv):.1f}"
    # the failure below the limit is concentrated in a few soundings
    res = P_all - Z_all
    Q_all = All["Q"]
    below = H_all < lim
    tail = below & (np.abs(res) > 0.15)
    i_big = int(np.argmax(np.abs(res)))
    m["ArcTailCut"] = "15"
    m["ArcTailN"] = f"{int(tail.sum()):d}"
    m["ArcTailRestN"] = f"{int((below & ~tail).sum()):d}"
    m["ArcTailRestRms"] = f"{srms(P_all, Z_all, below & ~tail):.1f}"
    m["ArcQCut"] = "0.85"
    m["ArcTailQN"] = f"{int((tail & (Q_all < 0.85)).sum()):d}"
    m["ArcSpikeRatio"] = f"{All['r1'][i_big] / All['r2'][i_big]:.1f}"
    m["ArcQSmallN"] = f"{int((Q_all < 0.85).sum()):d}"
    m["ArcQSmallRms"] = f"{srms(P_all, Z_all, Q_all < 0.85):.1f}"
    m["ArcQBigRms"] = f"{srms(P_all, Z_all, Q_all >= 0.85):.1f}"
    # what each term of the closed form buys, on the soundings
    hmj = json.load(open(os.path.join(OUT, "heff_model.json")))
    hw_v, H_v = hmj["hw"], hmj["H_PT"]

    def null_rms(H_, sg_):
        Ac_ = All["a2"] - All["gws"] * hw_v + All["gws"] * HM.closed_H(H_all, sg_, H_)
        return 100 * float(np.sqrt(np.mean(((All["a2"] - Ac_) / Ac_ - Z_all) ** 2)))

    m["ArcNullNoScale"] = f"{null_rms(1e12, All['sig']):.1f}"
    m["ArcNullNoSelf"] = f"{null_rms(H_v, 0 * All['sig']):.1f}"
    # zero crossing of the deviation at each station, from a line in 1/h_eff
    zc_ = []
    for s_ in ["01004"] + order:
        c_ = np.polyfit(1 / H_all[sid_all == s_], Z_all[sid_all == s_], 1)
        zc_.append(-c_[0] / c_[1])
    m["ArcZeroMin"], m["ArcZeroMax"] = f"{min(zc_):.2f}", f"{max(zc_):.2f}"
    grid = np.arange(8.0, 60.01, 0.1)
    m["ArcScaleBest"] = f"{grid[int(np.argmin([null_rms(g_, All['sig']) for g_ in grid]))]:.0f}"
    near = (H_all >= lim - 0.2) & (H_all < lim)
    low, mid, high = H_all < lim, (H_all >= lim) & (H_all < 5.0), H_all >= 5.0
    m["ArcHeffValidLo"] = f"{lim:.1f}"
    m["ArcHeffNear"] = f"{lim - 0.2:.1f}"
    m["ArcBinNearRms"] = f"{rms(near):.1f}"
    m["ArcBinNearN"] = f"{int(near.sum()):d}"
    for nm_, s in (("Low", low), ("Mid", mid), ("High", high)):
        m[f"ArcBin{nm_}N"] = f"{int(s.sum()):d}"
        m[f"ArcBin{nm_}Rms"] = f"{rms(s):.1f}"
        m[f"ArcBin{nm_}Bias"] = f"{100 * np.mean(P_all[s] - Z_all[s]):+.1f}"
    m["ArcResMax"] = f"{100 * np.max(np.abs(P_all - Z_all)):.0f}"
    m["ArcPoolN"] = f"{len(H_all):d}"
    m["ArcPoolRms"] = f"{100 * np.sqrt(np.mean((P_all - Z_all) ** 2)):.1f}"
    m["ArcPoolRival"] = f"{100 * np.sqrt(np.mean((Rv_all - Z_all) ** 2)):.1f}"
    # is (A1) the reason? Q on the six other stations
    Q = np.array(pts["Q"]); h_o = np.array(pts["heff"])
    m["ArcQLow"] = f"{np.median(Q[h_o < lim]):.2f}"
    m["ArcQHigh"] = f"{np.median(Q[h_o >= 5]):.2f}"
    Inv = np.array(pts["inv"] + ny["inv"])
    m["ArcInvCorr"] = f"{np.corrcoef(Inv, np.log(H_all))[0, 1]:.2f}"
    m["ArcNTot"] = f"{len(h_o):d}"
    m["ArcHeffMax"] = f"{H_all.max():.1f}"
    m["ArcZMin"] = f"{100 * Z_all.min():+.0f}"
    m["ArcZMax"] = f"{100 * Z_all.max():+.0f}"
    wz = [table[s[0]]["zwin"] for s in ST] + [table["01004"]["zwin"]]
    m["ArcZwinMin"] = f"{100 * min(wz):+.0f}"
    m["ArcZwinMax"] = f"{100 * max(wz):+.0f}"
    hw_ = [t["heffwin"] for t in six]
    m["ArcHeffWinMin"] = f"{min(hw_):.1f}"
    m["ArcHeffWinMax"] = f"{max(hw_):.1f}"
    with open(os.path.join(HERE, "..", "arctic_macros.tex"), "w") as f:
        f.write("% Generated by sim/arctic.py. Do not edit by hand.\n")
        for k in sorted(m):
            f.write(f"\\newcommand{{\\{k}}}{{{A.tex_minus(m[k])}}}\n")
    d["macros"] = m
    json.dump(d, open(os.path.join(OUT, "arctic.json"), "w"), indent=1)
    print(json.dumps(m, indent=1))


if __name__ == "__main__":
    if "--prep" in sys.argv:
        prep(sys.argv[sys.argv.index("--prep") + 1])
    elif "--macros-only" in sys.argv:
        summarize()
    else:
        main()
