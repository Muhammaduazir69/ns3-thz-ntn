#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""Data figures for the T-TST paper, from the caches written by analysis.py.

Style: figstyle.py (Times-like serif, CVD-checked palette of the
ieee-transactions figures.md; one color, marker and dash per method in every
figure; red for thresholds only). Figures are drawn at final size (3.5 in
column, 7.16 in page width), so 7 pt text is 7 pt.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figures")
import figstyle as FS
from figstyle import BLUE, ORANGE, GREEN, GOLD, PURPLE, TEAL, RED, GREY, INK, INK2, COL, PAGE
FS.setup()
SHADE = "#dae3f2"
WIN = [("W340", 325.152888, 380.197353), ("W410", 380.197353, 424.763020),
       ("W460", 448.001085, 474.689092), ("W670", 620.700807, 752.033113),
       ("W850", 752.033113, 916.171582)]
WLAB = ["340", "410", "460", "670", "850"]
ELEVS = [10, 15, 20, 30, 45, 60, 90]


def save(fig, name):
    fig.savefig(os.path.join(FIG, name))
    plt.close(fig)
    print("wrote", name)


def fig_spectrum():
    r = np.load(os.path.join(OUT, "ref.npz"))
    s = np.load(os.path.join(OUT, "short.npz"))
    f = r["f"]
    pw = list(r["pwv"])
    j = ELEVS.index(30)
    fig, ax = plt.subplots(figsize=(COL, 2.25))
    import analysis as AN
    for _, _, _, slo, shi in AN.WINDOWS:  # the range searched for each window center
        ax.axvspan(slo, shi, color=SHADE, alpha=0.7, lw=0)
    Z = np.load(os.path.join(OUT, "atm.npz"))
    ax.semilogy(Z["spec_f"], Z["spec_dry"], color=INK, lw=0.9, ls=FS.DASH["longdash"],
                label="Annex 1, high-lat. winter, 2 mm")
    ax.semilogy(f, r["ref"][pw.index(10.0), j], color=BLUE, lw=1.0, label="Annex 1, Sec. 1.1, 10 mm")
    ax.semilogy(f, s["a2"][pw.index(10.0), j], color=ORANGE, lw=1.4, ls=FS.DASH["dash"],
                label="Annex 2 (P.676-13), 10 mm")
    ax.semilogy(f, s["a2r12"][pw.index(10.0), j], color=PURPLE, lw=1.0, ls=FS.DASH["dashdot"],
                label="P.676-12 form, 10 mm")
    cs = json.load(open(os.path.join(OUT, "windows.json")))["centers"]
    for c, lab in zip(cs, WLAB):
        ax.text(c, 1.4, lab, ha="center", va="bottom", fontsize=7, color=INK)
    ax.set_xlim(300, 1000)
    ax.set_ylim(1, 3e3)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Slant gaseous loss at 30° (dB)")
    ax.legend(loc="lower center", ncol=2, fontsize=7, bbox_to_anchor=(0.5, 1.0), columnspacing=0.8)
    save(fig, "fig2_spectrum.pdf")


def atm_window_worst():
    """Largest |deviation| per window over the realizable six-atmosphere set
    (sea level, reference loss at most A_max), per shortcut row of atm.npz."""
    import analysis as AN
    Z = np.load(os.path.join(OUT, "atm.npz"))
    meta = json.load(open(os.path.join(OUT, "atm_meta.json")))
    out = {k: np.full(5, np.nan) for k in AN.ATM_SHORT}
    for key, m in meta.items():
        if key.startswith("native|") or key.split("|")[2] != "0":
            continue
        if m["rh_s_sea"] < AN.RH_FLOOR or m["rh_max"] > 1.0:
            continue
        X = Z[key]
        for si, k in enumerate(AN.ATM_SHORT, start=1):
            for w in range(5):
                D = (X[0][:, w] <= AN.AMAX) & np.isfinite(X[si][:, w])
                if D.any():
                    v = float(np.max(np.abs(X[si][:, w] - X[0][:, w])[D]))
                    out[k][w] = v if np.isnan(out[k][w]) else max(out[k][w], v)
    return out


def fig_errors():
    lines = json.load(open(os.path.join(OUT, "fig_extra.json")))
    aw = atm_window_worst()
    fig, ax = plt.subplots(1, 2, figsize=(PAGE, 2.8))
    x = np.arange(5)
    # (key, label, style key, hollow marker)
    rows = [("Module", "module quadrature", "module", False),
            ("Straight", "straight ray", "straight", True),
            ("Csc", "cosecant law", "csc", False),
            ("AnnexTwo", "Annex 2, heights", "a2", False),
            ("AnnexTwoVt", "Annex 2, column", "a2vt", True),
            ("AnnexTwoRtw", "P.676-12, heights", "a2r12", False),
            ("AnnexTwoVtRtw", "P.676-12, column", "a2vtr12", True)]
    for k, lab, st, hollow in rows:  # Annex 2 of P.676-13 exists only at 340 GHz: offset so it is not hidden
        dx = 0.12 if k in ("AnnexTwo", "AnnexTwoVt") else 0.0
        ax[0].semilogy(x + dx, aw[k], lw=0.9, ms=4, label=lab, **FS.style(st, hollow))
    ax[0].semilogy(x, lines["pseudo"], lw=0.9, ms=4.5, label="no 1780 GHz line", **FS.style("pseudo"))
    ax[0].semilogy(x, lines["prof"], lw=0.9, ms=4.5, label="column-only", **FS.style("prof"))
    ax[0].semilogy(x, lines["moddev"], lw=0.9, ms=5.5, label="vs. R22 model", **FS.style("moddev"))
    ax[0].axhline(0.5, **FS.THRESH)
    ax[0].text(4.15, 0.42, "0.5 dB", color=RED, fontsize=7, ha="right", va="top")
    ax[0].set_xticks(x, WLAB)
    ax[0].set_xlabel("Window center (GHz)")
    ax[0].set_ylabel(r"Max. $|\Delta A|$ (dB)")
    ax[0].set_ylim(1e-2, 100)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, fontsize=7, bbox_to_anchor=(0.5, 1.0),
               columnspacing=1.0, handletextpad=0.4)
    ax[0].set_title("(a) loss deviation", fontsize=8)
    bw = json.load(open(os.path.join(OUT, "bw.json")))
    for k, lab, st, hollow in rows:
        y = [max(abs(v) for v in bw[f"{k}|{w[0]}"]) if bw.get(f"{k}|{w[0]}") else np.nan for w in WIN]
        dx = {"AnnexTwoVt": 0.08, "AnnexTwoVtRtw": 0.08, "Straight": -0.08}.get(k, 0.0)
        ax[1].plot(x + dx, y, lw=0.9, ms=4, **FS.style(st, hollow))
    ax[1].set_xticks(x, WLAB)
    ax[1].set_xlabel("Window center (GHz)")
    ax[1].set_ylabel(r"Max. $|\Delta B_u|$ (GHz)")
    ax[1].set_title("(b) usable-bandwidth deviation", fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.80))
    save(fig, "fig3_errors.pdf")


def fig_heff():
    """Annex 2 zenith deviation at 341.5 GHz against h_eff = V / rho_s: a year at Ny-Alesund, six other
    Arctic stations and the six P.835-6 atmospheres at native column, with Corollary 1."""
    import analysis as AN
    ex = json.load(open(os.path.join(OUT, "fig_extra.json")))["sonde"]
    Z = np.load(os.path.join(OUT, "atm.npz"))
    meta = json.load(open(os.path.join(OUT, "atm_meta.json")))
    fig, ax = plt.subplots(figsize=(COL, 2.3))
    he, zr, gr = np.array(ex["heff"]), 100 * np.array(ex["zrel"]), np.array(ex["group"])
    arc = json.load(open(os.path.join(OUT, "arctic.json")))["points"]
    ax.plot(arc["heff"], 100 * np.array(arc["zrel"]), "x", ms=2.6, color=INK2, mew=0.5, ls="none", label="six other stations")
    ax.plot(he[gr == "mid"], zr[gr == "mid"], ".", ms=2, color="#999999", ls="none", label="Ny-\u00c5lesund, other months")
    ax.plot(he[gr == "win"], zr[gr == "win"], "o", ms=2.6, color=BLUE, mfc="none", mew=0.6, ls="none", label="Ny-\u00c5lesund, Dec. to Feb.")
    ax.plot(he[gr == "sum"], zr[gr == "sum"], "s", ms=2.6, color=ORANGE, mfc="none", mew=0.6, ls="none", label="Ny-\u00c5lesund, June to Aug.")
    first = True
    for pr in AN.ATM_PROFS:
        X = Z[f"native|{pr}"]
        y = 100 * (X[4][-1, 0] - X[0][-1, 0]) / X[0][-1, 0]
        hx = meta[f"native|{pr}"]["heff"]
        ax.plot(hx, y, "D", ms=3.4, color=INK, ls="none", label="P.835-6 atmospheres" if first else None)
        first = False
        if pr == "high-winter":
            ax.annotate("high-latitude winter", (hx, y), xytext=(7.5, 25), textcoords="data", fontsize=7, color=INK,
                        ha="center", va="center", arrowprops=dict(arrowstyle="-", lw=0.4, color=INK2, shrinkA=1, shrinkB=2))
    ax.axhline(0, color=INK, lw=0.5)
    hmj = json.load(open(os.path.join(OUT, "heff_model.json")))
    sg, wt = float(np.median(hmj["sondes"]["sig"])), float(np.median(hmj["sondes"]["wet"]))
    hg = np.geomspace(0.7, 30.0, 200)
    Hc = (1 - sg) / (1 / hg + 1 / hmj["H_PT"]) + sg * hg / 2
    ax.plot(hg, 100 * (hmj["hw"] / Hc - 1) * wt, color=INK, lw=0.9, ls="--", label="Corollary 1")
    ax.set_xscale("log")
    ax.set_xlim(0.7, 30)
    ax.set_xticks([1, 2, 5, 10, 20, 30])
    ax.set_xticklabels(["1", "2", "5", "10", "20", "30"])
    ax.set_xlabel(r"$h_{\mathrm{eff}} = V/\rho_s$ (km)")
    ax.set_ylabel("Annex 2 deviation, zenith (%)")
    ax.legend(loc="upper right", fontsize=6)
    save(fig, "fig_heff.pdf")


def fig_lines():
    cur = json.load(open(os.path.join(OUT, "fig_extra.json")))["linecurves"]
    fig, ax = plt.subplots(figsize=(COL, 2.2))
    sty = [(BLUE, "o", "-"), (ORANGE, "s", FS.DASH["dash"]), (GREEN, "^", FS.DASH["dashdot"]),
           (PURPLE, "D", FS.DASH["dot"]), (GOLD, "v", FS.DASH["longdash"])]
    for (w, *_), (c, m, ls), lab in zip(WIN, sty, WLAB):
        y = np.maximum(np.array(cur[w]), 1e-4)
        ax.semilogy(np.arange(len(y)), y, color=c, lw=0.9, marker=m, ms=2.8, ls=ls, label=f"{lab} GHz")
    ax.axhline(0.5, color=RED, lw=0.7, ls=(0, (1, 1.5)))
    ax.axhline(0.1, color=RED, lw=0.7, ls=(0, (4, 2)))
    n = max(len(v) for v in cur.values())
    ax.text(0.62 * n, 0.56, "0.5 dB", color=RED, fontsize=7, ha="center", va="bottom")
    ax.text(0.62 * n, 0.082, "0.1 dB", color=RED, fontsize=7, ha="center", va="top")
    ax.set_xlim(0, n + 1)
    ax.set_ylim(0.01, 100)
    ax.set_xlabel("Line units kept (greedy order)")
    ax.set_ylabel(r"Worst $|\Delta A|$ over operable domain (dB)")
    ax.legend(loc="upper right", ncol=2)
    save(fig, "fig4_lines.pdf")


def fig_models():
    S = np.load(os.path.join(OUT, "models_slant.npz"))
    M = np.load(os.path.join(OUT, "models.npz"))
    f = S["spec_f"]
    fig, ax = plt.subplots(2, 1, figsize=(COL, 3.6))
    ref = S["spec_R22"]
    seg = np.split(np.arange(len(f)), np.where(np.diff(f) > 3)[0] + 1)
    sty = [("p676", "P.676-13 Annex 1", BLUE, "-"), ("R98", "Rosenkranz 1998", GREY, FS.DASH["dash"]),
           ("R17", "Rosenkranz 2017", GOLD, FS.DASH["dashdot"]), ("R24", "Rosenkranz 2024", PURPLE, FS.DASH["dot"]),
           ("MWL24", "MWL24 continuum", GREEN, FS.DASH["longdash"])]
    for key, lab, c, ls in sty:
        y = 100 * (S[f"spec_{key}"] - ref) / ref
        for k, sg in enumerate(seg):
            ax[0].plot(f[sg], y[sg], color=c, ls=ls, lw=0.9, label=lab if k == 0 else None)
    ax[0].axhline(0, color=INK, lw=0.5)
    ax[0].set_ylabel("vs. R22 (%)")
    ax[0].set_xlabel("Frequency (GHz)")
    ax[0].set_ylim(-64, 48)
    ax[0].set_yticks([-40, -20, 0, 20, 40])
    ax[0].set_xlim(320, 920)
    # inside the empty lower-left part of the axes, clear of every curve (round 3, item 14)
    ax[0].legend(loc="lower left", ncol=2, fontsize=6.8, columnspacing=0.8, borderaxespad=0.3)
    ax[0].set_title("(a) slant path, 30°, 10 mm", fontsize=8)
    fh = M["horiz_f"]
    for key, rh, lab, c, ls in [("p676", 20, "P.676-13, 20% RH", BLUE, FS.DASH["dot"]),
                                ("p676", 50, "P.676-13, 50% RH", BLUE, "-"),
                                ("p676", 90, "P.676-13, 90% RH", BLUE, FS.DASH["dash"]),
                                ("MWL24", 50, "MWL24, 50% RH", GREEN, FS.DASH["longdash"])]:
        y = 100 * (M[f"horiz_{key}_{rh}"] - M[f"horiz_R22_{rh}"]) / M[f"horiz_R22_{rh}"]
        ax[1].plot(fh, y, color=c, ls=ls, lw=0.9, label=lab)
    ax[1].axhline(0, color=INK, lw=0.5)
    ax[1].set_xlabel("Frequency (GHz)")
    ax[1].set_ylabel("vs. R22 (%)")
    ax[1].set_ylim(-15, 15)
    ax[1].legend(loc="lower left", ncol=2, fontsize=7)
    ax[1].set_title("(b) horizontal path, sea level, 20 °C", fontsize=8)
    fig.tight_layout(h_pad=0.4)
    save(fig, "fig5_models.pdf")


def fig_avail():
    import analysis as AN
    from matplotlib.lines import Line2D
    res = dict(np.load(os.path.join(OUT, "avail.npz")))
    cs = json.load(open(os.path.join(OUT, "windows.json")))["centers"]
    fig, ax = plt.subplots(1, 2, figsize=(PAGE, 2.55))
    cols = {"Svalbard": BLUE, "Madrid": ORANGE}
    styles = (("Aref", "-", "reference (Annex 1)"), ("Aa2", FS.DASH["dash"], "Annex 2 (P.676-13)"),
              ("Aa2r12", FS.DASH["dashdot"], "P.676-12 form"))
    thr = np.geomspace(0.1, 200, 200)
    for p, k in enumerate((0, 1)):
        for name in cols:  # Singapore stays at zero in both windows and is omitted
            for which, ls, tag in styles:
                if not np.isfinite(res[f"{name}_{which}"][:, :, k, :]).all():
                    continue
                rr, W = AN.rate_samples(res, name, which, k, cs[k])
                av = np.array([W[rr >= t].sum() for t in thr])
                ax[p].semilogx(thr, 100 * av, color=cols[name], ls=ls, lw=1.0)
        ax[p].axvline(AN.RATE_TARGET, color=RED, lw=0.7)
        ax[p].set_xlim(0.1, 200)
        ax[p].set_ylim(0, 100)
        ax[p].set_xlabel("Clear-air rate threshold (Gbit/s)")
        ax[p].set_ylabel("Visible time above threshold (%)")
        ax[p].set_title(f"({'ab'[p]}) {WLAB[k]} GHz window", fontsize=8)
    h = [Line2D([0], [0], color=c, lw=1.4) for c in cols.values()] + \
        [Line2D([0], [0], color=INK, ls=ls, lw=1.0) for _, ls, _ in styles]
    l = list(cols) + [tag for _, _, tag in styles]
    fig.legend(h, l, loc="upper center", ncol=5, fontsize=7, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, "fig6_avail.pdf")


if __name__ == "__main__":
    import sys
    for name in sys.argv[1:] or ["spectrum", "errors", "lines", "models", "avail", "heff"]:
        globals()[f"fig_{name}"]()
