#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""Fig. 1: slant geometry (left) and the model chain whose shortcuts are scored (right).

House diagram style (deng_style.py): sans-serif bold labels, rounded role boxes, a legend
strip across the top, orthogonal lines. Geometry is schematic and not to scale. Every box
names something the study runs; nothing here is a placeholder.
Run: sim/.venv/bin/python sim/make_fig1.py  ->  figures/fig1_stack.pdf
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Polygon, Rectangle, Wedge, Arc

import deng_style as D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "fig1_stack.pdf")

W, H = 7.16, 3.05
D.setup(base=6.6)
plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                     "mathtext.bf": "Arial:bold"})
D.LINES.update({
    "chain": ("#a6b4cc", "-", 3.0, "-|>"),
    "check": (D.LINE_ORAN, "-", 1.2, "-|>"),
    "short": (D.LINE_NEW, "-", 1.2, "-|>"),
    "ray": (D.NEARRT, "-", 1.6, "-"),
    "los": (D.DASH_NAVY, "--", 0.9, "-"),
})

fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")

# ==== legend strip
D.legend_strip(ax, 0.10, H - 0.11, [("chain", "Reference chain"), ("check", "Independent check"),
                                    ("short", "Shortcut or alternative scored"),
                                    ("ray", "Refracted ray"), ("los", "Straight line of sight")],
               gap=0.16, fs=6.6)

# ==== left: geometry
PX0, PY0, PW, PH = 0.05, 0.05, 2.22, H - 0.30
D.panel(ax, PX0, PY0, PW, PH, "Slant geometry (not to scale)", fs=6.8)
ax.add_patch(Rectangle((PX0 + 0.01, PY0 + 0.01), PW - 0.02, PH - 0.30, fc="white", ec="none", zorder=1))
cx, cy, R = 1.10, -13.3, 13.75          # Earth (drawn as a large circle; only its cap shows)
clip = Rectangle((PX0 + 0.01, PY0 + 0.01), PW - 0.02, PH - 0.30, transform=ax.transData)
for dr, a in ((0.95, 0.25), (0.62, 0.45), (0.30, 0.75)):   # graded atmosphere, 100 km drawn thick
    c = Circle((cx, cy), R + dr, fc=D.PALE_BLUE2 if a < 0.5 else D.PALE_BLUE, ec="none", alpha=a, zorder=1.2)
    ax.add_patch(c)
    c.set_clip_path(clip)
top = Circle((cx, cy), R + 0.95, fc="none", ec=D.PANEL_EDGE, lw=0.7, ls=(0, (3, 2)), zorder=1.3)
ax.add_patch(top)
top.set_clip_path(clip)
earth = Circle((cx, cy), R, fc=D.GREEN_PALE, ec="#8cc56a", lw=0.8, zorder=1.4)
ax.add_patch(earth)
earth.set_clip_path(clip)


def surf(x, dh=0.0):
    return cy + np.sqrt((R + dh) ** 2 - (x - cx) ** 2)


# ground terminal: dish on a mast
gx = 0.36
gy = surf(gx)
ax.add_patch(Rectangle((gx - 0.018, gy), 0.036, 0.13, fc=D.INK, ec="none", zorder=6))
ax.add_patch(Wedge((gx + 0.02, gy + 0.16), 0.085, 200, 20, fc=D.SMO, ec=D.SMO, zorder=6))
tx, ty = gx + 0.05, gy + 0.19
# elevated site: a mountain with a smaller dish
mx = 1.86
my = surf(mx)
ax.add_patch(Polygon([(mx - 0.20, my - 0.02), (mx, my + 0.20), (mx + 0.22, surf(mx + 0.22) - 0.02)],
                     closed=True, fc="#c9b28c", ec="#a68a5b", lw=0.6, zorder=5))
ax.add_patch(Wedge((mx + 0.01, my + 0.24), 0.05, 200, 20, fc=D.SMO, ec=D.SMO, zorder=6))
# HAPS at 20 km: a small airship
hx, hy = 1.08, surf(1.08, 0.36)
ax.add_patch(Ellipse((hx, hy), 0.26, 0.085, fc=D.NEARRT_PALE, ec=D.NEARRT_EDGE, lw=0.7, zorder=6))
ax.add_patch(Rectangle((hx - 0.03, hy - 0.07), 0.06, 0.03, fc=D.NEARRT_EDGE, ec="none", zorder=6))
# satellite at 550 km
sx, sy = 1.86, PY0 + PH - 0.62
ax.add_patch(Rectangle((sx - 0.055, sy - 0.05), 0.11, 0.10, fc=D.RT, ec=D.RT, zorder=7))
for s in (-1, 1):
    ax.add_patch(Rectangle((sx + s * 0.075 - (0.16 if s < 0 else 0), sy - 0.03), 0.16, 0.06,
                           fc=D.DU_PALE, ec=D.DU_EDGE, lw=0.6, zorder=7))
# rays: the refracted ray leaves above the line of sight and bends toward the ground
ax.plot([tx, sx - 0.07], [ty, sy - 0.05], color=D.DASH_NAVY, lw=0.9, ls=(0, (3, 2)), zorder=4)
t = np.linspace(0, 1, 60)
p0, p1, p2 = np.array([tx, ty]), np.array([tx + 0.40, ty + 0.98]), np.array([sx - 0.07, sy - 0.05])
bez = np.outer((1 - t) ** 2, p0) + np.outer(2 * (1 - t) * t, p1) + np.outer(t ** 2, p2)
ax.plot(bez[:, 0], bez[:, 1], color=D.NEARRT, lw=1.6, zorder=4.5, solid_capstyle="round")
# elevation angle
ax.plot([tx, tx + 0.48], [ty, ty], color=D.INK, lw=0.5, zorder=4)
ang = np.degrees(np.arctan2(sy - 0.05 - ty, sx - 0.07 - tx))
ax.add_patch(Arc((tx, ty), 0.62, 0.62, theta1=0, theta2=ang, color=D.INK, lw=0.5, zorder=4))
ax.text(tx + 0.34, ty + 0.07, r"$\theta$", fontsize=7.5, color=D.INK, ha="left", va="bottom", zorder=8)

kw = dict(fontsize=6.4, weight="bold", color=D.INK, zorder=9)
ax.text(sx, sy + 0.10, "LEO, 550 km", ha="center", va="bottom", **kw)
ax.text(hx + 0.17, hy, "HAPS, 20 km", ha="left", va="center", **kw)
ax.text(PX0 + PW - 0.05, my - 0.10, "2.5 and 5 km sites", ha="right", va="top", **kw)
ax.text(gx - 0.12, gy - 0.10, "Sea level", ha="left", va="top", **kw)
ax.text(PX0 + PW - 0.06, surf(PX0 + PW - 0.30, 0.95) + 0.05, "Top of layers, 100 km", ha="right",
        va="bottom", **kw)

# ==== right: model chain
X0 = 2.38
cw = {"lvl": 0.74, "chk": 1.04, "ref": 1.30, "sc": 1.40}
gap = 0.10
xl = X0
xc = xl + cw["lvl"] + 0.07
xr = xc + cw["chk"] + gap
xs = xr + cw["ref"] + gap
hdr_y = H - 0.36
hk = dict(fontsize=6.4, weight="bold", color=D.INK, ha="center", va="center")
ax.text(xl + cw["lvl"] / 2, hdr_y, "Level", **hk)
ax.text(xc + cw["chk"] / 2, hdr_y, "Independent checks", **hk)
ax.text(xr + cw["ref"] / 2, hdr_y, "Reference", **hk)
ax.text(xs + cw["sc"] / 2, hdr_y, "Scored against it", **hk)

rh, rg = 0.58, 0.12
r_top = hdr_y - 0.16
rows = [
    ("Atmosphere", ("Independent code", "ITU-Rpy P.835-6\nprofiles"),
     ("ITU-R P.835-6", "six atmospheres,\nrealizable columns"),
     ("Profile shortcut", "column-only profile")),
    ("Specific\nattenuation", ("Codes and data", "ITU-Rpy, am,\npublished measurements"),
     ("P.676-13 Annex 1", "44 O$_2$ and 35 H$_2$O lines,\ndry continuum"),
     ("Line sets and models", "no 1780 GHz line, pruned sets,\nRosenkranz models, MWL24")),
    ("Slant path", ("Independent codes", "ITU-Rpy ray tracing,\nam on slant paths"),
     ("Refracted ray", "922 layers,\napparent elevation"),
     ("Path shortcuts", "straight ray, module quadrature,\ncosecant law, Annex 2 (P.676-13),\nP.676-12 form")),
]
ys = []
for i, (lvl, chk, ref, sc) in enumerate(rows):
    y = r_top - rh - i * (rh + rg)
    ys.append(y)
    D.box(ax, xl, y, cw["lvl"], rh, lvl, role="smo", fs=7.0)
    D.box(ax, xc, y, cw["chk"], rh, chk[0], role="green_pale", sub=chk[1], fs=7.4, sfs=6.7)
    D.box(ax, xr, y, cw["ref"], rh, ref[0], role="du", sub=ref[1], fs=7.4, sfs=6.7)
    D.box(ax, xs, y, cw["sc"], rh, sc[0], role="card", sub=sc[1], fs=7.4, sfs=6.5)
    ym = y + rh / 2
    D.line(ax, [(xc + cw["chk"], ym), (xr, ym)], "check")
    D.line(ax, [(xr + cw["ref"], ym), (xs, ym)], "short")
# metrics row
my_ = ys[-1] - rg - 0.40
D.box(ax, xl, my_, cw["lvl"], 0.40, "Metrics", role="smo", fs=7.0)
D.box(ax, xc, my_, xs + cw["sc"] - xc, 0.40, "Deviation from the reference at the 340 to 850 GHz windows",
      role="nonrt_pale", sub=r"loss $\Delta A$ (dB and %), usable bandwidth $B_u$, clear-air rate at availability",
      fs=7.4, sfs=6.7)
# reference chain: block arrows down the reference column
xm = xr + cw["ref"] / 2
for i in range(len(ys)):
    y_from = ys[i]
    y_to = ys[i + 1] + rh if i + 1 < len(ys) else my_ + 0.40
    D.line(ax, [(xm, y_from), (xm, y_to)], "chain")

fig.savefig(OUT)
print("wrote", OUT)
