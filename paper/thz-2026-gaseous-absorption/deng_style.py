#!/usr/bin/env python3
"""House diagram style adopted from Deng et al., IEEE OJ-COMS 2026 (see references/style-deng2026.md).

Use:
    import deng_style as D
    D.setup()                               # sans-serif bold labels, embedded TrueType fonts
    D.box(ax, x, y, w, h, "Near-RT RIC", role="nearrt")
    D.line(ax, [(x0, y0), (x1, y0), (x1, y1)], kind="oran", label="E2")
    D.legend_strip(ax, x, y, [("oran", "Existing O-RAN interface"), ("new", "New interface")])
    print(D.latex_preamble())               # colors + \\cellH/\\cellM/\\cellL for the comparison table
Run `python deng_style.py demo.pdf` to render a sample.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# framework view (white bold text on saturated fills)
SMO = "#3a71ae"
AI_CTRL = "#69c8bd"
NONRT = "#75bd42"
NEARRT = "#ee822f"
RT = "#3d7199"
NF = "#91ace0"
TPN = "#9b78b8"
CN = "#cd7e8a"
# reference-architecture view (dark text on pastels)
SMO_REF = "#9b5eb3"
NONRT_PALE = "#fce6d5"
NEARRT_PALE, NEARRT_EDGE = "#f5b482", "#f0934b"
Y1_PALE = "#fee695"
GREEN_PALE = "#acd78e"
DU_PALE, DU_EDGE = "#b6c7ea", "#8fa9df"
# lines, borders, fills
LINE_ORAN = "#35b44a"
LINE_NEW = "#e54c5e"
LINE_FH = "#4771c3"
ATTACH = "#ee822f"
DASH_NAVY = "#2e54a1"
PANEL_EDGE = "#808080"
CREAM = "#fef1c7"
PALE_BLUE = "#dae3f2"
PALE_BLUE2 = "#ebf5fe"
INK = "#1a1a1a"
# comparison-table cells
CELL_H, CELL_M, CELL_L = "#ffe599", "#bdd7ee", "#c6efce"

ROLES = {  # role: (fill, edge, text color)
    "smo": (SMO, SMO, "white"), "ai": (AI_CTRL, AI_CTRL, "white"), "nonrt": (NONRT, NONRT, "white"),
    "nearrt": (NEARRT, NEARRT, "white"), "rt": (RT, RT, "white"), "nf": (NF, NF, INK),
    "tpn": (TPN, TPN, "white"), "cn": (CN, CN, "white"),
    "smo_ref": (SMO_REF, SMO_REF, "white"), "nonrt_pale": (NONRT_PALE, NEARRT_EDGE, INK),
    "nearrt_pale": (NEARRT_PALE, NEARRT_EDGE, INK), "y1": (Y1_PALE, "#e8c75a", INK),
    "green_pale": (GREEN_PALE, "#8cc56a", INK), "du": (DU_PALE, DU_EDGE, INK),
    "card": (CREAM, "#f3d98a", INK), "panel": ("white", PANEL_EDGE, INK), "pale": (PALE_BLUE, "#b9c8e3", INK),
}
LINES = {  # kind: (color, linestyle, width, arrowstyle)
    "oran": (LINE_ORAN, "-", 1.3, "-"), "new": (LINE_NEW, "-", 1.3, "-"), "fh": (LINE_FH, "-", 1.3, "-"),
    "attach": (ATTACH, "-", 1.2, "<|-|>"), "flow": ("#a6b4cc", "-", 3.0, "-|>"), "control": (INK, "--", 1.0, "-|>"),
}


def setup(base=7.0):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "font.size": base, "font.weight": "bold", "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": 0.8, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    })


def box(ax, x, y, w, h, title, role="nf", sub=None, dashed=False, fs=7.2, z=5, radius=None, sfs=None):
    """rounded flat box; title bold centered (raised when there is a subtitle)"""
    fc, ec, tc = ROLES[role]
    r = radius if radius is not None else 0.12 * min(w, h)
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=%g" % r, fc=fc, ec=DASH_NAVY if dashed else ec,
                       lw=1.2 if dashed else 0.8, ls=(0, (4, 2)) if dashed else "-", zorder=z)
    ax.add_patch(p)
    three = bool(sub) and sub.count("\n") >= 2      # three-line subtitle: lift the title, lower the block
    ty = y + h * (0.78 if three else 0.64 if sub else 0.5)
    ax.text(x + w / 2, ty, title, ha="center", va="center", color=tc, fontsize=fs, weight="bold", zorder=z + 1)
    if sub:
        ax.text(x + w / 2, y + h * (0.34 if three else 0.30), sub, ha="center", va="center", color=tc, fontsize=sfs or fs - 1.3, weight="normal",
                zorder=z + 1, linespacing=1.15)
    return p


def panel(ax, x, y, w, h, title, fs=7.4):
    """domain panel: white, thin grey border, bold title top-left"""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="square,pad=0", fc="white", ec=PANEL_EDGE, lw=0.8, zorder=1))
    ax.text(x + 0.06, y + h - 0.06, title, ha="left", va="top", fontsize=fs, weight="bold", color=INK, zorder=2)


def line(ax, pts, kind="oran", label=None, label_at=0, dx=0.04, dy=0.04, fs=6.6, z=3):
    """orthogonal polyline through pts; the last segment carries the arrow head of the kind"""
    col, ls, lw, arr = LINES[kind]
    for a, b in zip(pts[:-2], pts[1:-1]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=col, lw=lw, ls=ls, zorder=z, solid_capstyle="butt")
    ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle=arr, mutation_scale=7, color=col, lw=lw, ls=ls,
                                 shrinkA=0, shrinkB=0, zorder=z))
    if label:
        a, b = pts[label_at], pts[label_at + 1]
        ax.text((a[0] + b[0]) / 2 + dx, (a[1] + b[1]) / 2 + dy, label, fontsize=fs, weight="bold", color=INK,
                ha="left", va="bottom", zorder=z + 2)


def legend_strip(ax, x, y, items, gap=0.20, fs=6.8):
    """one row of line samples + bold labels; items = [(kind, text), ...]; returns end x"""
    from matplotlib.textpath import TextPath
    for kind, text in items:
        col, ls, lw, _ = LINES[kind]
        ax.plot([x, x + 0.28], [y, y], color=col, lw=lw + 0.3, ls=ls)
        t = ax.text(x + 0.34, y, text, fontsize=fs, weight="bold", va="center", ha="left", color=INK)
        width = TextPath((0, 0), text, size=fs).get_extents().width / 72.0
        x += 0.34 + width * 1.08 + gap
    return x


def latex_preamble():
    """colors and cell macros for the comparison table (needs \\usepackage[table]{xcolor})"""
    return r"""\definecolor{cellH}{HTML}{FFE599}\definecolor{cellM}{HTML}{BDD7EE}\definecolor{cellL}{HTML}{C6EFCE}
\newcommand{\cH}{\cellcolor{cellH}H}\newcommand{\cM}{\cellcolor{cellM}M}\newcommand{\cL}{\cellcolor{cellL}L}
\newcommand{\cY}{\cellcolor{cellL}\checkmark}\newcommand{\cP}{\cellcolor{cellM}P}\newcommand{\cN}{--}
\newcommand{\rot}[1]{\rotatebox{90}{\textbf{#1}\,}}"""


def demo(path):
    setup()
    fig = plt.figure(figsize=(3.5, 2.6))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 3.5); ax.set_ylim(0, 2.6); ax.axis("off")
    legend_strip(ax, 0.08, 2.48, [("oran", "O-RAN interface"), ("new", "New interface")], fs=6.2)
    panel(ax, 0.05, 0.05, 3.4, 2.3, "RAN Domain")
    box(ax, 0.20, 1.80, 3.10, 0.36, "SMO / Non-RT RIC", role="smo")
    box(ax, 0.90, 1.10, 1.70, 0.40, "Near-RT RIC", role="nearrt", sub="xApps")
    box(ax, 0.30, 0.25, 1.20, 0.40, "O-DU", role="du")
    box(ax, 2.00, 0.25, 1.20, 0.40, "O-CU", role="nf", dashed=True)
    line(ax, [(1.75, 1.80), (1.75, 1.50)], "oran", label="A1")
    line(ax, [(1.40, 1.10), (1.40, 0.85), (0.90, 0.85), (0.90, 0.65)], "oran", label="E2", label_at=1)
    line(ax, [(2.10, 1.10), (2.10, 0.85), (2.60, 0.85), (2.60, 0.65)], "new")
    fig.savefig(path)


if __name__ == "__main__":
    demo(sys.argv[1] if len(sys.argv) > 1 else "deng_demo.pdf")
