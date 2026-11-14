#!/usr/bin/env python3
"""Block-bootstrap interval for the out-of-sample rms of the closed form (Corollary 1).

Read-only on cached outputs: out/heff_model.json (pred, zrel, written by heff_model.py),
out/sonde.json (dates, same order) and data/sondes_nyalesund.npz (the 109 soundings seen
before the closed form was written). Daily soundings are autocorrelated, so the 247
out-of-sample soundings are resampled in blocks of consecutive calendar days (moving-block
bootstrap). Writes ../oosci_macros.tex.
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
BLOCK_DAYS = 7
DRAWS = 5000
SEED = 20261110

h = json.load(open(os.path.join(OUT, "heff_model.json")))["sondes"]
rows = json.load(open(os.path.join(OUT, "sonde.json")))
pred, zrel = np.array(h["pred"]), np.array(h["zrel"])
dates = np.array([r["date"] for r in rows], dtype="datetime64[D]")
seen = {str(x) for x in np.load(os.path.join(HERE, "data", "sondes_nyalesund.npz"))["date"]}
oos = np.array([r["date"] not in seen for r in rows])
e2 = ((pred - zrel) ** 2)[oos]
d = dates[oos]
order = np.argsort(d)
e2, d = e2[order], d[order]
day = (d - d[0]).astype(int)

# blocks: every window of BLOCK_DAYS consecutive calendar days that holds at least one sounding
starts = np.unique(day)
blocks = [np.where((day >= s) & (day < s + BLOCK_DAYS))[0] for s in starts]
rng = np.random.default_rng(SEED)
n = len(e2)
stats = np.empty(DRAWS)
for k in range(DRAWS):
    idx = []
    while len(idx) < n:
        idx.extend(blocks[rng.integers(len(blocks))])
    idx = np.array(idx[:n])
    stats[k] = 100 * np.sqrt(e2[idx].mean())
lo, hi = np.percentile(stats, [2.5, 97.5])
point = 100 * np.sqrt(e2.mean())
print(f"n={n} rms={point:.2f} 95% {BLOCK_DAYS}-day block bootstrap {lo:.2f} to {hi:.2f} ({DRAWS} draws)")
with open(os.path.join(HERE, "..", "oosci_macros.tex"), "w") as f:
    f.write("%% written by sim/oos_block_bootstrap.py from cached outputs; do not edit\n")
    f.write(f"\\newcommand{{\\HmOosRmsLo}}{{{lo:.1f}}}\n")
    f.write(f"\\newcommand{{\\HmOosRmsHi}}{{{hi:.1f}}}\n")
    f.write(f"\\newcommand{{\\HmOosBlockDays}}{{{BLOCK_DAYS}}}\n")
