#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""Specific-attenuation grids from the independent spectroscopic models.

PyRTlib (Larosa et al. 2024) supplies the Rosenkranz model family and the
IAP RAS MWL24 water-vapor continuum. Each model is evaluated on a coarsened
copy of the P.676-13 Annex 1 layer grid (every 9th boundary) for the P.835-6
Section 1.1 atmosphere scaled to each precipitable water vapor (PWV) column,
and the P.676-13 Annex 1 kernel is evaluated on the same grid so that grid
effects cancel in model-to-model differences. Output: out/models.npz.
"""
import sys
import time

import numpy as np

import thzatm as A
from pyrtlib.absorption_model import H2OAbsModel, O2AbsModel, N2AbsModel
from pyrtlib.rt_equation import RTEquation
from pyrtlib.utils import import_lineshape


def setup(h2o, o2):
    """Select a PyRTlib model and load its line lists once."""
    H2OAbsModel.model = h2o
    O2AbsModel.model = o2
    N2AbsModel.model = "MWL24" if h2o == "MWL24" else o2
    H2OAbsModel.h2oll = import_lineshape("h2oll")
    O2AbsModel.o2ll = import_lineshape("o2ll")
    H2OAbsModel.set_ll()
    O2AbsModel.set_ll()


def grid(fs, T, P, rho):
    """(n_layers, n_freq) dB/km for the currently selected model."""
    e = rho * T / 216.7
    cols = []
    for fi in fs:
        aw, ad = RTEquation.clearsky_absorption(np.asarray(P, float), np.asarray(T, float),
                                                np.asarray(e, float), float(fi))
        cols.append((np.asarray(aw) + np.asarray(ad)).ravel() * A.NP_TO_DB)
    return np.stack(cols, 1)

MODELS = [  # (label, H2O model, O2 model)
    ("R98", "R98", "R98"),
    ("R17", "R17", "R17"),
    ("R22", "R22SD", "R22"),
    ("R24", "R24", "R24"),
    ("MWL24", "MWL24", "R24"),
]
PWVS = [1, 2, 5, 10, 15, 20, 30, 45, 60]
CENTERS = None  # filled from out/windows.json by analysis.py; defaults below


def coarse_layers(step=6):
    h, d = A.annex1_layers()
    top = h + d
    idx = np.arange(0, len(h), step)
    hl = h[idx]
    ht = np.append(hl[1:], top[-1])
    return hl, ht - hl


def main():
    import json
    centers = np.array(json.load(open("out/windows.json"))["centers"])
    spec_f = np.concatenate([np.arange(326.0, 425.0, 2.0), np.arange(448.0, 475.0, 2.0),
                             np.arange(620.0, 917.0, 2.0)])
    hl, dl = coarse_layers(9)
    hm = hl + 0.5 * dl
    out = {"h_lo": hl, "delta": dl, "centers": centers,
           "spec_f": spec_f}
    t0 = time.time()
    # every realizable atmosphere of the six-atmosphere study (round 2), plus
    # the Sec. 1.1 profile at 10 mm for the spectrum of Fig. 5(a)
    import analysis as AN
    cases = AN.realizable_cases()
    out["cases"] = np.array([f"{pr}|{p:g}" for pr, p in cases])
    states = {i: A.state(hm, pwv_mm=p, profile=pr) for i, (pr, p) in enumerate(cases)}
    spec_state = A.state(hm, pwv_mm=10)
    for i, (T, P, rho) in states.items():
        out[f"p676_c_{i}"] = A.gamma_grid(centers, hm, T, P, rho)
    out["p676_s_10"] = A.gamma_grid(spec_f, hm, *spec_state)
    for lab, hm2o, o2 in MODELS:
        setup(hm2o, o2)
        for i, (T, P, rho) in states.items():
            out[f"{lab}_c_{i}"] = grid(centers, T, P, rho)
        out[f"{lab}_s_10"] = grid(spec_f, *spec_state)
        print(f"{lab} done {time.time() - t0:.0f}s", flush=True)
    H = np.load("out/horiz.npz")
    out.update({k: H[k] for k in H.files})
    np.savez_compressed(sys.argv[1] if len(sys.argv) > 1 else "out/models.npz", **out)
    print("saved", time.time() - t0)


def _unused():
    # Horizontal sea-level path at 20 C, 20/50/90 % relative humidity, the
    # conditions of the multipass-cell measurement in the review.
    fh = np.arange(330.0, 490.01, 1.0)
    T0, P0 = 293.15, 1013.25
    es = 6.1121 * np.exp((18.678 - 20.0 / 234.5) * 20.0 / (20.0 + 257.14))  # P.453 eq. (9)
    out["horiz_f"] = fh
    for rh in (20, 50, 90):
        e = es * rh / 100.0
        rho = e * 216.7 / T0
        out[f"horiz_p676_{rh}"] = A.gamma_p676(fh, T0, P0, rho)
        for lab, hm2o, o2 in MODELS:
            out[f"horiz_{lab}_{rh}"] = np.array(
                [float(np.squeeze(A.gamma_pyrtlib(fi, np.array([T0]), np.array([P0]),
                                                  np.array([rho]), hm2o, o2))) for fi in fh])
    np.savez_compressed(sys.argv[1] if len(sys.argv) > 1 else "out/models.npz", **out)
    print("saved", time.time() - t0)


if __name__ == "__main__":
    main()
