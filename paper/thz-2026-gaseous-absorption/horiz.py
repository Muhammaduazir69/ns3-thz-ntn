#!/usr/bin/env python3
"""Horizontal sea-level specific attenuation, 330-490 GHz, 20/50/90 % RH at
20 C: P.676-13 and the PyRTlib models (line lists loaded once per model)."""
import numpy as np
from pyrtlib.absorption_model import H2OAbsModel, O2AbsModel, N2AbsModel
from pyrtlib.rt_equation import RTEquation
from pyrtlib.utils import import_lineshape
import thzatm as A
fh = np.arange(330.0, 490.01, 1.0)
T0, P0 = 293.15, 1013.25
es = 6.1121 * np.exp((18.678 - 20.0 / 234.5) * 20.0 / (20.0 + 257.14))
out = {"horiz_f": fh}
for rh in (20, 50, 90):
    rho = es * rh / 100.0 * 216.7 / T0
    out[f"horiz_p676_{rh}"] = A.gamma_p676(fh, T0, P0, rho)
for lab, h2o, o2 in [("R98", "R98", "R98"), ("R17", "R17", "R17"), ("R22", "R22SD", "R22"),
                     ("R24", "R24", "R24"), ("MWL24", "MWL24", "R24")]:
    H2OAbsModel.model = h2o; O2AbsModel.model = o2
    N2AbsModel.model = "MWL24" if h2o == "MWL24" else o2
    H2OAbsModel.h2oll = import_lineshape("h2oll"); O2AbsModel.o2ll = import_lineshape("o2ll")
    H2OAbsModel.set_ll(); O2AbsModel.set_ll()
    for rh in (20, 50, 90):
        rho = es * rh / 100.0 * 216.7 / T0
        e = rho * T0 / 216.7
        v = []
        for fi in fh:
            aw, ad = RTEquation.clearsky_absorption(np.array([P0]), np.array([T0]), np.array([e]), fi)
            v.append(float(np.squeeze(aw[0] + ad[0])) * A.NP_TO_DB)
        out[f"horiz_{lab}_{rh}"] = np.array(v)
    print(lab, flush=True)
np.savez_compressed("out/horiz.npz", **out)
