#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""Ny-Alesund (WMO 01004) radiosondes, 12 UTC, from the University of Wyoming
archive (https://weather.uwyo.edu/wsgi/sounding, TEXT:CSV, BUFR source),
resampled to a fixed height grid and stored in data/sondes_nyalesund.npz.

Usage: python -I prep_sondes.py <dir with nya_YYYY-MM-DD_12.csv files> [--year]
Without --year, the December to February and July subset goes to
data/sondes_nyalesund.npz (the held-out test of the absorbing-height analysis);
with --year, all files go to data/sondes_nyalesund_year.npz (the sounding rule of
the availability study).
Only numeric columns are read. A sounding is kept if it reaches 300 hPa with
humidity reported to that level.
"""
import csv
import glob
import os
import sys

import numpy as np

ZG = np.concatenate([np.arange(0.0, 3.0, 0.025), np.arange(3.0, 10.0, 0.1), np.arange(10.0, 30.01, 0.5)])


def surface_rules(z, rho):
    """Alternatives to the first record for rho_s: the second record, and the mean over the lowest 100 m."""
    zz = z[0] + np.arange(0.0, 0.1001, 0.005)
    return float(rho[1]), float(np.mean(np.interp(zz, z, rho)))


def read(path):
    rows = []
    with open(path, newline="") as f:
        r = csv.reader(f)
        head = next(r, None)
        if not head or "pressure_hPa" not in head:
            return None
        ip, iz, it, iw = (head.index(k) for k in ("pressure_hPa", "geopotential height_m",
                                                   "temperature_C", "mixing ratio_g/kg"))
        for row in r:
            try:
                rows.append([float(row[ip]), float(row[iz]), float(row[it]), float(row[iw])])
            except (ValueError, IndexError):
                continue
    if len(rows) < 50:
        return None
    a = np.array(rows)
    a = a[np.argsort(-a[:, 0])]                   # by decreasing pressure
    a = a[np.concatenate([[True], np.diff(a[:, 1]) > 0])]  # strictly rising
    if a[-1, 0] > 300.0:
        return None
    P, zgp, T, w = a[:, 0], a[:, 1] / 1000.0, a[:, 2] + 273.15, a[:, 3]
    z = 6356.766 * zgp / (6356.766 - zgp)        # geopotential -> geometric km
    e = P * w / (621.97 + w)                       # hPa
    rho = e * 216.7 / T
    return z, P, T, rho


def main(src, year=False):
    out = {"zg": ZG}
    dates, Z0, Ps, Ts, Rs, tops, R2, R100 = [], [], [], [], [], [], [], []
    for pth in sorted(glob.glob(os.path.join(src, "nya_*_12.csv"))):
        if not year and os.path.basename(pth)[9:11] not in ("12", "01", "02", "07"):
            continue
        s = read(pth)
        if s is None:
            continue
        z, P, T, rho = s
        zz = z[0] + ZG
        m = zz <= z[-1]
        Ps.append(np.where(m, np.exp(np.interp(zz, z, np.log(P))), np.nan))
        Ts.append(np.where(m, np.interp(zz, z, T), np.nan))
        Rs.append(np.where(m, np.interp(zz, z, rho), np.nan))
        Z0.append(z[0])
        a_, b_ = surface_rules(z, rho)
        R2.append(a_)
        R100.append(b_)
        tops.append(z[-1])
        dates.append(os.path.basename(pth)[4:14])
    out.update(date=np.array(dates), z0=np.array(Z0), P=np.array(Ps), T=np.array(Ts),
               rho=np.array(Rs), top=np.array(tops), rs2=np.array(R2), rs100=np.array(R100))
    dst = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                       "sondes_nyalesund_year.npz" if year else "sondes_nyalesund.npz")
    np.savez_compressed(dst, **out)
    print(len(dates), "soundings ->", dst)


if __name__ == "__main__":
    main(sys.argv[1], "--year" in sys.argv)
