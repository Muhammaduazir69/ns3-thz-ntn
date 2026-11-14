#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""Slant-path gaseous absorption engines for the T-TST paper (paper 4).

Everything here is a numerical model; nothing is fitted to the paper's
results. Four groups of code:

1. ITU-R P.835-6 Section 1.1 reference atmosphere, an exact port of
   ThzNtnMolecularAbsorption::GetAtmosphericConditions in contrib/thz-ntn,
   plus a surface water-vapor density parameter so the column can be swept.
2. ITU-R P.676-13 Annex 1 specific attenuation (44 O2 + 35 H2O lines and the
   dry continuum), vectorized over frequency, with optional line masks.
3. Slant-path integrators: (a) the Annex 1 Section 2.2 refracted ray on the
   exponentially spaced layers of eqs. (14)-(15); (b) the same layers with the
   refractive index forced to 1 (straight ray); (c) an exact mirror of the
   thz-ntn slab integrator; (d) the flat-Earth cosecant law.
4. Wrappers around two independent codes: ITU-Rpy (P.676 Annex 1 and Annex 2)
   and PyRTlib (the Rosenkranz model family and the IAP RAS MWL24 continuum).
"""
from __future__ import annotations

import math
import os
import warnings

import numpy as np

R_E = 6371.0  # km, P.676 Annex 1 Section 2.2 and the thz-ntn module
NP_TO_DB = 10.0 / math.log(10.0)

_HERE = os.path.dirname(os.path.abspath(__file__))
# P.676-13 Annex 1 Tables 1 and 2, vendored in data/ (see data/README.md for
# the source and SHA-256); identical to the tables compiled into thz-ntn and
# to the ITU-Rpy 0.4.0 revision-12 files.
_LINES = os.path.join(_HERE, "data", "p676_13_annex1_lines.json")


def _load_lines():
    import json
    d = json.load(open(_LINES))
    return np.array(d["O2"], float), np.array(d["H2O"], float)


O2, H2O = _load_lines()  # same tables the C++ model compiles in


_TEX_HYPHEN_MINUS = None


def tex_minus(s):
    """Typeset the sign of every negative number in a macro value as a true
    minus, \\ensuremath{-}, so a table cell or sentence never prints a hyphen
    ("-12 to -11"). Math segments ($...$), en-dash ranges (7.5--15) and
    hyphenated words (p835-std) are left alone. Only the sign glyph changes."""
    global _TEX_HYPHEN_MINUS
    import re
    if _TEX_HYPHEN_MINUS is None:
        _TEX_HYPHEN_MINUS = re.compile(r"(?<![\w\-.{^_])-(?=\.?\d)")
    seg = str(s).split("$")
    seg[0::2] = [_TEX_HYPHEN_MINUS.sub(r"\\ensuremath{-}", x) for x in seg[0::2]]
    return "$".join(seg)

# ---------------------------------------------------------------------------
# 1. ITU-R P.835-6 Section 1.1 (port of the C++ closed form)
# ---------------------------------------------------------------------------


def p835_TP(h_km):
    """Temperature (K) and total pressure (hPa) at geometric height h (km)."""
    h = np.maximum(np.asarray(h_km, float), 0.0)
    hp = 6356.766 * h / (6356.766 + h)
    T = np.select(
        [hp <= 11, hp <= 20, hp <= 32, hp <= 47, hp <= 51, hp <= 71, hp <= 84.852],
        [288.15 - 6.5 * hp, 216.65 + 0 * hp, 216.65 + (hp - 20), 228.65 + 2.8 * (hp - 32),
         270.65 + 0 * hp, 270.65 - 2.8 * (hp - 51), 214.65 - 2.0 * (hp - 71)],
        186.946)
    with np.errstate(all="ignore"):
        P = np.select(
            [hp <= 11, hp <= 20, hp <= 32, hp <= 47, hp <= 51, hp <= 71, hp <= 84.852],
            [1013.25 * (288.15 / (288.15 - 6.5 * hp)) ** (-34.1632 / 6.5),
             226.3226 * np.exp(-34.1632 * (hp - 11.0) / 216.65),
             54.74980 * (216.65 / (216.65 + (hp - 20.0))) ** 34.1632,
             8.680422 * (228.65 / (228.65 + 2.8 * (hp - 32.0))) ** (34.1632 / 2.8),
             1.109106 * np.exp(-34.1632 * (hp - 47.0) / 270.65),
             0.6694167 * (270.65 / (270.65 - 2.8 * (hp - 51.0))) ** (-34.1632 / 2.8),
             0.03956649 * (214.65 / (214.65 - 2.0 * (hp - 71.0))) ** (-34.1632 / 2.0)],
            0.00373405)  # continuity with the 71-84.852 km segment
    return T, P


def p835_rho(h_km, rho0=7.5, scale_km=2.0):
    """Water-vapor density (g/m3) with the 2e-6 mixing-ratio floor."""
    h = np.maximum(np.asarray(h_km, float), 0.0)
    T, P = p835_TP(h)
    rho = rho0 * np.exp(-h / scale_km)
    e = rho * T / 216.7
    # P.835-6 Sec. 1.1: the volume mixing ratio e/P is floored at 2e-6
    # (thz-ntn THZ-15).
    return np.where(e / P < 2e-6, 2e-6 * P * 216.7 / T, rho)


def p835_seasonal(h_km, kind):
    """ITU-R P.835-6 Section 2 seasonal profiles via ITU-Rpy (T, P, rho)."""
    import itur.models.itu835 as m
    # ITU-Rpy selects low (<22 deg), mid (<45 deg) or high latitude by |lat|.
    lat = {"low": 10.0, "mid-summer": 30.0, "mid-winter": 30.0,
           "high-summer": 60.0, "high-winter": 60.0}[kind]
    season = "winter" if kind.endswith("winter") else "summer"
    h = np.asarray(h_km, float)
    shape = np.shape(h)
    T = np.asarray(m.temperature(lat, h, season=season).value, float).reshape(shape)
    P = np.asarray(m.pressure(lat, h, season=season).value, float).reshape(shape)
    rho = np.asarray(m.water_vapour_density(lat, h, season=season).value, float).reshape(shape)
    return T, P, rho


# ---------------------------------------------------------------------------
# Layer grid of P.676-13 Annex 1 Section 2.2, eqs. (14)-(15)
# ---------------------------------------------------------------------------


def annex1_layers(n=922, h0=0.0):
    i = np.arange(n)
    delta = 1e-4 * np.exp(i / 100.0)                      # layer thickness, km
    h_lo = 1e-4 * (np.exp(i / 100.0) - 1.0) / (np.exp(0.01) - 1.0)
    keep = h_lo + h0 < 100.0
    return h_lo[keep] + h0, delta[keep]


def annex1_layers_from(h_lower, h_upper=100.0):
    """Layer grid of P.676-13 Annex 1 eqs. (16a)-(16d) for a path that starts
    at h_lower (km), e.g. an elevated ground station or a HAPS terminal."""
    if h_lower <= 0.0:
        return annex1_layers()
    e1, e2 = math.exp(0.01), math.exp(0.02)
    il = math.floor(100 * math.log(1e4 * h_lower * (e1 - 1) + 1) + 1)
    iu = math.ceil(100 * math.log(1e4 * h_upper * (e1 - 1) + 1) + 1)
    m = (e2 - e1) / (math.exp(iu / 100) - math.exp(il / 100)) * (h_upper - h_lower)
    i = np.arange(il, iu)                                  # layers il .. iu-1
    delta = m * np.exp((i - 1) / 100.0)
    h_lo = h_lower + m * (np.exp((i - 1) / 100.0) - np.exp((il - 1) / 100.0)) / (e1 - 1)
    return h_lo, delta


def column_mm(h_lo, delta, rho):
    """Integrated water vapor (kg/m2 = mm), rho in g/m3 at layer midpoints."""
    return float(np.sum(rho * delta))  # g/m3 * km = kg/m2


def rho0_for_pwv(pwv_mm, scale_km=2.0):
    """Surface density giving the requested column on the Annex 1 layers."""
    h_lo, d = annex1_layers()
    hm = h_lo + 0.5 * d
    lo, hi = 1e-4, 200.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if column_mm(h_lo, d, p835_rho(hm, mid, scale_km)) < pwv_mm:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# 2. ITU-R P.676-13 Annex 1 specific attenuation
# ---------------------------------------------------------------------------


def gamma_p676(f, T, Ptot, rho, o2_mask=None, h2o_mask=None, dry_cont=True,
               split=False):
    """Specific attenuation (dB/km), f (GHz) array; scalar T, Ptot, rho.

    Ptot is TOTAL pressure; the dry partial pressure p = Ptot - e enters the
    oxygen terms, as the recommendation defines (and as thz-ntn does).
    """
    f = np.atleast_1d(np.asarray(f, float))[:, None]
    th = 300.0 / T
    e = max(rho, 0.0) * T / 216.7
    p = max(Ptot - e, 0.0)
    o = O2 if o2_mask is None else O2[o2_mask]
    w = H2O if h2o_mask is None else H2O[h2o_mask]
    # oxygen
    f0, a1, a2, a3, a4, a5, a6 = (o[:, k][None, :] for k in range(7))
    S = a1 * 1e-7 * p * th ** 3 * np.exp(a2 * (1 - th))
    df = a3 * 1e-4 * (p * th ** (0.8 - a4) + 1.1 * e * th)
    df = np.sqrt(df * df + 2.25e-6)
    dl = (a5 + a6 * th) * 1e-4 * (p + e) * th ** 0.8
    F = f / f0 * ((df - dl * (f0 - f)) / ((f0 - f) ** 2 + df ** 2)
                  + (df - dl * (f0 + f)) / ((f0 + f) ** 2 + df ** 2))
    Nox = (S * F).sum(axis=1)
    # water vapor
    f0, b1, b2, b3, b4, b5, b6 = (w[:, k][None, :] for k in range(7))
    S = b1 * 1e-1 * e * th ** 3.5 * np.exp(b2 * (1 - th))
    df = b3 * 1e-4 * (p * th ** b4 + b5 * e * th ** b6)
    df = 0.535 * df + np.sqrt(0.217 * df * df + 2.1316e-12 * f0 * f0 / th)
    F = f / f0 * (df / ((f0 - f) ** 2 + df ** 2) + df / ((f0 + f) ** 2 + df ** 2))
    Nw = (S * F).sum(axis=1)
    fv = f[:, 0]
    Nd = 0.0
    if dry_cont:
        d0 = 5.6e-4 * (p + e) * th ** 0.8
        Nd = fv * p * th * th * (6.14e-5 / (d0 * (1 + (fv / d0) ** 2))
                                 + 1.4e-12 * p * th ** 1.5 / (1 + 1.9e-5 * fv ** 1.5))
    if split:
        return 0.1820 * fv * (Nox + Nd), 0.1820 * fv * Nw
    return np.maximum(0.1820 * fv * (Nox + Nw + Nd), 0.0)


def gamma_per_line(fc, T, Ptot, rho):
    """(n_layers, 44 + 35 + 1) specific attenuation (dB/km) at one carrier fc,
    split into each O2 line, each H2O line and the dry continuum. Their sum
    over the last axis equals gamma_p676 (the recommendation's N'' is a sum)."""
    T = np.asarray(T, float)[:, None]
    th = 300.0 / T
    e = np.maximum(np.asarray(rho, float)[:, None], 0.0) * T / 216.7
    p = np.maximum(np.asarray(Ptot, float)[:, None] - e, 0.0)
    f = float(fc)
    f0, a1, a2, a3, a4, a5, a6 = (O2[:, k][None, :] for k in range(7))
    S = a1 * 1e-7 * p * th ** 3 * np.exp(a2 * (1 - th))
    df = a3 * 1e-4 * (p * th ** (0.8 - a4) + 1.1 * e * th)
    df = np.sqrt(df * df + 2.25e-6)
    dl = (a5 + a6 * th) * 1e-4 * (p + e) * th ** 0.8
    F = f / f0 * ((df - dl * (f0 - f)) / ((f0 - f) ** 2 + df ** 2)
                  + (df - dl * (f0 + f)) / ((f0 + f) ** 2 + df ** 2))
    go = 0.1820 * f * S * F
    f0, b1, b2, b3, b4, b5, b6 = (H2O[:, k][None, :] for k in range(7))
    S = b1 * 1e-1 * e * th ** 3.5 * np.exp(b2 * (1 - th))
    df = b3 * 1e-4 * (p * th ** b4 + b5 * e * th ** b6)
    df = 0.535 * df + np.sqrt(0.217 * df * df + 2.1316e-12 * f0 * f0 / th)
    F = f / f0 * (df / ((f0 - f) ** 2 + df ** 2) + df / ((f0 + f) ** 2 + df ** 2))
    gw = 0.1820 * f * S * F
    d0 = 5.6e-4 * (p + e) * th ** 0.8
    Nd = f * p * th * th * (6.14e-5 / (d0 * (1 + (f / d0) ** 2))
                            + 1.4e-12 * p * th ** 1.5 / (1 + 1.9e-5 * f ** 1.5))
    return np.concatenate([go, gw, 0.1820 * f * Nd], axis=1)


def sat_vapor_pressure(T):
    """Saturation vapor pressure over water (hPa), ITU-R P.453 eq. (9)."""
    t = np.asarray(T, float) - 273.15
    return 6.1121 * np.exp((18.678 - t / 234.5) * t / (t + 257.14))


def max_rh(h_mid, T, rho, top_km=10.0):
    """Largest relative humidity (over water) below top_km."""
    m = np.asarray(h_mid) <= top_km
    e = np.asarray(rho)[m] * np.asarray(T)[m] / 216.7
    return float(np.max(e / sat_vapor_pressure(np.asarray(T)[m])))


def sky_temperature(A_db, T_surface):
    """Sky noise temperature (K) seen through a gaseous loss A_db, ITU-R
    P.618-14 Sec. 3: T_sky = T_mr (1 - 10^(-A/10)) + 2.7 * 10^(-A/10), with the
    mean radiating temperature T_mr = 37.34 + 0.81 T_s."""
    t = 10 ** (-np.asarray(A_db, float) / 10)
    return (37.34 + 0.81 * T_surface) * (1 - t) + 2.7 * t


def planck_temperature(T, f_ghz):
    """Noise-equivalent (Planck) temperature (K) of a blackbody at T for the
    available noise power per hertz: (hf/k) / (exp(hf/kT) - 1)."""
    x = 0.0479924 * np.asarray(f_ghz, float)  # hf/k in K for f in GHz
    T = np.asarray(T, float)
    return x / np.expm1(x / T)


def sky_temperature_rt(paths, G, T_layer, f_ghz, T_bg=2.725):
    """Downwelling sky noise temperature (K) at the terminal by radiative
    transfer along the ray (round 3, D25). paths: (n_el, n_layer) path
    lengths a_n (km); G: (n_layer, n_freq) specific attenuation (dB/km);
    T_layer: (n_layer,) layer temperature (K). Each layer emits at its own
    Planck temperature and is attenuated by every layer below it; the
    cosmic background enters through the whole column. Returns
    (T_sky, T_mr) of shape (n_el, n_freq), with T_mr = (T_sky - J(T_bg) t) /
    (1 - t) the mean radiating temperature that P.618 approximates."""
    tau = paths[:, :, None] * G[None, :, :] / (10 * np.log10(np.e))  # nepers
    below = np.cumsum(tau, axis=1) - tau  # optical depth between ground and layer
    J = planck_temperature(T_layer[:, None], np.asarray(f_ghz, float)[None, :])  # (n_layer, n_freq)
    em = (J[None] * (-np.expm1(-tau)) * np.exp(-below)).sum(axis=1)
    t = np.exp(-tau.sum(axis=1))
    Jb = planck_temperature(T_bg, np.asarray(f_ghz, float))[None, :]
    T_sky = em + Jb * t
    with np.errstate(invalid="ignore", divide="ignore"):
        T_mr = np.where(t < 1 - 1e-12, em / (1 - t), J[0][None, :] * np.ones_like(em))
    return T_sky, T_mr


def gamma_grid(f, h_mid, T, P, rho, **kw):
    """(n_layers, n_freq) specific attenuation, dB/km."""
    return np.stack([gamma_p676(f, T[i], P[i], rho[i], **kw) for i in range(len(h_mid))])


# ---------------------------------------------------------------------------
# 3. Slant-path geometry
# ---------------------------------------------------------------------------


def refractive_index(T, P, rho):
    """ITU-R P.453 eqs. (1)-(2) with the dry partial pressure."""
    e = rho * T / 216.7
    pd = P - e
    N = 77.6 * pd / T + 72.0 * e / T + 3.75e5 * e / T ** 2
    return 1.0 + N * 1e-6


def annex1_path_lengths(el_deg, h_lo, delta, n=None):
    """Path length a_n (km) in each layer, P.676-13 Annex 1 eqs. (17)-(19).

    n: refractive index at each layer (None -> straight ray, n = 1).
    Returns (a_n, bending_ok) where layers above the path apex are zero (the
    ray is Earth-space and rising, so every layer is crossed once).
    """
    r = R_E + h_lo
    beta = math.pi / 2 - math.radians(el_deg)  # incidence angle at layer 1
    a = np.zeros_like(h_lo)
    nr = np.ones_like(h_lo) if n is None else n
    b = beta
    for i in range(len(h_lo)):
        ri, di = r[i], delta[i]
        cb = math.cos(b)
        ai = -ri * cb + 0.5 * math.sqrt(4 * ri * ri * cb * cb + 8 * ri * di + 4 * di * di)
        a[i] = ai
        # eq. (18b); eq. (18a) is deprecated in P.676-13 for degraded accuracy
        alpha = math.asin(max(-1.0, min(1.0, ri / (ri + di) * math.sin(b))))
        if i + 1 < len(h_lo):
            s = nr[i] / nr[i + 1] * math.sin(alpha)
            b = math.asin(max(-1.0, min(1.0, s)))
    return a


def geometric_elevation_reached(launch_deg, h_lo, delta, n, h_sat=550.0):
    """Geometric elevation (deg), seen from the terminal, of the point at
    altitude h_sat that a refracted ray launched at launch_deg reaches.
    Layers follow eqs. (17), (18b), (19); above the top layer the ray is straight."""
    r = R_E + h_lo
    b = math.pi / 2 - math.radians(launch_deg)
    nr = np.ones_like(h_lo) if n is None else n
    phi = 0.0
    alpha = b
    for i in range(len(h_lo)):
        ri, di = r[i], delta[i]
        cb = math.cos(b)
        ai = -ri * cb + 0.5 * math.sqrt(4 * ri * ri * cb * cb + 8 * ri * di + 4 * di * di)
        phi += math.asin(max(-1.0, min(1.0, ai * math.sin(b) / (ri + di))))
        alpha = math.asin(max(-1.0, min(1.0, ri / (ri + di) * math.sin(b))))
        if i + 1 < len(h_lo):
            b = math.asin(max(-1.0, min(1.0, nr[i] / nr[i + 1] * math.sin(alpha))))
    r_top = r[-1] + delta[-1]
    r_sat = R_E + h_sat
    phi += alpha - math.asin(r_top * math.sin(alpha) / r_sat)
    r_obs = r[0]  # terminal radius: R_E at sea level, R_E + h_lower when elevated
    return math.degrees(math.atan2(r_sat * math.cos(phi) - r_obs, r_sat * math.sin(phi)))


def apparent_elevation(theta_geo, h_lo, delta, n, h_sat=550.0):
    """Launch (apparent) elevation whose refracted ray reaches a satellite at
    altitude h_sat seen at geometric elevation theta_geo (P.676-13 Annex 1
    Sec. 2.2.3 relates the two; here solved exactly by ray tracing)."""
    from scipy.optimize import brentq
    if theta_geo >= 89.999 or n is None:
        return theta_geo
    g = lambda x: geometric_elevation_reached(x, h_lo, delta, n, h_sat) - theta_geo
    return brentq(g, theta_geo - 0.05, min(90.0, theta_geo + 2.0), xtol=1e-7)


def annex1_paths_geo(theta_geo, h_lo, delta, n, h_sat=550.0):
    """Refracted-ray layer path lengths toward a satellite at geometric elevation theta_geo."""
    return annex1_path_lengths(apparent_elevation(theta_geo, h_lo, delta, n, h_sat), h_lo, delta, n)


def module_slab_weights(el_deg, h_ground=0.0, h_sat=550.0, steps=10):
    """Exact mirror of ThzNtnMolecularAbsorption::ComputeSlantPathAbsorption.

    Returns (h_mid, ds_km) of the midpoint quadrature the C++ model uses.
    """
    bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 100.0]
    th = max(math.radians(el_deg), math.radians(5.0))
    ct = math.cos(th)
    hlo, htop = min(h_ground, h_sat), min(max(h_ground, h_sat), 100.0)
    Rg = R_E + hlo
    hm, ds = [], []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        bot, top = max(lo, hlo), min(hi, htop)
        if bot >= top:
            continue
        dh = (top - bot) / steps
        for k in range(steps):
            m = bot + (k + 0.5) * dh
            arg = 1.0 - (Rg / (R_E + m)) ** 2 * ct * ct
            s = math.sqrt(arg) if arg > 0 else math.sin(math.radians(5.0))
            s = max(s, math.sin(math.radians(1.0)))
            hm.append(m)
            ds.append(dh / s)
    return np.array(hm), np.array(ds)


# ---------------------------------------------------------------------------
# Atmosphere state on a grid
# ---------------------------------------------------------------------------


def state(h_mid, pwv_mm=None, rho0=None, profile="p835-std", scale_km=2.0):
    """(T, P, rho) at heights h_mid for a profile, scaled to pwv_mm if given."""
    if profile == "p835-std":
        T, P = p835_TP(h_mid)
        if rho0 is None:
            rho0 = 7.5 if pwv_mm is None else rho0_for_pwv(pwv_mm, scale_km)
        return T, P, p835_rho(h_mid, rho0, scale_km)
    T, P, rho = p835_seasonal(h_mid, profile)
    if pwv_mm is not None:
        h_lo, d = annex1_layers()
        _, _, r_ref = p835_seasonal(h_lo + 0.5 * d, profile)
        rho = rho * pwv_mm / column_mm(h_lo, d, r_ref)
    return T, P, rho


def profile_for_pwv(pwv_mm):
    """Base atmosphere used for a column: the P.835-6 Sec. 1.1 profile up to
    15 mm (relative humidity at most 80 percent below 10 km), the Sec. 2
    low-latitude profile from 20 to 50 mm (at most 92 percent); scaling either
    further would supersaturate part of the column."""
    return "p835-std" if pwv_mm <= 15.5 else "low"


def atm_state(h_mid, profile, scale):
    """(T, P, rho) of a base profile with water vapor multiplied by scale,
    exactly what the thz-ntn WaterVapourScale attribute does."""
    if profile == "p835-std":
        return state(h_mid, rho0=7.5 * scale)
    T, P, rho = p835_seasonal(h_mid, profile)
    return T, P, rho * scale


def scale_for_pwv(pwv_mm):
    """Scale giving the column on the Annex 1 layers."""
    prof = profile_for_pwv(pwv_mm)
    if prof == "p835-std":
        return rho0_for_pwv(pwv_mm) / 7.5, prof
    h_lo, d = annex1_layers()
    _, _, r = p835_seasonal(h_lo + 0.5 * d, prof)
    return pwv_mm / column_mm(h_lo, d, r), prof


# ---------------------------------------------------------------------------
# 4. Independent codes
# ---------------------------------------------------------------------------


def itur_annex1(f, el, rho0):
    """ITU-Rpy P.676 Annex 1 slant path (its own P.835 and ray tracing)."""
    import itur.models.itu676 as m
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.array([m.gaseous_attenuation_slant_path(fi, el, rho0, 1013.25, 288.15,
                                                          mode="exact").value for fi in np.atleast_1d(f)])


_A2_P1 = _A2_P2 = None


def p676_13_annex2(f, el, rho_s, P_s=1013.25, T_s=288.15, V_t=None):
    """P.676-13 Annex 2, the method in force (Secs. 1.1 and 2.1, or 2.2 with V_t).

    Ps is TOTAL surface pressure (hPa). h_o from data Part 1 (eq. 31), h_w from
    eq. (37) and Table 4, K_V from data Part 2 (eq. 39); gamma_o and gamma_w
    are Annex 1 at the surface (eqs. 30, 36). The coefficient files cover
    1-350 GHz only, so the method is undefined above 350 GHz: NaN is returned
    there. Validated against ITU-R SG3 workbook CG-3M3J-13-ValEx-Rev8.3.0
    (test_annex2_rev13.py).
    """
    global _A2_P1, _A2_P2
    if _A2_P1 is None:
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        _A2_P1 = np.loadtxt(os.path.join(d, "p676_13_annex2_part1.txt"))
        _A2_P1 = _A2_P1[np.argsort(_A2_P1[:, 0])]
        _A2_P2 = np.loadtxt(os.path.join(d, "p676_13_annex2_part2.txt"))
        _A2_P2 = _A2_P2[np.argsort(_A2_P2[:, 0])]
    f = np.atleast_1d(np.asarray(f, float))
    go, gw = gamma_p676(f, T_s, P_s, rho_s, split=True)
    ao, bo, co, do = (np.interp(f, _A2_P1[:, 0], _A2_P1[:, k]) for k in range(1, 5))
    ho = ao + bo * T_s + co * P_s + do * rho_s
    if V_t is None:
        hw = 5.6585e-5 * f + 1.8348 + sum(ai / ((f - fi) ** 2 + bi) for fi, ai, bi in
                                          ((22.235080, 2.6846, 2.7649), (183.310087, 5.8905, 4.9219),
                                           (325.152888, 2.9810, 3.0748)))
        A = (go * ho + gw * hw) / np.sin(np.radians(el))
    else:
        aV, bV, cV, dV = (np.interp(f, _A2_P2[:, 0], _A2_P2[:, k]) for k in range(1, 5))
        KV = aV + bV * rho_s + cV * T_s + dV * P_s
        A = (go * ho + KV * V_t) / np.sin(np.radians(el))
    return np.where(f <= 350.0, A, np.nan)


def itur_annex2(f, el, rho_s, P_s=1013.25, T_s=288.15, V_t=None, h=0.0):
    """ITU-Rpy 0.4.0 P.676 Annex 2. NOTE: this is the superseded P.676-12
    closed-form method (equivalent heights; zenith water vapor if V_t), kept
    only to quantify what tools still built on it report."""
    import itur.models.itu676 as m
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if V_t is None:
            return np.array([m.gaseous_attenuation_slant_path(fi, el, rho_s, P_s, T_s,
                                                              mode="approx").value
                             for fi in np.atleast_1d(f)])
        return np.array([m.gaseous_attenuation_slant_path(fi, el, rho_s, P_s, T_s, V_t=V_t, h=h,
                                                          mode="approx").value
                         for fi in np.atleast_1d(f)])


def gamma_pyrtlib(f, T, P, rho, model, o2_model=None):
    """(n_layers, n_freq) dB/km from PyRTlib: H2O(model) + O2(o2_model) + N2."""
    from pyrtlib.absorption_model import H2OAbsModel, O2AbsModel, N2AbsModel
    from pyrtlib.rt_equation import RTEquation
    from pyrtlib.utils import import_lineshape
    o2m = o2_model or model
    H2OAbsModel.model = model
    O2AbsModel.model = o2m
    N2AbsModel.model = o2m if model != "MWL24" else "MWL24"
    H2OAbsModel.h2oll = import_lineshape("h2oll")
    O2AbsModel.o2ll = import_lineshape("o2ll")
    H2OAbsModel.set_ll()
    O2AbsModel.set_ll()
    e = rho * T / 216.7
    awet, adry = RTEquation.clearsky_absorption(np.asarray(P, float), np.asarray(T, float),
                                                np.asarray(e, float), np.asarray(f, float))
    return (np.asarray(awet) + np.asarray(adry)) * NP_TO_DB  # Np/km -> dB/km
