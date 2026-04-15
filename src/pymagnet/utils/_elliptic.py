# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2026 Peter Dunne
"""Elliptic integral implementations: Bulirsch cel and Carlson symmetric forms.

Provides two complete elliptic integral implementations:

- ``cel``: Bulirsch's AGM iteration (production default) — fast single-pass
  algorithm with errtol=1e-12 for double-precision accuracy.
- ``cel_carlson``: Carlson symmetric decomposition via RF, RD, RJ — available
  for cases where the individual symmetric integrals are needed directly.

The Carlson primitives (RF, RD, RJ, RC) are also exported for direct use.

References:
    - DLMF §19.2: https://dlmf.nist.gov/19.2 (Bulirsch)
    - DLMF §19.36: https://dlmf.nist.gov/19.36 (Carlson)
    - Carlson, B. C. (1995). Numerical computation of real or complex
      elliptic integrals. *Numerical Algorithms*, 10, 13–26.
"""

from math import fabs, nan, sqrt

from numba import float64, njit, vectorize


@njit(cache=True)
def _elliprf(x, y, z, tol=1e-12):
    """Carlson RF(x, y, z) via the duplication algorithm (DLMF 19.36.i).

    Computes the symmetric elliptic integral of the first kind:
        RF(x,y,z) = (1/2) ∫_0^∞ dt / sqrt((t+x)(t+y)(t+z))

    Args:
        x, y, z: Non-negative real arguments (at most one may be zero).
        tol: Convergence tolerance.

    Returns:
        float: Value of RF(x, y, z).
    """
    while True:
        xl = sqrt(x)
        yl = sqrt(y)
        zl = sqrt(z)
        lam = xl * yl + yl * zl + zl * xl
        x = 0.25 * (x + lam)
        y = 0.25 * (y + lam)
        z = 0.25 * (z + lam)
        A = (x + y + z) / 3.0
        dx = 1.0 - x / A
        dy = 1.0 - y / A
        dz = 1.0 - z / A
        eps = max(abs(dx), abs(dy), abs(dz))
        if eps < tol:
            break

    # Taylor series truncated at 7th order
    e2 = dx * dy - dz * dz
    e3 = dx * dy * dz
    return (1.0 + e2 * (-1.0 / 10.0 + e2 / 14.0 + e3 / 24.0)
            + e3 / 6.0) / sqrt(A)


@njit(cache=True)
def _elliprd(x, y, z, tol=1e-12):
    """Carlson RD(x, y, z) via the duplication algorithm (DLMF 19.36.ii).

    Computes the symmetric elliptic integral of the second kind:
        RD(x,y,z) = (3/2) ∫_0^∞ dt / ((t+z)*sqrt((t+x)(t+y)(t+z)))

    Note: RD(x, y, z) = RJ(x, y, z, z), but this dedicated routine is
    faster because it avoids the extra complexity of the RJ algorithm.

    Args:
        x, y: Non-negative real arguments.
        z: Positive real argument.
        tol: Convergence tolerance.

    Returns:
        float: Value of RD(x, y, z).
    """
    fac = 1.0
    sigma = 0.0

    while True:
        xl = sqrt(x)
        yl = sqrt(y)
        zl = sqrt(z)
        lam = xl * yl + yl * zl + zl * xl
        sigma += fac / (zl * (z + lam))
        fac *= 0.25
        x = 0.25 * (x + lam)
        y = 0.25 * (y + lam)
        z = 0.25 * (z + lam)
        A = (x + y + 3.0 * z) / 5.0
        dx = 1.0 - x / A
        dy = 1.0 - y / A
        dz = 1.0 - z / A
        eps = max(abs(dx), abs(dy), abs(dz))
        if eps < tol:
            break

    # Taylor series
    e2 = dx * dy
    e3 = e2 * dz
    e4 = dz * dz
    e5 = e2 * e4
    return (3.0 * sigma
            + fac * (1.0
                     - 3.0 * e2 / 14.0
                     + e3 / 6.0
                     + 9.0 * e4 / 22.0
                     - 3.0 * e5 / 26.0)
            / (A * sqrt(A)))


@njit(cache=True)
def _elliprj(x, y, z, p, tol=1e-12):
    """Carlson RJ(x, y, z, p) via the duplication algorithm (DLMF 19.36.ii).

    Computes the symmetric elliptic integral of the third kind:
        RJ(x,y,z,p) = (3/2) ∫_0^∞ dt / ((t+p)*sqrt((t+x)(t+y)(t+z)))

    Args:
        x, y, z: Non-negative real arguments (at most one may be zero).
        p: Positive real parameter.
        tol: Convergence tolerance.

    Returns:
        float: Value of RJ(x, y, z, p).
    """
    fac = 1.0
    sigma = 0.0

    while True:
        xl = sqrt(x)
        yl = sqrt(y)
        zl = sqrt(z)
        lam = xl * yl + yl * zl + zl * xl

        # RC accumulation (Carlson 1995, eq. 4.4)
        d = p * (xl + yl + zl) + xl * yl * zl
        alpha = d * d
        beta = p * (p + lam) ** 2
        sigma += fac * _elliprc(alpha, beta)

        fac *= 0.25
        x = 0.25 * (x + lam)
        y = 0.25 * (y + lam)
        z = 0.25 * (z + lam)
        p = 0.25 * (p + lam)
        A = (x + y + z + 2.0 * p) / 5.0
        dx = 1.0 - x / A
        dy = 1.0 - y / A
        dz = 1.0 - z / A
        dp = 1.0 - p / A
        eps = max(abs(dx), abs(dy), abs(dz), abs(dp))
        if eps < tol:
            break

    # Elementary symmetric polynomials of {dx, dy, dz, dp, dp}
    # (dp appears twice because A = (x+y+z+2p)/5)
    # Using dx+dy+dz = -2*dp (since e1 = 0 by construction)
    xyz = dx * dy + dy * dz + dz * dx  # sum of pairwise products of x,y,z devs
    E2 = xyz - 3.0 * dp * dp
    E3 = dx * dy * dz + 2.0 * dp * E2 + 4.0 * dp * dp * dp
    E4 = (2.0 * dx * dy * dz + dp * (E2 + 3.0 * dp * dp)) * dp
    E5 = dp * dp * dx * dy * dz

    return (3.0 * sigma
            + fac * (1.0
                     - 3.0 * E2 / 14.0
                     + E3 / 6.0
                     + 9.0 * E2 * E2 / 88.0
                     - 3.0 * E4 / 22.0
                     - 9.0 * E2 * E3 / 52.0
                     + 3.0 * E5 / 26.0)
            / (A * sqrt(A)))


@njit(cache=True)
def _elliprc(x, y, tol=1e-12):
    """Carlson RC(x, y) — degenerate case of RF.

    RC(x, y) = RF(x, y, y) = (1/2) ∫_0^∞ dt / ((t+y)*sqrt(t+x))

    Args:
        x: Non-negative real argument.
        y: Positive real argument.
        tol: Convergence tolerance.

    Returns:
        float: Value of RC(x, y).
    """
    while True:
        xl = sqrt(x)
        yl = sqrt(y)
        lam = 2.0 * xl * yl + y
        x = 0.25 * (x + lam)
        y = 0.25 * (y + lam)
        A = (x + 2.0 * y) / 3.0
        s = (y - A) / A
        if abs(s) < tol:
            break

    # Taylor series
    return (1.0 + s * s * (3.0 / 10.0 + s * (1.0 / 7.0
            + s * (3.0 / 8.0 + s * 9.0 / 22.0)))) / sqrt(A)


@vectorize([float64(float64, float64, float64, float64)], target="parallel")
def cel(kc, p, c, s):
    """Bulirsch's complete elliptic integral cel(kc, p, c, s).

    Production implementation using the AGM iteration with errtol=1e-12
    for full double-precision accuracy. This is faster than the Carlson
    decomposition while achieving comparable accuracy.

    See DLMF §19.2, https://dlmf.nist.gov/19.2

    Args:
        kc: Complementary modulus (kc = k').
        p: Real parameter.
        c: Real parameter.
        s: Real parameter.

    Returns:
        float: Value of the complete elliptic integral.
    """
    if kc == 0.0:
        return nan

    errtol = 1e-12
    k = fabs(kc)
    pp = p
    cc = c
    ss = s
    em = 1.0

    if p > 0:
        pp = sqrt(p)
        ss = s / pp
    else:
        f = kc * kc
        q = 1.0 - f
        g = 1.0 - pp
        f = f - pp
        q = q * (ss - c * pp)
        pp = sqrt(f / g)
        cc = (c - ss) / g
        ss = -q / (g * g * pp) + cc * pp

    f = cc
    cc = cc + ss / pp
    g = k / pp
    ss = 2.0 * (ss + f * g)
    pp = g + pp
    g = em
    em = k + em
    kk = k

    while fabs(g - k) > g * errtol:
        k = 2.0 * sqrt(kk)
        kk = k * em
        f = cc
        cc = cc + ss / pp
        g = kk / pp
        ss = 2.0 * (ss + f * g)
        pp = g + pp
        g = em
        em = k + em

    return (3.14159265358979323846 / 2.0) * (ss + cc * em) / (em * (em + pp))


@vectorize([float64(float64, float64, float64, float64)], target="parallel")
def cel_carlson(kc, p, c, s):
    """Bulirsch cel(kc, p, c, s) via Carlson symmetric elliptic integrals.

    Alternative implementation using the Carlson decomposition:
        cel(kc, p, c, s) = c * RF(0, kc², 1) + (s - c*p)/3 * RJ(0, kc², 1, p)

    When p = 1, uses the cheaper RD instead of RJ:
        cel(kc, 1, c, s) = c * RF(0, kc², 1) + (s - c)/3 * RD(0, kc², 1)

    This is slower than ``cel`` but provides the Carlson decomposition
    if individual symmetric integrals are needed.

    Args:
        kc: Complementary modulus (kc = k').
        p: Real parameter.
        c: Real parameter.
        s: Real parameter.

    Returns:
        float: Value of the complete elliptic integral.
    """
    if kc == 0.0:
        return nan

    kc2 = kc * kc
    rf = _elliprf(0.0, kc2, 1.0)

    if p == 1.0:
        return c * rf + (s - c) / 3.0 * _elliprd(0.0, kc2, 1.0)
    else:
        return c * rf + (s - c * p) / 3.0 * _elliprj(0.0, kc2, 1.0, p)
