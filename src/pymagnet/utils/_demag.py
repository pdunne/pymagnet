# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Self-consistent demagnetization solver.

Solves H_int = H_ext - N * M with M = f_MH(H_int) via 1D root-finding.
"""

import warnings
from collections.abc import Callable
from dataclasses import dataclass

import numpy as _np
from numba import njit, prange
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, fsolve


@dataclass
class DemagResult:
    """Result of the demagnetization solver."""

    M_solution: float
    H_int: float
    converged: bool
    solver_used: str  # "brentq" or "fsolve"


def tanh_MH_model(Ms: float, chi: float):
    """Return an analytical M(H) function: M = Ms * tanh(chi * H / Ms).

    At low fields this gives M ≈ chi * H (linear regime), and saturates
    smoothly to Ms at high fields.

    Args:
        Ms: Saturation magnetization.
        chi: Dimensionless initial susceptibility (slope dM/dH at H=0).

    Returns:
        Callable mapping H (float or array) -> M.
    """

    def f_MH(H):
        return Ms * _np.tanh(chi * _np.asarray(H, dtype=float) / Ms)

    return f_MH


def build_MH_interpolator(
    H_data: _np.ndarray, M_data: _np.ndarray
) -> PchipInterpolator:
    """Build a monotonicity-preserving interpolator from an MH curve.

    Args:
        H_data: Applied field values, strictly monotonically increasing.
        M_data: Magnetization values, monotonically non-decreasing.

    Returns:
        PchipInterpolator mapping H -> M.
    """
    H_data = _np.asarray(H_data, dtype=float)
    M_data = _np.asarray(M_data, dtype=float)

    if H_data.shape != M_data.shape:
        raise ValueError(
            f"H_data and M_data must have the same shape, "
            f"got {H_data.shape} and {M_data.shape}"
        )

    if H_data.ndim != 1:
        raise ValueError(f"Expected 1D arrays, got ndim={H_data.ndim}")

    dH = _np.diff(H_data)
    if _np.any(dH <= 0):
        raise ValueError("H_data must be strictly monotonically increasing")

    dM = _np.diff(M_data)
    if _np.any(dM < 0):
        raise ValueError("M_data must be monotonically non-decreasing")

    return PchipInterpolator(H_data, M_data)


def solve_demagnetization(
    H_ext: float,
    MH_interp: Callable,
    M_sat: float,
    N: float = 0.5,
) -> DemagResult:
    """Solve the self-consistent demagnetization equation.

    Finds M such that M = f_MH(H_ext - N * M) by solving
    g(M) = M - f_MH(H_ext - N * M) = 0.

    Args:
        H_ext: External applied field magnitude.
        MH_interp: Callable mapping H -> M. Can be a PchipInterpolator from
            build_MH_interpolator or an analytical model from tanh_MH_model.
        M_sat: Saturation magnetization (upper bracket bound).
        N: Demagnetizing factor. Default 0.5 (2D circle / infinite cylinder).

    Returns:
        DemagResult with M_solution, H_int, convergence flag, and solver used.
    """
    if isinstance(MH_interp, PchipInterpolator):
        H_max = MH_interp.x[-1]
        if H_ext > H_max:
            warnings.warn(
                f"H_ext={H_ext} exceeds interpolation domain max={H_max}; "
                f"extrapolation will be used",
                stacklevel=2,
            )

    def residual(M):
        M_val = _np.asarray(M, dtype=float)
        result = M_val - MH_interp(H_ext - N * M_val)
        return float(result) if result.ndim == 0 else result

    g0 = residual(0.0)
    g_sat = residual(M_sat)

    # g0 == 0 means M=0 is the solution; g_sat == 0 means M=M_sat is the solution
    if abs(g0) < 1e-12:
        return DemagResult(
            M_solution=0.0, H_int=H_ext, converged=True, solver_used="exact"
        )
    if abs(g_sat) < 1e-12:
        return DemagResult(
            M_solution=M_sat,
            H_int=H_ext - N * M_sat,
            converged=True,
            solver_used="exact",
        )

    if g0 * g_sat < 0:
        M_solution = brentq(residual, 0.0, M_sat, xtol=1e-6, rtol=1e-6)
        converged = True
        solver_used = "brentq"
    else:
        warnings.warn(
            "Bracket [0, M_sat] invalid for brentq; falling back to fsolve",
            stacklevel=2,
        )
        sol, _info, ier, _msg = fsolve(residual, x0=M_sat / 2, full_output=True)
        M_solution = float(sol[0])
        converged = ier == 1
        solver_used = "fsolve"

    H_int = H_ext - N * M_solution
    return DemagResult(
        M_solution=M_solution,
        H_int=H_int,
        converged=converged,
        solver_used=solver_used,
    )


# ------------------------------------------------------------------
# Numba-accelerated solver (tanh model only)
# ------------------------------------------------------------------
@njit(cache=True)
def _tanh_residual(M, H_ext, N, Ms, chi):
    """g(M) = M - Ms * tanh(chi * (H_ext - N*M) / Ms)"""
    H_int = H_ext - N * M
    return M - Ms * _np.tanh(chi * H_int / Ms)


@njit(cache=True)
def _brentq_njit(H_ext, N, Ms, chi, a, b, xtol=1e-10, maxiter=100):
    """Brent's method root-finder compiled with numba.

    Finds the root of _tanh_residual in [a, b].
    Returns (root, converged) tuple.
    """
    fa = _tanh_residual(a, H_ext, N, Ms, chi)
    fb = _tanh_residual(b, H_ext, N, Ms, chi)

    if abs(fa) < 1e-12:
        return a, True
    if abs(fb) < 1e-12:
        return b, True
    if fa * fb > 0:
        return (a + b) * 0.5, False

    # Ensure |f(b)| <= |f(a)| (b is the better guess)
    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa

    c = a
    fc = fa
    d = b - a
    mflag = True

    for _ in range(maxiter):
        if abs(fb) < xtol:
            return b, True
        if abs(b - a) < xtol:
            return b, True

        # Inverse quadratic interpolation or secant
        if fa != fc and fb != fc:
            s = (
                a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
            )
        else:
            s = b - fb * (b - a) / (fb - fa)

        # Conditions for bisection
        cond1 = not ((3 * a + b) / 4 < s < b or b < s < (3 * a + b) / 4)
        cond2 = mflag and abs(s - b) >= abs(b - c) / 2
        cond3 = (not mflag) and abs(s - b) >= abs(c - d) / 2
        cond4 = mflag and abs(b - c) < xtol
        cond5 = (not mflag) and abs(c - d) < xtol

        if cond1 or cond2 or cond3 or cond4 or cond5:
            s = (a + b) / 2
            mflag = True
        else:
            mflag = False

        fs = _tanh_residual(s, H_ext, N, Ms, chi)
        d = c
        c = b
        fc = fb

        if fa * fs < 0:
            b = s
            fb = fs
        else:
            a = s
            fa = fs

        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa

    return b, True


@njit(cache=True)
def solve_demag_tanh(H_ext, Ms, chi, N=0.5):
    """Solve demagnetization for the tanh model (single scalar H_ext).

    Args:
        H_ext: External applied field magnitude.
        Ms: Saturation magnetization.
        chi: Dimensionless initial susceptibility.
        N: Demagnetizing factor (default 0.5).

    Returns:
        (M_solution, H_int, converged) tuple.
    """
    M_sol, converged = _brentq_njit(H_ext, N, Ms, chi, 0.0, Ms)
    H_int = H_ext - N * M_sol
    return M_sol, H_int, converged


@njit(parallel=True, cache=True)
def solve_demag_tanh_batch(H_ext_array, Ms, chi, N=0.5):
    """Solve demagnetization for an array of H_ext values in parallel.

    Args:
        H_ext_array: 1D array of external field values.
        Ms: Saturation magnetization.
        chi: Dimensionless initial susceptibility.
        N: Demagnetizing factor (default 0.5).

    Returns:
        (M_solutions, H_int_solutions) arrays of same shape as H_ext_array.
    """
    n = H_ext_array.shape[0]
    M_out = _np.empty(n)
    H_int_out = _np.empty(n)

    for i in prange(n):  # ty:ignore[not-iterable]
        M_sol, H_int, _ = solve_demag_tanh(H_ext_array[i], Ms, chi, N)
        M_out[i] = M_sol
        H_int_out[i] = H_int

    return M_out, H_int_out
