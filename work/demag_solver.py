#!/usr/bin/env python3
"""Self-consistent demagnetization solver for 2D circular magnetic elements.

Solves H_int = H_ext - N * M with M = f_MH(H_int) via 1D root-finding.
"""

import warnings
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
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
        return Ms * np.tanh(chi * np.asarray(H, dtype=float) / Ms)

    return f_MH


def build_MH_interpolator(
    H_data: np.ndarray, M_data: np.ndarray
) -> PchipInterpolator:
    """Build a monotonicity-preserving interpolator from an MH curve.

    Args:
        H_data: Applied field values, strictly monotonically increasing.
        M_data: Magnetization values, monotonically non-decreasing.

    Returns:
        PchipInterpolator mapping H -> M.
    """
    H_data = np.asarray(H_data, dtype=float)
    M_data = np.asarray(M_data, dtype=float)

    if H_data.shape != M_data.shape:
        raise ValueError(
            f"H_data and M_data must have the same shape, "
            f"got {H_data.shape} and {M_data.shape}"
        )

    if H_data.ndim != 1:
        raise ValueError(f"Expected 1D arrays, got ndim={H_data.ndim}")

    dH = np.diff(H_data)
    if np.any(dH <= 0):
        raise ValueError("H_data must be strictly monotonically increasing")

    dM = np.diff(M_data)
    if np.any(dM < 0):
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
        M_val = np.asarray(M, dtype=float)
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
        sol, info, ier, msg = fsolve(
            residual, x0=M_sat / 2, full_output=True
        )
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
    return M - Ms * np.tanh(chi * H_int / Ms)


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
            s = (a * fb * fc / ((fa - fb) * (fa - fc))
                 + b * fa * fc / ((fb - fa) * (fb - fc))
                 + c * fa * fb / ((fc - fa) * (fc - fb)))
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
    M_out = np.empty(n)
    H_int_out = np.empty(n)

    for i in prange(n):
        M_sol, H_int, _ = solve_demag_tanh(H_ext_array[i], Ms, chi, N)
        M_out[i] = M_sol
        H_int_out[i] = H_int

    return M_out, H_int_out


# ------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------
if __name__ == "__main__":

    def _run_tests():
        tol = 1e-4
        passed = 0
        total = 0

        # --- Test 1: Linear material ---
        total += 1
        chi = 10.0
        H = np.linspace(0, 1e6, 500)
        M = chi * H
        interp = build_MH_interpolator(H, M)
        M_sat = M[-1]

        H_ext = 1e5
        N = 0.5
        result = solve_demagnetization(H_ext, interp, M_sat, N=N)
        M_analytical = chi * H_ext / (1 + N * chi)

        if abs(result.M_solution - M_analytical) < tol * M_analytical:
            print(f"Test 1 (linear material):    PASS  "
                  f"M={result.M_solution:.2f}, expected={M_analytical:.2f}")
            passed += 1
        else:
            print(f"Test 1 (linear material):    FAIL  "
                  f"M={result.M_solution:.2f}, expected={M_analytical:.2f}")

        # --- Test 2: Saturated input (nonlinear curve with saturation) ---
        total += 1
        M_s = 1e6
        H_nl = np.linspace(0, 1e6, 500)
        M_nl = M_s * np.tanh(H_nl / 1e4)
        interp_nl = build_MH_interpolator(H_nl, M_nl)

        H_ext_sat = 1e6
        result_sat = solve_demagnetization(H_ext_sat, interp_nl, M_s, N=N)

        if abs(result_sat.M_solution - M_s) / M_s < 0.01:
            print(f"Test 2 (saturated input):    PASS  "
                  f"M={result_sat.M_solution:.2f}, M_sat={M_s:.2f}")
            passed += 1
        else:
            print(f"Test 2 (saturated input):    FAIL  "
                  f"M={result_sat.M_solution:.2f}, M_sat={M_s:.2f}")

        # --- Test 3: Zero field ---
        total += 1
        result_zero = solve_demagnetization(0.0, interp, M_sat, N=N)

        if abs(result_zero.M_solution) < 1e-6:
            print(f"Test 3 (zero field):         PASS  "
                  f"M={result_zero.M_solution:.6f}")
            passed += 1
        else:
            print(f"Test 3 (zero field):         FAIL  "
                  f"M={result_zero.M_solution:.6f}")

        # --- Test 4: N generalisation ---
        total += 1
        all_ok = True
        for N_val, label in [(1 / 3, "sphere"), (1.0, "thin film")]:
            result_N = solve_demagnetization(H_ext, interp, M_sat, N=N_val)
            M_expected = chi * H_ext / (1 + N_val * chi)
            H_int_expected = H_ext - N_val * M_expected

            if abs(result_N.M_solution - M_expected) > tol * M_expected:
                print(f"Test 4 ({label} N={N_val}):  FAIL  "
                      f"M={result_N.M_solution:.2f}, expected={M_expected:.2f}")
                all_ok = False
            if abs(result_N.H_int - H_int_expected) > tol * abs(H_int_expected):
                print(f"Test 4 ({label} N={N_val}):  FAIL  "
                      f"H_int={result_N.H_int:.2f}, expected={H_int_expected:.2f}")
                all_ok = False

        if all_ok:
            print(f"Test 4 (N generalisation):   PASS")
            passed += 1

        # --- Test 5: tanh_MH_model analytical function ---
        total += 1
        Ms_t = 1e6
        chi_t = 10.0
        f_mh = tanh_MH_model(Ms=Ms_t, chi=chi_t)

        H_ext_t = 1e5
        N_t = 0.5
        result_t = solve_demagnetization(H_ext_t, f_mh, Ms_t, N=N_t)

        # At low field the tanh model is approximately linear: M ≈ chi*H,
        # so the analytical solution M = chi*H_ext/(1+N*chi) should be close.
        M_approx = chi_t * H_ext_t / (1 + N_t * chi_t)

        if result_t.converged and abs(result_t.M_solution - M_approx) / M_approx < 0.05:
            print(f"Test 5 (tanh_MH_model):      PASS  "
                  f"M={result_t.M_solution:.2f}, linear_approx={M_approx:.2f}")
            passed += 1
        else:
            print(f"Test 5 (tanh_MH_model):      FAIL  "
                  f"M={result_t.M_solution:.2f}, linear_approx={M_approx:.2f}, "
                  f"converged={result_t.converged}")

        # --- Test 6: numba scalar solver matches scipy ---
        total += 1
        M_nb, _, conv_nb = solve_demag_tanh(H_ext_t, Ms_t, chi_t, N=N_t)

        if conv_nb and abs(M_nb - result_t.M_solution) / result_t.M_solution < 1e-6:
            print(f"Test 6 (numba scalar):       PASS  "
                  f"M={M_nb:.2f}, scipy={result_t.M_solution:.2f}")
            passed += 1
        else:
            print(f"Test 6 (numba scalar):       FAIL  "
                  f"M={M_nb:.2f}, scipy={result_t.M_solution:.2f}")

        # --- Test 7: numba batch solver ---
        total += 1
        H_batch = np.linspace(0, 5e5, 1000)
        M_batch, _ = solve_demag_tanh_batch(H_batch, Ms_t, chi_t, N=N_t)

        # Verify endpoints
        batch_ok = True
        if abs(M_batch[0]) > 1e-6:
            print(f"Test 7 (numba batch):        FAIL  M[0]={M_batch[0]:.6f} != 0")
            batch_ok = False
        if not np.all(np.diff(M_batch) >= 0):
            print("Test 7 (numba batch):        FAIL  M not monotonic")
            batch_ok = False
        # Check against scalar solver at a few points
        for idx in [0, 250, 500, 750, 999]:
            M_ref, _, _ = solve_demag_tanh(H_batch[idx], Ms_t, chi_t, N=N_t)
            if abs(M_batch[idx] - M_ref) > 1e-4:
                print(f"Test 7 (numba batch):        FAIL  "
                      f"idx={idx}, batch={M_batch[idx]:.4f}, ref={M_ref:.4f}")
                batch_ok = False
                break
        if batch_ok:
            print("Test 7 (numba batch):        PASS  1000 points, all consistent")
            passed += 1

        print(f"\n{passed}/{total} tests passed")

        # --- Benchmark: scipy vs numba ---
        import time

        n_bench = 10_000
        H_bench = np.linspace(0, 5e5, n_bench)
        f_mh_bench = tanh_MH_model(Ms=Ms_t, chi=chi_t)

        # Warmup numba (first call compiles)
        solve_demag_tanh_batch(H_bench[:10], Ms_t, chi_t, N=N_t)

        t0 = time.perf_counter()
        for h in H_bench:
            solve_demagnetization(h, f_mh_bench, Ms_t, N=N_t)
        t_scipy = time.perf_counter() - t0

        t0 = time.perf_counter()
        solve_demag_tanh_batch(H_bench, Ms_t, chi_t, N=N_t)
        t_numba = time.perf_counter() - t0

        print(f"\nBenchmark ({n_bench} solves):")
        print(f"  scipy:  {t_scipy:.4f}s")
        print(f"  numba:  {t_numba:.4f}s")
        print(f"  speedup: {t_scipy / t_numba:.1f}x")

    _run_tests()
