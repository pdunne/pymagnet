#!/usr/bin/env python3
"""Self-consistent demagnetization solver for 2D circular magnetic elements.

Solves H_int = H_ext - N * M with M = f_MH(H_int) via 1D root-finding.
"""

import warnings
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, fsolve


@dataclass
class DemagResult:
    """Result of the demagnetization solver."""

    M_solution: float
    H_int: float
    converged: bool
    solver_used: str  # "brentq" or "fsolve"


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
    MH_interp: PchipInterpolator,
    M_sat: float,
    N: float = 0.5,
) -> DemagResult:
    """Solve the self-consistent demagnetization equation.

    Finds M such that M = f_MH(H_ext - N * M) by solving
    g(M) = M - f_MH(H_ext - N * M) = 0.

    Args:
        H_ext: External applied field magnitude.
        MH_interp: Prebuilt PchipInterpolator mapping H -> M.
        M_sat: Saturation magnetization (upper bracket bound).
        N: Demagnetizing factor. Default 0.5 (2D circle / infinite cylinder).

    Returns:
        DemagResult with M_solution, H_int, convergence flag, and solver used.
    """
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

        print(f"\n{passed}/{total} tests passed")

    _run_tests()
