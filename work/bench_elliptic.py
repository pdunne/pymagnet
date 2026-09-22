#!/usr/bin/env python3
"""Benchmark: scipy vs Bulirsch cel vs Carlson cel_carlson."""

import time

import numpy as np
from scipy.special import ellipe, ellipk

from pymagnet.utils._elliptic import cel as cel_bulirsch
from pymagnet.utils._elliptic import cel_carlson


# ---------------------------------------------------------------------------
# Scipy wrappers for K(k) and E(k) via cel parameterisation
# ---------------------------------------------------------------------------
def cel_scipy_K(kc_arr):
    """K(k) via scipy — equivalent to cel(kc, 1, 1, 1)."""
    return ellipk(1.0 - kc_arr**2)


def cel_scipy_E(kc_arr):
    """E(k) via scipy — equivalent to cel(kc, 1, 1, kc²)."""
    return ellipe(1.0 - kc_arr**2)


def bench(label, func, *args, warmup=2, repeats=10):
    """Run func(*args) with warmup, then time `repeats` iterations."""
    for _ in range(warmup):
        func(*args)

    t0 = time.perf_counter()
    for _ in range(repeats):
        func(*args)
    elapsed = (time.perf_counter() - t0) / repeats
    print(f"  {label:30s}  {elapsed*1e3:8.2f} ms")
    return elapsed


def main():
    rng = np.random.default_rng(42)

    # ===================================================================
    # Part 1: Three-way performance — K(k) and E(k)
    # ===================================================================
    print("=" * 70)
    print("Part 1: K(k) and E(k)  —  scipy vs Bulirsch vs Carlson")
    print("=" * 70)
    print("  (scipy uses optimised C; Bulirsch/Carlson use numba @vectorize)")

    for n in [1_000, 10_000, 100_000]:
        kc = rng.uniform(0.01, 0.99, size=n)
        ones = np.ones(n)

        # --- K(k): cel(kc, 1, 1, 1) ---
        print(f"\n  K(k)  —  array size: {n:,}")
        t_sp = bench("scipy ellipk", cel_scipy_K, kc)
        t_bul = bench("Bulirsch cel", cel_bulirsch, kc, ones, ones, ones)
        t_car = bench("Carlson cel_carlson", cel_carlson, kc, ones, ones, ones)
        print(f"  {'Bulirsch / scipy':30s}  {t_sp/t_bul:8.2f}x")
        print(f"  {'Bulirsch / Carlson':30s}  {t_car/t_bul:8.2f}x")

        # --- E(k): cel(kc, 1, 1, kc²) ---
        kc2 = kc**2
        print(f"\n  E(k)  —  array size: {n:,}")
        t_sp = bench("scipy ellipe", cel_scipy_E, kc)
        t_bul = bench("Bulirsch cel", cel_bulirsch, kc, ones, ones, kc2)
        t_car = bench("Carlson cel_carlson", cel_carlson, kc, ones, ones, kc2)
        print(f"  {'Bulirsch / scipy':30s}  {t_sp/t_bul:8.2f}x")
        print(f"  {'Bulirsch / Carlson':30s}  {t_car/t_bul:8.2f}x")

    # ===================================================================
    # Part 2: Bulirsch vs Carlson — general cel(kc, p, c, s)
    # ===================================================================
    print("\n" + "=" * 70)
    print("Part 2: General cel(kc, p, c, s)  —  Bulirsch vs Carlson")
    print("=" * 70)
    print("  (scipy has no general cel equivalent)")

    for n in [1_000, 10_000, 100_000]:
        kc = rng.uniform(0.01, 0.99, size=n)

        # Brho case: cel(kc, 1, 1, -1)
        p1 = np.ones(n)
        c1 = np.ones(n)
        s1 = -np.ones(n)
        print(f"\n  Brho cel(kc, 1, 1, -1)  —  array size: {n:,}")
        t_bul = bench("Bulirsch", cel_bulirsch, kc, p1, c1, s1)
        t_car = bench("Carlson", cel_carlson, kc, p1, c1, s1)
        print(f"  {'Bulirsch / Carlson':30s}  {t_car/t_bul:8.2f}x")

        # Bz case: cel(kc, γ², 1, γ)
        gamma = rng.uniform(-0.9, 0.9, size=n)
        p2 = gamma**2
        c2 = np.ones(n)
        s2 = gamma
        print(f"\n  Bz cel(kc, γ², 1, γ)  —  array size: {n:,}")
        t_bul = bench("Bulirsch", cel_bulirsch, kc, p2, c2, s2)
        t_car = bench("Carlson", cel_carlson, kc, p2, c2, s2)
        print(f"  {'Bulirsch / Carlson':30s}  {t_car/t_bul:8.2f}x")

    # ===================================================================
    # Part 3: Accuracy comparison
    # ===================================================================
    print("\n" + "=" * 70)
    print("Part 3: Accuracy  —  max |error| vs scipy reference")
    print("=" * 70)

    N = 10_000
    kc = rng.uniform(0.01, 0.99, size=N)
    ones = np.ones(N)
    k_sq = 1.0 - kc**2

    ref_K = ellipk(k_sq)
    ref_E = ellipe(k_sq)

    bul_K = cel_bulirsch(kc, ones, ones, ones)
    car_K = cel_carlson(kc, ones, ones, ones)
    bul_E = cel_bulirsch(kc, ones, ones, kc**2)
    car_E = cel_carlson(kc, ones, ones, kc**2)

    print(f"\n  {'':30s}  {'max |err|':>12s}  {'rel err':>12s}")
    print(f"  {'-'*30}  {'-'*12}  {'-'*12}")
    for label, result, ref in [
        ("Bulirsch K(k)", bul_K, ref_K),
        ("Carlson  K(k)", car_K, ref_K),
        ("Bulirsch E(k)", bul_E, ref_E),
        ("Carlson  E(k)", car_E, ref_E),
    ]:
        abs_err = np.max(np.abs(result - ref))
        rel_err = np.max(np.abs((result - ref) / ref))
        print(f"  {label:30s}  {abs_err:12.2e}  {rel_err:12.2e}")


if __name__ == "__main__":
    main()
