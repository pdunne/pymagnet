#!/usr/bin/env python3
"""Benchmark: numba finite difference kernels vs np.gradient."""

import time

import numpy as np

from pymagnet.utils._gradient_kernels import _gradient_2d, _gradient_3d


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


def np_gradient_2d(F, dx, dy):
    return np.gradient(F, dx, dy)


def np_gradient_3d(F, dx, dy, dz):
    return np.gradient(F, dx, dy, dz)


def numba_gradient_2d(F, dx, dy):
    return _gradient_2d(F, dx, dy)


def numba_gradient_3d(F, dx, dy, dz):
    return _gradient_3d(F, dx, dy, dz)


def main():
    rng = np.random.default_rng(42)

    print("=" * 60)
    print("2D Gradient Benchmark")
    print("=" * 60)

    for n in [100, 500, 1000, 2000]:
        F = np.ascontiguousarray(rng.standard_normal((n, n)))
        dx, dy = 0.5, 0.3
        print(f"\n  Grid: {n}x{n} ({F.size:,} points)")
        t_np = bench("np.gradient", np_gradient_2d, F, dx, dy)
        t_nb = bench("numba _gradient_2d", numba_gradient_2d, F, dx, dy)
        print(f"  {'Speedup':30s}  {t_np/t_nb:8.2f}x")

    print("\n" + "=" * 60)
    print("3D Gradient Benchmark")
    print("=" * 60)

    for n in [20, 50, 100, 150]:
        F = np.ascontiguousarray(rng.standard_normal((n, n, n)))
        dx, dy, dz = 0.5, 0.3, 0.8
        print(f"\n  Grid: {n}x{n}x{n} ({F.size:,} points)")
        t_np = bench("np.gradient", np_gradient_3d, F, dx, dy, dz)
        t_nb = bench("numba _gradient_3d", numba_gradient_3d, F, dx, dy, dz)
        print(f"  {'Speedup':30s}  {t_np/t_nb:8.2f}x")


if __name__ == "__main__":
    main()
