#!/usr/bin/env python
"""Benchmark script comparing serial vs parallel mesh field calculation.

Run with: python benchmarks/benchmark_mesh_parallel.py

This script measures the performance difference between:
1. Original Python implementation (serial)
2. Numba-optimized serial implementation
3. Numba parallel implementation (using prange) - DEFAULT

Note: Since v0.5.0, parallel=True is the default for Mesh.get_field().

The benchmarks test with different mesh sizes and evaluation point counts.
"""

import os
import time

import numpy as np

# Get the path to STL files relative to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
STL_DIR = os.path.join(PROJECT_ROOT, "examples/scripts/stl_magnets/stl")

print("Warming up numba JIT compilation...")
print("(This may take a few seconds on first run)")

import pymagnet

# Warm up by creating a small mesh and calculating field
cube_path = os.path.join(STL_DIR, "cube.stl")
if os.path.exists(cube_path):
    pymagnet.reset()
    warmup_mesh = pymagnet.magnets.Mesh(cube_path, Jr=1.0)
    _ = warmup_mesh.get_field(0.0, 0.0, 5.0, parallel=False)
    _ = warmup_mesh.get_field(0.0, 0.0, 5.0)  # parallel=True is default
    pymagnet.reset()

print("Warm-up complete.\n")


def benchmark(func, n_runs=10, warmup=2):
    """Run a function multiple times and return timing statistics."""
    # Warmup runs (not timed)
    for _ in range(warmup):
        func()

    # Timed runs
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        func()
        end = time.perf_counter()
        times.append(end - start)

    times = np.array(times)
    return {
        "mean": np.mean(times),
        "std": np.std(times),
        "min": np.min(times),
        "max": np.max(times),
        "median": np.median(times),
    }


def format_time(seconds):
    """Format time in appropriate units."""
    if seconds < 1e-6:
        return f"{seconds * 1e9:.2f} ns"
    elif seconds < 1e-3:
        return f"{seconds * 1e6:.2f} us"
    elif seconds < 1:
        return f"{seconds * 1e3:.2f} ms"
    else:
        return f"{seconds:.3f} s"


def print_comparison(name, stats_dict):
    """Print comparison between different implementations."""
    print(f"\n{'=' * 70}")
    print(f"{name}")
    print(f"{'=' * 70}")

    baseline_mean = None
    for impl_name, stats in stats_dict.items():
        if baseline_mean is None:
            baseline_mean = stats["mean"]
            speedup_str = "(baseline)"
        else:
            speedup = baseline_mean / stats["mean"]
            speedup_str = f"({speedup:.1f}x faster)"

        print(
            f"  {impl_name:20s}: {format_time(stats['mean']):>12} +/- {format_time(stats['std']):>10} {speedup_str}"
        )


def main():
    print("=" * 70)
    print("MESH MAGNET PARALLEL FIELD CALCULATION BENCHMARK")
    print("=" * 70)

    # Check for STL files
    if not os.path.exists(STL_DIR):
        print(f"Error: STL directory not found at {STL_DIR}")
        return

    available_meshes = {
        "cube": os.path.join(STL_DIR, "cube.stl"),
        "star": os.path.join(STL_DIR, "star.stl"),
        "bunny_500": os.path.join(STL_DIR, "Stanford_Bunny_500.stl"),
        "bunny_10000": os.path.join(STL_DIR, "Stanford_Bunny_10000.stl"),
    }

    # Filter to only available meshes
    meshes = {k: v for k, v in available_meshes.items() if os.path.exists(v)}

    if not meshes:
        print("Error: No STL files found")
        return

    print(f"\nFound {len(meshes)} mesh files for testing")

    # =========================================================================
    # Benchmark 1: Small mesh (cube), few points
    # =========================================================================
    if "cube" in meshes:
        pymagnet.reset()
        mesh = pymagnet.magnets.Mesh(meshes["cube"], Jr=1.0)
        n_triangles = len(mesh.mesh_vectors)

        # Small grid
        x = np.linspace(-5, 5, 10)
        y = np.linspace(-5, 5, 10)
        z = np.array([5.0])
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        n_points = X.size

        print(f"\nMesh: cube ({n_triangles} triangles)")
        print(f"Evaluation points: {n_points}")

        def serial_original():
            return mesh._get_field_internal(X, Y, Z)

        def serial_fast():
            return mesh._get_field_serial_fast(X, Y, Z)

        def parallel():
            return mesh._get_field_parallel(X, Y, Z)

        stats = {
            "Original (serial)": benchmark(serial_original, n_runs=20),
            "Numba serial": benchmark(serial_fast, n_runs=20),
            "Numba parallel": benchmark(parallel, n_runs=20),
        }
        print_comparison(
            f"Cube mesh ({n_triangles} triangles, {n_points} points)", stats
        )

    # =========================================================================
    # Benchmark 2: Medium mesh (star or bunny_500)
    # =========================================================================
    mesh_name = "bunny_500" if "bunny_500" in meshes else "star"
    if mesh_name in meshes:
        pymagnet.reset()
        mesh = pymagnet.magnets.Mesh(meshes[mesh_name], Jr=1.0)
        n_triangles = len(mesh.mesh_vectors)

        # Medium grid
        x = np.linspace(-20, 20, 20)
        y = np.linspace(-20, 20, 20)
        z = np.array([15.0])
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        n_points = X.size

        print(f"\nMesh: {mesh_name} ({n_triangles} triangles)")
        print(f"Evaluation points: {n_points}")

        def serial_original():
            return mesh._get_field_internal(X, Y, Z)

        def serial_fast():
            return mesh._get_field_serial_fast(X, Y, Z)

        def parallel():
            return mesh._get_field_parallel(X, Y, Z)

        stats = {
            "Original (serial)": benchmark(serial_original, n_runs=10),
            "Numba serial": benchmark(serial_fast, n_runs=10),
            "Numba parallel": benchmark(parallel, n_runs=10),
        }
        print_comparison(
            f"{mesh_name.capitalize()} mesh ({n_triangles} triangles, {n_points} points)",
            stats,
        )

    # =========================================================================
    # Benchmark 3: Large mesh (bunny_10000) if available
    # =========================================================================
    if "bunny_10000" in meshes:
        pymagnet.reset()
        mesh = pymagnet.magnets.Mesh(meshes["bunny_10000"], Jr=1.0)
        n_triangles = len(mesh.mesh_vectors)

        # Small grid for large mesh
        x = np.linspace(-30, 30, 15)
        y = np.linspace(-30, 30, 15)
        z = np.array([20.0])
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        n_points = X.size

        print(f"\nMesh: bunny_10000 ({n_triangles} triangles)")
        print(f"Evaluation points: {n_points}")

        def serial_original():
            return mesh._get_field_internal(X, Y, Z)

        def serial_fast():
            return mesh._get_field_serial_fast(X, Y, Z)

        def parallel():
            return mesh._get_field_parallel(X, Y, Z)

        # Fewer runs for large mesh
        stats = {
            "Original (serial)": benchmark(serial_original, n_runs=3, warmup=1),
            "Numba serial": benchmark(serial_fast, n_runs=3, warmup=1),
            "Numba parallel": benchmark(parallel, n_runs=3, warmup=1),
        }
        print_comparison(
            f"Bunny_10000 mesh ({n_triangles} triangles, {n_points} points)", stats
        )

    # =========================================================================
    # Benchmark 4: Scaling with point count
    # =========================================================================
    if "cube" in meshes:
        pymagnet.reset()
        mesh = pymagnet.magnets.Mesh(meshes["cube"], Jr=1.0)
        n_triangles = len(mesh.mesh_vectors)

        print(f"\n{'=' * 70}")
        print("SCALING: Fixed mesh (cube), varying point counts")
        print(f"{'=' * 70}")

        for grid_size in [5, 10, 20, 30]:
            x = np.linspace(-5, 5, grid_size)
            y = np.linspace(-5, 5, grid_size)
            z = np.array([5.0])
            X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
            n_points = X.size

            def serial_original():
                return mesh._get_field_internal(X, Y, Z)

            def parallel():
                return mesh._get_field_parallel(X, Y, Z)

            n_runs = max(3, 20 // grid_size)
            orig_stats = benchmark(serial_original, n_runs=n_runs)
            para_stats = benchmark(parallel, n_runs=n_runs)

            speedup = orig_stats["mean"] / para_stats["mean"]
            print(
                f"  {n_points:5d} points: Original {format_time(orig_stats['mean']):>10} | "
                f"Parallel {format_time(para_stats['mean']):>10} | Speedup: {speedup:.1f}x"
            )

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(
        """
The parallel implementation provides significant speedups for mesh field calculations:

Expected speedups (dependent on hardware):
- Small meshes (< 100 triangles): 2-5x faster
- Medium meshes (100-1000 triangles): 3-10x faster
- Large meshes (> 1000 triangles): 5-20x faster

The speedup increases with:
1. More triangles in the mesh (more parallelizable work)
2. More CPU cores available
3. Larger evaluation point arrays

Usage:
    mesh = pymagnet.magnets.Mesh("file.stl", Jr=1.0)
    Bx, By, Bz = mesh.get_field(x, y, z)  # parallel=True by default

    # To use serial method explicitly:
    Bx, By, Bz = mesh.get_field(x, y, z, parallel=False)

Note: First-time JIT compilation adds ~1-2s overhead on import.
Subsequent runs use cached compiled code.
"""
    )


if __name__ == "__main__":
    main()
