"""Compares serial vs parallel magnetic field calculation for STL mesh magnets.

This example demonstrates the performance improvement from using numba-parallelized
field calculation for mesh magnets. Parallel processing is enabled by default
and provides significant speedups for meshes with many triangles.

Usage:
    python parallel_field_comparison.py

Note:
    Since v0.5.0, parallel=True is the default for Mesh.get_field().
    Use parallel=False to explicitly use the serial method.
"""

import os
import time

import numpy as np

import pymagnet as pm

# Get path to STL files
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
STL_DIR = os.path.join(SCRIPT_DIR, "stl")


def format_time(seconds):
    """Format time in appropriate units."""
    if seconds < 1e-3:
        return f"{seconds * 1e6:.1f} us"
    elif seconds < 1:
        return f"{seconds * 1e3:.1f} ms"
    else:
        return f"{seconds:.2f} s"


def benchmark_field_calculation(mesh, x, y, z, n_runs=5):
    """Benchmark serial vs parallel field calculation.

    Args:
        mesh: Mesh magnet object
        x, y, z: Coordinate arrays
        n_runs: Number of runs for averaging

    Returns:
        dict: Timing results for serial and parallel methods
    """
    # Warmup run for JIT compilation (parallel=True is default)
    _ = mesh.get_field(x, y, z)

    # Benchmark serial
    serial_times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        Bx_s, By_s, Bz_s = mesh.get_field(x, y, z, parallel=False)
        serial_times.append(time.perf_counter() - start)

    # Benchmark parallel (default)
    parallel_times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        Bx_p, By_p, Bz_p = mesh.get_field(x, y, z)  # parallel=True is default
        parallel_times.append(time.perf_counter() - start)

    # Verify results match (calculate relative difference)
    Bn_s = np.sqrt(Bx_s**2 + By_s**2 + Bz_s**2)
    Bn_p = np.sqrt(Bx_p**2 + By_p**2 + Bz_p**2)
    rel_diff = np.abs(Bn_s - Bn_p) / np.maximum(Bn_s, 1e-10)
    mean_rel_diff = np.mean(rel_diff)

    return {
        "serial_mean": np.mean(serial_times),
        "serial_std": np.std(serial_times),
        "parallel_mean": np.mean(parallel_times),
        "parallel_std": np.std(parallel_times),
        "mean_rel_diff": mean_rel_diff,
        "Bx": Bx_p,
        "By": By_p,
        "Bz": Bz_p,
    }


def main():
    print("=" * 60)
    print("SERIAL vs PARALLEL MESH FIELD CALCULATION")
    print("=" * 60)

    # Check for STL files
    stl_files = {
        "cube": os.path.join(STL_DIR, "cube.stl"),
        "star": os.path.join(STL_DIR, "star.stl"),
        "bunny_500": os.path.join(STL_DIR, "Stanford_Bunny_500.stl"),
        "bunny_10000": os.path.join(STL_DIR, "Stanford_Bunny_10000.stl"),
    }

    # Filter to available files
    available = {k: v for k, v in stl_files.items() if os.path.exists(v)}

    if not available:
        print(f"No STL files found in {STL_DIR}")
        return

    print(f"\nFound {len(available)} STL files")
    print("\nWarming up numba JIT compilation...")

    # Warmup with smallest mesh
    first_file = list(available.values())[0]
    pm.reset()
    warmup_mesh = pm.magnets.Mesh(first_file, Jr=1.0)
    _ = warmup_mesh.get_field(0.0, 0.0, 5.0)  # parallel=True is default
    pm.reset()

    print("Warmup complete.\n")

    # Run benchmarks
    results = {}

    for name, filepath in available.items():
        pm.reset()
        mesh = pm.magnets.Mesh(filepath, Jr=1.0)
        n_triangles = len(mesh.mesh_vectors)

        # Create evaluation grid
        grid_size = 20
        x = np.linspace(-20, 20, grid_size)
        y = np.linspace(-20, 20, grid_size)
        z = np.array([15.0])
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        n_points = X.size

        print(f"\n{'-' * 60}")
        print(f"Mesh: {name}")
        print(f"Triangles: {n_triangles}")
        print(f"Evaluation points: {n_points}")
        print(f"{'-' * 60}")

        # Run fewer iterations for large meshes
        n_runs = 3 if n_triangles > 1000 else 5

        result = benchmark_field_calculation(mesh, X, Y, Z, n_runs=n_runs)
        results[name] = result

        speedup = result["serial_mean"] / result["parallel_mean"]

        print(
            f"  Serial:   {format_time(result['serial_mean']):>10} +/- {format_time(result['serial_std'])}"
        )
        print(
            f"  Parallel: {format_time(result['parallel_mean']):>10} +/- {format_time(result['parallel_std'])}"
        )
        print(f"  Speedup:  {speedup:.1f}x")
        print(
            f"  Mean rel diff: {result['mean_rel_diff']:.2e} (due to parallel accumulation order)"
        )

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(
        f"\n{'Mesh':<15} {'Triangles':>10} {'Serial':>12} {'Parallel':>12} {'Speedup':>10}"
    )
    print("-" * 60)

    for name, result in results.items():
        pm.reset()
        mesh = pm.magnets.Mesh(available[name], Jr=1.0)
        n_tri = len(mesh.mesh_vectors)
        speedup = result["serial_mean"] / result["parallel_mean"]
        print(
            f"{name:<15} {n_tri:>10} {format_time(result['serial_mean']):>12} "
            f"{format_time(result['parallel_mean']):>12} {speedup:>9.1f}x"
        )

    print("\n" + "=" * 60)
    print("USAGE")
    print("=" * 60)
    print("""
Parallel processing is now the default for mesh field calculations:

    import pymagnet as pm

    # Create mesh magnet
    mesh = pm.magnets.Mesh("your_mesh.stl", Jr=1.0)

    # Calculate field (parallel=True by default)
    Bx, By, Bz = mesh.get_field(x, y, z)

    # To use serial method explicitly:
    Bx, By, Bz = mesh.get_field(x, y, z, parallel=False)

Parallel processing is optimal for:
- Meshes with > 100 triangles
- Evaluation grids with many points
- Batch calculations where startup overhead is amortized

For small meshes (< 50 triangles), you may want to use parallel=False
to avoid parallel overhead.
""")


if __name__ == "__main__":
    main()
