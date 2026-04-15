#!/usr/bin/env python
"""Profiling script for pymagnet.

Runs the main computational hot-paths and reports wall-clock timing.
Optionally generates cProfile dumps or line_profiler output.

Usage:
  python benchmarks/profile_codebase.py                  # wall-clock summary table
  python benchmarks/profile_codebase.py --cprofile       # save .prof files (open with snakeviz)
  python benchmarks/profile_codebase.py --line           # line_profiler on hot functions
  python benchmarks/profile_codebase.py --snakeviz       # --cprofile + open snakeviz afterwards

Environment:
  NUMBA_DISABLE_JIT=1  Disables JIT so line_profiler can instrument @njit functions.
                        Use this together with --line for line-level timing inside
                        Numba kernels (runs slower, but fully visible).

  OMP_NUM_THREADS=N    Control thread count for the parallel mesh path.

Scenarios:
  1. analytic_cube       - Cube.get_field() over a 20×20×20 grid
  2. analytic_cylinder   - Cylinder.get_field() over a 20×20×20 grid
  3. mesh_small          - Mesh.get_field() on a cube STL (~12 triangles, 10×10 grid)
  4. mesh_large          - Mesh.get_field() on bunny_500 STL (~500 triangles, 20×20 grid)
  5. prism_force         - calc_force_prism() between two offset cubes
  6. mesh_force          - Mesh.get_force_torque() between cube STL and a Prism
"""

import argparse
import cProfile
import io
import os
import pstats
import statistics
import subprocess
import sys
import time

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
STL_DIR = os.path.join(PROJECT_ROOT, "examples/scripts/stl_magnets/stl")
PROFILES_DIR = os.path.join(SCRIPT_DIR, "profiles")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _format_time(seconds):
    if seconds < 1e-6:
        return f"{seconds * 1e9:.1f} ns"
    elif seconds < 1e-3:
        return f"{seconds * 1e6:.1f} µs"
    elif seconds < 1:
        return f"{seconds * 1e3:.1f} ms"
    else:
        return f"{seconds:.3f} s"


def _time_scenario(func, n_runs=5):
    """Run *func* n_runs times, return median wall-clock time (seconds)."""
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        func()
        times.append(time.perf_counter() - t0)
    return statistics.median(times)


def _stl(name):
    return os.path.join(STL_DIR, name)


# ---------------------------------------------------------------------------
# Scenario factories
# ---------------------------------------------------------------------------

def _make_scenarios(pm):
    """Build and return a list of (name, description, callable) tuples.

    Each callable is safe to call multiple times (state is reset via pm.reset()).
    """
    # --- shared grids ---
    x20 = np.linspace(-20, 20, 20)
    y20 = np.linspace(-20, 20, 20)
    z20 = np.linspace(-20, 20, 20)
    X3d, Y3d, Z3d = np.meshgrid(x20, y20, z20, indexing="ij")

    x10 = np.linspace(-10, 10, 10)
    y10 = np.linspace(-10, 10, 10)
    Z_plane = np.array([15.0])
    X2d, Y2d, Z2d_plane = np.meshgrid(x10, y10, Z_plane, indexing="ij")

    x_mesh_lg = np.linspace(-30, 30, 20)
    y_mesh_lg = np.linspace(-30, 30, 20)
    Z_plane_lg = np.array([20.0])
    X_lg, Y_lg, Z_lg = np.meshgrid(x_mesh_lg, y_mesh_lg, Z_plane_lg, indexing="ij")

    scenarios = []

    # ------------------------------------------------------------------
    # 1. Analytic Cube
    # ------------------------------------------------------------------
    def _build_cube_scenario():
        pm.reset()
        cube = pm.magnets.Cube(Jr=1.0, width=10.0)
        return cube

    cube = _build_cube_scenario()

    def scenario_analytic_cube():
        cube.get_field(X3d, Y3d, Z3d)

    scenarios.append((
        "analytic_cube",
        f"Cube.get_field() — 20×20×20 grid ({X3d.size} pts)",
        scenario_analytic_cube,
    ))

    # ------------------------------------------------------------------
    # 2. Analytic Cylinder
    # ------------------------------------------------------------------
    pm.reset()
    cylinder = pm.magnets.Cylinder(Jr=1.0, radius=5.0, length=10.0)

    def scenario_analytic_cylinder():
        cylinder.get_field(X3d, Y3d, Z3d)

    scenarios.append((
        "analytic_cylinder",
        f"Cylinder.get_field() — 20×20×20 grid ({X3d.size} pts)",
        scenario_analytic_cylinder,
    ))

    # ------------------------------------------------------------------
    # 3. Mesh — small (cube STL)
    # ------------------------------------------------------------------
    cube_stl = _stl("cube.stl")
    if os.path.exists(cube_stl):
        pm.reset()
        mesh_small = pm.magnets.Mesh(cube_stl, Jr=1.0)
        n_tri_small = len(mesh_small.mesh_vectors)

        def scenario_mesh_small():
            mesh_small.get_field(X2d, Y2d, Z2d_plane)

        scenarios.append((
            "mesh_small",
            f"Mesh.get_field() — cube STL ({n_tri_small} tri), 10×10 grid ({X2d.size} pts)",
            scenario_mesh_small,
        ))
    else:
        print(f"  [skip] mesh_small — {cube_stl} not found")

    # ------------------------------------------------------------------
    # 4. Mesh — large (bunny_500 or star as fallback)
    # ------------------------------------------------------------------
    for mesh_name, label in [("Stanford_Bunny_500.stl", "bunny_500"), ("star.stl", "star")]:
        large_stl = _stl(mesh_name)
        if os.path.exists(large_stl):
            pm.reset()
            mesh_large = pm.magnets.Mesh(large_stl, Jr=1.0)
            n_tri_large = len(mesh_large.mesh_vectors)

            def scenario_mesh_large(m=mesh_large):
                m.get_field(X_lg, Y_lg, Z_lg)

            scenarios.append((
                "mesh_large",
                f"Mesh.get_field() — {label} ({n_tri_large} tri), 20×20 grid ({X_lg.size} pts)",
                scenario_mesh_large,
            ))
            break
    else:
        print("  [skip] mesh_large — no suitable STL found")

    # ------------------------------------------------------------------
    # 5. Prism force (two cubes offset along z)
    # ------------------------------------------------------------------
    pm.reset()
    pm.magnets.Cube(Jr=1.0, width=10.0, center=[0.0, 0.0, 0.0])
    cube_force_target = pm.magnets.Cube(Jr=1.0, width=10.0, center=[0.0, 0.0, 15.0])

    def scenario_prism_force():
        pm.forces.calc_force_prism(cube_force_target, num_samples=10)

    scenarios.append((
        "prism_force",
        "calc_force_prism() — two 10mm cubes, 10 samples/face",
        scenario_prism_force,
    ))

    # ------------------------------------------------------------------
    # 6. Mesh force (cube STL + Prism source)
    # ------------------------------------------------------------------
    if os.path.exists(cube_stl):
        pm.reset()
        pm.magnets.Prism(width=10.0, depth=10.0, height=10.0, Jr=1.0, center=[0.0, 0.0, 0.0])
        mesh_force_target = pm.magnets.Mesh(cube_stl, Jr=1.0, center=[0.0, 0.0, 15.0])

        def scenario_mesh_force(m=mesh_force_target):
            m.get_force_torque(depth=2)

        scenarios.append((
            "mesh_force",
            "Mesh.get_force_torque() — cube STL vs Prism source, depth=2",
            scenario_mesh_force,
        ))

    return scenarios


# ---------------------------------------------------------------------------
# Warm-up
# ---------------------------------------------------------------------------

def _warmup(pm):
    """Trigger Numba JIT compilation before any timed run."""
    print("Warming up Numba JIT (first run only)...")
    cube_stl = _stl("cube.stl")
    if os.path.exists(cube_stl):
        pm.reset()
        m = pm.magnets.Mesh(cube_stl, Jr=1.0)
        x0, y0, z0 = np.array([0.0]), np.array([0.0]), np.array([5.0])
        m.get_field(x0, y0, z0, parallel=False)
        m.get_field(x0, y0, z0)   # parallel=True
    else:
        # Fallback: Cube field to warm up analytic kernels
        pm.reset()
        c = pm.magnets.Cube(Jr=1.0, width=10.0)
        c.get_field(np.array([0.0]), np.array([0.0]), np.array([5.0]))
    pm.reset()
    print("Warm-up complete.\n")


# ---------------------------------------------------------------------------
# Wall-clock summary
# ---------------------------------------------------------------------------

def run_summary(pm, scenarios, n_runs=5):
    col_name = 40
    col_time = 12
    col_desc = 60
    header = f"{'Scenario':<{col_name}} {'Median time':>{col_time}}  Description"
    print("=" * (col_name + col_time + col_desc + 4))
    print("WALL-CLOCK SUMMARY")
    print(f"  (median of {n_runs} runs, Numba JIT already warmed up)")
    print("=" * (col_name + col_time + col_desc + 4))
    print(header)
    print("-" * (col_name + col_time + col_desc + 4))

    for name, desc, func in scenarios:
        t = _time_scenario(func, n_runs=n_runs)
        print(f"  {name:<{col_name - 2}} {_format_time(t):>{col_time}}  {desc}")

    print("=" * (col_name + col_time + col_desc + 4))


# ---------------------------------------------------------------------------
# cProfile
# ---------------------------------------------------------------------------

def run_cprofile(pm, scenarios, top_n=20):
    os.makedirs(PROFILES_DIR, exist_ok=True)
    print(f"\nRunning cProfile — saving .prof files to: {PROFILES_DIR}")
    print(f"(Showing top {top_n} cumulative-time entries per scenario)\n")

    for name, desc, func in scenarios:
        prof_path = os.path.join(PROFILES_DIR, f"{name}.prof")
        pr = cProfile.Profile()
        pr.enable()
        func()
        pr.disable()
        pr.dump_stats(prof_path)

        # Pretty-print inline summary
        stream = io.StringIO()
        ps = pstats.Stats(pr, stream=stream).sort_stats("cumulative")
        ps.print_stats(top_n)
        summary = stream.getvalue()

        print(f"{'─' * 70}")
        print(f"Scenario: {name}  ({desc})")
        print(f"Profile saved: {prof_path}")
        print(summary)

    print(f"\nTo explore interactively:")
    print(f"  pip install snakeviz")
    print(f"  snakeviz {PROFILES_DIR}/<scenario>.prof")


# ---------------------------------------------------------------------------
# line_profiler
# ---------------------------------------------------------------------------

def run_line_profiler(pm, scenarios):
    """Instrument known hot functions with line_profiler and print results."""
    try:
        from line_profiler import LineProfiler
    except ImportError:
        print("line_profiler not installed — run: pip install line_profiler")
        sys.exit(1)

    # Import function objects for instrumentation.
    # Note: @njit functions (_get_field_parallel_njit etc.) are opaque to
    # line_profiler unless NUMBA_DISABLE_JIT=1 is set — instrument the
    # Python-level wrappers instead, which dispatch into them.
    from pymagnet.utils._elliptic import cel
    from pymagnet.magnets._polygon3D import Mesh
    from pymagnet.magnets._magnet3D import Prism, Cube, Cylinder
    from pymagnet.forces._prism_force import calc_force_prism

    hotspots = [
        # (label, function_object)
        ("Mesh._get_field_parallel", Mesh._get_field_parallel),
        ("Mesh._get_field_serial_fast", Mesh._get_field_serial_fast),
        ("Mesh._get_field_internal", Mesh._get_field_internal),
        ("Prism._get_field_internal", Prism._get_field_internal),
        ("Cylinder._get_field_internal", Cylinder._get_field_internal),
        ("calc_force_prism", calc_force_prism),
        ("cel", cel),
    ]

    # Also try to add njit wrappers — only meaningful with NUMBA_DISABLE_JIT=1
    jit_disabled = os.environ.get("NUMBA_DISABLE_JIT", "0") == "1"
    if not jit_disabled:
        print(
            "Note: NUMBA_DISABLE_JIT is not set.\n"
            "  @njit functions will appear as single-line opaque calls.\n"
            "  For line-level breakdown inside Numba kernels, re-run with:\n"
            "    NUMBA_DISABLE_JIT=1 python benchmarks/profile_codebase.py --line\n"
        )

    lp = LineProfiler()
    for label, fn in hotspots:
        try:
            lp.add_function(fn)
        except Exception:
            print(f"  [skip] cannot instrument {label} (not a pure-Python function)")

    print("Running line_profiler on all scenarios...\n")
    for name, desc, func in scenarios:
        print(f"  Profiling: {name}")
        lp.enable_by_count()
        func()
        lp.disable_by_count()

    stream = io.StringIO()
    lp.print_stats(stream=stream)
    print(stream.getvalue())


# ---------------------------------------------------------------------------
# snakeviz launcher
# ---------------------------------------------------------------------------

def launch_snakeviz():
    os.makedirs(PROFILES_DIR, exist_ok=True)
    profs = [f for f in os.listdir(PROFILES_DIR) if f.endswith(".prof")]
    if not profs:
        print("No .prof files found — run with --cprofile first.")
        return
    first = os.path.join(PROFILES_DIR, sorted(profs)[0])
    print(f"Opening snakeviz for: {first}")
    subprocess.run(["snakeviz", first])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description="Profile pymagnet computational hot-paths.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--cprofile",
        action="store_true",
        help="Run cProfile and save .prof files under benchmarks/profiles/",
    )
    parser.add_argument(
        "--line",
        action="store_true",
        help="Run line_profiler on known hotspot functions",
    )
    parser.add_argument(
        "--snakeviz",
        action="store_true",
        help="Run cProfile then open the first .prof file in snakeviz",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=5,
        help="Number of timed repetitions for the wall-clock summary (default: 5)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Top-N entries to print in cProfile output (default: 20)",
    )
    return parser.parse_args()


def main():
    args = _parse_args()

    # Import pymagnet after arg parsing so --help is fast
    print("Importing pymagnet...")
    import pymagnet as pm

    _warmup(pm)

    print("Building scenarios...")
    scenarios = _make_scenarios(pm)
    print(f"  {len(scenarios)} scenarios ready.\n")

    # Always print wall-clock summary
    run_summary(pm, scenarios, n_runs=args.runs)

    if args.cprofile or args.snakeviz:
        run_cprofile(pm, scenarios, top_n=args.top)

    if args.line:
        run_line_profiler(pm, scenarios)

    if args.snakeviz:
        launch_snakeviz()


if __name__ == "__main__":
    main()
