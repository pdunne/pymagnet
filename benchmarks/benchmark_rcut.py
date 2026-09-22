#!/usr/bin/env python
"""Benchmark: r_cut accuracy vs speed trade-off for Mesh.get_field().

Evaluates the field of Stanford_Bunny_10000 on a 30×30×30 volume grid
(±70 mm per axis) at a sweep of distance-cutoff values.  For each r_cut the
script records:

  * median wall-clock time (ms)
  * max absolute error  vs the no-cutoff reference  (T)
  * RMS error           vs the no-cutoff reference  (T)
  * nRMSE              — RMS error / RMS of reference field  (dimensionless)
  * pointwise rel. RMS — RMS of (err / |B_ref|) per point   (dimensionless)

A two-panel matplotlib figure is saved to reports/profiling/rcut_benchmark.png
and also shown interactively.

Usage
-----
    uv run python benchmarks/benchmark_rcut.py
"""

import os
import time

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
STL_DIR = os.path.join(PROJECT_ROOT, "examples", "scripts", "stl_magnets", "stl")
REPORT_DIR = os.path.join(PROJECT_ROOT, "reports", "profiling")
os.makedirs(REPORT_DIR, exist_ok=True)

BUNNY_STL = os.path.join(STL_DIR, "Stanford_Bunny_10000.stl")

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------
XMAX = YMAX = ZMAX = 70.0   # mm
N_PER_AXIS = 30

# ---------------------------------------------------------------------------
# r_cut sweep (mm).  np.inf = no cutoff (reference).
# ---------------------------------------------------------------------------
R_CUT_VALUES = [10, 20, 30, 40, 50, 60, 80, 100, 130, 160, 200, np.inf]

# Timed runs per r_cut (+ 1 warm-up run not included in timing)
N_RUNS = 3

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt(seconds):
    if seconds < 1e-3:
        return f"{seconds * 1e6:.1f} µs"
    if seconds < 1:
        return f"{seconds * 1e3:.1f} ms"
    return f"{seconds:.3f} s"


def _median_time(func, n_runs, warmup=1):
    for _ in range(warmup):
        func()
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        func()
        times.append(time.perf_counter() - t0)
    return float(np.median(times))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if not os.path.exists(BUNNY_STL):
        raise FileNotFoundError(
            f"Stanford_Bunny_10000.stl not found at:\n  {BUNNY_STL}"
        )

    print("Importing pymagnet …")
    import pymagnet as pm
    from pymagnet.magnets import Mesh

    # -- JIT warm-up ---------------------------------------------------------
    cube_stl = os.path.join(STL_DIR, "cube.stl")
    if os.path.exists(cube_stl):
        print("Warming up Numba JIT …")
        pm.reset()
        _m = Mesh(cube_stl, Jr=1.0)
        _x = np.linspace(-10, 10, 5)
        _X, _Y, _Z = np.meshgrid(_x, _x, _x, indexing="ij")
        _m.get_field(_X, _Y, _Z)                       # no cutoff
        _m.get_field(_X, _Y, _Z, r_cut=5.0)            # with cutoff
        pm.reset()
        print("Warm-up complete.\n")

    # -- Load bunny ----------------------------------------------------------
    pm.reset()
    print(f"Loading {os.path.basename(BUNNY_STL)} …")
    magnet = Mesh(BUNNY_STL, Jr=1.0)
    n_tri = len(magnet.mesh_vectors)
    print(f"  {n_tri} triangles loaded.")

    # -- Build evaluation grid -----------------------------------------------
    axis = np.linspace(-XMAX, XMAX, N_PER_AXIS)
    X, Y, Z = np.meshgrid(axis, axis, axis, indexing="ij")
    n_pts = X.size
    print(f"\nEvaluation grid: {N_PER_AXIS}³ = {n_pts} points  "
          f"(±{XMAX:.0f} mm per axis)")
    print(f"Timed runs per r_cut: {N_RUNS}  (+1 warm-up)\n")

    # -- Reference field (no cutoff) -----------------------------------------
    print("Computing reference field (r_cut = ∞) …")
    ref_time = _median_time(
        lambda: magnet.get_field(X, Y, Z),
        n_runs=N_RUNS,
    )
    Bx_ref, By_ref, Bz_ref = magnet.get_field(X, Y, Z)
    B_ref_mag = np.sqrt(Bx_ref**2 + By_ref**2 + Bz_ref**2)
    print(f"  Reference median: {_fmt(ref_time)}\n")

    # Pre-compute denominators for relative error metrics
    B_ref_rms = float(np.sqrt(np.mean(B_ref_mag**2)))   # for nRMSE
    # guard against near-zero field (points very close to a node)
    B_ref_safe = np.maximum(B_ref_mag, 1e-12 * B_ref_mag.max())

    # -- Sweep r_cut ---------------------------------------------------------
    header = (
        f"{'r_cut':>10}  {'median':>10}  {'speedup':>8}"
        f"  {'max |err|':>12}  {'RMS err':>12}  {'nRMSE':>10}  {'pw rel RMS':>12}"
    )
    print(header)
    print("-" * len(header))

    results = []   # (r_cut, time_s, max_err, rms_err, nrmse, pw_rel_rms)

    for r_cut in R_CUT_VALUES:
        label = f"{r_cut:.0f} mm" if np.isfinite(r_cut) else "∞  (ref)"

        t_med = _median_time(
            lambda rc=r_cut: magnet.get_field(X, Y, Z, r_cut=rc),
            n_runs=N_RUNS,
        )

        Bx, By, Bz = magnet.get_field(X, Y, Z, r_cut=r_cut)
        err = np.sqrt((Bx - Bx_ref)**2 + (By - By_ref)**2 + (Bz - Bz_ref)**2)
        max_err   = float(err.max())
        rms_err   = float(np.sqrt(np.mean(err**2)))
        nrmse     = rms_err / B_ref_rms
        pw_rel_rms = float(np.sqrt(np.mean((err / B_ref_safe)**2)))

        speedup = ref_time / t_med

        results.append((r_cut, t_med, max_err, rms_err, nrmse, pw_rel_rms))

        print(f"  {label:>8}  {_fmt(t_med):>10}  {speedup:>7.2f}×"
              f"  {max_err:>12.3e}  {rms_err:>12.3e}"
              f"  {nrmse:>10.4f}  {pw_rel_rms:>12.4f}")

    # -- Plot ----------------------------------------------------------------
    finite_mask = [np.isfinite(r[0]) for r in results]
    r_vals      = np.array([r[0] for r, m in zip(results, finite_mask) if m])
    t_vals      = np.array([r[1] * 1e3 for r, m in zip(results, finite_mask) if m])  # ms
    max_errs    = np.array([r[2] for r, m in zip(results, finite_mask) if m])
    rms_errs    = np.array([r[3] for r, m in zip(results, finite_mask) if m])
    nrmses      = np.array([r[4] for r, m in zip(results, finite_mask) if m])
    pw_rel_rmss = np.array([r[5] for r, m in zip(results, finite_mask) if m])

    ref_ms = ref_time * 1e3

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    fig.suptitle(
        f"r_cut trade-off — Stanford Bunny 10 000 triangles\n"
        f"{N_PER_AXIS}³ = {n_pts} eval points  (±{XMAX:.0f} mm)",
        fontsize=11,
    )

    # --- top panel: duration ------------------------------------------------
    ax1.plot(r_vals, t_vals, "o-", color="steelblue", label="with r_cut")
    ax1.axhline(ref_ms, color="steelblue", linestyle="--", alpha=0.6,
                label=f"no cutoff ({ref_ms:.0f} ms)")
    ax1.set_ylabel("Median time (ms)")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # secondary y-axis: speedup
    ax1b = ax1.twinx()
    ax1b.plot(r_vals, ref_ms / t_vals, "s--", color="darkorange", alpha=0.7,
              label="speedup")
    ax1b.axhline(1.0, color="darkorange", linestyle=":", alpha=0.4)
    ax1b.set_ylabel("Speedup  (×)", color="darkorange")
    ax1b.tick_params(axis="y", labelcolor="darkorange")
    ax1b.legend(fontsize=9, loc="upper right")

    # --- middle panel: absolute errors --------------------------------------
    ax2.semilogy(r_vals, max_errs, "o-", color="crimson",  label="max |error| (T)")
    ax2.semilogy(r_vals, rms_errs, "s-", color="darkgreen", label="RMS error (T)")
    ax2.set_ylabel("Absolute error (T)")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, which="both")

    # --- bottom panel: relative errors --------------------------------------
    ax3.semilogy(r_vals, nrmses,      "o-", color="royalblue",
                 label="nRMSE  =  RMS(err) / RMS(|B_ref|)")
    ax3.semilogy(r_vals, pw_rel_rmss, "s-", color="darkorchid",
                 label="pointwise rel. RMS  =  RMS(err / |B_ref|)")
    ax3.set_xlabel("r_cut  (mm)")
    ax3.set_ylabel("Relative error  (dimensionless)")
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    out_path = os.path.join(REPORT_DIR, "rcut_benchmark.png")
    plt.savefig(out_path, dpi=150)
    print(f"\nFigure saved → {out_path}")
    plt.show()


if __name__ == "__main__":
    main()
