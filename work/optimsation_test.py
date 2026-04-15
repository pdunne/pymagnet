#!/usr/bin/env python3

import json
import multiprocessing as mp
from itertools import product
from pathlib import Path

import numpy as np
import polars as pl
from tqdm import tqdm

import pymagnet as pm

import matplotlib.pyplot as plt
import seaborn as sns


def gen_alternating_magnets(N=5, width=10, height=10, Jr=1.0):
    """Generate N rectangular magnets with alternating remanence (+/-Jr),
    keeping phi=90.0 for all magnets, centered along x."""

    # Center the array along x
    x0 = -(N - 1) * width / 2
    y0 = -height / 2 - 0.01  # slight offset to avoid perfect grid alignment issues

    for i in range(N):
        center = (x0, y0)

        # Smart alternation of Jr sign
        Jr_i = Jr * (-1) ** i

        _ = pm.magnets.Rectangle(
            width=width, height=height, Jr=Jr_i, center=center, phi=90.0
        )

        x0 += width


def gen_lin_halbach(N=5, width=10, height=10, Jr=1.0):
    """Generate a linear Halbach array of N magnets (N >= 1),
    keeping the same orientation change pattern as the original:
    [LEFT(-Jr,0), DOWN(-Jr,90), RIGHT(+Jr,0), UP(+Jr,90)] repeating."""

    # Starting position and row (same as original code)
    x0 = -(N - 1) * width / 2
    y0 = -height / 2 - 0.01  # slight offset to avoid perfect grid alignment issues

    for i in range(N):
        center = (x0, y0)
        idx = i % 4

        if idx == 0:  # LEFT
            phi = 0.0
            Jr_i = -Jr
        elif idx == 1:  # DOWN
            phi = 90.0
            Jr_i = -Jr
        elif idx == 2:  # RIGHT
            phi = 0.0
            Jr_i = Jr
        else:  # UP
            phi = 90.0
            Jr_i = Jr

        _ = pm.magnets.Rectangle(
            width=width, height=height, Jr=Jr_i, center=center, phi=phi
        )

        x0 += width  # move to the next slot along x


# ------------------------------------------------------------
# 1. Parameter grid generator
# ------------------------------------------------------------
def generate_grid(**params):
    keys = list(params.keys())
    for values in product(*params.values()):
        yield dict(zip(keys, values))


# ------------------------------------------------------------
# 2. Field statistics (Polars internal)
# ------------------------------------------------------------
def compute_stats(field_n):

    # field_n is a 2D array (e.g. shape (N, 801))
    flat = np.asarray(field_n).ravel()  # or .flatten()

    s = pl.Series("field", flat)
    finite = s.filter(s.is_finite())

    return {"min": finite.min(), "mean": finite.mean(), "max": finite.max()}


# ------------------------------------------------------------
# 3. One single simulation run
# ------------------------------------------------------------
def run_simulation(run_input):
    """
    run_input -> (run_id, params_dict)
    This function is run in PARALLEL by multiprocessing.
    """

    run_id, params = run_input
    h = params["height"]
    w = params["width"]
    cfg = params["configuration"]
    N = params["N"]

    # Reset magnets
    pm.reset()

    # Generate magnets
    if cfg == "alternating":
        gen_alternating_magnets(N=N, height=h, width=w, Jr=1.26)
    else:
        gen_lin_halbach(N=N, height=h, width=w, Jr=1.26)

    # Compute field
    points = pm.grid2D(4 * w, 2 * h, unit="mm", ymin=0, num_points=1001)
    field = pm.get_field_2D(points)

    # Stats
    stats = compute_stats(field.n)

    # Save everything to its own folder
    run_dir = Path("results") / f"run_{run_id:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Save metadata
    with open(run_dir / "metadata.json", "w") as f:
        json.dump(params, f, indent=2)

    # Save stats
    pl.DataFrame([stats]).write_parquet(run_dir / "stats.parquet")

    # Save raw field
    pl.DataFrame({"field": field.n.flatten()}).write_parquet(run_dir / "field.parquet")

    return {"run_id": run_id, **params, **stats}


# ------------------------------------------------------------
# 4. Automatic best‑configuration finder
# ------------------------------------------------------------
def find_best(df, key="max"):
    """Return the row where `key` is maximized."""
    return df.sort(key, descending=True).head(1)


# ------------------------------------------------------------
# 5. Pareto‑front analysis
# ------------------------------------------------------------
def pareto_front(df, objectives):
    """
    objectives: list of column names to maximize.
    Returns the Pareto‑optimal subset.
    """
    data = df.to_dicts()
    pareto = []

    for i, a in enumerate(data):
        dominated = False
        for j, b in enumerate(data):
            if i == j:
                continue
            # b dominates a if b >= a on all objectives and > on at least one
            if all(b[o] >= a[o] for o in objectives) and any(
                b[o] > a[o] for o in objectives
            ):
                dominated = True
                break
        if not dominated:
            pareto.append(a)

    return pl.DataFrame(pareto)


def plot_pareto_min_mean(
    df: pl.DataFrame, pf: pl.DataFrame, out="results/pareto_min_mean.pdf"
):
    pdf = df.to_pandas()
    ppf = pf.to_pandas()

    plt.figure(figsize=(7, 6))
    sns.scatterplot(
        data=pdf,
        x="min",
        y="mean",
        hue="configuration",
        alpha=0.35,
        s=60,
        edgecolor=None,
    )

    # Sort the front for a nice connecting line
    ppf_line = ppf.sort_values(["min", "mean"], ascending=[False, False])
    plt.plot(
        ppf_line["min"],
        ppf_line["mean"],
        color="red",
        linewidth=2,
        marker="o",
        label="Pareto front",
    )

    plt.title("Pareto Front (maximize min & mean)")
    plt.xlabel("min (finite-only)")
    plt.ylabel("mean (finite-only)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


# ------------------------------------------------------------
# 6. Main optimisation loop
# ------------------------------------------------------------
if __name__ == "__main__":
    # Build full parameter grid
    param_grid = list(
        generate_grid(
            height=[6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0],
            width=[6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0],
            configuration=["alternating", "linear_halbach"],
            N=[20],  # number of magnets
        )
    )

    # Assign run IDs
    tasks = [(i + 1, params) for i, params in enumerate(param_grid)]

    # Parallel execution
    print(f"Running {len(tasks)} simulations in parallel...")

    with mp.Pool(mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap(run_simulation, tasks), total=len(tasks)))

    # Final results dataframe
    df = pl.DataFrame(results)

    # Save summary
    Path("results").mkdir(exist_ok=True)
    df.write_parquet("results/summary.parquet")
    df.write_csv("results/summary.csv")

    print("\nSummary saved to results/summary.parquet")

    best = find_best(df, key="max")
    print("\nBest configuration (max field):")
    print(best)

    pdf = df.to_pandas()

    for cfg_name, cfg_df in pdf.groupby("configuration"):
        heat = cfg_df.pivot_table(index="height", columns="width", values="max")
        plt.figure(figsize=(8, 6))
        sns.heatmap(heat, annot=True, fmt=".3f", cmap="viridis")
        plt.title(f"Max Field vs Height x Width ({cfg_name})")
        plt.tight_layout()
        plt.savefig(f"results/summary_heatmap_max_{cfg_name}.pdf")
        plt.close()

        heat = cfg_df.pivot_table(index="height", columns="width", values="mean")
        plt.figure(figsize=(8, 6))
        sns.heatmap(heat, annot=True, fmt=".3f", cmap="plasma")
        plt.title(f"Mean Field vs Height x Width ({cfg_name})")
        plt.tight_layout()
        plt.savefig(f"results/summary_heatmap_mean_{cfg_name}.pdf")
        plt.close()

    plt.figure(figsize=(8, 5))
    sns.lineplot(data=pdf, x="width", y="max", hue="configuration", marker="o")
    plt.title("Max Field: Alternating vs Linear Halbach")
    plt.tight_layout()
    plt.savefig("results/config_comparison_max.pdf")
    plt.close()

    plt.figure(figsize=(8, 5))
    sns.lineplot(data=pdf, x="width", y="mean", hue="configuration", marker="o")
    plt.title("Mean Field: Alternating vs Linear Halbach")
    plt.tight_layout()
    plt.savefig("results/config_comparison_mean.pdf")
    plt.close()

    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=pdf, x="min", y="max", hue="configuration", s=100)
    plt.title("Min vs Max Field (Pareto Structure)")
    plt.tight_layout()
    plt.savefig("results/min_vs_max_scatter.pdf")
    plt.close()

    pf = pareto_front(df, objectives=["min", "mean"])
    print(pf)

    pf_pdf = pf.to_pandas()

    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=pdf, x="min", y="max", hue="configuration", alpha=0.3)
    sns.scatterplot(
        data=pf_pdf, x="min", y="max", s=200, color="red", label="Pareto Front"
    )
    plt.title("Pareto Front: Min vs Max")
    plt.legend()
    plt.tight_layout()
    plt.savefig("results/pareto_min_vs_max.pdf")
    plt.close()

    for cfg_name, cfg_df in pdf.groupby("configuration"):
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection="3d")

        H = sorted(cfg_df["height"].unique())
        W = sorted(cfg_df["width"].unique())
        Z = cfg_df.pivot_table(
            index="height", columns="width", values="max"
        ).values

        H_mesh, W_mesh = np.meshgrid(H, W, indexing="ij")

        ax.plot_surface(H_mesh, W_mesh, Z, cmap="viridis")
        ax.set_xlabel("Height")
        ax.set_ylabel("Width")
        ax.set_zlabel("Max Field")
        plt.title(f"3D Surface: Max Field ({cfg_name})")
        plt.tight_layout()
        plt.savefig(f"results/surface_max_field_{cfg_name}.pdf")
        plt.close()

    plot_pareto_min_mean(df, pf)


# ------------------------------------------------------------
