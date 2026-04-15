# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Runner that orchestrates magnet building, field calculation, and plotting."""

from pathlib import Path
from typing import Any

from ._builder import build_magnets
from ._loader import load
from ._schema import GridConfig, PlotConfig, SimulationConfig


def run(
    toml_path: str | Path, *, output_dir: str | Path | None = None
) -> dict[str, Any]:
    """Load config, create magnets, calculate fields, and generate plots.

    Each [[grid]]/[[plot]] pair produces one entry in the results lists.

    Args:
        toml_path: Path to TOML configuration file.
        output_dir: Output directory for saved figures. If None, figures are
            saved relative to the current working directory.

    Returns:
        dict with keys:
            config: SimulationConfig
            magnets: list of magnet instances
            results: list of dicts, one per grid/plot pair, each containing:
                points, field, figure (if plot enabled), saved_to (if saved)
            force: dict with force/torque (if force enabled)
    """
    config = load(toml_path)

    result: dict[str, Any] = {"config": config}

    # Resolve output directory
    out_path: Path | None = None
    if output_dir is not None:
        out_path = Path(output_dir)

    # Build magnets
    magnets = build_magnets(config)
    result["magnets"] = magnets

    # Process each grid/plot pair, grouping by figure_group where set
    results_list: list[dict[str, Any]] = []

    # Identify figure groups for 3D slice plots
    group_indices: dict[str, list[int]] = {}
    for i, pc in enumerate(config.plots):
        if pc.figure_group:
            group_indices.setdefault(pc.figure_group, []).append(i)

    processed: set[int] = set()

    for i, (gc, pc) in enumerate(zip(config.grids, config.plots)):
        if i in processed:
            continue

        entry: dict[str, Any] = {}

        if config.dimension == "2D":
            points, field = _calculate_field_2d(gc)
            entry["points"] = points
            entry["field"] = field

            if pc.enabled and field is not None:
                fig = _plot_2d(pc, points, field)
                entry["figure"] = fig
                if pc.save_fig:
                    path = _save_figure_2d(fig, pc, i, config.name, out_path)
                    entry["saved_to"] = str(path)
        elif pc.figure_group and len(group_indices.get(pc.figure_group, [])) > 1:
            # Grouped 3D slices: render multiple grid/plot pairs onto one figure
            member_indices = group_indices[pc.figure_group]
            processed.update(member_indices)
            fig = _plot_3d_grouped(config, member_indices)
            entry["figure"] = fig
            # Save using first member's plot config
            first_pc = config.plots[member_indices[0]]
            if first_pc.save_fig:
                path = _save_figure_3d(
                    fig, first_pc, member_indices[0], config.name, out_path
                )
                entry["saved_to"] = str(path)
        else:
            if pc.enabled:
                fig, points, field = _plot_3d(gc, pc)
                entry["figure"] = fig
                entry["points"] = points
                entry["field"] = field
                if pc.save_fig:
                    path = _save_figure_3d(fig, pc, i, config.name, out_path)
                    entry["saved_to"] = str(path)
            else:
                points, field = _calculate_field_3d(gc)
                entry["points"] = points
                entry["field"] = field

        results_list.append(entry)

    result["results"] = results_list

    # Force calculation
    if config.force.enabled:
        result["force"] = _calculate_force(config, magnets)

    return result


def _resolve_save_path(
    pc: PlotConfig,
    plot_index: int,
    config_name: str,
    output_dir: Path | None,
    default_ext: str = ".png",
) -> Path:
    """Determine the file path for saving a figure."""
    if pc.filename:
        filename = pc.filename
    else:
        base = config_name.replace(" ", "_").lower() if config_name else "plot"
        filename = f"{base}_{plot_index}{default_ext}"

    path = Path(filename)
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        path = output_dir / path

    return path


def _save_figure_2d(
    fig,
    pc: PlotConfig,
    plot_index: int,
    config_name: str,
    output_dir: Path | None,
) -> Path:
    """Save a matplotlib figure to disk."""
    path = _resolve_save_path(pc, plot_index, config_name, output_dir)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    return path


def _save_figure_3d(
    fig,
    pc: PlotConfig,
    plot_index: int,
    config_name: str,
    output_dir: Path | None,
) -> Path:
    """Save a plotly figure to disk (png/pdf via kaleido, or html fallback)."""
    path = _resolve_save_path(pc, plot_index, config_name, output_dir)

    if path.suffix.lower() in (".png", ".pdf", ".svg", ".jpeg", ".jpg"):
        try:
            fig.write_image(str(path))
        except (ValueError, ImportError):
            # kaleido not installed — fall back to HTML
            html_path = path.with_suffix(".html")
            fig.write_html(str(html_path))
            return html_path
    elif path.suffix.lower() == ".html":
        fig.write_html(str(path))
    else:
        html_path = path.with_suffix(".html")
        fig.write_html(str(html_path))
        return html_path

    return path


def _calculate_field_2d(gc: GridConfig):
    """Create 2D grid and compute field."""
    import pymagnet as pm

    kwargs: dict[str, Any] = {
        "num_points": gc.num_points,
        "unit": gc.unit,
    }
    if gc.xmin is not None:
        kwargs["xmin"] = gc.xmin
    if gc.ymin is not None:
        kwargs["ymin"] = gc.ymin

    points = pm.grid2D(gc.xmax or 50.0, gc.ymax or 50.0, **kwargs)
    field = pm.get_field_2D(points)
    return points, field


def _calculate_field_3d(gc: GridConfig):
    """Create 3D grid/slice and compute field."""
    import pymagnet as pm
    from pymagnet.utils._routines3D import grid3D, slice3D

    if gc.type == "slice":
        points = slice3D(
            plane=gc.plane,
            max1=gc.max1 or 30.0,
            max2=gc.max2 or 30.0,
            slice_value=gc.slice_value,
            unit=gc.unit,
            num_points=gc.num_points,
        )
    else:
        points = grid3D(
            gc.xmax or 30.0,
            gc.ymax or 30.0,
            gc.zmax or 30.0,
            num_points=gc.num_points,
            unit=gc.unit,
        )

    field = pm.get_field_3D(points)
    return points, field


def _plot_2d(pc: PlotConfig, points, field):
    """Generate 2D contour or stream plot."""
    import pymagnet.plots as plots

    kwargs: dict[str, Any] = {
        "plot_type": pc.type,
        "cmap": pc.cmap,
        "num_levels": pc.num_levels,
        "show_magnets": pc.show_magnets,
    }
    if pc.cmin is not None:
        kwargs["cmin"] = pc.cmin
    if pc.cmax is not None:
        kwargs["cmax"] = pc.cmax
    if pc.num_arrows is not None:
        kwargs["num_arrows"] = pc.num_arrows

    fig, _ax = plots.plot_2D_contour(points, field, **kwargs)
    return fig


def _plot_3d(gc: GridConfig, pc: PlotConfig):
    """Generate 3D slice or volume plot using quickplot functions."""
    import pymagnet.plots as plots

    if pc.type == "slice":
        kwargs: dict[str, Any] = {
            "max1": gc.max1 or 30.0,
            "max2": gc.max2 or 30.0,
            "unit": gc.unit,
            "num_points": gc.num_points,
            "opacity": pc.opacity,
            "magnet_opacity": pc.magnet_opacity,
            "cone_opacity": pc.cone_opacity,
            "planes": pc.planes,
            "show_magnets": pc.show_magnets,
        }
        if gc.min1 is not None:
            kwargs["min1"] = gc.min1
        if gc.min2 is not None:
            kwargs["min2"] = gc.min2
        if pc.cmin is not None:
            kwargs["cmin"] = pc.cmin
        if pc.cmax is not None:
            kwargs["cmax"] = pc.cmax
        if pc.num_levels is not None:
            kwargs["num_levels"] = pc.num_levels
        if pc.num_arrows is not None:
            kwargs["num_arrows"] = pc.num_arrows
        if gc.slice_value != 0.0:
            kwargs["slice_value"] = gc.slice_value

        fig, cache, _data = plots.slice_quickplot(**kwargs)
        points = None
        field = None
        if cache:
            first_plane = next(iter(cache.values()))
            points = first_plane.get("points")
            field = first_plane.get("field")
        return fig, points, field

    elif pc.type == "volume":
        kwargs = {
            "num_points": gc.num_points,
            "unit": gc.unit,
            "xmax": gc.xmax or 30.0,
            "ymax": gc.ymax or 30.0,
            "zmax": gc.zmax or 30.0,
            "show_magnets": pc.show_magnets,
        }
        if pc.cmin is not None:
            kwargs["cmin"] = pc.cmin
        if pc.cmax is not None:
            kwargs["cmax"] = pc.cmax
        if pc.opacity is not None:
            kwargs["opacity"] = pc.opacity
        if pc.magnet_opacity is not None:
            kwargs["magnet_opacity"] = pc.magnet_opacity

        fig, cache, _data = plots.volume_quickplot(**kwargs)
        points = cache.get("points")
        field = cache.get("field")
        return fig, points, field

    return None, None, None


def _plot_3d_grouped(config: SimulationConfig, indices: list[int]):
    """Render multiple grid/plot pairs onto a single plotly 3D figure.

    Each pair contributes its own slice (with independent bounds) to a shared
    figure, using the lower-level plotly helpers.
    """
    import plotly.graph_objects as go

    import pymagnet as pm
    from pymagnet.plots._plotly3D import (
        _draw_cones,
        _draw_surface_slice,
        _generate_all_meshes,
        reset_polyhedra,
    )
    from pymagnet.utils._routines3D import slice3D

    reset_polyhedra()
    data_objects: list = []

    # Use first plot config for shared settings
    first_pc = config.plots[indices[0]]
    if first_pc.show_magnets:
        data_objects.extend(
            _generate_all_meshes(magnet_opacity=first_pc.magnet_opacity)
        )

    colorscale = first_pc.cmap or "viridis"
    unit = config.grids[indices[0]].unit

    for idx in indices:
        gc = config.grids[idx]
        pc = config.plots[idx]

        for plane in pc.planes:
            max1 = gc.max1 or 30.0
            max2 = gc.max2 or 30.0
            min1 = gc.min1 if gc.min1 is not None else -max1
            min2 = gc.min2 if gc.min2 is not None else -max2

            points = slice3D(
                plane=plane,
                max1=max1,
                min1=min1,
                max2=max2,
                min2=min2,
                slice_value=gc.slice_value,
                unit=gc.unit,
                num_points=gc.num_points,
            )
            field = pm.get_field_3D(points)

            cmin = pc.cmin if pc.cmin is not None else 0
            cmax = pc.cmax if pc.cmax is not None else 0.5
            opacity = pc.opacity

            data_objects.append(
                _draw_surface_slice(
                    points,
                    field,
                    colorscale,
                    opacity=opacity,
                    cmin=cmin,
                    cmax=cmax,
                    showscale=True,
                )
            )

            if pc.num_arrows is not None:
                NA = gc.num_points // pc.num_arrows
                if NA < 1:
                    NA = 1
                data_objects.append(
                    _draw_cones(
                        points, field, NA=NA, cone_opacity=pc.cone_opacity
                    )
                )

    fig = go.Figure(data=data_objects)
    fig.update_layout(
        scene=dict(
            xaxis_title=f"x ({unit})",
            yaxis_title=f"y ({unit})",
            zaxis_title=f"z ({unit})",
        ),
        width=700,
        margin=dict(r=20, b=10, l=10, t=10),
    )
    fig.update_layout(scene_aspectmode="data")
    fig.show()
    return fig


def _calculate_force(config: SimulationConfig, magnets: list) -> dict[str, Any]:
    """Calculate force and torque on the target magnet."""
    fc = config.force
    target = magnets[fc.target_magnet]
    force, torque = target.get_force_torque(
        num_samples=fc.num_samples, unit=fc.unit
    )
    return {"force": force, "torque": torque}
