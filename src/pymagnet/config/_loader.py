# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""TOML configuration file loader."""

import tomllib
from pathlib import Path

from ._schema import (
    ForceConfig,
    GridConfig,
    MagnetConfig,
    PlotConfig,
    SimulationConfig,
    validate_config,
)


def load(toml_path: str | Path) -> SimulationConfig:
    """Load and validate a TOML configuration file.

    Args:
        toml_path: Path to the TOML file.

    Returns:
        SimulationConfig dataclass.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If validation fails.
    """
    toml_path = Path(toml_path)
    if not toml_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {toml_path}")

    with open(toml_path, "rb") as f:
        raw = tomllib.load(f)

    config = _build_config(raw, toml_path.parent)

    errors = validate_config(config)
    if errors:
        msg = "Configuration validation failed:\n" + "\n".join(
            f"  - {e}" for e in errors
        )
        raise ValueError(msg)

    return config


def _build_config(raw: dict, base_dir: Path) -> SimulationConfig:
    """Convert raw TOML dict into a SimulationConfig."""
    meta = raw.get("meta", {})

    config = SimulationConfig(
        name=meta.get("name", ""),
        description=meta.get("description", ""),
        dimension=meta.get("dimension", "3D"),
    )

    # Parse magnets
    for m in raw.get("magnets", []):
        mc = _parse_magnet(m, base_dir)
        config.magnets.append(mc)

    # Parse grids: support both [grid] (single dict) and [[grid]] (list of dicts)
    if "grid" in raw:
        raw_grids = raw["grid"]
        if isinstance(raw_grids, list):
            for g in raw_grids:
                config.grids.append(_parse_grid(g))
        else:
            config.grids.append(_parse_grid(raw_grids))

    # Parse plots: support both [plot] (single dict) and [[plot]] (list of dicts)
    if "plot" in raw:
        raw_plots = raw["plot"]
        if isinstance(raw_plots, list):
            for p in raw_plots:
                config.plots.append(_parse_plot(p))
        else:
            config.plots.append(_parse_plot(raw_plots))

    # If grids defined but no plots, add a default plot per grid
    if config.grids and not config.plots:
        for _ in config.grids:
            config.plots.append(PlotConfig())

    # If plots defined but no grids, add a default grid per plot
    if config.plots and not config.grids:
        for _ in config.plots:
            config.grids.append(GridConfig())

    # Parse force
    if "force" in raw:
        config.force = _parse_force(raw["force"])

    return config


def _parse_magnet(data: dict, base_dir: Path) -> MagnetConfig:
    """Parse a single magnet table into a MagnetConfig."""
    mc = MagnetConfig(
        type=data["type"],
        Jr=data["Jr"],
    )

    # Optional fields with defaults from dataclass
    for key in (
        "center",
        "phi",
        "theta",
        "alpha",
        "beta",
        "gamma",
        "mask_magnet",
        "width",
        "height",
        "depth",
        "radius",
        "length",
        "mesh_scale",
    ):
        if key in data:
            setattr(mc, key, data[key])

    # Resolve mesh filename relative to TOML file location
    if "filename" in data:
        filepath = Path(data["filename"])
        if not filepath.is_absolute():
            filepath = base_dir / filepath
        mc.filename = str(filepath)

    return mc


def _parse_grid(data: dict) -> GridConfig:
    """Parse grid table into a GridConfig."""
    gc = GridConfig()
    for key in (
        "type",
        "num_points",
        "unit",
        "xmax",
        "ymax",
        "zmax",
        "xmin",
        "ymin",
        "zmin",
        "plane",
        "max1",
        "max2",
        "min1",
        "min2",
        "slice_value",
    ):
        if key in data:
            setattr(gc, key, data[key])
    return gc


def _parse_plot(data: dict) -> PlotConfig:
    """Parse plot table into a PlotConfig."""
    pc = PlotConfig()
    for key in (
        "enabled",
        "type",
        "cmap",
        "cmin",
        "cmax",
        "num_levels",
        "num_arrows",
        "show_magnets",
        "vector_plot",
        "opacity",
        "magnet_opacity",
        "cone_opacity",
        "planes",
        "figure_group",
        "save_fig",
        "filename",
    ):
        if key in data:
            setattr(pc, key, data[key])
    return pc


def _parse_force(data: dict) -> ForceConfig:
    """Parse force table into a ForceConfig."""
    fc = ForceConfig()
    for key in ("enabled", "num_samples", "unit", "target_magnet"):
        if key in data:
            setattr(fc, key, data[key])
    return fc
