# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Configuration schema dataclasses and validation."""

from dataclasses import dataclass, field

VALID_MAGNET_TYPES = {
    "Rectangle",
    "Square",
    "Circle",
    "Prism",
    "Cube",
    "Cylinder",
    "Sphere",
    "Mesh",
}

REQUIRED_DIMENSIONS: dict[str, list[str]] = {
    "Rectangle": ["width", "height"],
    "Square": ["width"],
    "Circle": ["radius"],
    "Prism": ["width", "depth", "height"],
    "Cube": ["width"],
    "Cylinder": ["radius", "length"],
    "Sphere": ["radius"],
    "Mesh": ["filename"],
}

DIMENSION_2D_TYPES = {"Rectangle", "Square", "Circle"}
DIMENSION_3D_TYPES = {"Prism", "Cube", "Cylinder", "Sphere", "Mesh"}

VALID_GRID_TYPES = {"grid2D", "slice", "grid3D", "volume"}
VALID_PLOT_TYPES_2D = {"contour", "streamplot"}
VALID_PLOT_TYPES_3D = {"slice", "volume"}
VALID_PLANES = {"xy", "xz", "yz"}


@dataclass
class MagnetConfig:
    type: str
    Jr: float
    center: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    phi: float = 90.0
    theta: float = 0.0
    alpha: float = 0.0
    beta: float = 0.0
    gamma: float = 0.0
    mask_magnet: bool = False
    # Dimension fields (type-dependent)
    width: float | None = None
    height: float | None = None
    depth: float | None = None
    radius: float | None = None
    length: float | None = None
    # Mesh-specific
    filename: str | None = None
    mesh_scale: float = 1.0


@dataclass
class GridConfig:
    type: str = "grid2D"
    num_points: int = 100
    unit: str = "mm"
    xmax: float | None = None
    ymax: float | None = None
    zmax: float | None = None
    xmin: float | None = None
    ymin: float | None = None
    zmin: float | None = None
    plane: str = "xy"
    max1: float | None = None
    max2: float | None = None
    min1: float | None = None
    min2: float | None = None
    slice_value: float = 0.0


@dataclass
class PlotConfig:
    enabled: bool = True
    type: str = "contour"
    cmap: str = "viridis"
    cmin: float | None = None
    cmax: float | None = None
    num_levels: int = 11
    num_arrows: int | None = None
    show_magnets: bool = True
    vector_plot: bool = True
    # 3D-specific
    opacity: float = 0.8
    magnet_opacity: float = 1.0
    cone_opacity: float = 1.0
    planes: list[str] = field(default_factory=lambda: ["xy", "xz", "yz"])
    # Figure grouping: plots with the same non-empty figure_group render together
    figure_group: str = ""
    # Save to file
    save_fig: bool = False
    filename: str = ""


@dataclass
class ForceConfig:
    enabled: bool = False
    num_samples: int = 20
    unit: str = "mm"
    target_magnet: int = 0


@dataclass
class SimulationConfig:
    name: str = ""
    description: str = ""
    dimension: str = "3D"
    magnets: list[MagnetConfig] = field(default_factory=list)
    grids: list[GridConfig] = field(default_factory=list)
    plots: list[PlotConfig] = field(default_factory=list)
    force: ForceConfig = field(default_factory=ForceConfig)


def validate_config(config: SimulationConfig) -> list[str]:
    """Validate a SimulationConfig and return a list of error messages.

    Returns:
        list[str]: Empty list means valid configuration.
    """
    errors: list[str] = []

    # Validate dimension
    if config.dimension not in ("2D", "3D"):
        errors.append(
            f"Invalid dimension '{config.dimension}'. Must be '2D' or '3D'."
        )

    # Validate magnets
    if not config.magnets:
        errors.append("At least one magnet must be defined in [[magnets]].")

    for i, mc in enumerate(config.magnets):
        prefix = f"magnets[{i}]"

        if mc.type not in VALID_MAGNET_TYPES:
            errors.append(
                f"{prefix}: Invalid type '{mc.type}'. "
                f"Must be one of: {', '.join(sorted(VALID_MAGNET_TYPES))}"
            )
            continue

        # Check dimension compatibility
        if config.dimension == "2D" and mc.type in DIMENSION_3D_TYPES:
            errors.append(
                f"{prefix}: Type '{mc.type}' is a 3D magnet "
                f"but dimension is '2D'."
            )
        elif config.dimension == "3D" and mc.type in DIMENSION_2D_TYPES:
            errors.append(
                f"{prefix}: Type '{mc.type}' is a 2D magnet "
                f"but dimension is '3D'."
            )

        # Check required dimensions
        for dim_key in REQUIRED_DIMENSIONS.get(mc.type, []):
            val = getattr(mc, dim_key)
            if val is None:
                errors.append(
                    f"{prefix}: Missing required dimension '{dim_key}' "
                    f"for type '{mc.type}'."
                )

        # Check center length
        expected_len = 2 if config.dimension == "2D" else 3
        if len(mc.center) != expected_len:
            errors.append(
                f"{prefix}: center must have {expected_len} elements "
                f"for {config.dimension}, got {len(mc.center)}."
            )

    # Validate grids and plots (paired by index)
    if not config.grids:
        errors.append("At least one [[grid]] must be defined.")

    if config.plots and len(config.plots) != len(config.grids):
        errors.append(
            f"Number of [[plot]] entries ({len(config.plots)}) must match "
            f"number of [[grid]] entries ({len(config.grids)})."
        )

    for i, gc in enumerate(config.grids):
        prefix = f"grid[{i}]"
        if gc.type not in VALID_GRID_TYPES:
            errors.append(
                f"{prefix}: type '{gc.type}' invalid. "
                f"Must be one of: {', '.join(sorted(VALID_GRID_TYPES))}"
            )

        if config.dimension == "2D" and gc.type not in ("grid2D",):
            errors.append(
                f"{prefix}: type '{gc.type}' is not compatible with "
                f"dimension '2D'. Use 'grid2D'."
            )
        if config.dimension == "3D" and gc.type == "grid2D":
            errors.append(
                f"{prefix}: type 'grid2D' is not compatible with "
                f"dimension '3D'. Use 'slice', 'grid3D', or 'volume'."
            )

    for i, pc in enumerate(config.plots):
        prefix = f"plot[{i}]"
        if pc.enabled:
            if config.dimension == "2D" and pc.type not in VALID_PLOT_TYPES_2D:
                errors.append(
                    f"{prefix}: type '{pc.type}' invalid for 2D. "
                    f"Must be one of: {', '.join(sorted(VALID_PLOT_TYPES_2D))}"
                )
            if config.dimension == "3D" and pc.type not in VALID_PLOT_TYPES_3D:
                errors.append(
                    f"{prefix}: type '{pc.type}' invalid for 3D. "
                    f"Must be one of: {', '.join(sorted(VALID_PLOT_TYPES_3D))}"
                )
            for plane in pc.planes:
                if plane not in VALID_PLANES:
                    errors.append(
                        f"{prefix}: planes contains invalid plane '{plane}'. "
                        f"Must be one of: {', '.join(sorted(VALID_PLANES))}"
                    )

    # Validate force
    fc = config.force
    if fc.enabled:
        if not config.magnets:
            errors.append("force.enabled but no magnets defined.")
        elif fc.target_magnet < 0 or fc.target_magnet >= len(config.magnets):
            errors.append(
                f"force.target_magnet={fc.target_magnet} is out of range. "
                f"Must be 0..{len(config.magnets) - 1}."
            )

    return errors
