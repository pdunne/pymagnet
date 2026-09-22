# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Magnet instantiation from configuration dataclasses."""

import numpy as np

import pymagnet as pm
from pymagnet.magnets import Circle, Cube, Cylinder, Prism, Rectangle, Sphere, Square
from pymagnet.magnets._polygon3D import Mesh

from ._schema import MagnetConfig, SimulationConfig

# Maps type string to (class, dimension_arg_names)
MAGNET_REGISTRY: dict[str, tuple[type, list[str]]] = {
    "Rectangle": (Rectangle, ["width", "height"]),
    "Square": (Square, ["width"]),
    "Circle": (Circle, ["radius"]),
    "Prism": (Prism, ["width", "depth", "height"]),
    "Cube": (Cube, ["width"]),
    "Cylinder": (Cylinder, ["radius", "length"]),
    "Sphere": (Sphere, ["radius"]),
    "Mesh": (Mesh, ["filename"]),
}


def build_magnets(config: SimulationConfig) -> list:
    """Reset the magnet registry and instantiate all magnets from config.

    Args:
        config: Validated SimulationConfig.

    Returns:
        List of instantiated magnet objects.
    """
    pm.reset()
    magnets = []
    for mc in config.magnets:
        magnet = _create_magnet(mc, config.dimension)
        magnets.append(magnet)
    return magnets


def _create_magnet(mc: MagnetConfig, dimension: str):
    """Create a single magnet instance from its MagnetConfig."""
    cls, dim_keys = MAGNET_REGISTRY[mc.type]

    kwargs: dict = {}

    # Dimension arguments
    for key in dim_keys:
        val = getattr(mc, key)
        if val is not None:
            kwargs[key] = val

    # Common arguments
    kwargs["Jr"] = mc.Jr
    kwargs["center"] = np.array(mc.center)
    kwargs["phi"] = mc.phi
    kwargs["alpha"] = mc.alpha

    if dimension == "3D":
        kwargs["theta"] = mc.theta
        kwargs["beta"] = mc.beta
        kwargs["gamma"] = mc.gamma
        kwargs["mask_magnet"] = mc.mask_magnet

    # Mesh-specific
    if mc.type == "Mesh":
        kwargs["mesh_scale"] = mc.mesh_scale

    return cls(**kwargs)
