# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnet
User friendly magnetic field calculations

Top level module for exposing the public API of usable modules
"""

__version__ = "0.5.2"

from . import config, forces, magnets, plots, utils
from .magnets._magnet_base import list, reset
from .magnets._polygon3D import get_total_field_mesh
from .utils._demag import (
    solve_demag_tanh,
    solve_demag_tanh_batch,
    solve_demagnetization,
)
from .utils._routines2D import (
    BdotgradB_2D,
    FgradB_2D,
    get_field_2D,
    gradB_2D,
    grid2D,
    jacobian_B_2D,
)
from .utils._routines3D import (
    BdotgradB_3D,
    FgradB_3D,
    get_field_3D,
    gradB_3D,
    grid3D,
    jacobian_B_3D,
    line3D,
    slice3D,
)
from .utils._vector_structs import Jacobian2, Jacobian3

__all__ = [
    "BdotgradB_2D",
    "BdotgradB_3D",
    "FgradB_2D",
    "FgradB_3D",
    "Jacobian2",
    "Jacobian3",
    "config",
    "forces",
    "get_field_2D",
    "get_field_3D",
    "get_total_field_mesh",
    "gradB_2D",
    "gradB_3D",
    "grid2D",
    "grid3D",
    "jacobian_B_2D",
    "jacobian_B_3D",
    "line3D",
    "list",
    "magnets",
    "plots",
    "reset",
    "slice3D",
    "solve_demag_tanh",
    "solve_demag_tanh_batch",
    "solve_demagnetization",
    "utils",
]
