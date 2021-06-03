# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnet
User friendly magnetic field calculations

Top level module for exposing the public API of usable modules
"""
from . import plots
from . import magnets
from . import utils
from . import forces
from .magnets._magnet_base import reset_magnets, list_magnets
from .utils._routines2D import grid2D, get_field_2D
from .utils._routines3D import grid3D, slice3D, get_field_3D

__all__ = ["magnets", "plots", "utils"]
