# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnet
User friendly magnetic field calculations

Top level module for exposing the public API of usable modules
"""
from . import forces, magnets, plots, utils
from .magnets._magnet_base import list_magnets, reset_magnets
from .utils._routines2D import get_field_2D, grid2D
from .utils._routines3D import *

__all__ = ["magnets", "plots", "utils"]
