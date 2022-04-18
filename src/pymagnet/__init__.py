# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnet
User friendly magnetic field calculations

Top level module for exposing the public API of usable modules
"""
__version__ = "0.4.0"

from . import forces, magnets, plots, utils
from .magnets._magnet_base import list, reset
from .utils._routines2D import get_field_2D, grid2D
from .utils._routines3D import get_field_3D

__all__ = [forces, magnets, plots, utils, list, reset, get_field_2D, grid2D, get_field_3D]
