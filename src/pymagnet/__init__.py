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
from .magnets._magnet import reset_magnets, list_magnets
from .magnets._routines2 import grid2D, B_calc_2D
from .magnets._routines3 import grid3D, B_calc_3D
from math import pi as PI

"""float: Module level PI.
"""

PI_2 = PI / 2.0
"""float: Module level PI/2.
"""
PI_4 = PI / 4.0
"""float: Module level PI/4.
"""

u0 = 4e-7 * PI
"""float: Module level u0, permittivity of free space.
"""

__all__ = ["magnets", "plots"]
