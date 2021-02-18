# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Imports all magnet classes and routines
"""
from . import plot
from . import magnets
from .magnets._magnet import reset_magnets, list_magnets
from .magnets._routines2 import grid2D, B_calc_2D
from .magnets._routines3 import grid3D, B_calc_3D
from math import pi as PI

PI_2 = PI / 2.0
PI_4 = PI / 4.0

u0 = 4E-7 * PI

__all__ = ['magnets', 'plot']
