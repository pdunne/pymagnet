# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnet.magnets

This module imports the classes and functions in the private modules to create a public API.

"""
from ._magnet1D import magnetic_field_cylinder_1D, magnetic_field_prism_1D
from ._magnet2D import Circle, Magnet2D, Rectangle, Square
from ._magnet3D import Cube, Cylinder, Magnet3D, Prism, Sphere
from ._magnet_base import Magnet
from ._polygon2D import Line, LineUtils, Polygon, PolyMagnet
from ._polygon3D import Mesh

__all__ = [magnetic_field_prism_1D, magnetic_field_cylinder_1D, Magnet, Magnet2D, Rectangle, Square, Circle, Magnet3D, Prism, Cube, Cylinder, Sphere, Polygon, Line, LineUtils, PolyMagnet, Mesh]
