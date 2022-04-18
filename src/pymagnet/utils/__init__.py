# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnets.utils

This module imports the classes and functions in the private modules to create a public API, including:

- Quaternion()
- Global Constants
- Point structures
- Vector structures

"""

from ._point_structs import Point2, Point3
from ._quaternion import Quaternion
from ._routines2D import get_field_2D, grid2D, rotate_points_2D
from ._routines3D import get_field_3D, grid3D, line3D, slice3D
from ._trigonometry3D import altitude, norm_plane, rotate_points, signed_area
from ._vector_structs import Field1, Field2, Field3, Point_Array2, Point_Array3
from .global_const import ALIGN_CUTOFF, FP_CUTOFF, MAG_TOL, MU0, PI, PI_2, PI_4

__all__ = [Point2, Point3, Quaternion, grid2D, get_field_2D, rotate_points_2D, get_field_3D, grid3D, slice3D, line3D, signed_area, norm_plane, rotate_points, altitude, Field1, Field2, Field3, Point_Array2, Point_Array3, PI, PI_2, PI_4, MU0, FP_CUTOFF, ALIGN_CUTOFF, MAG_TOL]

# Not implemented yet:
# from ._assemblies2D import *
# from ._fit import *
