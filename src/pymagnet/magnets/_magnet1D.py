# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""1D High symmetry methods

This private module implements magnetic field calculations in z along the
symmetry centre of a cylindrical or cuboidal magnet.

"""
import numpy as _np
from ..utils.global_const import PI
from ..utils._vector_structs import Field1

__all__ = ["magnetic_field_prism_1D", "magnetic_field_cylinder_1D"]


def magnetic_field_prism_1D(magnet, z):
    """Calculates the magnetic field z-component due to a cuboid along its
    axial symmetry center.

    Args:
        magnet (Magnet3D): magnet object, one of Prism, or Cube
        z (ndarray): Array of points along the symmetry axis

    Returns:
        Field1: z-component of the magnetic field and associated unit ('T')
    """
    from ._magnet3D import Prism

    if issubclass(magnet.__class__, Prism):
        a = magnet.a
        b = magnet.b
        c = magnet.c
        Jr = magnet.Jr

        z_local = _np.asarray(z) - c - magnet.center[2]

        ab = a * b
        a_sq = _np.power(a, 2)
        b_sq = _np.power(b, 2)
        z_sq = _np.power(z_local, 2)
        zc = z_local + 2 * c
        zc_sq = _np.power(zc, 2)

        Bz = _np.arctan2(zc * _np.sqrt(a_sq + b_sq + zc_sq), ab) - _np.arctan2(
            z_local * _np.sqrt(a_sq + b_sq + z_sq), ab
        )
        Bz *= Jr / PI
        field = Field1(Bz)
        return field
    else:
        print(f"Error, the magnet should be a 3D magnet not {magnet.__class__}")
        return None


def magnetic_field_cylinder_1D(magnet, z):
    """Calculates the magnetic field z-component due to a cuboid along its
    axial symmetry center.

    Args:
        magnet (Magnet3D): magnet object, Cylinder
        z (ndarray): Array of points along the symmetry axis

    Returns:
        Field1: z-component of the magnetic field and associated unit ('T')
    """
    from ._magnet3D import Cylinder

    if issubclass(magnet.__class__, Cylinder):
        L = magnet.length
        R = magnet.radius
        Jr = magnet.Jr

        z_local = _np.asarray(z) - magnet.length / 2 - magnet.center[2]

        zL = z_local + L
        R_sq = _np.power(R, 2)
        z_sq = _np.power(z_local, 2)
        zL_sq = _np.power(zL, 2)

        Bz = (zL / _np.sqrt(zL_sq + R_sq)) - (z_local / _np.sqrt(z_sq + R_sq))
        Bz *= Jr / 2
        data = Field1(Bz)
        return data

    else:
        print(f"Error, the magnet should be a 3D magnet not {magnet.__class__}")
        return None
