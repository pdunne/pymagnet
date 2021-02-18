# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Calculate magnetic fields along symmetry axis
of a cylindrical or prismatic magnet
"""
import numpy as _np

__all__ = ['magnetic_field_prism_1D',
           'magnetic_field_cylinder_1D']

def magnetic_field_prism_1D(magnet, z):
    from ._magnet3 import Magnet_3D
    """Magnetic

    Args:
        magnet_object ([type]): [description]
        Jr ([type]): [description]
        z ([type]): [description]

    Returns:
        [type]: [description]
    """
    if issubclass(magnet.__class__, Magnet_3D):
        a = magnet.a
        b = magnet.b
        c = magnet.c
        Jr = magnet.Jr

        if type(z) is not _np.ndarray:
            z_local = _np.asarray(z) - c - magnet.zc
        else:
            z_local = z - c - magnet.zc

        ab = a * b
        a_sq = _np.power(a, 2)
        b_sq = _np.power(b, 2)
        z_sq = _np.power(z_local, 2)
        zc = z_local + 2 * c
        zc_sq = _np.power(zc, 2)

        data = (_np.arctan2(zc * _np.sqrt(a_sq + b_sq + zc_sq), ab) -
                _np.arctan2(z_local * _np.sqrt(a_sq + b_sq + z_sq), ab))
        data *= (Jr / _np.pi)

        return data
    else:
        print(
            f"Error, the magnet should be a 3D magnet not {magnet.__class__}")
        return _np.full_like(z, _np.nan)


def magnetic_field_cylinder_1D(magnet, z):
    from ._magnet3 import Cylinder
    """Magnetic

    Args:
        magnet_object ([type]): [description]
        Jr ([type]): [description]
        z ([type]): [description]

    Returns:
        [type]: [description]
    """
    if issubclass(magnet.__class__, Cylinder):
        L = magnet.length
        R = magnet.radius
        Jr = magnet.Jr
        if type(z) is not _np.ndarray:
            # z_local = z
            z_local = _np.asarray(z) - magnet.length / 2 - magnet.zc
        else:
            # z_local = z
            z_local = z - magnet.length / 2 - magnet.zc

        zL = z_local + L
        R_sq = _np.power(R, 2)
        z_sq = _np.power(z_local, 2)
        zL_sq = _np.power(zL, 2)

        data = ((zL / _np.sqrt(zL_sq + R_sq)) - 
            (z_local / _np.sqrt(z_sq + R_sq)))
        data *= (Jr / 2)
        return data

    else:
        print(
            f"Error, the magnet should be a 3D magnet not {magnet.__class__}")
        return _np.NaN
