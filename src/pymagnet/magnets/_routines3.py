# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Routines for Three Dimensional Magnet Classes
"""
__all__ = ['B_calc_3D', 'grid3D']

import numpy as _np
from ._fields import Vector3


def grid3D(ux, uy, uz, **kwargs):    
    """Grid of x,y, z points.
    Args:
        
        ux ([float]): [x upper limit]
        uy ([float]): [y upper limit]
        uz ([float]): [z upper limit]
        lx ([float, optional]): [x lower limit defaults to -ux]
        ly ([optional]): [y lower limit defaults to -uy]
        lz ([optional]): [x lower limit defaults to -ux]
        NP (int, optional): [number of points in each direction]. Defaults to 100.

    Returns:
        grid: numpy array
    """
    NP = kwargs.pop('NP', None)

    if NP is None:
        NPx = kwargs.pop('NPx', 100)
        NPy = kwargs.pop('NPy', 100)
        NPz = kwargs.pop('NPz', 100)
    else:
        NPx = NP
        NPy = NP
        NPz = NP

    lx = kwargs.pop('lx', -1*ux)
    ly = kwargs.pop('ly', -1*uy)
    lz = kwargs.pop('lz', -1 * uz)
    return(_np.mgrid[lx:ux:NPx * 1j, ly:uy:NPy * 1j, lz:uz:NPz * 1j])


def B_calc_3D(x, y, z):
    from ._magnet3 import Magnet_3D, Prism, Cube, Cylinder
    """Function to calculate magnetic field due to any array of points
       It sums the magnetic field B over each component of the magnetisation
       J = mu_0 M
    """
    if (_np.isscalar(x)):
        x = _np.atleast_1d(x)
    if (_np.isscalar(y)):
        y = _np.atleast_1d(y)
    if (_np.isscalar(z)):
        z = _np.atleast_1d(z)

    if (_np.ndim(x) == 3):  # volume meshgrid
        B = Vector3(_np.zeros_like(x), _np.zeros_like(x), _np.zeros_like(x))

    elif (_np.ndim(x) == 2):  # planar slice
        B = Vector3(_np.zeros_like(x), _np.zeros_like(x), _np.zeros_like(x))

    elif (_np.ndim(y) == 2):  # planar slice
        B = Vector3(_np.zeros_like(y), _np.zeros_like(y), _np.zeros_like(y))

    else:
        B = Vector3(_np.zeros(max(x.size, y.size, z.size)),
                    _np.zeros(max(x.size, y.size, z.size)),
                    _np.zeros(max(x.size, y.size, z.size)))  # line or single point

    for magnet in Magnet_3D.instances:
        if issubclass(magnet.__class__, Prism):
        
            if _np.fabs(magnet.Jx) > magnet.tol:
                # print('Calc Jx')
                Bx, By, Bz = magnet._calcB_prism_x(x - magnet.xc,
                                        y - magnet.yc,
                                        z - magnet.zc)
                B.x += Bx
                B.y += By
                B.z += Bz

            if _np.fabs(magnet.Jy) > magnet.tol:
                # print('Calc Jy')
                Bx, By, Bz = magnet._calcB_prism_y(x - magnet.xc,
                                        y - magnet.yc,
                                        z - magnet.zc)
                B.x += Bx
                B.y += By
                B.z += Bz

            if _np.fabs(magnet.Jz) > magnet.tol:
                # print('Calc Jz')
                Bx, By, Bz = magnet._calcB_prism_z(x - magnet.xc,
                                        y - magnet.yc,
                                        z - magnet.zc)
                B.x += Bx
                B.y += By
                B.z += Bz
            
        elif issubclass(magnet.__class__, Cylinder):
            print("Cylinder not yet supported")

    B.calc_norm()
    return B


def _rotvec3(x, y, z, theta, phi, xc, yc, zc):
    pass


def _inv_rotvec3(x, y, z, theta, phi, xc, yc, zc):
    pass
