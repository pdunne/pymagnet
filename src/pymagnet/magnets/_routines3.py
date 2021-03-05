# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Routines for Three Dimensional Magnet Classes
"""
__all__ = ["B_calc_3D", "grid3D"]

import numpy as _np
from ._fields import Vector3

# TODO:
# - Default grid size in calc_3D for x, y, z slices


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
    NP = kwargs.pop("NP", None)

    if NP is None:
        NPx = kwargs.pop("NPx", 100)
        NPy = kwargs.pop("NPy", 100)
        NPz = kwargs.pop("NPz", 100)
    else:
        NPx = NP
        NPy = NP
        NPz = NP

    lx = kwargs.pop("lx", -1 * ux)
    ly = kwargs.pop("ly", -1 * uy)
    lz = kwargs.pop("lz", -1 * uz)
    return _np.mgrid[lx : ux : NPx * 1j, ly : uy : NPy * 1j, lz : uz : NPz * 1j]


def B_calc_3D(x, y, z):
    from ._magnet3 import Magnet_3D, Sphere

    """Function to calculate magnetic field due to any array of points
       It sums the magnetic field B over each component of the magnetisation
       J = mu_0 M
    """
    B = _allocate_field_array3(x, y, z)

    for magnet in Magnet_3D.instances:
        Bx, By, Bz = magnet.calcB(x, y, z)
        B.x += Bx.reshape(B.x.shape)
        B.y += By.reshape(B.y.shape)
        B.z += Bz.reshape(B.z.shape)

    B.calc_norm()
    return B


def _allocate_field_array3(x, y, z):
    """Allocates empty Vector3 data structure

    Args:
        x (float/array): x co-ordinates
        y (float/array): y co-ordinates
        z (float/array): z co-ordinates

    Returns:
        Vector2: Empty data structure
    """

    # Ensure x,y,z are numpy arrays (even of element 1)
    if _np.isscalar(x):
        x = _np.atleast_1d(x)
    if _np.isscalar(y):
        y = _np.atleast_1d(y)
    if _np.isscalar(z):
        z = _np.atleast_1d(z)

    # Determine array shape:
    if _np.ndim(x) == 3:  # Volume meshgrid
        B = Vector3(_np.zeros_like(x), _np.zeros_like(x), _np.zeros_like(x))

    elif _np.ndim(x) == 2:  # planar slice
        B = Vector3(_np.zeros_like(x), _np.zeros_like(x), _np.zeros_like(x))

    elif _np.ndim(y) == 2:  # planar slice
        B = Vector3(_np.zeros_like(y), _np.zeros_like(y), _np.zeros_like(y))

    else:  # line or single point
        B = Vector3(
            _np.zeros(max(x.size, y.size, z.size)),
            _np.zeros(max(x.size, y.size, z.size)),
            _np.zeros(max(x.size, y.size, z.size)),
        )
    return B