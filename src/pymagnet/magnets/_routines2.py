# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne

# __all__ = ['B_calc_2D']
"""Routines for Two Dimensional Magnet Classes
"""
import numpy as _np
from ._fields import Vector2


def grid2D(ux, uy, **kwargs):
    """Grid of x,y points.
    Args:

        ux ([float]): [x upper limit]
        uy ([float]): [y upper limit]
        lx ([float, optional]): [x lower limit defaults to -ux]
        ly ([optional]): [x lower limit defaults to -ux]
        NP (int, optional): [number of points in each direction]. Defaults to 100.
    """
    NP = kwargs.pop("NP", 100)
    lx = kwargs.pop("lx", -1 * ux)
    ly = kwargs.pop("ly", -1 * uy)
    NPJ = NP * 1j
    return _np.mgrid[lx:ux:NPJ, ly:uy:NPJ]


def B_calc_2D(x, y):
    """Function to calculate magnetic field due to any array of points
    It sums the magnetic field B over each component of the magnetisation
    J = mu_0 M
    """
    from ._magnet2 import Magnet_2D

    # Empty data structure
    B = _allocate_field_array2(x, y)

    for magnet in Magnet_2D.instances:
        Bx, By = magnet.calcB(x, y)
        B.x += Bx
        B.y += By

    B.calc_norm()
    return B


def rotate_points_2D(x, y, alpha):
    """Counter-clockwise rotation of points x,y

    Rotates 2D coordinates using a rotation matrix

    Args:
        x (float/array): array of x coordinates
        y (float/array): array of x coordinates
        alpha (float): rotation angle w.r.t. x-axis

    Returns:
        tuple: (x', y') rotated array of points
    """
    if len(x) != len(y):
        raise Exception("Must have same number of points in x and y")

    rot_matrix = _np.array(
        [[_np.cos(alpha), -_np.sin(alpha)], [_np.sin(alpha), _np.cos(alpha)]]
    )
    stacked_points = _np.column_stack((_np.ravel(x), _np.ravel(y)))
    rotated_points = _np.dot(rot_matrix, stacked_points.T)
    x_rotated = rotated_points[0, :]
    y_rotated = rotated_points[1, :]

    return _np.reshape(x_rotated, x.shape), _np.reshape(y_rotated, y.shape)


def gradB_2D(B, x, y):
    """Calculates the magnetic field gradient for a 2D field.

    Args:
        B (Vector2): Magnetic field vector
        x (float/array): x coordinates
        y (float/array): y coordinates

    Returns:
        Vector2: Magnetic field gradient vector
    """
    from ._fields import Vector2

    dB = Vector2(_np.zeros_like(B), _np.zeros_like(B))
    Nx = x.shape[0]
    Ny = x.shape[1]
    dx = (x.max() - x.min()) / Nx
    dy = (y.max() - y.min()) / Ny
    dB.x, dB.y = _np.gradient(B, dx, dy)
    dB.calc_norm()
    return dB


def FgradB_2D(B, x, y, chi_m, c):
    """Calculates the magnetic field gradient force for a 2D field.

    Args:
        B (Vector2): Magnetic field vector
        x (float/array): x coordinates
        y (float/array): y coordinates

    Returns:
        Vector2: Magnetic field gradient force vector
    """
    from ._fields import Vector2
    from .. import u0

    BgB = Vector2(_np.zeros_like(B.n), _np.zeros_like(B.n))
    FB = Vector2(_np.zeros_like(B.n), _np.zeros_like(B.n))

    dB = gradB_2D(B, x, y)
    BgB.n = dB.n * B.n
    BgB.x = dB.x * B.n
    BgB.y = dB.y * B.n
    FB.n = (1 / u0) * chi_m * c * BgB.n
    FB.x = (1 / u0) * chi_m * c * BgB.x
    FB.y = (1 / u0) * chi_m * c * BgB.y
    return FB


def _allocate_field_array2(x, y):
    """Allocates empty Vector2 data structure

    Args:
        x (float/array): x co-ordinates
        y (float/array): y co-ordinates

    Returns:
        Vector2: Empty data structure
    """

    # Ensure x,y,z are numpy arrays (even of element 1)
    if _np.isscalar(x):
        x = _np.atleast_1d(x)
    if _np.isscalar(y):
        y = _np.atleast_1d(y)

    # Determine array shape:
    if _np.ndim(x) == 2:  # planar slice
        B = Vector2(_np.zeros_like(x), _np.zeros_like(x))
    else:  # line or single point
        B = Vector2(_np.zeros(max(x.size, y.size)), _np.zeros(max(x.size, y.size)))
    return B


def _get_field_array_shape2(x, y):
    """Allocates empty Vector2 data structure

    Args:
        x (float/array): x co-ordinates
        y (float/array): y co-ordinates

    Returns:
        Vector2: Empty data structure
    """

    # Ensure x,y,z are numpy arrays (even of element 1)
    if _np.isscalar(x):
        x = _np.atleast_1d(x)
    if _np.isscalar(y):
        y = _np.atleast_1d(y)

    # Determine array shape:
    if _np.ndim(x) == 2:  # planar slice
        array_shape = x.shape
    else:  # line or single point
        array_shape = max(x.size, y.size)
    return array_shape