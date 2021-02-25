# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
# __all__ = ['B_calc_2D']
"""Routines for Two Dimensional Magnet Classes
"""
import numpy as _np


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
    from ._fields import Vector2
    from ._magnet2 import Magnet_2D

    """Function to calculate magnetic field due to any array of points
       It sums the magnetic field B over each component of the magnetisation
       J = mu_0 M
    """
    if _np.isscalar(x):
        x = _np.atleast_1d(x)
    if _np.isscalar(y):
        y = _np.atleast_1d(y)
    if _np.ndim(x) == 2:  # planar slice
        B = Vector2(_np.zeros_like(x), _np.zeros_like(x))
    else:
        B = Vector2(
            _np.zeros(max(x.size, y.size)), _np.zeros(max(x.size, y.size))
        )  # line or single point

    for magnet in Magnet_2D.instances:
        if magnet.alpha_radians > Magnet_2D.tol:
            xi, yi = rotate_points_2D(x, y, magnet.alpha_radians)
            xci, yci = rotate_points_2D(
                _np.array([magnet.xc]), _np.array([magnet.yc]), magnet.alpha_radians
            )

        if _np.fabs(magnet.Jx) > Magnet_2D.tol:
            if magnet.alpha_radians > Magnet_2D.tol:

                # Calculate fields in local frame
                Btx = magnet._calcBx_mag_x(xi - xci[0], yi - yci[0])
                Bty = magnet._calcBy_mag_x(xi - xci[0], yi - yci[0])

                # Rotate fields to global frame
                Bxt, Byt = rotate_points_2D(Btx, Bty, 2 * _np.pi - magnet.alpha_radians)

                B.x += Bxt
                B.y += Byt
            else:
                B.x += magnet._calcBx_mag_x(x - magnet.xc, y - magnet.yc)
                B.y += magnet._calcBy_mag_x(x - magnet.xc, y - magnet.yc)

        if _np.fabs(magnet.Jy) > magnet.tol:
            if magnet.alpha_radians > Magnet_2D.tol:

                Btx = magnet._calcBx_mag_y(xi - magnet.xc - xci[0], yi - yci[0])
                Bty = magnet._calcBy_mag_y(xi - magnet.xc - xci[0], yi - yci[0])

                Bxt, Byt = rotate_points_2D(Btx, Bty, 2 * _np.pi - magnet.alpha_radians)
                B.x += Bxt
                B.y += Byt
            else:

                B.x += magnet._calcBx_mag_y(x - magnet.xc, y - magnet.yc)
                B.y += magnet._calcBy_mag_y(x - magnet.xc, y - magnet.yc)
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
