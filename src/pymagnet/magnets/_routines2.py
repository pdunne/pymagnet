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
    NP = kwargs.pop('NP', 100)
    lx = kwargs.pop('lx', -1*ux)
    ly = kwargs.pop('ly', -1*uy)
    NPJ = NP * 1j
    return(_np.mgrid[lx:ux:NPJ, ly:uy:NPJ])


def B_calc_2D(x, y):
    from ._fields import Vector2
    from ._magnet2 import Magnet_2D
    """Function to calculate magnetic field due to any array of points
       It sums the magnetic field B over each component of the magnetisation
       J = mu_0 M
    """
    if (_np.isscalar(x)):
        x = _np.atleast_1d(x)
    if (_np.isscalar(y)):
        y = _np.atleast_1d(y)
    if (_np.ndim(x) == 2):  # planar slice
        B = Vector2(_np.zeros_like(x), _np.zeros_like(x))
    else:
        B = Vector2(_np.zeros(max(x.size, y.size)), _np.zeros(
            max(x.size, y.size)))  # line or single point

    for magnet in Magnet_2D.instances:

        if _np.fabs(magnet.Jx) > magnet.tol:
            B.x += magnet._calcBx_mag_x(x - magnet.xc, y - magnet.yc)
            B.y += magnet._calcBy_mag_x(x - magnet.xc, y - magnet.yc)

        if _np.fabs(magnet.Jy) > magnet.tol:
            B.x += magnet._calcBx_mag_y(x - magnet.xc, y - magnet.yc)
            B.y += magnet._calcBy_mag_y(x - magnet.xc, y - magnet.yc)
    B.calc_norm()
    return B




def _rotvec2(x, y, theta, xc, yc):
    """Rotation of points x,y

    Args:
        x ([ndarray]): array of x points
        y ([ndarray]): array of y points
        theta ([float]): Rotation angle with respect to x-axis
        xc ([float]): translation x co-orinate
        yc ([float]): translation y co-orinate

    Returns:
        x' [ndarray], y' [ndarry]: tuple of array of rotated points
    """
    theta = _np.deg2rad(theta)

    # Rotation matrix, defined for x-y rotation
    R = _np.array([[_np.cos(theta), -_np.sin(theta), 0],
                   [_np.sin(theta), _np.cos(theta), 0], [0, 0, 1]])

    #Define affine transformation for translation
    a = _np.array([[1, 0, xc], [0, 1, yc], [0, 0, 1]])
    c = _np.array([[1, 0, -xc], [0, 1, -yc], [0, 0, 1]])

    # Transformation matrix
    M = _np.linalg.multi_dot([a, R, c])
    b = _np.column_stack((_np.ravel(x),
                          _np.ravel(y),
                          _np.ones_like(_np.ravel(x))))
    rot = _np.dot(M, b.T)

    # pick out the vectors of rotated x- and y-data
    x_rot = rot[0, :]
    y_rot = rot[1, :]
    return _np.reshape(x_rot, x.shape), _np.reshape(y_rot, y.shape)


def _inv_rotvec2(x, y, theta, xc, yc):
    """Inverse rotation of points x,y

    Args:
        x ([ndarray]): array of x points
        y ([ndarray]): array of y points
        theta ([float]): Rotation angle with respect to x-axis
        xc ([float]): translation x co-orinate
        yc ([float]): translation y co-orinate

    Returns:
        x' [ndarray], y' [ndarry]: tuple of array of rotated points
    """    

    theta = _np.deg2rad(theta)
    # Inverse Rotation matrix, defined for x-y rotation
    R = _np.array([[_np.cos(theta), _np.sin(theta), 0],
                   [-_np.sin(theta), _np.cos(theta), 0], [0, 0, 1]]).T

    #Define affine transformation for translation
    a = _np.array([[1, 0, xc], [0, 1, yc], [0, 0, 1]])
    c = _np.array([[1, 0, -xc], [0, 1, -yc], [0, 0, 1]])

    # Transformation matrix
    M = _np.linalg.multi_dot([a, R, c])
    b = _np.column_stack((_np.ravel(x),
                          _np.ravel(y),
                          _np.ones_like(_np.ravel(x))))
    rot = _np.dot(M, b.T)

    # pick out the vectors of rotated x- and y-data
    x_rot = rot[0, :]
    y_rot = rot[1, :]
    return _np.reshape(x_rot, x.shape), _np.reshape(y_rot, y.shape)
