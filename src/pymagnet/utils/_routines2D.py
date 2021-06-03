# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Routines for Two Dimensional Magnet Classes
"""
import numpy as _np
from ._vector_structs import Point_Array2, Field2
from .global_const import MU0

__all__ = ["grid2D", "get_field_2D", "rotate_points_2D"]


def grid2D(xmax, ymax, **kwargs):
    """Generates grid of x and y points

    Args:
        xmax (float): maximum x value
        ymax (float): maximum y value

    Kwargs:
        num_points (int): Number of points in each direction. Defaults to 100
        xmin (float): minimum x value. Defaults to -xmax
        ymin (float): minimum y value. Defaults to -ymax
        unit (str): unit length. Defaults to 'mm'

    Returns:
        Point_Array2: array of x and y values of shape (num_points, num_points) and associated unit
    """
    num_points = kwargs.pop("num_points", 100)
    xmin = kwargs.pop("xmin", -1 * xmax)
    ymin = kwargs.pop("ymin", -1 * ymax)
    unit = kwargs.pop("unit", "mm")
    NPJ = num_points * 1j
    x, y = _np.mgrid[xmin:xmax:NPJ, ymin:ymax:NPJ]
    return Point_Array2(x, y, unit=unit)


def get_field_2D(Point_Array2):
    """Calculates magnetic field at an array of points due to every instantated
    `Magnet2D` magnet.

    Args:
        Point_Array2 (Point_Array2): array of x,y points and associated unit, defaults to 'mm'

    Returns:
        Field2: array of Bx,By,|B| values and associated unit (defaults to 'T')
    """
    from ..magnets import Magnet2D

    # Empty data structure
    B = _allocate_field_array2(Point_Array2.x, Point_Array2.y)

    for magnet in Magnet2D.instances:
        Bx, By = magnet.get_field(Point_Array2.x, Point_Array2.y)
        B.x += Bx
        B.y += By

    B.calc_norm()
    return B


def rotate_points_2D(x, y, alpha):
    """Counter-clockwise rotation of points x,y

    Rotates 2D coordinates using a rotation matrix

    Args:
        x (ndarray): array of x coordinates
        y (ndarray): array of x coordinates
        alpha (float): rotation angle w.r.t. x-axis

    Returns:
        tuple: (x', y') rotated array of points
    """
    x = _np.atleast_1d(x)
    y = _np.atleast_1d(y)
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
        B (Field2): Magnetic field vector
        x (ndarray): x coordinates
        y (ndarray): y coordinates

    Returns:
        Field2: Magnetic field gradient vector
    """

    dB = Field2(_np.zeros_like(B), _np.zeros_like(B))
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
        B (Field2): Magnetic field vector
        x (ndarray): x coordinates
        y (ndarray): y coordinates

    Returns:
        Field2: Magnetic field gradient force vector
    """

    BgB = Field2(_np.zeros_like(B.n), _np.zeros_like(B.n))
    FB = Field2(_np.zeros_like(B.n), _np.zeros_like(B.n))

    dB = gradB_2D(B, x, y)
    BgB.n = dB.n * B.n
    BgB.x = dB.x * B.n
    BgB.y = dB.y * B.n
    FB.n = (1 / MU0) * chi_m * c * BgB.n
    FB.x = (1 / MU0) * chi_m * c * BgB.x
    FB.y = (1 / MU0) * chi_m * c * BgB.y
    return FB


def _allocate_field_array2(x, y):
    """Allocates empty Field2 data structure

    Args:
        x (ndarray): x co-ordinates
        y (ndarray): y co-ordinates

    Returns:
        Field2: Empty data structure
    """

    # Ensure x,y,z are numpy arrays (even of element 1)
    if _np.isscalar(x):
        x = _np.atleast_1d(x)
    if _np.isscalar(y):
        y = _np.atleast_1d(y)

    # Determine array shape:
    if _np.ndim(x) == 2:  # planar slice
        B = Field2(_np.zeros_like(x), _np.zeros_like(x))
    else:  # line or single point
        B = Field2(_np.zeros(max(x.size, y.size)), _np.zeros(max(x.size, y.size)))
    return B


def _get_field_array_shape2(x, y):
    """Allocates empty Field2 data structure

    Args:
        x (ndarray): x co-ordinates
        y (ndarray): y co-ordinates

    Returns:
        Field2: Empty data structure
    """

    # Ensure x,y,z are numpy arrays (even of element 1)
    x = _np.atleast_1d(x)
    y = _np.atleast_1d(y)

    # Determine array shape:
    if _np.ndim(x) == 2:  # planar slice
        array_shape = x.shape
    else:  # line or single point
        array_shape = max(x.size, y.size)
    return array_shape