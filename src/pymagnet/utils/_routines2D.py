# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Routines for Two Dimensional Magnet Classes"""

import numpy as _np

from ._gradient_kernels import _gradient_2d
from ._vector_structs import Field2, Jacobian2, Point_Array2
from .global_const import MU0


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
        Point_Array2: array of x and y values of shape (num_points, num_points)
        and associated unit
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
        Point_Array2 (Point_Array2): array of x,y points and associated unit,
        defaults to 'mm'

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


_NUMBA_THRESHOLD = 10_000  # use numba for arrays with >= 10k points


def _grid_spacing_2d(x, y):
    """Compute uniform grid spacings from coordinate arrays."""
    Nx, Ny = x.shape
    dx = (x.max() - x.min()) / Nx
    dy = (y.max() - y.min()) / Ny
    return dx, dy


def gradB_2D(B, x, y):
    """Calculates the spatial gradient of the magnetic field magnitude.

    Computes grad(|B|) using finite differences. Uses numba-accelerated
    kernels for grids >= 10k points, np.gradient for smaller grids.

    Args:
        B (ndarray): Magnetic field magnitude |B| (2D array)
        x (ndarray): x coordinates
        y (ndarray): y coordinates

    Returns:
        Field2: Gradient vector (d|B|/dx, d|B|/dy) and its norm
    """
    dx, dy = _grid_spacing_2d(x, y)
    dB = Field2(_np.zeros_like(B), _np.zeros_like(B))
    if B.size >= _NUMBA_THRESHOLD:
        dB.x, dB.y = _gradient_2d(_np.ascontiguousarray(B), dx, dy)
    else:
        dB.x, dB.y = _np.gradient(B, dx, dy)
    dB.calc_norm()
    return dB


def FgradB_2D(B, x, y, chi_m, c):
    """Calculates the magnetic field gradient force for a 2D field.

    Computes F = (chi_m / mu_0) * c * |B| * grad(|B|).
    This is a scalar approximation of the gradient force.

    Args:
        B (Field2): Magnetic field vector (must have .n for magnitude)
        x (ndarray): x coordinates
        y (ndarray): y coordinates
        chi_m (float): Magnetic susceptibility
        c (float): Material constant

    Returns:
        Field2: Magnetic field gradient force vector
    """
    dB = gradB_2D(B.n, x, y)
    scale = (1 / MU0) * chi_m * c
    FB = Field2(_np.zeros_like(B.n), _np.zeros_like(B.n))
    FB.x = scale * dB.x * B.n
    FB.y = scale * dB.y * B.n
    FB.n = scale * dB.n * B.n
    return FB


def jacobian_B_2D(B, x, y):
    """Calculates the Jacobian of the 2D magnetic field vector.

    Computes J_ij = dB_i/dx_j, the full tensor gradient of the vector field.

    Args:
        B (Field2): Magnetic field vector with .x, .y components
        x (ndarray): x coordinates (2D grid)
        y (ndarray): y coordinates (2D grid)

    Returns:
        Jacobian2: Dataclass with components dBx_dx, dBx_dy, dBy_dx, dBy_dy
    """
    dx, dy = _grid_spacing_2d(x, y)
    if B.x.size >= _NUMBA_THRESHOLD:
        dBx_dx, dBx_dy = _gradient_2d(_np.ascontiguousarray(B.x), dx, dy)
        dBy_dx, dBy_dy = _gradient_2d(_np.ascontiguousarray(B.y), dx, dy)
    else:
        dBx_dx, dBx_dy = _np.gradient(B.x, dx, dy)
        dBy_dx, dBy_dy = _np.gradient(B.y, dx, dy)
    return Jacobian2(dBx_dx=dBx_dx, dBx_dy=dBx_dy, dBy_dx=dBy_dx, dBy_dy=dBy_dy)


def BdotgradB_2D(B, x, y):
    """Computes (B . grad)B using the full Jacobian tensor.

    Calculates the directional derivative of B along B:
        [(B·∇)B]_x = Bx * dBx/dx + By * dBx/dy
        [(B·∇)B]_y = Bx * dBy/dx + By * dBy/dy

    This is the force-relevant quantity for magnetic gradient forces
    on paramagnetic materials.

    Args:
        B (Field2): Magnetic field vector with .x, .y components
        x (ndarray): x coordinates (2D grid)
        y (ndarray): y coordinates (2D grid)

    Returns:
        Field2: (B·∇)B vector and its norm
    """
    J = jacobian_B_2D(B, x, y)
    F = Field2(_np.zeros_like(B.x), _np.zeros_like(B.y))
    F.x = B.x * J.dBx_dx + B.y * J.dBx_dy
    F.y = B.x * J.dBy_dx + B.y * J.dBy_dy
    F.calc_norm()
    return F


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
