# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Routines for Three Dimensional Magnet Classes
"""
__all__ = ["get_field_3D", "grid3D", "slice3D"]

import numpy as _np
from ._vector_structs import Field3, Point_Array3


def grid3D(xmax, ymax, zmax, **kwargs):
    """Generates grid of x and y points

    Args:
        xmax (float): maximum x value
        ymax (float): maximum y value
        zmax (float): maximum y value

    Kwargs:
        num_points (int): Number of points in each direction. Defaults to 100
        xmin (float): minimum x value. Defaults to -xmax
        ymin (float): minimum y value. Defaults to -ymax
        zmin (float): minimum y value. Defaults to -zmax
        unit (string): unit length. Defaults to 'mm'

    Returns:
        Point_Array2: array of x and y values of shape (num_points, num_points) and associated unit
    """
    num_points = kwargs.pop("num_points", None)

    xmin = kwargs.pop("xmin", -1 * xmax)
    ymin = kwargs.pop("ymin", -1 * ymax)
    zmin = kwargs.pop("zmin", -1 * zmax)
    unit = kwargs.pop("unit", "mm")
    NPJ = num_points * 1j

    if num_points is None:
        num_points_x = kwargs.pop("num_points_x", 100)
        num_points_y = kwargs.pop("num_points_y", 100)
        num_points_z = kwargs.pop("num_points_z", 100)
    else:
        num_points_x = num_points
        num_points_y = num_points
        num_points_z = num_points

    x, y, z = _np.mgrid[
        xmin : xmax : num_points_x * 1j,
        ymin : ymax : num_points_y * 1j,
        zmin : zmax : num_points_z * 1j,
    ]

    return Point_Array3(x, y, z, unit=unit)


def slice3D(plane="xy", max1=1.0, max2=1.0, slice_value=0.0, unit="mm", **kwargs):
    """Generates a planar slice of values

    Args:
        plane (str, optional): plane. Defaults to "xy".
        max1 (float, optional): maximum along axis 1. Defaults to 1.0.
        max2 (float, optional): maximum along axis 2. Defaults to 1.0.
        slice_value (float, optional): constant value for third axis. Defaults to 0.0.
        unit (str, optional): length scale units. Defaults to "mm".

    Kwargs:
        num_points (int): Number of points in each direction. Defaults to 100
        min1 (float): minimum along axis 1. Defaults to -min1
        min2 (float): minimum along axis 2. Defaults to -min2

    Raises:
        Exception: plane type, must be one of 'xy', 'xz, 'yz', or 'custom'

    Returns:
        Point_Array3: array of x, y, and z values of shape (num_points, num_points) and associated unit
    """
    num_points = kwargs.pop("num_points", 100)
    min1 = kwargs.pop("min1", -1 * max1)
    min2 = kwargs.pop("min2", -1 * max2)
    NPj = num_points * 1j

    if plane.lower() == "xy":
        x, y = _np.mgrid[min1:max1:NPj, min2:max2:NPj]
        z = _np.asarray([slice_value])
        z = _np.tile(z, x.shape)

    elif plane.lower() == "xz":
        x, z = _np.mgrid[min1:max1:NPj, min2:max2:NPj]
        y = _np.asarray([slice_value])
        y = _np.tile(y, x.shape)

    elif plane.lower() == "yz":
        y, z = _np.mgrid[min1:max1:NPj, min2:max2:NPj]
        x = _np.asarray([slice_value])
        x = _np.tile(x, y.shape)

    elif plane.lower() == "custom":
        x = kwargs.pop("custom_x", _np.array([0.0]))
        y = kwargs.pop("custom_y", _np.array([0.0]))
        z = kwargs.pop("custom_z", _np.array([0.0]))

    else:
        raise Exception("plane must be one of 'xy', 'xz, 'yz', or 'custom'")

    return Point_Array3(x, y, z, unit=unit)


def get_field_3D(points):
    """Calculates magnetic field at an array of points due to every instantated
    `Magnet3D` magnet.

    Args:
        Point_Array3 (Point_Array3): array of x,y,z points and associated unit, defaults to 'mm'

    Returns:
        Field3: array of Bx,By,Bz,|B| values and associated unit (defaults to 'T')
    """
    from ..magnets import Magnet3D

    B = _allocate_field_array3(points.x, points.y, points.z)

    for magnet in Magnet3D.instances:
        Bx, By, Bz = magnet.get_field(points.x, points.y, points.z)
        B.x += Bx.reshape(B.x.shape)
        B.y += By.reshape(B.y.shape)
        B.z += Bz.reshape(B.z.shape)

    B.calc_norm()
    return B


def _allocate_field_array3(x, y, z):
    """Allocates empty Field3 data structure

    Args:
        x (ndarray): x co-ordinates
        y (ndarray): y co-ordinates
        z (ndarray): z co-ordinates

    Returns:
        Field3: Empty data structure
    """

    # Ensure x,y,z are numpy arrays (even of element 1)
    x = _np.atleast_1d(x)
    y = _np.atleast_1d(y)
    z = _np.atleast_1d(z)

    # Determine array shape:
    if _np.ndim(x) == 3:  # Volume meshgrid
        B = Field3(_np.zeros_like(x), _np.zeros_like(x), _np.zeros_like(x))

    elif _np.ndim(x) == 2:  # planar slice
        B = Field3(_np.zeros_like(x), _np.zeros_like(x), _np.zeros_like(x))

    elif _np.ndim(y) == 2:  # planar slice
        B = Field3(_np.zeros_like(y), _np.zeros_like(y), _np.zeros_like(y))

    else:  # line or single point
        B = Field3(
            _np.zeros(max(x.size, y.size, z.size)),
            _np.zeros(max(x.size, y.size, z.size)),
            _np.zeros(max(x.size, y.size, z.size)),
        )
    return B


def _get_max_array(list_arrays):
    """Gets the shape and size of the largest array in a list

    Args:
        list_arrays (list): list of numpy ndarrays

    Returns:
        tuple: max_shape (tuple), max_size (int)
    """
    max_shape = (1, 1)
    max_size = 0
    for array in list_arrays:
        if array.size > max_size:
            max_shape = array.shape
            max_size = array.size
    return max_shape, max_size


def _check_tile_array(array, max_size, max_shape):
    if array.size < max_size:
        return _np.tile(array, max_shape)
    else:
        return array


def _tile_arrays(x, y, z):
    """Tiles arrays to match the largest array of x, y, z

    Args:
        x (ndarray/float): x-coordinates
        y (ndarray/float): y-coordinates
        z (ndarray/float): z-coordinates

    Returns:
        tuple: x, y, z ndarrays
    """
    x = _np.asarray(x)
    y = _np.asarray(y)
    z = _np.asarray(z)

    max_shape, max_size = _get_max_array([x, y, z])

    x = _check_tile_array(x, max_size, max_shape)
    y = _check_tile_array(y, max_size, max_shape)
    z = _check_tile_array(z, max_size, max_shape)

    return x, y, z


def _apply_mask(magnet, field, mask):
    """Calculates B, or applies Nan, inside magnets

    Args:
        magnet (Magnet3D): instantiated magnet
        field (Field3): magnetic field vector
        mask (array): boolean array where True marks coordinates inside the magnet

    Returns:
        Field3: masked magnetic field vector
    """
    from ..magnets import Prism, Cylinder, Sphere

    J = magnet.get_Jr()
    mask_magnet = magnet._mask_magnet

    if mask_magnet:
        field.x[mask] = _np.NaN
        field.y[mask] = _np.NaN
        field.z[mask] = _np.NaN

    else:
        if issubclass(magnet.__class__, Prism):
            field.x[mask] += J[0]
            field.y[mask] += J[1]
            field.z[mask] += J[2]

        elif issubclass(magnet.__class__, Cylinder):
            pass

        elif issubclass(magnet.__class__, Sphere):
            field.x[mask] = 0.0
            field.y[mask] = 0.0
            field.z[mask] = magnet.Jr * 2 / 3

    return field