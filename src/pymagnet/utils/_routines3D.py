# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Routines for Three Dimensional Magnet Classes"""

import numpy as _np

from ._gradient_kernels import _gradient_2d, _gradient_3d
from ._vector_structs import Field3, Jacobian3, Point_Array3
from .global_const import MU0


def plane3D(origin, point_a, point_b, num_points=100):
    """Generates a plane of points defined by

    Args:
        origin (tuple): Origin of plane
        point_a (tuple): Point defining end of vector in x-direction
        point_b (tuple): Point defining end of vector in y-direction
        num_points (int, optional): Number of points in each direction, for a total of num_points^2. Defaults to 100.

    Returns:
        Point_Array3: Struct containing the x,y,z coordinates
    """

    vector_a = point_a - origin
    vector_b = point_b - origin
    num_points_j = num_points * 1j
    x_elements, y_elements = _np.mgrid[0:1:num_points_j, 0:1:num_points_j]
    x = origin[0] + x_elements * vector_a[0] + y_elements * vector_b[0]
    y = origin[1] + x_elements * vector_a[1] + y_elements * vector_b[1]
    z = origin[2] + x_elements * vector_a[2] + y_elements * vector_b[2]

    return Point_Array3(x=x, y=y, z=z)


def grid3D(xmax, ymax, zmax, **kwargs):
    """Generates grid of x, y, z points

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
    # NPJ = num_points * 1j

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


def point3D(point, **kwargs):
    """Returns a single point

    Args:
        point (tuple):
        unit (str, optional): length scale units. Defaults to "mm".

    Returns:
        Point_Array3: struct of x, y, and z values of shape (1) and associated unit
    """

    unit = kwargs.pop("unit", "mm")

    return Point_Array3(
        x=point[0],
        y=point[1],
        z=point[2],
        unit=unit,
    )


def line3D(start, end, num_points=100, **kwargs):
    """Generates a line of points

    Args:
        start (tuple): Starting point (x1,y1,z1)
        end (tuple): End point (x2,y2,z2)
        num_points (int): number of points to generate. Defaults to 100
        unit (str, optional): length scale units. Defaults to "mm".

    Returns:
        Point_Array3: array of x, y, and z values of shape (num_points) and associated unit
    """

    unit = kwargs.pop("unit", "mm")

    return Point_Array3(
        x=_np.linspace(start[0], end[0], num_points),
        y=_np.linspace(start[1], end[1], num_points),
        z=_np.linspace(start[2], end[2], num_points),
        unit=unit,
    )


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

    # Track which points are inside any magnet (for masking after summation)
    mask_magnet = _np.zeros(B.x.shape, dtype=bool)
    needs_masking = False

    for magnet in Magnet3D.instances:
        Bx, By, Bz = magnet.get_field(points.x, points.y, points.z)
        Bx = Bx.reshape(B.x.shape)
        By = By.reshape(B.y.shape)
        Bz = Bz.reshape(B.z.shape)

        # Replace NaN with zero before accumulation to prevent NaN
        # poisoning the sum (NaN + valid = NaN). NaN can arise from
        # mask_magnet flagging points inside a magnet, or from
        # numerical singularities at magnet edges.
        nan_mask = _np.isnan(Bx) | _np.isnan(By) | _np.isnan(Bz)
        if _np.any(nan_mask):
            Bx = _np.where(nan_mask, 0.0, Bx)
            By = _np.where(nan_mask, 0.0, By)
            Bz = _np.where(nan_mask, 0.0, Bz)

            if magnet._mask_magnet:
                mask_magnet |= nan_mask
                needs_masking = True

        B.x += Bx
        B.y += By
        B.z += Bz

    # Re-apply NaN mask for points inside any magnet
    if needs_masking:
        B.x[mask_magnet] = _np.nan
        B.y[mask_magnet] = _np.nan
        B.z[mask_magnet] = _np.nan

    B.calc_norm()
    return B


def _detect_slice_axes(x, y, z):
    """Identify which axes vary in a 2D slice from slice3D().

    Returns:
        tuple: (varying_axes, constant_axis) where varying_axes is a list
        of two axis labels ('x', 'y', 'z') and constant_axis is the fixed one.
    """
    ranges = {
        "x": float(x.max() - x.min()),
        "y": float(y.max() - y.min()),
        "z": float(z.max() - z.min()),
    }
    coords = {"x": x, "y": y, "z": z}
    constant = min(ranges, key=ranges.get)
    varying = [k for k in ("x", "y", "z") if k != constant]
    return varying, constant, coords


_NUMBA_THRESHOLD = 10_000  # use numba for arrays with >= 10k points


def _grid_spacing_3d(x, y, z):
    """Compute uniform grid spacings from 3D coordinate arrays."""
    Nx, Ny, Nz = x.shape
    dx = (x.max() - x.min()) / Nx
    dy = (y.max() - y.min()) / Ny
    dz = (z.max() - z.min()) / Nz
    return dx, dy, dz


def gradB_3D(B, x, y, z):
    """Calculates the spatial gradient of the magnetic field magnitude.

    Computes grad(|B|) using finite differences. Uses numba-accelerated
    kernels for grids >= 10k points, np.gradient for smaller grids.
    Supports both full 3D grids (from grid3D) and 2D slices (from slice3D).

    Args:
        B (ndarray): Magnetic field magnitude |B|
        x (ndarray): x coordinates
        y (ndarray): y coordinates
        z (ndarray): z coordinates

    Returns:
        Field3: Gradient vector (d|B|/dx, d|B|/dy, d|B|/dz) and its norm
    """
    B_arr = _np.ascontiguousarray(B)
    use_numba = B_arr.size >= _NUMBA_THRESHOLD

    if B_arr.ndim == 3:
        dx, dy, dz = _grid_spacing_3d(x, y, z)
        if use_numba:
            dBdx, dBdy, dBdz = _gradient_3d(B_arr, dx, dy, dz)
        else:
            dBdx, dBdy, dBdz = _np.gradient(B_arr, dx, dy, dz)
    elif B_arr.ndim == 2:
        varying, constant, coords = _detect_slice_axes(x, y, z)
        c1, c2 = coords[varying[0]], coords[varying[1]]
        N1, N2 = B_arr.shape
        d1 = (c1.max() - c1.min()) / N1
        d2 = (c2.max() - c2.min()) / N2
        if use_numba:
            grad1, grad2 = _gradient_2d(B_arr, d1, d2)
        else:
            grad1, grad2 = _np.gradient(B_arr, d1, d2)

        components = {"x": _np.zeros_like(B_arr), "y": _np.zeros_like(B_arr), "z": _np.zeros_like(B_arr)}
        components[varying[0]] = grad1
        components[varying[1]] = grad2
        dBdx, dBdy, dBdz = components["x"], components["y"], components["z"]
    else:
        raise ValueError(f"B must be 2D (slice) or 3D (grid), got ndim={B_arr.ndim}")

    dB = Field3(dBdx, dBdy, dBdz)
    dB.calc_norm()
    return dB


def FgradB_3D(B, x, y, z, chi_m, c):
    """Calculates the magnetic field gradient force for a 3D field.

    Computes F = (chi_m / mu_0) * c * |B| * grad(|B|).
    This is a scalar approximation of the gradient force.

    Args:
        B (Field3): Magnetic field vector (must have .n for magnitude)
        x (ndarray): x coordinates
        y (ndarray): y coordinates
        z (ndarray): z coordinates
        chi_m (float): Magnetic susceptibility
        c (float): Material constant

    Returns:
        Field3: Magnetic field gradient force vector
    """
    dB = gradB_3D(B.n, x, y, z)
    scale = (1 / MU0) * chi_m * c
    FB = Field3(_np.zeros_like(B.n), _np.zeros_like(B.n), _np.zeros_like(B.n))
    FB.x = scale * dB.x * B.n
    FB.y = scale * dB.y * B.n
    FB.z = scale * dB.z * B.n
    FB.n = scale * dB.n * B.n
    return FB


def jacobian_B_3D(B, x, y, z):
    """Calculates the Jacobian of the 3D magnetic field vector.

    Computes J_ij = dB_i/dx_j, the full tensor gradient of the vector field.
    Supports both full 3D grids (from grid3D) and 2D slices (from slice3D).

    Args:
        B (Field3): Magnetic field vector with .x, .y, .z components
        x (ndarray): x coordinates
        y (ndarray): y coordinates
        z (ndarray): z coordinates

    Returns:
        Jacobian3: Dataclass with all 9 partial derivative components
    """
    Bx = _np.ascontiguousarray(B.x)
    By = _np.ascontiguousarray(B.y)
    Bz = _np.ascontiguousarray(B.z)
    use_numba = Bx.size >= _NUMBA_THRESHOLD

    if Bx.ndim == 3:
        dx, dy, dz = _grid_spacing_3d(x, y, z)
        if use_numba:
            dBx_dx, dBx_dy, dBx_dz = _gradient_3d(Bx, dx, dy, dz)
            dBy_dx, dBy_dy, dBy_dz = _gradient_3d(By, dx, dy, dz)
            dBz_dx, dBz_dy, dBz_dz = _gradient_3d(Bz, dx, dy, dz)
        else:
            dBx_dx, dBx_dy, dBx_dz = _np.gradient(Bx, dx, dy, dz)
            dBy_dx, dBy_dy, dBy_dz = _np.gradient(By, dx, dy, dz)
            dBz_dx, dBz_dy, dBz_dz = _np.gradient(Bz, dx, dy, dz)
    elif Bx.ndim == 2:
        varying, constant, coords = _detect_slice_axes(x, y, z)
        c1, c2 = coords[varying[0]], coords[varying[1]]
        N1, N2 = Bx.shape
        d1 = (c1.max() - c1.min()) / N1
        d2 = (c2.max() - c2.min()) / N2

        zeros = _np.zeros_like(Bx)
        grad_fn = _gradient_2d if use_numba else lambda F, d1, d2: _np.gradient(F, d1, d2)
        # Compute gradients for each B component on the 2D slice
        result = {}
        for comp_name, comp_arr in [("Bx", Bx), ("By", By), ("Bz", Bz)]:
            g1, g2 = grad_fn(comp_arr, d1, d2)
            grads = {"x": zeros.copy(), "y": zeros.copy(), "z": zeros.copy()}
            grads[varying[0]] = g1
            grads[varying[1]] = g2
            result[comp_name] = grads

        dBx_dx, dBx_dy, dBx_dz = result["Bx"]["x"], result["Bx"]["y"], result["Bx"]["z"]
        dBy_dx, dBy_dy, dBy_dz = result["By"]["x"], result["By"]["y"], result["By"]["z"]
        dBz_dx, dBz_dy, dBz_dz = result["Bz"]["x"], result["Bz"]["y"], result["Bz"]["z"]
    else:
        raise ValueError(f"B components must be 2D (slice) or 3D (grid), got ndim={Bx.ndim}")

    return Jacobian3(
        dBx_dx=dBx_dx, dBx_dy=dBx_dy, dBx_dz=dBx_dz,
        dBy_dx=dBy_dx, dBy_dy=dBy_dy, dBy_dz=dBy_dz,
        dBz_dx=dBz_dx, dBz_dy=dBz_dy, dBz_dz=dBz_dz,
    )


def BdotgradB_3D(B, x, y, z):
    """Computes (B . grad)B using the full Jacobian tensor.

    Calculates the directional derivative of B along B:
        [(B·∇)B]_x = Bx * dBx/dx + By * dBx/dy + Bz * dBx/dz
        [(B·∇)B]_y = Bx * dBy/dx + By * dBy/dy + Bz * dBy/dz
        [(B·∇)B]_z = Bx * dBz/dx + By * dBz/dy + Bz * dBz/dz

    Supports both full 3D grids and 2D slices.

    Args:
        B (Field3): Magnetic field vector with .x, .y, .z components
        x (ndarray): x coordinates
        y (ndarray): y coordinates
        z (ndarray): z coordinates

    Returns:
        Field3: (B·∇)B vector and its norm
    """
    J = jacobian_B_3D(B, x, y, z)
    F = Field3(_np.zeros_like(B.x), _np.zeros_like(B.y), _np.zeros_like(B.z))
    F.x = B.x * J.dBx_dx + B.y * J.dBx_dy + B.z * J.dBx_dz
    F.y = B.x * J.dBy_dx + B.y * J.dBy_dy + B.z * J.dBy_dz
    F.z = B.x * J.dBz_dx + B.y * J.dBz_dy + B.z * J.dBz_dz
    F.calc_norm()
    return F


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
    if _np.ndim(x) == 3 or _np.ndim(x) == 2:  # Volume meshgrid
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
    from ..magnets import Cylinder, Prism, Sphere

    J = magnet.get_Jr()
    mask_magnet = magnet._mask_magnet

    if mask_magnet:
        field.x[mask] = _np.nan
        field.y[mask] = _np.nan
        field.z[mask] = _np.nan

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
