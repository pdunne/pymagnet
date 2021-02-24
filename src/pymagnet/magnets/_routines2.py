# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
# __all__ = ['B_calc_2D']
"""Routines for Two Dimensional Magnet Classes
"""
import numpy as _np
from ..magnets import u0
from ..magnets import PI
from matplotlib.path import Path as _Path


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
            xi, yi = rotate_points(x, y, magnet.alpha_radians)
            xci, yci = rotate_points(
                _np.array([magnet.xc]), _np.array([magnet.yc]), magnet.alpha_radians
            )

        if _np.fabs(magnet.Jx) > Magnet_2D.tol:
            if magnet.alpha_radians > Magnet_2D.tol:

                # Calculate fields in local frame
                Btx = magnet._calcBx_mag_x(xi - xci[0], yi - yci[0])
                Bty = magnet._calcBy_mag_x(xi - xci[0], yi - yci[0])

                # Rotate fields to global frame
                Bxt, Byt = rotate_points(Btx, Bty, 2 * _np.pi - magnet.alpha_radians)

                B.x += Bxt
                B.y += Byt
            else:
                B.x += magnet._calcBx_mag_x(x - magnet.xc, y - magnet.yc)
                B.y += magnet._calcBy_mag_x(x - magnet.xc, y - magnet.yc)

        if _np.fabs(magnet.Jy) > magnet.tol:
            if magnet.alpha_radians > Magnet_2D.tol:

                Btx = magnet._calcBx_mag_y(xi - magnet.xc - xci[0], yi - yci[0])
                Bty = magnet._calcBy_mag_y(xi - magnet.xc - xci[0], yi - yci[0])

                Bxt, Byt = rotate_points(Btx, Bty, 2 * _np.pi - magnet.alpha_radians)
                B.x += Bxt
                B.y += Byt
            else:

                B.x += magnet._calcBx_mag_y(x - magnet.xc, y - magnet.yc)
                B.y += magnet._calcBy_mag_y(x - magnet.xc, y - magnet.yc)
    B.calc_norm()
    return B


def rotate_points(x, y, alpha):
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


# def gradB(B, x, y):
#     """ compute the gradient of Magnetic Field """
#     dB = field_vec(_np.zeros_like(B), _np.zeros_like(B))
#     Nx = x.shape[0]
#     Ny = x.shape[1]
#     dx = (x.max() - x.min()) / Nx
#     dy = (y.max() - y.min()) / Ny
#     dB.x, dB.y = _np.gradient(B, dx, dy)
#     dB.abmag()
#     return dB


# def FgradB(B, x, y, chi_m, c):
#     """ compute the field gradient force """
#     BgB = field_vec(_np.zeros_like(B.n), _np.zeros_like(B.n))
#     FB = field_vec(_np.zeros_like(B.n), _np.zeros_like(B.n))
#     # Nx = x.shape[0]
#     # Ny = x.shape[1]
#     # dx = (x.max()-x.min())/Nx
#     # dy = (y.max()-y.min())/Ny
#     dB = gradB(B.n, x, y)
#     BgB.n = dB.n * B.n
#     BgB.x = dB.x * B.n
#     BgB.y = dB.y * B.n
#     FB.n = (1 / u0) * chi_m * c * BgB.n
#     FB.x = (1 / u0) * chi_m * c * BgB.x
#     FB.y = (1 / u0) * chi_m * c * BgB.y
#     return FB


# def mask_poly(N, apo, xc, yc, x, y):
#     k = _np.arange(0, N, 1)
#     # Case Structure for setting angle offset

#     def f(N):
#         if N % 2 == 0:
#             return {
#                 # 3: PI / N.
#                 # 4: PI / N,
#                 # 5: PI / N,
#                 6: 0,
#                 # 8: PI / N,
#             }.get(N, PI / N)
#         else:
#             return {
#                 3: PI / 2,
#                 7: PI / 2,
#                 # 5: PI / N,
#                 # 6: 0,
#                 # 8: PI / N,
#             }.get(N, PI / N / 2)

#     r = apo / _np.around(_np.cos(PI / N), 4)

#     xv = xc + r * _np.sin(2 * PI * k / N + f(N))
#     yv = yc + r * _np.cos(2 * PI * k / N + f(N))
#     poly_verts = _np.vstack((xv, yv)).T.tolist()
#     points = _np.vstack((x.flatten(), y.flatten())).T

#     path = _Path(poly_verts)
#     grid = path.contains_points(points)
#     grid = grid.reshape(x.shape)
#     return grid


# def mask_data(a, N, R1, B, x, y):
#     apo = R1 - a
#     xp = 0
#     yp = 0
#     grid = mask_poly(N, apo, xp, yp, x, y)
#     Bm = field_vec(B.x, B.y)
#     # Bmx = B.x.copy()
#     # Bmy = B.y.copy
#     Bm.x[~grid] = _np.nan
#     Bm.y[~grid] = _np.nan
#     # Bm = mag.field_vec(Bmx,Bmy)
#     Bm.abmag()
#     return Bm


# def radius_mask_data(x, y, B, R):
#     """Circular Data Mask

#     Args:
#         x ([array]): array of x coordinates
#         y ([array]): array of x coordinates
#         B ([field_vec object]): field_vec object
#                                 contains B.x, B.y, B.n
#         R ([float]): Radius of circular mask
#     """
#     r_data = _np.sqrt(x ** 2 + y ** 2)
#     B_radius_mask = field_vec(B.x, B.y)
#     B_radius_mask.x[((r_data) > R)] = _np.NaN
#     B_radius_mask.y[((r_data) > R)] = _np.NaN
#     B_radius_mask.abmag()
#     return B_radius_mask


# def radial_profile(data, center):
#     """[summary]

#     Args:
#         data ([type]): [description]
#         center ([type]): [description]

#     Returns:
#         [type]: [description]
#     """
#     y, x = _np.indices((data.shape))
#     r = _np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
#     r = r.astype(_np.int)

#     tbin = _np.bincount(r.ravel(), data.ravel())
#     nr = _np.bincount(r.ravel())
#     radialprofile = tbin / nr
#     return radialprofile


# def radial_extraction(x, y, Bm, R, a):
#     index_x = x.shape[0] // 2
#     index_y = y.shape[1] // 2

#     radial_length = _np.sqrt(
#         (x[0, -1] - x[index_x, index_y]) ** 2 + (y[0, -1] - y[index_x, index_y]) ** 2
#     )
#     ravg = radial_profile(Bm.n, [index_x, index_y])
#     xpl = _np.linspace(0, radial_length, len(ravg))
#     ravg = ravg[xpl <= R - a]
#     xpl = xpl[xpl <= R - a]
#     return xpl, ravg