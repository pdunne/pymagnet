# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Plotting routines for calculating along symmetry lines of cubes, cuboids, and cylinders

"""
import matplotlib.pyplot as _plt
import numpy as _np
from ..utils._vector_structs import Point_Array1

from ..magnets import (
    Cylinder,
    Prism,
    magnetic_field_prism_1D,
    magnetic_field_cylinder_1D,
)

__all__ = ["plot_1D_field"]


def plot_1D_field(magnet, unit="mm", **kwargs):
    """Calculates and plots the magnetic field along the central symmetry axis
    of a cylinder or cuboid magnet, assuming the magnetic field is collinear

    Args:
        magnet (Magnet3D): Must be a Magnet3D type of magnet, either Prism, Cube,or Cylinder.

    Kwargs:
        num_points (int): Number of points to calculate. Defaults to 101.

    Returns:
        tuple: Point_Array1, Field1: point array struct containing z and the unit (e/g. 'mm'), vector array containing Bz and the field unit (e.g. 'T').
    """

    num_points = kwargs.pop("num_points", 101)
    return_data = kwargs.pop("return_data", False)
    points = Point_Array1(_np.zeros(num_points), unit=unit)

    if issubclass(magnet.__class__, Cylinder):
        mag_boundary = magnet.length / 2
        points.z = _np.linspace(
            -2 * magnet.length + magnet.center[2],
            2 * magnet.length + magnet.center[2],
            num_points,
        )
        field = magnetic_field_cylinder_1D(magnet, points.z)

        # if true, apply NaNs to inside the magnet
        if magnet._mask_magnet:
            mask = _generate_mask_1D(mag_boundary, magnet.center[2], points.z)
            field.z[mask] = _np.NaN

    elif issubclass(magnet.__class__, Prism):
        mag_boundary = magnet.height / 2
        points.z = _np.linspace(
            -2 * magnet.height + magnet.center[2],
            2 * magnet.height + magnet.center[2],
            num_points,
        )
        field = magnetic_field_prism_1D(magnet, points.z)

        # if true, apply NaNs to inside the magnet
        if magnet._mask_magnet:
            mask = _generate_mask_1D(mag_boundary, magnet.center[2], points.z)
            field.z[mask] = _np.NaN

    else:
        print("Error")
        return None

    fig, ax = _plt.subplots(figsize=(8, 8))
    unit_length = "(" + points.unit + ")"
    field_unit = "(" + field.unit + ")"
    _plt.xlabel(r"$z$ " + unit_length)
    _plt.ylabel(r"$B_z$ " + field_unit)
    _plt.plot(points.z, field.z)
    _plt.axvline(x=-mag_boundary + magnet.center[2], c="blue", ls="--")
    _plt.axvline(x=mag_boundary + magnet.center[2], c="red", ls="--")
    _plt.axvline(x=0.0, c="k", ls="-")
    _plt.show()

    if return_data:
        return points, field


def _generate_mask_1D(mag_boundary, zc, z):
    """Generates mask for points inside magnet

    Args:
        mag_boundary (float): half the height/length of a magnet
        zc (float): center of the magnet
        z (ndarray): z-coordinates

    Returns:
        ndarray: boolean mask array
    """
    zn = -mag_boundary + zc
    zp = mag_boundary + zc
    mask = _np.logical_and(z > zn, z < zp)

    return mask