# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Plotting routines for calculating along symmetry lines of 3D magnets

"""
import matplotlib.pyplot as _plt
import numpy as _np
from pymagnet import magnets as _mag


def plot_1D_field(magnet, **kwargs):
    """Calculates and plots the magnetic field along the central symmetry axis
    of a cylinder or cuboid magnet, assuming the magnetic field is collinear

    Args:
        magnet (magnet object): Must be a Magnet_3D type of magnet, either Prism, Cube,
        or Cylinder.

    Kwargs:
        NP (int): Number of points to calculate. Defaults to 101.

    Returns:
        tuple: z,Bz: arrays of z the distance from the magnet surface, and Bz the
        magnetic field.
    """

    NP = kwargs.pop("NP", 101)
    return_data = kwargs.pop("return_data", False)

    if issubclass(magnet.__class__, _mag.Cylinder):
        mag_boundary = magnet.length / 2
        z = _np.linspace(
            -2 * magnet.length + magnet.zc, 2 * magnet.length + magnet.zc, NP
        )
        Bz = _mag.magnetic_field_cylinder_1D(magnet, z)

        # if true, apply NaNs to inside the magnet
        if magnet._mask_magnet:
            mask = _generate_mask_1D(mag_boundary, magnet.zc, z)
            Bz[mask] = _np.NaN

    elif issubclass(magnet.__class__, _mag.Prism):
        mag_boundary = magnet.height / 2
        z = _np.linspace(
            -2 * magnet.height + magnet.zc, 2 * magnet.height + magnet.zc, NP
        )
        Bz = _mag.magnetic_field_prism_1D(magnet, z)

        # if true, apply NaNs to inside the magnet
        if magnet._mask_magnet:
            mask = _generate_mask_1D(mag_boundary, magnet.zc, z)
            Bz[mask] = _np.NaN

    else:
        print("Error")
        return None

    _, _ = _plt.subplots()
    _plt.xlabel(r"$z$ (mm)")
    _plt.ylabel(r"$B_z$ (mT)")
    _plt.plot(z * 1e3, Bz * 1e3, label="Cube")
    _plt.axvline(x=-mag_boundary * 1e3 + magnet.zc * 1e3, c="blue", ls="--")
    _plt.axvline(x=mag_boundary * 1e3 + magnet.zc * 1e3, c="red", ls="--")
    _plt.axvline(x=0.0, c="k", ls="-")
    _plt.show()

    if return_data:
        return z, Bz


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