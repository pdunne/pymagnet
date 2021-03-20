# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Routines for converting between coordinate systems
"""

import numpy as _np


def cart2pol(x, y):
    """Converts from cartesian to polar coordinates

    Args:
        x (float/array): x coordinates
        y (float/array): y coordinates

    Returns:
        tuple: rho, phi
    """
    rho = _np.sqrt(x ** 2 + y ** 2)
    phi = _np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    """Converts from polar to cartesian coordinates

    Args:
        rho (float/array): radial coordinates
        phi (float/array): azimuthal coordinates

    Returns:
        tuple: x,y
    """
    x = rho * _np.cos(phi)
    y = rho * _np.sin(phi)
    return (x, y)


def vector_pol2cart(Brho, Bphi, phi):
    """Converts Vectors from polar to cartesian coordinates

    Args:
        Brho (float/array): radial vector component
        Bphi (float/array): azimuthal vector component
        phi (float/array): azimuthal coordinates

    Returns:
        tuple: x,y
    """
    Bx = Brho * _np.cos(phi) - Bphi * _np.sin(phi)
    By = Brho * _np.sin(phi) + Bphi * _np.cos(phi)
    return Bx, By


def cart2sph(x, y, z):
    """Converts from cartesian to spherical coordinates

    Args:
        x (float/array): x coordinates
        y (float/array): y coordinates
        z (float/array): z coordinates

    Returns:
        tuple: r, theta, phi
    """
    r = _np.sqrt(x ** 2 + y ** 2 + z ** 2)
    phi = _np.arctan2(y, x)

    # Hide the warning for situtations where there is a divide by zero.
    # This returns a NaN in the array, which is ignored for plotting.
    with _np.errstate(divide="ignore", invalid="ignore"):
        theta = _np.arccos(z / r)
    return (r, theta, phi)


def sph2cart(r, theta, phi):
    """Converts from spherical to cartesian coordinates

    Args:
        r (float/array): radial coordinates
        theta (float/array): azimuthal angles
        phi (float/array): polar angle

    Returns:
        tuple: x,y,z
    """
    x = r * _np.sin(theta) * _np.cos(phi)
    y = r * _np.sin(theta) * _np.sin(phi)
    z = r * _np.cos(theta)
    return x, y, z


def vector_sph2cart(Br, Btheta, Bphi, theta, phi):
    """Converts Vectors from spherical to cartesian coordinates

    Args:
        Br (float/array): radial vector component
        Btheta (float/array): polar vector component
        Bphi (float/array): azimuthal vector component
        theta (float/array): azimuthal angles
        phi (float/array): polar angle

    Returns:
        tuple: Bx,By,Bz
    """
    Bx = (
        Br * _np.sin(theta) * _np.cos(phi)
        + Btheta * _np.cos(theta) * _np.cos(phi)
        - Bphi * _np.sin(phi)
    )

    By = (
        Br * _np.sin(theta) * _np.sin(phi)
        + Btheta * _np.cos(theta) * _np.sin(phi)
        + Bphi * _np.cos(phi)
    )

    Bz = Br * _np.cos(theta) - Btheta * _np.sin(theta)
    return Bx, By, Bz


def sphere_sph2cart(Br, Btheta, theta, phi):
    """Converts magnetic field of a sphere from spherical to cartesian coordinates

    Args:
        Br (float/array): radial vector component
        Btheta (float/array): polar vector component
        theta (float/array): azimuthal angles
        phi (float/array): polar angle

    Returns:
        tuple: Bx,By,Bz
    """
    Bx = Br * _np.sin(theta) * _np.cos(phi) + Btheta * _np.cos(theta) * _np.cos(phi)

    By = Br * _np.sin(theta) * _np.sin(phi) + Btheta * _np.cos(theta) * _np.sin(phi)

    Bz = Br * _np.cos(theta) - Btheta * _np.sin(theta)
    return Bx, By, Bz