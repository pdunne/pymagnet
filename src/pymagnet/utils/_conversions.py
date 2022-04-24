# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Routines for converting between coordinate systems, between 2D cartesian and polar,
as well as 3D cartesian, cylindrical, and spherical.

"""
from typing import Optional

import numpy as _np

from .typeconf import Numeric64, Tuple2_64, Tuple3_64


def cart2pol(x: Numeric64, y: Numeric64) -> Tuple2_64:
    """Converts from cartesian to polar coordinates

    Args:
        x (ndarray): x coordinates
        y (ndarray): y coordinates

    Returns:
        tuple: rho, phi
    """
    rho: Numeric64 = _np.sqrt(x**2 + y**2)
    phi: Numeric64 = _np.arctan2(y, x)
    return rho, phi


def pol2cart(rho: Numeric64, phi: Numeric64) -> Tuple2_64:
    """Converts from polar to cartesian coordinates

    Args:
        rho (ndarray): radial coordinates
        phi (ndarray): azimuthal coordinates

    Returns:
        tuple: x,y
    """
    x: Numeric64 = rho * _np.cos(phi)
    y: Numeric64 = rho * _np.sin(phi)
    return x, y


def vector_pol2cart(Brho: Numeric64, Bphi: Numeric64, phi: Numeric64) -> Tuple2_64:
    """Converts Vectors from polar to cartesian coordinates

    Args:
        Brho (ndarray): radial vector component
        Bphi (ndarray): azimuthal vector component
        phi (ndarray): azimuthal coordinates

    Returns:
        tuple: Bx, By
    """
    Bx: Numeric64 = Brho * _np.cos(phi) - Bphi * _np.sin(phi)
    By: Numeric64 = Brho * _np.sin(phi) + Bphi * _np.cos(phi)
    return Bx, By


def cart2sph(x: Numeric64, y: Numeric64, z: Numeric64) -> Tuple3_64:
    """Converts from cartesian to spherical coordinates

    Args:
        x (ndarray): x coordinates
        y (ndarray): y coordinates
        z (ndarray): z coordinates

    Returns:
        tuple: r, theta, phi
    """
    r: Numeric64 = _np.sqrt(x**2 + y**2 + z**2)
    phi: Numeric64 = _np.arctan2(y, x)

    # Hide the warning for situtations where there is a divide by zero.
    # This returns a NaN in the array, which is ignored for plotting.
    with _np.errstate(divide="ignore", invalid="ignore"):
        theta: Numeric64 = _np.arccos(z / r)
    return (r, theta, phi)


def sph2cart(r: Numeric64, theta: Numeric64, phi: Numeric64) -> Tuple3_64:
    """Converts from spherical to cartesian coordinates

    Args:
        r (ndarray): radial coordinates
        theta (ndarray): azimuthal angles
        phi (ndarray): polar angle

    Returns:
        tuple: x,y,z
    """
    x: Numeric64 = r * _np.sin(theta) * _np.cos(phi)
    y: Numeric64 = r * _np.sin(theta) * _np.sin(phi)
    z: Numeric64 = r * _np.cos(theta)
    return x, y, z


def vector_sph2cart(
    Br: Numeric64, Btheta: Numeric64, Bphi: Numeric64, theta: Numeric64, phi: Numeric64
) -> Tuple3_64:
    """Converts Vectors from spherical to cartesian coordinates

    Args:
        Br (ndarray): radial vector component
        Btheta (ndarray): polar vector component
        Bphi (ndarray): azimuthal vector component
        theta (ndarray): azimuthal angles
        phi (ndarray): polar angle

    Returns:
        tuple: Bx,By,Bz
    """
    Bx: Numeric64 = (
        Br * _np.sin(theta) * _np.cos(phi)
        + Btheta * _np.cos(theta) * _np.cos(phi)
        - Bphi * _np.sin(phi)
    )

    By: Numeric64 = (
        Br * _np.sin(theta) * _np.sin(phi)
        + Btheta * _np.cos(theta) * _np.sin(phi)
        + Bphi * _np.cos(phi)
    )

    Bz: Numeric64 = Br * _np.cos(theta) - Btheta * _np.sin(theta)
    return Bx, By, Bz


def sphere_sph2cart(
    Br: Numeric64, Btheta: Numeric64, theta: Numeric64, phi: Numeric64
) -> Tuple3_64:
    """Converts magnetic field of a sphere from spherical to cartesian coordinates

    Args:
        Br (ndarray): radial vector component
        Btheta (ndarray): polar vector component
        theta (ndarray): azimuthal angles
        phi (ndarray): polar angle

    Returns:
        tuple: Bx,By,Bz
    """
    Bx: Numeric64 = Br * _np.sin(theta) * _np.cos(phi) + Btheta * _np.cos(
        theta
    ) * _np.cos(phi)

    By: Numeric64 = Br * _np.sin(theta) * _np.sin(phi) + Btheta * _np.cos(
        theta
    ) * _np.sin(phi)

    Bz: Numeric64 = Br * _np.cos(theta) - Btheta * _np.sin(theta)
    return Bx, By, Bz


def get_unit_value_meter(unit: str) -> Optional[float]:
    """Returns a queried metre unit as a number
    Example:
        factor = get_unit_value_meter('cm')
        print(f"factor for 'cm' is {factor}")

    Args:
        unit (string): SI length unit

    Returns:
        float: SI prefix factor
    """
    si_prefixes = {
        "Ym": 1e24,
        "Zm": 1e21,
        "Em": 1e18,
        "Pm": 1e15,
        "Tm": 1e12,
        "Gm": 1e9,
        "Mm": 1e6,
        "km": 1e3,
        "hm": 1e2,
        "dam": 1e1,
        "m": 1.0,
        "dm": 1e-1,
        "cm": 1e-2,
        "mm": 1e-3,
        "µm": 1e-6,
        "um": 1e-6,
        "nm": 1e-9,
        "Ang": 1e-10,
        "pm": 1e-12,
        "fm": 1e-15,
        "am": 1e-18,
        "zm": 1e-21,
        "ym": 1e-24,
    }

    return si_prefixes.get(unit, None)


def get_unit_value_tesla(unit: str) -> Optional[float]:
    """Returns a queried magnetic flux density unit as a number
    Example:
        factor = get_unit_value_meter('mT')
        print(f"factor for 'mT' is {factor}")

    Args:
        unit (string): SI length unit

    Returns:
        float: SI prefix factor
    """

    si_prefixes = {
        "YT": 1e24,
        "ZT": 1e21,
        "ET": 1e18,
        "PT": 1e15,
        "TT": 1e12,
        "GT": 1e9,
        "MT": 1e6,
        "kT": 1e3,
        "hT": 1e2,
        "daT": 1e1,
        "T": 1,
        "dT": 1e-1,
        "cT": 1e-2,
        "mT": 1e-3,
        "µT": 1e-6,
        "uT": 1e-6,
        "nT": 1e-9,
        "pT": 1e-12,
        "fT": 1e-15,
        "aT": 1e-18,
        "zT": 1e-21,
        "yT": 1e-24,
    }

    return si_prefixes.get(unit, None)
