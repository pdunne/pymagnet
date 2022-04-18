# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Unit and Integration Tests

TODO:
    - All tests need to be updated and expanded to reflect the new API.
"""
from numpy import round as rnd
from numpy.linalg import norm

from pymagnet import get_field_2D, get_field_3D, magnets, reset
from pymagnet.utils import Point2, Point3


def test_create_cylinder():
    """Creating a cylinder magnet"""
    R, L, Jr = 1.0, 1.0, 1.0
    center = (0.0, 0.0, 0.0)
    magnet = magnets.Cylinder(R, L, Jr, center=center)
    assert isinstance(magnet, magnets.Cylinder)


def test_create_prism():
    """Creating a cuboid magnet"""
    a, b, c, Jr = 1.0, 1.0, 1.0, 1.0
    center = (0.0, 0.0, 0.0)
    magnet = magnets.Prism(a, b, c, Jr, center=center)
    assert isinstance(magnet, magnets.Prism)


def test_1D_cylinder_surface_max():
    """Cylinder: Central magnetic field on surface"""
    R, L = 1.0, 2.0
    Jr = 1.0
    zc = -L / 2
    m_cyl = magnets.Cylinder(radius=R, length=L, Jr=Jr, center=(0.0, 0.0, zc))
    Bz = magnets.magnetic_field_cylinder_1D(m_cyl, 0.0).z
    assert rnd(Bz, 5) == 0.44721


def test_1D_cylinder_center_value():
    """Cylinder: Magnetic field at magnet center"""
    R, L = 1.0, 2.0
    Jr = 1.0
    zc = -L / 2
    m_cyl = magnets.Cylinder(radius=R, length=L, Jr=Jr, center=(0.0, 0.0, zc))

    Bz = magnets.magnetic_field_cylinder_1D(m_cyl, zc).z
    assert rnd(Bz, 5) == 0.70711


def test_1D_prism_surface_max():
    """Cuboid: Central magnetic field on surface"""
    width = 1.0
    depth = 2.0
    height = 3.0
    Jr = 1.0
    zc = -height / 2
    m_quad = magnets.Prism(
        width=width, depth=depth, height=height, Jr=Jr, center=(0.0, 0.0, zc)
    )
    Bz = magnets.magnetic_field_prism_1D(m_quad, 0).z
    assert rnd(Bz, 5) == 0.48344


def test_1D_prism_center_value():
    """Cuboid: Magnetic field at magnet center\n"""
    width = 1.0
    depth = 2.0
    height = 3.0
    Jr = 1.0
    zc = -height / 2
    m_quad = magnets.Prism(
        width=width, depth=depth, height=height, Jr=Jr, center=(0.0, 0.0, zc)
    )
    Bz = magnets.magnetic_field_prism_1D(m_quad, zc).z
    assert rnd(Bz, 5) == 0.88775


def test_2D_two_magnets():
    """Test 2 Rectangular Magnets\n"""

    reset()
    width = 20e-3
    height = 40e-3
    center = (-0.75 * width, 0)
    _ = magnets.Rectangle(
        width=width, height=height, Jr=1.0, center=center, theta=45.0
    )
    center = (0.75 * width, 0)
    _ = magnets.Rectangle(
        width=width, height=height, Jr=1.0, center=center, theta=-45.0
    )

    x, y = 0.0, height / 2
    B = get_field_2D(Point2(x, y))
    result = rnd(B.n[0], 5)
    assert result == 0.13822


def test_3D_four_cubes():
    """Test 4 Cubes in Pseudo-Quadrupolar Arrangement"""

    reset()
    width = 10e-3
    a = width / 2
    hGap = a
    theta, phi = 0.0, 90.0

    # Add top left magnet
    _ = magnets.Cube(
        width=width, Jr=1.0, center=(-a - hGap, 0, a), theta=theta, phi=phi
    )
    # Add bottom left magnet
    _ = magnets.Cube(
        a=a, Jr=-1.0, center=(-a - hGap, 0, -a), theta=theta, phi=phi
    )
    # Add top right magnet
    _ = magnets.Cube(a=a, Jr=1.0, center=(a + hGap, 0, a), theta=theta, phi=phi)
    # Add bottom right magnet
    _ = magnets.Cube(
        a=a, Jr=-1.0, center=(a + hGap, 0, -a), theta=theta, phi=phi
    )

    x, y, z = a / 2, a, a / 3
    result = rnd(get_field_3D(Point3(x, y, z)).n[0], 5)
    assert result == 0.62385


def test_solenoid_off_axis():
    """Off-axis test of solenoid magnetic field"""

    R = 5e-3
    L = 20e-3
    m_cyl = magnets.Cylinder(radius=R, length=L, Jr=1.0, center=(0.0, 0.0, 0))

    result = round(
        norm(m_cyl._calcB_cyl(m_cyl.radius * 0.4, m_cyl.length * 0.2)), 5
    )
    assert result == 0.86351
