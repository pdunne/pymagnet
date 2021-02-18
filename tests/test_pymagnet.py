# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne

import pkg_resources
from pathlib import Path
import unittest
import os
import sys
import pymagnet as pm  # Magnetic library
import numpy as np  # we use numpy for handing vectors and matrices
import matplotlib.pyplot as plt  # Matplotlib is used to generate our plots


_REQUIREMENTS_PATH = ('tests/requirements.txt')

class MagnetTest(unittest.TestCase):

    def test_python_major_version(self):
        """Requires at least Python 3.6"""
        err_msg = "Requires at least Python 3.6\nYou are using Python {}.{}.".format(
            sys.version_info.major, sys.version_info.minor)
        self.assertEqual(sys.version_info.major == 3 and sys.version_info.minor >= 6,
                         True, msg=err_msg)

    def test_requirements(self):
        """Test that each required package is available."""
        # Ref: https://stackoverflow.com/a/45474387/
        with open(_REQUIREMENTS_PATH) as f:
            requirements = pkg_resources.parse_requirements(f)
            for requirement in requirements:
                requirement = str(requirement)
                with self.subTest(requirement=requirement):
                    pkg_resources.require(requirement)

    def test_create_cylinder(self):
        """Creating a cylinder magnet"""
        R, L, Jr = 1.0, 1.0, 1.0
        center = (0.0, 0.0, 0.0)
        self.assertEqual(pm.magnets.Cylinder(R, L, Jr, center=center).mag_type,
                         'Cylinder')

    def test_create_prism(self):
        """Creating a cuboid magnet"""
        a, b, c, Jr = 1.0, 1.0, 1.0, 1.0
        center = (0.0, 0.0, 0.0)
        self.assertEqual(pm.magnets.Prism(a, b, c, Jr, center=center).mag_type,
                         'Prism')

    def test_create_cube(self):
        """Creating a cube magnet"""
        a, Jr = 1.0, 1.0
        center = (0.0, 0.0, 0.0)
        self.assertEqual(pm.magnets.Cube(a, Jr, center=center).mag_type,
                         'Cube')

    def test_create_rectangle(self):
        """Creating a 2D rectangular magnet"""
        a, b, Jr = 1.0, 1.0, 1.0
        center = (0.0, 0.0)
        self.assertEqual(pm.magnets.Rectangle(a, b, Jr, center=center).mag_type,
                         'Rectangle')

    def test_create_square(self):
        """Creating a 2D square magnet"""
        a, Jr = 1.0, 1.0,
        center = (0.0, 0.0)
        self.assertEqual(pm.magnets.Square(a, Jr, center=center).mag_type,
                         'Square')

    def test_cylinder_class_instance(self):
        """Creating another Cylinder"""
        R, L, Jr = 1.0, 1.0, 1.0
        center = (0.0, 0.0, 0.0)
        magnet = pm.magnets.Cylinder(radius=R, length=L, Jr=Jr, center=center)
        self.assertIsInstance(magnet, pm.magnets.Cylinder)

    def test_1D_cylinder_surface_max(self):
        """Cylinder: Central magnetic field on surface"""
        R, L = 1.0, 2.0
        Jr = 1.0
        zc = -L/2
        m_cyl = pm.magnets.Cylinder(
            radius=R, length=L, Jr=Jr, center=(0.0, 0.0, zc))
        Bz = pm.magnets.magnetic_field_cylinder_1D(m_cyl, 0.0)
        self.assertEqual(round(Bz, 5), 0.44721)

    def test_1D_cylinder_center_value(self):
        """Cylinder: Magnetic field at magnet center"""
        R, L = 1.0, 2.0
        Jr = 1.0
        zc = -L/2
        m_cyl = pm.magnets.Cylinder(
            radius=R, length=L, Jr=Jr, center=(0.0, 0.0, zc))

        Bz = pm.magnets.magnetic_field_cylinder_1D(m_cyl, zc)
        self.assertEqual(round(Bz, 5), 0.70711)

    def test_1D_prism_surface_max(self):
        """Cuboid: Central magnetic field on surface"""
        width = 1.0
        depth = 2.0
        height = 3.0
        Jr = 1.0
        zc = -height / 2
        m_quad = pm.magnets.Prism(width=width, depth=depth, height=height,
                                    Jr=Jr, center=(0.0, 0.0, zc))
        Bz = pm.magnets.magnetic_field_prism_1D(m_quad, 0)
        self.assertEqual(round(Bz, 5), 0.48344)

    def test_1D_prism_center_value(self):
        """Cuboid: Magnetic field at magnet center\n"""
        width = 1.0
        depth = 2.0
        height = 3.0
        Jr = 1.0
        zc = -height / 2
        m_quad = pm.magnets.Prism(width=width, depth=depth, height=height,
                                   Jr=Jr, center=(0.0, 0.0, zc))
        Bz = pm.magnets.magnetic_field_prism_1D(m_quad, zc)
        self.assertEqual(round(Bz, 5), 0.88775)

    def test_2D_two_magnets(self):
        """Test 2 Rectangular Magnets\n"""

        pm.reset_magnets()
        width = 20e-3
        height = 40e-3
        center = (-0.75*width, 0)
        mag1 = pm.magnets.Rectangle(
            width=width, height=height, Jr=1.0, center=center, theta=45.0)
        center = (0.75*width, 0)
        mag2 = pm.magnets.Rectangle(
            width=width, height=height, Jr=1.0, center=center, theta=-45.0)

        x, y = 0.0, height/2
        B = pm.magnets.B_calc_2D(x, y)
        result = round(B.n[0], 5)
        self.assertEqual(result, 0.42462)

    def test_3D_four_cubes(self):
        """Test 4 Cubes in Pseudo-Quadrupolar Arrangement"""

        pm.reset_magnets()
        width = 10e-3
        a = width/2
        hGap = a
        theta, phi = 0.0, 90.0

        # Add top left magnet
        _ = pm.magnets.Cube(width=width, Jr=1.0, center=(-a - hGap,
                                                 0, a), theta=theta, phi=phi)
        # Add bottom left magnet
        _ = pm.magnets.Cube(
            a=a, Jr=-1.0, center=(-a - hGap, 0, -a), theta=theta, phi=phi)
        # Add top right magnet
        _ = pm.magnets.Cube(a=a, Jr=1.0, center=(
            a + hGap, 0, a), theta=theta, phi=phi)
        # Add bottom right magnet
        _ = pm.magnets.Cube(a=a, Jr=-1.0, center=(a + hGap,
                                                  0, -a), theta=theta, phi=phi)

        x, y, z = a/2, a, a/3
        result = round(pm.magnets.B_calc_3D(x, y, z).n[0], 5)
        self.assertEqual(result, 0.15011)
    
    def test_solenoid_off_axis(self):
        """Off-axis test of solenoid magnetic field"""

        R = 5e-3
        L = 20e-3
        m_cyl = pm.magnets.Cylinder(radius=R, length=L, Jr=1.0,
                            center=(0.0, 0.0, 0))

        result = round(np.linalg.norm(
            m_cyl._calcB_cyl(m_cyl.radius * .4, m_cyl.length * .2)), 5)
        self.assertEqual(result, 0.86351)


if __name__ == "__main__":
    unittest.main(argv=[''], verbosity=2,)
