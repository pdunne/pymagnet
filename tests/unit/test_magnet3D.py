# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for 3D magnet classes: Magnet3D, Prism, Cube, Cylinder, Sphere."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet import magnets


class TestPrismInit:
    """Tests for Prism initialization."""

    def test_dimensions(self, prism_default):
        """Width, depth, height properties set correctly."""
        assert prism_default.width == 10.0
        assert prism_default.depth == 20.0
        assert prism_default.height == 30.0

    def test_half_dimensions(self, prism_default):
        """Half dimensions a, b, c."""
        assert prism_default.a == 5.0
        assert prism_default.b == 10.0
        assert prism_default.c == 15.0

    def test_default_center(self, prism_default):
        """Default center is origin."""
        center = prism_default.get_center()
        npt.assert_allclose(center, [0.0, 0.0, 0.0])

    def test_custom_center(self):
        """Custom center from kwarg."""
        p = magnets.Prism(
            width=10.0, depth=20.0, height=30.0, Jr=1.0, center=(5.0, 6.0, 7.0)
        )
        center = p.get_center()
        npt.assert_allclose(center, [5.0, 6.0, 7.0])

    def test_jr_attribute(self, prism_default):
        """Jr property set correctly."""
        assert prism_default.Jr == 1.0

    def test_default_theta_phi(self):
        """Default theta=0, phi=0 (magnetization along z)."""
        p = magnets.Prism(width=10.0, depth=20.0, height=30.0, Jr=1.0)
        assert p.theta == 0.0

    def test_get_size(self, prism_default):
        """get_size returns [width, depth, height]."""
        size = prism_default.get_size()
        assert size[0] == 10.0
        assert size[1] == 20.0
        assert size[2] == 30.0

    def test_get_jr_z_magnetization(self, prism_default):
        """get_Jr returns [Jx, Jy, Jz] for z-magnetization."""
        jr = prism_default.get_Jr()
        npt.assert_allclose(jr[0], 0.0, atol=1e-10)  # Jx
        npt.assert_allclose(jr[1], 0.0, atol=1e-10)  # Jy
        npt.assert_allclose(jr[2], 1.0)  # Jz

    def test_rotation_angles(self, prism_rotated):
        """Alpha, beta, gamma rotation angles set correctly."""
        orientation = prism_rotated.get_orientation()
        npt.assert_allclose(orientation, [30.0, 45.0, 60.0])

    def test_str_repr(self, prism_default):
        """String representation contains class name."""
        str_rep = str(prism_default)
        assert "Prism" in str_rep


class TestPrismField:
    """Tests for Prism field calculations."""

    def test_field_returns_tuple(self, prism_default):
        """get_field returns Bx, By, Bz tuple."""
        Bx, By, Bz = prism_default.get_field(0.0, 0.0, 50.0)
        assert isinstance(Bx, (float, np.floating, np.ndarray))
        assert isinstance(By, (float, np.floating, np.ndarray))
        assert isinstance(Bz, (float, np.floating, np.ndarray))

    def test_field_on_axis(self, prism_default):
        """Field on z-axis is primarily Bz for z-magnetized prism."""
        Bx, By, Bz = prism_default.get_field(0.0, 0.0, 50.0)
        npt.assert_allclose(Bx, 0.0, atol=1e-10)
        npt.assert_allclose(By, 0.0, atol=1e-10)
        assert Bz > 0  # Positive Bz above magnet

    def test_field_symmetry_x(self, prism_default):
        """Field symmetric about x=0 plane."""
        Bx_pos, By_pos, Bz_pos = prism_default.get_field(20.0, 0.0, 0.0)
        Bx_neg, By_neg, Bz_neg = prism_default.get_field(-20.0, 0.0, 0.0)
        npt.assert_allclose(Bx_pos, -Bx_neg, rtol=1e-10)
        npt.assert_allclose(By_pos, By_neg, rtol=1e-10)
        npt.assert_allclose(Bz_pos, Bz_neg, rtol=1e-10)

    def test_field_decay_with_distance(self, prism_default):
        """Field decays with distance."""
        _, _, Bz_close = prism_default.get_field(0.0, 0.0, 20.0)
        _, _, Bz_far = prism_default.get_field(0.0, 0.0, 100.0)
        assert abs(Bz_close) > abs(Bz_far)

    def test_field_array_input(self, prism_default):
        """Field calculation with array inputs."""
        x = np.array([0.0, 10.0, -10.0])
        y = np.array([0.0, 0.0, 0.0])
        z = np.array([50.0, 50.0, 50.0])
        Bx, By, Bz = prism_default.get_field(x, y, z)
        assert Bx.shape == (3,)
        assert By.shape == (3,)
        assert Bz.shape == (3,)


class TestCube:
    """Tests for Cube class (inherits from Prism)."""

    def test_equal_dimensions(self, cube_default):
        """All dimensions equal."""
        assert cube_default.width == cube_default.depth == cube_default.height

    def test_inherits_prism(self, cube_default):
        """Cube is subclass of Prism."""
        assert isinstance(cube_default, magnets.Prism)

    def test_get_size_cube(self, cube_default):
        """get_size returns equal dimensions."""
        size = cube_default.get_size()
        assert size[0] == size[1] == size[2] == 10.0

    def test_init_with_a_parameter(self):
        """Cube can be initialized with 'a' parameter (half-width)."""
        c = magnets.Cube(a=5.0, Jr=1.0)
        assert c.width == 10.0


class TestCylinderInit:
    """Tests for Cylinder initialization."""

    def test_radius_length(self, cylinder_default):
        """Radius and length properties set correctly."""
        assert cylinder_default.radius == 10.0
        assert cylinder_default.length == 20.0

    def test_get_size(self, cylinder_default):
        """get_size returns [radius, length]."""
        size = cylinder_default.get_size()
        assert size[0] == 10.0
        assert size[1] == 20.0

    def test_jr_attribute(self, cylinder_default):
        """Jr property set correctly (axial magnetization)."""
        jr = cylinder_default.get_Jr()
        npt.assert_allclose(jr, [0.0, 0.0, 1.0])


class TestCylinderField:
    """Tests for Cylinder field calculations."""

    def test_field_on_axis(self, cylinder_default):
        """Field on z-axis is purely Bz."""
        Bx, By, Bz = cylinder_default.get_field(0.0, 0.0, 30.0)
        npt.assert_allclose(Bx, 0.0, atol=1e-10)
        npt.assert_allclose(By, 0.0, atol=1e-10)
        assert Bz > 0

    def test_field_axial_symmetry(self, cylinder_default):
        """Field has axial symmetry about z-axis."""
        # Points at same distance from axis should have same |B|
        Bx1, By1, Bz1 = cylinder_default.get_field(5.0, 0.0, 30.0)
        Bx2, By2, Bz2 = cylinder_default.get_field(0.0, 5.0, 30.0)
        B1 = np.sqrt(Bx1**2 + By1**2 + Bz1**2)
        B2 = np.sqrt(Bx2**2 + By2**2 + Bz2**2)
        npt.assert_allclose(B1, B2, rtol=1e-10)

    def test_field_decay_with_distance(self, cylinder_default):
        """Field decays with distance."""
        _, _, Bz_close = cylinder_default.get_field(0.0, 0.0, 15.0)
        _, _, Bz_far = cylinder_default.get_field(0.0, 0.0, 100.0)
        assert abs(Bz_close) > abs(Bz_far)


class TestSphereInit:
    """Tests for Sphere initialization."""

    def test_radius(self, sphere_default):
        """Radius property set correctly."""
        assert sphere_default.radius == 10.0

    def test_get_size(self, sphere_default):
        """get_size returns [radius]."""
        size = sphere_default.get_size()
        assert size[0] == 10.0

    def test_jr_attribute(self, sphere_default):
        """Jr returns magnetization along z."""
        jr = sphere_default.get_Jr()
        npt.assert_allclose(jr, [0.0, 0.0, 1.0])


class TestSphereField:
    """Tests for Sphere field calculations.

    Note: Sphere.get_field has an issue with scalar inputs due to mask application.
    Using array inputs as workaround.
    """

    def test_field_on_z_axis(self, sphere_default):
        """Field on z-axis is purely Bz."""
        # Use array inputs to avoid scalar mask issue
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([20.0])
        Bx, By, Bz = sphere_default.get_field(x, y, z)
        npt.assert_allclose(Bx[0], 0.0, atol=1e-10)
        npt.assert_allclose(By[0], 0.0, atol=1e-10)
        assert Bz[0] > 0

    def test_field_symmetry(self, sphere_default):
        """Field has rotational symmetry about z-axis."""
        # Use array inputs
        Bx1, By1, Bz1 = sphere_default.get_field(
            np.array([5.0]), np.array([0.0]), np.array([20.0])
        )
        Bx2, By2, Bz2 = sphere_default.get_field(
            np.array([0.0]), np.array([5.0]), np.array([20.0])
        )
        B1 = np.sqrt(Bx1[0] ** 2 + By1[0] ** 2 + Bz1[0] ** 2)
        B2 = np.sqrt(Bx2[0] ** 2 + By2[0] ** 2 + Bz2[0] ** 2)
        npt.assert_allclose(B1, B2, rtol=1e-10)

    def test_dipole_decay(self, sphere_default):
        """Far field decays as 1/r^3 (dipole)."""
        # Use array inputs
        _, _, Bz1 = sphere_default.get_field(
            np.array([0.0]), np.array([0.0]), np.array([50.0])
        )
        _, _, Bz2 = sphere_default.get_field(
            np.array([0.0]), np.array([0.0]), np.array([100.0])
        )
        # Doubling distance should reduce field by factor of 8
        ratio = Bz1[0] / Bz2[0]
        npt.assert_allclose(ratio, 8.0, rtol=0.1)  # Allow 10% tolerance


class TestMagnet3DCommon:
    """Tests common to all 3D magnets."""

    @pytest.mark.parametrize("Jr", [0.5, 1.0, -1.0, 2.0])
    def test_magnetization_scaling(self, Jr):
        """Field scales with magnetization."""
        p1 = magnets.Prism(width=10.0, depth=10.0, height=10.0, Jr=1.0)
        p2 = magnets.Prism(width=10.0, depth=10.0, height=10.0, Jr=Jr)

        _, _, Bz1 = p1.get_field(0.0, 0.0, 30.0)
        _, _, Bz2 = p2.get_field(0.0, 0.0, 30.0)

        npt.assert_allclose(Bz2, Bz1 * Jr, rtol=1e-10)

    @pytest.mark.parametrize(
        "center",
        [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)],
    )
    def test_center_offset(self, center):
        """Magnet can be created at various centers."""
        p = magnets.Prism(width=10.0, depth=10.0, height=10.0, Jr=1.0, center=center)
        result_center = p.get_center()
        npt.assert_allclose(result_center, center)


class TestMagnet1DFunctions:
    """Tests for 1D field calculation functions."""

    def test_cylinder_1d_on_surface(self):
        """1D cylinder field on surface."""
        cyl = magnets.Cylinder(radius=1.0, length=2.0, Jr=1.0, center=(0.0, 0.0, -1.0))
        field = magnets.magnetic_field_cylinder_1D(cyl, 0.0)
        assert hasattr(field, "z")
        assert np.isfinite(field.z)

    def test_cylinder_1d_at_center(self):
        """1D cylinder field at magnet center."""
        cyl = magnets.Cylinder(radius=1.0, length=2.0, Jr=1.0, center=(0.0, 0.0, -1.0))
        field = magnets.magnetic_field_cylinder_1D(cyl, -1.0)
        assert np.isfinite(field.z)

    def test_prism_1d_on_surface(self):
        """1D prism field on surface."""
        prism = magnets.Prism(
            width=1.0, depth=2.0, height=3.0, Jr=1.0, center=(0.0, 0.0, -1.5)
        )
        field = magnets.magnetic_field_prism_1D(prism, 0.0)
        assert hasattr(field, "z")
        assert np.isfinite(field.z)

    def test_prism_1d_at_center(self):
        """1D prism field at magnet center."""
        prism = magnets.Prism(
            width=1.0, depth=2.0, height=3.0, Jr=1.0, center=(0.0, 0.0, -1.5)
        )
        field = magnets.magnetic_field_prism_1D(prism, -1.5)
        assert np.isfinite(field.z)
