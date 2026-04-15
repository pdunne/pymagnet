# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for numerical stability and edge cases."""

import numpy as np
import numpy.testing as npt

from pymagnet import magnets


class TestArccosDomain:
    """Tests for arccos domain safety."""

    def test_quaternion_handles_near_unity_cos(self):
        """Quaternion handles cos values near 1."""
        from pymagnet.utils._quaternion import q_angle_from_axis

        # Create a small but not too small rotation
        q = q_angle_from_axis(0.01, (0, 0, 1))
        # Should not raise ValueError for arccos
        angle, axis = q.get_axisangle()
        assert np.isfinite(angle)
        assert abs(angle - 0.01) < 0.001  # Close to 0.01 radians

    def test_quaternion_handles_negative_near_unity(self):
        """Quaternion handles cos values near -1."""
        from pymagnet.utils._quaternion import Quaternion

        # This is a 180-degree rotation
        q = Quaternion(0.0, 1.0, 0.0, 0.0)
        angle, axis = q.get_axisangle()
        assert np.isfinite(angle)


class TestDivisionByZero:
    """Tests for division by zero handling."""

    def test_rectangle_field_at_corner(self):
        """Rectangle field at corner is handled."""
        m = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        # Point at exact corner
        Bx, By = m.get_field(5.0, 10.0)
        # May be NaN, inf, or array at singularity, but shouldn't crash
        # Just check it returns something
        assert Bx is not None

    def test_prism_field_at_edge(self):
        """Prism field at edge is handled."""
        m = magnets.Prism(width=10.0, depth=10.0, height=20.0, Jr=1.0)
        # Point on edge
        Bx, By, Bz = m.get_field(5.0, 5.0, 10.0)
        # May be NaN at singularity, but shouldn't crash
        assert isinstance(Bx, (float, np.floating, np.ndarray))

    def test_cylinder_field_on_axis_at_pole(self):
        """Cylinder field on axis at pole is finite."""
        m = magnets.Cylinder(radius=5.0, length=10.0, Jr=1.0)
        # Exactly on axis at top surface
        Bx, By, Bz = m.get_field(0.0, 0.0, 5.0)
        # Should be finite (analytically defined)
        assert np.isfinite(Bx)
        assert np.isfinite(By)
        assert np.isfinite(Bz)


class TestDegenerateGeometry:
    """Tests for degenerate geometric cases."""

    def test_zero_magnetization(self):
        """Zero magnetization produces zero field."""
        m = magnets.Rectangle(width=10.0, height=20.0, Jr=0.0)
        Bx, By = m.get_field(20.0, 20.0)
        npt.assert_allclose(Bx, 0.0, atol=1e-15)
        npt.assert_allclose(By, 0.0, atol=1e-15)

    def test_negative_magnetization(self):
        """Negative magnetization works correctly."""
        m1 = magnets.Prism(width=10.0, depth=10.0, height=20.0, Jr=1.0)
        Bx1, By1, Bz1 = m1.get_field(0.0, 0.0, 30.0)

        from pymagnet import reset

        reset()

        m2 = magnets.Prism(width=10.0, depth=10.0, height=20.0, Jr=-1.0)
        Bx2, By2, Bz2 = m2.get_field(0.0, 0.0, 30.0)

        # Field should be negated
        npt.assert_allclose(Bx2, -Bx1, rtol=1e-10)
        npt.assert_allclose(By2, -By1, rtol=1e-10)
        npt.assert_allclose(Bz2, -Bz1, rtol=1e-10)


class TestLargeCoordinates:
    """Tests for numerical stability with large coordinates."""

    def test_field_at_large_distance(self):
        """Field at large distance is finite and small."""
        m = magnets.Prism(width=10.0, depth=10.0, height=20.0, Jr=1.0)
        Bx, By, Bz = m.get_field(1000.0, 0.0, 0.0)
        assert np.isfinite(Bx)
        assert np.isfinite(By)
        assert np.isfinite(Bz)
        # Field should be very small at large distance
        B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
        assert B_mag < 0.001  # Should be much smaller than 1

    def test_field_at_very_large_distance(self):
        """Field at very large distance doesn't overflow."""
        m = magnets.Cylinder(radius=5.0, length=10.0, Jr=1.0)
        Bx, By, Bz = m.get_field(10000.0, 0.0, 0.0)
        assert np.isfinite(Bx)
        assert np.isfinite(By)
        assert np.isfinite(Bz)


class TestSmallCoordinates:
    """Tests for numerical stability with small coordinates."""

    def test_field_very_close_to_surface(self):
        """Field very close to surface is finite."""
        m = magnets.Prism(width=10.0, depth=10.0, height=20.0, Jr=1.0)
        # Just above the surface
        Bx, By, Bz = m.get_field(0.0, 0.0, 10.001)
        assert np.isfinite(Bx)
        assert np.isfinite(By)
        assert np.isfinite(Bz)

    def test_small_magnet_dimensions(self):
        """Small magnet dimensions work correctly."""
        m = magnets.Prism(width=0.001, depth=0.001, height=0.001, Jr=1.0)
        Bx, By, Bz = m.get_field(0.0, 0.0, 0.01)
        assert np.isfinite(Bx)
        assert np.isfinite(By)
        assert np.isfinite(Bz)


class TestArrayOperations:
    """Tests for numerical stability with array operations."""

    def test_large_array_computation(self):
        """Large array field computation is stable."""
        m = magnets.Prism(width=10.0, depth=10.0, height=20.0, Jr=1.0)
        x = np.linspace(-50, 50, 50)
        y = np.zeros(50)
        z = np.full(50, 30.0)
        Bx, By, Bz = m.get_field(x, y, z)
        assert np.all(np.isfinite(Bx))
        assert np.all(np.isfinite(By))
        assert np.all(np.isfinite(Bz))

    def test_meshgrid_computation(self):
        """Meshgrid field computation is stable."""
        m = magnets.Prism(width=10.0, depth=10.0, height=20.0, Jr=1.0)
        x, y = np.mgrid[-30:30:10j, -30:30:10j]
        z = np.full_like(x, 30.0)
        Bx, By, Bz = m.get_field(x, y, z)
        # Check that most values are finite (some may be NaN at singularities)
        finite_fraction = np.sum(np.isfinite(Bz)) / Bz.size
        assert finite_fraction > 0.9


class TestQuaternionStability:
    """Tests for quaternion numerical stability."""

    def test_double_rotation_returns_original(self):
        """Double 180-degree rotation returns original."""
        from pymagnet.utils._quaternion import q_angle_from_axis

        # 180-degree rotation about z
        q = q_angle_from_axis(np.pi, (0, 0, 1))
        # Apply twice
        v = np.array([1.0, 0.0, 0.0])
        v_rot = q * v
        v_rot2 = q * v_rot
        npt.assert_allclose(v_rot2, v, atol=1e-10)

    def test_small_angle_rotation(self):
        """Small angle rotation is accurate."""
        from pymagnet.utils._quaternion import q_angle_from_axis

        q = q_angle_from_axis(1e-6, (0, 0, 1))
        v = np.array([1.0, 0.0, 0.0])
        v_rot = q * v
        # For small angle, v_rot ≈ v + (1e-6) * (0,0,1) × (1,0,0) = (1, 1e-6, 0)
        npt.assert_allclose(v_rot[0], 1.0, atol=1e-8)
        npt.assert_allclose(v_rot[1], 1e-6, atol=1e-8)
        npt.assert_allclose(v_rot[2], 0.0, atol=1e-10)


class TestTrigonometry3DStability:
    """Tests for trigonometry3D numerical stability."""

    def test_nearly_flat_triangle(self):
        """Nearly flat triangle is handled."""
        from pymagnet.utils._trigonometry3D import signed_area

        # Triangle with very small z deviation
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 1e-10]])
        area = signed_area(triangle)
        assert np.isfinite(area)

    def test_very_small_triangle(self):
        """Very small triangle is handled."""
        from pymagnet.utils._trigonometry3D import signed_area

        triangle = np.array([[0.0, 0.0, 0.0], [1e-8, 0.0, 0.0], [0.5e-8, 1e-8, 0.0]])
        area = signed_area(triangle)
        assert np.isfinite(area)
        assert area >= 0
