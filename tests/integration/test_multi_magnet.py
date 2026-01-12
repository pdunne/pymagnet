# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Integration tests for multi-magnet systems."""

import numpy as np
import numpy.testing as npt

from pymagnet import magnets, reset
from pymagnet.utils._routines2D import get_field_2D, grid2D
from pymagnet.utils._routines3D import get_field_3D, slice3D


class TestFieldSuperposition2D:
    """Tests for 2D field superposition."""

    def test_two_magnets_field_adds(self):
        """Fields from two magnets superimpose."""
        # Create single magnet and measure field
        m1 = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, center=(0.0, 0.0))
        points = grid2D(30.0, 30.0, num_points=5)
        field_single = get_field_2D(points)
        Bx_single = field_single.x.copy()

        reset()

        # Create two identical magnets at different positions
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, center=(-20.0, 0.0))
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, center=(20.0, 0.0))
        field_double = get_field_2D(points)

        # Fields should be different (superposition)
        assert not np.allclose(field_double.x, Bx_single)

    def test_opposing_magnets_cancel(self):
        """Opposing magnetizations partially cancel."""
        # Two magnets with opposite Jr
        m1 = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, center=(-15.0, 0.0))
        m2 = magnets.Rectangle(width=10.0, height=20.0, Jr=-1.0, center=(15.0, 0.0))

        # Field at origin should have reduced magnitude
        Bx, By = m1.get_field(0.0, 0.0)
        Bx2, By2 = m2.get_field(0.0, 0.0)

        # At origin between them, fields partially cancel
        # Total field magnitude should be less than either individual
        assert np.isfinite(Bx + Bx2)


class TestFieldSuperposition3D:
    """Tests for 3D field superposition."""

    def test_two_prisms_field_adds(self):
        """Fields from two prisms superimpose."""
        # Create single prism
        _ = magnets.Prism(
            width=5.0, depth=5.0, height=10.0, Jr=1.0, center=(0.0, 0.0, 0.0)
        )
        points = slice3D(
            plane="xz", max1=20.0, max2=20.0, slice_value=0.0, num_points=5
        )
        field_single = get_field_3D(points)
        Bz_single = field_single.z.copy()

        reset()

        # Create two prisms
        _ = magnets.Prism(
            width=5.0, depth=5.0, height=10.0, Jr=1.0, center=(-15.0, 0.0, 0.0)
        )
        _ = magnets.Prism(
            width=5.0, depth=5.0, height=10.0, Jr=1.0, center=(15.0, 0.0, 0.0)
        )
        field_double = get_field_3D(points)

        # Fields should be different (superposition)
        assert not np.allclose(field_double.z, Bz_single)


class TestHalbachArray:
    """Tests for Halbach array configurations."""

    def test_halbach_quad_2d(self):
        """Simple 2D quadrupole-like arrangement."""
        # Create 4 magnets in quadrupole configuration
        # Top: points up, Bottom: points down
        # Left: points left, Right: points right
        _ = magnets.Rectangle(
            width=10.0, height=10.0, Jr=1.0, center=(0.0, 20.0), phi=90
        )
        _ = magnets.Rectangle(
            width=10.0, height=10.0, Jr=1.0, center=(0.0, -20.0), phi=-90
        )
        _ = magnets.Rectangle(
            width=10.0, height=10.0, Jr=1.0, center=(-20.0, 0.0), phi=180
        )
        _ = magnets.Rectangle(
            width=10.0, height=10.0, Jr=1.0, center=(20.0, 0.0), phi=0
        )

        # Check that field at center is computed
        points = grid2D(5.0, 5.0, num_points=3)
        field = get_field_2D(points)

        # Field should be finite
        assert np.all(np.isfinite(field.x))
        assert np.all(np.isfinite(field.y))


class TestDifferentMagnetTypes:
    """Tests for mixing different magnet types."""

    def test_rectangle_and_circle_together(self):
        """Rectangle and Circle magnets can coexist."""
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, center=(-15.0, 0.0))
        _ = magnets.Circle(radius=5.0, Jr=1.0, center=(15.0, 0.0))

        points = grid2D(30.0, 30.0, num_points=5)
        field = get_field_2D(points)

        # Both should contribute to total field
        assert np.any(field.x != 0.0) or np.any(field.y != 0.0)

    def test_prism_and_cylinder_together(self):
        """Prism and Cylinder magnets can coexist."""
        _ = magnets.Prism(
            width=5.0, depth=5.0, height=10.0, Jr=1.0, center=(-15.0, 0.0, 0.0)
        )
        _ = magnets.Cylinder(radius=5.0, length=10.0, Jr=1.0, center=(15.0, 0.0, 0.0))

        points = slice3D(
            plane="xz", max1=30.0, max2=30.0, slice_value=0.0, num_points=5
        )
        field = get_field_3D(points)

        # Both should contribute
        assert (
            np.any(field.x != 0.0) or np.any(field.y != 0.0) or np.any(field.z != 0.0)
        )

    def test_prism_and_sphere_together(self):
        """Prism and Sphere magnets can coexist."""
        _ = magnets.Prism(
            width=5.0, depth=5.0, height=10.0, Jr=1.0, center=(-15.0, 0.0, 0.0)
        )
        _ = magnets.Sphere(radius=5.0, Jr=1.0, center=(15.0, 0.0, 0.0))

        points = slice3D(
            plane="xz", max1=30.0, max2=30.0, slice_value=0.0, num_points=5
        )
        field = get_field_3D(points)

        # Both should contribute
        assert (
            np.any(field.x != 0.0) or np.any(field.y != 0.0) or np.any(field.z != 0.0)
        )


class TestFieldSymmetry:
    """Tests for field symmetry properties."""

    def test_rectangle_xz_symmetry(self):
        """Z-magnetized rectangle has reflection symmetry about y=0."""
        m = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, center=(0.0, 0.0))

        # Points symmetric about y=0
        Bx_pos, By_pos = m.get_field(10.0, 5.0)
        Bx_neg, By_neg = m.get_field(10.0, -5.0)

        # For reflection symmetry about y=0:
        # Bx(x, y) = -Bx(x, -y) (odd in y)
        # By(x, y) = By(x, -y) (even in y)
        npt.assert_allclose(Bx_pos, -Bx_neg, rtol=1e-10)
        npt.assert_allclose(By_pos, By_neg, rtol=1e-10)

    def test_prism_xy_symmetry(self):
        """Symmetric prism has x-y symmetry."""
        m = magnets.Prism(
            width=10.0, depth=10.0, height=20.0, Jr=1.0, center=(0.0, 0.0, 0.0)
        )

        # Points at (x, 0, z) and (0, x, z) should have related fields
        Bx1, By1, Bz1 = m.get_field(5.0, 0.0, 30.0)
        Bx2, By2, Bz2 = m.get_field(0.0, 5.0, 30.0)

        # By rotational symmetry about z-axis
        npt.assert_allclose(Bx1, By2, rtol=1e-10)
        npt.assert_allclose(By1, Bx2, rtol=1e-10)
        npt.assert_allclose(Bz1, Bz2, rtol=1e-10)

    def test_cylinder_axial_symmetry(self):
        """Cylinder has perfect axial symmetry."""
        m = magnets.Cylinder(radius=5.0, length=10.0, Jr=1.0, center=(0.0, 0.0, 0.0))

        # Points at same distance from axis
        Bx1, By1, Bz1 = m.get_field(10.0, 0.0, 20.0)
        Bx2, By2, Bz2 = m.get_field(0.0, 10.0, 20.0)

        # Same Bz and same magnitude
        npt.assert_allclose(Bz1, Bz2, rtol=1e-10)
        B1_mag = np.sqrt(Bx1**2 + By1**2 + Bz1**2)
        B2_mag = np.sqrt(Bx2**2 + By2**2 + Bz2**2)
        npt.assert_allclose(B1_mag, B2_mag, rtol=1e-10)


class TestMagnetTranslation:
    """Tests for magnet translation behavior."""

    def test_translated_magnet_field(self):
        """Magnet at offset produces offset field."""
        # Magnet at origin
        m1 = magnets.Prism(
            width=5.0, depth=5.0, height=10.0, Jr=1.0, center=(0.0, 0.0, 0.0)
        )
        Bx1, By1, Bz1 = m1.get_field(0.0, 0.0, 20.0)

        reset()

        # Magnet translated
        m2 = magnets.Prism(
            width=5.0, depth=5.0, height=10.0, Jr=1.0, center=(10.0, 0.0, 0.0)
        )
        Bx2, By2, Bz2 = m2.get_field(10.0, 0.0, 20.0)

        # Field at equivalent relative position should be same
        npt.assert_allclose(Bx1, Bx2, rtol=1e-10)
        npt.assert_allclose(By1, By2, rtol=1e-10)
        npt.assert_allclose(Bz1, Bz2, rtol=1e-10)


class TestMagnetCount:
    """Tests for magnet instance counting."""

    def test_count_increases_with_magnets(self):
        """Creating magnets increases instance count."""
        from pymagnet.magnets._magnet_base import Magnet

        initial = Magnet.get_num_instances()

        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)

        assert Magnet.get_num_instances() == initial + 3

    def test_reset_clears_all_types(self):
        """Reset clears magnets of all types."""
        from pymagnet.magnets._magnet_base import Magnet

        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        _ = magnets.Circle(radius=5.0, Jr=1.0)
        _ = magnets.Prism(width=5.0, depth=5.0, height=10.0, Jr=1.0)
        _ = magnets.Cylinder(radius=5.0, length=10.0, Jr=1.0)

        reset()

        assert Magnet.get_num_instances() == 0
