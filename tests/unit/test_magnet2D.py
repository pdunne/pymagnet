# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for 2D magnet classes: Magnet2D, Rectangle, Square, Circle."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet import magnets
from pymagnet.utils import Point2


class TestRectangleInit:
    """Tests for Rectangle initialization."""

    def test_dimensions(self, rectangle_default):
        """Width and height properties set correctly."""
        assert rectangle_default.width == 20.0
        assert rectangle_default.height == 40.0

    def test_half_dimensions(self, rectangle_default):
        """Half dimensions a=width/2, b=height/2."""
        assert rectangle_default.a == 10.0
        assert rectangle_default.b == 20.0

    def test_default_center(self, rectangle_default):
        """Default center is origin."""
        center = rectangle_default.get_center()
        npt.assert_allclose(center, [0.0, 0.0])

    def test_custom_center(self, rectangle_offset):
        """Custom center from kwarg."""
        center = rectangle_offset.get_center()
        npt.assert_allclose(center, [10.0, 5.0])

    def test_jr_attribute(self, rectangle_default):
        """Jr property set correctly."""
        assert rectangle_default.Jr == 1.0

    def test_default_phi(self):
        """Default phi=90 (magnetized in y direction)."""
        r = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        assert r.phi == 90.0

    def test_custom_phi(self):
        """Custom phi angle."""
        r = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, phi=45.0)
        assert r.phi == 45.0

    def test_jx_jy_components_default(self, rectangle_default):
        """Jx, Jy computed from Jr and phi (default phi=90)."""
        # phi=90 means magnetization in y direction
        npt.assert_allclose(rectangle_default.Jx, 0.0, atol=1e-10)
        npt.assert_allclose(rectangle_default.Jy, 1.0)

    def test_jx_jy_components_phi_0(self):
        """Jx, Jy with phi=0 (magnetization in x direction)."""
        r = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, phi=0.0)
        npt.assert_allclose(r.Jx, 1.0)
        npt.assert_allclose(r.Jy, 0.0, atol=1e-10)

    def test_get_size(self, rectangle_default):
        """get_size returns [width, height]."""
        size = rectangle_default.get_size()
        assert size[0] == 20.0
        assert size[1] == 40.0

    def test_get_jr(self, rectangle_default):
        """get_Jr returns [Jx, Jy]."""
        jr = rectangle_default.get_Jr()
        npt.assert_allclose(jr[0], 0.0, atol=1e-10)  # Jx
        npt.assert_allclose(jr[1], 1.0)  # Jy

    def test_alpha_rotation(self, rectangle_rotated):
        """Alpha rotation angle set correctly."""
        orientation = rectangle_rotated.get_orientation()
        assert orientation == 45.0

    def test_str_repr(self, rectangle_default):
        """String representations contain class name."""
        str_rep = str(rectangle_default)
        assert "Rectangle" in str_rep


class TestRectangleField:
    """Tests for Rectangle field calculations."""

    def test_field_returns_tuple(self, rectangle_default):
        """get_field returns Bx, By tuple."""
        Bx, By = rectangle_default.get_field(0.0, 50.0)
        assert isinstance(Bx, (float, np.floating, np.ndarray))
        assert isinstance(By, (float, np.floating, np.ndarray))

    def test_field_outside_magnet(self, rectangle_default):
        """Field outside magnet is finite and non-zero."""
        Bx, By = rectangle_default.get_field(0.0, 50.0)
        assert np.isfinite(Bx)
        assert np.isfinite(By)
        # Field should be non-zero outside
        assert abs(By) > 0

    def test_field_decay_with_distance(self, rectangle_default):
        """Field magnitude decreases with distance."""
        _, By_close = rectangle_default.get_field(0.0, 25.0)
        _, By_far = rectangle_default.get_field(0.0, 100.0)
        assert abs(By_close) > abs(By_far)

    def test_field_symmetry_y_axis(self, rectangle_default):
        """Field is symmetric about y-axis for centered magnet."""
        Bx_pos, By_pos = rectangle_default.get_field(10.0, 50.0)
        Bx_neg, By_neg = rectangle_default.get_field(-10.0, 50.0)
        # By should be same, Bx should be opposite
        npt.assert_allclose(By_pos, By_neg, rtol=1e-10)
        npt.assert_allclose(Bx_pos, -Bx_neg, rtol=1e-10)

    def test_field_array_input(self, rectangle_default):
        """Field calculation with array inputs."""
        x = np.array([0.0, 10.0, -10.0])
        y = np.array([50.0, 50.0, 50.0])
        Bx, By = rectangle_default.get_field(x, y)
        assert Bx.shape == (3,)
        assert By.shape == (3,)

    @pytest.mark.parametrize("phi", [0.0, 45.0, 90.0, 135.0, 180.0])
    def test_magnetization_directions(self, phi):
        """Field calculation works for various magnetization directions."""
        r = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, phi=phi)
        Bx, By = r.get_field(0.0, 30.0)
        assert np.isfinite(Bx)
        assert np.isfinite(By)


class TestSquare:
    """Tests for Square class (inherits from Rectangle)."""

    def test_equal_dimensions(self, square_default):
        """Width equals height."""
        assert square_default.width == square_default.height

    def test_inherits_rectangle(self, square_default):
        """Square is subclass of Rectangle."""
        assert isinstance(square_default, magnets.Rectangle)

    def test_get_size_square(self, square_default):
        """get_size returns equal dimensions."""
        size = square_default.get_size()
        assert size[0] == size[1] == 20.0


class TestCircleInit:
    """Tests for Circle initialization."""

    def test_radius(self, circle_default):
        """Radius property set correctly."""
        assert circle_default.radius == 10.0

    def test_default_phi(self):
        """Default phi=0 (magnetized in x direction)."""
        c = magnets.Circle(radius=10.0, Jr=1.0)
        assert c.phi == 0.0

    def test_get_size(self, circle_default):
        """get_size returns [radius]."""
        size = circle_default.get_size()
        assert size[0] == 10.0

    def test_jr_attribute(self, circle_default):
        """Jr property set correctly."""
        assert circle_default.Jr == 1.0


class TestCircleField:
    """Tests for Circle field calculations."""

    def test_field_returns_tuple(self, circle_default):
        """get_field returns Bx, By tuple."""
        Bx, By = circle_default.get_field(20.0, 0.0)
        assert isinstance(Bx, (float, np.floating, np.ndarray))
        assert isinstance(By, (float, np.floating, np.ndarray))

    def test_field_outside_magnet(self, circle_default):
        """Field outside magnet is finite."""
        Bx, By = circle_default.get_field(20.0, 0.0)
        assert np.isfinite(Bx)
        assert np.isfinite(By)

    def test_field_decay_with_distance(self, circle_default):
        """Field decays with distance."""
        Bx_close, _ = circle_default.get_field(15.0, 0.0)
        Bx_far, _ = circle_default.get_field(50.0, 0.0)
        assert abs(Bx_close) > abs(Bx_far)


class TestMagnet2DCommon:
    """Tests common to all 2D magnets."""

    @pytest.mark.parametrize("Jr", [0.5, 1.0, -1.0, 2.0])
    def test_magnetization_scaling(self, Jr):
        """Field scales with magnetization."""
        r1 = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        r2 = magnets.Rectangle(width=10.0, height=20.0, Jr=Jr)

        _, By1 = r1.get_field(0.0, 30.0)
        _, By2 = r2.get_field(0.0, 30.0)

        npt.assert_allclose(By2, By1 * Jr, rtol=1e-10)

    @pytest.mark.parametrize(
        "center", [(0.0, 0.0), (10.0, 0.0), (0.0, 10.0), (-5.0, -5.0)]
    )
    def test_center_offset(self, center):
        """Magnet can be created at various centers."""
        r = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0, center=center)
        result_center = r.get_center()
        npt.assert_allclose(result_center, center)
