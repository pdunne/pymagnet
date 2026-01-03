# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for Point2 and Point3 classes."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.utils import Point2, Point3


class TestPoint2Init:
    """Tests for Point2 initialization."""

    def test_init_scalar(self):
        """Initialize with scalar values."""
        p = Point2(1.0, 2.0)
        assert p.x == 1.0
        assert p.y == 2.0

    def test_init_array(self):
        """Initialize with numpy arrays."""
        x = np.array([1.0, 2.0])
        y = np.array([3.0, 4.0])
        p = Point2(x, y)
        npt.assert_array_equal(p.x, x)
        npt.assert_array_equal(p.y, y)

    def test_init_integers(self):
        """Initialize with integers."""
        p = Point2(1, 2)
        assert p.x == 1
        assert p.y == 2


class TestPoint2Arithmetic:
    """Tests for Point2 arithmetic operations."""

    def test_addition(self, point2_origin, point2_unit_x):
        """Point2 + Point2."""
        result = point2_origin + point2_unit_x
        assert result.x == 1.0
        assert result.y == 0.0

    def test_subtraction(self):
        """Point2 - Point2."""
        p1 = Point2(3.0, 4.0)
        p2 = Point2(1.0, 2.0)
        result = p1 - p2
        assert result.x == 2.0
        assert result.y == 2.0

    def test_multiplication(self):
        """Elementwise multiplication."""
        p1 = Point2(2.0, 3.0)
        p2 = Point2(4.0, 5.0)
        result = p1 * p2
        assert result.x == 8.0
        assert result.y == 15.0

    def test_division(self):
        """Elementwise division."""
        p1 = Point2(6.0, 8.0)
        p2 = Point2(2.0, 4.0)
        result = p1.__div__(p2)
        assert result.x == 3.0
        assert result.y == 2.0

    def test_division_by_zero(self):
        """Division by zero returns 0."""
        p1 = Point2(6.0, 8.0)
        p2 = Point2(0.0, 4.0)
        result = p1.__div__(p2)
        assert result.x == 0  # Division by zero returns 0
        assert result.y == 2.0


class TestPoint2Comparison:
    """Tests for Point2 comparison operations (magnitude-based)."""

    def test_lt(self):
        """Less than comparison (magnitude-based)."""
        p1 = Point2(1.0, 0.0)  # mag^2 = 1
        p2 = Point2(2.0, 0.0)  # mag^2 = 4
        assert p1 < p2
        assert not p2 < p1

    def test_le(self):
        """Less than or equal."""
        p1 = Point2(1.0, 0.0)
        p2 = Point2(1.0, 0.0)
        assert p1 <= p2
        assert p2 <= p1

    def test_gt(self):
        """Greater than comparison."""
        p1 = Point2(2.0, 0.0)
        p2 = Point2(1.0, 0.0)
        assert p1 > p2
        assert not p2 > p1

    def test_ge(self):
        """Greater than or equal."""
        p1 = Point2(2.0, 0.0)
        p2 = Point2(2.0, 0.0)
        assert p1 >= p2
        assert p2 >= p1

    def test_eq(self):
        """Equality comparison (exact match)."""
        p1 = Point2(1.0, 2.0)
        p2 = Point2(1.0, 2.0)
        assert p1 == p2

    def test_ne(self):
        """Not equal comparison."""
        p1 = Point2(1.0, 2.0)
        p2 = Point2(1.0, 3.0)
        assert p1 != p2


class TestPoint2Distance:
    """Tests for Point2 distance calculations."""

    def test_distance_to(self, point2_origin, point2_unit_x):
        """Distance between two points."""
        dist = point2_origin.distance_to(point2_unit_x)
        assert dist == 1.0

    def test_distance_to_diagonal(self):
        """Distance to diagonal point."""
        p1 = Point2(0.0, 0.0)
        p2 = Point2(3.0, 4.0)
        dist = p1.distance_to(p2)
        assert dist == 5.0  # 3-4-5 triangle

    def test_distance_to_origin(self):
        """Distance from origin."""
        p = Point2(3.0, 4.0)
        dist = p.distance_to_origin()
        assert dist == 5.0

    def test_distance_to_origin_from_origin(self, point2_origin):
        """Distance from origin to origin is 0."""
        dist = point2_origin.distance_to_origin()
        assert dist == 0.0


class TestPoint2Norm:
    """Tests for Point2 norm calculation."""

    def test_norm_unit_vector(self, point2_unit_x):
        """Norm of unit vector is 1."""
        norm = point2_unit_x._norm()
        assert norm == 1.0

    def test_norm_345(self):
        """Norm of 3-4 vector is 5."""
        p = Point2(3.0, 4.0)
        norm = p._norm()
        assert norm == 5.0

    def test_norm_zero(self, point2_origin):
        """Norm of zero vector is 0."""
        norm = point2_origin._norm()
        assert norm == 0.0


class TestPoint2Repr:
    """Tests for Point2 string representations."""

    def test_repr(self):
        """repr returns string format."""
        p = Point2(1.0, 2.0)
        assert repr(p) == "(1.0, 2.0)"

    def test_str(self):
        """str returns same format."""
        p = Point2(1.0, 2.0)
        assert str(p) == "(1.0, 2.0)"


# ==================== Point3 Tests ====================


class TestPoint3Init:
    """Tests for Point3 initialization."""

    def test_init_scalar(self):
        """Initialize with scalar values."""
        p = Point3(1.0, 2.0, 3.0)
        assert p.x == 1.0
        assert p.y == 2.0
        assert p.z == 3.0

    def test_init_array(self):
        """Initialize with numpy arrays."""
        x = np.array([1.0, 2.0])
        y = np.array([3.0, 4.0])
        z = np.array([5.0, 6.0])
        p = Point3(x, y, z)
        npt.assert_array_equal(p.x, x)
        npt.assert_array_equal(p.y, y)
        npt.assert_array_equal(p.z, z)

    def test_inherits_from_point2(self):
        """Point3 should inherit from Point2."""
        p = Point3(1.0, 2.0, 3.0)
        assert isinstance(p, Point2)


class TestPoint3Arithmetic:
    """Tests for Point3 arithmetic operations."""

    def test_addition(self):
        """Point3 + Point3."""
        p1 = Point3(1.0, 2.0, 3.0)
        p2 = Point3(4.0, 5.0, 6.0)
        result = p1 + p2
        assert result.x == 5.0
        assert result.y == 7.0
        assert result.z == 9.0

    def test_subtraction(self):
        """Point3 - Point3."""
        p1 = Point3(4.0, 5.0, 6.0)
        p2 = Point3(1.0, 2.0, 3.0)
        result = p1 - p2
        assert result.x == 3.0
        assert result.y == 3.0
        assert result.z == 3.0

    def test_multiplication(self):
        """Elementwise multiplication."""
        p1 = Point3(2.0, 3.0, 4.0)
        p2 = Point3(5.0, 6.0, 7.0)
        result = p1 * p2
        assert result.x == 10.0
        assert result.y == 18.0
        assert result.z == 28.0

    def test_division(self):
        """Elementwise division."""
        p1 = Point3(6.0, 8.0, 10.0)
        p2 = Point3(2.0, 4.0, 5.0)
        result = p1.__div__(p2)
        assert result.x == 3.0
        assert result.y == 2.0
        assert result.z == 2.0

    def test_division_by_zero_z(self):
        """Division by zero in z returns 0."""
        p1 = Point3(6.0, 8.0, 10.0)
        p2 = Point3(2.0, 4.0, 0.0)
        result = p1.__div__(p2)
        assert result.z == 0


class TestPoint3Distance:
    """Tests for Point3 distance calculations."""

    def test_distance_to(self, point3_origin, point3_unit_x):
        """Distance between two points."""
        dist = point3_origin.distance_to(point3_unit_x)
        assert dist == 1.0

    def test_distance_to_3d(self):
        """3D distance calculation."""
        p1 = Point3(0.0, 0.0, 0.0)
        p2 = Point3(1.0, 2.0, 2.0)
        dist = p1.distance_to(p2)
        assert dist == 3.0  # sqrt(1 + 4 + 4) = 3

    def test_distance_to_origin(self):
        """Distance from origin in 3D."""
        p = Point3(1.0, 2.0, 2.0)
        dist = p.distance_to_origin()
        assert dist == 3.0


class TestPoint3Norm:
    """Tests for Point3 norm calculation."""

    def test_norm_unit_vector(self, point3_unit_x):
        """Norm of unit vector is 1."""
        norm = point3_unit_x._norm()
        assert norm == 1.0

    def test_norm_3d(self):
        """3D norm calculation."""
        p = Point3(1.0, 2.0, 2.0)
        norm = p._norm()
        assert norm == 3.0

    def test_norm_zero(self, point3_origin):
        """Norm of zero vector is 0."""
        norm = point3_origin._norm()
        assert norm == 0.0


class TestPoint3Repr:
    """Tests for Point3 string representations."""

    def test_repr(self):
        """repr returns 3D string format."""
        p = Point3(1.0, 2.0, 3.0)
        assert repr(p) == "(1.0, 2.0, 3.0)"

    def test_str(self):
        """str returns same format."""
        p = Point3(1.0, 2.0, 3.0)
        assert str(p) == "(1.0, 2.0, 3.0)"


class TestPoint3Comparison:
    """Tests for Point3 comparison operations.

    Note: Point3 comparisons use only x,y (inherited from Point2).
    """

    def test_eq_same_xy_different_z(self):
        """Equality based on x,y only (z ignored in current impl)."""
        p1 = Point3(1.0, 2.0, 3.0)
        p2 = Point3(1.0, 2.0, 5.0)
        # Current implementation only compares x,y
        assert p1 == p2  # This documents current behavior

    def test_ne(self):
        """Not equal when x or y differ."""
        p1 = Point3(1.0, 2.0, 3.0)
        p2 = Point3(1.0, 3.0, 3.0)
        assert p1 != p2
