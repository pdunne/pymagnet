# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for triangle rotation and manipulation functions.

CRITICAL: These are Numba JIT compiled functions.
Tests must verify:
- Correct compilation of @jit decorated functions
- Numerical accuracy
- Edge cases (degenerate triangles)
"""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.utils._trigonometry3D import (
    _largest_side_RA,
    _rotate_triangle,
    align_triangle_to_y,
    align_triangle_xz,
    altitude,
    check_sign,
    norm_plane,
    return_axis_vector,
    return_z_vector,
    rotate_points,
    signed_area,
)
from pymagnet.utils._quaternion import Quaternion, q_angle_from_axis
from pymagnet.utils.global_const import PI


class TestSignedArea:
    """Tests for signed_area function (Numba JIT)."""

    @pytest.mark.numba
    def test_numba_compiles(self, equilateral_triangle):
        """signed_area should compile on first call."""
        result = signed_area(equilateral_triangle)
        assert np.isfinite(result)

    def test_equilateral_area(self, equilateral_triangle):
        """Equilateral triangle with unit side has area sqrt(3)/4."""
        area = signed_area(equilateral_triangle)
        expected_area = np.sqrt(3) / 4
        # Note: signed_area assumes triangle is in xz plane
        npt.assert_allclose(abs(area), expected_area, rtol=1e-10)

    def test_right_angle_triangle_xz(self, right_angle_triangle_xz):
        """Right triangle in xz plane has area 0.5."""
        area = signed_area(right_angle_triangle_xz)
        npt.assert_allclose(abs(area), 0.5, rtol=1e-10)

    def test_clockwise_negative_area(self):
        """Clockwise ordering gives negative area."""
        # Counter-clockwise triangle in xz
        ccw = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        # Clockwise (reversed)
        cw = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        area_ccw = signed_area(ccw)
        area_cw = signed_area(cw)

        # Areas should have opposite signs
        assert area_ccw * area_cw < 0

    def test_degenerate_zero_area(self, degenerate_triangle):
        """Collinear points should give zero area."""
        area = signed_area(degenerate_triangle)
        npt.assert_allclose(area, 0.0, atol=1e-10)


class TestNormPlane:
    """Tests for norm_plane function."""

    def test_xy_plane_normal(self, right_angle_triangle):
        """Triangle in xy plane has normal along z."""
        normal = norm_plane(right_angle_triangle)
        # Normal should be along z (positive or negative)
        npt.assert_allclose(abs(normal[2]), 1.0, rtol=1e-10)
        npt.assert_allclose(normal[0], 0.0, atol=1e-10)
        npt.assert_allclose(normal[1], 0.0, atol=1e-10)

    def test_xz_plane_normal(self, right_angle_triangle_xz):
        """Triangle in xz plane has normal along y."""
        normal = norm_plane(right_angle_triangle_xz)
        # Normal should be along y (positive or negative)
        npt.assert_allclose(abs(normal[1]), 1.0, rtol=1e-10)
        npt.assert_allclose(normal[0], 0.0, atol=1e-10)
        npt.assert_allclose(normal[2], 0.0, atol=1e-10)

    def test_unit_normal(self, equilateral_triangle):
        """Returned normal should be unit vector."""
        normal = norm_plane(equilateral_triangle)
        norm_length = np.linalg.norm(normal)
        npt.assert_allclose(norm_length, 1.0, rtol=1e-10)

    def test_arbitrary_triangle(self, triangle_arbitrary):
        """Normal of arbitrary triangle is unit vector."""
        normal = norm_plane(triangle_arbitrary)
        norm_length = np.linalg.norm(normal)
        npt.assert_allclose(norm_length, 1.0, rtol=1e-10)


class TestAltitude:
    """Tests for altitude calculation using Heron's formula."""

    def test_equilateral_triangle(self):
        """Altitude of equilateral triangle = sqrt(3)/2 * side."""
        # Unit equilateral: all sides = 1
        alt = altitude(1.0, 1.0, 1.0)
        expected = np.sqrt(3) / 2
        npt.assert_allclose(alt, expected, rtol=1e-10)

    def test_right_triangle_345(self):
        """Altitude of 3-4-5 right triangle to hypotenuse."""
        # 3-4-5 right triangle, altitude to hypotenuse (5)
        alt = altitude(5.0, 3.0, 4.0)
        # Area = 0.5 * 3 * 4 = 6, altitude = 2 * Area / base = 2.4
        expected = 2.4
        npt.assert_allclose(alt, expected, rtol=1e-10)

    def test_isoceles_triangle(self):
        """Altitude of isoceles triangle."""
        # Isoceles with base=2, sides=sqrt(2)
        alt = altitude(2.0, np.sqrt(2), np.sqrt(2))
        # Height should be 1
        npt.assert_allclose(alt, 1.0, rtol=1e-10)

    def test_very_thin_triangle(self):
        """Very thin triangle with small altitude."""
        # Long thin triangle (base=10, sides slightly longer than 5)
        # For a valid triangle: a < b + c, so 10 < 5.1 + 5.1 = 10.2 (barely valid)
        alt = altitude(10.0, 5.1, 5.1)
        assert alt > 0
        assert alt < 2  # Should be relatively small


class TestLargestSideRA:
    """Tests for _largest_side_RA function."""

    def test_identifies_longest_side(self):
        """Identifies longest side correctly."""
        # Triangle with side 0 being longest
        triangle = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [5.0, 0.0, 3.0]])
        longest_side, RA1, RA2 = _largest_side_RA(triangle)
        assert longest_side == 0  # Side from vertex 0 to vertex 1

    def test_returns_right_angle_triangles(self, equilateral_triangle):
        """Returns two right-angled triangle dimensions."""
        longest_side, RA1, RA2 = _largest_side_RA(equilateral_triangle)
        # RA1 and RA2 should have [base, height] format
        assert len(RA1) == 2
        assert len(RA2) == 2
        assert RA1[0] > 0  # base
        assert RA1[1] > 0  # height

    def test_total_base_equals_longest_side(self, equilateral_triangle):
        """Sum of RA triangle bases should equal longest side."""
        longest_side, RA1, RA2 = _largest_side_RA(equilateral_triangle)
        total_base = RA1[0] + RA2[0]
        # For equilateral, longest side = 1
        npt.assert_allclose(total_base, 1.0, rtol=1e-10)

    def test_altitude_preserved(self):
        """Both RA triangles have same altitude."""
        triangle = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.0, 0.0, 2.0]])
        longest_side, RA1, RA2 = _largest_side_RA(triangle)
        npt.assert_allclose(RA1[1], RA2[1], rtol=1e-10)


class TestCheckSign:
    """Tests for check_sign function."""

    def test_same_positive_signs(self):
        """Same positive signs returns True."""
        v1 = np.array([1.0, 2.0, 3.0])
        v2 = np.array([4.0, 5.0, 6.0])
        assert check_sign(v1, v2) is True

    def test_same_negative_signs(self):
        """Same negative signs returns True."""
        v1 = np.array([-1.0, -2.0, -3.0])
        v2 = np.array([-4.0, -5.0, -6.0])
        assert check_sign(v1, v2) is True

    def test_opposite_signs(self):
        """Opposite signs returns False."""
        v1 = np.array([1.0, 2.0, 3.0])
        v2 = np.array([-1.0, -2.0, -3.0])
        assert check_sign(v1, v2) is False

    def test_mixed_signs(self):
        """Mixed signs returns False."""
        v1 = np.array([1.0, -2.0, 3.0])
        v2 = np.array([1.0, 2.0, 3.0])
        assert check_sign(v1, v2) is False

    def test_zero_elements(self):
        """Zero elements are treated as same sign."""
        v1 = np.array([0.0, 1.0, 2.0])
        v2 = np.array([0.0, 3.0, 4.0])
        assert check_sign(v1, v2) is True


class TestReturnAxisVector:
    """Tests for return_axis_vector function."""

    def test_returns_unit_vector(self, equilateral_triangle):
        """Returns normalized vector."""
        vec = return_axis_vector(equilateral_triangle, 0)
        norm = np.linalg.norm(vec)
        npt.assert_allclose(norm, 1.0, rtol=1e-10)

    def test_correct_side_selection(self):
        """Returns vector for correct side."""
        triangle = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]])

        # Side 0: from vertex 0 to vertex 1 -> along x
        vec0 = return_axis_vector(triangle, 0)
        npt.assert_allclose(abs(vec0[0]), 1.0, rtol=1e-10)

        # Side 2: from vertex 0 to vertex 2 -> along y
        vec2 = return_axis_vector(triangle, 2)
        npt.assert_allclose(abs(vec2[1]), 1.0, rtol=1e-10)


class TestReturnZVector:
    """Tests for return_z_vector function."""

    def test_returns_unit_vector(self, equilateral_triangle):
        """Returns normalized vector."""
        vec = return_z_vector(equilateral_triangle, 0)
        norm = np.linalg.norm(vec)
        npt.assert_allclose(norm, 1.0, rtol=1e-10)


class TestAlignTriangleToY:
    """Tests for align_triangle_to_y function."""

    def test_xy_plane_triangle_rotates_to_xz(self, right_angle_triangle):
        """Triangle in xy plane gets rotated to xz plane (normal -> y)."""
        norm_vec = norm_plane(right_angle_triangle)
        y_axis = np.array([0, 1, 0])
        rot_axis = np.cross(y_axis, norm_vec)

        aligned, rotation = align_triangle_to_y(right_angle_triangle, rot_axis, norm_vec)

        # Check that aligned triangle has normal along y
        new_normal = norm_plane(aligned)
        npt.assert_allclose(abs(new_normal[1]), 1.0, rtol=1e-10)

    def test_already_aligned_unchanged(self, right_angle_triangle_xz):
        """Triangle already in xz plane stays unchanged or gets 180 rotation."""
        norm_vec = norm_plane(right_angle_triangle_xz)
        y_axis = np.array([0, 1, 0])
        rot_axis = np.cross(y_axis, norm_vec)

        aligned, rotation = align_triangle_to_y(
            right_angle_triangle_xz, rot_axis, norm_vec
        )

        # Rotation should be identity OR 180 degree rotation (parallel or anti-parallel)
        # Check that it's either identity (w=1) or 180 rotation (w=0)
        assert abs(rotation.w) > 0.99 or abs(rotation.w) < 0.01


class TestAlignTriangleXZ:
    """Tests for align_triangle_xz function."""

    def test_returns_quaternions(self, right_angle_triangle_xz):
        """Returns two quaternion rotations."""
        second_rot, third_rot = align_triangle_xz(right_angle_triangle_xz, 0)
        assert isinstance(second_rot, Quaternion)
        assert isinstance(third_rot, Quaternion)


class TestRotateTriangle:
    """Tests for _rotate_triangle function - main rotation pipeline."""

    def test_returns_correct_tuple(self, equilateral_triangle):
        """Returns (rotation, triangle, offset, RA1, RA2)."""
        Jr = np.array([0.0, 1.0, 0.0])
        result = _rotate_triangle(equilateral_triangle, Jr)
        assert len(result) == 5

        total_rotation, rotated_tri, offset, RA1, RA2 = result
        assert isinstance(total_rotation, Quaternion)
        assert rotated_tri.shape == (3, 3)
        assert offset.shape == (3,)
        assert len(RA1) == 2
        assert len(RA2) == 2

    def test_shape_preserved(self, equilateral_triangle):
        """Triangle shape (side lengths) preserved after rotation."""
        Jr = np.array([0.0, 1.0, 0.0])
        _, rotated_tri, _, _, _ = _rotate_triangle(equilateral_triangle, Jr)

        # Calculate side lengths of original
        orig_sides = [
            np.linalg.norm(equilateral_triangle[1] - equilateral_triangle[0]),
            np.linalg.norm(equilateral_triangle[2] - equilateral_triangle[1]),
            np.linalg.norm(equilateral_triangle[2] - equilateral_triangle[0]),
        ]

        # Calculate side lengths of rotated
        rot_sides = [
            np.linalg.norm(rotated_tri[1] - rotated_tri[0]),
            np.linalg.norm(rotated_tri[2] - rotated_tri[1]),
            np.linalg.norm(rotated_tri[2] - rotated_tri[0]),
        ]

        # Sides should match (order might differ due to rotation)
        orig_sorted = np.sort(orig_sides)
        rot_sorted = np.sort(rot_sides)
        npt.assert_allclose(orig_sorted, rot_sorted, rtol=1e-10)

    def test_area_preserved(self, triangle_arbitrary):
        """Area preserved after rotation."""
        Jr = np.array([0.0, 1.0, 0.0])
        _, rotated_tri, _, _, _ = _rotate_triangle(triangle_arbitrary, Jr)

        # Area via cross product
        def triangle_area(tri):
            v1 = tri[1] - tri[0]
            v2 = tri[2] - tri[0]
            return 0.5 * np.linalg.norm(np.cross(v1, v2))

        orig_area = triangle_area(triangle_arbitrary)
        rot_area = triangle_area(rotated_tri)
        npt.assert_allclose(orig_area, rot_area, rtol=1e-10)

    def test_ra_triangles_sum_to_original(self, equilateral_triangle):
        """RA triangle bases sum to original longest side."""
        Jr = np.array([0.0, 1.0, 0.0])
        _, _, _, RA1, RA2 = _rotate_triangle(equilateral_triangle, Jr)

        # For equilateral, all sides = 1, so bases should sum to 1
        total_base = RA1[0] + RA2[0]
        npt.assert_allclose(total_base, 1.0, rtol=1e-10)


class TestRotatePoints:
    """Tests for rotate_points function."""

    def test_identity_rotation(self, equilateral_triangle, identity_quaternion):
        """Identity rotation returns same points."""
        rotated = rotate_points(equilateral_triangle, identity_quaternion)
        npt.assert_allclose(rotated, equilateral_triangle, rtol=1e-10)

    def test_90_degree_rotation(self, right_angle_triangle_xz, rotation_90_y):
        """90 degree rotation about y-axis."""
        rotated = rotate_points(right_angle_triangle_xz, rotation_90_y)
        # Point at (1,0,0) should go to (0,0,-1) after 90 deg about y
        # Check that rotated points are different from original
        assert not np.allclose(rotated, right_angle_triangle_xz)


class TestNumericalStability:
    """Tests for numerical edge cases in trigonometry functions."""

    def test_nearly_degenerate_triangle(self):
        """Nearly degenerate triangle (very thin)."""
        # Very thin triangle
        triangle = np.array(
            [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [5.0, 0.0, 0.001]]
        )
        Jr = np.array([0.0, 1.0, 0.0])

        # Should not raise an error
        result = _rotate_triangle(triangle, Jr)
        assert result is not None

    def test_triangle_with_large_coordinates(self):
        """Triangle with large coordinate values."""
        triangle = np.array(
            [[1000.0, 0.0, 1000.0], [1001.0, 0.0, 1000.0], [1000.5, 0.0, 1000.866]]
        )
        Jr = np.array([0.0, 1.0, 0.0])

        result = _rotate_triangle(triangle, Jr)
        assert result is not None

    def test_triangle_with_small_coordinates(self):
        """Triangle with small coordinate values."""
        triangle = np.array(
            [[0.0, 0.0, 0.0], [1e-6, 0.0, 0.0], [0.5e-6, 0.0, 0.866e-6]]
        )
        Jr = np.array([0.0, 1.0, 0.0])

        result = _rotate_triangle(triangle, Jr)
        assert result is not None
