# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for force calculation functions."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.forces._mesh_force import (
    get_centroid,
    triangle_area,
    get_midpoints,
    divide_triangle_centroid,
    divide_triangle_regular,
)


class TestGetCentroid:
    """Tests for get_centroid function."""

    def test_equilateral_centroid(self):
        """Centroid of equilateral triangle."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, np.sqrt(3) / 2, 0.0]])
        centroid = get_centroid(triangle)
        expected = np.array([0.5, np.sqrt(3) / 6, 0.0])
        npt.assert_allclose(centroid, expected, rtol=1e-10)

    def test_right_angle_centroid(self):
        """Centroid of right-angle triangle."""
        triangle = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        centroid = get_centroid(triangle)
        expected = np.array([1.0, 4.0 / 3, 0.0])
        npt.assert_allclose(centroid, expected, rtol=1e-10)

    def test_3d_triangle_centroid(self):
        """Centroid of triangle in 3D space."""
        triangle = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        centroid = get_centroid(triangle)
        expected = np.array([4.0, 5.0, 6.0])
        npt.assert_allclose(centroid, expected, rtol=1e-10)

    def test_centroid_returns_3d_array(self):
        """Centroid returns (3,) array."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        centroid = get_centroid(triangle)
        assert centroid.shape == (3,)


class TestTriangleArea:
    """Tests for triangle_area function."""

    def test_unit_right_triangle_area(self):
        """Area of unit right triangle is 0.5."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        area = triangle_area(triangle)
        npt.assert_allclose(area, 0.5, rtol=1e-10)

    def test_equilateral_triangle_area(self):
        """Area of equilateral triangle with side 1."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, np.sqrt(3) / 2, 0.0]])
        area = triangle_area(triangle)
        expected = np.sqrt(3) / 4
        npt.assert_allclose(area, expected, rtol=1e-10)

    def test_345_triangle_area(self):
        """Area of 3-4-5 right triangle."""
        triangle = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        area = triangle_area(triangle)
        npt.assert_allclose(area, 6.0, rtol=1e-10)

    def test_degenerate_triangle_area(self):
        """Degenerate (collinear) triangle has zero area."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        area = triangle_area(triangle)
        npt.assert_allclose(area, 0.0, atol=1e-10)

    def test_3d_triangle_area(self):
        """Area of triangle in 3D space."""
        # Triangle in xz plane
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        area = triangle_area(triangle)
        npt.assert_allclose(area, 0.5, rtol=1e-10)


class TestGetMidpoints:
    """Tests for get_midpoints function."""

    def test_midpoints_count(self):
        """get_midpoints returns 3 midpoints."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        midpoints = get_midpoints(triangle)
        assert midpoints.shape == (3, 3)

    def test_midpoints_locations(self):
        """Midpoints are at correct locations."""
        triangle = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        midpoints = get_midpoints(triangle)

        # Midpoint of edge 0-1
        npt.assert_allclose(midpoints[0], [1.0, 0.0, 0.0], rtol=1e-10)
        # Midpoint of edge 1-2
        npt.assert_allclose(midpoints[1], [1.0, 1.0, 0.0], rtol=1e-10)
        # Midpoint of edge 2-0
        npt.assert_allclose(midpoints[2], [0.0, 1.0, 0.0], rtol=1e-10)

    def test_midpoints_on_edges(self):
        """All midpoints lie on triangle edges."""
        triangle = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.5, 3.0, 0.0]])
        midpoints = get_midpoints(triangle)

        # Check that each midpoint is equidistant from two vertices
        for i in range(3):
            v1 = triangle[i]
            v2 = triangle[(i + 1) % 3]
            dist1 = np.linalg.norm(midpoints[i] - v1)
            dist2 = np.linalg.norm(midpoints[i] - v2)
            npt.assert_allclose(dist1, dist2, rtol=1e-10)


class TestDivideTriangleCentroid:
    """Tests for divide_triangle_centroid function."""

    def test_depth_1_creates_3_triangles(self):
        """Depth 1 creates 3 sub-triangles."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = divide_triangle_centroid(triangle, depth=1)
        assert result.shape[0] == 3
        assert result.shape[1] == 3
        assert result.shape[2] == 3

    def test_depth_2_creates_9_triangles(self):
        """Depth 2 creates 9 (3^2) sub-triangles."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = divide_triangle_centroid(triangle, depth=2)
        assert result.shape[0] == 9

    def test_depth_3_creates_27_triangles(self):
        """Depth 3 creates 27 (3^3) sub-triangles."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = divide_triangle_centroid(triangle, depth=3)
        assert result.shape[0] == 27

    def test_total_area_preserved(self):
        """Total area of sub-triangles equals original area."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        original_area = triangle_area(triangle)

        result = divide_triangle_centroid(triangle, depth=2)
        total_area = sum(triangle_area(t) for t in result)

        npt.assert_allclose(total_area, original_area, rtol=1e-10)


class TestDivideTriangleRegular:
    """Tests for divide_triangle_regular function."""

    def test_depth_1_creates_4_triangles(self):
        """Depth 1 creates 4 sub-triangles."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = divide_triangle_regular(triangle, depth=1)
        assert result.shape[0] == 4
        assert result.shape[1] == 3
        assert result.shape[2] == 3

    def test_depth_2_creates_16_triangles(self):
        """Depth 2 creates 16 (4^2) sub-triangles."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = divide_triangle_regular(triangle, depth=2)
        assert result.shape[0] == 16

    def test_depth_3_creates_64_triangles(self):
        """Depth 3 creates 64 (4^3) sub-triangles."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = divide_triangle_regular(triangle, depth=3)
        assert result.shape[0] == 64

    def test_total_area_preserved(self):
        """Total area of sub-triangles equals original area."""
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        original_area = triangle_area(triangle)

        result = divide_triangle_regular(triangle, depth=2)
        total_area = sum(triangle_area(t) for t in result)

        npt.assert_allclose(total_area, original_area, rtol=1e-10)

    def test_all_subtriangles_congruent_depth_1(self):
        """At depth 1, all 4 sub-triangles have equal area."""
        triangle = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        result = divide_triangle_regular(triangle, depth=1)

        areas = [triangle_area(t) for t in result]
        for area in areas:
            npt.assert_allclose(area, areas[0], rtol=1e-10)


class TestGetAreaTriangles:
    """Tests for get_area_triangles function."""

    def test_computes_multiple_areas(self):
        """Computes areas for multiple triangles."""
        from pymagnet.forces._mesh_force import get_area_triangles

        triangles = np.array(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            ]
        )
        areas = np.zeros(2)
        get_area_triangles(triangles, areas)

        npt.assert_allclose(areas[0], 0.5, rtol=1e-10)
        npt.assert_allclose(areas[1], 2.0, rtol=1e-10)


class TestNumericalStability:
    """Tests for numerical stability of force calculations."""

    def test_small_triangle_centroid(self):
        """Centroid works for small triangles."""
        triangle = np.array(
            [[0.0, 0.0, 0.0], [1e-6, 0.0, 0.0], [0.5e-6, 1e-6, 0.0]]
        )
        centroid = get_centroid(triangle)
        assert np.all(np.isfinite(centroid))

    def test_large_triangle_centroid(self):
        """Centroid works for large triangles."""
        triangle = np.array([[0.0, 0.0, 0.0], [1e6, 0.0, 0.0], [0.5e6, 1e6, 0.0]])
        centroid = get_centroid(triangle)
        assert np.all(np.isfinite(centroid))

    def test_small_triangle_area(self):
        """Area works for small triangles."""
        triangle = np.array(
            [[0.0, 0.0, 0.0], [1e-6, 0.0, 0.0], [0.5e-6, 1e-6, 0.0]]
        )
        area = triangle_area(triangle)
        assert np.isfinite(area)
        assert area >= 0

    def test_division_preserves_centroid(self):
        """Centroid of divided triangles average to original centroid."""
        triangle = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.5, 3.0, 0.0]])
        original_centroid = get_centroid(triangle)

        divided = divide_triangle_regular(triangle, depth=1)
        centroids = np.array([get_centroid(t) for t in divided])
        areas = np.array([triangle_area(t) for t in divided])

        # Weighted average by area
        weighted_centroid = np.sum(
            centroids * areas[:, np.newaxis], axis=0
        ) / np.sum(areas)

        npt.assert_allclose(weighted_centroid, original_centroid, rtol=1e-10)
