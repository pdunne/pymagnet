# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for 2D routines."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.utils._routines2D import (
    grid2D,
    get_field_2D,
    rotate_points_2D,
    gradB_2D,
    _allocate_field_array2,
)


class TestGrid2D:
    """Tests for grid2D function."""

    def test_default_returns_point_array2(self):
        """grid2D returns Point_Array2."""
        from pymagnet.utils._vector_structs import Point_Array2

        result = grid2D(10.0, 10.0)
        assert isinstance(result, Point_Array2)

    def test_default_num_points(self):
        """Default num_points is 100."""
        result = grid2D(10.0, 10.0)
        assert result.x.shape == (100, 100)
        assert result.y.shape == (100, 100)

    def test_custom_num_points(self):
        """Custom num_points is respected."""
        result = grid2D(10.0, 10.0, num_points=50)
        assert result.x.shape == (50, 50)

    def test_symmetric_limits(self):
        """Default limits are symmetric about origin."""
        result = grid2D(10.0, 5.0, num_points=3)
        # Should go from -10 to 10 in x, -5 to 5 in y
        npt.assert_allclose(result.x.min(), -10.0)
        npt.assert_allclose(result.x.max(), 10.0)
        npt.assert_allclose(result.y.min(), -5.0)
        npt.assert_allclose(result.y.max(), 5.0)

    def test_custom_min_values(self):
        """Custom xmin, ymin are respected."""
        result = grid2D(10.0, 10.0, xmin=0.0, ymin=0.0, num_points=3)
        npt.assert_allclose(result.x.min(), 0.0)
        npt.assert_allclose(result.y.min(), 0.0)

    def test_default_unit(self):
        """Default unit is mm."""
        result = grid2D(10.0, 10.0)
        assert result.unit == "mm"

    def test_custom_unit(self):
        """Custom unit is stored."""
        result = grid2D(10.0, 10.0, unit="cm")
        assert result.unit == "cm"


class TestGetField2D:
    """Tests for get_field_2D function."""

    def test_returns_field2(self):
        """get_field_2D returns Field2."""
        from pymagnet.utils._vector_structs import Field2

        points = grid2D(10.0, 10.0, num_points=5)
        result = get_field_2D(points)
        assert isinstance(result, Field2)

    def test_empty_registry_returns_zero_field(self):
        """With no magnets, field is zero everywhere."""
        points = grid2D(10.0, 10.0, num_points=5)
        result = get_field_2D(points)
        npt.assert_allclose(result.x, 0.0)
        npt.assert_allclose(result.y, 0.0)

    def test_single_magnet_field(self):
        """Single magnet produces non-zero field."""
        from pymagnet import magnets

        _ = magnets.Rectangle(width=5.0, height=10.0, Jr=1.0)
        points = grid2D(20.0, 20.0, num_points=5)
        result = get_field_2D(points)

        # Field should have non-zero values
        assert np.any(result.x != 0.0) or np.any(result.y != 0.0)

    def test_norm_computed(self):
        """Field norm is computed."""
        from pymagnet import magnets

        _ = magnets.Rectangle(width=5.0, height=10.0, Jr=1.0)
        points = grid2D(20.0, 20.0, num_points=5)
        result = get_field_2D(points)

        # Norm should be non-zero where field is non-zero
        assert result.n.shape == result.x.shape


class TestRotatePoints2D:
    """Tests for rotate_points_2D function."""

    def test_zero_rotation_identity(self):
        """Zero rotation returns original points."""
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.0, 0.0, 0.0])
        x_rot, y_rot = rotate_points_2D(x, y, 0.0)
        npt.assert_allclose(x_rot, x, atol=1e-10)
        npt.assert_allclose(y_rot, y, atol=1e-10)

    def test_90_degree_rotation(self):
        """90° rotation moves x to y."""
        x = np.array([1.0])
        y = np.array([0.0])
        x_rot, y_rot = rotate_points_2D(x, y, np.pi / 2)
        npt.assert_allclose(x_rot, [0.0], atol=1e-10)
        npt.assert_allclose(y_rot, [1.0], atol=1e-10)

    def test_180_degree_rotation(self):
        """180° rotation negates coordinates."""
        x = np.array([1.0])
        y = np.array([2.0])
        x_rot, y_rot = rotate_points_2D(x, y, np.pi)
        npt.assert_allclose(x_rot, [-1.0], atol=1e-10)
        npt.assert_allclose(y_rot, [-2.0], atol=1e-10)

    def test_preserves_distance_from_origin(self):
        """Rotation preserves distance from origin."""
        x = np.array([3.0])
        y = np.array([4.0])
        original_dist = np.sqrt(x**2 + y**2)

        x_rot, y_rot = rotate_points_2D(x, y, np.pi / 6)
        rotated_dist = np.sqrt(x_rot**2 + y_rot**2)

        npt.assert_allclose(rotated_dist, original_dist, rtol=1e-10)

    def test_preserves_shape(self):
        """Rotation preserves array shape."""
        x = np.array([[1.0, 2.0], [3.0, 4.0]])
        y = np.array([[0.0, 0.0], [0.0, 0.0]])
        x_rot, y_rot = rotate_points_2D(x, y, np.pi / 4)
        assert x_rot.shape == x.shape
        assert y_rot.shape == y.shape

    def test_mismatched_lengths_raises(self):
        """Mismatched x, y lengths raise exception."""
        x = np.array([1.0, 2.0])
        y = np.array([0.0])
        with pytest.raises(Exception):
            rotate_points_2D(x, y, 0.0)


class TestAllocateFieldArray2:
    """Tests for _allocate_field_array2 function."""

    def test_scalar_inputs(self):
        """Handles scalar inputs."""
        result = _allocate_field_array2(1.0, 2.0)
        assert result.x.shape == (1,)
        assert result.y.shape == (1,)

    def test_1d_array_inputs(self):
        """Handles 1D array inputs."""
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.0, 0.0, 0.0])
        result = _allocate_field_array2(x, y)
        assert result.x.shape == (3,)

    def test_2d_array_inputs(self):
        """Handles 2D grid inputs."""
        x, y = np.mgrid[0:5, 0:5]
        result = _allocate_field_array2(x, y)
        assert result.x.shape == (5, 5)

    def test_initializes_to_zero(self):
        """Allocated array is initialized to zero."""
        x = np.array([1.0, 2.0])
        y = np.array([3.0, 4.0])
        result = _allocate_field_array2(x, y)
        npt.assert_allclose(result.x, 0.0)
        npt.assert_allclose(result.y, 0.0)


class TestGradB2D:
    """Tests for gradB_2D function."""

    def test_uniform_field_zero_gradient(self):
        """Uniform field has zero gradient."""
        from pymagnet.utils._vector_structs import Field2

        x, y = np.mgrid[-5:5:10j, -5:5:10j]
        B = Field2(np.ones_like(x), np.zeros_like(x))
        B.n = np.ones_like(x)

        dB = gradB_2D(B.n, x, y)
        # Gradient of uniform field should be ~0
        npt.assert_allclose(dB.x, 0.0, atol=0.1)
        npt.assert_allclose(dB.y, 0.0, atol=0.1)

    def test_linear_field_constant_gradient(self):
        """Linear field has constant gradient."""
        from pymagnet.utils._vector_structs import Field2

        x, y = np.mgrid[0:10:11j, 0:10:11j]
        B = x  # Field increases linearly with x

        dB = gradB_2D(B, x, y)
        # Gradient in x should be ~1, in y should be ~0
        npt.assert_allclose(dB.x[1:-1, 1:-1], 1.0, atol=0.1)
        npt.assert_allclose(dB.y[1:-1, 1:-1], 0.0, atol=0.1)
