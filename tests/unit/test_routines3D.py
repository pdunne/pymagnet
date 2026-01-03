# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for 3D routines."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.utils._routines3D import (
    grid3D,
    line3D,
    slice3D,
    point3D,
    plane3D,
    get_field_3D,
    _allocate_field_array3,
    _get_max_array,
    _tile_arrays,
)


class TestGrid3D:
    """Tests for grid3D function."""

    def test_returns_point_array3(self):
        """grid3D returns Point_Array3."""
        from pymagnet.utils._vector_structs import Point_Array3

        result = grid3D(10.0, 10.0, 10.0, num_points=5)
        assert isinstance(result, Point_Array3)

    def test_default_num_points(self):
        """Default num_points is 100."""
        result = grid3D(10.0, 10.0, 10.0)
        assert result.x.shape == (100, 100, 100)

    def test_custom_num_points(self):
        """Custom num_points is respected."""
        result = grid3D(10.0, 10.0, 10.0, num_points=10)
        assert result.x.shape == (10, 10, 10)

    def test_asymmetric_num_points(self):
        """Different num_points per axis."""
        result = grid3D(10.0, 10.0, 10.0, num_points_x=5, num_points_y=10, num_points_z=15)
        assert result.x.shape == (5, 10, 15)

    def test_symmetric_limits(self):
        """Default limits are symmetric about origin."""
        result = grid3D(10.0, 5.0, 2.0, num_points=3)
        npt.assert_allclose(result.x.min(), -10.0)
        npt.assert_allclose(result.x.max(), 10.0)
        npt.assert_allclose(result.y.min(), -5.0)
        npt.assert_allclose(result.y.max(), 5.0)
        npt.assert_allclose(result.z.min(), -2.0)
        npt.assert_allclose(result.z.max(), 2.0)

    def test_custom_min_values(self):
        """Custom xmin, ymin, zmin are respected."""
        result = grid3D(10.0, 10.0, 10.0, xmin=0.0, ymin=0.0, zmin=0.0, num_points=3)
        npt.assert_allclose(result.x.min(), 0.0)
        npt.assert_allclose(result.y.min(), 0.0)
        npt.assert_allclose(result.z.min(), 0.0)

    def test_default_unit(self):
        """Default unit is mm."""
        result = grid3D(10.0, 10.0, 10.0, num_points=3)
        assert result.unit == "mm"

    def test_custom_unit(self):
        """Custom unit is stored."""
        result = grid3D(10.0, 10.0, 10.0, unit="cm", num_points=3)
        assert result.unit == "cm"


class TestLine3D:
    """Tests for line3D function."""

    def test_returns_point_array3(self):
        """line3D returns Point_Array3."""
        from pymagnet.utils._vector_structs import Point_Array3

        result = line3D((0, 0, 0), (1, 1, 1))
        assert isinstance(result, Point_Array3)

    def test_default_num_points(self):
        """Default num_points is 100."""
        result = line3D((0, 0, 0), (1, 1, 1))
        assert result.x.shape == (100,)

    def test_custom_num_points(self):
        """Custom num_points is respected."""
        result = line3D((0, 0, 0), (1, 1, 1), num_points=50)
        assert result.x.shape == (50,)

    def test_start_end_points(self):
        """Line starts and ends at correct points."""
        result = line3D((1, 2, 3), (4, 5, 6), num_points=10)
        npt.assert_allclose(result.x[0], 1.0)
        npt.assert_allclose(result.y[0], 2.0)
        npt.assert_allclose(result.z[0], 3.0)
        npt.assert_allclose(result.x[-1], 4.0)
        npt.assert_allclose(result.y[-1], 5.0)
        npt.assert_allclose(result.z[-1], 6.0)

    def test_custom_unit(self):
        """Custom unit is stored."""
        result = line3D((0, 0, 0), (1, 1, 1), unit="um")
        assert result.unit == "um"


class TestSlice3D:
    """Tests for slice3D function."""

    def test_xy_plane(self):
        """XY plane slice has constant z."""
        result = slice3D(plane="xy", max1=10.0, max2=10.0, slice_value=5.0, num_points=3)
        npt.assert_allclose(result.z, 5.0)
        assert result.x.shape == (3, 3)

    def test_xz_plane(self):
        """XZ plane slice has constant y."""
        result = slice3D(plane="xz", max1=10.0, max2=10.0, slice_value=5.0, num_points=3)
        npt.assert_allclose(result.y, 5.0)
        assert result.x.shape == (3, 3)

    def test_yz_plane(self):
        """YZ plane slice has constant x."""
        result = slice3D(plane="yz", max1=10.0, max2=10.0, slice_value=5.0, num_points=3)
        npt.assert_allclose(result.x, 5.0)
        assert result.y.shape == (3, 3)

    def test_custom_plane(self):
        """Custom plane with specified arrays."""
        custom_x = np.array([1.0, 2.0])
        custom_y = np.array([3.0, 4.0])
        custom_z = np.array([5.0, 6.0])
        result = slice3D(
            plane="custom",
            custom_x=custom_x,
            custom_y=custom_y,
            custom_z=custom_z,
        )
        npt.assert_allclose(result.x, custom_x)
        npt.assert_allclose(result.y, custom_y)
        npt.assert_allclose(result.z, custom_z)

    def test_invalid_plane_raises(self):
        """Invalid plane name raises exception."""
        with pytest.raises(Exception):
            slice3D(plane="invalid")

    def test_custom_unit(self):
        """Custom unit is stored."""
        result = slice3D(plane="xy", max1=10.0, max2=10.0, unit="m")
        assert result.unit == "m"


class TestPoint3D:
    """Tests for point3D function."""

    def test_returns_point_array3(self):
        """point3D returns Point_Array3."""
        from pymagnet.utils._vector_structs import Point_Array3

        result = point3D((1, 2, 3))
        assert isinstance(result, Point_Array3)

    def test_correct_coordinates(self):
        """Point has correct coordinates."""
        result = point3D((1.5, 2.5, 3.5))
        npt.assert_allclose(result.x, 1.5)
        npt.assert_allclose(result.y, 2.5)
        npt.assert_allclose(result.z, 3.5)

    def test_custom_unit(self):
        """Custom unit is stored."""
        result = point3D((1, 2, 3), unit="cm")
        assert result.unit == "cm"


class TestPlane3D:
    """Tests for plane3D function."""

    def test_returns_point_array3(self):
        """plane3D returns Point_Array3."""
        from pymagnet.utils._vector_structs import Point_Array3

        origin = np.array([0, 0, 0])
        point_a = np.array([1, 0, 0])
        point_b = np.array([0, 1, 0])
        result = plane3D(origin, point_a, point_b, num_points=5)
        assert isinstance(result, Point_Array3)

    def test_shape(self):
        """Plane has correct shape."""
        origin = np.array([0, 0, 0])
        point_a = np.array([1, 0, 0])
        point_b = np.array([0, 1, 0])
        result = plane3D(origin, point_a, point_b, num_points=10)
        assert result.x.shape == (10, 10)

    def test_plane_contains_origin(self):
        """Plane contains the origin point."""
        origin = np.array([1, 2, 3])
        point_a = np.array([2, 2, 3])
        point_b = np.array([1, 3, 3])
        result = plane3D(origin, point_a, point_b, num_points=3)
        # Origin should be at corner (0, 0)
        npt.assert_allclose(result.x[0, 0], 1.0)
        npt.assert_allclose(result.y[0, 0], 2.0)
        npt.assert_allclose(result.z[0, 0], 3.0)


class TestGetField3D:
    """Tests for get_field_3D function."""

    def test_returns_field3(self):
        """get_field_3D returns Field3."""
        from pymagnet.utils._vector_structs import Field3

        points = grid3D(10.0, 10.0, 10.0, num_points=3)
        result = get_field_3D(points)
        assert isinstance(result, Field3)

    def test_empty_registry_returns_zero_field(self):
        """With no magnets, field is zero everywhere."""
        points = grid3D(10.0, 10.0, 10.0, num_points=3)
        result = get_field_3D(points)
        npt.assert_allclose(result.x, 0.0)
        npt.assert_allclose(result.y, 0.0)
        npt.assert_allclose(result.z, 0.0)

    def test_single_magnet_field(self):
        """Single magnet produces non-zero field."""
        from pymagnet import magnets

        _ = magnets.Prism(width=5.0, depth=5.0, height=10.0, Jr=1.0)
        points = slice3D(plane="xz", max1=20.0, max2=20.0, slice_value=0.0, num_points=5)
        result = get_field_3D(points)

        # Field should have non-zero values
        assert np.any(result.x != 0.0) or np.any(result.y != 0.0) or np.any(result.z != 0.0)

    def test_norm_computed(self):
        """Field norm is computed."""
        from pymagnet import magnets

        _ = magnets.Prism(width=5.0, depth=5.0, height=10.0, Jr=1.0)
        points = slice3D(plane="xz", max1=20.0, max2=20.0, slice_value=0.0, num_points=5)
        result = get_field_3D(points)

        assert result.n.shape == result.x.shape


class TestAllocateFieldArray3:
    """Tests for _allocate_field_array3 function."""

    def test_scalar_inputs(self):
        """Handles scalar inputs."""
        result = _allocate_field_array3(1.0, 2.0, 3.0)
        assert result.x.shape == (1,)
        assert result.y.shape == (1,)
        assert result.z.shape == (1,)

    def test_1d_array_inputs(self):
        """Handles 1D array inputs."""
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.0, 0.0, 0.0])
        z = np.array([0.0, 0.0, 0.0])
        result = _allocate_field_array3(x, y, z)
        assert result.x.shape == (3,)

    def test_2d_array_inputs(self):
        """Handles 2D planar inputs."""
        x, y = np.mgrid[0:5, 0:5]
        z = np.zeros_like(x)
        result = _allocate_field_array3(x, y, z)
        assert result.x.shape == (5, 5)

    def test_3d_array_inputs(self):
        """Handles 3D grid inputs."""
        x, y, z = np.mgrid[0:3, 0:3, 0:3]
        result = _allocate_field_array3(x, y, z)
        assert result.x.shape == (3, 3, 3)

    def test_initializes_to_zero(self):
        """Allocated array is initialized to zero."""
        x = np.array([1.0, 2.0])
        y = np.array([3.0, 4.0])
        z = np.array([5.0, 6.0])
        result = _allocate_field_array3(x, y, z)
        npt.assert_allclose(result.x, 0.0)
        npt.assert_allclose(result.y, 0.0)
        npt.assert_allclose(result.z, 0.0)


class TestGetMaxArray:
    """Tests for _get_max_array function."""

    def test_finds_largest(self):
        """Finds largest array in list."""
        arrays = [np.array([1, 2]), np.array([1, 2, 3, 4, 5]), np.array([1])]
        max_shape, max_size = _get_max_array(arrays)
        assert max_size == 5
        assert max_shape == (5,)

    def test_2d_arrays(self):
        """Handles 2D arrays."""
        arrays = [np.ones((3, 3)), np.ones((5, 5)), np.ones((2, 2))]
        max_shape, max_size = _get_max_array(arrays)
        assert max_size == 25
        assert max_shape == (5, 5)


class TestTileArrays:
    """Tests for _tile_arrays function."""

    def test_tiles_scalar(self):
        """Tiles scalar to match array."""
        x, y, z = _tile_arrays(1.0, np.array([1, 2, 3]), 0.0)
        assert x.shape == (3,)
        assert z.shape == (3,)

    def test_no_change_for_equal_sizes(self):
        """Arrays of equal size unchanged."""
        x_in = np.array([1.0, 2.0, 3.0])
        y_in = np.array([4.0, 5.0, 6.0])
        z_in = np.array([7.0, 8.0, 9.0])
        x, y, z = _tile_arrays(x_in, y_in, z_in)
        npt.assert_allclose(x, x_in)
        npt.assert_allclose(y, y_in)
        npt.assert_allclose(z, z_in)
