# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for vector and point array structures."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.utils._vector_structs import (
    Point_Array1,
    Point_Array2,
    Point_Array3,
    Field1,
    Field2,
    Field3,
)


class TestPointArray1:
    """Tests for Point_Array1 class."""

    def test_init_scalar(self):
        """Initialize with scalar value."""
        p = Point_Array1(5.0)
        npt.assert_allclose(p.z, [5.0])

    def test_init_array(self):
        """Initialize with array."""
        p = Point_Array1([1.0, 2.0, 3.0])
        npt.assert_allclose(p.z, [1.0, 2.0, 3.0])

    def test_default_unit(self):
        """Default unit is mm."""
        p = Point_Array1(1.0)
        assert p.unit == "mm"

    def test_custom_unit(self):
        """Custom unit is stored."""
        p = Point_Array1(1.0, unit="cm")
        assert p.unit == "cm"

    def test_invalid_unit_raises(self):
        """Invalid unit raises ValueError."""
        with pytest.raises(ValueError):
            Point_Array1(1.0, unit="invalid")

    def test_get_unit(self):
        """get_unit returns current unit."""
        p = Point_Array1(1.0, unit="m")
        assert p.get_unit() == "m"

    def test_str_repr(self):
        """String representation contains unit."""
        p = Point_Array1(1.0, unit="cm")
        str_rep = str(p)
        assert "cm" in str_rep

    def test_repr_contains_unit(self):
        """repr contains unit."""
        p = Point_Array1(1.0, unit="cm")
        repr_str = repr(p)
        assert "cm" in repr_str


class TestPointArray2:
    """Tests for Point_Array2 class."""

    def test_init_scalars(self):
        """Initialize with scalar values."""
        p = Point_Array2(1.0, 2.0)
        npt.assert_allclose(p.x, [1.0])
        npt.assert_allclose(p.y, [2.0])

    def test_init_arrays(self):
        """Initialize with arrays."""
        p = Point_Array2([1.0, 2.0], [3.0, 4.0])
        npt.assert_allclose(p.x, [1.0, 2.0])
        npt.assert_allclose(p.y, [3.0, 4.0])

    def test_default_unit(self):
        """Default unit is mm."""
        p = Point_Array2(1.0, 2.0)
        assert p.unit == "mm"

    def test_custom_unit(self):
        """Custom unit is stored."""
        p = Point_Array2(1.0, 2.0, unit="um")
        assert p.unit == "um"

    def test_invalid_unit_raises(self):
        """Invalid unit raises ValueError."""
        with pytest.raises(ValueError):
            Point_Array2(1.0, 2.0, unit="invalid")

    def test_get_unit(self):
        """get_unit returns current unit."""
        p = Point_Array2(1.0, 2.0, unit="m")
        assert p.get_unit() == "m"

    def test_str_contains_xy(self):
        """String representation contains x and y info."""
        p = Point_Array2(1.0, 2.0)
        str_rep = str(p)
        assert "x:" in str_rep
        assert "y:" in str_rep


class TestPointArray3:
    """Tests for Point_Array3 class."""

    def test_init_scalars(self):
        """Initialize with scalar values."""
        p = Point_Array3(1.0, 2.0, 3.0)
        npt.assert_allclose(p.x, [1.0])
        npt.assert_allclose(p.y, [2.0])
        npt.assert_allclose(p.z, [3.0])

    def test_init_arrays(self):
        """Initialize with arrays."""
        p = Point_Array3([1.0, 2.0], [3.0, 4.0], [5.0, 6.0])
        npt.assert_allclose(p.x, [1.0, 2.0])
        npt.assert_allclose(p.y, [3.0, 4.0])
        npt.assert_allclose(p.z, [5.0, 6.0])

    def test_inherits_point_array2(self):
        """Point_Array3 is subclass of Point_Array2."""
        p = Point_Array3(1.0, 2.0, 3.0)
        assert isinstance(p, Point_Array2)

    def test_rotate_identity(self):
        """Rotation by zero angles is identity."""
        p = Point_Array3([1.0], [0.0], [0.0])
        rotated = p.rotate(0, 0, 0)
        npt.assert_allclose(rotated.x, p.x, atol=1e-10)
        npt.assert_allclose(rotated.y, p.y, atol=1e-10)
        npt.assert_allclose(rotated.z, p.z, atol=1e-10)

    def test_rotate_90_about_z(self):
        """90° rotation about z moves x to y."""
        p = Point_Array3([1.0], [0.0], [0.0])
        rotated = p.rotate(90, 0, 0)  # alpha is about z
        npt.assert_allclose(rotated.x, [0.0], atol=1e-10)
        npt.assert_allclose(rotated.y, [1.0], atol=1e-10)
        npt.assert_allclose(rotated.z, [0.0], atol=1e-10)

    def test_rotate_preserves_magnitude(self):
        """Rotation preserves vector magnitude."""
        p = Point_Array3([1.0], [2.0], [3.0])
        original_mag = np.sqrt(p.x**2 + p.y**2 + p.z**2)

        rotated = p.rotate(30, 45, 60)
        rotated_mag = np.sqrt(rotated.x**2 + rotated.y**2 + rotated.z**2)

        npt.assert_allclose(rotated_mag, original_mag, rtol=1e-10)

    def test_str_contains_xyz(self):
        """String representation contains x, y, z info."""
        p = Point_Array3(1.0, 2.0, 3.0)
        str_rep = str(p)
        assert "x:" in str_rep
        assert "y:" in str_rep
        assert "z:" in str_rep


class TestField1:
    """Tests for Field1 class."""

    def test_init_scalar(self):
        """Initialize with scalar value."""
        f = Field1(0.5)
        npt.assert_allclose(f.z, [0.5])

    def test_init_array(self):
        """Initialize with array."""
        f = Field1([0.1, 0.2, 0.3])
        npt.assert_allclose(f.z, [0.1, 0.2, 0.3])

    def test_default_unit(self):
        """Default unit is T (Tesla)."""
        f = Field1(0.5)
        assert f.unit == "T"

    def test_custom_unit(self):
        """Custom unit is stored."""
        f = Field1(500.0, unit="mT")
        assert f.unit == "mT"

    def test_invalid_unit_raises(self):
        """Invalid unit raises ValueError."""
        with pytest.raises(ValueError):
            Field1(0.5, unit="invalid")

    def test_str_contains_bz(self):
        """String representation contains Bz."""
        f = Field1(0.5)
        str_rep = str(f)
        assert "Bz" in str_rep


class TestField2:
    """Tests for Field2 class."""

    def test_init_scalars(self):
        """Initialize with scalar values."""
        f = Field2(0.1, 0.2)
        npt.assert_allclose(f.x, [0.1])
        npt.assert_allclose(f.y, [0.2])

    def test_init_arrays(self):
        """Initialize with arrays."""
        f = Field2([0.1, 0.2], [0.3, 0.4])
        npt.assert_allclose(f.x, [0.1, 0.2])
        npt.assert_allclose(f.y, [0.3, 0.4])

    def test_default_unit(self):
        """Default unit is T (Tesla)."""
        f = Field2(0.1, 0.2)
        assert f.unit == "T"

    def test_custom_unit(self):
        """Custom unit is stored."""
        f = Field2(100.0, 200.0, unit="mT")
        assert f.unit == "mT"

    def test_invalid_unit_raises(self):
        """Invalid unit raises ValueError."""
        with pytest.raises(ValueError):
            Field2(0.1, 0.2, unit="invalid")

    def test_has_norm_attribute(self):
        """Field2 has norm attribute n."""
        f = Field2(0.1, 0.2)
        assert hasattr(f, "n")

    def test_calc_norm(self):
        """calc_norm computes vector magnitude."""
        f = Field2([3.0], [4.0])
        f.calc_norm()
        npt.assert_allclose(f.n, [5.0], rtol=1e-10)

    def test_calc_norm_array(self):
        """calc_norm works with arrays."""
        f = Field2([3.0, 0.0], [4.0, 5.0])
        f.calc_norm()
        npt.assert_allclose(f.n, [5.0, 5.0], rtol=1e-10)

    def test_str_contains_bx_by(self):
        """String representation contains Bx, By."""
        f = Field2(0.1, 0.2)
        str_rep = str(f)
        assert "Bx" in str_rep
        assert "By" in str_rep


class TestField3:
    """Tests for Field3 class."""

    def test_init_scalars(self):
        """Initialize with scalar values."""
        f = Field3(0.1, 0.2, 0.3)
        npt.assert_allclose(f.x, [0.1])
        npt.assert_allclose(f.y, [0.2])
        npt.assert_allclose(f.z, [0.3])

    def test_init_arrays(self):
        """Initialize with arrays."""
        f = Field3([0.1, 0.2], [0.3, 0.4], [0.5, 0.6])
        npt.assert_allclose(f.x, [0.1, 0.2])
        npt.assert_allclose(f.y, [0.3, 0.4])
        npt.assert_allclose(f.z, [0.5, 0.6])

    def test_inherits_point_array3(self):
        """Field3 is subclass of Point_Array3."""
        f = Field3(0.1, 0.2, 0.3)
        assert isinstance(f, Point_Array3)

    def test_default_unit(self):
        """Default unit is T (Tesla)."""
        f = Field3(0.1, 0.2, 0.3)
        assert f.unit == "T"

    def test_has_norm_attribute(self):
        """Field3 has norm attribute n."""
        f = Field3(0.1, 0.2, 0.3)
        assert hasattr(f, "n")

    def test_calc_norm(self):
        """calc_norm computes 3D vector magnitude."""
        f = Field3([2.0], [3.0], [6.0])
        f.calc_norm()
        npt.assert_allclose(f.n, [7.0], rtol=1e-10)

    def test_calc_norm_array(self):
        """calc_norm works with arrays."""
        f = Field3([1.0, 2.0], [2.0, 3.0], [2.0, 6.0])
        f.calc_norm()
        npt.assert_allclose(f.n, [3.0, 7.0], rtol=1e-10)

    def test_str_contains_bx_by_bz(self):
        """String representation contains Bx, By, Bz."""
        f = Field3(0.1, 0.2, 0.3)
        str_rep = str(f)
        assert "Bx" in str_rep
        assert "By" in str_rep
        assert "Bz" in str_rep


class TestUnitConversion:
    """Tests for unit conversion in vector structs."""

    def test_field1_unit_attribute(self):
        """Field1 stores unit correctly."""
        f = Field1(0.5, unit="uT")
        assert f.unit == "uT"

    def test_field2_unit_attribute(self):
        """Field2 stores unit correctly."""
        f = Field2(0.1, 0.2, unit="nT")
        assert f.unit == "nT"

    def test_point_array_accepts_valid_length_units(self):
        """Point arrays accept valid length units."""
        valid_units = ["km", "m", "cm", "mm", "um"]
        for unit in valid_units:
            p = Point_Array1(1.0, unit=unit)
            assert p.unit == unit

    def test_field_accepts_valid_field_units(self):
        """Field structs accept valid field units."""
        valid_units = ["T", "mT", "uT", "nT"]
        for unit in valid_units:
            f = Field1(1.0, unit=unit)
            assert f.unit == unit


class TestArrayShapes:
    """Tests for handling various array shapes."""

    def test_point_array2_2d_grid(self):
        """Point_Array2 handles 2D grid arrays."""
        x, y = np.mgrid[0:5, 0:5]
        p = Point_Array2(x, y)
        assert p.x.shape == (5, 5)
        assert p.y.shape == (5, 5)

    def test_point_array3_3d_grid(self):
        """Point_Array3 handles 3D grid arrays."""
        x, y, z = np.mgrid[0:3, 0:3, 0:3]
        p = Point_Array3(x, y, z)
        assert p.x.shape == (3, 3, 3)
        assert p.y.shape == (3, 3, 3)
        assert p.z.shape == (3, 3, 3)

    def test_field2_2d_grid(self):
        """Field2 handles 2D grid arrays."""
        bx = np.ones((5, 5))
        by = np.ones((5, 5))
        f = Field2(bx, by)
        assert f.x.shape == (5, 5)
        assert f.n.shape == (5, 5)

    def test_field3_3d_grid(self):
        """Field3 handles 3D grid arrays."""
        bx = np.ones((3, 3, 3))
        by = np.ones((3, 3, 3))
        bz = np.ones((3, 3, 3))
        f = Field3(bx, by, bz)
        assert f.x.shape == (3, 3, 3)
        assert f.n.shape == (3, 3, 3)
