# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Shared pytest fixtures for pymagnet test suite.

This module provides reusable fixtures for:
- Magnet instances (2D and 3D)
- Point grids and arrays
- Triangle geometries
- Registry state management
"""

import numpy as np
import pytest

import pymagnet
from pymagnet import magnets, reset
from pymagnet.utils import Point2, Point3
from pymagnet.utils._quaternion import Quaternion, q_angle_from_axis
from pymagnet.utils._vector_structs import Field2, Field3, Point_Array2, Point_Array3
from pymagnet.utils.global_const import PI


# ==================== Session/Module Setup ====================


@pytest.fixture(autouse=True)
def reset_magnet_registry():
    """Reset magnet registry before and after each test."""
    reset()
    yield
    reset()


# ==================== 2D Magnet Fixtures ====================


@pytest.fixture
def rectangle_default():
    """Default rectangle magnet: 20x40mm, Jr=1.0T, center at origin."""
    return magnets.Rectangle(width=20.0, height=40.0, Jr=1.0)


@pytest.fixture
def rectangle_rotated():
    """Rectangle rotated 45 degrees."""
    return magnets.Rectangle(width=20.0, height=40.0, Jr=1.0, alpha=45.0)


@pytest.fixture
def rectangle_offset():
    """Rectangle offset from origin."""
    return magnets.Rectangle(width=20.0, height=40.0, Jr=1.0, center=(10.0, 5.0))


@pytest.fixture
def square_default():
    """Default square magnet: 20x20mm, Jr=1.0T."""
    return magnets.Square(width=20.0, Jr=1.0)


@pytest.fixture
def circle_default():
    """Default circle magnet: radius=10mm, Jr=1.0T."""
    return magnets.Circle(radius=10.0, Jr=1.0)


# ==================== 3D Magnet Fixtures ====================


@pytest.fixture
def prism_default():
    """Default prism magnet: 10x20x30mm, Jr=1.0T."""
    return magnets.Prism(width=10.0, depth=20.0, height=30.0, Jr=1.0)


@pytest.fixture
def prism_rotated():
    """Prism with rotation about all axes."""
    return magnets.Prism(
        width=10.0,
        depth=20.0,
        height=30.0,
        Jr=1.0,
        alpha=30.0,
        beta=45.0,
        gamma=60.0,
    )


@pytest.fixture
def cube_default():
    """Default cube magnet: 10mm side, Jr=1.0T."""
    return magnets.Cube(width=10.0, Jr=1.0)


@pytest.fixture
def cylinder_default():
    """Default cylinder magnet: radius=10mm, length=20mm, Jr=1.0T."""
    return magnets.Cylinder(radius=10.0, length=20.0, Jr=1.0)


@pytest.fixture
def sphere_default():
    """Default sphere magnet: radius=10mm, Jr=1.0T."""
    return magnets.Sphere(radius=10.0, Jr=1.0)


# ==================== Point Fixtures ====================


@pytest.fixture
def point2_origin():
    """Point2 at origin."""
    return Point2(0.0, 0.0)


@pytest.fixture
def point2_unit_x():
    """Point2 at (1, 0)."""
    return Point2(1.0, 0.0)


@pytest.fixture
def point2_unit_y():
    """Point2 at (0, 1)."""
    return Point2(0.0, 1.0)


@pytest.fixture
def point3_origin():
    """Point3 at origin."""
    return Point3(0.0, 0.0, 0.0)


@pytest.fixture
def point3_unit_x():
    """Point3 at (1, 0, 0)."""
    return Point3(1.0, 0.0, 0.0)


@pytest.fixture
def point3_unit_y():
    """Point3 at (0, 1, 0)."""
    return Point3(0.0, 1.0, 0.0)


@pytest.fixture
def point3_unit_z():
    """Point3 at (0, 0, 1)."""
    return Point3(0.0, 0.0, 1.0)


# ==================== Point Array Fixtures ====================


@pytest.fixture
def grid_2d_small():
    """Small 2D grid for testing: 10x10 points."""
    x, y = np.mgrid[-1:1:10j, -1:1:10j]
    return Point_Array2(x, y, unit="mm")


@pytest.fixture
def grid_3d_small():
    """Small 3D grid for testing: 5x5x5 points."""
    x, y, z = np.mgrid[-1:1:5j, -1:1:5j, -1:1:5j]
    return Point_Array3(x, y, z, unit="mm")


@pytest.fixture
def line_along_z():
    """Line of points along z-axis."""
    z = np.linspace(-10, 10, 50)
    return Point_Array3(np.zeros_like(z), np.zeros_like(z), z, unit="mm")


# ==================== Quaternion Fixtures ====================


@pytest.fixture
def identity_quaternion():
    """Identity quaternion (no rotation)."""
    return Quaternion()


@pytest.fixture
def rotation_90_z():
    """90 degree rotation about z-axis."""
    return q_angle_from_axis(PI / 2, (0, 0, 1))


@pytest.fixture
def rotation_90_x():
    """90 degree rotation about x-axis."""
    return q_angle_from_axis(PI / 2, (1, 0, 0))


@pytest.fixture
def rotation_90_y():
    """90 degree rotation about y-axis."""
    return q_angle_from_axis(PI / 2, (0, 1, 0))


@pytest.fixture
def rotation_180_z():
    """180 degree rotation about z-axis."""
    return q_angle_from_axis(PI, (0, 0, 1))


# ==================== Triangle Fixtures ====================


@pytest.fixture
def equilateral_triangle():
    """Equilateral triangle in xz plane with unit side length."""
    return np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.0, np.sqrt(3) / 2]], dtype=np.float64
    )


@pytest.fixture
def right_angle_triangle():
    """Right-angled triangle in xy plane."""
    return np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64
    )


@pytest.fixture
def right_angle_triangle_xz():
    """Right-angled triangle in xz plane."""
    return np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64
    )


@pytest.fixture
def degenerate_triangle():
    """Degenerate (collinear points) triangle."""
    return np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64
    )


@pytest.fixture
def triangle_arbitrary():
    """Arbitrary triangle not aligned to any axis."""
    return np.array(
        [[1.0, 2.0, 3.0], [4.0, 1.0, 2.0], [2.0, 5.0, 1.0]], dtype=np.float64
    )


@pytest.fixture
def triangle_345():
    """3-4-5 right triangle in xz plane."""
    return np.array(
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, 4.0]], dtype=np.float64
    )


# ==================== STL Test Files ====================


@pytest.fixture
def stl_cube_path(tmp_path):
    """Create a simple cube STL file for testing."""
    from stl import mesh as stl_mesh

    # Define vertices of a unit cube
    vertices = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    )

    # Define faces (12 triangles for 6 faces)
    faces = np.array(
        [
            [0, 3, 1],
            [1, 3, 2],  # bottom
            [4, 5, 7],
            [5, 6, 7],  # top
            [0, 1, 5],
            [0, 5, 4],  # front
            [2, 3, 7],
            [2, 7, 6],  # back
            [0, 4, 7],
            [0, 7, 3],  # left
            [1, 2, 6],
            [1, 6, 5],  # right
        ]
    )

    cube = stl_mesh.Mesh(np.zeros(faces.shape[0], dtype=stl_mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            cube.vectors[i][j] = vertices[f[j]]

    stl_path = tmp_path / "test_cube.stl"
    cube.save(str(stl_path))
    return stl_path


# ==================== Numerical Tolerances ====================


@pytest.fixture
def float_tolerance():
    """Standard floating point comparison tolerance."""
    return {"rtol": 1e-5, "atol": 1e-8}


@pytest.fixture
def loose_tolerance():
    """Looser tolerance for numerical algorithms."""
    return {"rtol": 1e-3, "atol": 1e-6}


# ==================== Parametrized Values ====================

MAGNETIZATION_VALUES = [0.1, 0.5, 1.0, -1.0, 2.5]
ROTATION_ANGLES = [0, 30, 45, 90, 180, 270]
POSITION_OFFSETS_2D = [(0, 0), (1, 0), (0, 1), (1, 1), (-1, -1)]
POSITION_OFFSETS_3D = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]
SI_LENGTH_UNITS = ["m", "cm", "mm", "um"]
SI_FIELD_UNITS = ["T", "mT", "uT", "nT"]
