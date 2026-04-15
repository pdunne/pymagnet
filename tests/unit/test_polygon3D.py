# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for Mesh class (3D polygon magnets from STL files)."""

import os
import tempfile

import numpy as np
import numpy.testing as npt
import pytest

# Try to import stl for creating test meshes
try:
    from stl import mesh as stl_mesh

    HAS_STL = True
except ImportError:
    HAS_STL = False


@pytest.fixture
def cube_stl_file():
    """Create a temporary STL file with a simple cube."""
    if not HAS_STL:
        pytest.skip("numpy-stl not available")

    # Define the 8 vertices of a unit cube centered at origin
    vertices = np.array(
        [
            [-0.5, -0.5, -0.5],
            [0.5, -0.5, -0.5],
            [0.5, 0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
            [0.5, -0.5, 0.5],
            [0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
        ]
    )

    # Define the 12 triangles composing the cube
    faces = np.array(
        [
            [0, 3, 1],
            [1, 3, 2],  # Bottom
            [0, 4, 7],
            [0, 7, 3],  # Left
            [4, 5, 6],
            [4, 6, 7],  # Top
            [5, 1, 2],
            [5, 2, 6],  # Right
            [2, 3, 6],
            [3, 7, 6],  # Back
            [0, 1, 5],
            [0, 5, 4],  # Front
        ]
    )

    # Create the mesh
    cube = stl_mesh.Mesh(np.zeros(faces.shape[0], dtype=stl_mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            cube.vectors[i][j] = vertices[f[j], :]

    # Write to temporary file
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        cube.save(f.name)
        filepath = f.name

    yield filepath

    # Cleanup
    if os.path.exists(filepath):
        os.unlink(filepath)


@pytest.fixture
def tetrahedron_stl_file():
    """Create a temporary STL file with a simple tetrahedron."""
    if not HAS_STL:
        pytest.skip("numpy-stl not available")

    # Define the 4 vertices of a tetrahedron
    vertices = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0.5, np.sqrt(3) / 2, 0],
            [0.5, np.sqrt(3) / 6, np.sqrt(2 / 3)],
        ]
    )

    # Define the 4 triangular faces
    faces = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]])

    # Create the mesh
    tetra = stl_mesh.Mesh(np.zeros(faces.shape[0], dtype=stl_mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            tetra.vectors[i][j] = vertices[f[j], :]

    # Write to temporary file
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        tetra.save(f.name)
        filepath = f.name

    yield filepath

    # Cleanup
    if os.path.exists(filepath):
        os.unlink(filepath)


@pytest.mark.skipif(not HAS_STL, reason="numpy-stl not available")
class TestMeshInit:
    """Tests for Mesh initialization."""

    def test_mesh_loads_stl(self, cube_stl_file):
        """Mesh loads STL file successfully."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        assert m is not None
        assert m.Jr == 1.0

    def test_mesh_has_vectors(self, cube_stl_file):
        """Mesh has mesh_vectors attribute."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        assert hasattr(m, "mesh_vectors")
        assert m.mesh_vectors.shape[0] == 12  # Cube has 12 triangles
        assert m.mesh_vectors.shape[1] == 3  # Each triangle has 3 vertices
        assert m.mesh_vectors.shape[2] == 3  # Each vertex is 3D

    def test_mesh_has_normals(self, cube_stl_file):
        """Mesh has mesh_normals attribute with unit normals."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        assert hasattr(m, "mesh_normals")
        assert m.mesh_normals.shape[0] == 12
        # Check normals are unit vectors
        norms = np.linalg.norm(m.mesh_normals, axis=1)
        npt.assert_allclose(norms, 1.0, rtol=1e-10)

    def test_mesh_has_volume(self, cube_stl_file):
        """Mesh has volume attribute."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        assert hasattr(m, "volume")
        # Unit cube has volume 1.0
        npt.assert_allclose(abs(m.volume), 1.0, rtol=0.01)

    def test_mesh_has_centroid(self, cube_stl_file):
        """Mesh has centroid attribute."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        assert hasattr(m, "centroid")
        # Centroid should be at origin for centered cube
        npt.assert_allclose(m.centroid, [0.0, 0.0, 0.0], atol=1e-10)

    def test_mesh_jr_components(self, cube_stl_file):
        """Mesh has Jx, Jy, Jz components."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        # Default theta=0, phi=90 means magnetization along z
        assert hasattr(m, "Jx")
        assert hasattr(m, "Jy")
        assert hasattr(m, "Jz")

    def test_mesh_get_jr(self, cube_stl_file):
        """get_Jr returns magnetization vector."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        J = m.get_Jr()
        assert len(J) == 3

    def test_mesh_custom_center(self, cube_stl_file):
        """Mesh can be created with custom center."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0, center=(5.0, 6.0, 7.0))
        center = m.get_center()
        npt.assert_allclose(center, [5.0, 6.0, 7.0])

    def test_mesh_scale(self, cube_stl_file):
        """Mesh can be scaled."""
        from pymagnet import magnets

        # Scale by 10
        m = magnets.Mesh(cube_stl_file, Jr=1.0, mesh_scale=10.0)
        # Volume should be 10^3 = 1000
        npt.assert_allclose(abs(m.volume), 1000.0, rtol=0.01)

    def test_mesh_str_repr(self, cube_stl_file):
        """String representation contains class name."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        str_rep = str(m)
        assert "Mesh" in str_rep

    def test_mesh_mag_type(self, cube_stl_file):
        """Mesh has correct mag_type."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        assert m.mag_type == "Mesh"


@pytest.mark.skipif(not HAS_STL, reason="numpy-stl not available")
class TestMeshField:
    """Tests for Mesh field calculations."""

    def test_field_returns_tuple(self, cube_stl_file):
        """get_field returns Bx, By, Bz tuple."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        x = np.array([0.0, 0.0, 2.0])
        y = np.array([0.0, 0.0, 0.0])
        z = np.array([2.0, 0.0, 0.0])
        Bx, By, Bz = m.get_field(x, y, z)
        assert Bx.shape == (3,)
        assert By.shape == (3,)
        assert Bz.shape == (3,)

    def test_field_finite(self, cube_stl_file):
        """Field values are finite away from magnet."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([5.0])
        Bx, By, Bz = m.get_field(x, y, z)
        assert np.isfinite(Bx[0])
        assert np.isfinite(By[0])
        assert np.isfinite(Bz[0])


@pytest.mark.skipif(not HAS_STL, reason="numpy-stl not available")
class TestChargeSheetFunctions:
    """Tests for charge sheet calculation functions."""

    def test_charge_sheet_x_finite(self, cube_stl_file):
        """_charge_sheet_x returns finite values."""
        from pymagnet.magnets._polygon3D import _charge_sheet_x

        result = _charge_sheet_x(1.0, 1.0, 1.0, 2.0, 0.0, 0.0)
        assert np.isfinite(result)

    def test_charge_sheet_y_finite(self, cube_stl_file):
        """_charge_sheet_y returns finite values."""
        from pymagnet.magnets._polygon3D import _charge_sheet_y

        result = _charge_sheet_y(1.0, 1.0, 1.0, 0.0, 2.0, 0.0)
        assert np.isfinite(result)

    def test_charge_sheet_z_finite(self, cube_stl_file):
        """_charge_sheet_z returns finite values."""
        from pymagnet.magnets._polygon3D import _charge_sheet_z

        result = _charge_sheet_z(1.0, 1.0, 1.0, 0.0, 0.0, 2.0)
        assert np.isfinite(result)

    def test_charge_sheet_y_at_zero(self, cube_stl_file):
        """_charge_sheet_y returns 0 at y=0 (singularity handled)."""
        from pymagnet.magnets._polygon3D import _charge_sheet_y

        result = _charge_sheet_y(1.0, 1.0, 1.0, 1.0, 0.0, 1.0)
        assert result == 0.0

    def test_charge_sheet_vectorized(self):
        """Charge sheet functions work with array inputs."""
        from pymagnet.magnets._polygon3D import _charge_sheet_x

        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.0, 0.0, 0.0])
        z = np.array([1.0, 1.0, 1.0])
        result = _charge_sheet_x(1.0, 1.0, 1.0, x, y, z)
        assert result.shape == (3,)


@pytest.mark.skipif(not HAS_STL, reason="numpy-stl not available")
class TestCalcBTriangle:
    """Tests for single triangle field calculation."""

    def test_calcB_triangle_returns_components(self, cube_stl_file):
        """calcB_triangle returns field components and rotation info."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0)
        triangle = m.mesh_vectors[0]
        Jr = m.Jnorm[0]

        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([2.0])

        result = m.calcB_triangle(triangle, Jr, x, y, z, 0)
        assert len(result) == 6  # Bx, By, Bz, rotated_triangle, offset, rotation


@pytest.mark.skipif(not HAS_STL, reason="numpy-stl not available")
class TestMeshRotation:
    """Tests for mesh rotation functionality."""

    def test_mesh_with_rotation(self, cube_stl_file):
        """Mesh can be created with rotation angles."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0, alpha=45.0, beta=0.0, gamma=0.0)
        orientation = m.get_orientation()
        npt.assert_allclose(orientation, [45.0, 0.0, 0.0])

    def test_mesh_normals_unit_after_rotation(self, cube_stl_file):
        """Mesh normals remain unit vectors after rotation."""
        from pymagnet import magnets

        m = magnets.Mesh(cube_stl_file, Jr=1.0, alpha=30.0, beta=45.0, gamma=60.0)
        norms = np.linalg.norm(m.mesh_normals, axis=1)
        npt.assert_allclose(norms, 1.0, rtol=1e-10)
