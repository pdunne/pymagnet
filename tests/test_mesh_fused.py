# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2026 Peter Dunne
"""Tests for the multi-magnet fused mesh field kernel.

Verifies that ``get_total_field_mesh`` produces results numerically identical
to summing per-magnet ``get_field()`` calls (to float64 precision), and that
the distance-cutoff variant is accurate when the cutoff encompasses all
contributing triangles.
"""

import os

import numpy as np
import pytest

import pymagnet as pm
from pymagnet.magnets import Mesh, get_total_field_mesh

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_STL_DIR = os.path.join(
    os.path.dirname(__file__), "..", "examples", "scripts", "stl_magnets", "stl"
)


def _stl(name):
    return os.path.join(_STL_DIR, name)


_CUBE_STL = _stl("cube.stl")
_HAS_CUBE = os.path.exists(_CUBE_STL)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def reset_registry():
    """Ensure a clean magnet registry for each test."""
    pm.reset()
    yield
    pm.reset()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _HAS_CUBE, reason="cube.stl not found")
def test_single_magnet_fused_matches_get_field():
    """get_total_field_mesh with one magnet == Mesh.get_field()."""
    m = Mesh(_CUBE_STL, Jr=1.0)
    x = np.linspace(-15, 15, 10)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    # get_field() returns (Bx, By, Bz) tuple
    Bx_ref, By_ref, Bz_ref = m.get_field(X, Y, Z)
    B_fused = get_total_field_mesh([m], X, Y, Z)

    assert np.allclose(B_fused.x, Bx_ref, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.y, By_ref, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.z, Bz_ref, rtol=1e-10, atol=1e-15)


@pytest.mark.skipif(not _HAS_CUBE, reason="cube.stl not found")
def test_two_magnets_fused_matches_sum():
    """get_total_field_mesh([m1, m2]) == m1.get_field() + m2.get_field()."""
    m1 = Mesh(_CUBE_STL, Jr=1.0, center=[0.0, 0.0, 0.0])
    m2 = Mesh(_CUBE_STL, Jr=0.8, center=[0.0, 0.0, 30.0])

    x = np.linspace(-20, 20, 12)
    z = np.linspace(-10, 40, 12)
    X, Y, Z = np.meshgrid(x, x, z, indexing="ij")

    Bx1, By1, Bz1 = m1.get_field(X, Y, Z)
    Bx2, By2, Bz2 = m2.get_field(X, Y, Z)
    B_fused = get_total_field_mesh([m1, m2], X, Y, Z)

    assert np.allclose(B_fused.x, Bx1 + Bx2, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.y, By1 + By2, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.z, Bz1 + Bz2, rtol=1e-10, atol=1e-15)


@pytest.mark.skipif(not _HAS_CUBE, reason="cube.stl not found")
def test_three_magnets_fused_matches_sum():
    """Three-magnet fused == sum of three individual calls."""
    m1 = Mesh(_CUBE_STL, Jr=1.0, center=[0.0, 0.0,  0.0])
    m2 = Mesh(_CUBE_STL, Jr=1.0, center=[20.0, 0.0,  0.0])
    m3 = Mesh(_CUBE_STL, Jr=0.5, center=[0.0, 20.0,  0.0])

    x = np.linspace(-25, 25, 10)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    Bx1, By1, Bz1 = m1.get_field(X, Y, Z)
    Bx2, By2, Bz2 = m2.get_field(X, Y, Z)
    Bx3, By3, Bz3 = m3.get_field(X, Y, Z)
    B_fused = get_total_field_mesh([m1, m2, m3], X, Y, Z)

    assert np.allclose(B_fused.x, Bx1 + Bx2 + Bx3, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.y, By1 + By2 + By3, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.z, Bz1 + Bz2 + Bz3, rtol=1e-10, atol=1e-15)


@pytest.mark.skipif(not _HAS_CUBE, reason="cube.stl not found")
def test_fused_with_generous_rcut_matches_sum():
    """Distance cutoff larger than the mesh extent gives the exact answer."""
    m1 = Mesh(_CUBE_STL, Jr=1.0, center=[0.0, 0.0,  0.0])
    m2 = Mesh(_CUBE_STL, Jr=1.0, center=[0.0, 0.0, 30.0])

    x = np.linspace(-20, 20, 8)
    z = np.linspace(-10, 40, 8)
    X, Y, Z = np.meshgrid(x, x, z, indexing="ij")

    Bx1, By1, Bz1 = m1.get_field(X, Y, Z)
    Bx2, By2, Bz2 = m2.get_field(X, Y, Z)
    # r_cut=200 mm encompasses the whole scene (cube is ~5 mm)
    B_fused = get_total_field_mesh([m1, m2], X, Y, Z, r_cut=200.0)

    assert np.allclose(B_fused.x, Bx1 + Bx2, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.y, By1 + By2, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.z, Bz1 + Bz2, rtol=1e-10, atol=1e-15)


@pytest.mark.skipif(not _HAS_CUBE, reason="cube.stl not found")
def test_mesh_instances_iterable():
    """Passing Mesh.instances (WeakSet) works as the meshes argument."""
    m1 = Mesh(_CUBE_STL, Jr=1.0, center=[0.0, 0.0,  0.0])
    m2 = Mesh(_CUBE_STL, Jr=1.0, center=[0.0, 0.0, 20.0])

    x = np.linspace(-10, 10, 6)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    Bx1, By1, Bz1 = m1.get_field(X, Y, Z)
    Bx2, By2, Bz2 = m2.get_field(X, Y, Z)

    # Use the class-level registry directly
    B_fused = get_total_field_mesh(Mesh.instances, X, Y, Z)

    assert np.allclose(B_fused.x, Bx1 + Bx2, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.y, By1 + By2, rtol=1e-10, atol=1e-15)
    assert np.allclose(B_fused.z, Bz1 + Bz2, rtol=1e-10, atol=1e-15)
