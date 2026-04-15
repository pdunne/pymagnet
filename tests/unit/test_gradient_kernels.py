# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests validating numba gradient kernels against np.gradient."""

import numpy as np
import numpy.testing as npt

from pymagnet.utils._gradient_kernels import _gradient_2d, _gradient_3d


class TestKernelVsNumpy2D:
    """Validate _gradient_2d matches np.gradient on 2D arrays."""

    def test_random_square_array(self):
        """Random 20x20 array matches np.gradient."""
        rng = np.random.default_rng(42)
        F = rng.standard_normal((20, 20))
        dx, dy = 0.5, 0.3

        numba_dx, numba_dy = _gradient_2d(np.ascontiguousarray(F), dx, dy)
        np_dx, np_dy = np.gradient(F, dx, dy)

        npt.assert_allclose(numba_dx, np_dx, atol=1e-12)
        npt.assert_allclose(numba_dy, np_dy, atol=1e-12)

    def test_random_rectangular_array(self):
        """Non-square 15x30 array matches np.gradient."""
        rng = np.random.default_rng(123)
        F = rng.standard_normal((15, 30))
        dx, dy = 1.0, 0.7

        numba_dx, numba_dy = _gradient_2d(np.ascontiguousarray(F), dx, dy)
        np_dx, np_dy = np.gradient(F, dx, dy)

        npt.assert_allclose(numba_dx, np_dx, atol=1e-12)
        npt.assert_allclose(numba_dy, np_dy, atol=1e-12)

    def test_small_3x3_array(self):
        """Minimum-size 3x3 array matches np.gradient."""
        F = np.array([[1.0, 2.0, 4.0], [3.0, 5.0, 7.0], [6.0, 8.0, 9.0]])
        dx, dy = 1.0, 1.0

        numba_dx, numba_dy = _gradient_2d(np.ascontiguousarray(F), dx, dy)
        np_dx, np_dy = np.gradient(F, dx, dy)

        npt.assert_allclose(numba_dx, np_dx, atol=1e-12)
        npt.assert_allclose(numba_dy, np_dy, atol=1e-12)

    def test_different_spacings(self):
        """Works with very different dx and dy values."""
        rng = np.random.default_rng(99)
        F = rng.standard_normal((10, 10))
        dx, dy = 0.001, 100.0

        numba_dx, numba_dy = _gradient_2d(np.ascontiguousarray(F), dx, dy)
        np_dx, np_dy = np.gradient(F, dx, dy)

        npt.assert_allclose(numba_dx, np_dx, atol=1e-8)
        npt.assert_allclose(numba_dy, np_dy, atol=1e-8)


class TestKernelVsNumpy3D:
    """Validate _gradient_3d matches np.gradient on 3D arrays."""

    def test_random_cube_array(self):
        """Random 10x10x10 array matches np.gradient."""
        rng = np.random.default_rng(42)
        F = rng.standard_normal((10, 10, 10))
        dx, dy, dz = 0.5, 0.3, 0.8

        numba_dx, numba_dy, numba_dz = _gradient_3d(
            np.ascontiguousarray(F), dx, dy, dz
        )
        np_dx, np_dy, np_dz = np.gradient(F, dx, dy, dz)

        npt.assert_allclose(numba_dx, np_dx, atol=1e-12)
        npt.assert_allclose(numba_dy, np_dy, atol=1e-12)
        npt.assert_allclose(numba_dz, np_dz, atol=1e-12)

    def test_random_rectangular_array(self):
        """Non-cube 8x12x6 array matches np.gradient."""
        rng = np.random.default_rng(123)
        F = rng.standard_normal((8, 12, 6))
        dx, dy, dz = 1.0, 0.5, 2.0

        numba_dx, numba_dy, numba_dz = _gradient_3d(
            np.ascontiguousarray(F), dx, dy, dz
        )
        np_dx, np_dy, np_dz = np.gradient(F, dx, dy, dz)

        npt.assert_allclose(numba_dx, np_dx, atol=1e-12)
        npt.assert_allclose(numba_dy, np_dy, atol=1e-12)
        npt.assert_allclose(numba_dz, np_dz, atol=1e-12)

    def test_small_3x3x3_array(self):
        """Minimum-size 3x3x3 array matches np.gradient."""
        rng = np.random.default_rng(7)
        F = rng.standard_normal((3, 3, 3))
        dx, dy, dz = 1.0, 1.0, 1.0

        numba_dx, numba_dy, numba_dz = _gradient_3d(
            np.ascontiguousarray(F), dx, dy, dz
        )
        np_dx, np_dy, np_dz = np.gradient(F, dx, dy, dz)

        npt.assert_allclose(numba_dx, np_dx, atol=1e-12)
        npt.assert_allclose(numba_dy, np_dy, atol=1e-12)
        npt.assert_allclose(numba_dz, np_dz, atol=1e-12)
