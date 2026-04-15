# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2026 Peter Dunne
"""Numba-accelerated finite difference gradient kernels.

These replicate the behaviour of ``np.gradient`` (second-order central
differences, first-order forward/backward at edges) but can be called
from ``@njit`` code and benefit from parallel execution via ``prange``.
"""

import numpy as np
from numba import njit, prange


@njit(cache=True, parallel=True)
def _gradient_2d(F, dx, dy):
    """Compute the gradient of a 2D scalar field using finite differences.

    Equivalent to ``np.gradient(F, dx, dy)`` but numba-compatible.

    Args:
        F: 2D array of shape (Nx, Ny).
        dx: Grid spacing along axis 0.
        dy: Grid spacing along axis 1.

    Returns:
        Tuple (dF_dx, dF_dy), each of shape (Nx, Ny).
    """
    Nx, Ny = F.shape
    dF_dx = np.empty_like(F)
    dF_dy = np.empty_like(F)

    for i in prange(Nx):
        for j in range(Ny):
            # Gradient along axis 0 (x)
            if i == 0:
                dF_dx[i, j] = (F[1, j] - F[0, j]) / dx
            elif i == Nx - 1:
                dF_dx[i, j] = (F[Nx - 1, j] - F[Nx - 2, j]) / dx
            else:
                dF_dx[i, j] = (F[i + 1, j] - F[i - 1, j]) / (2.0 * dx)

            # Gradient along axis 1 (y)
            if j == 0:
                dF_dy[i, j] = (F[i, 1] - F[i, 0]) / dy
            elif j == Ny - 1:
                dF_dy[i, j] = (F[i, Ny - 1] - F[i, Ny - 2]) / dy
            else:
                dF_dy[i, j] = (F[i, j + 1] - F[i, j - 1]) / (2.0 * dy)

    return dF_dx, dF_dy


@njit(cache=True, parallel=True)
def _gradient_3d(F, dx, dy, dz):
    """Compute the gradient of a 3D scalar field using finite differences.

    Equivalent to ``np.gradient(F, dx, dy, dz)`` but numba-compatible.

    Args:
        F: 3D array of shape (Nx, Ny, Nz).
        dx: Grid spacing along axis 0.
        dy: Grid spacing along axis 1.
        dz: Grid spacing along axis 2.

    Returns:
        Tuple (dF_dx, dF_dy, dF_dz), each of shape (Nx, Ny, Nz).
    """
    Nx, Ny, Nz = F.shape
    dF_dx = np.empty_like(F)
    dF_dy = np.empty_like(F)
    dF_dz = np.empty_like(F)

    for i in prange(Nx):
        for j in range(Ny):
            for k in range(Nz):
                # Gradient along axis 0 (x)
                if i == 0:
                    dF_dx[i, j, k] = (F[1, j, k] - F[0, j, k]) / dx
                elif i == Nx - 1:
                    dF_dx[i, j, k] = (F[Nx - 1, j, k] - F[Nx - 2, j, k]) / dx
                else:
                    dF_dx[i, j, k] = (F[i + 1, j, k] - F[i - 1, j, k]) / (2.0 * dx)

                # Gradient along axis 1 (y)
                if j == 0:
                    dF_dy[i, j, k] = (F[i, 1, k] - F[i, 0, k]) / dy
                elif j == Ny - 1:
                    dF_dy[i, j, k] = (F[i, Ny - 1, k] - F[i, Ny - 2, k]) / dy
                else:
                    dF_dy[i, j, k] = (F[i, j + 1, k] - F[i, j - 1, k]) / (2.0 * dy)

                # Gradient along axis 2 (z)
                if k == 0:
                    dF_dz[i, j, k] = (F[i, j, 1] - F[i, j, 0]) / dz
                elif k == Nz - 1:
                    dF_dz[i, j, k] = (F[i, j, Nz - 1] - F[i, j, Nz - 2]) / dz
                else:
                    dF_dz[i, j, k] = (F[i, j, k + 1] - F[i, j, k - 1]) / (2.0 * dz)

    return dF_dx, dF_dy, dF_dz
