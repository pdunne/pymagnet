# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Numba-compatible quaternion operations.

This module provides quaternion operations as pure functions operating on
numpy arrays, compatible with Numba's nopython mode. Quaternions are
represented as (4,) arrays in [w, x, y, z] order.

These functions can be used in @njit decorated functions and with prange
for parallel execution.

Example:
    from pymagnet.utils._quaternion_numba import (
        quat_from_axis_angle, quat_rotate_vector, quat_multiply
    )

    # Create rotation quaternion for 90° about z-axis
    q = quat_from_axis_angle(np.pi/2, np.array([0., 0., 1.]))

    # Rotate a vector
    v = np.array([1., 0., 0.])
    v_rotated = quat_rotate_vector(q, v)  # Returns [0, 1, 0]
"""

import numpy as np
from numba import njit

from .global_const import FP_CUTOFF

# =============================================================================
# Core Quaternion Operations
# =============================================================================


@njit(cache=True)
def quat_identity():
    """Return identity quaternion [1, 0, 0, 0].

    Returns:
        ndarray: (4,) identity quaternion
    """
    return np.array([1.0, 0.0, 0.0, 0.0])


@njit(cache=True)
def quat_conjugate(q):
    """Return conjugate of quaternion.

    Args:
        q (ndarray): (4,) quaternion [w, x, y, z]

    Returns:
        ndarray: (4,) conjugate quaternion [w, -x, -y, -z]
    """
    return np.array([q[0], -q[1], -q[2], -q[3]])


@njit(cache=True)
def quat_multiply(q1, q2):
    """Compute Hamilton product of two quaternions.

    Args:
        q1 (ndarray): (4,) first quaternion [w, x, y, z]
        q2 (ndarray): (4,) second quaternion [w, x, y, z]

    Returns:
        ndarray: (4,) product quaternion q1 * q2
    """
    w1, x1, y1, z1 = q1[0], q1[1], q1[2], q1[3]
    w2, x2, y2, z2 = q2[0], q2[1], q2[2], q2[3]

    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


@njit(cache=True)
def quat_multiply3(q1, q2, q3):
    """Compute Hamilton product of three quaternions: q1 * q2 * q3.

    Args:
        q1 (ndarray): (4,) first quaternion
        q2 (ndarray): (4,) second quaternion
        q3 (ndarray): (4,) third quaternion

    Returns:
        ndarray: (4,) product quaternion q1 * q2 * q3
    """
    return quat_multiply(quat_multiply(q1, q2), q3)


# =============================================================================
# Quaternion Creation
# =============================================================================


@njit(cache=True)
def quat_from_axis_angle(theta, axis):
    """Create rotation quaternion from axis-angle representation.

    Args:
        theta (float): rotation angle in radians
        axis (ndarray): (3,) rotation axis (will be normalized)

    Returns:
        ndarray: (4,) unit quaternion [w, x, y, z]
    """
    # Normalize axis
    norm = np.sqrt(axis[0] ** 2 + axis[1] ** 2 + axis[2] ** 2)
    if norm < FP_CUTOFF:
        return quat_identity()

    ax = axis[0] / norm
    ay = axis[1] / norm
    az = axis[2] / norm

    half_theta = theta / 2.0
    s = np.sin(half_theta)
    c = np.cos(half_theta)

    return np.array([c, ax * s, ay * s, az * s])


@njit(cache=True)
def quat_from_euler(alpha, beta, gamma):
    """Create quaternion from Euler angles (z-y-x convention).

    Args:
        alpha (float): rotation about z-axis in radians
        beta (float): rotation about y-axis in radians
        gamma (float): rotation about x-axis in radians

    Returns:
        ndarray: (4,) quaternion [w, x, y, z]
    """
    ca = np.cos(alpha / 2)
    sa = np.sin(alpha / 2)
    cb = np.cos(beta / 2)
    sb = np.sin(beta / 2)
    cg = np.cos(gamma / 2)
    sg = np.sin(gamma / 2)

    w = ca * cb * cg + sa * sb * sg
    x = sa * cb * cg - ca * sb * sg
    y = ca * sb * cg + sa * cb * sg
    z = ca * cb * sg - sa * sb * cg

    return np.array([w, x, y, z])


# =============================================================================
# Vector Rotation
# =============================================================================


@njit(cache=True)
def quat_rotate_vector(q, v):
    """Rotate 3D vector by quaternion.

    Uses optimized formula: v' = v + 2*w*(q_xyz × v) + 2*(q_xyz × (q_xyz × v))
    which is equivalent to q * v * q' but faster.

    Args:
        q (ndarray): (4,) unit quaternion [w, x, y, z]
        v (ndarray): (3,) vector to rotate

    Returns:
        ndarray: (3,) rotated vector
    """
    w, x, y, z = q[0], q[1], q[2], q[3]
    vx, vy, vz = v[0], v[1], v[2]

    # t = 2 * cross(q.xyz, v)
    tx = 2.0 * (y * vz - z * vy)
    ty = 2.0 * (z * vx - x * vz)
    tz = 2.0 * (x * vy - y * vx)

    # result = v + w*t + cross(q.xyz, t)
    return np.array([
        vx + w * tx + (y * tz - z * ty),
        vy + w * ty + (z * tx - x * tz),
        vz + w * tz + (x * ty - y * tx),
    ])


@njit(cache=True)
def quat_rotate_vector_inverse(q, v):
    """Rotate 3D vector by inverse (conjugate) of quaternion.

    Equivalent to q' * v * q, i.e., the reverse rotation.

    Args:
        q (ndarray): (4,) unit quaternion [w, x, y, z]
        v (ndarray): (3,) vector to rotate

    Returns:
        ndarray: (3,) rotated vector
    """
    # For unit quaternion, inverse = conjugate, so negate x,y,z
    w, x, y, z = q[0], -q[1], -q[2], -q[3]
    vx, vy, vz = v[0], v[1], v[2]

    tx = 2.0 * (y * vz - z * vy)
    ty = 2.0 * (z * vx - x * vz)
    tz = 2.0 * (x * vy - y * vx)

    return np.array([
        vx + w * tx + (y * tz - z * ty),
        vy + w * ty + (z * tx - x * tz),
        vz + w * tz + (x * ty - y * tx),
    ])


@njit(cache=True)
def quat_rotate_points(q, points):
    """Rotate array of 3D points by quaternion.

    Args:
        q (ndarray): (4,) unit quaternion [w, x, y, z]
        points (ndarray): (N, 3) array of points to rotate

    Returns:
        ndarray: (N, 3) array of rotated points
    """
    n = points.shape[0]
    result = np.empty((n, 3))

    w, x, y, z = q[0], q[1], q[2], q[3]

    for i in range(n):
        vx, vy, vz = points[i, 0], points[i, 1], points[i, 2]

        tx = 2.0 * (y * vz - z * vy)
        ty = 2.0 * (z * vx - x * vz)
        tz = 2.0 * (x * vy - y * vx)

        result[i, 0] = vx + w * tx + (y * tz - z * ty)
        result[i, 1] = vy + w * ty + (z * tx - x * tz)
        result[i, 2] = vz + w * tz + (x * ty - y * tx)

    return result


@njit(cache=True)
def quat_rotate_points_inverse(q, points):
    """Rotate array of 3D points by inverse of quaternion.

    Args:
        q (ndarray): (4,) unit quaternion [w, x, y, z]
        points (ndarray): (N, 3) array of points to rotate

    Returns:
        ndarray: (N, 3) array of rotated points
    """
    n = points.shape[0]
    result = np.empty((n, 3))

    # Conjugate for inverse rotation
    w, x, y, z = q[0], -q[1], -q[2], -q[3]

    for i in range(n):
        vx, vy, vz = points[i, 0], points[i, 1], points[i, 2]

        tx = 2.0 * (y * vz - z * vy)
        ty = 2.0 * (z * vx - x * vz)
        tz = 2.0 * (x * vy - y * vx)

        result[i, 0] = vx + w * tx + (y * tz - z * ty)
        result[i, 1] = vy + w * ty + (z * tx - x * tz)
        result[i, 2] = vz + w * tz + (x * ty - y * tx)

    return result


# =============================================================================
# Utility Functions
# =============================================================================


@njit(cache=True)
def quat_norm(q):
    """Compute norm (magnitude) of quaternion.

    Args:
        q (ndarray): (4,) quaternion

    Returns:
        float: quaternion norm
    """
    return np.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3] ** 2)


@njit(cache=True)
def quat_normalize(q):
    """Normalize quaternion to unit length.

    Args:
        q (ndarray): (4,) quaternion

    Returns:
        ndarray: (4,) unit quaternion
    """
    norm = quat_norm(q)
    if norm < FP_CUTOFF:
        return quat_identity()
    return q / norm


@njit(cache=True)
def quat_to_axis_angle(q):
    """Extract axis-angle representation from quaternion.

    Args:
        q (ndarray): (4,) unit quaternion [w, x, y, z]

    Returns:
        tuple: (theta, axis) where theta is angle and axis is (3,) unit vector
    """
    # Clamp w to [-1, 1] to handle numerical errors
    w = q[0]
    if w > 1.0:
        w = 1.0
    elif w < -1.0:
        w = -1.0

    theta = 2.0 * np.arccos(w)

    # Compute axis
    s = np.sqrt(1.0 - w * w)
    if s < FP_CUTOFF:
        # Angle is ~0, axis is arbitrary
        axis = np.array([1.0, 0.0, 0.0])
    else:
        axis = np.array([q[1] / s, q[2] / s, q[3] / s])

    return theta, axis


@njit(cache=True)
def quat_to_rotation_matrix(q):
    """Convert quaternion to 3x3 rotation matrix.

    Args:
        q (ndarray): (4,) unit quaternion [w, x, y, z]

    Returns:
        ndarray: (3, 3) rotation matrix
    """
    w, x, y, z = q[0], q[1], q[2], q[3]

    # Pre-compute products
    xx = x * x
    yy = y * y
    zz = z * z
    xy = x * y
    xz = x * z
    yz = y * z
    wx = w * x
    wy = w * y
    wz = w * z

    R = np.empty((3, 3))

    R[0, 0] = 1.0 - 2.0 * (yy + zz)
    R[0, 1] = 2.0 * (xy - wz)
    R[0, 2] = 2.0 * (xz + wy)

    R[1, 0] = 2.0 * (xy + wz)
    R[1, 1] = 1.0 - 2.0 * (xx + zz)
    R[1, 2] = 2.0 * (yz - wx)

    R[2, 0] = 2.0 * (xz - wy)
    R[2, 1] = 2.0 * (yz + wx)
    R[2, 2] = 1.0 - 2.0 * (xx + yy)

    return R


# =============================================================================
# Comparison and Checks
# =============================================================================


@njit(cache=True)
def quat_is_identity(q, tol=1e-10):
    """Check if quaternion is approximately identity.

    Args:
        q (ndarray): (4,) quaternion
        tol (float): tolerance for comparison

    Returns:
        bool: True if quaternion is approximately [1, 0, 0, 0]
    """
    return (
        np.abs(q[0] - 1.0) < tol
        and np.abs(q[1]) < tol
        and np.abs(q[2]) < tol
        and np.abs(q[3]) < tol
    )


@njit(cache=True)
def vectors_parallel(v1, v2, tol=1e-6):
    """Check if two vectors are parallel (same or opposite direction).

    Args:
        v1 (ndarray): (3,) first vector
        v2 (ndarray): (3,) second vector
        tol (float): tolerance for cross product magnitude

    Returns:
        bool: True if vectors are parallel
    """
    cross = np.array([
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    ])
    return np.sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2) < tol


@njit(cache=True)
def vectors_same_direction(v1, v2):
    """Check if two vectors point in the same direction.

    Args:
        v1 (ndarray): (3,) first vector (should be normalized)
        v2 (ndarray): (3,) second vector (should be normalized)

    Returns:
        bool: True if dot product > 0
    """
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] > 0


# =============================================================================
# Safe Math Operations
# =============================================================================


@njit(cache=True)
def safe_arccos(x):
    """Compute arccos with domain clamping to prevent NaN.

    Args:
        x (float): input value

    Returns:
        float: arccos of clamped value
    """
    if x > 1.0:
        return 0.0
    elif x < -1.0:
        return np.pi
    return np.arccos(x)


@njit(cache=True)
def vec3_normalize(v):
    """Normalize a 3D vector.

    Args:
        v (ndarray): (3,) vector

    Returns:
        ndarray: (3,) unit vector (or zero vector if input is zero)
    """
    norm = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    if norm < FP_CUTOFF:
        return np.array([0.0, 0.0, 0.0])
    return v / norm


@njit(cache=True)
def vec3_cross(a, b):
    """Compute cross product of two 3D vectors.

    Args:
        a (ndarray): (3,) first vector
        b (ndarray): (3,) second vector

    Returns:
        ndarray: (3,) cross product a × b
    """
    return np.array([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ])


@njit(cache=True)
def vec3_dot(a, b):
    """Compute dot product of two 3D vectors.

    Args:
        a (ndarray): (3,) first vector
        b (ndarray): (3,) second vector

    Returns:
        float: dot product a · b
    """
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


@njit(cache=True)
def vec3_norm(v):
    """Compute norm (magnitude) of a 3D vector.

    Args:
        v (ndarray): (3,) vector

    Returns:
        float: vector norm
    """
    return np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
