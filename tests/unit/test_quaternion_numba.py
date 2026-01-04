# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for numba-compatible quaternion operations."""

import numpy as np
import pytest

from pymagnet.utils._quaternion import Quaternion, q_angle_from_axis
from pymagnet.utils._quaternion_numba import (
    quat_conjugate,
    quat_from_axis_angle,
    quat_from_euler,
    quat_identity,
    quat_is_identity,
    quat_multiply,
    quat_multiply3,
    quat_norm,
    quat_normalize,
    quat_rotate_points,
    quat_rotate_points_inverse,
    quat_rotate_vector,
    quat_rotate_vector_inverse,
    quat_to_axis_angle,
    quat_to_rotation_matrix,
    safe_arccos,
    vec3_cross,
    vec3_dot,
    vec3_norm,
    vec3_normalize,
    vectors_parallel,
    vectors_same_direction,
)


class TestQuaternionBasics:
    """Test basic quaternion operations."""

    def test_identity(self):
        """Test identity quaternion creation."""
        q = quat_identity()
        assert np.allclose(q, [1, 0, 0, 0])

    def test_conjugate(self):
        """Test quaternion conjugate."""
        q = np.array([0.5, 0.5, 0.5, 0.5])
        qc = quat_conjugate(q)
        assert np.allclose(qc, [0.5, -0.5, -0.5, -0.5])

    def test_norm(self):
        """Test quaternion norm calculation."""
        q = np.array([1.0, 2.0, 3.0, 4.0])
        expected = np.sqrt(1 + 4 + 9 + 16)
        assert np.isclose(quat_norm(q), expected)

    def test_normalize(self):
        """Test quaternion normalization."""
        q = np.array([1.0, 2.0, 3.0, 4.0])
        qn = quat_normalize(q)
        assert np.isclose(quat_norm(qn), 1.0)

    def test_is_identity(self):
        """Test identity check."""
        assert quat_is_identity(quat_identity())
        assert not quat_is_identity(np.array([0.0, 1.0, 0.0, 0.0]))


class TestQuaternionCreation:
    """Test quaternion creation functions."""

    def test_from_axis_angle_x(self):
        """Test axis-angle creation for x-axis rotation."""
        angle = np.pi / 2
        axis = np.array([1.0, 0.0, 0.0])
        q = quat_from_axis_angle(angle, axis)

        # Should be [cos(π/4), sin(π/4), 0, 0]
        expected = np.array([np.cos(np.pi / 4), np.sin(np.pi / 4), 0, 0])
        assert np.allclose(q, expected)

    def test_from_axis_angle_y(self):
        """Test axis-angle creation for y-axis rotation."""
        angle = np.pi / 2
        axis = np.array([0.0, 1.0, 0.0])
        q = quat_from_axis_angle(angle, axis)

        expected = np.array([np.cos(np.pi / 4), 0, np.sin(np.pi / 4), 0])
        assert np.allclose(q, expected)

    def test_from_axis_angle_z(self):
        """Test axis-angle creation for z-axis rotation."""
        angle = np.pi / 2
        axis = np.array([0.0, 0.0, 1.0])
        q = quat_from_axis_angle(angle, axis)

        expected = np.array([np.cos(np.pi / 4), 0, 0, np.sin(np.pi / 4)])
        assert np.allclose(q, expected)

    def test_from_axis_angle_zero(self):
        """Test that zero angle gives identity."""
        q = quat_from_axis_angle(0.0, np.array([1.0, 0.0, 0.0]))
        assert np.allclose(q, quat_identity())

    def test_from_axis_angle_normalizes(self):
        """Test that axis is normalized."""
        axis = np.array([3.0, 0.0, 0.0])  # Not unit length
        q = quat_from_axis_angle(np.pi / 2, axis)

        # Result should be same as normalized axis
        q_expected = quat_from_axis_angle(np.pi / 2, np.array([1.0, 0.0, 0.0]))
        assert np.allclose(q, q_expected)

    def test_to_axis_angle_roundtrip(self):
        """Test axis-angle roundtrip conversion."""
        angle = np.pi / 3
        axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)

        q = quat_from_axis_angle(angle, axis)
        theta, axis_out = quat_to_axis_angle(q)

        assert np.isclose(theta, angle)
        assert np.allclose(axis_out, axis)


class TestQuaternionMultiplication:
    """Test quaternion multiplication."""

    def test_multiply_identity_left(self):
        """Test multiplication with identity on left."""
        q = np.array([0.5, 0.5, 0.5, 0.5])
        result = quat_multiply(quat_identity(), q)
        assert np.allclose(result, q)

    def test_multiply_identity_right(self):
        """Test multiplication with identity on right."""
        q = np.array([0.5, 0.5, 0.5, 0.5])
        result = quat_multiply(q, quat_identity())
        assert np.allclose(result, q)

    def test_multiply_inverse(self):
        """Test that q * q_conjugate = identity (for unit quaternion)."""
        q = quat_from_axis_angle(np.pi / 4, np.array([1.0, 0.0, 0.0]))
        result = quat_multiply(q, quat_conjugate(q))
        assert np.allclose(result, quat_identity())

    def test_multiply3(self):
        """Test triple multiplication."""
        q1 = quat_from_axis_angle(np.pi / 6, np.array([1.0, 0.0, 0.0]))
        q2 = quat_from_axis_angle(np.pi / 4, np.array([0.0, 1.0, 0.0]))
        q3 = quat_from_axis_angle(np.pi / 3, np.array([0.0, 0.0, 1.0]))

        result = quat_multiply3(q1, q2, q3)
        expected = quat_multiply(quat_multiply(q1, q2), q3)

        assert np.allclose(result, expected)

    def test_multiply_matches_class(self):
        """Test that numba multiply matches Quaternion class."""
        # Create random quaternions
        q1_nb = quat_from_axis_angle(np.pi / 3, np.array([1.0, 2.0, 3.0]))
        q2_nb = quat_from_axis_angle(np.pi / 4, np.array([0.0, 1.0, 1.0]))

        # Create equivalent Quaternion objects
        q1_py = q_angle_from_axis(np.pi / 3, (1.0, 2.0, 3.0))
        q2_py = q_angle_from_axis(np.pi / 4, (0.0, 1.0, 1.0))

        # Multiply using both methods
        result_nb = quat_multiply(q1_nb, q2_nb)
        result_py = q1_py * q2_py

        # Compare
        assert np.allclose(result_nb[0], result_py.w)
        assert np.allclose(result_nb[1], result_py.x)
        assert np.allclose(result_nb[2], result_py.y)
        assert np.allclose(result_nb[3], result_py.z)


class TestVectorRotation:
    """Test vector rotation by quaternion."""

    def test_rotate_x_by_z_90(self):
        """Test rotating x-axis vector 90° about z-axis."""
        q = quat_from_axis_angle(np.pi / 2, np.array([0.0, 0.0, 1.0]))
        v = np.array([1.0, 0.0, 0.0])

        result = quat_rotate_vector(q, v)

        # x -> y
        assert np.allclose(result, [0, 1, 0], atol=1e-10)

    def test_rotate_y_by_x_90(self):
        """Test rotating y-axis vector 90° about x-axis."""
        q = quat_from_axis_angle(np.pi / 2, np.array([1.0, 0.0, 0.0]))
        v = np.array([0.0, 1.0, 0.0])

        result = quat_rotate_vector(q, v)

        # y -> z
        assert np.allclose(result, [0, 0, 1], atol=1e-10)

    def test_rotate_z_by_y_90(self):
        """Test rotating z-axis vector 90° about y-axis."""
        q = quat_from_axis_angle(np.pi / 2, np.array([0.0, 1.0, 0.0]))
        v = np.array([0.0, 0.0, 1.0])

        result = quat_rotate_vector(q, v)

        # z -> x
        assert np.allclose(result, [1, 0, 0], atol=1e-10)

    def test_rotate_identity(self):
        """Test that identity quaternion doesn't change vector."""
        v = np.array([1.0, 2.0, 3.0])
        result = quat_rotate_vector(quat_identity(), v)
        assert np.allclose(result, v)

    def test_rotate_inverse(self):
        """Test that inverse rotation undoes rotation."""
        q = quat_from_axis_angle(np.pi / 3, np.array([1.0, 1.0, 0.0]))
        v = np.array([1.0, 2.0, 3.0])

        v_rot = quat_rotate_vector(q, v)
        v_back = quat_rotate_vector_inverse(q, v_rot)

        assert np.allclose(v_back, v)

    def test_rotate_points(self):
        """Test rotating multiple points."""
        q = quat_from_axis_angle(np.pi / 2, np.array([0.0, 0.0, 1.0]))
        points = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ])

        result = quat_rotate_points(q, points)

        expected = np.array([
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-1.0, 1.0, 0.0],
        ])
        assert np.allclose(result, expected, atol=1e-10)

    def test_rotate_points_inverse(self):
        """Test inverse rotation of multiple points."""
        q = quat_from_axis_angle(np.pi / 4, np.array([1.0, 1.0, 1.0]))
        points = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])

        points_rot = quat_rotate_points(q, points)
        points_back = quat_rotate_points_inverse(q, points_rot)

        assert np.allclose(points_back, points)

    def test_rotation_matches_class(self):
        """Test that numba rotation matches Quaternion class."""
        q_nb = quat_from_axis_angle(np.pi / 3, np.array([1.0, 2.0, 3.0]))
        q_py = q_angle_from_axis(np.pi / 3, (1.0, 2.0, 3.0))

        v = np.array([1.0, 2.0, 3.0])

        result_nb = quat_rotate_vector(q_nb, v)
        result_py = q_py * v

        assert np.allclose(result_nb, result_py)


class TestRotationMatrix:
    """Test rotation matrix conversion."""

    def test_identity_matrix(self):
        """Test that identity quaternion gives identity matrix."""
        R = quat_to_rotation_matrix(quat_identity())
        assert np.allclose(R, np.eye(3))

    def test_rotation_matrix_orthogonal(self):
        """Test that rotation matrix is orthogonal."""
        q = quat_from_axis_angle(np.pi / 5, np.array([1.0, 2.0, 3.0]))
        R = quat_to_rotation_matrix(q)

        # R * R^T should be identity
        assert np.allclose(R @ R.T, np.eye(3))

        # det(R) should be 1
        assert np.isclose(np.linalg.det(R), 1.0)

    def test_rotation_matrix_matches_vector_rotation(self):
        """Test that matrix rotation matches quaternion rotation."""
        q = quat_from_axis_angle(np.pi / 4, np.array([0.0, 1.0, 0.0]))
        R = quat_to_rotation_matrix(q)

        v = np.array([1.0, 2.0, 3.0])

        result_quat = quat_rotate_vector(q, v)
        result_matrix = R @ v

        assert np.allclose(result_quat, result_matrix)


class TestVectorUtilities:
    """Test vector utility functions."""

    def test_vec3_dot(self):
        """Test dot product."""
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([4.0, 5.0, 6.0])
        assert np.isclose(vec3_dot(a, b), 32.0)

    def test_vec3_cross(self):
        """Test cross product."""
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 1.0, 0.0])
        result = vec3_cross(a, b)
        assert np.allclose(result, [0, 0, 1])

    def test_vec3_norm(self):
        """Test vector norm."""
        v = np.array([3.0, 4.0, 0.0])
        assert np.isclose(vec3_norm(v), 5.0)

    def test_vec3_normalize(self):
        """Test vector normalization."""
        v = np.array([3.0, 4.0, 0.0])
        vn = vec3_normalize(v)
        assert np.isclose(vec3_norm(vn), 1.0)
        assert np.allclose(vn, [0.6, 0.8, 0])

    def test_vectors_parallel(self):
        """Test parallel vector detection."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([2.0, 0.0, 0.0])
        v3 = np.array([-1.0, 0.0, 0.0])
        v4 = np.array([0.0, 1.0, 0.0])

        assert vectors_parallel(v1, v2)  # Same direction
        assert vectors_parallel(v1, v3)  # Opposite direction
        assert not vectors_parallel(v1, v4)  # Perpendicular

    def test_vectors_same_direction(self):
        """Test same direction detection."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([2.0, 0.0, 0.0])
        v3 = np.array([-1.0, 0.0, 0.0])

        assert vectors_same_direction(v1, v2)
        assert not vectors_same_direction(v1, v3)

    def test_safe_arccos(self):
        """Test safe arccos with out-of-range values."""
        assert np.isclose(safe_arccos(0.5), np.arccos(0.5))
        assert np.isclose(safe_arccos(1.0), 0.0)
        assert np.isclose(safe_arccos(-1.0), np.pi)
        assert np.isclose(safe_arccos(1.1), 0.0)  # Clamped to 1
        assert np.isclose(safe_arccos(-1.1), np.pi)  # Clamped to -1


class TestTrigonometryNumbaFunctions:
    """Test numba-compatible trigonometry functions."""

    def test_rotate_triangle_njit_basic(self):
        """Test basic triangle rotation."""
        from pymagnet.utils._trigonometry3D import _rotate_triangle_njit

        # Simple triangle in xy plane
        triangle = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.866, 0.0],
        ])

        q, rotated, offset, ra1, ra2 = _rotate_triangle_njit(triangle)

        # Check quaternion is unit length
        assert np.isclose(quat_norm(q), 1.0)

        # Check RA triangles have valid dimensions
        assert ra1[0] >= 0 and ra1[1] >= 0
        assert ra2[0] >= 0 and ra2[1] >= 0

    def test_rotate_triangle_njit_matches_original(self):
        """Test that numba version matches original."""
        from pymagnet.utils._trigonometry3D import _rotate_triangle, _rotate_triangle_njit

        # Test multiple triangles
        triangles = [
            np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]),
            np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 10.]]),
            np.array([[0., 0., 0.], [1., 0., 0.], [0.5, 0., 0.866]]),
        ]

        for triangle in triangles:
            # Original function
            q_orig, rot_orig, off_orig, ra1_orig, ra2_orig = _rotate_triangle(
                triangle, 1.0
            )

            # Numba function
            q_njit, rot_njit, off_njit, ra1_njit, ra2_njit = _rotate_triangle_njit(
                triangle
            )

            # Compare RA triangles (should be identical)
            assert np.allclose(ra1_orig, ra1_njit, rtol=1e-5)
            assert np.allclose(ra2_orig, ra2_njit, rtol=1e-5)

            # Compare offsets (may differ in sign due to different vertex selection)
            # but the offset should be on the rotated triangle
            assert np.allclose(np.abs(off_orig), np.abs(off_njit), rtol=1e-5)

    def test_norm_plane_njit(self):
        """Test numba plane normal calculation."""
        from pymagnet.utils._trigonometry3D import norm_plane, norm_plane_njit

        triangle = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])

        norm_orig = norm_plane(triangle)
        norm_njit = norm_plane_njit(triangle)

        assert np.allclose(norm_orig, norm_njit)

    def test_rotate_vector_by_quat_njit(self):
        """Test coordinate array rotation."""
        from pymagnet.utils._trigonometry3D import rotate_vector_by_quat_njit

        q = quat_from_axis_angle(np.pi / 2, np.array([0.0, 0.0, 1.0]))

        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.0, 0.0, 0.0])
        z = np.array([0.0, 0.0, 0.0])

        x_rot, y_rot, z_rot = rotate_vector_by_quat_njit(q, x, y, z)

        # Each x value should become y
        assert np.allclose(x_rot, [0, 0, 0], atol=1e-10)
        assert np.allclose(y_rot, [1, 2, 3], atol=1e-10)
        assert np.allclose(z_rot, [0, 0, 0], atol=1e-10)
