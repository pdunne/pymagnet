# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for Quaternion class and rotation operations.

Critical tests for:
- Quaternion initialization and properties
- Hamilton product (quaternion multiplication)
- Vector rotation
- Euler angle conversion
- Axis-angle generation
"""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.utils._quaternion import Quaternion, q_angle_from_axis
from pymagnet.utils.global_const import PI


class TestQuaternionInit:
    """Tests for Quaternion initialization."""

    def test_default_identity(self):
        """Default quaternion should be identity (1; 0,0,0)."""
        q = Quaternion()
        assert q.w == 1.0
        assert q.x == 0.0
        assert q.y == 0.0
        assert q.z == 0.0

    def test_custom_components(self):
        """Quaternion with custom w,x,y,z components."""
        q = Quaternion(0.5, 0.5, 0.5, 0.5)
        assert q.w == 0.5
        assert q.x == 0.5
        assert q.y == 0.5
        assert q.z == 0.5

    def test_array_components(self):
        """Quaternion with numpy array components."""
        w = np.array([1.0, 0.0])
        x = np.array([0.0, 1.0])
        y = np.array([0.0, 0.0])
        z = np.array([0.0, 0.0])
        q = Quaternion(w, x, y, z)
        npt.assert_array_equal(q.w, w)
        npt.assert_array_equal(q.x, x)

    def test_repr(self):
        """String representation of quaternion."""
        q = Quaternion(1.0, 0.0, 0.0, 0.0)
        repr_str = repr(q)
        assert "1.0" in repr_str


class TestQuaternionConjugate:
    """Tests for quaternion conjugate operation."""

    def test_conjugate_signs(self):
        """Conjugate negates x,y,z but preserves w."""
        q = Quaternion(1.0, 2.0, 3.0, 4.0)
        qc = q.get_conjugate()
        assert qc.w == 1.0
        assert qc.x == -2.0
        assert qc.y == -3.0
        assert qc.z == -4.0

    def test_double_conjugate_identity(self):
        """Double conjugate returns original."""
        q = Quaternion(1.0, 2.0, 3.0, 4.0)
        qcc = q.get_conjugate().get_conjugate()
        assert qcc.w == q.w
        assert qcc.x == q.x
        assert qcc.y == q.y
        assert qcc.z == q.z

    def test_identity_conjugate(self, identity_quaternion):
        """Conjugate of identity is identity."""
        qc = identity_quaternion.get_conjugate()
        assert qc.w == 1.0
        assert qc.x == 0.0
        assert qc.y == 0.0
        assert qc.z == 0.0


class TestQuaternionMultiplication:
    """Tests for Hamilton product."""

    def test_identity_multiplication_left(self, identity_quaternion):
        """identity * q = q."""
        q = Quaternion(0.5, 0.5, 0.5, 0.5)
        result = identity_quaternion * q
        npt.assert_allclose(result.w, q.w)
        npt.assert_allclose(result.x, q.x)
        npt.assert_allclose(result.y, q.y)
        npt.assert_allclose(result.z, q.z)

    def test_identity_multiplication_right(self, identity_quaternion):
        """q * identity = q."""
        q = Quaternion(0.5, 0.5, 0.5, 0.5)
        result = q * identity_quaternion
        npt.assert_allclose(result.w, q.w)
        npt.assert_allclose(result.x, q.x)
        npt.assert_allclose(result.y, q.y)
        npt.assert_allclose(result.z, q.z)

    def test_multiplication_associativity(self):
        """(q1 * q2) * q3 = q1 * (q2 * q3)."""
        q1 = Quaternion(0.5, 0.5, 0.5, 0.5)
        q2 = Quaternion(0.7, 0.1, 0.3, 0.6)
        q3 = Quaternion(0.2, 0.8, 0.4, 0.1)

        left = (q1 * q2) * q3
        right = q1 * (q2 * q3)

        npt.assert_allclose(left.w, right.w, rtol=1e-10)
        npt.assert_allclose(left.x, right.x, rtol=1e-10)
        npt.assert_allclose(left.y, right.y, rtol=1e-10)
        npt.assert_allclose(left.z, right.z, rtol=1e-10)

    def test_multiplication_non_commutativity(self):
        """q1 * q2 != q2 * q1 in general."""
        q1 = q_angle_from_axis(PI / 2, (1, 0, 0))
        q2 = q_angle_from_axis(PI / 2, (0, 1, 0))

        r1 = q1 * q2
        r2 = q2 * q1

        # They should NOT be equal
        assert not np.allclose([r1.w, r1.x, r1.y, r1.z], [r2.w, r2.x, r2.y, r2.z])

    def test_q_times_conjugate_gives_scalar(self):
        """q * q_conjugate gives |q|^2 as scalar quaternion."""
        q = Quaternion(0.5, 0.5, 0.5, 0.5)
        result = q * q.get_conjugate()
        # Result should be scalar (w != 0, x=y=z=0)
        npt.assert_allclose(result.x, 0.0, atol=1e-10)
        npt.assert_allclose(result.y, 0.0, atol=1e-10)
        npt.assert_allclose(result.z, 0.0, atol=1e-10)
        # w should be |q|^2
        norm_sq = 0.5**2 + 0.5**2 + 0.5**2 + 0.5**2
        npt.assert_allclose(result.w, norm_sq, rtol=1e-10)


class TestVectorRotation:
    """Tests for rotating vectors with quaternions."""

    def test_rotation_90_about_z(self, rotation_90_z):
        """90 degree rotation about z: (1,0,0) -> (0,1,0)."""
        v = np.array([1.0, 0.0, 0.0])
        result = rotation_90_z * v
        npt.assert_allclose(result[0], 0.0, atol=1e-10)
        npt.assert_allclose(result[1], 1.0, atol=1e-10)
        npt.assert_allclose(result[2], 0.0, atol=1e-10)

    def test_rotation_90_about_x(self, rotation_90_x):
        """90 degree rotation about x: (0,1,0) -> (0,0,1)."""
        v = np.array([0.0, 1.0, 0.0])
        result = rotation_90_x * v
        npt.assert_allclose(result[0], 0.0, atol=1e-10)
        npt.assert_allclose(result[1], 0.0, atol=1e-10)
        npt.assert_allclose(result[2], 1.0, atol=1e-10)

    def test_rotation_90_about_y(self, rotation_90_y):
        """90 degree rotation about y: (0,0,1) -> (1,0,0)."""
        v = np.array([0.0, 0.0, 1.0])
        result = rotation_90_y * v
        npt.assert_allclose(result[0], 1.0, atol=1e-10)
        npt.assert_allclose(result[1], 0.0, atol=1e-10)
        npt.assert_allclose(result[2], 0.0, atol=1e-10)

    def test_rotation_preserves_magnitude(self, rotation_90_z):
        """Rotation should preserve vector magnitude."""
        v = np.array([3.0, 4.0, 5.0])
        original_mag = np.linalg.norm(v)
        result = rotation_90_z * v
        result_mag = np.sqrt(result[0] ** 2 + result[1] ** 2 + result[2] ** 2)
        npt.assert_allclose(result_mag, original_mag, rtol=1e-10)

    @pytest.mark.parametrize("angle", [0, PI / 4, PI / 2, PI, 3 * PI / 2, 2 * PI])
    def test_rotation_angles_preserve_magnitude(self, angle):
        """Test various rotation angles preserve magnitude."""
        q = q_angle_from_axis(angle, (0, 0, 1))
        v = np.array([1.0, 0.0, 0.0])
        result = q * v
        result_mag = np.sqrt(result[0] ** 2 + result[1] ** 2 + result[2] ** 2)
        npt.assert_allclose(result_mag, 1.0, rtol=1e-10)

    def test_identity_rotation_unchanged(self, identity_quaternion):
        """Identity quaternion should not change vector."""
        v = np.array([1.0, 2.0, 3.0])
        result = identity_quaternion * v
        npt.assert_allclose(result[0], 1.0, atol=1e-10)
        npt.assert_allclose(result[1], 2.0, atol=1e-10)
        npt.assert_allclose(result[2], 3.0, atol=1e-10)

    def test_rotation_180_about_z(self, rotation_180_z):
        """180 degree rotation about z: (1,0,0) -> (-1,0,0)."""
        v = np.array([1.0, 0.0, 0.0])
        result = rotation_180_z * v
        npt.assert_allclose(result[0], -1.0, atol=1e-10)
        npt.assert_allclose(result[1], 0.0, atol=1e-10)
        npt.assert_allclose(result[2], 0.0, atol=1e-10)

    def test_rotation_with_list_input(self, rotation_90_z):
        """Rotation with list input instead of array."""
        v = [1.0, 0.0, 0.0]
        result = rotation_90_z * v
        npt.assert_allclose(result[0], 0.0, atol=1e-10)
        npt.assert_allclose(result[1], 1.0, atol=1e-10)

    def test_rotation_with_tuple_input(self, rotation_90_z):
        """Rotation with tuple input."""
        v = (1.0, 0.0, 0.0)
        result = rotation_90_z * v
        npt.assert_allclose(result[0], 0.0, atol=1e-10)
        npt.assert_allclose(result[1], 1.0, atol=1e-10)


class TestAxisAngle:
    """Tests for q_angle_from_axis function."""

    def test_zero_angle_identity(self):
        """Zero angle should give identity quaternion."""
        q = q_angle_from_axis(0.0, (0, 0, 1))
        npt.assert_allclose(q.w, 1.0, atol=1e-10)
        npt.assert_allclose(q.x, 0.0, atol=1e-10)
        npt.assert_allclose(q.y, 0.0, atol=1e-10)
        npt.assert_allclose(q.z, 0.0, atol=1e-10)

    def test_axis_normalization(self):
        """Non-unit axis vectors should be normalized."""
        # Use non-unit axis (2,0,0) instead of (1,0,0)
        q = q_angle_from_axis(PI / 2, (2, 0, 0))
        # Should behave same as unit axis
        v = np.array([0.0, 1.0, 0.0])
        result = q * v
        npt.assert_allclose(result[0], 0.0, atol=1e-10)
        npt.assert_allclose(result[1], 0.0, atol=1e-10)
        npt.assert_allclose(result[2], 1.0, atol=1e-10)

    def test_zero_norm_axis_raises(self):
        """Zero-norm axis should raise ValueError."""
        with pytest.raises(ValueError):
            q_angle_from_axis(PI / 2, (0, 0, 0))

    @pytest.mark.parametrize(
        "axis", [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]
    )
    def test_principal_axes(self, axis):
        """Test rotation about principal axes."""
        q = q_angle_from_axis(PI / 2, axis)
        # Quaternion should be normalized
        norm = np.sqrt(q.w**2 + q.x**2 + q.y**2 + q.z**2)
        npt.assert_allclose(norm, 1.0, rtol=1e-10)

    def test_2pi_rotation_identity(self):
        """2*PI rotation should return to original."""
        q = q_angle_from_axis(2 * PI, (0, 0, 1))
        v = np.array([1.0, 0.0, 0.0])
        result = q * v
        npt.assert_allclose(result[0], 1.0, atol=1e-10)
        npt.assert_allclose(result[1], 0.0, atol=1e-10)
        npt.assert_allclose(result[2], 0.0, atol=1e-10)


class TestGenRotationQuaternion:
    """Tests for gen_rotation_quaternion static method."""

    def test_zero_angles_identity(self):
        """All zero Euler angles should give identity."""
        q = Quaternion.gen_rotation_quaternion(0.0, 0.0, 0.0)
        npt.assert_allclose(q.w, 1.0, atol=1e-10)
        npt.assert_allclose(q.x, 0.0, atol=1e-10)
        npt.assert_allclose(q.y, 0.0, atol=1e-10)
        npt.assert_allclose(q.z, 0.0, atol=1e-10)

    def test_single_axis_alpha(self):
        """Rotation about z-axis only (alpha)."""
        q = Quaternion.gen_rotation_quaternion(alpha_rad=PI / 2)
        v = np.array([1.0, 0.0, 0.0])
        result = q * v
        npt.assert_allclose(result[0], 0.0, atol=1e-10)
        npt.assert_allclose(result[1], 1.0, atol=1e-10)
        npt.assert_allclose(result[2], 0.0, atol=1e-10)

    def test_single_axis_beta(self):
        """Rotation about y-axis only (beta)."""
        q = Quaternion.gen_rotation_quaternion(beta_rad=PI / 2)
        v = np.array([0.0, 0.0, 1.0])
        result = q * v
        npt.assert_allclose(result[0], 1.0, atol=1e-10)
        npt.assert_allclose(result[1], 0.0, atol=1e-10)
        npt.assert_allclose(result[2], 0.0, atol=1e-10)

    def test_single_axis_gamma(self):
        """Rotation about x-axis only (gamma)."""
        q = Quaternion.gen_rotation_quaternion(gamma_rad=PI / 2)
        v = np.array([0.0, 1.0, 0.0])
        result = q * v
        npt.assert_allclose(result[0], 0.0, atol=1e-10)
        npt.assert_allclose(result[1], 0.0, atol=1e-10)
        npt.assert_allclose(result[2], 1.0, atol=1e-10)

    def test_combined_rotation(self):
        """Combined rotations about multiple axes."""
        q = Quaternion.gen_rotation_quaternion(
            alpha_rad=PI / 4, beta_rad=PI / 6, gamma_rad=PI / 3
        )
        # Just verify it creates a normalized quaternion
        norm = np.sqrt(q.w**2 + q.x**2 + q.y**2 + q.z**2)
        npt.assert_allclose(norm, 1.0, rtol=1e-10)


class TestEulerToQuaternion:
    """Tests for euler_to_quaternion static method."""

    def test_zero_angles(self):
        """Zero Euler angles give identity."""
        q = Quaternion.euler_to_quaternion(0, 0, 0)
        npt.assert_allclose(q.w, 1.0, atol=1e-10)
        npt.assert_allclose(q.x, 0.0, atol=1e-10)
        npt.assert_allclose(q.y, 0.0, atol=1e-10)
        npt.assert_allclose(q.z, 0.0, atol=1e-10)

    def test_normalized_result(self):
        """Result should be a unit quaternion."""
        q = Quaternion.euler_to_quaternion(PI / 4, PI / 3, PI / 6)
        norm = np.sqrt(q.w**2 + q.x**2 + q.y**2 + q.z**2)
        npt.assert_allclose(norm, 1.0, rtol=1e-10)


class TestQuaternionMultiplyErrors:
    """Test error handling in quaternion multiplication."""

    def test_invalid_vector_length(self):
        """Vector with wrong length should raise exception."""
        q = Quaternion()
        with pytest.raises(Exception):
            q * [1, 2]  # Too short

    def test_invalid_type(self):
        """Invalid type should raise exception."""
        q = Quaternion()
        with pytest.raises(Exception):
            q * "invalid"

    def test_invalid_type_number(self):
        """Scalar number should raise exception."""
        q = Quaternion()
        with pytest.raises(Exception):
            q * 5.0


class TestPrepareVector:
    """Tests for _prepare_vector static method."""

    def test_scalar_inputs(self):
        """Scalar inputs converted to arrays."""
        result = Quaternion._prepare_vector(1.0, 2.0, 3.0)
        assert result.shape == (3, 1)
        npt.assert_allclose(result[:, 0], [1.0, 2.0, 3.0])

    def test_array_inputs_same_length(self):
        """Arrays of same length."""
        x = np.array([1.0, 2.0])
        y = np.array([3.0, 4.0])
        z = np.array([5.0, 6.0])
        result = Quaternion._prepare_vector(x, y, z)
        assert result.shape == (3, 2)

    def test_array_extension(self):
        """Shorter arrays should be extended."""
        x = np.array([1.0, 2.0])
        y = np.array([3.0])  # Will be tiled
        z = np.array([5.0])  # Will be tiled
        result = Quaternion._prepare_vector(x, y, z)
        assert result.shape == (3, 2)
        # y should be [3.0, 3.0]
        npt.assert_allclose(result[1, :], [3.0, 3.0])


class TestVecNorm:
    """Tests for vec_norm static method."""

    def test_unit_vector(self):
        """Unit vector has norm 1."""
        result = Quaternion.vec_norm(1.0, 0.0, 0.0)
        npt.assert_allclose(result, 1.0)

    def test_345_vector(self):
        """3-4-0 vector has norm 5 (no z component)."""
        result = Quaternion.vec_norm(3.0, 4.0, 0.0)
        npt.assert_allclose(result, 5.0)

    def test_array_norms(self):
        """Multiple vectors return multiple norms."""
        x = np.array([1.0, 3.0])
        y = np.array([0.0, 4.0])
        z = np.array([0.0, 0.0])
        result = Quaternion.vec_norm(x, y, z)
        npt.assert_allclose(result, [1.0, 5.0])


class TestGetAxisAngle:
    """Tests for get_axisangle method."""

    def test_identity_angle(self):
        """Near-identity quaternion has near-zero angle.

        Note: Pure identity quaternion (1,0,0,0) raises ValueError when
        normalizing (0,0,0) axis. This is expected behavior.
        """
        # Test with a very small rotation instead
        q = q_angle_from_axis(0.001, (0, 0, 1))
        theta, axis = q.get_axisangle()
        npt.assert_allclose(theta, 0.001, atol=1e-6)

    def test_90_degree_rotation(self, rotation_90_z):
        """90 degree rotation about z."""
        theta, axis = rotation_90_z.get_axisangle()
        npt.assert_allclose(theta, PI / 2, rtol=1e-10)
        # Axis should be along z
        npt.assert_allclose(axis[2], 1.0, atol=1e-10)
