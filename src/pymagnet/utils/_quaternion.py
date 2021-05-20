# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Quaternion module

Implements quaternion multiplication for convenient rotation of vectors in 3D.

Example:
    Rotation of a vector about the x-axis:

        import numpy as np
        import pymagnet as pm
        vector1 = np.array([1,0,0])
        rotate_about_z = pm.magnets.Quaternion.q_angle_from_axis(np.pi/2, (0, 0, 1))
        vector2 = rotate_about_z * vector1

"""
__all__ = ["Quaternion"]

import numpy as _np
from ..utils.global_const import FP_CUTOFF, MAG_TOL


class Quaternion:
    """Quaternion class.
    overloading of multiplication symbol allows easy quaternion multiplications

    """

    def __init__(self, w=1.0, x=0.0, y=0.0, z=0.0):
        """Initialse a pure quaternion (1; 0, 0, 0)

        Args:
            w (ndarray, optional): scalar quaternion. Defaults to 1.0.
            x (ndarray, optional): vector component. Defaults to 0.0.
            y (ndarray, optional): vector component. Defaults to 0.0.
            z (ndarray, optional): vector component. Defaults to 0.0.
        """
        self.w = _np.asarray(w)
        self.x = _np.asarray(x)
        self.y = _np.asarray(y)
        self.z = _np.asarray(z)

    def q_angle_from_axis(theta, vec):
        """Generates a rotation quaternion for an angle `theta` about an axis `vec`
        This is a normailsed, i.e. unit quaternion.

        Args:
            theta (float): angle of rotation
            vec (tuple/array): axis vector

        Example:
            90 degree rotation about the x axis:

                rotation_quaternion = Quaternion.q_angle_from_axis(np.pi/2, (1, 0, 0) )

        Returns:
            Quaternion: rotation quaternion
        """
        vec = Quaternion._normalise_axis(vec)
        w = _np.cos(theta / 2.0)
        vec *= _np.sin(theta / 2.0)
        x = vec[0]
        y = vec[1]
        z = vec[2]
        rotation_quaternion = Quaternion(w, x, y, z)
        return rotation_quaternion

    def get_conjugate(self):
        """Returns quaternion conjugate

        Returns:
            quaternion: quaternion conjugate
        """
        return Quaternion(self.w, -self.x, -self.y, -self.z)

    @staticmethod
    def gen_rotation_quaternion(alpha_rad=0.0, beta_rad=0.0, gamma_rad=0.0):
        """Generates quaternion for rotation around z, y, x axes

        Args:
            alpha_rad (float): angle to z-axis. Defaults to 0.0.
            beta_rad (float): angle to y-axis. Defaults to 0.0.
            gamma_rad (float): angle to x-axis. Defaults to 0.0.

        Returns:
            Quaternion: rotation quaternion
        """

        rotate_about_x = Quaternion()
        rotate_about_y = Quaternion()
        rotate_about_z = Quaternion()

        forward_rotation = Quaternion()

        if _np.fabs(alpha_rad) > MAG_TOL:
            rotate_about_z = Quaternion.q_angle_from_axis(alpha_rad, (0, 0, 1))

        if _np.fabs(beta_rad) > MAG_TOL:
            rotate_about_y = Quaternion.q_angle_from_axis(beta_rad, (0, 1, 0))

        if _np.fabs(gamma_rad) > MAG_TOL:
            rotate_about_x = Quaternion.q_angle_from_axis(gamma_rad, (1, 0, 0))

        forward_rotation = rotate_about_x * rotate_about_z * rotate_about_y

        return forward_rotation

    @staticmethod
    def _prepare_vector(x, y, z):
        """Creates 3xN array where each column corresponds to x, y, z vectors
        all of the same length. If some vectors are shorter, their values are repeated
        to ensure and array of (x,y,z) points for quaternion rotation.

        Args:
            x (ndarray): x coordinates
            y (ndarray): y coordinates
            z (ndarray): z coordinates

        Returns:
            ndarray: 3xN numpy ndarray
        """
        # Check input x,y,z are numpy arrays and if not, convert to numpy arrays
        x = _np.asarray(x)
        y = _np.asarray(y)
        z = _np.asarray(z)

        longest_array_length = _np.max([x.size, y.size, z.size])

        # Ensure the unravelled arrays are all of the same length
        x = Quaternion._check_extend_array(x.ravel(), longest_array_length)
        y = Quaternion._check_extend_array(y.ravel(), longest_array_length)
        z = Quaternion._check_extend_array(z.ravel(), longest_array_length)

        return _np.array([x, y, z])

    @staticmethod
    def _check_extend_array(array, max):
        """Extend a 1D vector to length 'max' if it is shorter than 'max'

        Args:
            array (ndarray): numpy array
            max (int): array length to compare to

        Returns:
            ndarray: numpy array of length 'max'
        """

        if max % array.size != 0:
            raise ValueError("Incorrect gridding of x,y,z")
        # If array is smaller than longest array length 'max', extend it by tiling
        if array.size < max:
            extended_array = _np.tile(array, (max // array.size))

            # # If there is a remainder, 'n' append the first n elements of the array
            # # to ensure the final length is 'max'.
            # if max % array.size != 0:
            #     extended_array = _np.append(
            #         extended_array, array[0 : (max % array.size)]
            #     )
            return extended_array
        else:
            return array

    @staticmethod
    def _normalise_axis(vec):
        """Normalise

        Args:
            v (array): [description]

        Returns:
            [type]: [description]
        """
        vec = _np.asarray(vec)
        if _np.fabs(_np.linalg.norm(vec, axis=0)) < FP_CUTOFF:
            raise ValueError("Vec norm should be non-zero")
        return vec / _np.linalg.norm(vec, axis=0)

    @staticmethod
    def vec_norm(x, y, z):
        """Normalises each x,y,z vector

        Args:
            x (ndarray): x array
            y (ndarray): y array
            z (ndarray): z array

        Returns:
            array: 3xN array of normalised vectors
        """
        vec = Quaternion._prepare_vector(x, y, z)
        return _np.linalg.norm(vec, axis=0)

    def __repr__(self):
        str = f"({self.w}; {self.x}, {self.y}, {self.z})"
        return str

    def __mul__(self, b):
        if isinstance(b, Quaternion):
            return self._multiply_with_quaternion(b)
        elif isinstance(b, (list, tuple, _np.ndarray)):
            if len(b) != 3:
                raise Exception(f"Input vector has invalid length {len(b)}")
            return self._multiply_with_vector(b)
        else:
            raise Exception(f"Multiplication with unknown type {type(b)}")

    def as_tuple(self):
        """Returns quaternion as tuple of arrays

        Returns:
            tuple: w (array), x (array), y (array), z (array)
        """
        return self.w, self.x, self.y, self.z

    def _multiply_with_quaternion(self, q2):
        """Computes the Hamilton product of two quaternions

        Args:
            q2 (quaternion): quaternion

        Returns:
            quaternion: Hamilton product
        """
        w1, x1, y1, z1 = self.as_tuple()
        w2, x2, y2, z2 = q2.as_tuple()
        w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
        x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
        y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
        z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
        result = Quaternion(w, x, y, z)
        return result

    def _multiply_with_vector(self, v):
        """Converts vector ``v`` to quaternion ``q2`` and then performs rotation by
        multiplying `q * q2 * q'`


        Args:
            v (vector/array): vector to be rotated

        Returns:
            tuple: multiplied vector x', y', z'
        """
        q2 = Quaternion(_np.zeros_like(v[0]), v[0], v[1], v[2])
        result = self * q2 * self.get_conjugate()
        return result.x, result.y, result.z

    def get_axisangle(self):
        theta = _np.arccos(self.w) * 2.0
        vec = _np.hstack([self.x, self.y, self.z])
        return theta, self._normalise_axis(vec)

    @staticmethod
    def euler_to_quaternion(alpha, beta, gamma):
        """Converts Euler angles to quaternion

        Args:
            alpha (float): angle to z-axis
            beta (float): angle to y-axis
            gamma (float): angle to x-axis

        Returns:
            Quaternion: Euler angles as a Quaternion
        """

        qw = _np.cos(alpha / 2) * _np.cos(beta / 2) * _np.cos(gamma / 2) + _np.sin(
            alpha / 2
        ) * _np.sin(beta / 2) * _np.sin(gamma / 2)
        qx = _np.sin(alpha / 2) * _np.cos(beta / 2) * _np.cos(gamma / 2) - _np.cos(
            alpha / 2
        ) * _np.sin(beta / 2) * _np.sin(gamma / 2)
        qy = _np.cos(alpha / 2) * _np.sin(beta / 2) * _np.cos(gamma / 2) + _np.sin(
            alpha / 2
        ) * _np.cos(beta / 2) * _np.sin(gamma / 2)
        qz = _np.cos(alpha / 2) * _np.cos(beta / 2) * _np.sin(gamma / 2) - _np.sin(
            alpha / 2
        ) * _np.sin(beta / 2) * _np.cos(gamma / 2)

        return Quaternion(qw, qx, qy, qz)