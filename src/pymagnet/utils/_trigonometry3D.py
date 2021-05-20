# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Contains functions needed to rotate and translate a triangle to lie in the xz plane
and to divide it into two right angled triangles
"""
from ._quaternion import Quaternion
from .global_const import PI, ALIGN_CUTOFF

import numpy as _np
from numba import jit

__all__ = ["signed_area", "norm_plane", "rotate_points", "altitude"]


@jit
def signed_area(triangle):
    """Calculates signed area of a triangle. Area area < 0 for clockwise ordering.
    Assumes the triangle is in the xz plane (i.e. with the normal parallel to y).

    Args:
        triangle (ndarray): 3x3 array of vertices

    Returns:
        float: signed area
    """

    j = 1
    NP = 3
    area = 0.0

    for i in range(NP):
        j = j % NP

        area += (triangle[j][0] - triangle[i][0]) * (triangle[j][2] + triangle[i][2])
        j += 1

    # check winding order of polygon, area < 0 for clockwise ordering of points
    area /= 2.0

    return area


def norm_plane(vec):
    """Calculates the normal to a triangular plane

    Args:
        vec (ndarray/list/tuple): (N,1) array

    Returns:
        ndarray: normal vector (N,)
    """
    norm = _np.cross(vec[1] - vec[0], vec[2] - vec[0])
    norm = norm / _np.linalg.norm(norm)
    return norm


def rotate_points(points, rotation_quaternion):
    """Rotates a set of points

    Args:
        points ([type]): [description]
        rotation_quaternion ([type]): [description]

    Returns:
        [type]: [description]
    """
    x_rot, y_rot, z_rot = rotation_quaternion * points.T
    rotate_points = _np.vstack([x_rot, y_rot, z_rot]).T

    return rotate_points


def altitude(a, b, c):
    """Gets altitude to side `a` of a triangle.

    Args:
        a (float): longest side
        b (float): triangle side
        c (float): triangle side

    Returns:
        float: altitude to side `a`
    """
    s = (a + b + c) / 2
    return 2 * _np.sqrt(s * (s - a) * (s - b) * (s - c)) / a


def _largest_side_RA(triangle):
    """Determines largest side of triangle

    Args:
        triangle (ndarray): triangle vertices

    Returns:
        tuple: longest_side (int), length (float), right_angle (bool)
    """
    A = _np.linalg.norm([triangle[1] - triangle[0]])
    B = _np.linalg.norm([triangle[2] - triangle[1]])
    C = _np.linalg.norm([triangle[2] - triangle[0]])
    sides = _np.array([[0, 1, 2], [A, B, C]]).T
    sorted_sides = sides[sides[:, 1].argsort()][::-1]
    longest_side = int(sorted_sides[0, 0])

    # Get sizes of two right angled triangles

    # Get Altitude to longest side
    alt_side = altitude(sorted_sides[0, 1], sorted_sides[1, 1], sorted_sides[2, 1])
    left_side = sides[(longest_side + 1) % 3, 1]
    right_side = sides[(longest_side - 1) % 3, 1]

    p = _np.sqrt(left_side ** 2 - alt_side ** 2)
    q = _np.sqrt(right_side ** 2 - alt_side ** 2)

    RA_triangle1 = _np.array([p, alt_side])
    RA_triangle2 = _np.array([q, alt_side])

    return longest_side, RA_triangle1, RA_triangle2


def check_sign(vector_1, vector_2):
    """Returns true if the signs of all elements of two arrays are the same

    Args:
        vector_1 (ndarray): input array 2
        vector_2 (ndarray): input array 2

    Returns:
        boolean: True if elements in two arrays have the same sign
    """
    sign_comp_1 = _np.fabs(vector_1 + vector_2)
    sign_comp_2 = _np.fabs(vector_1) + _np.fabs(vector_2)

    return _np.allclose(sign_comp_1, sign_comp_2, atol=1e-6)


def return_axis_vector(triangle, longest_side):
    """Returns vector collinear to the longest side of a triangle

    Args:
        triangle (ndarray): triangle vertices
        longest_side (int): longest side index

    Returns:
        ndarray: vector corresponding to longest side
    """
    vector_A = triangle[1] - triangle[0]
    vector_A = vector_A / _np.linalg.norm(vector_A)

    vector_B = triangle[2] - triangle[1]
    vector_B = vector_B / _np.linalg.norm(vector_B)
    vector_C = triangle[2] - triangle[0]
    vector_C = vector_C / _np.linalg.norm(vector_C)

    vec_dict = {
        0: vector_A,
        1: vector_B,
        2: vector_C,
    }

    vec = vec_dict[longest_side]

    return vec


def return_z_vector(triangle, longest_side):
    """Returns altitude vector from longest side

    Args:
        triangle (ndarray): triangle vertices
        longest_side (int): longest side index

    Returns:
        ndarray: altitude vector
    """
    new_vertex_dict = {
        0: [triangle[2, 0], triangle[0, 1], triangle[0, 2]],
        1: [triangle[0, 0], triangle[0, 1], triangle[1, 2]],
        2: [triangle[1, 0], triangle[2, 1], triangle[2, 2]],
    }
    new_vertex = _np.array(new_vertex_dict[longest_side])
    vec_z_dict = {
        0: triangle[2] - new_vertex,
        1: triangle[0] - new_vertex,
        2: triangle[1] - new_vertex,
    }
    vec_z = vec_z_dict[longest_side]
    vec_z = vec_z / _np.linalg.norm(vec_z)

    return vec_z


def align_triangle_to_y(triangle, rot_axis, norm_vec):
    """Rotates and translates a triangle in lie in the xz plane

    Args:
        triangle (ndarray): vertices of a triangle
        rot_axis (ndarray): axis about which to rotate triangle
        norm_vec (ndarray): normal to triangle

    Returns:
        tuple: aligned_triangle (ndarray, rotated triangle), first_rotation (quaternion, align to y-axis)
    """
    y_axis = _np.array([0, 1, 0])

    if _np.linalg.norm(rot_axis) < ALIGN_CUTOFF:

        # Check if parallel or anti-parallel
        if check_sign(y_axis, norm_vec):
            # Parallel
            first_rotation = Quaternion()
            aligned_triangle = triangle

        else:
            # Anti-parallel
            first_rotation = Quaternion.q_angle_from_axis(PI, y_axis)
            aligned_triangle = rotate_points(triangle, first_rotation)

    else:
        angle = -_np.arccos(_np.dot(y_axis, norm_vec))
        first_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)
        aligned_triangle = rotate_points(triangle, first_rotation)
    return aligned_triangle, first_rotation


def align_triangle_xz(triangle, longest_side):
    """Returns quaternions needed to rotate a triangle lying in the xz plane to be
    aligned with its longest side along x and altitude along z

    Args:
        triangle (ndarray): vertices of triangle
        longest_side (int): index of the longest side of triangle

    Returns:
        tuple: second_rotation (quaternion, align to x-axis), third_rotation (quaternion, align to z-axis):
    """
    x_axis = _np.array([1, 0, 0])
    y_axis = _np.array([0, 1, 0])
    z_axis = _np.array([0, 0, 1])

    # side_list = [0, 1, 2]
    # side_list.pop(longest_side)

    vec_x = return_axis_vector(triangle, longest_side)
    rot_axis = _np.cross(x_axis, vec_x)

    # Check aligment of base of triangle with x-axis
    if _np.linalg.norm(rot_axis) < ALIGN_CUTOFF:

        # Check if parallel or anti-parallel
        if check_sign(x_axis, vec_x):
            # Parallel
            second_rotation = Quaternion()
            tri_x = triangle
        else:
            # Anti-parallel
            second_rotation = Quaternion.q_angle_from_axis(PI, y_axis)
            tri_x = rotate_points(triangle, second_rotation)

    else:
        angle = -_np.arccos(_np.dot(x_axis, vec_x))
        second_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)
        tri_x = rotate_points(triangle, second_rotation)

    vec_z = return_z_vector(tri_x, longest_side)
    rot_axis = _np.cross(z_axis, vec_z)

    # Check aligment of triangle altitude with z-axis
    if _np.all(_np.fabs([rot_axis]) < ALIGN_CUTOFF):

        # Check if parallel anti-parallel
        if check_sign(z_axis, vec_z):
            # Parallel
            third_rotation = Quaternion()
        else:
            # Anti-parallel
            third_rotation = Quaternion.q_angle_from_axis(PI, y_axis)

    else:
        angle = -_np.arccos(_np.dot(z_axis, vec_z))
        third_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)

    return second_rotation, third_rotation


def _rotate_triangle(triangle, Jr):
    """Gets rotation angles needed for transformation of coordinates from
    global frame to a local frame

    Args:
        triangle (ndarray): vertices of a triangular plane

    Returns:
        tuple: rotated_triangle (ndarray), right_angled (Boolean), total_rotation (Quaternion)
    """

    loc_triangle = triangle - triangle[0]

    y_axis = _np.array([0, 1, 0])

    longest_side, RA_triangle1, RA_triangle2 = _largest_side_RA(loc_triangle)

    # align norm to y axis
    norm_vec_y = norm_plane(loc_triangle)
    rot_axis = _np.cross(y_axis, norm_vec_y)

    aligned_triangle, first_rotation = align_triangle_to_y(
        loc_triangle, rot_axis, norm_vec_y
    )

    second_rotation, third_rotation = align_triangle_xz(aligned_triangle, longest_side)
    total_rotation = third_rotation * second_rotation * first_rotation

    rotated_triangle = rotate_points(triangle, total_rotation)

    origin_vertex = _np.argwhere(
        rotated_triangle[:, 0] == rotated_triangle[:, 0].min()
    ).ravel()[0]

    offset = rotated_triangle[origin_vertex]

    return total_rotation, rotated_triangle, offset, RA_triangle1, RA_triangle2