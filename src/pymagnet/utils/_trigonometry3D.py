# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Contains functions needed to rotate and translate a triangle to lie in the xz plane
and to divide it into two right angled triangles
"""

import numpy as _np
from numba import jit

from ._quaternion import Quaternion, q_angle_from_axis
from .global_const import ALIGN_CUTOFF, PI

# Module-level constants to avoid repeated array creation
_X_AXIS = _np.array([1.0, 0.0, 0.0])
_Y_AXIS = _np.array([0.0, 1.0, 0.0])
_Z_AXIS = _np.array([0.0, 0.0, 1.0])

# Vertex pair indices for each edge: edge 0 = v0→v1, edge 1 = v1→v2, edge 2 = v0→v2
_EDGE_VERTICES = ((0, 1), (1, 2), (0, 2))

# Opposite vertex for each edge
_OPPOSITE_VERTEX = {0: 2, 1: 0, 2: 1}


def _safe_arccos(x):
    """Compute arccos with domain clamping to prevent NaN from floating-point errors.

    Args:
        x (float or ndarray): Input value(s), should be in [-1, 1]

    Returns:
        float or ndarray: arccos(clip(x, -1, 1))
    """
    return _np.arccos(_np.clip(x, -1.0, 1.0))


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
    """Calculates the normal to a triangular plane.

    Args:
        vec (ndarray): (3,3) array of triangle vertices

    Returns:
        ndarray: unit normal vector (3,)

    Raises:
        ValueError: If triangle is degenerate (collinear vertices)
    """
    norm = _np.cross(vec[1] - vec[0], vec[2] - vec[0])
    length = _np.linalg.norm(norm)

    if length < 1e-10:
        raise ValueError("Degenerate triangle: vertices are collinear")

    return norm / length


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
    """Gets altitude to side `a` of a triangle using Heron's formula.

    Args:
        a (float): base side (altitude is perpendicular to this)
        b (float): second side
        c (float): third side

    Returns:
        float: altitude to side `a`

    Raises:
        ValueError: If triangle is degenerate (sides violate triangle inequality)
    """
    if a < 1e-10:
        raise ValueError("Degenerate triangle: base side has zero length")

    s = (a + b + c) / 2
    radicand = s * (s - a) * (s - b) * (s - c)

    if radicand < 0:
        raise ValueError(
            f"Degenerate triangle: sides ({a:.6f}, {b:.6f}, {c:.6f}) "
            "violate triangle inequality"
        )

    return 2 * _np.sqrt(radicand) / a


def _largest_side_RA(triangle):
    """Determines largest side of triangle and decomposes into right-angled triangles.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices

    Returns:
        tuple: (longest_side, RA_triangle1, RA_triangle2)
            - longest_side (int): index of longest edge (0, 1, or 2)
            - RA_triangle1 (ndarray): [base, altitude] of first right triangle
            - RA_triangle2 (ndarray): [base, altitude] of second right triangle
    """
    # Compute edge lengths efficiently
    edge_lengths = _np.array([
        _np.linalg.norm(triangle[1] - triangle[0]),  # Edge 0: v0 → v1
        _np.linalg.norm(triangle[2] - triangle[1]),  # Edge 1: v1 → v2
        _np.linalg.norm(triangle[2] - triangle[0]),  # Edge 2: v0 → v2
    ])

    # Use argmax (O(n)) instead of argsort (O(n log n))
    longest_side = int(_np.argmax(edge_lengths))

    # Get lengths for altitude calculation
    a = edge_lengths[longest_side]  # longest side
    b = edge_lengths[(longest_side + 1) % 3]
    c = edge_lengths[(longest_side + 2) % 3]

    # Altitude from opposite vertex to longest side (Heron's formula)
    alt_side = altitude(a, b, c)

    # Decompose: drop altitude from apex to base
    # Creates two right triangles with legs (p, h) and (q, h) where p + q = a
    # Using Pythagorean theorem: p² + h² = b², q² + h² = c²
    left_side = edge_lengths[(longest_side + 1) % 3]
    right_side = edge_lengths[(longest_side - 1) % 3]

    p = _np.sqrt(_np.maximum(left_side**2 - alt_side**2, 0))
    q = _np.sqrt(_np.maximum(right_side**2 - alt_side**2, 0))

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
    """Returns unit vector along the longest side of a triangle.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices
        longest_side (int): index of longest edge (0, 1, or 2)

    Returns:
        ndarray: unit vector (3,) along longest edge
    """
    i, j = _EDGE_VERTICES[longest_side]
    vec = triangle[j] - triangle[i]
    return vec / _np.linalg.norm(vec)


def return_z_vector(triangle, longest_side):
    """Returns unit altitude vector from longest side toward opposite vertex.

    For a triangle aligned in the xz plane with longest side along x,
    this returns the direction toward the apex (along z).

    Args:
        triangle (ndarray): (3,3) array of triangle vertices
        longest_side (int): index of longest edge (0, 1, or 2)

    Returns:
        ndarray: unit altitude vector (3,)

    Raises:
        ValueError: If altitude has zero length (degenerate triangle)
    """
    # Get the vertex opposite to the longest side
    opposite_idx = _OPPOSITE_VERTEX[longest_side]

    # Get endpoints of the longest side
    start_idx, end_idx = _EDGE_VERTICES[longest_side]

    # Project opposite vertex onto the line of longest side
    edge_vec = triangle[end_idx] - triangle[start_idx]
    edge_unit = edge_vec / _np.linalg.norm(edge_vec)

    to_opposite = triangle[opposite_idx] - triangle[start_idx]
    projection_length = _np.dot(to_opposite, edge_unit)
    foot_of_altitude = triangle[start_idx] + projection_length * edge_unit

    # Altitude vector: from foot to opposite vertex
    altitude_vec = triangle[opposite_idx] - foot_of_altitude
    length = _np.linalg.norm(altitude_vec)

    if length < 1e-10:
        raise ValueError("Degenerate triangle: altitude has zero length")

    return altitude_vec / length


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
            first_rotation = q_angle_from_axis(PI, y_axis)
            aligned_triangle = rotate_points(triangle, first_rotation)

    else:
        angle = -_safe_arccos(_np.dot(y_axis, norm_vec))
        first_rotation = q_angle_from_axis(angle, rot_axis)
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
            second_rotation = q_angle_from_axis(PI, y_axis)
            tri_x = rotate_points(triangle, second_rotation)

    else:
        angle = -_safe_arccos(_np.dot(x_axis, vec_x))
        second_rotation = q_angle_from_axis(angle, rot_axis)
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
            third_rotation = q_angle_from_axis(PI, y_axis)

    else:
        angle = -_safe_arccos(_np.dot(z_axis, vec_z))
        third_rotation = q_angle_from_axis(angle, rot_axis)

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
