import numpy as _np
from numpy.lib.twodim_base import tri

# from ._magnet import Magnet
from ._magnet3 import Magnet_3D
from ._routines3 import Vector3
from ._quaternion import Quaternion

PI = _np.pi
u0 = PI * 4e-7


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
    area = 0

    for i in range(NP):
        j = j % NP

        area += (triangle[j][0] - triangle[i][0]) * (triangle[j][2] + triangle[i][2])
        j += 1

    # check winding order of polygon, area < 0 for clockwise ordering of points
    area /= 2.0

    return area


def _charge_sheet(a, b, Jr, x, y, z):
    """Calculates the magnetic field of a right angled charge sheet

    Args:
        a ([type]): [description]
        b ([type]): [description]
        Jr ([type]): [description]
        x ([type]): [description]
        y ([type]): [description]
        z ([type]): [description]

    Returns:
        tuple: Bx, By, Bz
    """

    sig = Jr / (4 * PI)

    # ----Bx, Bz

    S_ab = 1 / _np.sqrt(a ** 2 + b ** 2)
    a_ab = _np.sqrt(1 + b ** 2 / a ** 2)
    r1 = _np.sqrt(x ** 2 + y ** 2 + z ** 2)

    ax = a - x
    bz = b - z

    r2 = _np.sqrt(ax ** 2 + y ** 2 + bz ** 2)
    r3 = _np.sqrt(ax ** 2 + y ** 2 + z ** 2)

    t1 = _np.log(r1 - (a * x + b * z) * S_ab)
    t2 = _np.log(r2 + (a * ax + b * bz) * S_ab)

    Bx = sig * (b * S_ab * (t1 - t2) + _np.log((r2 + bz) / (r3 - z)))
    Bx[(r2 + bz) == 0] = 0
    Bx[_np.isnan(Bx)] = 0

    Bz = sig * (a * S_ab * (t2 - t1) + _np.log((r1 - x) / (r3 + ax)))
    Bz[(r1 - x) == 0] = 0
    Bz[_np.isnan(Bz)] = 0

    # ----Bz
    g = b / a
    si_alpha = 1.0 / (a_ab ** 2)
    beta = -(x + z * g) * si_alpha
    gamma = _np.sqrt((r1 ** 2 * si_alpha) - (beta ** 2))

    A = -gamma * g
    B = gamma * a_ab
    C = z + beta * g

    t3 = (a + beta) / gamma
    t4 = beta / gamma

    t1 = _np.sqrt(B ** 2 - A ** 2 - C ** 2)
    t2 = 1.0 / t1
    atan_1 = t2 * (C + (A + B) * (_np.sqrt(1 + t3 * t3) + t3))
    atan_2 = t2 * (C + (A + B) * (_np.sqrt(1 + t4 * t4) + t4))
    By = (y / (a_ab * t1)) * _np.arctan2(atan_1 - atan_2, 1 + atan_1 * atan_2)

    t2 = 1.0 / y
    atan_1 = t2 * (r3 + z + (x - a))
    atan_2 = t2 * (r1 + z + x)
    By = By + _np.arctan2(atan_1 - atan_2, 1 + atan_1 * atan_2)
    # By[y == 0] = 0
    By[(r2 + bz) == 0] = 0
    By = sig * By
    By[_np.isnan(By)] = 0

    # Bn = _np.sqrt(Bx ** 2 + By ** 2 + Bz ** 2)
    return Bx, By, Bz


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


def gen_rotation_quaternion(alpha, beta, gamma):
    """Generates quaternion for rotation around z, y, x axes

    Args:
        alpha ([type]): [description]
        beta ([type]): [description]
        gamma ([type]): [description]

    Returns:
        [type]: [description]
    """
    from functools import reduce

    q_z = Quaternion.q_angle_from_axis(_np.deg2rad(alpha), (0, 0, 1))
    q_y = Quaternion.q_angle_from_axis(_np.deg2rad(beta), (0, 1, 0))
    q_x = Quaternion.q_angle_from_axis(_np.deg2rad(gamma), (1, 0, 0))

    q_forward = [q for q in [q_z, q_y, q_x] if q is not None]

    forward_rotation = reduce(lambda x, y: x * y, q_forward)

    q_reverse = [q.get_conjugate() for q in [q_x, q_y, q_z] if q is not None]

    reverse_rotation = reduce(lambda x, y: x * y, q_reverse)
    return forward_rotation, reverse_rotation


def rot_points(points, rotation_quaternion):
    """Rotates a set of points

    Args:
        points ([type]): [description]
        rotation_quaternion ([type]): [description]

    Returns:
        [type]: [description]
    """
    x_rot, y_rot, z_rot = rotation_quaternion * points.T
    rot_points = _np.vstack([x_rot, y_rot, z_rot]).T

    return rot_points


def _gen_two_RA_triangles(triangle):
    """Splits a triangle into two right angled triangles

    Args:
        triangle (ndarray): 3x3 array of vertices

    Returns:
        tuple: tuple of the lengths a and b of two sides of each triangle
    """
    A = _np.linalg.norm([triangle[1] - triangle[0]])
    B = _np.linalg.norm([triangle[2] - triangle[1]])
    C = _np.linalg.norm([triangle[2] - triangle[0]])

    new_vertex = [triangle[1, 0], triangle[0, 1], triangle[0, 2]]

    A1 = A

    B1 = _np.linalg.norm(triangle[1] - new_vertex)
    C1 = _np.linalg.norm(triangle[0] - new_vertex)
    A2 = B
    B2 = B1
    C2 = C - C1

    triangle_1 = _np.array([C1, B1])
    triangle_2 = _np.array([C2, B2])

    return triangle_1, triangle_2


def _calcB_2_triangles(triangle1, triangle2, Jr, x, y, z):
    """Calculates the magnetic field due to two split right angled triangles
    in their local frame.

    Args:
        triangle1 (ndarray): Vertices of triangle 1
        triangle2 (ndarray): Vertices of triangle 2
        Jr (float): normal remnant magnetisation
        x (ndarray): x coordinates
        y (ndarray): y coordinates
        z (ndarray): z coordinates

    Returns:
        tuple: Bx, By, Bz magnetic field components
    """

    # Calc RA1 Field
    Btx, Bty, Btz = _charge_sheet(triangle1[0], triangle1[1], Jr, x, y, z)

    # Rotate into local of RA2
    rotate_about_z = Quaternion.q_angle_from_axis(PI, (0, 0, 1))
    pos_vec_RA2 = Quaternion._prepare_vector(x - triangle1[0], y, z)
    x_local, y_local, z_local = rotate_about_z * pos_vec_RA2

    # Calc RA2 Field
    Btx2, Bty2, Btz2 = _charge_sheet(
        triangle2[0], triangle2[1], Jr, x_local + triangle1[0], y_local, z_local
    )

    # Inverse Rot of RA2 Field
    Bvec = Quaternion._prepare_vector(Btx2, Bty2, Btz2)
    Btx2, Bty2, Btz2 = rotate_about_z.get_conjugate() * Bvec

    Btx += Btx2
    Bty += Bty2
    Btz += Btz2

    return Btx, Bty, Btz


def largest_side_RA(triangle):
    """Determines largest side of triangle and if it is a right angled triangle

    Args:
        triangle (ndarray): triangle vertices

    Returns:
        tuple: longest_side (int), length (float), right_angle (bool)
    """
    A = _np.linalg.norm([triangle[1] - triangle[0]])
    B = _np.linalg.norm([triangle[2] - triangle[1]])
    C = _np.linalg.norm([triangle[2] - triangle[0]])
    sides = _np.array([[0, 1, 2], [A, B, C]]).T
    sides = sides[sides[:, 1].argsort()]
    longest_side = int(sides[-1, 0])
    length = sides[-1, 1]
    right_angle = _np.allclose(sides[0, 1] ** 2 + sides[1, 1] ** 2, sides[2, 1] ** 2)
    return longest_side, length, right_angle


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

    return _np.allclose(sign_comp_1, sign_comp_2)


def return_axis_vector(triangle, longest_side):
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
    new_vertex_dict = {
        0: [triangle[2, 0], triangle[0, 1], triangle[0, 2]],
        1: [triangle[0, 0], triangle[1, 1], triangle[1, 2]],
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

    return vec_z, new_vertex


def align_triangle_to_y(triangle, rot_axis, norm_vec):

    y_axis = _np.array([0, 1, 0])
    null_vector = _np.array([0, 0, 0])

    if _np.allclose(rot_axis, null_vector):

        # Check if anti-parallel
        if check_sign(y_axis, norm_vec):
            # Parallel
            first_rotation = Quaternion()
            aligned_triangle = triangle.copy()

        else:
            # Anti-parallel
            first_rotation = Quaternion.q_angle_from_axis(PI, y_axis)
            aligned_triangle = rot_points(triangle, first_rotation)

    else:
        angle = -_np.arccos(_np.dot(y_axis, norm_vec))
        first_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)
        aligned_triangle = rot_points(triangle, first_rotation)
    return aligned_triangle, first_rotation


def align_triangle_xz(triangle, longest_side):

    x_axis = _np.array([1, 0, 0])
    y_axis = _np.array([0, 1, 0])
    z_axis = _np.array([0, 0, 1])
    null_vector = _np.array([0, 0, 0])
    side_list = [0, 1, 2]
    side_list.pop(longest_side)

    vec_x = return_axis_vector(triangle, longest_side)
    rot_axis = _np.cross(x_axis, vec_x)

    if _np.allclose(rot_axis, null_vector):
        # Check if anti-parallel
        if check_sign(x_axis, vec_x):
            # Parallel
            second_rotation = Quaternion()
            tri_x = triangle.copy()
        else:
            # Anti-parallel
            second_rotation = Quaternion.q_angle_from_axis(PI, x_axis)
            tri_x = rot_points(triangle, second_rotation)

    else:
        angle = -_np.arccos(_np.dot(x_axis, vec_x))
        second_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)
        tri_x = rot_points(triangle, second_rotation)

    vec_z, new_vertex = return_z_vector(tri_x.copy(), longest_side)

    rot_axis = _np.cross(z_axis, vec_z)
    #
    if _np.allclose(rot_axis, null_vector):
        # Check if anti-parallel
        if check_sign(z_axis, vec_z):
            # Parallel
            third_rotation = Quaternion()
            tri_z = tri_x.copy()
        else:
            # Anti-parallel
            third_rotation = Quaternion.q_angle_from_axis(PI, y_axis)
            tri_z = rot_points(tri_x.copy(), third_rotation)

    else:
        angle = -_np.arccos(_np.dot(z_axis, vec_z))
        third_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)
        tri_z = rot_points(tri_x.copy(), third_rotation)
    print(tri_z)
    return second_rotation, third_rotation, new_vertex


def _rotate_triangle(triangle):
    """Gets rotation angles needed for transformation of coordinates from
    global frame to a local frame

    Args:
        triangle (ndarray): vertices of a triangular plane

    Returns:
        tuple: rotated_triangle (ndarray), right_angled (Boolean), total_rotation (Quaternion)
    """

    triangle -= triangle[0]

    x_axis = _np.array([1, 0, 0])
    y_axis = _np.array([0, 1, 0])
    z_axis = _np.array([0, 0, 1])

    null_vector = _np.array([0, 0, 0])

    longest_side, _, _ = largest_side_RA(triangle)

    # align norm to y axis
    norm_vec = norm_plane(triangle)
    rot_axis = _np.cross(y_axis, norm_vec)

    aligned_triangle, first_rotation = align_triangle_to_y(triangle, rot_axis, norm_vec)

    second_rotation, third_rotation, _ = align_triangle_xz(
        aligned_triangle, longest_side
    )
    total_rotation = third_rotation * second_rotation * first_rotation

    rotated_triangle = rot_points(triangle, total_rotation)

    # Reorder vertices so origin vertex is first
    origin_vertex = (longest_side + 1) % 3
    rotated_triangle = rotated_triangle[
        [(origin_vertex + 0) % 3, (origin_vertex + 1) % 3, (origin_vertex + 2) % 3], :
    ]

    return rotated_triangle, total_rotation


def calcB_triangle(triangle, Jr, x, y, z):
    """Calculates the magnetic field due to a triangle

    Args:
        triangle2(ndarray): Vertices of triangle
        Jr (float): normal remnant magnetisation
        x (ndarray): x coordinates
        y (ndarray): y coordinates
        z (ndarray): z coordinates

    Returns:
        tuple: Bx, By, Bz magnetic field components
    """

    loc_triangle, q_fwd = _rotate_triangle(triangle)

    # Prepare points and quaternion
    pos_vec = Quaternion._prepare_vector(
        x - triangle[0, 0], y - triangle[0, 1], z - triangle[0, 2]
    )
    # Rotate points
    x_rot, y_rot, z_rot = q_fwd * pos_vec

    # loc_triangle, q_fwd = _rotate_triangle(triangle)
    loc_RA_1, loc_RA_2 = _gen_two_RA_triangles(loc_triangle)

    Btx, Bty, Btz = _calcB_2_triangles(loc_RA_1, loc_RA_2, Jr, x_rot, y_rot, z_rot)

    Bvec = Quaternion._prepare_vector(Btx, Bty, Btz)
    Bx, By, Bz = q_fwd.get_conjugate() * Bvec

    return Bx, By, Bz