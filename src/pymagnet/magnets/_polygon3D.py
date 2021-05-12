import numpy as _np

from ._magnet3 import Magnet_3D
from ..utils._routines3D import Vector3
from ..utils._quaternion import Quaternion

from stl import mesh
import plotly.graph_objects as _go

from numba import jit, vectorize, float64
from math import sqrt, log, fabs, atan2

from ..utils.global_const import PI, FP_CUTOFF, ALIGN_CUTOFF


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


def largest_side_RA(triangle):
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

    return vec_z, new_vertex


def align_triangle_to_y(triangle, rot_axis, norm_vec):

    y_axis = _np.array([0, 1, 0])
    # null_vector = _np.array([0, 0, 0])

    # if _np.all(_np.fabs([rot_axis]) < ALIGN_CUTOFF):
    if _np.linalg.norm(rot_axis) < ALIGN_CUTOFF:

        # Check if anti-parallel
        if check_sign(y_axis, norm_vec):
            # Parallel
            first_rotation = Quaternion()
            aligned_triangle = triangle

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

    side_list = [0, 1, 2]
    side_list.pop(longest_side)

    vec_x = return_axis_vector(triangle, longest_side)
    rot_axis = _np.cross(x_axis, vec_x)

    # if _np.allclose(rot_axis, null_vector, atol=1e-4):
    # if _np.all(_np.fabs([rot_axis]) < ALIGN_CUTOFF):
    if _np.linalg.norm(rot_axis) < ALIGN_CUTOFF:

        # Check if anti-parallel
        if check_sign(x_axis, vec_x):
            # Parallel
            second_rotation = Quaternion()
            tri_x = triangle
        else:
            # Anti-parallel
            second_rotation = Quaternion.q_angle_from_axis(PI, y_axis)
            tri_x = rot_points(triangle, second_rotation)

    else:
        angle = -_np.arccos(_np.dot(x_axis, vec_x))
        second_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)
        tri_x = rot_points(triangle, second_rotation)

    vec_z, new_vertex = return_z_vector(tri_x, longest_side)
    rot_axis = _np.cross(z_axis, vec_z)

    if _np.all(_np.fabs([rot_axis]) < ALIGN_CUTOFF):
        # Check if anti-parallel
        if check_sign(z_axis, vec_z):
            # Parallel
            third_rotation = Quaternion()
        else:
            # Anti-parallel
            third_rotation = Quaternion.q_angle_from_axis(PI, y_axis)

    else:
        angle = -_np.arccos(_np.dot(z_axis, vec_z))
        third_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)

    return second_rotation, third_rotation, new_vertex


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

    longest_side, RA_triangle1, RA_triangle2 = largest_side_RA(loc_triangle)

    # align norm to y axis
    norm_vec_y = norm_plane(loc_triangle)
    rot_axis = _np.cross(y_axis, norm_vec_y)

    aligned_triangle, first_rotation = align_triangle_to_y(
        loc_triangle, rot_axis, norm_vec_y
    )

    second_rotation, third_rotation, _ = align_triangle_xz(
        aligned_triangle, longest_side
    )
    total_rotation = third_rotation * second_rotation * first_rotation

    rotated_triangle = rot_points(triangle, total_rotation)
    origin_vertex = _np.argwhere(
        rotated_triangle[:, 0] == rotated_triangle[:, 0].min()
    ).ravel()[0]

    offset = rotated_triangle[origin_vertex]

    return total_rotation, rotated_triangle, offset, RA_triangle1, RA_triangle2


class Mesh(Magnet_3D):
    """Mesh Magnet Class."""

    mag_type = "Mesh"

    def __init__(
        self,
        filename,
        Jr=1.0,  # local magnetisation
        **kwargs,
    ):

        super().__init__(Jr, **kwargs)

        self.phi = kwargs.pop("theta", 90)
        self.phi_rad = _np.deg2rad(self.phi)
        self.theta = kwargs.pop("phi", 0)
        self.theta_rad = _np.deg2rad(self.theta)

        self.mesh_scale = kwargs.pop("mesh_scale", 1.0)
        self._filename = filename

        self.mesh_vectors, self.mesh_normals = self._import_mesh()

        self.Jx = _np.around(
            Jr * _np.cos(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        )
        self.Jy = _np.around(
            Jr * _np.sin(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        )
        self.Jz = _np.around(Jr * _np.cos(self.theta_rad), decimals=6)
        self.tol = 1e-4  # sufficient for 0.01 degree accuracy

        self.J = _np.array([self.Jx, self.Jy, self.Jz])

        self.Jnorm = _np.dot(self.J, self.mesh_normals.T)

        self.start = kwargs.pop("start", 0)
        self.start = _np.min([self.start, len(self.mesh_vectors) - 1])
        self.stop = kwargs.pop("stop", len(self.mesh_vectors))
        self.stop = _np.min([self.stop, len(self.mesh_vectors)])

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def get_Jr(self):
        return self.J

    def size(self):
        """Returns magnet dimesions

        Returns:
            size[ndarray]: numpy array [width, depth, height]
        """
        # return _np.array([self.width, self.depth, self.height])
        pass

    def calcB(self, x, y, z):
        """Calculates the magnetic field at point(s) x,y,z due to a 3D magnet
        The calculations are always performed in local coordinates with the centre of the magnet at origin and z magnetisation pointing along the local z' axis.

        The rotations and translations are performed first, and the internal field calculation functions are called.

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates

        Returns:
            tuple: Bx(ndarray), By(ndarray), Bz(ndarray)  field vector
        """
        # from ._routines3 import _tile_arrays, _apply_mask

        B = self._calcB_local(x, y, z)

        # xloc, yloc, zloc = _tile_arrays(x - self.xc, y - self.yc, z - self.zc)
        # mask = self._generate_mask(xloc, yloc, zloc)
        # B = _apply_mask(self, B, mask)

        return B.x, B.y, B.z

    def _calcB_local(self, x, y, z):
        """Internal magnetic field calculation methods.
        Iterates over each component of the prism magnetised in x, y, and z
        (in local coordinates).

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates

        Returns:
            Vector3: Magnetic field array
        """
        from ..utils._routines3D import _allocate_field_array3

        B = _allocate_field_array3(x, y, z)
        vec_shape = B.x.shape
        B.x = B.x.ravel()
        B.y = B.y.ravel()
        B.z = B.z.ravel()

        for i in range(self.start, self.stop):
            # for i in range(len(self.mesh_vectors)):
            if _np.fabs(self.Jnorm[i] / self.Jr) > 1e-4:
                # triangle = self.mesh_vectors[i].copy()
                Btx, Bty, Btz, _, _, _ = self.calcB_triangle(
                    self.mesh_vectors[i],
                    # 1.0,  # FIXME: Temporary change for testing values
                    self.Jnorm[i],
                    x,
                    y,
                    z,
                    i,
                )

                B.x += Btx
                B.y += Bty
                B.z += Btz

        B.x = _np.reshape(B.x, vec_shape)
        B.y = _np.reshape(B.y, vec_shape)
        B.z = _np.reshape(B.z, vec_shape)

        B.n = _np.linalg.norm([B.x, B.y, B.z], axis=0)
        return B

    #     def _import_mesh(self):
    #         plot_data, magnet_mesh = self.stl2mesh3d()
    #         magnet_mesh.vectors *= self.mesh_scale
    #
    #         magnet_mesh.normals = magnet_mesh.normals / _np.linalg.norm(
    #             magnet_mesh.normals, axis=1, keepdims=True
    #         )
    #         return plot_data, magnet_mesh

    def _import_mesh(self):
        stl_mesh = mesh.Mesh.from_file(self._filename)

        offset = self.center()

        stl_mesh.translate(offset / self.mesh_scale)

        if _np.any(_np.fabs([self.alpha_rad, self.beta_rad, self.gamma_rad]) > 1e-5):
            mesh_rotation = Quaternion.gen_rotation_quaternion(
                self.alpha_rad, self.beta_rad, self.gamma_rad
            )

            angle, axis = mesh_rotation.get_axisangle()
            stl_mesh.rotate(axis, angle)

        mesh_vectors = stl_mesh.vectors.astype(_np.float64)
        mesh_normals = stl_mesh.normals.astype(_np.float64)

        mesh_vectors *= self.mesh_scale

        mesh_normals = mesh_normals / _np.linalg.norm(
            mesh_normals, axis=1, keepdims=True
        )

        return mesh_vectors, mesh_normals

    def _generate_mask(self, x, y, z):
        """Generates mask of points inside a magnet

        Args:
            x (ndarray/float): x-coordinates
            y (ndarray/float): y-coordinates
            z (ndarray/float): z-coordinates
        """
        pass

    def calcB_triangle(self, triangle, Jr, x, y, z, i):
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

        (
            total_rotation,
            rotated_triangle,
            offset,
            RA_triangle1,
            RA_triangle2,
        ) = _rotate_triangle(triangle, Jr)

        # Prepare points and quaternion
        pos_vec = Quaternion._prepare_vector(x, y, z)

        # Rotate points
        x_rot, y_rot, z_rot = total_rotation * pos_vec

        # FIXME:
        norm1 = norm_plane(triangle)

        # if (
        #     _np.fabs(norm1[1]) < 1e-4 > 0
        #     and _np.all(norm1[0::2] > ALIGN_CUTOFF)
        #     and Jr > 0
        # ):
        #     print("Trig 1", i)
        # RA_triangle1, RA_triangle2 = RA_triangle2, RA_triangle1
        #     if (
        #         not _np.allclose(_np.fabs(norm1), [0, 0, 1], atol=1e-5)
        #         and not _np.allclose(_np.fabs(norm1), [1, 0, 0], atol=1e-5)
        #         and not _np.allclose(_np.fabs(norm1), [1, 0, 1], atol=1e-5)
        #     ):
        #         print("Trig 2")

        if _np.allclose(norm1, [0, -1, 0], atol=1e-5) and Jr < 0:
            # print("Trig 2", i)
            RA_triangle1, RA_triangle2 = RA_triangle2, RA_triangle1

        # pass
        # print(f"{i},{norm1[0]},{norm1[1]},{norm1[2]},{Jr}")

        Btx, Bty, Btz = self._calcB_2_triangles(
            RA_triangle1,
            RA_triangle2,
            Jr,
            x_rot - offset[0],
            y_rot - offset[1],
            z_rot - offset[2],
        )

        Bvec = Quaternion._prepare_vector(Btx, Bty, Btz)
        Bx, By, Bz = total_rotation.get_conjugate() * Bvec

        return Bx, By, Bz, rotated_triangle, offset, total_rotation

    def _calcB_2_triangles(self, triangle1, triangle2, Jr, x, y, z):
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
        Btx, Bty, Btz = self._charge_sheet(triangle1[0], triangle1[1], Jr, x, y, z)

        # Rotate into local of RA2
        rotate_about_z = Quaternion.q_angle_from_axis(PI, (0, 0, 1))
        pos_vec_RA2 = Quaternion._prepare_vector(x - triangle1[0], y, z)

        x_local, y_local, z_local = rotate_about_z * pos_vec_RA2
        # print("Jr:", Jr)
        # Calc RA2 Field
        Btx2, Bty2, Btz2 = self._charge_sheet(
            triangle2[0], triangle2[1], Jr, x_local + triangle2[0], y_local, z_local
        )

        # Inverse Rot of RA2 Field
        Bvec = Quaternion._prepare_vector(Btx2, Bty2, Btz2)
        Btx2, Bty2, Btz2 = rotate_about_z.get_conjugate() * Bvec

        Btx += Btx2
        Bty += Bty2
        Btz += Btz2

        return Btx, Bty, Btz

    @staticmethod
    def _charge_sheet(a, b, Jr, x, y, z):
        sigma = Jr / (4 * PI)
        with _np.errstate(all="ignore"):
            Bx = _charge_sheet_x(a, b, sigma, x, y, z)
            By = _charge_sheet_y(a, b, sigma, x, y, z)
            Bz = _charge_sheet_z(a, b, sigma, x, y, z)
        return Bx, By, Bz


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64)], target="parallel"
)
def _charge_sheet_x(a, b, sigma, x, y, z):
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

    S_ab = 1.0 / sqrt(a ** 2 + b ** 2)

    r1 = sqrt(x ** 2 + y ** 2 + z ** 2)

    ax = a - x
    bz = b - z

    r2 = sqrt(ax ** 2 + y ** 2 + bz ** 2)
    r3 = sqrt(ax ** 2 + y ** 2 + z ** 2)
    t1_log = r1 - (a * x + b * z) * S_ab

    t2_log = r2 + (a * ax + b * bz) * S_ab

    if fabs(t1_log) > FP_CUTOFF:
        t1 = log(t1_log)
    else:
        t1 = 0.0

    if fabs(t2_log) > FP_CUTOFF:
        t2 = log(t2_log)
    else:
        t2 = 0.0

    dt = t1 - t2
    Bx = b * S_ab * dt

    r3_minus_z = r3 - z
    r2_plus_bz = r2 + bz
    if fabs(r3_minus_z) > 0.0:
        r2_over_r3 = r2_plus_bz / r3_minus_z
        if fabs(r2_over_r3) > FP_CUTOFF:
            Bx += log(r2_plus_bz / r3_minus_z)

    Bx *= sigma

    return Bx


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64)], target="parallel"
)
def _charge_sheet_z(a, b, sigma, x, y, z):
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

    S_ab = 1.0 / sqrt(a ** 2 + b ** 2)

    r1 = sqrt(x ** 2 + y ** 2 + z ** 2)

    ax = a - x
    bz = b - z

    r2 = sqrt(ax ** 2 + y ** 2 + bz ** 2)
    r3 = sqrt(ax ** 2 + y ** 2 + z ** 2)

    t1_log = r1 - (a * x + b * z) * S_ab
    t2_log = r2 + (a * ax + b * bz) * S_ab

    if fabs(t1_log) > FP_CUTOFF:
        t1 = log(t1_log)
    else:
        t1 = 0.0

    if fabs(t2_log) > FP_CUTOFF:
        t2 = log(t2_log)
    else:
        t2 = 0.0

    Bz = a * S_ab * (t2 - t1)

    r3_plus_ax = r3 + ax
    r1_minus_x = r1 - x

    if fabs(r3_plus_ax) > 0.0:
        r1_over_r3 = r1_minus_x / r3_plus_ax
        if fabs(r1_over_r3) > FP_CUTOFF:
            Bz += log(r1_over_r3)
    Bz *= sigma

    return Bz


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64)], target="parallel"
)
def _charge_sheet_y(a, b, sigma, x, y, z):
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
    if fabs(y) > FP_CUTOFF:
        a_ab = sqrt(1 + b ** 2 / a ** 2)
        r1 = sqrt(x ** 2 + y ** 2 + z ** 2)

        ax = a - x

        r3 = sqrt(ax ** 2 + y ** 2 + z ** 2)

        g = b / a

        si_alpha = 1.0 / (a_ab ** 2)

        beta = -(x + z * g) * si_alpha
        gamma_sq = (r1 ** 2 * si_alpha) - (beta ** 2)
        if gamma_sq > 0.0:
            gamma = sqrt(gamma_sq)
        else:
            gamma = 0.0

        A = -gamma * g
        B = gamma * a_ab
        C = z + beta * g
        ABC_diff = B ** 2 - A ** 2 - C ** 2

        if fabs(gamma) > 0.0:
            t3 = (a + beta) / gamma
            t4 = beta / gamma
        else:
            t3 = 0.0
            t4 = 0.0

        if ABC_diff > 0.0:
            t1 = sqrt(ABC_diff)
            t2 = 1.0 / t1
        else:
            t1 = 0.0
            t2 = 0.0

        atan_1 = t2 * (C + (A + B) * (sqrt(1 + t3 * t3) + t3))
        atan_2 = t2 * (C + (A + B) * (sqrt(1 + t4 * t4) + t4))
        atan_y = atan_1 - atan_2
        atan_x = 1 + atan_1 * atan_2

        a_ab_t1 = a_ab * t1

        if fabs(atan_y) < FP_CUTOFF and fabs(atan_x) < FP_CUTOFF:
            By = 0.0
        elif a_ab_t1 > FP_CUTOFF:
            By = (y / (a_ab_t1)) * atan2(atan_y, atan_x)
        else:
            By = 0.0

        t2 = 1 / y
        atan_1 = t2 * (r3 + z + (x - a))
        atan_2 = t2 * (r1 + z + x)
        atan_y = atan_1 - atan_2
        atan_x = 1 + atan_1 * atan_2

        if fabs(atan_y) < FP_CUTOFF and fabs(atan_x) < FP_CUTOFF:
            By = 0.0
        else:
            By = By + atan2(atan_y, atan_x)

        By = sigma * By

    else:
        By = 0.0

    return By
