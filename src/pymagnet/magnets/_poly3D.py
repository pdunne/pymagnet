import numpy as _np
from stl import mesh
import plotly.graph_objects as _go


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

    if _np.allclose(rot_axis, null_vector, rtol=1e-4):

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

    if _np.allclose(rot_axis, null_vector, atol=1e-4):
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

    if _np.allclose(rot_axis, null_vector, atol=1e-4):
        # Check if anti-parallel
        if check_sign(z_axis, vec_z):
            # Parallel
            third_rotation = Quaternion()
            # tri_z = tri_x.copy()
        else:
            # Anti-parallel
            third_rotation = Quaternion.q_angle_from_axis(PI, y_axis)
            # tri_z = rot_points(tri_x.copy(), third_rotation)

    else:
        angle = -_np.arccos(_np.dot(z_axis, vec_z))
        third_rotation = Quaternion.q_angle_from_axis(angle, rot_axis)
        # tri_z = rot_points(tri_x.copy(), third_rotation)
    return second_rotation, third_rotation, new_vertex


def _rotate_triangle(triangle):
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

    origin_vertex = (longest_side + 1) % 3

    if norm_vec_y[1] < 0:
        offset = rotated_triangle[(origin_vertex - 1) % 3]
    else:
        offset = rotated_triangle[(origin_vertex + 0) % 3]

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

        self.mesh_scale = kwargs.pop("mesh_scale", 1e-3)
        self._filename = filename

        self.plot_data, self.mesh = self._import_mesh()

        self.Jx = _np.around(
            Jr * _np.cos(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        )
        self.Jy = _np.around(
            Jr * _np.sin(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        )
        self.Jz = _np.around(Jr * _np.cos(self.theta_rad), decimals=6)
        self.tol = 1e-4  # sufficient for 0.01 degree accuracy

        self.J = _np.array([self.Jx, self.Jy, self.Jz])

        self.Jnorm = _np.dot(self.J, self.mesh.normals.T)

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
        from ._routines3 import _allocate_field_array3

        B = _allocate_field_array3(x, y, z)
        vec_shape = B.x.shape
        B.x = B.x.ravel()
        B.y = B.y.ravel()
        B.z = B.z.ravel()

        for i in range(len(self.mesh.vectors)):
            triangle = self.mesh.vectors[i].copy()
            Btx, Bty, Btz, _, _, _ = self.calcB_triangle(
                triangle, self.Jnorm[i], x, y, z
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
        stl_mesh.vectors *= self.mesh_scale
        _, cog, _ = stl_mesh.get_mass_properties()
        print("COG:", cog)
        offset = self.center()
        print(offset)
        stl_mesh.translate(offset)
        _, cog, _ = stl_mesh.get_mass_properties()
        print("COG:", cog)

        mesh_rotation = Quaternion.gen_rotation_quaternion(
            self.alpha_rad, self.beta_rad, self.gamma_rad
        )
        angle, axis = mesh_rotation.get_axisangle()
        if _np.fabs(angle) > self.tol:
            stl_mesh.rotate(axis, angle)

        stl_mesh.normals = stl_mesh.normals / _np.linalg.norm(
            stl_mesh.normals, axis=1, keepdims=True
        )

        p, q, r = stl_mesh.vectors.shape  # (p, 3, 3)

        # the array stl_mesh.vectors.reshape(p*q, r) can contain multiple copies of the same vertex;
        # extract unique vertices from all mesh triangles
        vertices, ixr = _np.unique(
            stl_mesh.vectors.reshape(p * q, r), return_inverse=True, axis=0
        )
        I = _np.take(ixr, [3 * k for k in range(p)])
        J = _np.take(ixr, [3 * k + 1 for k in range(p)])
        K = _np.take(ixr, [3 * k + 2 for k in range(p)])
        x, y, z = vertices.T
        trace = _go.Mesh3d(x=x, y=y, z=z, i=I, j=J, k=K)

        # optional parameters to make it look nicer
        trace.update(
            flatshading=True, lighting_facenormalsepsilon=0, lighting_ambient=0.7
        )
        return trace, stl_mesh

    def _generate_mask(self, x, y, z):
        """Generates mask of points inside a magnet

        Args:
            x (ndarray/float): x-coordinates
            y (ndarray/float): y-coordinates
            z (ndarray/float): z-coordinates
        """
        pass

    def calcB_triangle(self, triangle, Jr, x, y, z):
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
        ) = _rotate_triangle(triangle.copy())

        # Prepare points and quaternion
        pos_vec = Quaternion._prepare_vector(x, y, z)

        # Rotate points
        x_rot, y_rot, z_rot = total_rotation * pos_vec
        if Jr < 0:
            RA_triangle1, RA_triangle2 = RA_triangle2, RA_triangle1

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

        return Bx, By, Bz