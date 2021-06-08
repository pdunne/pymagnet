from ._magnet3D import Magnet3D
from ..utils._quaternion import Quaternion
from ..utils.global_const import MAG_TOL, PI, FP_CUTOFF, ALIGN_CUTOFF
from ..utils._trigonometry3D import norm_plane, _rotate_triangle

from stl import mesh
import numpy as _np
from numba import jit, vectorize, float64
from math import sqrt, log, fabs, atan2


class Mesh(Magnet3D):
    """Mesh Magnet Class."""

    mag_type = "Mesh"

    def __init__(
        self,
        filename,
        Jr=1.0,  # local magnetisation
        **kwargs,
    ):
        """Init Method

        Args:
            filename (string): path to stl file to be imported
            Jr (float, optional): Signed remnant magnetisation. Defaults to 1.0.

        Kwargs:
            phi (float):
            theta (float):
            mesh_scale (float): scaling factor if mesh needs to be resized. Defaults to 1.0
        """
        super().__init__(Jr, **kwargs)

        self.phi = kwargs.pop("theta", 90.0)
        self.phi_rad = _np.deg2rad(self.phi)
        self.theta = kwargs.pop("phi", 0.0)
        self.theta_rad = _np.deg2rad(self.theta)

        self.mesh_scale = kwargs.pop("mesh_scale", 1.0)
        self._filename = filename

        (
            self.mesh_vectors,
            self.mesh_normals,
            self.volume,
            self.centroid,
        ) = self._import_mesh()

        self.Jx = _np.around(
            Jr * _np.cos(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        )
        self.Jy = _np.around(
            Jr * _np.sin(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        )
        self.Jz = _np.around(Jr * _np.cos(self.theta_rad), decimals=6)
        self.tol = MAG_TOL  # sufficient for 0.01 degree accuracy
        self.J = _np.array([self.Jx, self.Jy, self.Jz])

        # FIXME: Sort out rotation of magnetisation with rotation of mesh
        # if _np.any(
        #     _np.fabs([self.alpha_rad, self.beta_rad, self.gamma_rad]) > ALIGN_CUTOFF
        # ):
        #     mag_rotation = Quaternion.gen_rotation_quaternion(
        #         self.alpha_rad, self.beta_rad, self.gamma_rad
        #     )
        #     Jrot = mag_rotation * self.J
        #     self.Jx = Jrot[0]
        #     self.Jy = Jrot[1]
        #     self.Jz = Jrot[2]
        #     self.J = _np.array([self.Jx, self.Jy, self.Jz])

        self.Jnorm = _np.dot(self.J, self.mesh_normals.T)

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Center {self.get_center()}\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Center {self.get_center()}\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def get_Jr(self):
        """Returns Magnetisation vector

        Returns:
            ndarray: [Jx, Jy, Jz]
        """
        return self.J

    def size(self):
        """Returns magnet dimesions

        Returns:
            size (ndarray): numpy array [width, depth, height]
        """
        pass

    def get_center(self):
        """Returns magnet center

        Returns:
            ndarray: magnet center
        """
        return self.center

    def get_field(self, x, y, z):
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
        B = self._get_field_internal(x, y, z)

        return B.x, B.y, B.z

    def get_force_torque(self, depth=4, unit="mm"):
        """Calculates the force and torque on a prism magnet due to all other magnets.

        Args:
            depth (int, optional): Number of recursions of division by 4 per simplex
            unit (str, optional): Length scale. Defaults to 'mm'.

        Returns:
            tuple: force (ndarray (3,) ) and torque (ndarray (3,) )
        """
        from ..forces._mesh_force import calc_force_mesh

        force, torque = calc_force_mesh(self, depth, unit)
        return force, torque

    def _get_field_internal(self, x, y, z):
        """Internal magnetic field calculation methods.
        Iterates over each triangle that makes up the mesh magnet and calculates the magnetic field

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates

        Returns:
            Field3: Magnetic field array
        """
        from ..utils._routines3D import _allocate_field_array3

        B = _allocate_field_array3(x, y, z)
        vec_shape = B.x.shape
        B.x = B.x.ravel()
        B.y = B.y.ravel()
        B.z = B.z.ravel()

        # debug for loop, used when needing to check certain triangles, or groups of triangles
        # for i in range(self.start, self.stop):
        for i in range(len(self.mesh_vectors)):
            if _np.fabs(self.Jnorm[i] / self.Jr) > 1e-4:
                Btx, Bty, Btz, _, _, _ = self.calcB_triangle(
                    self.mesh_vectors[i],
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

    def _import_mesh(self):
        """Imports mesh from STL file

        Returns:
            tuple: mesh_vectors (ndarray of mesh triangles), mesh_normals (ndarray of normals to each triangle)
        """
        stl_mesh = mesh.Mesh.from_file(self._filename)

        if _np.any(
            _np.fabs([self.alpha_rad, self.beta_rad, self.gamma_rad]) > ALIGN_CUTOFF
        ):
            mesh_rotation = Quaternion.gen_rotation_quaternion(
                self.alpha_rad, self.beta_rad, self.gamma_rad
            )

            angle, axis = mesh_rotation.get_axisangle()
            stl_mesh.rotate(axis, angle)

        # to ensure that the initial center is set to the centroid
        _, centroid, _ = stl_mesh.get_mass_properties()
        stl_mesh.translate(-centroid)

        offset = self.get_center()
        stl_mesh.translate(offset / self.mesh_scale)

        # get values after translation
        volume, centroid, _ = stl_mesh.get_mass_properties()

        mesh_vectors = stl_mesh.vectors.astype(_np.float64)
        mesh_normals = stl_mesh.normals.astype(_np.float64)

        # scale values
        volume *= self.mesh_scale ** 3
        centroid *= self.mesh_scale
        mesh_vectors *= self.mesh_scale

        mesh_normals = mesh_normals / _np.linalg.norm(
            mesh_normals, axis=1, keepdims=True
        )

        return mesh_vectors, mesh_normals, volume, centroid

    def _generate_mask(self, x, y, z):
        """Generates mask of points inside a magnet
        NOTE: not implemented for Mesh magnets.
        Args:
            x (ndarray/float): x-coordinates
            y (ndarray/float): y-coordinates
            z (ndarray/float): z-coordinates
        """
        pass

    def calcB_triangle(self, triangle, Jr, x, y, z, i):
        """Calculates the magnetic field due to a triangle

        Args:
            triangle (ndarray): Vertices of a triangle
            Jr (float): Remnant magnetisation component normal to triangle
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

        norm1 = norm_plane(triangle)

        if _np.allclose(norm1, [0, -1, 0], atol=ALIGN_CUTOFF) and Jr < 0:
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
        sigma = Jr
        with _np.errstate(all="ignore"):
            Bx = _charge_sheet_x(a, b, sigma, x, y, z)
            By = _charge_sheet_y(a, b, sigma, x, y, z)
            Bz = _charge_sheet_z(a, b, sigma, x, y, z)
        return Bx, By, Bz


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64)], target="parallel"
)
def _charge_sheet_x(a, b, sigma, x, y, z):
    """Calculates the x-component of the magnetic field of a right angled charge sheet

    Args:
        a (float): triangle base
        b (float): triangle altitude
        sigma (float): normal magnetic charge density in tesla
        x (ndarray): x-coordinates
        y (ndarray): y-coordinates
        z (ndarray): z-coordinates

    Returns:
        ndarray: Bx magnetic field component
    """
    S_ab = 1.0 / sqrt(a ** 2 + b ** 2)

    r1 = sqrt(x ** 2 + y ** 2 + z ** 2)

    ax = a - x
    bz = b - z

    r2 = sqrt(ax ** 2 + y ** 2 + bz ** 2)
    r3 = sqrt(ax ** 2 + y ** 2 + z ** 2)
    t1_log = r1 - (a * x + b * z) * S_ab

    t2_log = r2 + (a * ax + b * bz) * S_ab

    # To avoid introduction of noise when taking log of values
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

    # if statements are used to avoid singularities such as divide by zero and miminise noise
    if fabs(r3_minus_z) > 0.0:
        r2_over_r3 = r2_plus_bz / r3_minus_z
        if fabs(r2_over_r3) > FP_CUTOFF:
            Bx += log(r2_plus_bz / r3_minus_z)

    Bx *= sigma / PI / 4

    return Bx


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64)], target="parallel"
)
def _charge_sheet_z(a, b, sigma, x, y, z):
    """Calculates the z-component of the magnetic field of a right angled charge sheet

    Args:
    a (float): triangle base
    b (float): triangle altitude
    sigma (float): normal magnetic charge density in tesla
    x (ndarray): x-coordinates
    y (ndarray): y-coordinates
    z (ndarray): z-coordinates

    Returns:
    ndarray: Bz magnetic field component
    """

    S_ab = 1.0 / sqrt(a ** 2 + b ** 2)

    r1 = sqrt(x ** 2 + y ** 2 + z ** 2)

    ax = a - x
    bz = b - z

    r2 = sqrt(ax ** 2 + y ** 2 + bz ** 2)
    r3 = sqrt(ax ** 2 + y ** 2 + z ** 2)

    t1_log = r1 - (a * x + b * z) * S_ab
    t2_log = r2 + (a * ax + b * bz) * S_ab

    # To avoid introduction of noise when taking log of values
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

    # if statements are used to avoid singularities such as divide by zero and miminise noise
    if fabs(r3_plus_ax) > 0.0:
        r1_over_r3 = r1_minus_x / r3_plus_ax
        if fabs(r1_over_r3) > FP_CUTOFF:
            Bz += log(r1_over_r3)
    Bz *= sigma / PI / 4

    return Bz


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64)], target="parallel"
)
def _charge_sheet_y(a, b, sigma, x, y, z):
    """Calculates the y-component of the magnetic field of a right angled charge sheet

    Args:
    a (float): triangle base
    b (float): triangle altitude
    sigma (float): normal magnetic charge density in tesla
    x (ndarray): x-coordinates
    y (ndarray): y-coordinates
    z (ndarray): z-coordinates

    Returns:
    ndarray: By magnetic field component
    """

    # if statements are used to avoid singularities such as divide by zero
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

        # To avoid undefined atan2 behaviour
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

        By = sigma * By / PI / 2

    else:
        By = 0.0

    return By
