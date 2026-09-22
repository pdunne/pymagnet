from math import atan2, fabs, log, sqrt
from os import environ as _environ

import numpy as _np
from numba import float64, njit, prange, vectorize
from stl import mesh

from ..utils._quaternion import Quaternion, q_angle_from_axis
from ..utils._trigonometry3D import (
    _rotate_triangle,
    _rotate_triangle_njit,
    norm_plane,
    norm_plane_njit,
    rotate_vector_by_quat_inverse_njit,
    rotate_vector_by_quat_njit,
)
from ..utils.global_const import ALIGN_CUTOFF, FP_CUTOFF, MAG_TOL, PI
from ._magnet3D import Magnet3D

# Used to hide deprecated warning that doesn't affect the code
# Issue: https://github.com/numba/numba/issues/5520
_environ["KMP_WARNINGS"] = "0"


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

        self.phi = kwargs.pop("phi", 90.0)
        self.phi_rad = _np.deg2rad(self.phi)
        self.theta = kwargs.pop("theta", 0.0)
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

    def get_field(self, x, y, z, parallel=True, r_cut=_np.inf):
        """Calculates the magnetic field at point(s) x,y,z due to a 3D magnet
        The calculations are always performed in local coordinates with the centre of the magnet at origin and z magnetisation pointing along the local z' axis.

        The rotations and translations are performed first, and the internal field calculation functions are called.

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates
            parallel (bool): If True, use parallel numba implementation
            r_cut (float): Distance cutoff. Triangles whose centroid is farther
                than ``r_cut`` from an evaluation point are skipped.  Units must
                match the mesh coordinates (typically mm).  Default: no cutoff.

        Returns:
            tuple: Bx(ndarray), By(ndarray), Bz(ndarray)  field vector
        """
        if parallel:
            B = self._get_field_parallel(x, y, z, r_cut=r_cut)
        else:
            B = self._get_field_internal(x, y, z)

        return B.x, B.y, B.z

    def _get_field_parallel(self, x, y, z, r_cut=_np.inf):
        """Parallel magnetic field calculation using numba.

        Delegates to :meth:`_get_field_parallel_pts`, which parallelises over
        evaluation points (prange outer loop) rather than triangles.  This is
        the only correct parallel strategy: the original triangles-outer kernel
        had a race condition — multiple threads writing to the same output
        array indices without synchronisation — producing non-deterministic
        errors up to ~10% on large meshes.

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates
            r_cut (float): Distance cutoff passed to
                :meth:`_get_field_parallel_pts`.  Default: no cutoff.

        Returns:
            Field3: Magnetic field array
        """
        return self._get_field_parallel_pts(x, y, z, r_cut=r_cut)

    def _get_field_serial_fast(self, x, y, z):
        """Serial but numba-optimized magnetic field calculation.

        Uses the numba-compiled functions but processes triangles serially.
        Useful for comparison with parallel version.

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

        x_flat = _np.asarray(x).ravel().astype(_np.float64)
        y_flat = _np.asarray(y).ravel().astype(_np.float64)
        z_flat = _np.asarray(z).ravel().astype(_np.float64)

        mesh_vectors = _np.ascontiguousarray(self.mesh_vectors, dtype=_np.float64)
        Jnorm = _np.ascontiguousarray(self.Jnorm, dtype=_np.float64)

        Bx, By, Bz = _get_field_serial_njit(
            mesh_vectors, Jnorm, self.Jr, x_flat, y_flat, z_flat
        )

        Bx[~_np.isfinite(Bx)] = 0.0
        By[~_np.isfinite(By)] = 0.0
        Bz[~_np.isfinite(Bz)] = 0.0

        B.x = Bx.reshape(vec_shape)
        B.y = By.reshape(vec_shape)
        B.z = Bz.reshape(vec_shape)
        B.n = _np.linalg.norm([B.x, B.y, B.z], axis=0)

        return B

    def _get_field_parallel_pts(self, x, y, z, r_cut=_np.inf):
        """Transposed parallel field calculation: prange over evaluation points.

        Precomputes per-triangle rotation data once, then parallelises over
        the evaluation-point axis rather than the triangle axis.  This layout
        enables per-point triangle culling via an optional distance cutoff.

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates
            r_cut (float): distance cutoff in the same units as the mesh
                coordinates.  Triangles whose centroid is farther than r_cut
                from an evaluation point are skipped.  Default: np.inf
                (no culling — full accuracy).

        Returns:
            Field3: Magnetic field array
        """
        from ..utils._routines3D import _allocate_field_array3

        B = _allocate_field_array3(x, y, z)
        vec_shape = B.x.shape

        x_flat = _np.asarray(x).ravel().astype(_np.float64)
        y_flat = _np.asarray(y).ravel().astype(_np.float64)
        z_flat = _np.asarray(z).ravel().astype(_np.float64)

        mesh_vectors = _np.ascontiguousarray(self.mesh_vectors, dtype=_np.float64)
        Jnorm = _np.ascontiguousarray(self.Jnorm, dtype=_np.float64)

        rotations, offsets, RA_tris1, RA_tris2, swap_flags, active, centroids = (
            _precompute_triangle_data(mesh_vectors, Jnorm, self.Jr)
        )

        Bx, By, Bz = _get_field_parallel_pts_njit(
            rotations,
            offsets,
            RA_tris1,
            RA_tris2,
            swap_flags,
            active,
            centroids,
            Jnorm,
            x_flat,
            y_flat,
            z_flat,
            float(r_cut),
        )

        Bx[~_np.isfinite(Bx)] = 0.0
        By[~_np.isfinite(By)] = 0.0
        Bz[~_np.isfinite(Bz)] = 0.0

        B.x = Bx.reshape(vec_shape)
        B.y = By.reshape(vec_shape)
        B.z = Bz.reshape(vec_shape)
        B.n = _np.linalg.norm([B.x, B.y, B.z], axis=0)

        return B

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
        volume *= self.mesh_scale**3
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
        rotate_about_z = q_angle_from_axis(PI, (0, 0, 1))
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
    S_ab = 1.0 / sqrt(a**2 + b**2)

    r1 = sqrt(x**2 + y**2 + z**2)

    ax = a - x
    bz = b - z

    r2 = sqrt(ax**2 + y**2 + bz**2)
    r3 = sqrt(ax**2 + y**2 + z**2)
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

    S_ab = 1.0 / sqrt(a**2 + b**2)

    r1 = sqrt(x**2 + y**2 + z**2)

    ax = a - x
    bz = b - z

    r2 = sqrt(ax**2 + y**2 + bz**2)
    r3 = sqrt(ax**2 + y**2 + z**2)

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
        a_ab = sqrt(1 + b**2 / a**2)
        r1 = sqrt(x**2 + y**2 + z**2)

        ax = a - x

        r3 = sqrt(ax**2 + y**2 + z**2)

        g = b / a

        si_alpha = 1.0 / (a_ab**2)

        beta = -(x + z * g) * si_alpha
        gamma_sq = (r1**2 * si_alpha) - (beta**2)
        if gamma_sq > 0.0:
            gamma = sqrt(gamma_sq)
        else:
            gamma = 0.0

        A = -gamma * g
        B = gamma * a_ab
        C = z + beta * g
        ABC_diff = B**2 - A**2 - C**2

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


# =============================================================================
# Numba-compatible charge sheet functions for parallel processing
# =============================================================================


@njit(cache=True)
def _charge_sheet_x_scalar(a, b, sigma, x, y, z):
    """Scalar version of charge sheet x-component (numba-compatible)."""
    S_ab = 1.0 / _np.sqrt(a * a + b * b)

    r1 = _np.sqrt(x * x + y * y + z * z)

    ax = a - x
    bz = b - z

    r2 = _np.sqrt(ax * ax + y * y + bz * bz)
    r3 = _np.sqrt(ax * ax + y * y + z * z)
    t1_log = r1 - (a * x + b * z) * S_ab
    t2_log = r2 + (a * ax + b * bz) * S_ab

    if _np.fabs(t1_log) > FP_CUTOFF:
        t1 = _np.log(t1_log)
    else:
        t1 = 0.0

    if _np.fabs(t2_log) > FP_CUTOFF:
        t2 = _np.log(t2_log)
    else:
        t2 = 0.0

    dt = t1 - t2
    Bx = b * S_ab * dt

    r3_minus_z = r3 - z
    r2_plus_bz = r2 + bz

    if _np.fabs(r3_minus_z) > 0.0:
        r2_over_r3 = r2_plus_bz / r3_minus_z
        if _np.fabs(r2_over_r3) > FP_CUTOFF:
            Bx += _np.log(r2_plus_bz / r3_minus_z)

    Bx *= sigma / PI / 4.0

    return Bx


@njit(cache=True)
def _charge_sheet_y_scalar(a, b, sigma, x, y, z):
    """Scalar version of charge sheet y-component (numba-compatible)."""
    if _np.fabs(y) > FP_CUTOFF:
        a_ab = _np.sqrt(1.0 + b * b / (a * a))
        r1 = _np.sqrt(x * x + y * y + z * z)

        ax = a - x
        r3 = _np.sqrt(ax * ax + y * y + z * z)

        g = b / a
        si_alpha = 1.0 / (a_ab * a_ab)

        beta = -(x + z * g) * si_alpha
        gamma_sq = (r1 * r1 * si_alpha) - (beta * beta)
        if gamma_sq > 0.0:
            gamma = _np.sqrt(gamma_sq)
        else:
            gamma = 0.0

        A = -gamma * g
        B = gamma * a_ab
        C = z + beta * g
        ABC_diff = B * B - A * A - C * C

        if _np.fabs(gamma) > 0.0:
            t3 = (a + beta) / gamma
            t4 = beta / gamma
        else:
            t3 = 0.0
            t4 = 0.0

        if ABC_diff > 0.0:
            t1 = _np.sqrt(ABC_diff)
            t2 = 1.0 / t1
        else:
            t1 = 0.0
            t2 = 0.0

        atan_1 = t2 * (C + (A + B) * (_np.sqrt(1.0 + t3 * t3) + t3))
        atan_2 = t2 * (C + (A + B) * (_np.sqrt(1.0 + t4 * t4) + t4))
        atan_y = atan_1 - atan_2
        atan_x = 1.0 + atan_1 * atan_2

        a_ab_t1 = a_ab * t1

        if _np.fabs(atan_y) < FP_CUTOFF and _np.fabs(atan_x) < FP_CUTOFF:
            By = 0.0
        elif a_ab_t1 > FP_CUTOFF:
            By = (y / a_ab_t1) * _np.arctan2(atan_y, atan_x)
        else:
            By = 0.0

        t2 = 1.0 / y
        atan_1 = t2 * (r3 + z + (x - a))
        atan_2 = t2 * (r1 + z + x)
        atan_y = atan_1 - atan_2
        atan_x = 1.0 + atan_1 * atan_2

        if _np.fabs(atan_y) < FP_CUTOFF and _np.fabs(atan_x) < FP_CUTOFF:
            By = 0.0
        else:
            By = By + _np.arctan2(atan_y, atan_x)

        By = sigma * By / PI / 2.0
    else:
        By = 0.0

    return By


@njit(cache=True)
def _charge_sheet_z_scalar(a, b, sigma, x, y, z):
    """Scalar version of charge sheet z-component (numba-compatible)."""
    S_ab = 1.0 / _np.sqrt(a * a + b * b)

    r1 = _np.sqrt(x * x + y * y + z * z)

    ax = a - x
    bz = b - z

    r2 = _np.sqrt(ax * ax + y * y + bz * bz)
    r3 = _np.sqrt(ax * ax + y * y + z * z)

    t1_log = r1 - (a * x + b * z) * S_ab
    t2_log = r2 + (a * ax + b * bz) * S_ab

    if _np.fabs(t1_log) > FP_CUTOFF:
        t1 = _np.log(t1_log)
    else:
        t1 = 0.0

    if _np.fabs(t2_log) > FP_CUTOFF:
        t2 = _np.log(t2_log)
    else:
        t2 = 0.0

    Bz = a * S_ab * (t2 - t1)

    r3_plus_ax = r3 + ax
    r1_minus_x = r1 - x

    if _np.fabs(r3_plus_ax) > 0.0:
        r1_over_r3 = r1_minus_x / r3_plus_ax
        if _np.fabs(r1_over_r3) > FP_CUTOFF:
            Bz += _np.log(r1_over_r3)

    Bz *= sigma / PI / 4.0

    return Bz


@njit(cache=True)
def _charge_sheet_njit(a, b, sigma, x, y, z):
    """Compute all three components of charge sheet field for arrays.

    Args:
        a (float): triangle base
        b (float): triangle altitude
        sigma (float): charge density (Jr)
        x, y, z (ndarray): coordinate arrays (1D)

    Returns:
        tuple: (Bx, By, Bz) arrays
    """
    n = x.size
    Bx = _np.empty(n)
    By = _np.empty(n)
    Bz = _np.empty(n)

    for i in range(n):
        Bx[i] = _charge_sheet_x_scalar(a, b, sigma, x[i], y[i], z[i])
        By[i] = _charge_sheet_y_scalar(a, b, sigma, x[i], y[i], z[i])
        Bz[i] = _charge_sheet_z_scalar(a, b, sigma, x[i], y[i], z[i])

    return Bx, By, Bz


# Pre-computed rotation quaternion for 180° about z-axis
_ROTATE_Z_PI = _np.array([0.0, 0.0, 0.0, 1.0])  # [cos(π/2), 0, 0, sin(π/2)]


@njit(cache=True)
def _calcB_2_triangles_njit(triangle1, triangle2, Jr, x, y, z):
    """Calculate field from two right-angled triangles (numba-compatible).

    Args:
        triangle1 (ndarray): (2,) [base, altitude] of first RA triangle
        triangle2 (ndarray): (2,) [base, altitude] of second RA triangle
        Jr (float): normal magnetization component
        x, y, z (ndarray): coordinate arrays

    Returns:
        tuple: (Bx, By, Bz) field arrays
    """
    # Calculate field from first right-angled triangle
    Btx, Bty, Btz = _charge_sheet_njit(triangle1[0], triangle1[1], Jr, x, y, z)

    # Rotate coordinates into local frame of second triangle (180° about z)
    # For 180° rotation about z: x' = -x, y' = -y, z' = z
    x_offset = x - triangle1[0]
    x_local = -x_offset
    y_local = -y
    z_local = z.copy()

    # Calculate field from second right-angled triangle
    Btx2, Bty2, Btz2 = _charge_sheet_njit(
        triangle2[0], triangle2[1], Jr, x_local + triangle2[0], y_local, z_local
    )

    # Inverse rotation of field (180° about z): Bx' = -Bx, By' = -By, Bz' = Bz
    Btx2_rot = -Btx2
    Bty2_rot = -Bty2
    Btz2_rot = Btz2

    # Sum contributions
    Bx_out = Btx + Btx2_rot
    By_out = Bty + Bty2_rot
    Bz_out = Btz + Btz2_rot

    return Bx_out, By_out, Bz_out


@njit(cache=True)
def _calcB_triangle_njit(triangle, Jr, x_flat, y_flat, z_flat):
    """Calculate magnetic field from a single triangle (numba-compatible).

    Args:
        triangle (ndarray): (3, 3) triangle vertices
        Jr (float): normal magnetization component
        x_flat, y_flat, z_flat (ndarray): flattened coordinate arrays

    Returns:
        tuple: (Bx, By, Bz) field arrays
    """
    # Rotate triangle to standard orientation
    total_rotation, _rotated_triangle, offset, RA_triangle1, RA_triangle2 = (
        _rotate_triangle_njit(triangle)
    )

    # Rotate coordinate points
    x_rot, y_rot, z_rot = rotate_vector_by_quat_njit(
        total_rotation, x_flat, y_flat, z_flat
    )

    # Get triangle normal
    norm1 = norm_plane_njit(triangle)

    # Check if we need to swap triangles (anti-parallel to -y with negative Jr)
    if (
        _np.fabs(norm1[0]) < ALIGN_CUTOFF
        and norm1[1] < -0.99
        and _np.fabs(norm1[2]) < ALIGN_CUTOFF
        and Jr < 0
    ):
        RA_triangle1, RA_triangle2 = RA_triangle2, RA_triangle1

    # Calculate field in rotated frame
    Btx, Bty, Btz = _calcB_2_triangles_njit(
        RA_triangle1,
        RA_triangle2,
        Jr,
        x_rot - offset[0],
        y_rot - offset[1],
        z_rot - offset[2],
    )

    # Rotate field back to global frame (inverse rotation)
    Bx, By, Bz = rotate_vector_by_quat_inverse_njit(total_rotation, Btx, Bty, Btz)

    return Bx, By, Bz


@njit(cache=True)
def _get_field_serial_njit(
    mesh_vectors, Jnorm, Jr, x_flat, y_flat, z_flat, threshold=1e-4
):
    """Calculate magnetic field from all triangles serially (for comparison).

    Args:
        mesh_vectors (ndarray): (N, 3, 3) array of triangle vertices
        Jnorm (ndarray): (N,) normal magnetization for each triangle
        Jr (float): remnant magnetization magnitude
        x_flat, y_flat, z_flat (ndarray): flattened coordinate arrays
        threshold (float): minimum |Jnorm/Jr| to include triangle

    Returns:
        tuple: (Bx, By, Bz) total field arrays
    """
    n_triangles = mesh_vectors.shape[0]
    n_points = x_flat.size

    Bx_total = _np.zeros(n_points)
    By_total = _np.zeros(n_points)
    Bz_total = _np.zeros(n_points)

    for i in range(n_triangles):
        if _np.fabs(Jnorm[i] / Jr) > threshold:
            Btx, Bty, Btz = _calcB_triangle_njit(
                mesh_vectors[i], Jnorm[i], x_flat, y_flat, z_flat
            )

            Bx_total += Btx
            By_total += Bty
            Bz_total += Btz

    return Bx_total, By_total, Bz_total


# ---------------------------------------------------------------------------
# Transposed (points-outer) parallel kernel
# ---------------------------------------------------------------------------


def _precompute_triangle_data(mesh_vectors, Jnorm, Jr, threshold=1e-4):
    """Precompute per-triangle rotation data for the points-outer kernel.

    Calls _rotate_triangle_njit once per active triangle and stores the
    results in flat contiguous arrays so they can be passed to a @njit kernel.

    Args:
        mesh_vectors (ndarray): (N, 3, 3) triangle vertices
        Jnorm (ndarray): (N,) normal magnetisation per triangle
        Jr (float): remnant magnetisation magnitude
        threshold (float): minimum |Jnorm/Jr| to include a triangle

    Returns:
        rotations  (N, 4)  float64 — rotation quaternion per triangle
        offsets    (N, 3)  float64 — translation offset after rotation
        RA_tris1   (N, 2)  float64 — [base, altitude] of first RA sub-triangle
        RA_tris2   (N, 2)  float64 — [base, altitude] of second RA sub-triangle
        swap_flags (N,)    bool    — True when RA sub-triangles must be swapped
        active     (N,)    bool    — |Jnorm/Jr| > threshold
        centroids  (N, 3)  float64 — triangle centroids for distance culling
    """
    n_tri = mesh_vectors.shape[0]
    rotations = _np.zeros((n_tri, 4), dtype=_np.float64)
    offsets = _np.zeros((n_tri, 3), dtype=_np.float64)
    RA_tris1 = _np.zeros((n_tri, 2), dtype=_np.float64)
    RA_tris2 = _np.zeros((n_tri, 2), dtype=_np.float64)
    swap_flags = _np.zeros(n_tri, dtype=_np.bool_)
    active = _np.abs(Jnorm / Jr) > threshold
    centroids = mesh_vectors.mean(axis=1).astype(_np.float64)

    for i in range(n_tri):
        if active[i]:
            q, _, offset, RA1, RA2 = _rotate_triangle_njit(mesh_vectors[i])
            rotations[i] = q
            offsets[i] = offset
            RA_tris1[i] = RA1
            RA_tris2[i] = RA2
            norm1 = norm_plane_njit(mesh_vectors[i])
            swap_flags[i] = (
                _np.fabs(norm1[0]) < ALIGN_CUTOFF
                and norm1[1] < -0.99
                and _np.fabs(norm1[2]) < ALIGN_CUTOFF
                and Jnorm[i] < 0
            )

    return (
        _np.ascontiguousarray(rotations),
        _np.ascontiguousarray(offsets),
        _np.ascontiguousarray(RA_tris1),
        _np.ascontiguousarray(RA_tris2),
        _np.ascontiguousarray(swap_flags),
        _np.ascontiguousarray(active),
        _np.ascontiguousarray(centroids),
    )


@njit(parallel=True, cache=True)
def _get_field_parallel_pts_njit(
    rotations,
    offsets,
    RA_tris1,
    RA_tris2,
    swap_flags,
    active,
    centroids,
    Jnorm,
    x_flat,
    y_flat,
    z_flat,
    r_cut,
):
    """Transposed parallel field kernel: prange over evaluation points.

    For each evaluation point, iterates over all active triangles and
    accumulates field contributions using scalar arithmetic (no per-triangle
    array allocations).  Quaternion rotations are inlined.

    Args:
        rotations  (N_tri, 4): rotation quaternion per triangle
        offsets    (N_tri, 3): translation offset per triangle
        RA_tris1   (N_tri, 2): first RA sub-triangle [base, altitude]
        RA_tris2   (N_tri, 2): second RA sub-triangle [base, altitude]
        swap_flags (N_tri,):   swap RA sub-triangles when True
        active     (N_tri,):   include triangle when True
        centroids  (N_tri, 3): triangle centroids for distance culling
        Jnorm      (N_tri,):   normal magnetisation per triangle
        x_flat, y_flat, z_flat: flattened evaluation-point coordinates
        r_cut (float): distance cutoff — triangles whose centroid is farther
            than r_cut from the evaluation point are skipped.
            Pass np.inf to disable culling.

    Returns:
        tuple: (Bx, By, Bz) total field arrays
    """
    n_tri = rotations.shape[0]
    n_pts = x_flat.size
    r_cut_sq = r_cut * r_cut

    Bx_total = _np.zeros(n_pts)
    By_total = _np.zeros(n_pts)
    Bz_total = _np.zeros(n_pts)

    for j in prange(n_pts):  # ty:ignore[not-iterable]
        bx = 0.0
        by = 0.0
        bz = 0.0
        xj = x_flat[j]
        yj = y_flat[j]
        zj = z_flat[j]

        for i in range(n_tri):
            if not active[i]:
                continue

            # ---- distance culling ----
            ddx = centroids[i, 0] - xj
            ddy = centroids[i, 1] - yj
            ddz = centroids[i, 2] - zj
            if ddx * ddx + ddy * ddy + ddz * ddz > r_cut_sq:
                continue

            # ---- inline quaternion rotation of evaluation point ----
            qw = rotations[i, 0]
            qx = rotations[i, 1]
            qy = rotations[i, 2]
            qz = rotations[i, 3]

            tx = 2.0 * (qy * zj - qz * yj)
            ty = 2.0 * (qz * xj - qx * zj)
            tz = 2.0 * (qx * yj - qy * xj)
            xr = xj + qw * tx + (qy * tz - qz * ty)
            yr = yj + qw * ty + (qz * tx - qx * tz)
            zr = zj + qw * tz + (qx * ty - qy * tx)

            xr -= offsets[i, 0]
            yr -= offsets[i, 1]
            zr -= offsets[i, 2]

            # ---- RA sub-triangle selection ----
            if swap_flags[i]:
                a1 = RA_tris2[i, 0]
                b1 = RA_tris2[i, 1]
                a2 = RA_tris1[i, 0]
                b2 = RA_tris1[i, 1]
            else:
                a1 = RA_tris1[i, 0]
                b1 = RA_tris1[i, 1]
                a2 = RA_tris2[i, 0]
                b2 = RA_tris2[i, 1]

            jr = Jnorm[i]

            # ---- first RA triangle ----
            Bx1 = _charge_sheet_x_scalar(a1, b1, jr, xr, yr, zr)
            By1 = _charge_sheet_y_scalar(a1, b1, jr, xr, yr, zr)
            Bz1 = _charge_sheet_z_scalar(a1, b1, jr, xr, yr, zr)

            # ---- second RA triangle (180° z-rotation of coords) ----
            xl = -(xr - a1)
            yl = -yr
            # zl == zr (unchanged)
            Bx2 = -_charge_sheet_x_scalar(a2, b2, jr, xl + a2, yl, zr)
            By2 = -_charge_sheet_y_scalar(a2, b2, jr, xl + a2, yl, zr)
            Bz2 = _charge_sheet_z_scalar(a2, b2, jr, xl + a2, yl, zr)

            # Sum in rotated frame
            Bxr = Bx1 + Bx2
            Byr = By1 + By2
            Bzr = Bz1 + Bz2

            # ---- inline inverse quaternion rotation of field ----
            # Conjugate: negate xyz components
            tx = 2.0 * (-qy * Bzr - (-qz) * Byr)
            ty = 2.0 * (-qz * Bxr - (-qx) * Bzr)
            tz = 2.0 * (-qx * Byr - (-qy) * Bxr)
            bx += Bxr + qw * tx + (-qy * tz - (-qz) * ty)
            by += Byr + qw * ty + (-qz * tx - (-qx) * tz)
            bz += Bzr + qw * tz + (-qx * ty - (-qy) * tx)

        Bx_total[j] = bx
        By_total[j] = by
        Bz_total[j] = bz

    return Bx_total, By_total, Bz_total


# ---------------------------------------------------------------------------
# Multi-magnet fused path
# ---------------------------------------------------------------------------


def _precompute_all_meshes(meshes, threshold=1e-4):
    """Precompute and concatenate triangle data for a list of Mesh magnets.

    Calls ``_precompute_triangle_data`` once per magnet and concatenates the
    results so that all triangles from all magnets can be processed in a single
    ``_get_field_parallel_pts_njit`` call.

    Args:
        meshes: iterable of ``Mesh`` instances
        threshold (float): minimum |Jnorm/Jr| to include a triangle

    Returns:
        tuple: (rotations, offsets, RA_tris1, RA_tris2, swap_flags, active,
                centroids, Jnorm) — all contiguous float64/bool arrays with
                rows from every magnet concatenated in order.
    """
    mesh_list = list(meshes)
    parts = [
        _precompute_triangle_data(
            _np.ascontiguousarray(m.mesh_vectors, dtype=_np.float64),
            _np.ascontiguousarray(m.Jnorm, dtype=_np.float64),
            m.Jr,
            threshold,
        )
        for m in mesh_list
    ]
    Jnorm_parts = [_np.ascontiguousarray(m.Jnorm, dtype=_np.float64) for m in mesh_list]
    combined = tuple(_np.concatenate([p[k] for p in parts]) for k in range(7))
    return (*combined, _np.concatenate(Jnorm_parts))


def get_total_field_mesh(meshes, x, y, z, r_cut=_np.inf):
    """Compute the total magnetic field from multiple Mesh magnets in one pass.

    Concatenates triangle data from all meshes and evaluates the field at
    every point ``(x, y, z)`` using the points-outer parallel kernel
    ``_get_field_parallel_pts_njit``.  An optional distance cutoff skips
    triangles whose centroid is farther than ``r_cut`` from an evaluation
    point, which can give a large speedup for sparse or localised geometries.

    Args:
        meshes (list): list (or any iterable) of ``Mesh`` instances.  Pass
            ``pm.magnets.Mesh.instances`` to include all currently registered
            meshes.
        x (float or ndarray): x co-ordinates of the evaluation points.
        y (float or ndarray): y co-ordinates of the evaluation points.
        z (float or ndarray): z co-ordinates of the evaluation points.
        r_cut (float): distance cutoff in the same length units as the mesh
            coordinates.  Triangles farther than ``r_cut`` from a point are
            skipped.  Default: ``np.inf`` (no culling — full accuracy).

    Returns:
        Field3: total magnetic field array (attributes ``.x``, ``.y``, ``.z``,
            ``.n``).

    Example::

        import pymagnet as pm
        import numpy as np

        pm.reset()
        m1 = pm.magnets.Mesh("left.stl",  Jr=1.0, center=[-30, 0, 0])
        m2 = pm.magnets.Mesh("right.stl", Jr=1.0, center=[ 30, 0, 0])

        x = np.linspace(-60, 60, 40)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

        # Single fused pass — equivalent to summing m1.get_field() + m2.get_field()
        B = pm.magnets.get_total_field_mesh([m1, m2], X, Y, Z, r_cut=40.0)
    """
    from ..utils._routines3D import _allocate_field_array3

    B = _allocate_field_array3(x, y, z)
    vec_shape = B.x.shape

    x_flat = _np.asarray(x).ravel().astype(_np.float64)
    y_flat = _np.asarray(y).ravel().astype(_np.float64)
    z_flat = _np.asarray(z).ravel().astype(_np.float64)

    rotations, offsets, RA_tris1, RA_tris2, swap_flags, active, centroids, Jnorm = (
        _precompute_all_meshes(meshes)
    )

    Bx, By, Bz = _get_field_parallel_pts_njit(
        rotations,
        offsets,
        RA_tris1,
        RA_tris2,
        swap_flags,
        active,
        centroids,
        Jnorm,
        x_flat,
        y_flat,
        z_flat,
        float(r_cut),
    )

    Bx[~_np.isfinite(Bx)] = 0.0
    By[~_np.isfinite(By)] = 0.0
    Bz[~_np.isfinite(Bz)] = 0.0

    B.x = Bx.reshape(vec_shape)
    B.y = By.reshape(vec_shape)
    B.z = Bz.reshape(vec_shape)
    B.n = _np.linalg.norm([B.x, B.y, B.z], axis=0)
    return B
