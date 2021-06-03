__all__ = ["calc_force_prism"]

import numpy as _np
from ..utils.global_const import FP_CUTOFF, MU0
from ..utils._conversions import get_unit_value_meter
from ..utils._vector_structs import Point_Array3
from ..utils._routines3D import _allocate_field_array3
from ..utils._quaternion import Quaternion
from ..magnets import Magnet3D


def _get_ranges_prism(active_magnet):
    """Gets min and max x,y,z values for magnet faces.

    This is needed to generate a grid of points on each surface for calculating
    the force and torque on the magnet due to the other magnets.

    Args:
        active_magnet (Prism): Target magnet

    Returns:
        tuple: xlim (tuple), ylim (tuple), zlim (tuple), areas (ndarray) - areas of each face of the cuboid
    """
    size_x, size_y, size_z = active_magnet.get_size()

    xlim = (-size_x / 2, size_x / 2)
    ylim = (-size_y / 2, size_y / 2)
    zlim = (-size_z / 2, size_z / 2)

    areas = _np.array([size_y * size_z, size_x * size_z, size_x * size_y])
    return xlim, ylim, zlim, areas


def _gen_grid(xlim=(0, 10), ylim=(0, 10), num_samples=10):
    """Generates a 2D grid of N**2 points for a given pair of x and y limits.

    Args:
        xlim (tuple, optional): Min and max x values. Defaults to (0, 10).
        ylim (tuple, optional): Min and max y values. Defaults to (0, 10).
        num_samples (int, optional): Number of points to generate in each direction. Defaults to 10.

    Returns:
        ndarray: (num_samples**2, 3) array of points
    """
    xmin = xlim[0]
    xmax = xlim[1]
    ymin = ylim[0]
    ymax = ylim[1]

    dx = xmax - xmin
    dy = ymax - ymin

    xmin += dx / 2 / num_samples
    xmax -= dx / 2 / num_samples

    ymin += dy / 2 / num_samples
    ymax -= dy / 2 / num_samples

    NJ = num_samples * 1j
    x, y = _np.mgrid[xmin:xmax:NJ, ymin:ymax:NJ]
    return x.ravel(), y.ravel()


def _gen_planar_grid(
    xlim,
    ylim,
    zlim,
    center,
    face="x",
    num_samples=10,
    unit="mm",
    reverse_rotation=None,
):
    """Returns grids of points for two planes of a cuboid

    Args:
        xlim (tuple): min and max x values (float)
        ylim (tuple): min and max y values (float)
        zlim (tuple): min and max z values (float)
        face (str, optional): family of planes to generate (x, y, or z). Defaults to "x".
        num_samples (int, optional): Number of grid points to generate (10x10). Defaults to 10.
        unit (str, optional): Length scale. Defaults to "mm".

    Returns:
        tuple: points_lower, points_upper both (num_samples**2, 3) ndarrays
    """
    xc, yc, zc = center

    if face.lower() == "x":
        range_1 = ylim
        range_2 = zlim
        planes = xlim
    elif face.lower() == "y":
        range_1 = xlim
        range_2 = zlim
        planes = ylim

    elif face.lower() == "z":
        range_1 = xlim
        range_2 = ylim
        planes = zlim

    # Gets centroids of num_samples
    points_1, points_2 = _gen_grid(xlim=range_1, ylim=range_2, num_samples=num_samples)
    points_3_lower = _np.tile(planes[0], points_1.shape)
    points_3_upper = _np.tile(planes[1], points_1.shape)

    if reverse_rotation is not None:
        if face.lower() == "x":
            pos_vec = Quaternion._prepare_vector(points_3_lower, points_1, points_2)
            x, y, z = reverse_rotation * pos_vec
            points_lower = Point_Array3(x + xc, y + yc, z + zc, unit=unit)

            pos_vec = Quaternion._prepare_vector(points_3_upper, points_1, points_2)
            x, y, z = reverse_rotation * pos_vec
            points_upper = Point_Array3(x + xc, y + yc, z + zc, unit=unit)

        elif face.lower() == "y":
            pos_vec = Quaternion._prepare_vector(points_1, points_3_lower, points_2)
            x, y, z = reverse_rotation * pos_vec
            points_lower = Point_Array3(x + xc, y + yc, z + zc, unit=unit)

            pos_vec = Quaternion._prepare_vector(points_1, points_3_upper, points_2)
            x, y, z = reverse_rotation * pos_vec
            points_upper = Point_Array3(x + xc, y + yc, z + zc, unit=unit)

        elif face.lower() == "z":

            pos_vec = Quaternion._prepare_vector(points_1, points_2, points_3_lower)
            x, y, z = reverse_rotation * pos_vec
            points_lower = Point_Array3(x + xc, y + yc, z + zc, unit=unit)

            pos_vec = Quaternion._prepare_vector(points_1, points_2, points_3_upper)
            x, y, z = reverse_rotation * pos_vec
            points_upper = Point_Array3(x + xc, y + yc, z + zc, unit=unit)
    else:

        if face.lower() == "x":
            points_lower = Point_Array3(
                points_3_lower + xc, points_1 + yc, points_2 + zc, unit=unit
            )
            points_upper = Point_Array3(
                points_3_upper + xc, points_1 + yc, points_2 + zc, unit=unit
            )

        elif face.lower() == "y":
            points_lower = Point_Array3(
                points_1 + xc, points_3_lower + yc, points_2 + zc, unit=unit
            )
            points_upper = Point_Array3(
                points_1 + xc, points_3_upper + yc, points_2 + zc, unit=unit
            )

        elif face.lower() == "z":
            points_lower = Point_Array3(
                points_1 + xc, points_2 + yc, points_3_lower + zc, unit=unit
            )
            points_upper = Point_Array3(
                points_1 + xc, points_2 + yc, points_3_upper + zc, unit=unit
            )

    return points_lower, points_upper


def _calc_field_face(active_magnet, points):
    """Calculates the total force and torque acting on one face of a magnet due
    to all other instantiated magnets

    Args:
        active_magnet (Magnet3D): Target magnet
        points (Point_Array3): Grid of points (x,y,z)

    Returns:
        tuple: total_field (float), total_torque (float)
    """
    field = _allocate_field_array3(points.x, points.y, points.z)
    for magnet in Magnet3D.instances:
        if magnet is not active_magnet:
            Bx, By, Bz = magnet.get_field(points.x, points.y, points.z)
            field.x += Bx.reshape(field.x.shape)
            field.y += By.reshape(field.y.shape)
            field.z += Bz.reshape(field.z.shape)
    xc, yc, zc = active_magnet.get_center()
    field.x[~_np.isfinite(field.x)] = 0.0
    field.y[~_np.isfinite(field.y)] = 0.0
    field.z[~_np.isfinite(field.z)] = 0.0
    total_field = _np.array((_np.sum(field.x), _np.sum(field.y), _np.sum(field.z)))

    torques = _np.cross(
        _np.array(
            [points.x.ravel() - xc, points.y.ravel() - yc, points.z.ravel() - zc]
        ).T,
        _np.array([field.x.ravel(), field.y.ravel(), field.z.ravel()]).T,
    )
    total_torque = _np.sum(torques, axis=0)

    return total_field, total_torque


def calc_force_prism(active_magnet, num_samples=20, unit="mm"):
    """Calculates the total force on a cuboidal magnet due to all other instantiated magnets

    Args:
        active_magnet (Prism): Target Magnet
        num_samples (int, optional): Number of grid points per face. Defaults to 20.
        unit (str, optional): Length scale. Defaults to "mm".

    Returns:
        total_field (float), total_torque (float)
    """
    num_points_sq = num_samples * num_samples
    Jvec = active_magnet.get_Jr()
    Jnorm = active_magnet.Jr

    if _np.any(
        _np.fabs(
            _np.array(
                [
                    active_magnet.alpha_rad,
                    active_magnet.beta_rad,
                    active_magnet.gamma_rad,
                ]
            )
        )
        > active_magnet.tol
    ):

        _, reverse_rotation = active_magnet._generate_rotation_quaternions()
        Jrot = reverse_rotation * Jvec
    else:
        reverse_rotation = None
        Jrot = Jvec

    center = active_magnet.get_center()

    xlim, ylim, zlim, areas = _get_ranges_prism(active_magnet)
    force = _np.zeros(3)
    torque = _np.zeros(3)
    if _np.fabs(Jvec[0] / Jnorm) > FP_CUTOFF:
        points_lower, points_upper = _gen_planar_grid(
            xlim,
            ylim,
            zlim,
            center,
            face="x",
            num_samples=num_samples,
            unit=unit,
            reverse_rotation=reverse_rotation,
        )

        total_field, total_torque = _calc_field_face(active_magnet, points_lower)
        force += total_field * -Jrot[0] * areas[0] / num_points_sq
        torque += total_torque * -Jrot[0] * areas[0] / num_points_sq

        total_field, total_torque = _calc_field_face(active_magnet, points_upper)
        force += total_field * Jrot[0] * areas[0] / num_points_sq
        torque += total_torque * Jrot[0] * areas[0] / num_points_sq

    if _np.fabs(Jvec[1] / Jnorm) > FP_CUTOFF:
        points_lower, points_upper = _gen_planar_grid(
            xlim,
            ylim,
            zlim,
            center,
            face="y",
            num_samples=num_samples,
            unit=unit,
            reverse_rotation=reverse_rotation,
        )

        total_field, total_torque = _calc_field_face(active_magnet, points_lower)
        force += total_field * -Jrot[1] * areas[1] / num_points_sq
        torque += total_torque * -Jrot[1] * areas[1] / num_points_sq

        total_field, total_torque = _calc_field_face(active_magnet, points_upper)
        force += total_field * Jrot[1] * areas[1] / num_points_sq
        torque += total_torque * Jrot[1] * areas[1] / num_points_sq

    if _np.fabs(Jvec[2] / Jnorm) > FP_CUTOFF:
        points_lower, points_upper = _gen_planar_grid(
            xlim,
            ylim,
            zlim,
            center,
            face="z",
            num_samples=num_samples,
            unit=unit,
            reverse_rotation=reverse_rotation,
        )

        total_field, total_torque = _calc_field_face(active_magnet, points_lower)
        force += total_field * -Jrot[2] * areas[2] / num_points_sq
        torque += total_torque * -Jrot[2] * areas[2] / num_points_sq

        total_field, total_torque = _calc_field_face(active_magnet, points_upper)
        force += total_field * Jrot[2] * areas[2] / num_points_sq
        torque += total_torque * Jrot[2] * areas[2] / num_points_sq

    scaling_factor = get_unit_value_meter(points_lower.get_unit())
    force /= MU0 / scaling_factor / scaling_factor
    torque /= MU0 / scaling_factor / scaling_factor / scaling_factor

    return force, torque