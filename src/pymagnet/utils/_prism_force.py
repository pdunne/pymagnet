__all__ = ["calc_force_prism"]
import numpy as _np
from .global_const import FP_CUTOFF, MU0
from ._conversions import get_unit_value_meter
from ._vector_structs import Point_Array3
from ._routines3D import _allocate_field_array3
from ..magnets import Magnet3D

# from numba import njit, guvectorize


def _get_ranges_prism(active_magnet):
    size_x, size_y, size_z = active_magnet.get_size()
    xc, yc, zc = active_magnet.get_center()
    xlim = (-size_x / 2 + xc, size_x / 2 + xc)
    ylim = (-size_y / 2 + yc, size_y / 2 + yc)
    zlim = (-size_z / 2 + zc, size_z / 2 + zc)

    areas = _np.array([size_y * size_z, size_x * size_z, size_x * size_y])
    return xlim, ylim, zlim, areas


def _gen_grid(xlim=(0, 10), ylim=(0, 10), num_rectangles=10):
    xmin = xlim[0]
    xmax = xlim[1]
    ymin = ylim[0]
    ymax = ylim[1]

    dx = xmax - xmin
    dy = ymax - ymin

    xmin += dx / 2 / num_rectangles
    xmax -= dx / 2 / num_rectangles

    ymin += dy / 2 / num_rectangles
    ymax -= dy / 2 / num_rectangles

    NJ = num_rectangles * 1j
    x, y = _np.mgrid[xmin:xmax:NJ, ymin:ymax:NJ]
    return x, y


def _gen_planar_grid(xlim, ylim, zlim, face="x", num_rectangles=10, unit="mm"):
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

    # Gets centroids of num_rectangles
    points_1, points_2 = _gen_grid(
        xlim=range_1, ylim=range_2, num_rectangles=num_rectangles
    )
    points_3_lower = _np.tile(planes[0], points_1.shape)
    points_3_upper = _np.tile(planes[1], points_1.shape)

    if face.lower() == "x":
        points_lower = Point_Array3(points_3_lower, points_1, points_2)
        points_upper = Point_Array3(points_3_upper, points_1, points_2)

    elif face.lower() == "y":
        points_lower = Point_Array3(points_1, points_3_lower, points_2)
        points_upper = Point_Array3(points_1, points_3_upper, points_2)

    elif face.lower() == "z":
        points_lower = Point_Array3(points_1, points_2, points_3_lower, unit=unit)
        points_upper = Point_Array3(points_1, points_2, points_3_upper, unit=unit)
    return points_lower, points_upper


def _calc_field_face(active_magnet, points):
    field = _allocate_field_array3(points.x, points.y, points.z)
    for magnet in Magnet3D.instances:
        if magnet is not active_magnet:
            Bx, By, Bz = magnet.calcB(points.x, points.y, points.z)
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


def calc_force_prism(active_magnet, num_rectangles=20, unit="mm"):
    num_points_sq = num_rectangles * num_rectangles
    Jvec = active_magnet.get_Jr()
    Jnorm = active_magnet.Jr

    xlim, ylim, zlim, areas = _get_ranges_prism(active_magnet)
    force = _np.zeros(3)
    torque = _np.zeros(3)
    if Jvec[0] / Jnorm > FP_CUTOFF:
        points_lower, points_upper = _gen_planar_grid(
            xlim, ylim, zlim, face="x", num_rectangles=num_rectangles, unit=unit
        )

        total_field, total_torque = _calc_field_face(active_magnet, points_lower)
        force += total_field * -Jvec[0] * areas[0] / num_points_sq
        torque += total_torque * -Jvec[0] * areas[0] / num_points_sq

        total_field, total_torque = _calc_field_face(active_magnet, points_upper)
        force += total_field * Jvec[0] * areas[0] / num_points_sq
        torque += total_torque * Jvec[0] * areas[0] / num_points_sq

    if Jvec[1] / Jnorm > FP_CUTOFF:
        points_lower, points_upper = _gen_planar_grid(
            xlim, ylim, zlim, face="y", num_rectangles=num_rectangles, unit=unit
        )

        total_field, total_torque = _calc_field_face(active_magnet, points_lower)
        force += total_field * -Jvec[1] * areas[1] / num_points_sq
        torque += total_torque * -Jvec[1] * areas[1] / num_points_sq

        total_field, total_torque = _calc_field_face(active_magnet, points_upper)
        force += total_field * Jvec[1] * areas[1] / num_points_sq
        torque += total_torque * Jvec[1] * areas[1] / num_points_sq

    if Jvec[2] / Jnorm > FP_CUTOFF:
        points_lower, points_upper = _gen_planar_grid(
            xlim, ylim, zlim, face="z", num_rectangles=num_rectangles, unit=unit
        )

        total_field, total_torque = _calc_field_face(active_magnet, points_lower)
        force += total_field * -Jvec[2] * areas[2] / num_points_sq
        torque += total_torque * -Jvec[2] * areas[2] / num_points_sq

        total_field, total_torque = _calc_field_face(active_magnet, points_upper)
        force += total_field * Jvec[2] * areas[2] / num_points_sq
        torque += total_torque * Jvec[2] * areas[2] / num_points_sq

    scaling_factor = get_unit_value_meter(points_lower.get_unit())
    force /= MU0 / scaling_factor / scaling_factor
    torque /= MU0 / scaling_factor / scaling_factor

    return force, torque