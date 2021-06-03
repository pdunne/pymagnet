__all__ = ["calc_force_cylinder"]
import numpy as _np
from ..utils.global_const import MU0, PI
from ..utils._conversions import get_unit_value_meter
from ..utils._vector_structs import Point_Array3
from ..utils._quaternion import Quaternion


from ._prism_force import _calc_field_face


def _get_ranges_cylinder(active_magnet):
    """Gets min and max x,y,z values for cylinder magnet faces.

    This is needed to generate a grid of points on each surface for calculating
    the force and torque on the magnet due to the other magnets.

    Args:
        active_magnet (Cylinder): Target magnet

    Returns:
        tuple: xlim (tuple), ylim (tuple), zlim (tuple), areas (ndarray) - areas of each face of the cuboid
    """
    radius, length = active_magnet.get_size()
    xlim = (-radius, radius)
    ylim = (-radius, radius)
    zlim = (-length / 2, length / 2)

    areas = _np.array([1, 1]) * PI * radius ** 2
    return xlim, ylim, zlim, areas


def _gen_cylinder_grid(active_magnet, xlim, ylim, num_segments=10):
    """Generates a 2D grid of N**2 points for a given pair of x and y limits.

    Args:
        active_magnet (Cylinder): Target magnet.
        xlim (tuple, optional): Min and max x values. Defaults to (0, 10).
        ylim (tuple, optional): Min and max y values. Defaults to (0, 10).
        num_rectangles (int, optional): Number of points to generate in each direction. Defaults to 10.

    Returns:
        ndarray: (num_rectangles**2, 3) array of points
    """
    NS = num_segments * 1j
    radius, _ = active_magnet.get_size()
    xc, yc, _ = active_magnet.get_center()
    x, y = _np.mgrid[xlim[0] : xlim[1] : NS, ylim[0] : ylim[1] : NS]
    index = (x - xc) ** 2 + (y - yc) ** 2 < radius ** 2
    x = x[index]
    y = y[index]

    return x, y


def _gen_cylinder_face_grid(
    active_magnet, xlim, ylim, zlim, num_segments=20, unit="mm", reverse_rotation=None
):
    """Returns grids of points for the two circular faces of a cylinder

    Args:
        xlim (tuple): min and max x values (float)
        ylim (tuple): min and max y values (float)
        zlim (tuple): min and max z values (float)
        num_segments (int, optional): Number of grid points to generate (10x10). Defaults to 20.
        unit (str, optional): Length scale. Defaults to "mm".

    Returns:
        tuple: points_lower, points_upper both (num_rectangles**2, 3) ndarrays
    """
    xc, yc, zc = active_magnet.get_center()

    x, y = _gen_cylinder_grid(active_magnet, xlim, ylim, num_segments=num_segments)
    z = _np.tile([zlim[0]], x.shape)

    if reverse_rotation is not None:
        pos_vec = Quaternion._prepare_vector(x, y, z)
        x_rot, y_rot, z_rot = reverse_rotation * pos_vec
        points_lower = Point_Array3(
            x_rot.ravel() + xc, y_rot.ravel() + yc, z_rot.ravel() + zc, unit=unit
        )

        z = _np.tile([zlim[1]], x.shape)
        pos_vec = Quaternion._prepare_vector(x, y, z)
        x_rot, y_rot, z_rot = reverse_rotation * pos_vec
        points_upper = Point_Array3(
            x_rot.ravel() + xc, y_rot.ravel() + yc, z_rot.ravel() + zc, unit=unit
        )

    else:
        points_lower = Point_Array3(
            x.ravel() + xc, y.ravel() + yc, z.ravel() + zc, unit=unit
        )
        z = _np.tile([zlim[1]], x.shape)
        points_upper = Point_Array3(
            x.ravel() + xc, y.ravel() + yc, z.ravel() + zc, unit=unit
        )

    return points_lower, points_upper


# def _calc_field_face(active_magnet, points):
#     field = _allocate_field_array3(points.x, points.y, points.z)
#     for magnet in Magnet3D.instances:
#         if magnet is not active_magnet:
#             #             print(magnet, active_magnet)
#             Bx, By, Bz = magnet.get_field(points.x, points.y, points.z)
#             field.x += Bx.reshape(field.x.shape)
#             field.y += By.reshape(field.y.shape)
#             field.z += Bz.reshape(field.z.shape)
#     xc, yc, zc = active_magnet.get_center()
#     field.x[~_np.isfinite(field.x)] = 0.0
#     field.y[~_np.isfinite(field.y)] = 0.0
#     field.z[~_np.isfinite(field.z)] = 0.0
#     total_field = _np.array((_np.sum(field.x), _np.sum(field.y), _np.sum(field.z)))
#
#     torques = _np.cross(
#         _np.array(
#             [points.x.ravel() - xc, points.y.ravel() - yc, points.z.ravel() - zc]
#         ).T,
#         _np.narray([field.x.ravel(), field.y.ravel(), field.z.ravel()]).T,
#     )
#     total_torque = _np.sum(torques, axis=0)
#
#     return total_field, total_torque


def calc_force_cylinder(active_magnet, num_segments=20, unit="mm"):
    """Calculates the total force on a cylinder magnet due to all other instantiated magnets

    Args:
        active_magnet (Cylinder): Target Magnet
        num_rectangles (int, optional): Number of grid points per face. Defaults to 20.
        unit (str, optional): Length scale. Defaults to "mm".

    Returns:
        total_field (float), total_torque (float)
    """
    Jvec = active_magnet.get_Jr()

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

    xlim, ylim, zlim, areas = _get_ranges_cylinder(active_magnet)

    points_lower, points_upper = _gen_cylinder_face_grid(
        active_magnet,
        xlim,
        ylim,
        zlim,
        num_segments=num_segments,
        unit=unit,
        reverse_rotation=reverse_rotation,
    )
    num_points = _np.size(points_lower.x)

    force = _np.zeros(3)
    torque = _np.zeros(3)

    total_field, total_torque = _calc_field_face(active_magnet, points_lower)

    force += total_field * -Jrot[2] * areas[0] / num_points
    torque += total_torque * -Jrot[2] * areas[0] / num_points

    total_field, total_torque = _calc_field_face(active_magnet, points_upper)
    force += total_field * Jrot[2] * areas[1] / num_points
    torque += total_torque * Jrot[2] * areas[1] / num_points

    scaling_factor = get_unit_value_meter(points_lower.get_unit())
    force /= MU0 / scaling_factor / scaling_factor
    torque /= MU0 / scaling_factor / scaling_factor / scaling_factor
    return force, torque