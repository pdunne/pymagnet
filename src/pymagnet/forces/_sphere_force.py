__all__ = ["calc_force_sphere"]

import numpy as _np
from ..utils.global_const import FP_CUTOFF, MU0, PI
from ..utils._conversions import get_unit_value_meter
from ..utils._vector_structs import Point_Array3
from ..utils._routines3D import _allocate_field_array3
from ..utils._quaternion import Quaternion
from ..magnets import Magnet3D


# def _gen_sphere_grid_cartesian(active_magnet, num_samples=200, unit="mm", thresh=1e-4):
#
#     radius = active_magnet.get_size()
#     xc, yc, zc = active_magnet.get_center()
#     Jr = active_magnet.Jr
#
#     NP = num_samples * 1j
#     x, y, z = _np.mgrid[
#         -radius + xc : radius + xc : NP,
#         -radius + yc : radius + yc : NP,
#         -radius + zc : radius + zc : NP,
#     ]
#
#     index = _np.abs((x ** 2 + y ** 2 + z ** 2) - radius ** 2) / radius ** 2 <= thresh
#
#     x = x[index]
#     y = y[index]
#     z = z[index]
#
#     Jnorm = Jr * z / radius
#     points = Point_Array3(x.ravel(), y.ravel(), z.ravel(), unit=unit)
#
#     return points, Jnorm.ravel()
#


def _gen_sphere_grid(active_magnet, num_samples=10, unit="mm"):

    NS = num_samples * 1j
    radius = active_magnet.get_size()
    Jr = active_magnet.Jr

    u, v = _np.mgrid[0 : 2 * PI * (num_samples - 1) / num_samples : NS, 0:PI:NS]

    x = radius * _np.cos(u) * _np.sin(v)
    y = radius * _np.sin(u) * _np.sin(v)
    z = radius * _np.cos(v)
    Jnorm = Jr * _np.cos(v)

    # Delete origin duplicates
    Jnorm = _np.delete(Jnorm, _np.arange(num_samples, Jnorm.size, num_samples)).ravel()

    y = _np.delete(y, _np.arange(num_samples, x.size, num_samples))
    z = _np.delete(z, _np.arange(num_samples, x.size, num_samples))
    x = _np.delete(x, _np.arange(num_samples, x.size, num_samples))
    points = Point_Array3(x.ravel(), y.ravel(), z.ravel(), unit=unit)

    return points, Jnorm


def _calc_field_face_sphere(active_magnet, points):
    field = _allocate_field_array3(points.x, points.y, points.z)
    for magnet in Magnet3D.instances:
        if magnet is not active_magnet:
            Bx, By, Bz = magnet.calcB(points.x, points.y, points.z)
            field.x += Bx
            field.y += By
            field.z += Bz
    xc, yc, zc = active_magnet.get_center()
    field.x[~_np.isfinite(field.x)] = 0.0
    field.y[~_np.isfinite(field.y)] = 0.0
    field.z[~_np.isfinite(field.z)] = 0.0
    fields = _np.array((field.x, field.y, field.z))

    torques = _np.cross(
        _np.array([points.x - xc, points.y - yc, points.z - zc]).T,
        _np.array([field.x, field.y, field.z]).T,
    ).T
    return fields, torques


def calc_force_sphere(active_magnet, num_samples=200, unit="mm"):

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

    points, Jn = _gen_sphere_grid(
        active_magnet,
        num_samples=num_samples,
        unit=unit,
        Jvec=Jrot,
        reverse_rotation=reverse_rotation,
    )

    area = 4 * PI * active_magnet.radius ** 2

    fields, torques = _calc_field_face_sphere(active_magnet, points)
    forces = Jn * fields

    torques = Jn * torques

    num_samples = _np.size(points.x)

    force = _np.sum(forces, axis=1) * area / num_samples
    torque = _np.sum(torques, axis=1) * area / num_samples

    scaling_factor = get_unit_value_meter(points.get_unit())
    force /= MU0 / scaling_factor / scaling_factor
    torque /= MU0 / scaling_factor / scaling_factor / scaling_factor
    return force, torque