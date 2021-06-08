from ..magnets import Magnet3D
from ..utils.global_const import FP_CUTOFF, MU0, ALIGN_CUTOFF
from ..utils._conversions import get_unit_value_meter
from ..utils._routines3D import _allocate_field_array3
from ..utils._vector_structs import Point_Array3
import numpy as _np
from numba import njit, guvectorize

__all__ = [
    "get_centroid",
    "triangle_area",
    "get_area_triangles",
    "get_midpoints",
    "divide_triangle_centroid",
    "divide_triangle_regular",
]


@njit
def get_centroid(triangle):
    """Gets the centroid of a triangle

    Args:
        triangle (ndarray): (3,3) array of vertices of a triangle in 3D

    Returns:
        ndarray: (3,) array of centroid coordinates
    """
    result = _np.zeros(3)
    for i in range(3):
        result[i] = (triangle[0, i] + triangle[1, i] + triangle[2, i]) / 3.0
    return result


# @guvectorize(['void(f8[:,:], f8[:])'],
#               '(x, x)->(x)')
# def get_centroid(triangle, result):
#     for i in range(3):
#         result[i] = (triangle[0,i] + triangle[1,i] + triangle[2,i])/3.0


@njit
def triangle_area(triangle):
    """Gets the area of a triangle. Computes the cross product area.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D

    Returns:
        float: area
    """
    return _np.linalg.norm(
        _np.cross((triangle[1] - triangle[0]), (triangle[2] - triangle[0])) / 2
    )


# @guvectorize(["void(f8[:,:, :], f8[:])"], "(x, y, y)->(x)")
@njit
def get_area_triangles(triangles, area):
    """Computes the area for an array of triangles

    Args:
        triangles (ndarray): (N,3,3) array of N triangles in 3D
        area (ndarray): (N,) array of areas
    """
    for i in range(len(triangles)):
        area[i] = _np.linalg.norm(
            _np.cross(
                (triangles[i][1] - triangles[i][0]), (triangles[i][2] - triangles[i][0])
            )
            / 2
        )


@njit
def _divide_triangle_centroid(triangle, depth=1, memo=_np.array(())):
    """Recursively divides a triangle into 3 using the centroid

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int, optional): Number of recursions to do. Defaults to 1.
        memo (ndarray, optional): Tracked array of generated vertices. Defaults to _np.array(()).

    Returns:
        ndarray: N*3*3 array of vertices. where N = 3**depth
    """

    centroid = _np.array([0.0, 0.0, 0.0])
    trig_array = _np.zeros((3, 3, 3))
    centroid = get_centroid(triangle)

    for i in range(3):
        centroid[i] = (triangle[0, i] + triangle[1, i] + triangle[2, i]) / 3
        trig_array[i, :] = triangle[:]

    trig_array[0, 2, :] = centroid
    trig_array[1, 0, :] = centroid
    trig_array[2, 1, :] = centroid
    depth -= 1

    if depth < 1:
        memo = _np.append(memo, trig_array)
    else:
        memo = _divide_triangle_centroid(trig_array[0, ...], depth=depth, memo=memo)
        memo = _divide_triangle_centroid(trig_array[1, ...], depth=depth, memo=memo)
        memo = _divide_triangle_centroid(trig_array[2, ...], depth=depth, memo=memo)

    return memo


@njit
def get_midpoints(triangle):
    """Get midpoints of the faces of a triangle

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D

    Returns:
        ndarray: (3,3) array of the three midpoints
    """
    midpoints = _np.zeros_like(triangle)
    for i in range(3):
        midpoints[i, ...] = (triangle[i] + triangle[(i + 1) % 3]) / 2
    return midpoints


@njit
def _divide_triangle_regular(triangle, depth=1, memo=_np.array(())):
    """Recursively divides a triangle into 4 using the midpoint of each face.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int, optional): Number of recursions to perform. Defaults to 1.
        memo (ndarray, optional): Tracked array of generated vertices. Defaults to _np.array(()).

    Returns:
        ndarray: (N*3*3,) array of vertices, where N = 4**depth
    """
    midpoints = get_midpoints(triangle)
    trig_array = _np.zeros((4, 3, 3))
    trig_array[0, 0, :] = triangle[0]
    trig_array[0, 1, :] = midpoints[0]
    trig_array[0, 2, :] = midpoints[2]

    trig_array[1, 0, :] = midpoints[0]
    trig_array[1, 1, :] = midpoints[1]
    trig_array[1, 2, :] = midpoints[2]

    trig_array[2, 0, :] = midpoints[2]
    trig_array[2, 1, :] = midpoints[1]
    trig_array[2, 2, :] = triangle[2]

    trig_array[3, 0, :] = midpoints[0]
    trig_array[3, 1, :] = triangle[1]
    trig_array[3, 2, :] = midpoints[1]
    depth -= 1

    if depth < 1:
        memo = _np.append(memo, trig_array)
    else:
        memo = _divide_triangle_regular(trig_array[0, ...], depth=depth, memo=memo)
        memo = _divide_triangle_regular(trig_array[1, ...], depth=depth, memo=memo)
        memo = _divide_triangle_regular(trig_array[2, ...], depth=depth, memo=memo)
        memo = _divide_triangle_regular(trig_array[3, ...], depth=depth, memo=memo)

    return memo


def divide_triangle_centroid(triangle, depth=1):
    """Recursively divides a triangle into 3 using the centroid.
    A wrapper over `_divide_triangle_centroid` to  reshape the output from
    (N*3*3,) to (N,3,3)

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int, optional): Number of recursions to perform. Defaults to 1.

    Returns:
        ndarray: (N,3,3) array of vertices, where N = 3**depth
    """
    mesh = _divide_triangle_centroid(triangle, depth=depth)
    mesh = mesh.reshape((mesh.shape[0] // 9, 3, 3))
    return mesh


def divide_triangle_regular(triangle, depth=1):
    """Recursively divides a triangle into 4 using the midpoint of each face.
    A wrapper over `_divide_triangle_regular` to  reshape the output from
    (N*3*3,) to (N,3,3)

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int, optional): Number of recursions to perform. Defaults to 1.

    Returns:
        ndarray: (N,3,3) array of vertices, where N = 4**depth
    """
    mesh = _divide_triangle_regular(triangle, depth=depth)
    mesh = mesh.reshape((mesh.shape[0] // 9, 3, 3))
    return mesh


def _calc_field_simplex(active_magnet, points):
    """Calculates the total force and torque acting on one simplex of a magnet due
    to all other instantiated magnets

    Args:
        active_magnet (Magnet3D): Target magnet
        points (Point_Array3): Grid of points (x,y,z)

    Returns:
        tuple: total_field (float), total_torque (float)
    """
    field = _allocate_field_array3(points.x, points.y, points.z)
    xc, yc, zc = active_magnet.centroid

    for magnet in Magnet3D.instances:
        if magnet is not active_magnet:
            Bx, By, Bz = magnet.get_field(points.x, points.y, points.z)
            field.x += Bx.reshape(field.x.shape)
            field.y += By.reshape(field.y.shape)
            field.z += Bz.reshape(field.z.shape)

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


def calc_force_mesh(active_magnet, depth=3, unit="mm"):
    force = _np.zeros(3)
    torque = _np.zeros(3)

    for i in range(len(active_magnet.mesh_vectors)):
        if _np.fabs(active_magnet.Jnorm[i] / active_magnet.Jr) > ALIGN_CUTOFF:
            triangle = active_magnet.mesh_vectors[i]
            mesh = divide_triangle_regular(triangle, depth=depth)
            centroids = _np.mean(mesh, axis=1)
            num_sub_triangles = len(centroids)

            points = Point_Array3(
                centroids[:, 0], centroids[:, 1], centroids[:, 2], unit=unit
            )
            local_field, local_torque = _calc_field_simplex(active_magnet, points)
            area = triangle_area(triangle)
            force += local_field * -active_magnet.Jnorm[i] * area / num_sub_triangles
            torque += local_torque * -active_magnet.Jnorm[i] * area / num_sub_triangles

    scaling_factor = get_unit_value_meter(points.get_unit())
    force /= MU0 / scaling_factor / scaling_factor
    torque /= MU0 / scaling_factor / scaling_factor / scaling_factor

    return force, torque