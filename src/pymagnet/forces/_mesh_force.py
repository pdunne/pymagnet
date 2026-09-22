import numpy as _np
from numba import njit, prange

from ..magnets import Magnet3D
from ..utils._conversions import get_unit_value_meter
from ..utils._routines3D import _allocate_field_array3
from ..utils._vector_structs import Point_Array3
from ..utils.global_const import ALIGN_CUTOFF, MU0


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


@njit(cache=True)
def triangle_area(triangle):
    """Gets the area of a triangle. Computes the cross product area.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D

    Returns:
        float: area
    """
    # Edge vectors
    e1_x = triangle[1, 0] - triangle[0, 0]
    e1_y = triangle[1, 1] - triangle[0, 1]
    e1_z = triangle[1, 2] - triangle[0, 2]

    e2_x = triangle[2, 0] - triangle[0, 0]
    e2_y = triangle[2, 1] - triangle[0, 1]
    e2_z = triangle[2, 2] - triangle[0, 2]

    # Cross product
    cx = e1_y * e2_z - e1_z * e2_y
    cy = e1_z * e2_x - e1_x * e2_z
    cz = e1_x * e2_y - e1_y * e2_x

    # Area = |cross| / 2
    return _np.sqrt(cx * cx + cy * cy + cz * cz) / 2.0


# @guvectorize(["void(f8[:,:, :], f8[:])"], "(x, y, y)->(x)")
# @njit
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


def _divide_triangle_centroid_fast(triangle, depth=1):
    """Divides a triangle into 3**depth sub-triangles using centroid subdivision.

    Uses pre-allocated arrays and iterative processing for better performance.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int, optional): Number of subdivision levels. Defaults to 1.

    Returns:
        ndarray: (3**depth, 3, 3) array of sub-triangle vertices
    """
    # Use a working buffer: start with input triangle, expand each iteration
    current = _np.array([triangle])

    for _ in range(depth):
        n = len(current)
        next_level = _np.zeros((n * 3, 3, 3))

        for i in range(n):
            tri = current[i]
            centroid = (tri[0] + tri[1] + tri[2]) / 3.0

            # Sub-triangle 0: v0, v1, centroid
            next_level[i * 3, 0, :] = tri[0]
            next_level[i * 3, 1, :] = tri[1]
            next_level[i * 3, 2, :] = centroid

            # Sub-triangle 1: centroid, v1, v2
            next_level[i * 3 + 1, 0, :] = centroid
            next_level[i * 3 + 1, 1, :] = tri[1]
            next_level[i * 3 + 1, 2, :] = tri[2]

            # Sub-triangle 2: v0, centroid, v2
            next_level[i * 3 + 2, 0, :] = tri[0]
            next_level[i * 3 + 2, 1, :] = centroid
            next_level[i * 3 + 2, 2, :] = tri[2]

        current = next_level

    return current


def get_midpoints(triangle):
    """Get midpoints of the edges of a triangle.

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
def _get_midpoints_njit(triangle):
    """Get midpoints of the edges of a triangle (numba-compatible).

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D

    Returns:
        ndarray: (3,3) array of the three midpoints
    """
    midpoints = _np.zeros((3, 3))
    for i in range(3):
        j = (i + 1) % 3
        for k in range(3):
            midpoints[i, k] = (triangle[i, k] + triangle[j, k]) / 2.0
    return midpoints


@njit
def _subdivide_triangle_regular_njit(tri, out, out_idx):
    """Subdivide a single triangle into 4 sub-triangles (numba-compatible).

    Args:
        tri (ndarray): (3,3) input triangle
        out (ndarray): output array to write results to
        out_idx (int): starting index in output array

    Returns:
        None (modifies out in-place)
    """
    midpoints = _get_midpoints_njit(tri)

    # Sub-triangle 0: v0, m01, m20
    for k in range(3):
        out[out_idx, 0, k] = tri[0, k]
        out[out_idx, 1, k] = midpoints[0, k]
        out[out_idx, 2, k] = midpoints[2, k]

    # Sub-triangle 1: m01, m12, m20 (center triangle)
    for k in range(3):
        out[out_idx + 1, 0, k] = midpoints[0, k]
        out[out_idx + 1, 1, k] = midpoints[1, k]
        out[out_idx + 1, 2, k] = midpoints[2, k]

    # Sub-triangle 2: m20, m12, v2
    for k in range(3):
        out[out_idx + 2, 0, k] = midpoints[2, k]
        out[out_idx + 2, 1, k] = midpoints[1, k]
        out[out_idx + 2, 2, k] = tri[2, k]

    # Sub-triangle 3: m01, v1, m12
    for k in range(3):
        out[out_idx + 3, 0, k] = midpoints[0, k]
        out[out_idx + 3, 1, k] = tri[1, k]
        out[out_idx + 3, 2, k] = midpoints[1, k]


@njit
def _divide_triangle_regular_fast_njit(triangle, depth):
    """Divides a triangle into 4**depth sub-triangles using midpoint subdivision.

    Uses pre-allocated arrays and iterative processing for better performance.
    Compiled with numba for additional speedup.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int): Number of subdivision levels

    Returns:
        ndarray: (4**depth, 3, 3) array of sub-triangle vertices
    """
    # Start with the input triangle
    current_count = 1
    current = _np.zeros((1, 3, 3))
    for i in range(3):
        for j in range(3):
            current[0, i, j] = triangle[i, j]

    for _ in range(depth):
        next_count = current_count * 4
        next_level = _np.zeros((next_count, 3, 3))

        for i in range(current_count):
            _subdivide_triangle_regular_njit(current[i], next_level, i * 4)

        current = next_level
        current_count = next_count

    return current


def _divide_triangle_regular_fast(triangle, depth=1):
    """Divides a triangle into 4**depth sub-triangles using midpoint subdivision.

    Wrapper around the numba-compiled implementation.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int, optional): Number of subdivision levels. Defaults to 1.

    Returns:
        ndarray: (4**depth, 3, 3) array of sub-triangle vertices
    """
    return _divide_triangle_regular_fast_njit(triangle.astype(_np.float64), depth)


def divide_triangle_centroid(triangle, depth=1):
    """Divides a triangle into 3**depth sub-triangles using centroid subdivision.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int, optional): Number of subdivision levels. Defaults to 1.

    Returns:
        ndarray: (N,3,3) array of vertices, where N = 3**depth
    """
    return _divide_triangle_centroid_fast(triangle, depth=depth)


def divide_triangle_regular(triangle, depth=1):
    """Divides a triangle into 4**depth sub-triangles using midpoint subdivision.

    Args:
        triangle (ndarray): (3,3) array of triangle vertices in 3D
        depth (int, optional): Number of subdivision levels. Defaults to 1.

    Returns:
        ndarray: (N,3,3) array of vertices, where N = 4**depth
    """
    return _divide_triangle_regular_fast(triangle, depth=depth)


def _precompute_force_centroids(mesh_vectors, Jnorm, Jr, depth):
    """Pack sub-triangle centroids and metadata for all active triangles.

    Args:
        mesh_vectors (ndarray): (N, 3, 3) triangle vertices
        Jnorm (ndarray): (N,) normal component of magnetisation per triangle
        Jr (float): remanence magnitude
        depth (int): subdivision depth (4**depth sub-triangles per triangle)

    Returns:
        tuple: (centroids_flat, Jnorm_active, areas_active, n_sub)
            - centroids_flat (ndarray): (n_active * n_sub, 3) centroid coordinates
            - Jnorm_active (ndarray): (n_active,) Jnorm values for active triangles
            - areas_active (ndarray): (n_active,) triangle areas
            - n_sub (int): sub-triangles per active triangle
    """
    n_sub = 4 ** depth
    active = _np.fabs(Jnorm / Jr) > ALIGN_CUTOFF
    idx = _np.where(active)[0]
    n_active = len(idx)

    centroids_flat = _np.empty((n_active * n_sub, 3), dtype=_np.float64)
    Jnorm_active = _np.empty(n_active, dtype=_np.float64)
    areas_active = _np.empty(n_active, dtype=_np.float64)

    for k, i in enumerate(idx):
        sub_tris = _divide_triangle_regular_fast_njit(
            _np.ascontiguousarray(mesh_vectors[i], dtype=_np.float64), depth
        )
        centroids_flat[k * n_sub : (k + 1) * n_sub] = sub_tris.mean(axis=1)
        Jnorm_active[k] = Jnorm[i]
        areas_active[k] = triangle_area(mesh_vectors[i])

    return centroids_flat, Jnorm_active, areas_active, n_sub


@njit(parallel=True, cache=True)
def _accumulate_force_torque_njit(
    centroids, Bx, By, Bz, Jnorm_active, areas, n_sub, xc, yc, zc
):
    """Accumulate force and torque from pre-computed field values.

    Uses prange over active triangles with scalar reductions.

    Args:
        centroids (ndarray): (n_active * n_sub, 3) centroid positions
        Bx, By, Bz (ndarray): (n_active * n_sub,) field components
        Jnorm_active (ndarray): (n_active,) Jnorm values
        areas (ndarray): (n_active,) triangle areas
        n_sub (int): sub-triangles per active triangle
        xc, yc, zc (float): centroid of active magnet (torque reference point)

    Returns:
        tuple: (force (3,), torque (3,))
    """
    n_active = Jnorm_active.shape[0]
    fx = 0.0
    fy = 0.0
    fz = 0.0
    tx = 0.0
    ty = 0.0
    tz = 0.0

    for k in prange(n_active):
        s = k * n_sub
        scale = -Jnorm_active[k] * areas[k] / n_sub
        sfx = 0.0
        sfy = 0.0
        sfz = 0.0
        stx = 0.0
        sty = 0.0
        stz = 0.0
        for j in range(s, s + n_sub):
            bx = Bx[j]
            by = By[j]
            bz = Bz[j]
            rx = centroids[j, 0] - xc
            ry = centroids[j, 1] - yc
            rz = centroids[j, 2] - zc
            sfx += bx
            sfy += by
            sfz += bz
            stx += ry * bz - rz * by
            sty += rz * bx - rx * bz
            stz += rx * by - ry * bx
        fx += sfx * scale
        fy += sfy * scale
        fz += sfz * scale
        tx += stx * scale
        ty += sty * scale
        tz += stz * scale

    force = _np.array([fx, fy, fz])
    torque = _np.array([tx, ty, tz])
    return force, torque


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
    # 1. Pre-compute all sub-triangle centroids for active triangles
    centroids, Jnorm_active, areas, n_sub = _precompute_force_centroids(
        active_magnet.mesh_vectors, active_magnet.Jnorm, active_magnet.Jr, depth
    )

    if len(centroids) == 0:
        return _np.zeros(3), _np.zeros(3)

    # 2. Evaluate field from all source magnets at all centroids (single batched call each)
    cx = centroids[:, 0]
    cy = centroids[:, 1]
    cz = centroids[:, 2]

    Bx = _np.zeros(len(centroids))
    By = _np.zeros(len(centroids))
    Bz = _np.zeros(len(centroids))
    for magnet in Magnet3D.instances:
        if magnet is not active_magnet:
            bx, by, bz = magnet.get_field(cx, cy, cz)
            Bx += bx.ravel()
            By += by.ravel()
            Bz += bz.ravel()

    Bx[~_np.isfinite(Bx)] = 0.0
    By[~_np.isfinite(By)] = 0.0
    Bz[~_np.isfinite(Bz)] = 0.0

    # 3. Accumulate force/torque in parallel over active triangles
    xc, yc, zc = active_magnet.centroid
    force, torque = _accumulate_force_torque_njit(
        centroids, Bx, By, Bz, Jnorm_active, areas, n_sub,
        float(xc), float(yc), float(zc),
    )

    # 4. Apply unit scaling
    scaling_factor = get_unit_value_meter(unit)
    assert scaling_factor is not None
    force = force / (MU0 / scaling_factor / scaling_factor)
    torque = torque / (MU0 / scaling_factor / scaling_factor / scaling_factor)
    return force, torque
