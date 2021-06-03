# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
""" INCOMPLETE Routines for creation of known magnet assemblies consisting of Halbachs, quadrupoles, and others. This does not yet work.

"""
import numpy as _np
from ..utils.global_const import PI, PI_2
from matplotlib.path import Path as _Path


def mask_poly(num_sides, apothem, xc, yc, x, y):
    """A

    Args:
        num_sides (int): Number of sides in polygon
        apothem (float): Polygon apothem
        xc (float): Centre of polygon, x
        yc (float): Centre of polygon, y
        x (ndarray): array of x coordinates
        y (ndarray): array of y coordinates

    Returns:
        ndarray: boolean array of masked points
    """

    # Stack x and y arrays
    points = _np.vstack((x.flatten(), y.flatten())).T

    # Array of number of sides
    k = _np.arange(0, num_sides, 1)

    radius = apothem / _np.around(_np.cos(PI / num_sides), 4)

    # Generate polygon vertices
    xv = xc + radius * _np.sin(2 * PI * k / num_sides + _angle_offset(num_sides))
    yv = yc + radius * _np.cos(2 * PI * k / num_sides + _angle_offset(num_sides))

    # Needs to be a list for _Path function below
    poly_verts = _np.vstack((xv, yv)).T.tolist()

    # Create a _Path object following the vertices
    path = _Path(poly_verts)

    # generate a boolean array of points inside the polygon and reshape
    grid = path.contains_points(points)
    grid = grid.reshape(x.shape)
    return grid


def _angle_offset(num_sides):
    """Case structure for setting angle offset

    Args:
        num_sides (int): number of sides in polygon

    Returns:
        float: offset angle in radians
    """
    # Even number of sides
    if num_sides % 2 == 0:
        return {
            6: 0,
        }.get(num_sides, PI / num_sides)
    # odd number of sides
    else:
        return {
            3: PI_2,
            7: PI_2,
        }.get(num_sides, PI_2 / num_sides)


def mask_data_2D(magnet_half_width, num_sides, radius, B, x, y, **kwargs):
    """Applies mask to field vector.

    All field values outside the mask are set to NaN

    Args:
        magnet_half_width (float): Half the magnet width for an array
        num_sides (int): Number of sides
        radius (float): distance from origin to centre of each magnet
        B (Vector2): Magnetic field vector
        x (ndarray): array of x coordinates
        y (ndarray): array of y coordinates
        **kwargs: keyword arguments

    Kwargs:
        center (tuple): centre of mask, defaults to (.0, .0)
        mask_function (function): mask function, defaults to `mask_poly`, can also be
        `radial_mask`

    Returns:
        Vector2: new vector of masked points
    """
    from ..utils._point_structs import Vector2

    apothem = radius - magnet_half_width
    center = kwargs.pop("center", (0.0, 0.0))
    mask_function = kwargs.pop("mask_function", mask_poly)
    xp = center[0]
    yp = center[1]

    grid = mask_function(num_sides, apothem, xp, yp, x, y)
    Bm = Vector2(B.x, B.y)
    Bm.x[~grid] = _np.nan
    Bm.y[~grid] = _np.nan
    Bm.calc_norm()
    return Bm


def radial_mask(radius, xp, yp, x, y):
    """Applies circular mask

    Args:
        radius (float): mask radius
        xp (float): center of circle, x
        yp (float): center of circle, y
        x (ndarray): array of x points
        y (ndarray): array of y points

    Returns:
        array: boolean array
    """
    r_data = _np.sqrt(_np.power(x - xp, 2) + _np.power(y - yp, 2))
    grid = (r_data) <= radius
    return grid


def radial_profile(data, center):
    """Calculates radial profile

    Args:
        data (ndarray): Data to be profiled must be a 2D array
        center (tuple): [description]

    Returns:
        [type]: [description]
    """
    # Generate a grid of index points depending on the shape of the data
    y, x = _np.indices((data.shape))

    # Generate a radial function in units of index
    r = _np.sqrt(_np.power(x - center[0], 2) + _np.power(y - center[1], 2))

    # Cast to integer
    r = r.astype(_np.int)

    # Bin the data
    tbin = _np.bincount(r.ravel(), data.ravel())
    nr = _np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile


def radial_extraction_2D(magnet_half_width, radius, x, y, B):
    """Returns the radial average of the normal component of 2D vector field

    Args:
        magnet_half_width (float): half the magnet width
        radius (float): Radius
        x (ndarray): array of x points
        y (ndarray): array of y points
        B (Vector2): Field vector

    Returns:
        tuple: array of radial points, array of radial average
    """

    # Get center indices
    index_x = x.shape[0] // 2
    index_y = y.shape[1] // 2

    # Generate radial length in units of index
    radial_length = _np.sqrt(
        _np.power(x[0, -1] - x[index_x, index_y], 2)
        + _np.power(y[0, -1] - y[index_x, index_y], 2)
    )

    # Get radial average of the normalised field
    radial_average = radial_profile(B.n, [index_x, index_y])

    # Calculate the x value for plotting `xpl`
    xpl = _np.linspace(0, radial_length, len(radial_average))

    # Discard data outside the mask
    radial_average = radial_average[xpl <= radius - magnet_half_width]
    xpl = xpl[xpl <= radius - magnet_half_width]
    return xpl, radial_average


def init_magnets(num_magnets=4, b_scale=1, assem_type="halbach"):
    from .. import reset_magnets
    from ..magnets._magnet2D import Rectangle
    from ._routines2D import grid2D

    """Initialise 

    Args:
        N (int, optional): [Number of magnets]. Defaults to 4.
        b_scale (int, optional): [ratio of b/a]. Defaults to 1.
        assem_type (str, optional): [halbach, 2mag, pseudo_quad]. Defaults to 'halbach'.
    Returns:
        mag_prop (dictionary): Magnet properties
        grid_prop (dictionary): Calculation grid (x,y)
    """
    reset_magnets()
    NP = 201
    width = 2
    height = width * b_scale
    Jr = 1.0

    if assem_type == "halbach".lower():
        theta = 90.0
        alpha = 0.0

        # assembly radius
        radius = (width * 2) / (
            2 * _np.tan(PI / num_magnets)
        ) + width * 1.0  # Assembly radius
        Jr = _np.array([Jr * (-1) ** num for num in range(num_magnets)])
        theta = _np.tile(theta, num_magnets)
        alpha = _np.tile(alpha, num_magnets)

        # Halbach - 4 element
        steps = _np.arange(0, num_magnets, 1)

        alpha = steps * 360 / num_magnets
        xc = -radius * _np.cos(_np.deg2rad(steps * 360 / num_magnets))
        yc = radius * _np.sin(_np.deg2rad(steps * 360 / num_magnets))

        for i in range(num_magnets):
            _ = Rectangle(
                width=width,
                height=height,
                Jr=Jr[i],
                center=(xc[i], yc[i]),
                theta=theta[i],
                alpha=alpha[i],
            )
        off = -width / 2
        UP = radius + off
        UPx = UP
        UPy = UP

    else:
        hGap = width / 2
        xc = width / 2 + hGap
        yc = 0
        theta = 0.0
        alpha = 0.0
        N = 2

        #  2 element
        _ = Rectangle(
            width=width, height=height, Jr=Jr, center=(xc, yc), theta=theta, alpha=alpha
        )
        _ = Rectangle(
            width=width,
            height=height,
            Jr=-Jr,
            center=(-xc, yc),
            theta=theta,
            alpha=alpha,
        )
        UPx = hGap
        UPy = width / 2

    x, y = grid2D(-UPx, UPx, -UPy, UPy, NP)

    mag_prop = {
        "width": width,
        "height": height,
        "xc": xc,
        "yc": yc,
        "NP": NP,
        "Jr": Jr,
        "theta": theta,
        "alpha": alpha,
        "N": num_magnets,
        "radius": radius,
        "assem_type": assem_type.lower(),
    }

    grid_prop = {"x": x, "y": y}

    return mag_prop, grid_prop


def calc_magnetic_field(mag_prop, grid_prop):
    from ..utils._routines2D import get_field_2D

    """Run magnetic calculation and mask data to central bore

    Args:
        mag_prop (dictionary): Magnet properties
        grid_prop (dictionary): Calculation grid (x,y)
    
    Updates grid_prop with magnetic field object
    and masked magnetic field object
    """
    width = mag_prop["width"]
    num_magnets = mag_prop["num_magnets"]

    x = grid_prop["x"]
    y = grid_prop["y"]
    B = get_field_2D(x, y)
    if mag_prop["assem_type"].lower() == "halbach":
        mask_radius = mag_prop["radius"]
    else:
        mask_radius = mag_prop["hGap"]
        N = 4

    Bm = mask_data_2D(width, num_magnets, mask_radius, B, x, y)
    grid_prop.update({"B": B, "Bm": Bm})


# def fit_power(N):
#     """Switch case dictionary for returning powers needed
#     to fit radially averaged magnetic field inside a bore.
#     Defaults to 3 (i.e. for fitting third order polynomial)

#     Args:
#         N (int): [Number of magnets]

#     Returns:
#         [int]: [polynomial order for fit]
#     """
#     return{
#            2: 1,
#            4: 1,
#            6: 2
#            }.get(N, 3)

# def fit_radial_profile(mag_prop, grid_prop):
#     """Calculate radial average inside cavity and fit it with an N-order
#     polynomial

#     Args:
#         mag_prop (dict): magnetic properties
#         grid_prop (dict): Grid properties for calculations (x,y, B, Bm)

#     Returns:
#         radial_param (dict): Radial B, BgradB
#     """

#     x = grid_prop['x']
#     y = grid_prop['y']
#     B = grid_prop['B']
#     a = mag_prop['a']
#     N = mag_prop['N']
#     centre = 0


#     xpl,B_ravg = mag.radial_extraction(x,y,B,a,centre)
#     xpl = xpl[np.isfinite(B_ravg)]
#     B_ravg = B_ravg[np.isfinite(B_ravg)]
#     fit_res = np.polyfit(xpl,B_ravg, fit_power(N))
#     B_fit = np.polyval(fit_res, xpl)
#     c = fit_res[-1]
#     m = fit_res[-2]

#     # num_list = [str(num) for num in fit_res]


#     BgB_fit = B_ravg*m

#     radial_param = {"xpl": xpl,
#                   "B_ravg": B_ravg,
#                   "fit_res": fit_res,
#                   "B_fit": B_fit,
#                   'BgB_fit':BgB_fit}


#     if fit_res.shape[0] > 2:
#         m_approx = np.gradient(B_ravg,xpl).mean()
#         B_fit2 = np.polyval( [m_approx, fit_res[-1]], xpl)
#         if fit_res.shape[0] > 3:
#             BgB_fit = B_fit *( 3*xpl**2 * fit_res[0] + 2*xpl * fit_res[1] + fit_res[2])
#         else:
#             BgB_fit = B_fit *( 2*xpl* fit_res[0] + fit_res[1])
#         m_approx = np.gradient(BgB_fit,xpl).mean()
#         BgB_fit2 = B_fit2*m_approx
#         radial_param = {"xpl": xpl,
#                     "B_ravg": B_ravg,
#                     "fit_res": fit_res,
#                     "B_fit": B_fit,
#                     "B_fit2": B_fit2,
#                     'BgB_fit': BgB_fit,
#                     'BgB_fit2': BgB_fit2
#                     }
#     return radial_param