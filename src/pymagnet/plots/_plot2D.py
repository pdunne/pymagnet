# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Plotting routines

This module contains all functions needed to plot lines and contours for 2D
magnetic sources, and 

"""

__all__ = ["plot_2D_line", "plot_2D_contour", "plot_3D_contour", "plot_sub_contour_3D"]

from .. import magnets as _mag
from ..utils._conversions import get_unit_value_meter, get_unit_value_tesla
from ..utils._vector_structs import Field2, Point_Array2

import numpy as _np
import matplotlib.pyplot as _plt
import matplotlib.cm as _cm
from matplotlib.patches import Rectangle as _Rect
from matplotlib.patches import Circle as _Circ
from matplotlib.patches import Arrow as _Arrow
from matplotlib.transforms import Affine2D


class patch(object):
    """Encodes magnet dimensions for drawing on plots"""

    def __init__(self, x, y, width, height, transform, type):
        """Initialse a patch

        Args:
            x (float): centre, x
            y (float): center, y
            width (float): width
            height (float): width
            transform (matplotlib Affine2D): transform object, `rotate_deg_around`
        """
        super().__init__()
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.transform = transform
        self.type = type

    def __repr__(self) -> str:
        return f"(x: {self.x}, y: {self.y} w:{self.width}, h: {self.height})"

    def __str__(self) -> str:
        return f"(x: {self.x}, y: {self.y} w:{self.width}, h: {self.height})"


class arrow(object):
    """Encodes magnetisation vector for drawing on plots"""

    def __init__(self, x, y, dx, dy, transform, width=3):
        """Init Arrow

        Args:
            x (float): arrow tail, x
            y (float): arrow tail, y
            dx (float): arrow head displacement, x
            dy (float): arrow head displacement, y
            transform (matplotlib Affine2D): transformation object, `translate`
            width (int, optional): Arrow width. Defaults to 3.
        """
        super().__init__()
        self.x = x
        self.y = y
        self.dx = dx
        self.dy = dy
        self.transform = transform
        self.width = width

    def __repr__(self) -> str:
        return (
            f"(x: {self.x}, y: {self.y}, dx: {self.dx}, dy: {self.dy}, w:{self.width})"
        )

    def __str__(self) -> str:
        return (
            f"(x: {self.x}, y: {self.y}, dx: {self.dx}, dy: {self.dy}, w:{self.width})"
        )


class magnet_patch(object):
    """Magnet drawing class"""

    def __init__(self, patch, arrow) -> None:
        super().__init__()
        self.patch = patch
        self.arrow = arrow

    def __str__(self) -> str:
        return self.patch.__str__() + self.arrow.__str__()


def plot_2D_line(point_array, field, **kwargs):
    """Line Plot of field from 2D magnet

    Args:
        point_array (Point_Array2): coordinates
        field (Field2): Magnetic Field

    Kwargs:
        xlab (str): xlabel
        ylab (str): ylabel
        axis_scale (str): unused
        save_fig (bool): Save to png file. Defaults to False

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """

    xlab = kwargs.pop("xlab", f"x ({point_array.unit})")
    ylab = kwargs.pop("ylab", f"B ({field.unit})")
    axis_scale = kwargs.pop("axis_scale", "equal")

    SAVE = kwargs.pop("save_fig", False)

    fig, ax = _plt.subplots(figsize=(8, 8))
    _plt.plot(point_array.x, field.n, label=r"$|\mathbf{B}|$")
    _plt.plot(point_array.x, field.x, label=r"$B_x$")
    _plt.plot(point_array.x, field.y, label=r"$B_y$")
    _plt.legend(loc="best")
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)

    _plt.show()

    if SAVE:
        _plt.savefig("line_plot.png", dpi=300)
        # _plt.savefig('contour_plot.pdf', dpi=300)
    return fig, ax


def plot_2D_contour(point_array, field, **kwargs):
    """Contour plot of field

    Args:
        point_array (Point_Array2): coordinates
        field (Field2): Magnetic Field

    Kwargs:
        save_fig (bool): Save to png file. Defaults to False
        xlab (str): x axis label
        ylab (str): y axis label
        clab (str): label for colorbar
        axis_scale (str): axis aspect ratio Defaults to 'equal'.
        show_magnets (bool): Draw magnets. Defaults to True
        field_component (str): Defaults to 'n'.
        plot_type (str): Draw `contour` or `streamplot`. Defaults to 'contour'
        cmap (str): Colormap. Defaults to `viridis`
        vector_plot (bool): Draw vectors as arrows. Defaults to False.
        vector_color (str): Arrow color. Defaults to 'w'
        cmin (float): Color scale minimum. Defaults to 0.0
        cmax (float): Color scale minimum. Defaults to twice the mean field
        num_levels (int): Number of contour levels. Defaults to 11.

    Raises:
        Exception: plot_type must be 'contour' or 'streamplot'

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    import matplotlib.cm as _cm
    from ..magnets._polygon2D import PolyMagnet

    show_magnets = kwargs.pop("show_magnets", True)

    xlab = kwargs.pop("xlab", f"x ({point_array.unit})")
    ylab = kwargs.pop("ylab", f"y ({point_array.unit})")
    clab = kwargs.pop("clab", f"B ({field.unit})")
    axis_scale = kwargs.pop("axis_scale", "equal")

    SAVE = kwargs.pop("save_fig", False)

    field_component = kwargs.pop("field_component", "n")

    plot_type = kwargs.pop("plot_type", "contour")
    fig, ax = _plt.subplots(figsize=(8, 8))

    if plot_type.lower() == "contour":
        cmap = kwargs.pop("cmap", "viridis")
        vector_plot = kwargs.pop("vector_plot", False)
        vector_color = kwargs.pop("vector_color", "w")
        NQ = kwargs.pop("vector_arrows", 11)

        if field_component == "x":
            field_chosen = field.x
        elif field_component == "y":
            field_chosen = field.y
        else:
            field_chosen = field.n

        finite_field = field_chosen[_np.isfinite(field_chosen)]
        cmin = kwargs.pop("cmin", 0)
        cmax = kwargs.pop("cmax", round(finite_field.mean() * 2, 1))
        num_levels = kwargs.pop("num_levels", 11)

        lev2 = _np.linspace(cmin, cmax, 256, endpoint=True)

        CS = _plt.contourf(
            point_array.x,
            point_array.y,
            field_chosen,
            levels=lev2,
            cmap=_plt.get_cmap(cmap),
            extend="max",
        )

        # Draw contour lines
        if num_levels > 1:
            lev1 = _np.linspace(cmin, cmax, num_levels, endpoint=True)
            _ = _plt.contour(
                point_array.x,
                point_array.y,
                field_chosen,
                vmin=cmin,
                vmax=cmax,
                levels=lev1,
                linewidths=1.0,
                colors="k",
            )
            CB = _plt.colorbar(CS, ticks=lev1)
        else:
            CB = _plt.colorbar(CS)

        # Draw field vectors
        if vector_plot:
            _vector_plot2(point_array, field, NQ, vector_color)

    elif plot_type.lower() == "streamplot":

        xpl = point_array.x[:, 0]
        ypl = point_array.y[0, :]
        cmap = kwargs.pop("cmap", None)

        if cmap is not None:
            cmin = kwargs.pop("cmin", -round(_np.nanmean(field.n), 1))
            cmax = kwargs.pop("cmax", round(_np.nanmean(field.n), 1))
            stream_shading = kwargs.pop("stream_color", "vertical")
            norm = _cm.colors.Normalize(vmin=cmin, vmax=cmax)

            stream_dict = {
                "normal": field.n.T,
                "horizontal": field.x.T,
                "vertical": field.y.T,
            }

            CS = _plt.streamplot(
                xpl,
                ypl,
                field.x.T / field.n.T,
                field.y.T / field.n.T,
                color=stream_dict.get(stream_shading, "normal"),
                density=1.2,
                norm=norm,
                cmap=cmap,
                linewidth=0.5,
            )
            CB = _plt.colorbar(CS.lines)

        else:
            color = kwargs.pop("color", "k")
            CS = _plt.streamplot(
                xpl,
                ypl,
                field.x.T / field.n.T,
                field.y.T / field.n.T,
                density=1.2,
                linewidth=0.5,
                color=color,
            )
            CB = None

    else:
        raise Exception("plot_type must be 'contour' or 'streamplot'")

    # Draw magnets and magnetisation arrows
    if show_magnets:
        _draw_magnets2(ax)
        if len(PolyMagnet.instances) > 0:
            for magnet in PolyMagnet.instances:
                poly = _plt.Polygon(
                    _np.array(magnet.polygon.vertices),
                    ec="k",
                    fc="w",
                    zorder=5,
                )
                ax.add_patch(poly)

    if CB is not None:
        CB.ax.get_yaxis().labelpad = 15
        CB.ax.set_ylabel(clab, rotation=270)
    _plt.axis(axis_scale)
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)
    _plt.show()
    fig.tight_layout()
    ax.axis("scaled")
    if SAVE:
        _plt.savefig("contour_plot.png", dpi=300)

    return fig, ax


def _num_patch_2D():
    """Generates patches and arrows for drawing

    Returns:
        tuple: (list, list) lists of patch and arrow objects
    """
    from ..magnets._magnet2D import Magnet2D, Rectangle, Circle

    patch_array = []
    for magnet in Magnet2D.instances:
        if issubclass(magnet.__class__, Rectangle):
            magnet_patch_tmp = _gen_rect_patch(magnet)
            patch_array.append(magnet_patch_tmp)
        elif issubclass(magnet.__class__, Circle):
            magnet_patch_tmp = _gen_circle_patch(magnet)
            patch_array.append(magnet_patch_tmp)

    return patch_array


def _gen_rect_patch(magnet):
    """Generates rectangular patches and arrows to represent magnets for 2D plots

    Args:
        magnet (Magnet2D): instance of a magnet class
    Returns:
        magnet_patch: magnet_patch data structure
    """
    # from matplotlib.transforms import Affine2D

    patch_tmp = patch(
        x=(magnet.center[0] - magnet.a),
        y=(magnet.center[1] - magnet.b),
        width=(2 * magnet.a),
        height=(2 * magnet.b),
        transform=Affine2D().rotate_deg_around(
            (magnet.center[0]),
            (magnet.center[1]),
            -magnet.alpha,
        ),
        type="rectangle",
    )

    Jnorm = magnet.get_Jr() / _np.abs(magnet.Jr)
    offset = _np.multiply(Jnorm, magnet.get_size()) / 2

    arrow_tmp = arrow(
        x=(magnet.center[0] - offset[0]),
        y=(magnet.center[1] - offset[1]),
        dx=(2 * offset[0]),
        dy=(2 * offset[1]),
        transform=Affine2D().translate(0, 0),
    )

    magnet_patch_tmp = magnet_patch(patch_tmp, arrow_tmp)
    return magnet_patch_tmp


def _gen_circle_patch(magnet):
    """Generates circcmaxar patches and arrows to represent magnets for 2D plots

    Args:
        magnet (Magnet2D): instance of a magnet class

    Returns:
        magnet_patch: magnet_patch data structure
    """
    # from matplotlib.transforms import Affine2D

    patch_tmp = patch(
        x=(magnet.center[0]),
        y=(magnet.center[1]),
        width=(magnet.radius),
        height=(magnet.radius),
        transform=Affine2D().rotate_deg_around(
            (magnet.center[0]),
            (magnet.center[1]),
            -magnet.alpha + magnet.phi,
        ),
        type="circle",
    )

    Jnorm = magnet.get_Jr() / _np.abs(magnet.Jr)
    offset = magnet.radius * Jnorm[0] / 2
    arrow_tmp = arrow(
        x=(magnet.center[0] - offset),
        y=(magnet.center[1]),
        dx=(2 * offset),
        dy=(2 * 0),
        transform=Affine2D().translate(0, 0),
    )
    magnet_patch_tmp = magnet_patch(patch_tmp, arrow_tmp)
    return magnet_patch_tmp


def _draw_magnets2(ax):
    """Draws Magnets and magnetisation arrows on 2D plots

    Args:
        ax (axis): axis reference
    """
    patch_array = _num_patch_2D()

    # Need original axis transform data before making additional transformations
    axis_transform = ax.transData
    for p in patch_array:
        if p.patch.type == "rectangle":
            sq = _Rect(
                xy=(p.patch.x, p.patch.y),
                width=p.patch.width,
                height=p.patch.height,
                fill=True,
                facecolor="#F5F5F5",
                edgecolor="k",
                zorder=5,
                transform=p.patch.transform + axis_transform,
            )

            ax.add_patch(sq)

        if p.patch.type == "circle":
            circ = _Circ(
                xy=(p.patch.x, p.patch.y),
                radius=p.patch.width,
                fill=True,
                facecolor="#F5F5F5",
                edgecolor="k",
                zorder=5,
                transform=p.patch.transform + axis_transform,
            )

            ax.add_patch(circ)

        ar = _Arrow(
            p.arrow.x,
            p.arrow.y,
            p.arrow.dx,
            p.arrow.dy,
            width=p.arrow.width,
            zorder=6,
            color="k",
            # must translate arrow before rotating
            transform=p.arrow.transform + p.patch.transform + axis_transform,
        )

        ax.add_patch(ar)


def _vector_plot2(points, field, NQ, vector_color):
    """Draws quiver plot of field vectors

    Args:
        points (Point_Array2): coordinates
        field (Field2): Magnetic field
        NQ (int): Plot every NQth arrow in quiver/vector plot
        vector_color (str): Color of arrows
    """
    NPx, NPy = points.x.shape
    if NQ != 0:
        NSx, NSy = NPx // NQ, NPy // NQ
        with _np.errstate(divide="ignore", invalid="ignore"):
            _plt.quiver(
                points.x[::NSx, ::NSy],
                points.y[::NSx, ::NSy],
                field.x[::NSx, ::NSy] / field.n[::NSx, ::NSy],
                field.y[::NSx, ::NSy] / field.n[::NSx, ::NSy],
                color=vector_color,
                alpha=1,
            )


def plot_3D_contour(points, field, plane, **kwargs):
    """Contour plot of field

    Args:
        points (Point_Array2): coordinates
        field (Field2): Magnetic field
        plane (str): Plane to draw contour on. Can be 'xy', 'xz', or 'yz'

    Raises:
        Exception: plot_type must be 'contour' or 'streamplot

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    axis_scale = kwargs.pop("axis_scale", "equal")

    plot_type = kwargs.pop("plot_type", "contour")

    xlab = kwargs.pop("xlab", "x (" + points.unit + ")")
    ylab = kwargs.pop("ylab", "y (" + points.unit + ")")
    zlab = kwargs.pop("zlab", "z (" + points.unit + ")")
    clab = kwargs.pop("clab", "B (" + field.unit + ")")

    SAVE = kwargs.pop("save_fig", False)

    finite_field = field.n[_np.isfinite(field.n)]

    cmax = kwargs.pop("cmax", round(finite_field.mean() * 2, 1))
    num_levels = kwargs.pop("num_levels", 11)

    if plane.lower() == "xy":
        plot_x = points.x
        plot_y = points.y
        plot_xlab = xlab
        plot_ylab = ylab
        stream_x = field.x
        stream_y = field.z

    elif plane.lower() == "xz":
        stream_x = field.x
        stream_y = field.z
        plot_x = points.x
        plot_y = points.z
        plot_xlab = xlab
        plot_ylab = zlab

    else:
        stream_x = field.y
        stream_y = field.z
        plot_x = points.y
        plot_y = points.z
        plot_xlab = ylab
        plot_ylab = zlab

    fig, ax = _plt.subplots(figsize=(8, 8))

    # Generate Contour Plot
    if plot_type.lower() == "contour":
        vector_plot = kwargs.pop("vector_plot", False)
        vector_color = kwargs.pop("vector_color", "w")
        NQ = kwargs.pop("vector_arrows", 11)

        cmap = kwargs.pop("cmap", "viridis")
        cmin = kwargs.pop("cmin", 0)
        lev2 = _np.linspace(cmin, cmax, 256, endpoint=True)
        CS = _plt.contourf(
            plot_x,
            plot_y,
            field.n,
            levels=lev2,
            cmap=_plt.get_cmap(cmap),
            extend="max",
        )

        # Draw contour lines
        if num_levels > 1:
            lev1 = _np.linspace(cmin, cmax, num_levels, endpoint=True)
            _ = _plt.contour(
                plot_x,
                plot_y,
                field.n,
                vmin=cmin,
                vmax=cmax,
                levels=lev1,
                linewidths=1.0,
                colors="k",
            )
            CB = _plt.colorbar(CS, ticks=lev1)

        else:
            CB = _plt.colorbar(CS)

        if vector_plot:
            B_2D = Field2(stream_x, stream_y, unit=field.unit)
            B_2D.n = field.n
            points_2D = Point_Array2(plot_x, plot_y, unit=points.unit)
            _vector_plot2(points_2D, B_2D, NQ, vector_color)

    # Generates streamplot
    elif plot_type.lower() == "streamplot":
        xpl = plot_x[:, 0]
        ypl = plot_y[0, :]
        cmap = kwargs.pop("cmap", None)
        if cmap is not None:
            cmin = kwargs.pop("cmin", -round(finite_field.mean() * 2, 1))
            cmax = kwargs.pop("cmax", round(finite_field.mean() * 2, 1))

            stream_shading = kwargs.pop("stream_shading", "vertical")
            norm = _cm.colors.Normalize(vmin=cmin, vmax=cmax)

            stream_dict = {
                "normal": field.n.T,
                "horizontal": stream_x.T,
                "vertical": stream_y.T,
            }

            CS = _plt.streamplot(
                xpl,
                ypl,
                stream_x.T / field.n.T,
                stream_y.T / field.n.T,
                color=stream_dict.get(stream_shading, "normal"),
                density=1.2,
                norm=norm,
                cmap=cmap,
                linewidth=0.5,
            )
            CB = _plt.colorbar(CS.lines)
        else:
            color = kwargs.pop("color", "k")
            CS = _plt.streamplot(
                xpl,
                ypl,
                stream_x.T / field.n.T,
                stream_y.T / field.n.T,
                density=1.2,
                linewidth=0.5,
                color=color,
            )
            CB = None

    else:
        raise Exception("plot_type must be 'contour' or 'streamplot'")

    if CB is not None:
        CB.ax.get_yaxis().labelpad = 15
        CB.ax.set_ylabel(clab, rotation=270)
    _plt.xlabel(plot_xlab)
    _plt.ylabel(plot_ylab)
    _plt.axis(axis_scale)

    if SAVE:
        _plt.savefig("contour_plot.png", dpi=300)

    return fig, ax


def plot_sub_contour_3D(plot_x, plot_y, plot_B, **kwargs):
    """Contour plot of a single magnetic field component of a 3D simulation

    Args:
        plot_x (ndarray): coordinates for x-axis of plot
        plot_y (ndarray): coordinates for y-axis of plot
        plot_B (ndarray): Magnetic field component to plot

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    cmap = kwargs.pop("cmap", "seismic")
    xlab = kwargs.pop("xlab", f"x (m)")
    ylab = kwargs.pop("ylab", f"y (m)")
    clab = kwargs.pop("clab", f"B (T)")

    axis_scale = kwargs.pop("axis_scale", "equal")

    SAVE = kwargs.pop("save_fig", False)

    cmin = kwargs.pop("cmin", -0.5)
    cmax = kwargs.pop("cmax", 0.5)
    num_levels = kwargs.pop("num_levels", 11)

    lev2 = _np.linspace(cmin, cmax, 256, endpoint=True)
    fig, ax = _plt.subplots(figsize=(8, 8))
    CS = _plt.contourf(
        plot_x, plot_y, plot_B, levels=lev2, cmap=_plt.get_cmap(cmap), extend="both"
    )

    if num_levels > 1:
        lev1 = _np.linspace(cmin, cmax, num_levels, endpoint=True)
        _ = _plt.contour(
            plot_x,
            plot_y,
            plot_B,
            vmin=cmin,
            vmax=cmax,
            levels=lev1,
            linewidths=1.0,
            colors="k",
        )
        CB = _plt.colorbar(CS, ticks=lev1)
    else:
        CB = _plt.colorbar(CS)

    CB.ax.get_yaxis().labelpad = 15
    CB.ax.set_ylabel(clab, rotation=270)
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)
    _plt.axis("equal")
    _plt.show()

    return fig, ax


def line_plot_cylinder(magnet, **kwargs):
    """Calculates and plots the magnetic field along the central axis
    of a cylinder

    This is an example helper function.

    Args:
        magnet (Cylinder): instance of magnetic cylinder

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    rho = _np.linspace(-2 * magnet.radius, 2 * magnet.radius, 51)
    z = _np.array([magnet.length * 1.1 / 2])

    Br, Bz = magnet._calcB_cyl(rho, z)
    fig, ax = _plt.subplots(figsize=(8, 8))
    _plt.plot(rho * 1, Bz, label=r"$B_z$")
    _plt.plot(rho * 1, Br, label=r"$B_r$")
    _plt.legend(loc="best")
    _plt.show()
    return fig, ax


def contour_plot_cylinder(magnet, **kwargs):
    """Calculates and plots the magnetic field
    of a cylinder

    This is an example helper function.


    Args:
        magnet (Cylinder): instance of magnetic cylinder

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    NP = 101
    NPJ = NP * 1j
    rho, z = _np.mgrid[
        -3 * magnet.radius : 3 * magnet.radius : NPJ,
        -magnet.length : magnet.length : NPJ,
    ]
    Br, Bz = magnet._calcB_cyl(rho, z)
    Bn = _np.sqrt(Bz ** 2 + Br ** 2)

    xlab = f"r (m)"
    ylab = "z (m)"

    # plot_B = Bn
    clab = r"$|B|$ (T)"
    cmap = "viridis"
    fig, ax = plot_sub_contour_3D(
        rho * 1,
        z * 1,
        Bn,
        xlab=xlab,
        ylab=ylab,
        clab=clab,
        cmap=cmap,
        cmin=0,
        cmax=1.0,
    )
    return fig, ax
