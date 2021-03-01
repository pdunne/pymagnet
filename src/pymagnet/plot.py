# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Plotting routines

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Todo:
    * Reimplement all data fits


.. _Google Python Style Guide:
   https://google.github.io/styleguide/pyguide.html
"""
import matplotlib.pyplot as _plt
from matplotlib.patches import Rectangle as _Rect
from matplotlib.patches import Circle as _Circ
from matplotlib.patches import Arrow as _Arrow
import numpy as _np
from pymagnet import magnets as _mag
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
    """Magnet drawing class

    Args:
        object ([type]): [description]
    """

    def __init__(self, patch, arrow) -> None:
        super().__init__()
        self.patch = patch
        self.arrow = arrow

    def __str__(self) -> str:
        return self.patch.__str__() + self.arrow.__str__()


def plot_1D_field(magnet, **kwargs):

    NP = kwargs.pop("NP", 101)
    return_data = kwargs.pop("return_data", False)

    if issubclass(magnet.__class__, _mag.Cylinder):
        mag_boundary = magnet.length / 2
        z = _np.linspace(
            -2 * magnet.length + magnet.zc, 2 * magnet.length + magnet.zc, NP
        )
        Bz = _mag.magnetic_field_cylinder_1D(magnet, z)

    elif issubclass(magnet.__class__, _mag.Prism):
        mag_boundary = magnet.height / 2
        z = _np.linspace(
            -2 * magnet.height + magnet.zc, 2 * magnet.height + magnet.zc, NP
        )
        Bz = _mag.magnetic_field_prism_1D(magnet, z)

    else:
        print("Error")
        return None

    _, _ = _plt.subplots()
    _plt.xlabel(r"$z$ (mm)")
    _plt.ylabel(r"$B_z$ (mT)")
    _plt.plot(z * 1e3, Bz * 1e3, label="Cube")
    _plt.axvline(x=-mag_boundary * 1e3 + magnet.zc * 1e3, c="blue", ls="--")
    _plt.axvline(x=mag_boundary * 1e3 + magnet.zc * 1e3, c="red", ls="--")
    _plt.axvline(x=0.0, c="k", ls="-")
    _plt.show()
    if return_data:
        return z, Bz
    else:
        return None


def plot_2D_line(x, Field, **kwargs):
    """Line Plot of field from 2D magnet

    Args:
        x ([array]): [assumes in m]
        Field ([field vector_2D]): [field vector object]
    """
    scale_x = kwargs.pop("scale_x", 1e-3)
    scale_cb = kwargs.pop("scale_cb", 1)

    xlab = kwargs.pop("xlab", f"x (mm)")
    ylab = kwargs.pop("ylab", f"B (T)")
    axis_scale = kwargs.pop("axis_scale", "equal")

    SAVE = kwargs.pop("save_fig", False)

    _, ax = _plt.subplots()
    _plt.plot(x / scale_x, Field.n / scale_cb, label=r"$|\mathbf{B}|$")
    _plt.plot(x / scale_x, Field.x / scale_cb, label=r"$B_x$")
    _plt.plot(x / scale_x, Field.y / scale_cb, label=r"$B_y$")
    _plt.legend(loc="best")
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)

    _plt.show()

    if SAVE:
        _plt.savefig("line_plot.png", dpi=300)
        # _plt.savefig('contour_plot.pdf', dpi=300)


def plot_2D_contour(x, y, Field, **kwargs):
    """Contour Plot of field

    Args:
        x ([array]): [assumes in µm]
        y ([array]): [assumes in nm]
        Field ([field vector_2D]): [field vector object]
    """
    scale_x = kwargs.pop("scale_x", 1e-3)
    scale_y = kwargs.pop("scale_y", 1e-3)
    scale_cb = kwargs.pop("scale_cb", 1)

    show_magnets = kwargs.pop("show_magnets", True)

    cmap = kwargs.pop("cmap", "magma")
    xlab = kwargs.pop("xlab", f"x (mm)")
    ylab = kwargs.pop("ylab", f"y (mm)")
    clab = kwargs.pop("clab", f"B (T)")
    axis_scale = kwargs.pop("axis_scale", "equal")

    field_component = kwargs.pop("field_component", "n")

    vector_plot = kwargs.pop("vector_plot", False)
    vector_color = kwargs.pop("vector_color", "w")
    NQ = kwargs.pop("vector_arrows", 11)

    if field_component == "x":
        field_chosen = Field.x
    elif field_component == "y":
        field_chosen = Field.y
    else:
        field_chosen = Field.n

    SAVE = kwargs.pop("save_fig", False)

    finite_field = field_chosen[_np.isfinite(field_chosen)] / scale_cb
    LL = kwargs.pop("LL", 0)
    UL = kwargs.pop("UL", round(finite_field.mean() * 2, 1))
    NL = kwargs.pop("NL", 11)

    lev2 = _np.linspace(LL, UL, 256, endpoint=True)

    _, ax = _plt.subplots()
    CS = _plt.contourf(
        x / scale_x,
        y / scale_y,
        field_chosen / scale_cb,
        levels=lev2,
        cmap=_plt.get_cmap(cmap),
        extend="max",
    )

    # Draw contour lines
    if NL > 1:
        lev1 = _np.linspace(LL, UL, NL, endpoint=True)
        _ = _plt.contour(
            x / scale_x,
            y / scale_y,
            field_chosen / scale_cb,
            vmin=LL,
            vmax=UL,
            levels=lev1,
            linewidths=1.0,
            colors="k",
        )
        CB = _plt.colorbar(CS, ticks=lev1)
    else:
        CB = _plt.colorbar(CS)

    # Draw field vectors
    if vector_plot:
        _vector_plot2(x, y, Field, NQ, scale_x, scale_y, vector_color)

    # Draw magnets and magnetisation arrows
    if show_magnets:
        _draw_magnets2(ax, scale_x, scale_y)

    CB.ax.get_yaxis().labelpad = 15
    CB.ax.set_ylabel(clab, rotation=270)
    _plt.axis(axis_scale)
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)
    _plt.show()

    if SAVE:
        _plt.savefig("contour_plot.png", dpi=300)


def _num_patch_2D(scale_x, scale_y):
    """Generates patches and arrows for drawing

    Returns:
        [tuple(list, list)]: lists of patch and arrow objects
    """
    from .magnets._magnet2 import Magnet_2D, Rectangle, Circle

    patch_array = []
    for magnet in Magnet_2D.instances:
        if issubclass(magnet.__class__, Rectangle):
            magnet_patch_tmp = _gen_rect_patch(magnet, scale_x, scale_y)
        elif issubclass(magnet.__class__, Circle):
            magnet_patch_tmp = _gen_circle_patch(magnet, scale_x, scale_y)
        patch_array.append(magnet_patch_tmp)
    return patch_array


def _gen_rect_patch(magnet, scale_x, scale_y):
    """Generates rectangular patches and arrows to represent magnets for 2D plots

    Args:
        magnet (Magnet object): instance of a magnet class
        scale_x (float): unit scaling
        scale_y (float): unit scaling

    Returns:
        struct: magnet_patch data structure
    """
    from matplotlib.transforms import Affine2D

    patch_tmp = patch(
        x=(magnet.center()[0] - magnet.a) / scale_x,
        y=(magnet.center()[1] - magnet.b) / scale_y,
        width=(2 * magnet.a) / scale_x,
        height=(2 * magnet.b) / scale_y,
        transform=Affine2D().rotate_deg_around(
            (magnet.center()[0]) / scale_x,
            (magnet.center()[1]) / scale_y,
            -magnet.alpha,
        ),
        type="rectangle",
    )

    Jnorm = magnet.get_Jr() / _np.abs(magnet.Jr)
    offset = _np.multiply(Jnorm, magnet.size()) / 2

    arrow_tmp = arrow(
        x=(magnet.center()[0] - offset[0]) / scale_x,
        y=(magnet.center()[1] - offset[1] / scale_y),
        dx=(2 * offset[0]) / scale_x,
        dy=(2 * offset[1]) / scale_y,
        transform=Affine2D().translate(0, magnet.center()[1] / scale_y),
    )

    magnet_patch_tmp = magnet_patch(patch_tmp, arrow_tmp)
    return magnet_patch_tmp


def _gen_circle_patch(magnet, scale_x, scale_y):
    """Generates circular patches and arrows to represent magnets for 2D plots

    Args:
        magnet (Magnet object): instance of a magnet class
        scale_x (float): unit scaling
        scale_y (float): unit scaling

    Returns:
        struct: magnet_patch data structure
    """
    from matplotlib.transforms import Affine2D

    patch_tmp = patch(
        x=(magnet.center()[0]) / scale_x,
        y=(magnet.center()[1]) / scale_y,
        width=(magnet.radius) / scale_x,
        height=(magnet.radius) / scale_y,
        transform=Affine2D().rotate_deg_around(
            (magnet.center()[0]) / scale_x,
            (magnet.center()[1]) / scale_y,
            -magnet.alpha,
        ),
        type="circle",
    )

    Jnorm = magnet.get_Jr() / _np.abs(magnet.Jr)
    offset = (
        magnet.radius
        * Jnorm
        * _np.array([_np.cos(magnet.theta_rad), _np.sin(magnet.theta_rad)])
        / 2
    )
    print(f"offset: ({offset[0]*1e3}, {offset[1]*1e3})")
    print(f"Center{magnet.center()*1e3}")
    arrow_tmp = arrow(
        x=(magnet.center()[0] - offset[0]) / scale_x,
        y=(magnet.center()[1] - offset[1]) / scale_y,
        dx=(2 * offset[0]) / scale_x,
        dy=(2 * offset[1]) / scale_y,
        transform=Affine2D().translate(0, 0),
    )
    print(f"x,y({arrow_tmp.x}, {arrow_tmp.y})")
    print(f"dx,dy({arrow_tmp.dx}, {arrow_tmp.dy})")
    magnet_patch_tmp = magnet_patch(patch_tmp, arrow_tmp)
    return magnet_patch_tmp


def _draw_magnets2(ax, scale_x, scale_y):
    """Draws Magnets and magnetisation arrows on 2D plots

    Args:
        ax (axis?): axis reference
        scale_x (float): unit scaling
        scale_y (float): unit scaling
    """
    patch_array = _num_patch_2D(scale_x, scale_y)

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


def _vector_plot2(x, y, Field, NQ, scale_x, scale_y, vector_color):
    """Draws quiver plot of field vectors

    Args:
        x (float/array): x coordinates
        y (float/array): y coordinates
        Field (Vector2): Field structure
        NQ (int): Plot every NQth vector
        scale_x (float): unit scaling
        scale_y (float): unit scaling
        vector_color (string): quiver color
    """
    NPx, NPy = x.shape
    if NQ != 0:
        NSx, NSy = NPx // NQ, NPy // NQ
        _plt.quiver(
            x[::NSx, ::NSy] / scale_x,
            y[::NSx, ::NSy] / scale_y,
            Field.x[::NSx, ::NSy] / Field.n[::NSx, ::NSy],
            Field.y[::NSx, ::NSy] / Field.n[::NSx, ::NSy],
            color=vector_color,
            alpha=1,
        )


def plot_3D_contour(x, y, z, Field, **kwargs):
    """Contour Plot of field

    Args:
        x ([array]): [assumes in µm]
        y ([array]): [assumes in nm]
        Field ([field vector]): [field vector object]
    """
    scale_x = kwargs.pop("scale_x", 1e-3)
    scale_y = kwargs.pop("scale_y", 1e-3)
    scale_z = kwargs.pop("scale_z", 1e-3)
    scale_cb = kwargs.pop("scale_cb", 1)
    axis_scale = kwargs.pop("axis_scale", "equal")

    patch_array = kwargs.pop("patch_array", _np.array([]))

    cmap = kwargs.pop("cmap", "magma")
    xlab = kwargs.pop("xlab", f"x (mm)")
    ylab = kwargs.pop("ylab", f"y (mm)")
    zlab = kwargs.pop("zlab", f"z (mm)")
    clab = kwargs.pop("clab", f"B (T)")

    SAVE = kwargs.pop("save_fig", False)

    finite_field = Field.n[_np.isfinite(Field.n)] / scale_cb
    LL = kwargs.pop("LL", 0)
    UL = kwargs.pop("UL", round(finite_field.mean() * 2, 1))
    NL = kwargs.pop("NL", 11)

    lev2 = _np.linspace(LL, UL, 256, endpoint=True)

    if _np.ndim(z) < 2:
        orientation = "xy"
        plot_x = x / scale_x
        plot_y = y / scale_y
        plot_xlab = xlab
        plot_ylab = ylab

    elif _np.ndim(y) < 2:
        orientation = "xz"
        plot_x = x / scale_x
        plot_y = z / scale_z
        plot_xlab = xlab
        plot_ylab = zlab

    else:
        orientation = "yz"
        plot_x = y / scale_y
        plot_y = z / scale_z
        plot_xlab = ylab
        plot_ylab = zlab

    print(orientation)
    _, ax = _plt.subplots()
    CS = _plt.contourf(
        plot_x,
        plot_y,
        Field.n / scale_cb,
        levels=lev2,
        cmap=_plt.get_cmap(cmap),
        extend="max",
    )

    if NL > 1:
        lev1 = _np.linspace(LL, UL, NL, endpoint=True)
        _ = _plt.contour(
            plot_x,
            plot_y,
            Field.n / scale_cb,
            vmin=LL,
            vmax=UL,
            levels=lev1,
            linewidths=1.0,
            colors="k",
        )
        CB = _plt.colorbar(CS, ticks=lev1)
    else:
        CB = _plt.colorbar(CS)

    if len(patch_array) > 0:
        for ii in patch_array:
            ii /= scale_x

            sq = _Rect(
                xy=(ii[0], ii[1]),
                width=ii[2],
                height=ii[3],
                fill=True,
                color="#F5F5F5",
                zorder=5,
            )
            ax.add_patch(sq)

            ar = _Arrow(
                ii[0] + ii[2] / 2 - ii[4] / 2,
                ii[1] + ii[3] / 2 - ii[5] / 2,
                ii[4],
                ii[5],
                width=5,
                zorder=6,
                color="k",
            )
            ax.add_patch(ar)

    CB.ax.get_yaxis().labelpad = 15
    CB.ax.set_ylabel(clab, rotation=270)
    _plt.xlabel(plot_xlab)
    _plt.ylabel(plot_ylab)
    _plt.axis(axis_scale)

    if SAVE:
        _plt.savefig("contour_plot.png", dpi=300)


def plot_sub_contour_3D(plot_x, plot_y, plot_B, **kwargs):
    """Contour plot of 3D simulation

    Args:
        plot_x ([type]): [description]
        plot_y ([type]): [description]
        plot_B ([type]): [description]
    """

    cmap = kwargs.pop("cmap", "seismic")
    xlab = kwargs.pop("xlab", f"x (mm)")
    ylab = kwargs.pop("ylab", f"y (mm)")
    clab = kwargs.pop("clab", f"B (T)")

    axis_scale = kwargs.pop("axis_scale", "equal")

    SAVE = kwargs.pop("save_fig", False)

    LL = kwargs.pop("LL", -0.5)
    UL = kwargs.pop("UL", 0.5)
    NL = kwargs.pop("NL", 11)

    lev2 = _np.linspace(LL, UL, 256, endpoint=True)
    _, ax = _plt.subplots()
    CS = _plt.contourf(
        plot_x, plot_y, plot_B, levels=lev2, cmap=_plt.get_cmap(cmap), extend="both"
    )

    if NL > 1:
        lev1 = _np.linspace(LL, UL, NL, endpoint=True)
        _ = _plt.contour(
            plot_x,
            plot_y,
            plot_B,
            vmin=LL,
            vmax=UL,
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


def line_plot_cylinder(magnet, **kwargs):
    """Calculates and plots the magnetic field along the central axis
    of a cylinder

    This is an example helper function.

    Args:
        magnet ([magnet object]):
    """
    rho = _np.linspace(-2 * magnet.radius, 2 * magnet.radius, 51)
    z = _np.array([magnet.length * 1.1 / 2])

    Br, Bz = magnet._calcB_cyl(rho, z)
    _, _ = _plt.subplots()
    _plt.plot(rho * 1e3, Bz, label="Bz")
    _plt.plot(rho * 1e3, Br, label="Br")
    _plt.legend(loc="best")
    _plt.show()


def contour_plot_cylinder(magnet, **kwargs):
    """Calculates and plots the magnetic field
    of a cylinder

    This is an example helper function.


    Args:
        magnet ([type]): [description]
    """
    NP = 101
    NPJ = NP * 1j
    rho, z = _np.mgrid[
        -3 * magnet.radius : 3 * magnet.radius : NPJ,
        -magnet.length : magnet.length : NPJ,
    ]
    Br, Bz = magnet._calcB_cyl(rho, z)
    Bn = _np.sqrt(Bz ** 2 + Br ** 2)

    xlab = f"r (mm)"
    ylab = "z (mm)"

    # plot_B = Bn
    clab = r"$|B|$ (T)"
    cmap = "viridis"
    plot_sub_contour_3D(
        rho * 1e3, z * 1e3, Bn, xlab=xlab, ylab=ylab, clab=clab, cmap=cmap, LL=0, UL=1.0
    )


def param_test_2D(width, height):
    """Example plots while varying the size of a permanent magnet

    This is an example helper function.

    Args:
        width ([type]): [description]
        height ([type]): [description]
    """
    x = _np.linspace(-2 * width, 2 * width, 100)
    y = 1e-3 + height

    B = _mag._routines2.B_calc_2D(x, y)
    plot_2D_line(x, B)
    x, y = _mag._routines2.grid2D(1.5 * width, height)
    B = _mag._routines2.B_calc_2D(x, y)
    cmap = "viridis"
    plot_2D_contour(x, y, B, UL=1, NL=11, cmap=cmap)