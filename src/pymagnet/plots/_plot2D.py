# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Plotting routines

This module contains all functions needed to plot lines and contours for 2D
magnetic sources, and

"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

try:
    import matplotlib.pyplot as _plt

except ImportError:
    _has_matplotlib = False
    warnings.warn("matplotlib is not installed", UserWarning, stacklevel=2)

else:
    _has_matplotlib = True
    import matplotlib.cm as _cm
    from matplotlib.patches import Arrow as _Arrow
    from matplotlib.patches import Circle as _Circ
    from matplotlib.patches import Rectangle as _Rect
    from matplotlib.transforms import Affine2D

import numpy as _np
from numpy.typing import NDArray

from ..utils import Field2, Field3, Point_Array2, Point_Array3

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

# from ..utils._conversions import get_unit_value_meter, get_unit_value_tesla
# from .. import magnets as _mag


class patch:
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


class arrow:
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
        if not _has_matplotlib:
            raise ImportError("matplotlib is required to use this plot function.")

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


class magnet_patch:
    """Magnet drawing class"""

    def __init__(self, patch, arrow) -> None:
        super().__init__()
        self.patch = patch
        self.arrow = arrow

    def __str__(self) -> str:
        return self.patch.__str__() + self.arrow.__str__()


def plot_2D_line(
    point_array: Point_Array2, field: Field2, **kwargs: Any
) -> tuple[Figure, Axes]:
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
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

    xlab = kwargs.pop("xlab", f"x ({point_array.unit})")
    ylab = kwargs.pop("ylab", f"B ({field.unit})")
    # axis_scale = kwargs.pop("axis_scale", "equal")

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


def plot_2D_contour(
    point_array: Point_Array2, field: Field2, **kwargs: Any
) -> tuple[Figure, Axes]:
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
        num_arrows (int or None): Number of arrows per axis. Defaults to None.
        vector_color (str): Arrow color. Defaults to 'w'
        cmin (float): Color scale minimum. Defaults to 0.0
        cmax (float): Color scale minimum. Defaults to twice the mean field
        num_levels (int): Number of contour levels. Defaults to 11.

    Raises:
        Exception: plot_type must be 'contour' or 'streamplot'

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

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
        vector_color = kwargs.pop("vector_color", "w")
        NQ = kwargs.pop("num_arrows", None)

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
        CS.set_edgecolor("face")

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
        if NQ is not None:
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


def plot_2D_contour_gradient(
    point_array: Point_Array2, field: Field2, **kwargs: Any
) -> tuple[Figure, Axes]:
    """Contour plot of the magnetic field gradient magnitude.

    Computes grad(|B|) and plots as a filled contour. Delegates to
    plot_2D_contour with appropriate default labels.

    Args:
        point_array (Point_Array2): coordinates
        field (Field2): Magnetic field (must have .n populated via calc_norm)

    Kwargs:
        All kwargs from plot_2D_contour are supported.
        clab defaults to "∇|B| (T/{unit})"

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    from ..utils._routines2D import gradB_2D

    grad_field = gradB_2D(field.n, point_array.x, point_array.y)
    kwargs.setdefault("clab", f"\u2207|B| (T/{point_array.unit})")
    if "cmax" not in kwargs:
        finite = grad_field.n[_np.isfinite(grad_field.n)]
        kwargs["cmax"] = float(_np.max(finite)) if finite.size > 0 else 1.0
    return plot_2D_contour(point_array, grad_field, **kwargs)


def plot_2D_contour_force(
    point_array: Point_Array2,
    field: Field2,
    chi_m: float,
    c: float,
    **kwargs: Any,
) -> tuple[Figure, Axes]:
    """Contour plot of the scalar magnetic gradient force.

    Computes F = (chi_m / mu_0) * c * |B| * grad(|B|) and plots as a
    filled contour. Delegates to plot_2D_contour with appropriate labels.

    Args:
        point_array (Point_Array2): coordinates
        field (Field2): Magnetic field
        chi_m (float): Magnetic susceptibility
        c (float): Material constant

    Kwargs:
        All kwargs from plot_2D_contour are supported.
        clab defaults to "F∇B (T²/{unit})"

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    from ..utils._routines2D import FgradB_2D

    force_field = FgradB_2D(field, point_array.x, point_array.y, chi_m, c)
    kwargs.setdefault("clab", f"F\u2207B (T\u00b2/{point_array.unit})")
    if "cmax" not in kwargs:
        finite = force_field.n[_np.isfinite(force_field.n)]
        kwargs["cmax"] = float(_np.max(finite)) if finite.size > 0 else 1.0
    return plot_2D_contour(point_array, force_field, **kwargs)


def plot_2D_contour_BdotgradB(
    point_array: Point_Array2, field: Field2, **kwargs: Any
) -> tuple[Figure, Axes]:
    """Contour plot of B . grad(B) tensor force.

    Computes the exact tensor product F_i = sum_j B_j * dB_j/dx_i
    and plots as a filled contour. Delegates to plot_2D_contour.

    Args:
        point_array (Point_Array2): coordinates
        field (Field2): Magnetic field

    Kwargs:
        All kwargs from plot_2D_contour are supported.
        clab defaults to "B·∇B (T²/{unit})"

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    from ..utils._routines2D import BdotgradB_2D

    bdg = BdotgradB_2D(field, point_array.x, point_array.y)
    kwargs.setdefault("clab", f"B\u00b7\u2207B (T\u00b2/{point_array.unit})")
    if "cmax" not in kwargs:
        finite = bdg.n[_np.isfinite(bdg.n)]
        kwargs["cmax"] = float(_np.max(finite)) if finite.size > 0 else 1.0
    return plot_2D_contour(point_array, bdg, **kwargs)


def plot_2D_contour_gradB2(
    point_array: Point_Array2, field: Field2, **kwargs: Any
) -> tuple[Figure, Axes]:
    """Contour plot of |∇(B²)| with ∇(B²) vector overlay.

    Computes ∇(|B|²) directly by taking the gradient of field.n².
    The contour shows |∇(B²)| and quiver arrows show the ∇(B²) direction.

    Args:
        point_array (Point_Array2): coordinates
        field (Field2): Magnetic field

    Kwargs:
        All kwargs from plot_2D_contour are supported.
        clab defaults to "|∇(B²)| (T²/{unit})"
        num_arrows (int): Number of quiver arrows per axis side. Defaults to 10.

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    from ..utils._routines2D import gradB_2D

    # Compute ∇(|B|²) directly from the gradient of the squared magnitude
    grad_B2 = gradB_2D(field.n**2, point_array.x, point_array.y)

    kwargs.setdefault("clab", f"|\u2207(B\u00b2)| (T\u00b2/{point_array.unit})")
    kwargs.setdefault("num_arrows", 10)
    if "cmax" not in kwargs:
        finite = grad_B2.n[_np.isfinite(grad_B2.n)]
        kwargs["cmax"] = float(_np.max(finite)) if finite.size > 0 else 1.0
    return plot_2D_contour(point_array, grad_B2, **kwargs)


def plot_2D_contour_gradB(
    point_array: Point_Array2, field: Field2, **kwargs: Any
) -> tuple[Figure, Axes]:
    """Contour plot of the Jacobian Frobenius norm with grad(|B|) vector overlay.

    Computes the full Jacobian J_ij = dB_i/dx_j via jacobian_B_2D, then plots
    ||J||_F = sqrt(sum of J_ij^2) as a filled contour with quiver arrows
    showing the direction of grad(|B|) overlaid.

    Args:
        point_array (Point_Array2): coordinates
        field (Field2): Magnetic field

    Kwargs:
        All kwargs from plot_2D_contour are supported.
        clab defaults to "||∇B|| (T/{unit})"
        num_arrows (int): Number of quiver arrows per axis side. Defaults to 10.

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    from ..utils._routines2D import gradB_2D, jacobian_B_2D

    J = jacobian_B_2D(field, point_array.x, point_array.y)

    # Frobenius norm of the Jacobian: total rate of field variation
    frob = _np.sqrt(
        J.dBx_dx**2 + J.dBx_dy**2 + J.dBy_dx**2 + J.dBy_dy**2
    )

    # Build a Field2 with grad(|B|) as the vector components and ||J||_F as norm
    grad_field = gradB_2D(field.n, point_array.x, point_array.y)
    grad_field.n = frob

    kwargs.setdefault("clab", f"||\u2207B|| (T/{point_array.unit})")
    kwargs.setdefault("num_arrows", 10)
    if "cmax" not in kwargs:
        finite = frob[_np.isfinite(frob)]
        kwargs["cmax"] = float(_np.max(finite)) if finite.size > 0 else 1.0
    return plot_2D_contour(point_array, grad_field, **kwargs)


def plot_2D_contour_jacobian(
    point_array: Point_Array2, field: Field2, **kwargs: Any
) -> tuple[Figure, NDArray]:
    """2x2 contour plot of all Jacobian tensor components.

    Computes the full 2D Jacobian J_ij = dB_i/dx_j and plots each of the
    4 components (dBx/dx, dBx/dy, dBy/dx, dBy/dy) as a subplot panel.

    Args:
        point_array (Point_Array2): coordinates
        field (Field2): Magnetic field

    Kwargs:
        cmap (str): Colormap. Defaults to 'RdBu_r' (diverging).
        cmin (float): Color scale minimum. Defaults to -cmax (symmetric).
        cmax (float): Color scale maximum. Defaults to max(|data|).
        num_levels (int): Number of contour levels. Defaults to 11.
        show_magnets (bool): Draw magnets. Defaults to True.
        save_fig (bool): Save to png file. Defaults to False.

    Returns:
        tuple: fig, axes (2x2 array of matplotlib axes)
    """
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

    from ..magnets._polygon2D import PolyMagnet
    from ..utils._routines2D import jacobian_B_2D

    J = jacobian_B_2D(field, point_array.x, point_array.y)

    cmap = kwargs.pop("cmap", "RdBu_r")
    num_levels = kwargs.pop("num_levels", 11)
    show_magnets = kwargs.pop("show_magnets", True)
    SAVE = kwargs.pop("save_fig", False)
    unit = point_array.unit

    components = [
        (J.dBx_dx, f"\u2202Bx/\u2202x (T/{unit})"),
        (J.dBx_dy, f"\u2202Bx/\u2202y (T/{unit})"),
        (J.dBy_dx, f"\u2202By/\u2202x (T/{unit})"),
        (J.dBy_dy, f"\u2202By/\u2202y (T/{unit})"),
    ]

    # Determine symmetric color limits from all components
    all_finite = _np.concatenate(
        [comp[_np.isfinite(comp)].ravel() for comp, _ in components]
    )
    default_cmax = float(_np.max(_np.abs(all_finite))) if all_finite.size > 0 else 1.0
    cmax = kwargs.pop("cmax", round(default_cmax, 2))
    cmin = kwargs.pop("cmin", -cmax)

    fig, axes = _plt.subplots(2, 2, figsize=(14, 12))
    lev_fill = _np.linspace(cmin, cmax, 256, endpoint=True)
    lev_lines = _np.linspace(cmin, cmax, num_levels, endpoint=True)

    for ax, (data, label) in zip(axes.ravel(), components):
        CS = ax.contourf(
            point_array.x,
            point_array.y,
            data,
            levels=lev_fill,
            cmap=_plt.get_cmap(cmap),
            extend="both",
        )
        CS.set_edgecolor("face")

        if num_levels > 1:
            ax.contour(
                point_array.x,
                point_array.y,
                data,
                levels=lev_lines,
                linewidths=0.5,
                colors="k",
            )
            CB = fig.colorbar(CS, ax=ax, ticks=lev_lines)
        else:
            CB = fig.colorbar(CS, ax=ax)

        CB.ax.get_yaxis().labelpad = 15
        CB.ax.set_ylabel(label, rotation=270)
        ax.set_xlabel(f"x ({unit})")
        ax.set_ylabel(f"y ({unit})")
        ax.set_aspect("equal")

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

    fig.tight_layout()
    if SAVE:
        _plt.savefig("jacobian_plot.png", dpi=300)

    return fig, axes


def _num_patch_2D():
    """Generates patches and arrows for drawing

    Returns:
        tuple: (list, list) lists of patch and arrow objects
    """
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

    from ..magnets._magnet2D import Circle, Magnet2D, Rectangle

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
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

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
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

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
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

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


def plot_3D_contour(
    points: Point_Array2 | Point_Array3,
    field: Field2 | Field3,
    plane: str,
    **kwargs: Any,
) -> tuple[Figure, Axes]:
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
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

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
        vector_color = kwargs.pop("vector_color", "w")
        NQ = kwargs.pop("num_arrows", None)

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

        if NQ is not None:
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
            CS.set_edgecolor("face")
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


def plot_sub_contour_3D(
    plot_x: NDArray[_np.floating],
    plot_y: NDArray[_np.floating],
    plot_B: NDArray[_np.floating],
    **kwargs: Any,
) -> tuple[Figure, Axes]:
    """Contour plot of a single magnetic field component of a 3D simulation

    Args:
        plot_x (ndarray): coordinates for x-axis of plot
        plot_y (ndarray): coordinates for y-axis of plot
        plot_B (ndarray): Magnetic field component to plot

    Returns:
        tuple: fig, ax reference to matplotlib figure and axis objects
    """
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

    cmap = kwargs.pop("cmap", "seismic")
    xlab = kwargs.pop("xlab", "x (m)")
    ylab = kwargs.pop("ylab", "y (m)")
    clab = kwargs.pop("clab", "B (T)")

    # axis_scale = kwargs.pop("axis_scale", "equal")

    # SAVE = kwargs.pop("save_fig", False)

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
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

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
    if not _has_matplotlib:
        raise ImportError("matplotlib is required to use this plot function.")

    NP = 101
    NPJ = NP * 1j
    rho, z = _np.mgrid[
        -3 * magnet.radius : 3 * magnet.radius : NPJ,
        -magnet.length : magnet.length : NPJ,
    ]
    Br, Bz = magnet._calcB_cyl(rho, z)
    Bn = _np.sqrt(Bz**2 + Br**2)

    xlab = "r (m)"
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
