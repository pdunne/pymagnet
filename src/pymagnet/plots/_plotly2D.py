# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""2D plotting routines using Plotly as the backend: NOT IMPLEMENTED YET

This module contains all functions needed to plot lines and contours for 2D
magnetic sources.

"""
import plotly.graph_objects as _go
import plotly.figure_factory as _ff
import numpy as _np
from ..utils._conversions import get_unit_value_meter


def _plotly_vector_plot2(x, y, Field, NQ, scale_x, scale_y, vector_color):
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
    plot_object = []
    NPx, NPy = x.shape
    if NQ != 0:
        NSx, NSy = NPx // NQ, NPy // NQ
        with _np.errstate(divide="ignore", invalid="ignore"):
            plot_object = _ff.create_quiver(
                x[::NSx, ::NSy] / scale_x,
                y[::NSx, ::NSy] / scale_y,
                Field.x[::NSx, ::NSy] / Field.n[::NSx, ::NSy],
                Field.y[::NSx, ::NSy] / Field.n[::NSx, ::NSy],
                scale=1.5,
                arrow_scale=0.5,  # name='quiver',
                line_width=1,
                line=dict(color="white"),
            )

    return plot_object


def plotly_2D_contour(x, y, Field, **kwargs):
    """Contour Plot of field

    Args:
        x ([array]): [assumes in µm]
        y ([array]): [assumes in nm]
        Field ([field vector_2D]): [field vector object]
    """

    scale_x = kwargs.pop("scale_x", 1)
    scale_y = kwargs.pop("scale_y", 1)
    scale_cb = kwargs.pop("scale_cb", 1)

    colorscale = kwargs.pop("colorscale", "viridis")
    cmin = kwargs.pop("cmin", 0.0)
    cmax = kwargs.pop("cmax", round(_np.nanmean(Field.n[:]) * 2, 1))
    cstep = kwargs.pop("cstep", 0.1)

    show_magnets = kwargs.pop("show_magnets", True)
    title = kwargs.pop("title", None)
    xlab = kwargs.pop("xlab", "x (mm)")
    ylab = kwargs.pop("ylab", "y (mm)")
    clab = kwargs.pop("clab", "B (T)")

    # field_component = kwargs.pop("field_component", "n")

    # plot_type = kwargs.pop("plot_type", "contour")

    plot_objects = []
    # cache = {}

    plot_objects.append(
        _plotly_contour2(
            x,
            y,
            Field.n,
            cmin=cmin,
            cmax=cmax,
            cstep=cstep,
            colorscale=colorscale,
        )
    )

    plot_objects.append(_plotly_draw_circle(p1=(-10, -10), p2=(10, 10)))
    plot_objects.append(_plotly_draw_arrow(head=(6, -6), tail=(-6, 6)))

    fig = _go.Figure(data=plot_objects)

    fig.update_layout(
        title=title,
        xaxis_title=xlab,
        yaxis_title=ylab,
    )
    fig.update_yaxes(
        scaleanchor="x",
        scaleratio=1,
    )

    fig.update_xaxes(range=[x.min() * 1, x.max() * 1])
    fig.update_yaxes(range=[y.min() * 1, y.max() * 1])
    # fig.update_yaxes(automargin=True)

    fig.update_layout(
        autosize=False,
        width=500,
        height=500,
        margin=dict(l=50, r=50, b=50, t=100, pad=0),
    )

    fig.show()


def _plotly_contour2(x, y, z, **kwargs):
    colorscale = kwargs.pop("colorscale", "viridis")
    cmin = kwargs.pop("cmin", 0.0)
    cmax = kwargs.pop("cmax", round(z.mean() * 2, 1))
    cstep = kwargs.pop("cstep", 0.1)

    plot_object = [
        _go.Contour(
            z=z,
            x=x,  # horizontal axis
            y=y,  # vertical axis
            colorscale=colorscale,
            connectgaps=True,
            line_smoothing=0.85,
            zmin=cmin,
            zmax=cmax,
            contours=dict(
                # coloring="heatmap",
                showlabels=True,  # show labels on contours
                start=cmin,
                end=cmax,
                size=cstep,
                labelfont=dict(  # label font properties
                    size=12,
                    color="white",
                ),
            ),
            colorbar=dict(
                title="|B| (T)",
                titleside="right",
                titlefont=dict(size=14, family="Arial, sans-serif"),
            ),
        )
    ]
    return plot_object


def _plotly_draw_circle(p1=(-0.5, -0.5), p2=(0.5, 0.5), **kwargs):
    fillcolor = kwargs.pop("fillcolor", "white")
    line_color = kwargs.pop("line_color", "black")
    return dict(
        type="circle",
        xref="x",
        yref="y",
        fillcolor=fillcolor,
        x0=p1[0],
        y0=p1[1],
        x1=p2[0],
        y1=p2[1],
        line_color=line_color,
    )


def _plotly_draw_arrow(head=(-0.5, -0.5), tail=(0.5, 0.5), **kwargs):

    arrowcolor = kwargs.pop("arrowcolor", "black")
    return dict(
        x=head[0],  # arrows' head
        y=head[1],  # arrows' head
        ax=tail[0],  # arrows' tail
        ay=tail[1],  # arrows' tail
        xref="x",
        yref="y",
        axref="x",
        ayref="y",
        text="",  # if you want only the arrow
        showarrow=True,
        arrowhead=3,
        arrowsize=1,
        arrowwidth=5,
        arrowcolor=arrowcolor,
    )
