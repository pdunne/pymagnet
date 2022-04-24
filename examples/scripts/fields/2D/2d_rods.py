"""Calculates and plots the magnetic field due to 4 biaxial rods
"""
import pymagnet as pm
import pymagnet.plots as mplt

pm.reset()
units = "mm"
radius = 10

width = 20
height = 20
hgap_x = width / 2
hgap_y = height / 2


center = (0, -2 * radius)

_ = pm.magnets.Circle(radius=radius, Jr=1.0, center=center, alpha=-90)

center = (0, 2 * radius)
_ = pm.magnets.Circle(radius=radius, Jr=1.0, center=center, alpha=-90)

center = (2 * radius, 0)
_ = pm.magnets.Circle(radius=radius, Jr=1.0, center=center, alpha=90)

center = (-2 * radius, 0)
_ = pm.magnets.Circle(radius=radius, Jr=1.0, center=center, alpha=90)


points = pm.grid2D(4 * radius, 4 * radius)
field = pm.get_field_2D(points)

_, _ = mplt.plot_2D_contour(
    points, field, cmin=0.0, cmax=0.6, num_levels=7, vector_plot=True, vector_arrows=11
)


_, _ = mplt.plot_2D_contour(
    points,
    field,
    cmin=-0.3,
    cmax=0.3,
    cmap="coolwarm",
    plot_type="streamplot",
    stream_color="vertical",
)

_, _ = mplt.plot_2D_contour(points, field, plot_type="streamplot")
