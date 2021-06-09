"""Calculates and plots the magnetic field due to 2 rectangular magnets
"""

import pymagnet as pm

pm.reset_magnets()  # clear magnet registry


width = 20
height = 20
unit = "mm"

# Set the space between magnets to be the width of one
half_gap = width / 2

# Center of first magnet
center = (-width / 2 - half_gap, 0)

# Create first magnet
_ = pm.magnets.Rectangle(width=width, height=height, Jr=1.0, center=center, theta=0.0)

# Centre of second magnet
center = (width / 2 + half_gap, 0)

# Create second magnet
_ = pm.magnets.Rectangle(width=width, height=height, Jr=1.0, center=center, theta=90.0)

# Prepare 100x100 grid of x,y coordinates to calculate the field
points = pm.grid2D(2 * width, 2 * height, unit=unit)

# Calculate the magnetic field due to all magnets in the registry
field = pm.get_field_2D(points)

# Plot the result, vector_plot = True toggles on the vector field plot
_, _ = pm.plots.plot_2D_contour(
    points,
    field,
    cmin=0.0,  # minimum field value
    cmax=0.5,  # maximum field value
    vector_plot=True,  # plot the vector field
)

_, _ = pm.plots.plot_2D_contour(
    points,
    field,
    cmin=-0.1,  # minimum field value
    cmax=0.1,  # maximum field value
    cmap="coolwarm",  # set the colormap
    plot_type="streamplot",
    show_magnets=True,
)
_, _ = pm.plots.plot_2D_contour(
    points, field, plot_type="streamplot", show_magnets=True
)
