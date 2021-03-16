import sys


import pymagnet as pm
import matplotlib.pyplot as plt
import numpy as np


# def plot_test(theta=0.0, phi=0.0):
#     pm.reset_magnets()
#     cmap = "viridis"

#     radius = 10e-3
#     Jr = 1.0
#     xc = radius
#     yc = radius / 2
#     zc = -radius / 2
#     center = (xc, yc, zc)
#     center2 = (-xc, -yc, -zc)
#     NN = 4
#     Jx, Jy, Jz = pm.magnets._routines.sph2cart(
#         Jr,
#         np.deg2rad(theta),
#         np.deg2rad(phi),
#     )

#     print(np.around([Jx, Jy, Jz], 3))

#     s1 = pm.magnets.Sphere(radius=radius, Jr=1.0, theta=theta, phi=phi, center=center)
#     #     s2 = pm.magnets.Sphere(radius=radius, Jr=-1.0, theta = theta, phi = phi, center=center2)

#     x, z = pm.grid2D(NN * radius, NN * radius)
#     y = np.array([0])
#     B = pm.B_calc_3D(x, y, z)
#     B.calc_norm()

#     ax1 = pm.plot.plot_3D_contour(x, y, z, B, UL=0.6, NL=13, cmap=cmap)
#     c = plt.Circle(
#         xy=(xc * 1e3, zc * 1e3),
#         radius=radius * 1e3,
#         edgecolor="k",
#         facecolor="white",
#         zorder=5,
#     )
#     c2 = plt.Circle(
#         xy=(-xc * 1e3, -zc * 1e3),
#         radius=radius * 1e3,
#         edgecolor="k",
#         facecolor="white",
#         zorder=5,
#     )
#     ax1.add_patch(c)
#     #     ax1.add_patch(c2)

#     B2 = pm.magnets._fields.Vector2(B.x, B.z)
#     B2.calc_norm()
#     pm.plot._vector_plot2(x, z, B2, NQ=11, scale_x=1e-3, scale_y=1e-3, vector_color="w")

#     x, y = pm.grid2D(NN * radius, NN * radius)
#     z = np.array([0])
#     B = pm.B_calc_3D(x, y, z)
#     B.calc_norm()

#     ax1 = pm.plot.plot_3D_contour(x, y, z, B, UL=0.6, NL=13, cmap=cmap)
#     c = plt.Circle(
#         xy=(xc * 1e3, yc * 1e3),
#         radius=radius * 1e3,
#         edgecolor="k",
#         facecolor="white",
#         zorder=5,
#     )
#     c2 = plt.Circle(
#         xy=(-xc * 1e3, -yc * 1e3),
#         radius=radius * 1e3,
#         edgecolor="k",
#         facecolor="white",
#         zorder=5,
#     )
#     ax1.add_patch(c)
#     #     ax1.add_patch(c2)

#     B2 = pm.magnets._fields.Vector2(B.x, B.y)
#     B2.calc_norm()
#     pm.plot._vector_plot2(x, y, B2, NQ=11, scale_x=1e-3, scale_y=1e-3, vector_color="w")

#     y, z = pm.grid2D(NN * radius, NN * radius)
#     x = np.array([0])
#     B = pm.B_calc_3D(x, y, z)
#     B.calc_norm()

#     ax1 = pm.plot.plot_3D_contour(x, y, z, B, UL=0.6, NL=13, cmap=cmap)
#     c = plt.Circle(
#         xy=(yc * 1e3, zc * 1e3),
#         radius=radius * 1e3,
#         edgecolor="k",
#         facecolor="white",
#         zorder=5,
#     )
#     c2 = plt.Circle(
#         xy=(-yc * 1e3, -zc * 1e3),
#         radius=radius * 1e3,
#         edgecolor="k",
#         facecolor="white",
#         zorder=5,
#     )
#     ax1.add_patch(c)
#     #     ax1.add_patch(c2)

#     B2 = pm.magnets._fields.Vector2(B.y, B.z)
#     B2.calc_norm()
#     pm.plot._vector_plot2(y, z, B2, NQ=11, scale_x=1e-3, scale_y=1e-3, vector_color="w")


# plot_test(theta=0.0, phi=90.0)


pm.reset_magnets()
a = 20e-3
b = 40e-3
c = 60e-3
center = (0, 0, 0)

# magnetised in x:
# theta, phi = 0.0, 90.0

# magnetised in y:
# theta, phi = 90.0, 90.0


# magnetised in z:
theta, phi = 90.0, 0.0
# or
# theta, phi= 0.0, 0.0
cmap = "viridis"

magnet = pm.magnets.Prism(
    width=a,
    depth=b,
    height=c,
    Jr=1.0,
    center=center,
    theta=theta,
    phi=phi,
    alpha=0,
    beta=45,
    gamma=0,
)
x, y = pm.grid2D(c / 2, c / 2)
z = (c / 2) * 1.001

B = pm.B_calc_3D(x, y, z)

_ = pm.plots.plot_3D_contour(x, y, z, B, UL=0.6, NL=13, cmap=cmap)
