"""Calculates and plots the magnetic field due to rectangular magnet arrays."""

import matplotlib.pyplot as plt
import numpy as np

import pymagnet as pm
import pymagnet.plots as mplt

WIDTH = 6.0
HEIGHT = 6.0
UNIT = "mm"
JR = 1.26
NUM_POINTS = 201
Y_EPSILON = 0.001  # slight offset to avoid perfect grid alignment issues
SAMPLE_HEIGHT = 1.0  # y-coordinate (mm) at which to sample the line field
VGAP = 6.0  # vertical gap between the two Halbach arrays
mag_offset = 4


def gen_alternating_magnets(N=5, width=WIDTH, height=HEIGHT, Jr=JR):
    """Generate N rectangular magnets with alternating remanence (+/-Jr),
    keeping phi=90.0 for all magnets, centered along x."""
    x0 = -(N - 1) * width / 2
    y0 = -height / 2 - Y_EPSILON

    for i in range(N):
        pm.magnets.Rectangle(
            width=width, height=height, Jr=Jr * (-1) ** i, center=(x0, y0), phi=90.0
        )
        x0 += width


def gen_lin_halbach(
    N=5,
    width=WIDTH,
    height=HEIGHT,
    Jr=JR,
    y_offset=0.0,
    alpha=0.0,
    flip=False,
    mag_offset=0,
):
    """Generate a linear Halbach array of N magnets (N >= 1), centered along x.

    The base pattern cycles: RIGHT, UP, LEFT, DOWN.
    When flip=True the vertical components are swapped (UP<->DOWN),
    producing the mirrored bottom array.
    """
    x0 = -(N - 1) * width / 2
    y0 = -height / 2 - Y_EPSILON + y_offset

    # Base pattern: (phi, Jr_sign) indexed by i % 4
    #   RIGHT(+Jr, 0), UP(+Jr, 90), LEFT(-Jr, 0), DOWN(-Jr, 90)
    pattern = [
        (0.0, 1),  # RIGHT
        (90.0, 1),  # UP
        (0.0, -1),  # LEFT
        (90.0, -1),  # DOWN
    ]

    for i in range(N):
        i -= mag_offset
        phi, sign = pattern[i % 4]
        if flip and phi == 90.0:
            sign = -sign

        pm.magnets.Rectangle(
            width=width,
            height=height,
            Jr=Jr * sign,
            center=(x0, y0),
            phi=phi,
            alpha=alpha,
        )
        x0 += width


def main():
    N = 10
    xmax = (N + 1) * WIDTH / 2

    # --- Single Halbach: line plot ---
    # pm.reset()
    # gen_lin_halbach(N=N, width=WIDTH, height=HEIGHT, Jr=JR)

    #     x = np.linspace(-xmax, xmax, NUM_POINTS)
    #     y = np.full_like(x, SAMPLE_HEIGHT)
    #     line_points = pm.utils.Point_Array2(x, y)
    #
    #     field_hal = pm.get_field_2D(line_points)
    #     fig, ax = pm.plots.plot_2D_line(line_points, field_hal)
    #     fig.savefig(f"halbach_{N:02d}.pdf", dpi=200)
    #     plt.close(fig)

    # for mag_offset in range(1):
    # pm.reset()
    # gen_lin_halbach(N=N, width=WIDTH, height=HEIGHT, Jr=JR, mag_offset=mag_offset)

    # field_hal = pm.get_field_2D(line_points)
    # fig, ax = pm.plots.plot_2D_line(line_points, field_hal)
    # fig.savefig(f"halbach_{N:02d}_offset{mag_offset}.pdf", dpi=200)
    # plt.close(fig)
    # --- Double Halbach: contour + line plot ---
    pm.reset()
    gen_lin_halbach(
        N=N,
        width=WIDTH,
        height=HEIGHT,
        Jr=JR,
        y_offset=HEIGHT + VGAP / 2,
        flip=True,
        mag_offset=0,
    )
    gen_lin_halbach(
        N=N,
        width=WIDTH,
        height=HEIGHT,
        Jr=JR,
        y_offset=-VGAP / 2,
        flip=False,
        mag_offset=mag_offset,
    )

    gen_lin_halbach(
        N=N,
        width=WIDTH,
        height=HEIGHT,
        Jr=JR,
        y_offset=2 * HEIGHT + VGAP / 2,
        flip=False,
        mag_offset=0,
    )
    gen_lin_halbach(
        N=N,
        width=WIDTH,
        height=HEIGHT,
        Jr=JR,
        y_offset=-HEIGHT - VGAP / 2,
        flip=True,
        mag_offset=mag_offset,
    )

    grid_points = pm.grid2D(
        xmax=xmax, ymax=32.0, unit=UNIT, ymin=-32, num_points=NUM_POINTS
    )
    field_grid = pm.get_field_2D(grid_points)

    fig, _ = mplt.plot_2D_contour(
        grid_points,
        field_grid,
        cmin=0.0,
        cmax=1.0,
        vector_plot=True,
        show_magnets=True,
    )
    fig.savefig(f"double_halbach_{N:02d}_offset{mag_offset:01d}_contour.pdf", dpi=200)
    plt.close(fig)

    # B·∇B tensor force contour

    fig, ax = mplt.plot_2D_contour_BdotgradB(
        grid_points, field_grid, cmin=0.0, cmax=0.5
    )
    fig.savefig(f"double_halbach_{N:02d}_offset{mag_offset:01d}_BdotgradB.pdf", dpi=200)
    plt.close(fig)

    # Jacobian 2x2 panel (all 4 partial derivatives)
    fig, axes = mplt.plot_2D_contour_jacobian(
        grid_points, field_grid, cmap="RdBu_r", num_levels=15
    )
    fig.savefig(f"double_halbach_{N:02d}_offset{mag_offset:01d}_jacobian.pdf", dpi=200)
    plt.close(fig)

    fig, ax = mplt.plot_2D_contour_gradB(
        grid_points, field_grid, cmax=0.5, show_magnets=True
    )


#
#         line_points = pm.utils.Point_Array2(x, y)
#         field_double_hal = pm.get_field_2D(line_points)
#         fig, ax = pm.plots.plot_2D_line(line_points, field_double_hal)
#         fig.savefig(f"double_halbach_{N:02d}_offset{mag_offset:01d}.pdf", dpi=200)
#         plt.close(fig)

# --- Alternating magnets: line plot ---
# pm.reset()
# gen_alternating_magnets(N=N, width=WIDTH, height=HEIGHT, Jr=JR)

# line_points = pm.utils.Point_Array2(x, y)
# field_alt = pm.get_field_2D(line_points)
# fig, ax = pm.plots.plot_2D_line(line_points, field_alt)
# fig.savefig(f"alternating_{N:02d}.pdf", dpi=200)
# plt.close(fig)

# # --- Comparison plot ---
# fig, ax = plt.subplots(figsize=(8, 8))
# ax.set_title(f"Field along line y={SAMPLE_HEIGHT:.1f} {UNIT} for N={N} magnets")
# ax.plot(x, field_alt.n, label="Alternating")
# ax.plot(x, field_hal.n, label="Linear Halbach")
# ax.plot(x, field_double_hal.n, label="Double Halbach")
# ax.legend(loc="best")
# ax.set_xlabel(f"x ({UNIT})")
# ax.set_ylabel(r"$|\mathbf{B}|$")
# fig.savefig(f"comparison_{N:02d}_offset{mag_offset:01d}.pdf", dpi=200)
# plt.show()
# plt.close(fig)


if __name__ == "__main__":
    main()
