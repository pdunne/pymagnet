# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Example plots
"""
def main():
    print("Line Plot through centre of a cylinder")
    
    import pymagnet as pm
    R = 5e-3
    L = 25e-3
    m_cyl = pm.magnets.Cylinder(radius=R, length=L, Jr=1.2,
                                center=(0.0, 0.0, 0))
    pm.plot.plot_1D_field(m_cyl);

    print("Contour Plot of 2D Magnets")
    pm.reset_magnets()
    cmap = 'viridis' # colormap
    width = 20e-3
    height = 20e-3
    hgap_x = width / 2

    center = (-width / 2 - hgap_x, 0)
    _ = pm.magnets.Rectangle(width=width, height=height,
                            Jr=1.0, center=center, theta=0.0)
    center = (width / 2 + hgap_x, 0)
    _ = pm.magnets.Rectangle(width=width, height=height,
                            Jr=1.0, center=center, theta=90.0)

    x, y = pm.grid2D(2 * width, 2 * height)
    B = pm.B_calc_2D(x, y)


    pm.plot.plot_2D_contour(x, y, B, UL=.6, NL=7, vector_plot=True, cmap=cmap)


if __name__ == '__main__':
    main()
