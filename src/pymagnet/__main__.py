def main():
    print("Example Plot")
    R = 5e-3
    L = 25e-3
    # import magnets
    # import plot
    import pymagnet as mag
    # from . import magnets
    # from . import plot
    m_cyl = mag.magnets.Cylinder(R=R, L=L, Jr=1.2,
                                        center=(0.0, 0.0, 0))

    mag.plot.plot_1D_field(m_cyl);


if __name__ == '__main__':
    main()
