import matplotlib.pyplot as plt
import numpy as np

import pymagnet as pm

PI = np.pi
MU0 = 4e-7 * PI


def plot_force_result(offsets, f_total, spacing=None, compare=None, plot_type='force'):
    _, _ = plt.subplots(figsize=(8, 8))
    
    if plot_type.lower() == 'torque':
        labels = [r"$\tau_x$", r"$\tau_y$", r"$\tau_z$"]
    else:
        labels = [r"$F_x$", r"$F_y$", r"$F_z$"]
    
    plt.plot(offsets, f_total[:, 0], label=labels[0])
    plt.plot(offsets, f_total[:, 1], label=labels[1])
    plt.plot(offsets, f_total[:, 2], label=labels[2])
    
    if spacing is not None and compare is not None:
        plt.scatter(spacing, compare[0])
        plt.scatter(spacing, compare[1])
        plt.scatter(spacing, compare[2])

    plt.legend(loc='best')
    if plot_type.lower() == 'torque':
        plt.ylabel(r'$\tau$ (mN.m)')
    else:
        plt.ylabel(r'$F$ (N)')
    plt.xlabel(r"$d$ (mm)")
    plt.grid(True)
    plt.show()


def gen_mesh_test(path, file, offset=0, Jr=1.0, theta=0, phi=0, alpha=0, beta=0, gamma=0, mask_magnet=False):
    pm.reset()

    center = (0, 0, -10)

    m0 = pm.magnets.Mesh(
        path + file,
        Jr=Jr,
        center=center,
        mesh_scale=5,
        theta=theta,
        phi=phi,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )
    
    center = (0.0 + offset, 0.0, 10)

    m1 = pm.magnets.Mesh(
        path + file,
        Jr=Jr,
        center=center,
        mesh_scale=5,
        theta=theta,
        phi=phi,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    return m0, m1


path = './examples/scripts/stl_magnets/stl/'
file = "cube.stl"

force_points = 101
f_total = np.zeros((force_points, 3))
t_total = np.zeros((force_points, 3))
offsets = np.linspace(-20, 20, force_points)


for i in range(force_points):
    m0, m1 = gen_mesh_test(path, file, offset=offsets[i])
    f_total[i], t_total[i] = m0.get_force_torque(depth=4, unit='mm')
