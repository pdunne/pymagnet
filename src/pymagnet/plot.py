# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Plotting routines


"""
import matplotlib.pyplot as _plt
from matplotlib.patches import Rectangle as _Rect
from matplotlib.patches import Arrow as _Arrow
import numpy as _np
from pymagnet import magnets as _mag

def plot_1D_field(magnet, **kwargs):

    NP = kwargs.pop('NP', 101)
    return_data = kwargs.pop('return_data', False)
    
    if issubclass(magnet.__class__, _mag.Cylinder):
        mag_boundary = magnet.length/2
        z = _np.linspace(-2 * magnet.length +  magnet.zc, 2 * magnet.length + magnet.zc, NP)
        Bz = _mag.magnetic_field_cylinder_1D(magnet, z)
        
    elif issubclass(magnet.__class__, _mag.Prism):
        mag_boundary = magnet.height / 2
        z = _np.linspace(-2*magnet.height + magnet.zc,
                         2*magnet.height + magnet.zc, NP)
        Bz = _mag.magnetic_field_prism_1D(magnet, z)

    else:
        print("Error")
        return None

    _, _ = _plt.subplots()
    _plt.xlabel(r'$z$ (mm)')
    _plt.ylabel(r'$B_z$ (mT)')
    _plt.plot(z*1e3, Bz*1e3,  label='Cube')
    _plt.axvline(x=-mag_boundary*1e3 + magnet.zc*1e3, c='blue', ls='--')
    _plt.axvline(x=mag_boundary*1e3 + magnet.zc*1e3, c='red', ls='--')
    _plt.axvline(x=0.0, c='k', ls='-')
    _plt.show()
    if return_data:
        return z,Bz
    else:
        return None


def plot_2D_line(x, Field, **kwargs):
    """Line Plot of field from 2D magnet

    Args:
        x ([array]): [assumes in m]
        Field ([field vector_2D]): [field vector object]
    """
    scale_x = kwargs.pop('scale_x', 1e-3)
    scale_cb = kwargs.pop('scale_cb', 1)


    xlab = kwargs.pop('xlab', f'x (mm)')
    ylab = kwargs.pop('ylab', f'B (T)')
    axis_scale = kwargs.pop('axis_scale', 'equal')

    SAVE = kwargs.pop('save_fig', False)

    _, ax = _plt.subplots()
    _plt.plot(x/scale_x, Field.n/scale_cb, label=r'$|\mathbf{B}|$')
    _plt.plot(x/scale_x, Field.x/scale_cb, label=r'$B_x$')
    _plt.plot(x/scale_x, Field.y/scale_cb, label=r'$B_y$')
    _plt.legend(loc='best')
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)

    _plt.show()

    if SAVE:
        _plt.savefig('line_plot.png', dpi=300)
        # _plt.savefig('contour_plot.pdf', dpi=300)


def plot_2D_contour(x, y, Field, **kwargs):
    """Contour Plot of field

    Args:
        x ([array]): [assumes in µm]
        y ([array]): [assumes in nm]
        Field ([field vector_2D]): [field vector object]
    """
    scale_x = kwargs.pop('scale_x', 1e-3)
    scale_y = kwargs.pop('scale_y', 1e-3)
    scale_cb = kwargs.pop('scale_cb', 1)

    show_magnets = kwargs.pop('show_magnets', True)

    cmap = kwargs.pop('cmap', 'magma')
    xlab = kwargs.pop('xlab', f'x (mm)')
    ylab = kwargs.pop('ylab', f'y (mm)')
    clab = kwargs.pop('clab', f'B (T)')
    axis_scale = kwargs.pop('axis_scale', 'equal')

    field_component = kwargs.pop('field_component', 'n')

    vector_plot = kwargs.pop('vector_plot', False)
    vector_color = kwargs.pop('vector_color','w')
    NQ = kwargs.pop('vector_arrows', 11)

    if field_component == 'x':
        field_chosen = Field.x
    elif field_component == 'y':
        field_chosen = Field.y
    else:
        field_chosen = Field.n

    SAVE = kwargs.pop('save_fig', False)

    finite_field = field_chosen[_np.isfinite(field_chosen)]/scale_cb
    LL = kwargs.pop('LL', 0)
    UL = kwargs.pop('UL', round(finite_field.mean()*2, 1))
    NL = kwargs.pop('NL', 11)

    lev2 = _np.linspace(LL, UL, 256, endpoint=True)

    _, ax = _plt.subplots()
    CS = _plt.contourf(x/scale_x, y/scale_y, field_chosen/scale_cb,
                       levels=lev2, cmap=_plt.get_cmap(cmap), extend='max')

    if NL > 1:
        lev1 = _np.linspace(LL, UL, NL, endpoint=True)
        _ = _plt.contour(x/scale_x, y/scale_y, field_chosen/scale_cb,
                          vmin=LL, vmax=UL, levels=lev1,
                          linewidths=1.0, colors='k')
        CB = _plt.colorbar(CS, ticks=lev1)
    else:
        CB = _plt.colorbar(CS)


    if vector_plot:
        NPx, NPy = x.shape
        if NQ != 0:
            NSx, NSy = NPx//NQ, NPy//NQ
            _plt.quiver(x[::NSx, ::NSy]/scale_x,
                        y[::NSx, ::NSy]/scale_y,
                        Field.x[::NSx, ::NSy]/Field.n[::NSx, ::NSy],
                        Field.y[::NSx, ::NSy]/Field.n[::NSx, ::NSy],
                        color=vector_color,
                        alpha=1)
    
    if show_magnets:
        patch_array =_num_patch_2D()
        for p in patch_array:
        
            sq = _Rect(xy=(p.patch.x/scale_x, p.patch.y/scale_y),
                       width=p.patch.width/scale_x,
                       height=p.patch.height/scale_y,
                       fill=True, facecolor='#F5F5F5', edgecolor ='k', zorder=5)
            ax.add_patch(sq)

            ar = _Arrow(p.arrow.x/scale_x, p.arrow.y/scale_y,
                        p.arrow.dx/scale_x, p.arrow.dy/scale_y,
                        width = p.arrow.width,
                        zorder=6,
                        color='k')
            ax.add_patch(ar)

    CB.ax.get_yaxis().labelpad = 15
    CB.ax.set_ylabel(clab, rotation=270)
    _plt.axis(axis_scale)
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)
    _plt.show()

    if SAVE:
        _plt.savefig('contour_plot.png', dpi=300)
        # _plt.savefig('contour_plot.pdf', dpi=300)


class patch(object):
    """Encodes magnet dimensions for drawing on plots

    Args:
        object ([type]): [description]
    """    
    def __init__(self, x, y, width, height ):
        super().__init__()
        self.x = x
        self.y = y
        self.width = width
        self.height = height

    def __repr__(self) -> str:
        return f"(x: {self.x}, y: {self.y} w:{self.width}, h: {self.height})"

    def __str__(self) -> str:
        return f"(x: {self.x}, y: {self.y} w:{self.width}, h: {self.height})"

class arrow(object):
    """Encodes magnetisation vector for drawing on plots

    Args:
        object ([type]): [description]
    """    
    def __init__(self, x, y, dx, dy, width=3):
        super().__init__()
        self.x = x
        self.y = y
        self.dx = dx
        self.dy = dy
        self.width = width

    def __repr__(self) -> str:
        return f"(x: {self.x}, y: {self.y}, dx: {self.dx}, dy: {self.dy}, w:{self.width})"

    def __str__(self) -> str:
        return f"(x: {self.x}, y: {self.y}, dx: {self.dx}, dy: {self.dy}, w:{self.width})"

class magnet_patch(object):
    """Magnet drawing class

    Args:
        object ([type]): [description]
    """    
    def __init__(self, patch, arrow)    -> None:
        super().__init__()
        self.patch = patch
        self.arrow = arrow

    def __str__(self) -> str:
        return self.patch.__str__() + self.arrow.__str__()

def _num_patch_2D():
    """Generates patches and arrows for drawing 

    Returns:
        [tuple(list, list)]: lists of patch and arrow objects
    """    
    from .magnets._magnet2 import Magnet_2D
    patch_array = []
    
    for magnet in Magnet_2D.instances:
        patch_tmp = patch(x = magnet.center()[0] - magnet.a,
                          y = magnet.center()[1] - magnet.b,
                          width= 2*magnet.a,
                          height= 2*magnet.b)
        Jnorm = magnet.get_Jr()/_np.abs(magnet.Jr)
        offset = _np.multiply(Jnorm, magnet.size())/2

        arrow_tmp = arrow(x = magnet.center()[0] - offset[0],
                          y = magnet.center()[1] - offset[1],
                          dx= 2 * offset[0],
                          dy= 2 * offset[1])

        magnet_patch_tmp = magnet_patch(patch_tmp, arrow_tmp)

        patch_array.append(magnet_patch_tmp)
    return patch_array


def plot_3D_contour(x, y, z, Field, **kwargs):
    """Contour Plot of field

    Args:
        x ([array]): [assumes in µm]
        y ([array]): [assumes in nm]
        Field ([field vector]): [field vector object]
    """
    scale_x = kwargs.pop('scale_x', 1e-3)
    scale_y = kwargs.pop('scale_y', 1e-3)
    scale_z = kwargs.pop('scale_z', 1e-3)
    scale_cb = kwargs.pop('scale_cb', 1)
    axis_scale = kwargs.pop('axis_scale', 'equal')

    patch_array = kwargs.pop('patch_array', _np.array([]))

    cmap = kwargs.pop('cmap', 'magma')
    xlab = kwargs.pop('xlab', f'x (mm)')
    ylab = kwargs.pop('ylab', f'y (mm)')
    zlab = kwargs.pop('zlab', f'z (mm)')
    clab = kwargs.pop('clab', f'B (T)')

    SAVE = kwargs.pop('save_fig', False)

    finite_field = Field.n[_np.isfinite(Field.n)]/scale_cb
    LL = kwargs.pop('LL', 0)
    UL = kwargs.pop('UL', round(finite_field.mean()*2, 1))
    NL = kwargs.pop('NL', 11)

    lev2 = _np.linspace(LL, UL, 256, endpoint=True)


    if _np.ndim(z) < 2:
        orientation = 'xy'
        plot_x = x/scale_x
        plot_y = y/scale_y
        plot_xlab = xlab
        plot_ylab = ylab

    elif _np.ndim(y) < 2:
        orientation = 'xz'
        plot_x = x/scale_x
        plot_y = z/scale_z
        plot_xlab = xlab
        plot_ylab = zlab

    else:
        orientation = 'yz'
        plot_x = y/scale_y
        plot_y = z/scale_z
        plot_xlab = ylab
        plot_ylab = zlab

    print(orientation)
    _, ax = _plt.subplots()
    CS = _plt.contourf(plot_x, plot_y, Field.n/scale_cb,
                       levels=lev2, cmap=_plt.get_cmap(cmap), extend='max')

    if NL > 1:
        lev1 = _np.linspace(LL, UL, NL, endpoint=True)
        _ = _plt.contour(plot_x, plot_y, Field.n/scale_cb,
                         vmin=LL, vmax=UL, levels=lev1,
                         linewidths=1.0, colors='k')
        CB = _plt.colorbar(CS, ticks=lev1)
    else:
        CB = _plt.colorbar(CS)

    if len(patch_array) > 0:
        for ii in patch_array:
            ii /= scale_x

            sq = _Rect(xy=(ii[0], ii[1]), width=ii[2], height=ii[3],
                       fill=True, color='#F5F5F5', zorder=5)
            ax.add_patch(sq)

            ar = _Arrow(ii[0] + ii[2]/2 - ii[4]/2, ii[1] + ii[3]/2 - ii[5]/2,
                        ii[4], ii[5],
                        width=5, zorder=6,
                        color='k')
            ax.add_patch(ar)

    CB.ax.get_yaxis().labelpad = 15
    CB.ax.set_ylabel(clab, rotation=270)
    _plt.xlabel(plot_xlab)
    _plt.ylabel(plot_ylab)
    _plt.axis(axis_scale)

    if SAVE:
        _plt.savefig('contour_plot.png', dpi=300)
        # _plt.savefig('contour_plot.pdf', dpi=300)


def plot_sub_contour_3D(plot_x, plot_y, plot_B, **kwargs):
    """Contour plot of 3D simulation

    Args:
        plot_x ([type]): [description]
        plot_y ([type]): [description]
        plot_B ([type]): [description]
    """

    cmap = kwargs.pop('cmap', 'seismic')
    xlab = kwargs.pop('xlab', f'x (mm)')
    ylab = kwargs.pop('ylab', f'y (mm)')
    clab = kwargs.pop('clab', f'B (T)')

    axis_scale = kwargs.pop('axis_scale', 'equal')

    SAVE = kwargs.pop('save_fig', False)

    LL = kwargs.pop('LL', -0.5)
    UL = kwargs.pop('UL', 0.5)
    NL = kwargs.pop('NL', 11)

    lev2 = _np.linspace(LL, UL, 256, endpoint=True)
    _, ax = _plt.subplots()
    CS = _plt.contourf(plot_x, plot_y, plot_B,
                      levels=lev2, cmap=_plt.get_cmap(cmap), extend='both')

    if NL > 1:
        lev1 = _np.linspace(LL, UL, NL, endpoint=True)
        _ = _plt.contour(plot_x, plot_y, plot_B,
                        vmin=LL, vmax=UL, levels=lev1,
                        linewidths=1.0, colors='k')
        CB = _plt.colorbar(CS, ticks=lev1)
    else:
        CB = _plt.colorbar(CS)

    CB.ax.get_yaxis().labelpad = 15
    CB.ax.set_ylabel(clab, rotation=270)
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)
    _plt.axis('equal')
    _plt.show()


def line_plot_cylinder(magnet, **kwargs):
    """Calculates and plots the magnetic field along the central axis
    of a cylinder

    This is an example helper function.

    Args:
        magnet ([magnet object]):
    """
    rho = _np.linspace(-2*magnet.radius, 2*magnet.radius, 51)
    z = _np.array([magnet.length * 1.1 / 2])

    Bz, Br = magnet._calcB_cyl(rho, z)
    _, _ = _plt.subplots()
    _plt.plot(rho*1e3, Bz, label='Bz')
    _plt.plot(rho*1e3, Br, label='Br')
    _plt.legend(loc='best')
    _plt.show()

def contour_plot_cylinder(magnet, **kwargs):
    """Calculates and plots the magnetic field 
    of a cylinder

    This is an example helper function.


    Args:
        magnet ([type]): [description]
    """    
    NP = 101
    NPJ = NP*1j
    rho, z = _np.mgrid[-3 * magnet.radius:3 *
                       magnet.radius:NPJ, -magnet.length:magnet.length:NPJ]
    Bz, Br = magnet._calcB_cyl(rho, z)
    Bn = _np.sqrt(Bz**2 + Br**2)

    xlab = f'r (mm)'
    ylab = 'z (mm)'

    # plot_B = Bn
    clab = r'$|B|$ (T)'
    cmap = 'viridis'
    plot_sub_contour_3D(
        rho*1e3, z*1e3, Bn, xlab=xlab, ylab=ylab, clab=clab, cmap=cmap, LL=0, UL=1.0)


def param_test_2D(width, height):
    """Example plots while varying the size of a permanent magnet

    This is an example helper function.

    Args:
        width ([type]): [description]
        height ([type]): [description]
    """    
    x = _np.linspace(-2*width, 2*width, 100)
    y = 1e-3 + height

    B = _mag._routines2.B_calc_2D(x, y)
    plot_2D_line(x, B)
    x, y = _mag._routines2.grid2D(1.5*width, height)
    B = _mag._routines2.B_calc_2D(x, y)
    cmap = 'viridis'
    patch_array = _num_patch_2D()
    plot_2D_contour(
        x, y, B, UL=1, NL=11, cmap=cmap)
