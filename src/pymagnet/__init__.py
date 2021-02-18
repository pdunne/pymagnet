from . import plot
from . import magnets
from .magnets._magnet import reset_magnets, list_magnets
from .magnets._routines2 import grid2D, B_calc_2D
from .magnets._routines3 import grid3D, B_calc_3D
from math import pi as PI

PI_2 = PI / 2.0
PI_4 = PI / 4.0

u0 = 4E-7 * PI

__all__ = ['magnets', 'plot']
