# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Two Dimensional Magnet Classes
"""
import numpy as _np
from ._magnet import Magnet
from ._fields import Point2
# from typing import List, Any, Union

__all__ = ['Magnet_2D', 'Rectangle', 'Square']

class Magnet_2D(Magnet):    
    """2D Magnet class.

    Args:
        width [float]: magnet width
        height [float]: magnet height
        Jr [float]: Remnant magnetisation

    Kwargs:
        center [Tuple(float, float), or Point2]: center of magnet, defaults to
        Point2(0.0, 0.0)
    """    
    mag_type = 'Magnet_2D'


    def __init__(self, width, height, Jr, **kwargs) -> None:
        super().__init__()
        self.width = width
        self.height = height
        self.a = width/2
        self.b = height/2
        self.Jr = Jr
        

        center = kwargs.pop('center', Point2(0.0, 0.0))

        if type(center) is tuple:
            center = Point2(center[0], center[1])

        self.xc = center.x
        self.yc = center.y

    def center(self):
        """Returns magnet centre

        Returns:
            center[ndarray]: numpy array [xc, yc]
        """        
        return _np.array([self.xc, self.yc])

    def size(self):
        """Returns magnet dimesions

        Returns:
            size[ndarray]: numpy array [width, height]
        """        
        return _np.array([self.width, self.height])


class Rectangle(Magnet_2D):
    """Rectangle Magnet Class

    Args:
        width [float]: magnet width [m] (defaults to 20e-3)
        height [float]: magnet height [m] (defaults to 40e-3)
        Jr [float]: Remnant magnetisation [T] (defaults to 1.0)
         
    Optional Arguments:
         centre: magnet centre (Point2 object) [default Point2(0.0, 0.0)]
         theta: Angle between magnetisation and x-axis [default 90.0 degrees] 

    """
    mag_type = 'Rectangle'

    def __init__(self,
                 width=20e-3, height=40e-3,  # magnet dimensions
                 Jr=1.0,  # local magnetisation
                 **kwargs):
        super().__init__(width, height, Jr, **kwargs)

        self.theta = kwargs.pop('theta', 90)
        self.theta_rad = _np.deg2rad(self.theta)

        self.Jx = _np.around(Jr * _np.cos(self.theta_rad), decimals=6)
        self.Jy = _np.around(Jr * _np.sin(self.theta_rad), decimals=6)
        self.tol = 1e-4  # sufficient for 0.01 degree accuracy

    def __str__(self):
        str = f'{self.__class__.mag_type}\n' + \
              f'J: {self.get_Jr()} (T)\n' + \
              f'Size: {self.size() * 2} (m)\n' + \
              f'Center {self.center()} (m)\n'
        return str

    def __repr__(self):
        str = f'{self.__class__.mag_type}\n' + \
              f'J: {self.get_Jr()} (T)\n' + \
              f'Size: {self.size()} (m)\n' + \
              f'Center {self.center()} (m)\n'
        return str

    def get_Jr(self):
        """Returns Magnetisation components [Jx, Jy]

        Returns:
            J[ndarray]: [Jx, Jy]
        """        
        return _np.array([self.Jx, self.Jy])

    def _calcBx_mag_x(self, x, y):
        """ Bx using 2D Model for rectangular sheets magnetised in x-plane

        Args:
            magnet ([type]): [magnet object]
            x ([type]): [x array]
            y ([type]): [y array]

        Returns:
            Bx[ndarray]
        """
        a = self.a
        b = self.b
        J = self.Jx
        return (J / (2 * _np.pi)) * (
            _np.arctan2((2 * a * (b + y)), (x**2 - a**2 + (y + b)**2))
            + _np.arctan2((2 * a * (b - y)), (x**2 - a**2 + (y - b)**2))
        )

    def _calcBy_mag_x(self, x, y):
        """ By using 2D Model for rectangular sheets magnetised in x-plane

        Args:
            magnet ([type]): [magnet object]
            x ([type]): [x array]
            y ([type]): [y array]

        Returns:
            [array]: [By]
        """
        a = self.a
        b = self.b
        J = self.Jx
        return (-J / (4 * _np.pi)) * (
            _np.log(((x - a)**2 + (y - b)**2) / ((x + a)**2 + (y - b)**2))
            - _np.log(((x - a)**2 + (y + b)**2) / ((x + a)**2 + (y + b)**2))
        )

    def _calcBx_mag_y(self, x, y):
        """ Bx using 2D Model for rectangular sheets magnetised in y-plane

        Args:
            magnet ([type]): [magnet object]
            x ([type]): [x array]
            y ([type]): [y array]

        Returns:
            [array]: [Bx]
        """
        a = self.a
        b = self.b
        J = self.Jy
        return (J / (4 * _np.pi)) * (
            _np.log(((x + a)**2 + (y - b)**2) / ((x + a)**2 + (y + b)**2))
            - _np.log(((x - a)**2 + (y - b)**2) / ((x - a)**2 + (y + b)**2))
        )

    def _calcBy_mag_y(self, x:float, y: float) -> float:
        """ By using 2D Model for rectangular sheets magnetised in y-plane

        Args:
            magnet ([type]): [magnet object]
            x ([type]): [x array]
            y ([type]): [y array]

        Returns:
            [array]: [By]
        """
        a = self.a
        b = self.b
        J = self.Jy
        return (J / (2 * _np.pi)) * (
            _np.arctan2((2 * b * (x + a)), ((x + a)**2 + y**2 - b**2))
            - _np.arctan2((2 * b * (x - a)), ((x - a)**2 + y**2 - b**2))
        )


class Square(Rectangle):
    """Square Magnet Class

    Args:
        width [float]: magnet side length [m] (defaults to 20e-3)
        Jr [float]: Remnant magnetisation [T] (defaults to 1.0)
         
    Optional Arguments:
         centre: magnet centre (Point2 object) [default Point2(0.0, 0.0)]
         theta: Angle between magnetisation and x-axis [default 90.0 degrees] 

    """
    mag_type = 'Square'

    def __init__(self,
                 width=20e-3,   # magnet dimensions
                 Jr=1.0,  # local magnetisation direction
                 **kwargs):

        super().__init__(width=width, height=width, Jr=Jr, **kwargs)
