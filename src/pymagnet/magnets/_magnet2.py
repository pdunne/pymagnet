# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""2D Magnet classes

This private module implements the rectangle and square 2D magnet classes

TODO:
    * Add __del__ method for removing strong ref in class instance list

"""
import numpy as _np
from ._magnet import Magnet
from ._fields import Point2

# from typing import List, Any, Union

__all__ = ["Magnet_2D", "Rectangle", "Square", "Circle"]


class Magnet_2D(Magnet):
    """2D Magnet class.

    Args:
        width [float]: magnet width
        height [float]: magnet height
        Jr [float]: Remnant magnetisation
        **kwargs: Arbitrary keyword arguments.

    kwargs:
        center [Tuple(float, float), or Point2]: center of magnet, defaults to
        Point2(0.0, 0.0)
    """

    mag_type = "Magnet_2D"

    def __init__(self, Jr, **kwargs) -> None:
        super().__init__()
        self.Jr = Jr

        # Magnet rotation w.r.t. x-axis
        self.alpha = kwargs.pop("alpha", 0.0)
        self.alpha_radians = _np.deg2rad(self.alpha)

        center = kwargs.pop("center", Point2(0.0, 0.0))

        if type(center) is tuple:
            center = Point2(center[0], center[1])

        self.xc = center.x
        self.yc = center.y

    def center(self):
        """Returns magnet centre

        Returns:
            center (ndarray): numpy array [xc, yc]
        """
        return _np.array([self.xc, self.yc])

    def get_orientation(self):
        """Returns magnet orientation, `alpha` in degrees

        Returns:
            float: alpha, rotation angle w.r.t x-axis.
        """

        return self.alpha


class Rectangle(Magnet_2D):
    """Rectangle Magnet Class

    Args:
        width [float]: magnet width [m] (defaults to 20e-3)
        height [float]: magnet height [m] (defaults to 40e-3)
        Jr [float]: Remnant magnetisation [T] (defaults to 1.0)

    Optional Arguments:
         centre: magnet centre (Point2 object) [default Point2(0.0, 0.0)]
         phi: Angle between magnetisation and x-axis [default 90.0 degrees]

    """

    mag_type = "Rectangle"

    def __init__(self, width=20e-3, height=40e-3, Jr=1.0, **kwargs):
        super().__init__(Jr, **kwargs)
        self.width = width
        self.height = height

        self.a = width / 2
        self.b = height / 2

        self.phi = kwargs.pop("phi", 90)
        self.phi_rad = _np.deg2rad(self.phi)

        self.Jx = _np.around(Jr * _np.cos(self.phi_rad), decimals=6)
        self.Jy = _np.around(Jr * _np.sin(self.phi_rad), decimals=6)
        self.tol = 1e-4  # sufficient for 0.01 degree accuracy

    def size(self):
        """Returns magnet dimesions

        Returns:
            size[ndarray]: numpy array [width, height]
        """
        return _np.array([self.width, self.height])

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size()} (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation: alpha {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size()} (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation: alpha {self.get_orientation()}\n"
        )
        return str

    def get_Jr(self):
        """Returns Magnetisation components [Jx, Jy]

        Returns:
            J[ndarray]: [Jx, Jy]
        """
        return _np.array([self.Jx, self.Jy])

    def calcB(self, x, y):
        """Calculates the magnetic field at point(s) x,y due to a rectangular magnet

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates

        Returns:
            Vector2: magnetic field vector
        """
        from ._routines2 import rotate_points_2D, _get_field_array_shape2

        array_shape = _get_field_array_shape2(x, y)
        Bx, By = _np.zeros(array_shape), _np.zeros(array_shape)

        if _np.fabs(self.alpha_radians) > Magnet_2D.tol:
            xi, yi = rotate_points_2D(x, y, self.alpha_radians)
            xci, yci = rotate_points_2D(
                _np.array([self.xc]), _np.array([self.yc]), self.alpha_radians
            )

        # Calculate field due to x-component of magnetisation
        if _np.fabs(self.Jx) > Magnet_2D.tol:
            if _np.fabs(self.alpha_radians) > Magnet_2D.tol:

                # Calculate fields in local frame
                Btx = self._calcBx_mag_x(xi - xci[0], yi - yci[0])
                Bty = self._calcBy_mag_x(xi - xci[0], yi - yci[0])

                # Rotate fields to global frame
                Bx, By = rotate_points_2D(Btx, Bty, 2 * _np.pi - self.alpha_radians)

            else:
                Bx = self._calcBx_mag_x(x - self.xc, y - self.yc)
                By = self._calcBy_mag_x(x - self.xc, y - self.yc)

        # Calculate field due to y-component of magnetisation
        if _np.fabs(self.Jy) > Magnet_2D.tol:
            if _np.fabs(self.alpha_radians) > Magnet_2D.tol:

                Btx = self._calcBx_mag_y(xi - xci[0], yi - yci[0])
                Bty = self._calcBy_mag_y(xi - xci[0], yi - yci[0])

                Bxt, Byt = rotate_points_2D(Btx, Bty, 2 * _np.pi - self.alpha_radians)
                Bx += Bxt
                By += Byt
            else:

                Bx += self._calcBx_mag_y(x - self.xc, y - self.yc)
                By += self._calcBy_mag_y(x - self.xc, y - self.yc)
        return Bx, By

    def _calcBx_mag_x(self, x, y):
        """Bx using 2D Model for rectangular sheets magnetised in x-plane

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
        # Hide the warning for situtations where there is a divide by zero.
        # This returns a NaN in the array, which is ignored for plotting.
        with _np.errstate(divide="ignore", invalid="ignore"):
            return (J / (2 * _np.pi)) * (
                _np.arctan2((2 * a * (b + y)), (x ** 2 - a ** 2 + (y + b) ** 2))
                + _np.arctan2((2 * a * (b - y)), (x ** 2 - a ** 2 + (y - b) ** 2))
            )

    def _calcBy_mag_x(self, x, y):
        """By using 2D Model for rectangular sheets magnetised in x-plane

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
        # Hide the warning for situtations where there is a divide by zero.
        # This returns a NaN in the array, which is ignored for plotting.
        with _np.errstate(divide="ignore", invalid="ignore"):
            return (-J / (4 * _np.pi)) * (
                _np.log(((x - a) ** 2 + (y - b) ** 2) / ((x + a) ** 2 + (y - b) ** 2))
                - _np.log(((x - a) ** 2 + (y + b) ** 2) / ((x + a) ** 2 + (y + b) ** 2))
            )

    def _calcBx_mag_y(self, x, y):
        """Bx using 2D Model for rectangular sheets magnetised in y-plane

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
        # Hide the warning for situtations where there is a divide by zero.
        # This returns a NaN in the array, which is ignored for plotting.
        with _np.errstate(divide="ignore", invalid="ignore"):
            return (J / (4 * _np.pi)) * (
                _np.log(((x + a) ** 2 + (y - b) ** 2) / ((x + a) ** 2 + (y + b) ** 2))
                - _np.log(((x - a) ** 2 + (y - b) ** 2) / ((x - a) ** 2 + (y + b) ** 2))
            )

    def _calcBy_mag_y(self, x: float, y: float) -> float:
        """By using 2D Model for rectangular sheets magnetised in y-plane

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
            _np.arctan2((2 * b * (x + a)), ((x + a) ** 2 + y ** 2 - b ** 2))
            - _np.arctan2((2 * b * (x - a)), ((x - a) ** 2 + y ** 2 - b ** 2))
        )


class Square(Rectangle):
    """Square Magnet Class

    Args:
        width [float]: magnet side length [m] (defaults to 20e-3)
        Jr [float]: Remnant magnetisation [T] (defaults to 1.0)

    Optional Arguments:
         centre: magnet centre (Point2 object) [default Point2(0.0, 0.0)]
         phi: Angle between magnetisation and x-axis [default 90.0 degrees]

    """

    mag_type = "Square"

    def __init__(self, width=20e-3, Jr=1.0, **kwargs):
        super().__init__(width=width, height=width, Jr=Jr, **kwargs)


class Circle(Magnet_2D):
    """Circle Magnet Class

    Args:
        radius [float]: magnet radius [m] (defaults to 10e-3)
        Jr [float]: Remnant magnetisation [T] (defaults to 1.0)

    Optional Arguments:
         center: magnet centre (Point2 object) [default Point2(0.0, 0.0)]

    """

    mag_type = "Circle"

    def __init__(
        self,
        radius=10e-3,
        Jr=1.0,  # local magnetisation
        **kwargs,
    ):
        super().__init__(Jr, **kwargs)
        self.radius = radius

        center = kwargs.pop("center", Point2(0.0, 0.0))

        if type(center) is tuple:
            center = Point2(center[0], center[1])

        self.xc = center.x
        self.yc = center.y

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size()} (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation: alpha {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size()} (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation: alpha {self.get_orientation()}\n"
        )
        return str

    def get_Jr(self):
        """Returns Magnetisation components [Jr]

        Returns:
            ndarray: [Jr]
        """
        return _np.array([self.Jr])

    def calcB(self, x, y):
        """Calculates the magnetic field due to long bipolar cylinder

        Args:
            x (float/array): x coordinates
            y (float/array): y coordinates

        Returns:
            tuple: Bx, By magnetic field in cartesian coordinates
        """
        from ._routines import cart2pol, vector_pol2cart
        from ._routines2 import rotate_points_2D

        if _np.fabs(self.alpha_radians) > Magnet_2D.tol:
            xi, yi = rotate_points_2D(x, y, self.alpha_radians)
            xci, yci = rotate_points_2D(
                _np.array([self.xc]), _np.array([self.yc]), self.alpha_radians
            )
            rho, phi = cart2pol(xi - xci[0], yi - yci[0])

            Brho, Bphi = self._calcB_polar(rho, phi)

            # Convert magnetic fields from cylindrical to cartesian
            Bx, By = vector_pol2cart(Brho, Bphi, phi)
            Bx, By = rotate_points_2D(Bx, By, 2 * _np.pi - self.alpha_radians)
            return Bx, By

        rho, phi = cart2pol(x - self.xc, y - self.yc)

        Brho, Bphi = self._calcB_polar(rho, phi)

        # Convert magnetic fields from cylindrical to cartesian
        Bx, By = vector_pol2cart(Brho, Bphi, phi)

        return Bx, By

    def _calcB_polar(self, rho, phi):
        """Calculates the magnetic field due to long bipolar cylinder in polar
        coordinates

        Args:
            rho (float/array): radial values
            phi (float/array): azimuthal values

        Returns:
            tuple: Br, Bphi magnetic field in polar coordinates
        """
        prefac = self.Jr * (self.radius ** 2 / rho ** 2) / 2

        Brho = prefac * _np.cos(phi)
        Bphi = prefac * _np.sin(phi)

        return Brho, Bphi
