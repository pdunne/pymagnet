# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""2D Magnet Classes

This private module implements the `Rectangle` and `Square` and `Circle` 2D magnet classes. The parent class `Magnet2D` implements the location and orientation
methods, i.e. magnet center and quaternion methods for rotating the magnet with respect
to each principal axis.

"""
import numpy as _np
from ._magnet_base import Magnet
from ..utils.global_const import MAG_TOL, PI

__all__ = ["Magnet2D", "Rectangle", "Square", "Circle"]


class Magnet2D(Magnet):
    """2D Magnet Base Class"""

    mag_type = "Magnet2D"

    def __init__(self, Jr, **kwargs) -> None:
        """Init Method

        Args:
            Jr (float): signed magnetised of remnant magnetisationnega

        Kwargs:
            alpha (float): Magnetisation orientation angle (in degrees). Defaults to 0.
            center (tuple or ndarrray): magnet center (x, y). Defaults to (0,0).
        """
        super().__init__()
        self.Jr = Jr

        # Magnet rotation w.r.t. x-axis
        self.alpha = kwargs.pop("alpha", 0.0)
        self.alpha_radians = _np.deg2rad(self.alpha)

        self.center = kwargs.pop("center", _np.array([0.0, 0.0]))
        self.center = _np.asarray(self.center)

    def get_center(self):
        """Returns magnet centre

        Returns:
            center (ndarray): numpy array [xc, yc]
        """
        return self.center

    def get_orientation(self):
        """Returns magnet orientation, `alpha` in degrees

        Returns:
            float: alpha, rotation angle w.r.t x-axis.
        """

        return self.alpha

    def get_field(self):
        """Calculates the magnetic field.

        This is a template that needs to be implemented for each magnet
        """
        pass

    def get_force_torque(self):
        """Calculates the force and torque on a magnet due to all other magnets.

        This is a template that needs to be implemented for each magnet.
        """
        pass


class Rectangle(Magnet2D):
    """Rectangular 2D Magnet Class"""

    mag_type = "Rectangle"

    def __init__(self, width=20.0, height=40.0, Jr=1.0, **kwargs):
        """Init Method

        Args:
            width (float, optional): Magnet Width. Defaults to 20.0.
            height (float, optional): Magnet Height. Defaults to 40.0.
            Jr (float, optional): Remnant Magnetisation. Defaults to 1.0.

        Kwargs:
            alpha (float): Magnetisation orientation angle (in degrees). Defaults to 0.
            center (tuple or ndarrray): magnet center (x, y). Defaults to (0,0).
            phi (float): Rotation Angle (in degrees) of magnet w.r.t x-axis. Defaults to 90.
        """
        super().__init__(Jr, **kwargs)
        self.width = width
        self.height = height

        self.a = width / 2
        self.b = height / 2

        self.phi = kwargs.pop("phi", 90)
        self.phi_rad = _np.deg2rad(self.phi)

        self.Jx = _np.around(Jr * _np.cos(self.phi_rad), decimals=6)
        self.Jy = _np.around(Jr * _np.sin(self.phi_rad), decimals=6)
        self.tol = MAG_TOL  # sufficient for 0.01 degree accuracy

    def get_size(self):
        """Returns magnet dimesions

        Returns:
            ndarray: numpy array [width, height]
        """
        return _np.array([self.width, self.height])

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.get_size()}\n"
            + f"Center {self.get_center()}\n"
            + f"Orientation: alpha {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.get_size()}\n"
            + f"Center {self.get_center()}\n"
            + f"Orientation: alpha {self.get_orientation()}\n"
        )
        return str

    def get_Jr(self):
        """Returns Magnetisation vector

        Returns:
            ndarray: [Jx, Jy]
        """
        return _np.array([self.Jx, self.Jy])

    def get_field(self, x, y):
        """Calculates the magnetic field at point(s) x,y due to a rectangular magnet

        Args:
            x (ndarray): x co-ordinates
            y (ndarray): y co-ordinates

        Returns:
            tuple: magnetic field vector Bx (ndarray), By (ndarray)
        """
        from ..utils._routines2D import rotate_points_2D, _get_field_array_shape2

        array_shape = _get_field_array_shape2(x, y)
        Bx, By = _np.zeros(array_shape), _np.zeros(array_shape)

        if _np.fabs(self.alpha_radians) > Magnet2D.tol:
            xi, yi = rotate_points_2D(
                x - self.center[0], y - self.center[1], self.alpha_radians
            )

        # Calculate field due to x-component of magnetisation
        if _np.fabs(self.Jx / self.Jr) > Magnet2D.tol:
            if _np.fabs(self.alpha_radians) > Magnet2D.tol:

                # Calculate fields in local frame
                Btx = self._calcBx_mag_x(xi, yi)
                Bty = self._calcBy_mag_x(xi, yi)

                # Rotate fields to global frame
                Bx, By = rotate_points_2D(Btx, Bty, 2 * PI - self.alpha_radians)

            else:
                Bx = self._calcBx_mag_x(x - self.center[0], y - self.center[1])
                By = self._calcBy_mag_x(x - self.center[0], y - self.center[1])

        # Calculate field due to y-component of magnetisation
        if _np.fabs(self.Jy / self.Jr) > Magnet2D.tol:
            if _np.fabs(self.alpha_radians) > Magnet2D.tol:

                Btx = self._calcBx_mag_y(xi, yi)
                Bty = self._calcBy_mag_y(xi, yi)

                Bxt, Byt = rotate_points_2D(Btx, Bty, 2 * PI - self.alpha_radians)
                Bx += Bxt
                By += Byt
            else:

                Bx += self._calcBx_mag_y(x - self.center[0], y - self.center[1])
                By += self._calcBy_mag_y(x - self.center[0], y - self.center[1])
        return Bx, By

    def _calcBx_mag_x(self, x, y):
        """Bx using 2D Model for rectangular sheets magnetised in x-plane

        Args:
            x (ndarray): x-coordinates
            y (ndarray): y-coordinates

        Returns:
            ndarray: Bx, x component of magnetic field
        """
        a = self.a
        b = self.b
        J = self.Jx
        # Hide the warning for situtations where there is a divide by zero.
        # This returns a NaN in the array, which is ignored for plotting.
        with _np.errstate(divide="ignore", invalid="ignore"):
            return (J / (2 * PI)) * (
                _np.arctan2((2 * a * (b + y)), (x ** 2 - a ** 2 + (y + b) ** 2))
                + _np.arctan2((2 * a * (b - y)), (x ** 2 - a ** 2 + (y - b) ** 2))
            )

    def _calcBy_mag_x(self, x, y):
        """By using 2D Model for rectangular sheets magnetised in x-plane

        Args:
            x (ndarray): x-coordinates
            y (ndarray): y-coordinates

        Returns:
            ndarray: By, x component of magnetic field
        """
        a = self.a
        b = self.b
        J = self.Jx
        # Hide the warning for situtations where there is a divide by zero.
        # This returns a NaN in the array, which is ignored for plotting.
        with _np.errstate(divide="ignore", invalid="ignore"):
            return (-J / (4 * PI)) * (
                _np.log(((x - a) ** 2 + (y - b) ** 2) / ((x + a) ** 2 + (y - b) ** 2))
                - _np.log(((x - a) ** 2 + (y + b) ** 2) / ((x + a) ** 2 + (y + b) ** 2))
            )

    def _calcBx_mag_y(self, x, y):
        """Bx using 2D Model for rectangular sheets magnetised in y-plane

        Args:
            x (ndarray): x-coordinates
            y (ndarray): y-coordinates

        Returns:
            ndarray: Bx, x component of magnetic field
        """
        a = self.a
        b = self.b
        J = self.Jy
        # Hide the warning for situtations where there is a divide by zero.
        # This returns a NaN in the array, which is ignored for plotting.
        with _np.errstate(divide="ignore", invalid="ignore"):
            return (J / (4 * PI)) * (
                _np.log(((x + a) ** 2 + (y - b) ** 2) / ((x + a) ** 2 + (y + b) ** 2))
                - _np.log(((x - a) ** 2 + (y - b) ** 2) / ((x - a) ** 2 + (y + b) ** 2))
            )

    def _calcBy_mag_y(self, x: float, y: float) -> float:
        """Bx using 2D Model for rectangular sheets magnetised in y-plane

        Args:
            x (ndarray): x-coordinates
            y (ndarray): y-coordinates

        Returns:
            ndarray: By, x component of magnetic field
        """
        a = self.a
        b = self.b
        J = self.Jy
        return (J / (2 * PI)) * (
            _np.arctan2((2 * b * (x + a)), ((x + a) ** 2 + y ** 2 - b ** 2))
            - _np.arctan2((2 * b * (x - a)), ((x - a) ** 2 + y ** 2 - b ** 2))
        )


class Square(Rectangle):
    """Square 2D Magnet Class"""

    mag_type = "Square"

    def __init__(self, width=20, Jr=1.0, **kwargs):
        """Init Method

        Args:
            width (float, optional): Sqaure side length. Defaults to 20.0.
            Jr (float, optional): Remnant Magnetisation. Defaults to 1.0.

        Kwargs:
             alpha (float): Magnetisation orientation angle (in degrees). Defaults to 0.
            center (tuple or ndarrray): magnet center (x, y). Defaults to (0,0).
            phi (float): Rotation Angle (in degrees) of magnet w.r.t x-axis. Defaults to 90.
        """
        super().__init__(width=width, height=width, Jr=Jr, **kwargs)


class Circle(Magnet2D):
    """Circle 2D Magnet Class"""

    mag_type = "Circle"

    def __init__(
        self,
        radius=10,
        Jr=1.0,  # local magnetisation
        **kwargs,
    ):
        """Init Method

        Args:
            radius (float, optional): Radius. Defaults to 10.0.
            Jr (float, optional): Remnant magnetisation. Defaults to 1.0.

        Kwargs:
            alpha (float): Unused. For rotations use phi instead
            center (tuple or ndarrray): magnet center (x, y). Defaults to (0,0)
            phi (float): Rotation Angle (in degrees) of magnet w.r.t x-axis. Defaults to 90.
        """
        super().__init__(Jr, **kwargs)
        self.radius = radius
        self.phi = kwargs.pop("phi", 0)
        self.phi_rad = _np.deg2rad(self.phi)

        self.Jx = _np.around(Jr * _np.cos(self.phi_rad), decimals=6)
        self.Jy = _np.around(Jr * _np.sin(self.phi_rad), decimals=6)
        self.tol = MAG_TOL  # sufficient for 0.01 degree accuracy

        self.center = kwargs.pop("center", _np.array([0.0, 0.0]))
        self.center = _np.asarray(self.center)

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.get_size()}\n"
            + f"Center {self.get_center()}\n"
            + f"Orientation: alpha {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.get_size()}\n"
            + f"Center {self.get_center()}\n"
            + f"Orientation: alpha {self.get_orientation()}\n"
        )
        return str

    def get_size(self):
        """Returns radius

        Returns:
            ndarray: radius
        """
        return _np.array([self.radius])

    def get_Jr(self):
        """Returns signed remnant magnetisation

        Returns:
            ndarray: remnant magnetisation
        """
        return _np.array([self.Jx, self.Jy])

    def get_field(self, x, y):
        """Calculates the magnetic field due to long bipolar cylinder

        Args:
            x (ndarray): x coordinates
            y (narray): y coordinates

        Returns:
            tuple: Bx, By magnetic field in cartesian coordinates
        """
        from ..utils._conversions import cart2pol, vector_pol2cart
        from ..utils._routines2D import rotate_points_2D

        if _np.fabs(self.alpha_radians) > Magnet2D.tol:
            xi, yi = rotate_points_2D(
                x - self.center[0], y - self.center[1], self.alpha_radians
            )

            rho, phi = cart2pol(xi, yi)

            Brho, Bphi = self._calcB_polar(rho, phi - self.phi_rad)

            # Convert magnetic fields from cylindrical to cartesian
            Bx, By = vector_pol2cart(Brho, Bphi, phi)
            Bx, By = rotate_points_2D(Bx, By, 2 * PI - self.alpha_radians)
            return Bx, By

        rho, phi = cart2pol(x - self.center[0], y - self.center[1])

        Brho, Bphi = self._calcB_polar(rho, phi - self.phi_rad)

        # Convert magnetic fields from cylindrical to cartesian
        Bx, By = vector_pol2cart(Brho, Bphi, phi)

        return Bx, By

    def _calcB_polar(self, rho, phi):
        """Calculates the magnetic field due to long bipolar cylinder in polar
        coordinates

        Args:
            rho (ndarray): radial values
            phi (ndarray): azimuthal values

        Returns:
            tuple: Br, Bphi magnetic field in polar coordinates
        """
        prefac = self.Jr * (self.radius ** 2 / rho ** 2) / 2

        Brho = prefac * _np.cos(phi)
        Bphi = prefac * _np.sin(phi)

        return Brho, Bphi
