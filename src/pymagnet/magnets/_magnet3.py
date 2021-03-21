# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""3D Magnet classes

This private module implements the `Prism`, `Cube`, `Cylinder`, and `Sphere`
3D magnet classes. The parent class `Magnet_3D` implements the location and orientation
methods, i.e. magnet center and quaternion methods for rotating the magnet with respect
to each principal axis.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

TODO:
    * Add __del__ method for removing strong ref in class instance list
    * Update __str__ and __repr__ methods to show orientation and magnetisation
"""
from numba import vectorize, float64
import numpy as _np
from ._fields import Point3
from ._magnet import Magnet
from ._quaternion import Quaternion

# from pymagma import PI

__all__ = ["Magnet_3D", "Prism", "Cube", "Cylinder", "Sphere"]


class Magnet_3D(Magnet):
    """3D Magnet Class

    Args:
         width: magnet width [default 10e-3]
         depth: magnet depth [default 20e-3]
         height: magnet height [default 30e-3]

    Optional Arguments:
         centre: magnet centre (Point3_3D object) [default Point3(.0, .0, .0)]
         J: remnant magnetisation in Tesla [default 1.0]

    """

    mag_type = "Magnet_3D"

    def __init__(self, Jr, **kwargs) -> None:
        super().__init__()

        self.Jr = Jr

        center = kwargs.pop("center", Point3(0.0, 0.0, 0.0))

        self._mask_magnet = kwargs.pop("mask_magnet", False)

        if type(center) is tuple:
            center = Point3(center[0], center[1], center[2])

        self.xc = center.x
        self.yc = center.y
        self.zc = center.z

        self.alpha = kwargs.pop("alpha", 0.0)  # rotation angle about z
        self.beta = kwargs.pop("beta", 0.0)  # rotation angle about y
        self.gamma = kwargs.pop("gamma", 0.0)  # rotation angle about x

        self.alpha_rad = _np.deg2rad(self.alpha)
        self.beta_rad = _np.deg2rad(self.beta)
        self.gamma_rad = _np.deg2rad(self.gamma)

    def center(self):
        """Returns magnet center

        Returns:
            ndarray: [center_x, center_y, center_z]
        """
        return _np.array([self.xc, self.yc, self.zc])

    def get_Jr(self):
        """Returns local magnetisation orientation

        Must be implemented for all classes
        """
        pass

    def get_orientation(self):
        """Returns magnet orientation, `alpha`, `beta`, `gamma` in degrees

        Returns:
            ndarray: alpha, beta, gamma rotation angles w.r.t z, y, and x axes
        """

        return _np.array([self.alpha, self.beta, self.gamma])

    def _generate_rotation_quaternions(self):
        """Generates single rotation quaternion for all non-zero rotation angles,
        which are:

            alpha: angle in degrees around z-axis
            beta: angle in degrees around y-axis
            gamma: angle in degrees around x-axis

        Returns:
            Quaternion: total rotation quaternion
        """
        from functools import reduce

        # Initialise quaternions
        rotate_about_x, rotate_about_y, rotate_about_z = None, None, None
        forward_rotation, reverse_rotation = None, None

        if _np.fabs(self.alpha_rad) > 1e-4:
            rotate_about_z = Quaternion.q_angle_from_axis(self.alpha_rad, (0, 0, 1))

        if _np.fabs(self.beta_rad) > 1e-4:
            rotate_about_y = Quaternion.q_angle_from_axis(self.beta_rad, (0, 1, 0))

        if _np.fabs(self.gamma_rad) > 1e-4:
            rotate_about_x = Quaternion.q_angle_from_axis(self.gamma_rad, (1, 0, 0))

        # Generate compound rotations
        # Order of rotation: beta  about y, alpha about z, gamma about x
        q_forward = [
            q for q in [rotate_about_y, rotate_about_z, rotate_about_x] if q is not None
        ]

        forward_rotation = reduce(lambda x, y: x * y, q_forward)

        q_reverse = [
            q.get_conjugate()
            for q in [rotate_about_x, rotate_about_z, rotate_about_y]
            if q is not None
        ]

        reverse_rotation = reduce(lambda x, y: x * y, q_reverse)
        return forward_rotation, reverse_rotation

    def calcB(self, x, y, z):
        """Calculates the magnetic field at point(s) x,y,z due to a 3D magnet
        The calculations are always performed in local coordinates with the centre of the magnet at origin and z magnetisation pointing along the local z' axis.

        The rotations and translations are performed first, and the internal field calculation functions are called.

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates

        Returns:
            tuple: Bx(ndarray), By(ndarray), Bz(ndarray)  field vector
        """
        from ._routines3 import _tile_arrays, _apply_mask

        # If any rotation angle is set, transform the data
        if _np.any(
            _np.fabs(
                _np.array(
                    [
                        self.alpha_rad,
                        self.beta_rad,
                        self.gamma_rad,
                    ]
                )
            )
            > Magnet.tol
        ):

            forward_rotation, reverse_rotation = self._generate_rotation_quaternions()

            # Generate 3xN array for quaternion rotation
            pos_vec = Quaternion._prepare_vector(x - self.xc, y - self.yc, z - self.zc)

            # Rotate points
            x_rot, y_rot, z_rot = forward_rotation * pos_vec

            # Calls internal child method to calculate the field
            B_local = self._calcB_local(x_rot, y_rot, z_rot)
            mask = self._generate_mask(x_rot, y_rot, z_rot)

            B_local = _apply_mask(self, B_local, mask)

            # Rearrange the field vectors in a 3xN array for quaternion rotation
            Bvec = Quaternion._prepare_vector(B_local.x, B_local.y, B_local.z)

            # Rotate the local fields back into the global frame using quaternions
            Bx, By, Bz = reverse_rotation * Bvec

            # finally return the fields
            return Bx, By, Bz
        else:
            B = self._calcB_local(x - self.xc, y - self.yc, z - self.zc)

            xloc, yloc, zloc = _tile_arrays(x - self.xc, y - self.yc, z - self.zc)
            mask = self._generate_mask(xloc, yloc, zloc)
            B = _apply_mask(self, B, mask)

            return B.x, B.y, B.z

    def _calcB_local(x, y, z):
        """Internal magnetic field calculation method. This should be defined
        for each magnet type.

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates
        """
        pass

    def _generate_mask(self, x, y, z):
        """Generates mask of points inside a magnet

        Args:
            x (ndarray/float): x-coordinates
            y (ndarray/float): y-coordinates
            z (ndarray/float): z-coordinates
        """
        pass

    def size(self):
        """Returns magnet dimesions

        Must be implemented for each magnet
        """
        pass


class Prism(Magnet_3D):
    """Prism Magnet Class. Cuboid (width, depth height).

    Args:
          width (float): magnet width [default 10e-3]
          depth (float): magnet depth [default 20e-3]
          height (float): magnet height [default 30e-3]
          Jr (float): remnant magnetisation in Tesla [default 1.0]

     Optional Arguments:
          centre: magnet centre (tuple or Point3) [default Point3(.0, .0, .0)]
          theta: Angle between magnetisation and x-axis [default 90.0 degrees]
          phi: Angle between magnetisation and z-axis [default 0.0 degrees]

    """

    mag_type = "Prism"

    def __init__(
        self,
        width=10e-3,
        depth=20e-3,
        height=30e-3,  # magnet dimensions
        Jr=1.0,  # local magnetisation direction
        **kwargs,
    ):
        self.width = width
        self.depth = depth
        self.height = height

        self.a = width / 2
        self.b = depth / 2
        self.c = height / 2

        super().__init__(Jr, **kwargs)

        self.phi = kwargs.pop("theta", 90)
        self.phi_rad = _np.deg2rad(self.phi)
        self.theta = kwargs.pop("phi", 0)
        self.theta_rad = _np.deg2rad(self.theta)

        self.Jx = _np.around(
            Jr * _np.cos(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        )
        self.Jy = _np.around(
            Jr * _np.sin(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        )
        self.Jz = _np.around(Jr * _np.cos(self.theta_rad), decimals=6)
        self.tol = 1e-4  # sufficient for 0.01 degree accuracy

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size()} (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size()} (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def size(self):
        return _np.array([self.a, self.b, self.c])

    def get_Jr(self):
        return _np.array([self.Jx, self.Jy, self.Jz])

    def size(self):
        """Returns magnet dimesions

        Returns:
            size[ndarray]: numpy array [width, depth, height]
        """
        return _np.array([self.width, self.depth, self.height])

    @staticmethod
    def _F1(a, b, c, x, y, z):
        """Helper Function F1 for 3D prismatic magnets

        Args:
            x ([type]): [description]
            y ([type]): [description]
            z ([type]): [description]

        Returns:
            [array]:
        """

        try:
            # Hide the warning for situtations where there is a divide by zero.
            # This returns a NaN in the array, which is ignored for plotting.
            with _np.errstate(divide="ignore", invalid="ignore"):
                data = _np.arctan(
                    ((y + b) * (z + c))
                    / (
                        (x + a)
                        * _np.sqrt(
                            _np.power((x + a), 2)
                            + _np.power((y + b), 2)
                            + _np.power((z + c), 2)
                        )
                    )
                )
        except ValueError:
            data = _np.NaN
        return data

    @staticmethod
    def _F2(a, b, c, x, y, z):
        """Helper Function F2 for 3D prismatic magnets

        Args:
            x (ndarray): x-coordinates
            y (ndarray): y-coordinates
            z (ndarray): z-coordinates

        Returns:
            [array]: [description]
        """

        try:
            xa_sq = _np.power((x + a), 2)
            yb_sq = _np.power((y + b), 2)
            zc_sq = _np.power((z + c), 2)
            znc_sq = _np.power((z - c), 2)
            # Hide the warning for situtations where there is a divide by zero.
            # This returns a NaN in the array, which is ignored for plotting.
            with _np.errstate(divide="ignore", invalid="ignore"):
                data = (_np.sqrt(xa_sq + yb_sq + znc_sq) + c - z) / (
                    _np.sqrt(xa_sq + yb_sq + zc_sq) - c - z
                )
        except ValueError:
            data = _np.NaN
        return data

    def _calcB_local(self, x, y, z):
        """Internal magnetic field calculation methods.
        Iterates over each component of the prism magnetised in x, y, and z
        (in local coordinates).

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates

        Returns:
            Vector3: Magnetic field array
        """
        from ._routines3 import _allocate_field_array3

        B = _allocate_field_array3(x, y, z)
        # Magnetic field due to component of M magnetised in x
        if _np.fabs(self.Jx) > self.tol:
            Bx, By, Bz = self._calcB_prism_x(x, y, z)
            B.x += Bx
            B.y += By
            B.z += Bz

        # Magnetic field due to component of M magnetised in y
        if _np.fabs(self.Jy) > self.tol:
            Bx, By, Bz = self._calcB_prism_y(x, y, z)
            B.x += Bx
            B.y += By
            B.z += Bz

        # Magnetic field due to component of M magnetised in z
        if _np.fabs(self.Jz) > self.tol:
            Bx, By, Bz = self._calcB_prism_z(x, y, z)
            B.x += Bx
            B.y += By
            B.z += Bz
        return B

    def _calcBx_prism_x(self, a, b, c, Jr, x, y, z):
        """x component of magnetic field for prism magnet
        magnetised in x

        Args:
            a ([type]): [description]
            b ([type]): [description]
            c ([type]): [description]
            Jr ([type]): [description]
            x ([type]): [description]
            y ([type]): [description]
            z ([type]): [description]

        Returns:
            [array]: [description]
        """

        try:
            data = -(Jr / (4 * _np.pi)) * (
                self._F1(a, b, c, -x, y, z)
                + self._F1(a, b, c, -x, y, -z)
                + self._F1(a, b, c, -x, -y, z)
                + self._F1(a, b, c, -x, -y, -z)
                + self._F1(a, b, c, x, y, z)
                + self._F1(a, b, c, x, y, -z)
                + self._F1(a, b, c, x, -y, z)
                + self._F1(a, b, c, x, -y, -z)
            )
        except ValueError:
            data = _np.NaN
        return data

    def _calcBy_prism_x(self, a, b, c, Jr, x, y, z):
        """y component of magnetic field for prism magnet
        magnetised in x

        Args:
            a ([type]): [description]
            b ([type]): [description]
            c ([type]): [description]
            Jr ([type]): [description]
            x ([type]): [description]
            y ([type]): [description]
            z ([type]): [description]

        Returns:
            [array]: [description]
        """

        try:
            data = _np.log(
                self._F2(a, b, c, -x, -y, z)
                * self._F2(a, b, c, x, y, z)
                / (self._F2(a, b, c, -x, y, z) * self._F2(a, b, c, x, -y, z))
            )
            data *= Jr / (4 * _np.pi)
        except ValueError:
            data = _np.NaN
        return data

    def _calcBz_prism_x(self, a, b, c, Jr, x, y, z):
        """z component of magnetic field for prism magnet
        magnetised in x

        Args:
            x ([type]): [description]
            y ([type]): [description]
            z ([type]): [description]

        Returns:
            [array]: [description]
        """
        try:
            # Hide the warning for situtations where there is a divide by zero.
            # This returns a NaN in the array, which is ignored for plotting.
            with _np.errstate(divide="ignore", invalid="ignore"):
                data = _np.log(
                    self._F2(a, c, b, -x, -z, y)
                    * self._F2(a, c, b, x, z, y)
                    / (self._F2(a, c, b, -x, z, y) * self._F2(a, c, b, x, -z, y))
                )
            data *= Jr / (4 * _np.pi)
        except ValueError:
            data = _np.NaN
        return data

    def _calcB_prism_x(self, x, y, z):
        """Magnetic field vector due to magnet magnetised in x

        Args:
            a ([type]): [description]
            b ([type]): [description]
            c ([type]): [description]
            Jr ([type]): [description]
            x ([type]): [description]
            y ([type]): [description]
            z ([type]): [description]

        Returns:
            [tuple]: Bx, By, Bz
        """
        a = self.a
        b = self.b
        c = self.c
        Jr = self.Jr

        Bx = self._calcBx_prism_x(a, b, c, Jr, x, y, z)
        By = self._calcBy_prism_x(a, b, c, Jr, x, y, z)
        Bz = self._calcBz_prism_x(a, b, c, Jr, x, y, z)
        return Bx, By, Bz

    def _calcB_prism_z(self, x, y, z):
        """Magnetic field vector due to magnet magnetised in y

        Args:
            a ([type]): [description]
            b ([type]): [description]
            c ([type]): [description]
            Jr ([type]): [description]
            x ([type]): [description]
            y ([type]): [description]
            z ([type]): [description]

        Returns:
            [tuple]: Bx, By, Bz
        """
        a = self.a
        b = self.b
        c = self.c
        Jr = self.Jr

        Bx = self._calcBz_prism_x(-c, b, a, Jr, -z, y, x)
        By = self._calcBy_prism_x(-c, b, a, Jr, -z, y, x)
        Bz = -1 * self._calcBx_prism_x(-c, b, a, Jr, -z, y, x)
        return Bx, By, Bz

    def _calcB_prism_y(self, x, y, z):
        """Magnetic field vector due to magnet magnetised in z

        Args:
            a ([type]): [description]
            b ([type]): [description]
            c ([type]): [description]
            Jr ([type]): [description]
            x ([type]): [description]
            y ([type]): [description]
            z ([type]): [description]

        Returns:
            [tuple]: Bx, By, Bz
        """
        a = self.a
        b = self.b
        c = self.c
        Jr = self.Jr

        Bx = self._calcBy_prism_x(-b, a, c, Jr, -y, x, z)
        By = -1 * self._calcBx_prism_x(-b, a, c, Jr, -y, x, z)
        Bz = self._calcBz_prism_x(-b, a, c, Jr, -y, x, z)
        return Bx, By, Bz

    def _generate_mask(self, x, y, z):
        """Generates mask of points inside a magnet

        Args:
            x (ndarray/float): x-coordinates
            y (ndarray/float): y-coordinates
            z (ndarray/float): z-coordinates
        """
        w, d, h = self.size()
        xn = -w / 2
        xp = w / 2
        yn = -d / 2
        yp = d / 2
        zn = -h / 2
        zp = h / 2

        mask_x = _np.logical_and(x > xn, x < xp)
        mask_y = _np.logical_and(y > yn, y < yp)
        mask_z = _np.logical_and(z > zn, z < zp)

        # merge logical masks
        mask = _np.logical_and(mask_x, mask_z)
        mask = _np.logical_and(mask, mask_y)

        return mask


class Cube(Prism):
    """Cube Magnet Class. Cuboid width x depth x height.

    Args:
          width [float]: magnet side length [default 10e-3]
          J: remnant magnetisation in Tesla [default 1.0]

     Optional Arguments:
          centre: magnet centre (tuple or Point3) [default Point3(.0, .0, .0)]
          theta: Angle between magnetisation and x-axis [default 90.0 degrees]
          phi: Angle between magnetisation and z-axis [default 0.0 degrees]

    """

    mag_type = "Cube"

    def __init__(
        self,
        width=10e-3,  # magnet dimensions
        Jr=1.0,  # local magnetisation direction
        **kwargs,
    ):

        super().__init__(width=width, depth=width, height=width, Jr=Jr, **kwargs)


class Cylinder(Magnet_3D):
    """Cylinder Magnet Class

    Args:
         radius: radius [default 10e-3]
         length: length [default 10e-3]
         J: remnant magnetisation in Tesla [default 1.0]


    Optional Arguments:
         centre: magnet centre (tuple or Point3) [default Point3(.0, .0, .0)]

    """

    mag_type = "Cylinder"

    def __init__(
        self,
        radius=10e-3,
        length=10e-3,  # magnet dimensions
        Jr=1.0,  # local magnetisation direction
        **kwargs,
    ):
        super().__init__(Jr, **kwargs)
        self.radius = radius
        self.length = length

        center = kwargs.pop("center", Point3(0.0, 0.0, 0.0))

        if type(center) is tuple:
            center = Point3(center[0], center[1], center[2])

        self.xc = center.x
        self.yc = center.y
        self.zc = center.z

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size()} (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.Jr} (T)\n"
            + f"Size: {self.size() } (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def size(self):
        """Returns magnet dimesions

        Returns:
            size[ndarray]: numpy array [radius, length]
        """
        return _np.array([self.radius, self.length])

    def get_Jr(self):
        return _np.array([0.0, 0.0, self.Jr])

    # Use Numba to create compiled Numpy ufunc that accepts arrays
    @staticmethod
    @vectorize([float64(float64, float64, float64, float64)])
    def _cel(kc, p, c, s):
        """Bulirsch's complete elliptic integral
        See NIST Handbook of Mathematical Functions, http://dlmf.nist.gov/19.2
        """
        if kc == 0:
            data = _np.NaN
            return data
        else:
            errtol = 0.000001
            k = _np.abs(kc)
            pp = p
            cc = c
            ss = s
            em = 1.0

            if p > 0:
                pp = _np.sqrt(p)
                ss = s / pp
            else:
                f = kc * kc
                q = 1.0 - f
                g = 1.0 - pp
                f = f - pp
                q = q * (ss - c * pp)
                pp = _np.sqrt(f / g)
                cc = (c - ss) / g
                ss = -q / (g * g * pp) + cc * pp
            f = cc
            cc = cc + ss / pp
            g = k / pp
            ss = 2 * (ss + f * g)
            pp = g + pp
            g = em
            em = k + em
            kk = k

            while _np.abs(g - k) > g * errtol:
                k = 2 * _np.sqrt(kk)
                kk = k * em
                f = cc
                cc = cc + ss / pp
                g = kk / pp
                ss = 2 * (ss + f * g)
                pp = g + pp
                g = em
                em = k + em
            data = (_np.pi / 2.0) * (ss + cc * em) / (em * (em + pp))
            return data

    def _calcB_local(self, x, y, z):
        """Internal magnetic field calculation methods.
        Calculates the field due to a cylindrical magnet/solenoid magnetised along z
        (in local coordinates).

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates

        Returns:
            Vector3: Magnetic field array
        """
        from ._routines import cart2pol, pol2cart
        from ._routines3 import _allocate_field_array3

        B = _allocate_field_array3(x, y, z)

        # Convert cartesian coordinates to cylindrical
        rho, phi = cart2pol(x, y)

        Brho, B.z = self._calcB_cyl(rho, z)

        # Convert magnetic fields from cylindrical to cartesian
        # We use pol2cart because Bphi is zero
        # If Bphi was != 0, then would have to use
        # `._routines.vector_pol2cart(Brho, Bphi, phi)`

        B.x, B.y = pol2cart(Brho, phi)

        return B

    def _calcB_cyl(self, rho, z):
        """Calculates the magnetic field due to a solenoid/cylinder in
        polar cylindrical coordinates

        Args:
            rho (float/array): radial coordinates
            z (float/array): axial coordinates

        Returns:
            tuple: Brho, Bz - magnetic field components
        """
        a = self.radius
        b = self.length / 2
        B0 = self.Jr / _np.pi

        zp = z + b
        zn = z - b

        zp_sq = _np.power(zp, 2)
        zn_sq = _np.power(zn, 2)
        rho_a_sq = _np.power(rho + a, 2)
        nrho_a_sq = _np.power(a - rho, 2)

        alphap = a / _np.sqrt(zp_sq + rho_a_sq)
        alphan = a / _np.sqrt(zn_sq + rho_a_sq)

        betap = zp / _np.sqrt(zp_sq + rho_a_sq)
        betan = zn / _np.sqrt(zn_sq + rho_a_sq)

        gamma = (a - rho) / (a + rho)

        kp = _np.sqrt((zp_sq + nrho_a_sq) / (zp_sq + rho_a_sq))

        kn = _np.sqrt((zn_sq + nrho_a_sq) / (zn_sq + rho_a_sq))

        Brho = B0 * (
            alphap * self._cel(kp, 1, 1, -1) - alphan * self._cel(kn, 1, 1, -1)
        )

        Bz = (B0 * a / (a + rho)) * (
            betap * self._cel(kp, gamma ** 2, 1, gamma)
            - betan * self._cel(kn, gamma ** 2, 1, gamma)
        )
        return Brho, Bz

    def _generate_mask(self, x, y, z):
        """Generates mask of points inside a cylindrical magnet

        Args:
            x (ndarray/float): x-coordinates
            y (ndarray/float): y-coordinates
            z (ndarray/float): z-coordinates
        """

        radius, length = self.size()
        data_norm = x ** 2 + y ** 2
        zn = -length / 2
        zp = length / 2

        mask_rho = data_norm < radius ** 2
        mask_z = _np.logical_and(z > zn, z < zp)
        mask = _np.logical_and(mask_rho, mask_z)

        return mask


class Sphere(Magnet_3D):
    """Sphere Magnet Class

    Args:
         radius: radius [default 10e-3]
         J: remnant magnetisation in Tesla [default 1.0]


    Optional Arguments:
         centre: magnet centre (tuple or Point3) [default Point3(.0, .0, .0)]
         theta: Angle between magnetisation and x-axis [default 90.0 degrees]
         phi: Angle between magnetisation and z-axis [default 0.0 degrees]
    """

    mag_type = "Sphere"

    def __init__(
        self,
        radius=10e-3,
        Jr=1.0,  # local magnetisation direction
        **kwargs,
    ):

        super().__init__(Jr, **kwargs)
        self.radius = radius

        # self.phi = kwargs.pop("phi", 90)
        # self.phi_rad = _np.deg2rad(self.phi)
        # self.theta = kwargs.pop("theta", 0)
        # self.theta_rad = _np.deg2rad(self.theta)

        # self.Jx = _np.around(
        #     Jr * _np.cos(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        # )
        # self.Jy = _np.around(
        #     Jr * _np.sin(self.phi_rad) * _np.sin(self.theta_rad), decimals=6
        # )
        # self.Jz = _np.around(Jr * _np.cos(self.theta_rad), decimals=6)

        center = kwargs.pop("center", Point3(0.0, 0.0, 0.0))

        if type(center) is tuple:
            center = Point3(center[0], center[1], center[2])

        self.xc = center.x
        self.yc = center.y
        self.zc = center.z

    def __str__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size()} (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"{self.__class__.mag_type}\n"
            + f"J: {self.get_Jr()} (T)\n"
            + f"Size: {self.size() } (m)\n"
            + f"Center {self.center()} (m)\n"
            + f"Orientation alpha,beta,gamma: {self.get_orientation()}\n"
        )
        return str

    def size(self):
        """Returns magnet dimesions

        Returns:
            size[ndarray]: numpy array [radius, length]
        """
        return _np.array([self.radius])

    def get_Jr(self):
        return _np.array([0.0, 0.0, self.Jr])

    def _calcB_local(self, x, y, z):
        """Internal magnetic field calculation methods.
        Calculates the field due to a spherical magnet magnetised along z
        (in local coordinates). Returns a dipolar field outside the magnet

        Args:
            x (float/array): x co-ordinates
            y (float/array): y co-ordinates
            z (float/array): z co-ordinates

        Returns:
            Vector3: Magnetic field array
        """
        from ._routines import cart2sph, sphere_sph2cart
        from ._routines3 import _allocate_field_array3

        B = _allocate_field_array3(x, y, z)

        # Convert to spherical coordinates
        r, theta, phi = cart2sph(x, y, z)

        # Calculates field for sphere magnetised along z
        Br, Btheta = self._calcB_spherical(r, theta)

        # Convert magnetic fields from spherical to cartesian
        B.x, B.y, B.z = sphere_sph2cart(Br, Btheta, theta, phi)
        return B

    # FIXME: Phi not rotating field
    # def _rotate_and_calcB(self, x, y, z):
    #     """Translates and rotates global coordinates into local before calculating
    #     the magnetic field and then rotating magnetic field into global frame

    #     Args:
    #         x (array): x coordinates
    #         y (array): y coordinates
    #         z (array): z coordinates

    #     TODO:
    #         - Fix the rotation around z axis for phi (local magnetisation)

    #     Returns:
    #         tuple: Bx, By, Bz: calculated magnetic field
    #     """

    #     from ._routines import cart2sph, sphere_sph2cart
    #     from functools import reduce

    #     x = Quaternion._input_to_numpy(x)
    #     y = Quaternion._input_to_numpy(y)
    #     z = Quaternion._input_to_numpy(z)

    #     rotate_about_y, rotate_about_z = None, None
    #     forward_rotation, reverse_rotation = None, None

    #     if _np.fabs(self.theta_rad) > self.tol:
    #         rotate_about_y = Quaternion.q_angle_from_axis(self.theta_rad, (0, 1, 0))

    #     if _np.fabs(self.phi_rad) > self.tol:
    #         rotate_about_z = Quaternion.q_angle_from_axis(self.phi_rad, (0, 0, 1))

    #     # Generate compound rotations
    #     q_forward = [q for q in [rotate_about_z, rotate_about_y] if q is not None]
    #     forward_rotation = reduce(lambda x, y: x * y, q_forward)

    #     q_reverse = [
    #         q.get_conjugate() for q in [rotate_about_y, rotate_about_z] if q is not None
    #     ]
    #     reverse_rotation = reduce(lambda x, y: x * y, q_reverse)

    #     # Generate 3xN array for quaternion rotation
    #     pos_vec = Quaternion._prepare_vector(x, y, z)

    #     # Rotate points
    #     x_rot, y_rot, z_rot = forward_rotation * pos_vec

    #     # Convert points to spherical coordinates
    #     r, theta, phi = cart2sph(x_rot, y_rot, z_rot)

    #     # Calculates field for sphere magnetised along z
    #     Br, Btheta = self._calcB_spherical(r, theta)

    #     # Convert magnetic fields from spherical to cartesian
    #     Bx, By, Bz = sphere_sph2cart(Br, Btheta, theta, phi)

    #     # Rotate using quaternions to get correct fields
    #     Bvec = Quaternion._prepare_vector(Bx, By, Bz)
    #     Bx_rot, By_rot, Bz_rot = reverse_rotation * Bvec
    #     return Bx_rot, By_rot, Bz_rot

    def _calcB_spherical(self, r, theta):
        """Calculates the magnetic field due to due to a sphere at any point
        in spherical coordinates

        Args:
            r (float/array): radial coordinates
            theta (float/array): azimuthal coordinates

        Returns:
            tuple: Br, Btheta
        """

        # Hide the warning for situtations where there is a divide by zero.
        # This returns a NaN in the array, which is ignored for plotting.
        with _np.errstate(divide="ignore", invalid="ignore"):
            preFac = self.Jr * (self.radius ** 3 / r ** 3) / 3.0

        Br = preFac * 2.0 * _np.cos(theta)
        Btheta = preFac * _np.sin(theta)

        return Br, Btheta

    def _generate_mask(self, x, y, z):
        """Generates mask of points inside a spherical magnet

        Args:
            x (ndarray/float): x-coordinates
            y (ndarray/float): y-coordinates
            z (ndarray/float): z-coordinates
        """

        data_norm = x ** 2 + y ** 2 + z ** 2
        mask = data_norm < self.radius ** 2

        return mask
