# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Three Dimensional Magnet Classes
"""
import numpy as _np
from ._fields import Point3
from ._magnet import Magnet
# from pymagma import PI

__all__ = ['Magnet_3D', 'Prism', 'Cube', 'Cylinder']

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
    mag_type = 'Magnet_3D'

    def __init__(self, width, depth, height, Jr, **kwargs) -> None:
        super().__init__()
        self.width = width
        self.depth = depth
        self.height = height

        self.a = width/2
        self.b = depth/2
        self.c = height/2

        self.Jr = Jr 

        center = kwargs.pop('center', Point3(0.0, 0.0, 0.0))

        if type(center) is tuple:
            center = Point3(center[0], center[1], center[2])

        self.xc = center.x
        self.yc = center.y
        self.zc = center.z

    def center(self):
        return _np.array([self.xc, self.yc, self.zc])
    
    def size(self):
        """Returns magnet dimesions

        Returns:
            size[ndarray]: numpy array [width, depth, height]
        """
        return _np.array([self.width, self.depth, self.height])



class Prism(Magnet_3D):
    """Prism Magnet Class. Cuboid width x depth x height.

   Args:
         width: magnet width [default 10e-3]
         depth: magnet depth [default 20e-3]
         height: magnet height [default 30e-3]
         J: remnant magnetisation in Tesla [default 1.0]

    Optional Arguments:
         centre: magnet centre (tuple or Point3) [default Point3(.0, .0, .0)]
         theta: Angle between magnetisation and x-axis [default 90.0 degrees]
         phi: Angle between magnetisation and z-axis [default 0.0 degrees]  

    """
    mag_type = 'Prism'
    _instances = []
    counter = 0

    def __init__(self,
                 width=10e-3, depth=20e-3, height=30e-3,   # magnet dimensions
                 Jr=1.0,  # local magnetisation direction
                 **kwargs):

        super().__init__(width, depth, height, Jr, **kwargs)

        self.theta = kwargs.pop('theta', 90)
        self.theta_rad = _np.deg2rad(self.theta)
        self.phi = kwargs.pop('phi', 0)
        self.phi_rad = _np.deg2rad(self.phi)

        self.Jx = _np.around(Jr*_np.cos(self.theta_rad)*_np.sin(self.phi_rad),
                             decimals=6)
        self.Jy = _np.around(Jr*_np.sin(self.theta_rad)*_np.sin(self.phi_rad),
                             decimals=6)
        self.Jz = _np.around(Jr*_np.cos(self.phi_rad),
                             decimals=6)
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
              f'Size: {self.size() * 2} (m)\n' + \
              f'Center {self.center()} (m)\n'
        return str

    def size(self):
        return _np.array([self.a, self.b, self.c])

    def get_Jr(self):
        return _np.array([self.Jx, self.Jy, self.Jz])

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
            data = _np.arctan(
                ((y + b) * (z + c)) /
                ((x + a) * _np.sqrt(_np.power((x + a), 2) +
                                    _np.power((y + b), 2) +
                                    _np.power((z + c), 2)))
            )
        except ValueError:
            data = _np.NaN
        return data

    @staticmethod
    def _F2(a, b, c, x, y, z):
        """Helper Function F2 for 3D prismatic magnets

        Args:
            x ([type]): [description]
            y ([type]): [description]
            z ([type]): [description]

        Returns:
            [array]: [description]
        """

        try:
            xa_sq = _np.power((x + a), 2)
            yb_sq = _np.power((y + b), 2)
            zc_sq = _np.power((z + c), 2)
            znc_sq = _np.power((z - c), 2)

            data = (_np.sqrt(xa_sq + yb_sq + znc_sq) + c - z) / (
                _np.sqrt(xa_sq + yb_sq + zc_sq) - c - z)
        except ValueError:
            data = _np.NaN
        return data

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
            data = -(Jr/(4 * _np.pi)) * (
                self._F1(a, b, c, -x, y, z) + self._F1(a, b, c, -x, y, -z)
                + self._F1(a, b, c, -x, -y, z) + self._F1(a, b, c, -x, -y, -z)
                + self._F1(a, b, c, x, y, z) + self._F1(a, b, c, x, y, -z)
                + self._F1(a, b, c, x, -y, z) + self._F1(a, b, c, x, -y, -z)
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
        # a = self.a
        # b = self.b
        # c = self.c
        # Jr = self.Jr

        try:
            data = _np.log(self._F2(a, b, c, -x, -y, z) *
                           self._F2(a, b, c, x, y, z) /
                           (self._F2(a, b, c, -x, y, z) *
                            self._F2(a, b, c, x, -y, z)))
            data *= (Jr / (4 * _np.pi))
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
            data = _np.log(self._F2(a, c, b, -x, -z, y) *
                           self._F2(a, c, b, x, z, y) /
                           (self._F2(a, c, b, -x, z, y) *
                            self._F2(a, c, b, x, -z, y)))
            data *= (Jr/(4 * _np.pi))
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
        Bz = -1*self._calcBx_prism_x(-c, b, a, Jr, -z, y, x)
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
        By = -1*self._calcBx_prism_x(-b, a, c, Jr, -y, x, z)
        Bz = self._calcBz_prism_x(-b, a, c, Jr, -y, x, z)
        return Bx, By, Bz


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
    mag_type = 'Cube'

    def __init__(self,
                 width=10e-3,   # magnet dimensions
                 Jr=1.0,  # local magnetisation direction
                 **kwargs):

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
    mag_type = 'Cylinder'

    def __init__(self,
                 radius=10e-3, length=10e-3,   # magnet dimensions
                 Jr=1.0,  # local magnetisation direction
                 **kwargs):

        self.radius = radius
        self.length = length
        self.Jr = Jr
        self.Jx = 0
        self.Jy = 0
        self.Jz = Jr

        center = kwargs.pop('center', Point3(0.0, 0.0, 0.0))

        if type(center) is tuple:
            center = Point3(center[0], center[1], center[2])

        self.xc = center.x
        self.yc = center.y
        self.zc = center.z

    def __str__(self):
        str = f'{self.__class__.mag_type}\n' + \
              f'J: {self.Jr} (T)\n' + \
              f'Size: {self.size()} (m)\n' + \
              f'Center {self.center()} (m)\n'
        return str

    def __repr__(self):
        str = f'{self.__class__.mag_type}\n' + \
              f'J: {self.Jr} (T)\n' + \
              f'Size: {self.size() } (m)\n' + \
              f'Center {self.center()} (m)\n'
        return str

    def size(self):
        """Returns magnet dimesions

        Returns:
            size[ndarray]: numpy array [radius, length]
        """
        return _np.array([self.radius, self.length])

    @staticmethod
    @_np.vectorize
    def _cel(kc, p, c, s):
        """ Burlisch's complete elliptic integral
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
                ss = s/pp
            else:
                f = kc*kc
                q = 1.0 - f
                g = 1.0 - pp
                f = f - pp
                q = q*(ss - c*pp)
                pp = _np.sqrt(f/g)
                cc = (c - ss)/g
                ss = - q/(g*g*pp) + cc*pp
            f = cc
            cc = cc + ss/pp
            g = k/pp
            ss = 2*(ss + f*g)
            pp = g + pp
            g = em
            em = k + em
            kk = k

            while _np.abs(g - k) > g * errtol:
                k = 2*_np.sqrt(kk)
                kk = k*em
                f = cc
                cc = cc + ss/pp
                g = kk/pp
                ss = 2*(ss + f*g)
                pp = g + pp
                g = em
                em = k + em
            data = (_np.pi / 2.) * (ss + cc * em) / (em * (em + pp))
            return data

    def _calcB_cyl(self, rho, z):
        """ Calculates the magnetic field due to 
            due to a solenoid at any point
            returns Bz,Br
            Omit B0 for B0 = 1 mT
        """
        a = self.radius
        b = self.length / 2
        B0 = self.Jr

        # rho -= self.xc
        # z -= self.zc

        zp = z + b
        zn = z - b

        zp_sq = _np.power(zp, 2)
        zn_sq = _np.power(zn, 2)
        rho_a_sq = _np.power(rho + a, 2)
        nrho_a_sq = _np.power(a - rho, 2)

        alphap = a/_np.sqrt(zp_sq + rho_a_sq)
        alphan = a/_np.sqrt(zn_sq + rho_a_sq)

        betap = zp/_np.sqrt(zp_sq + rho_a_sq)
        betan = zn/_np.sqrt(zn_sq + rho_a_sq)

        gamma = (a - rho) / (a + rho)

        kp = _np.sqrt((zp_sq + nrho_a_sq) / (zp_sq + rho_a_sq))

        kn = _np.sqrt((zn_sq + nrho_a_sq) / (zn_sq + rho_a_sq))

        Br = B0 * (alphap * self._cel(kp, 1, 1, -1)
                   - alphan * self._cel(kn, 1, 1, -1))

        Bz = (B0*a/(a + rho)) * (betap * self._cel(kp, gamma**2, 1, gamma)
                                 - betan * self._cel(kn, gamma**2, 1, gamma))

        return Bz / _np.pi, Br / _np.pi

