# This Source Code Form is subject to the terms of the Mozilla Public
# License, v.2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnets.utils._vector_structs

Private module consiting of vector and point array classes and their methods.
"""

import numpy as _np
from ._conversions import get_unit_value_meter, get_unit_value_tesla

__all__ = ["Field1", "Field2", "Field3", "Point_Array2", "Point_Array3"]


class Point_Array1(object):
    """1D point structure
    This is used to contain one position array (z), and the units
    ('mm', 'cm', etc)
    """

    def __init__(self, z, unit="mm"):
        """Init method

        Args:
            z (ndarray): z coordinates
            unit (str, optional): Unit of length. Defaults to "mm".

        Raises:
            ValueError: Unit must an SI prefix, e.g. km, m, cm, mm
        """
        self.z = _np.asarray(z)
        if get_unit_value_meter(unit) is not None:
            self.unit = unit
        else:
            raise ValueError("Error, not an SI prefix, e.g. km, m, cm, mm ")

    def __repr__(self) -> str:
        return f"[(Unit:{self.unit}) Array:{self.z}]"

    def __str__(self) -> str:
        return f"[(Unit:{self.unit}) Array:{self.z}]"

    def get_unit(self):
        """Gets unit

        Returns:
            str: unit
        """
        return self.unit

    def change_unit(self, new_unit, get_unit_value=get_unit_value_meter):
        """Converts point array to a different unit. e.g from 'cm' to 'mm'

        Args:
            new_unit (str): unit to be converted to
            get_unit_value (function, optional): Function for checking unit type. Defaults to get_unit_value_meter.
        """
        from ..magnets import Magnet, Prism, Cylinder

        current_unit_val = get_unit_value(self.get_unit())
        new_unit_val = get_unit_value(new_unit)
        scale_val = current_unit_val / new_unit_val

        self.z *= scale_val

        self.unit = new_unit
        for magnet in Magnet.instances:
            magnet.center = scale_val * magnet.center
            if issubclass(magnet.__class__, Prism):
                magnet.a = magnet.a * scale_val
                magnet.b = magnet.b * scale_val
                magnet.c = magnet.c * scale_val
                magnet.width = magnet.width * scale_val
                magnet.depth = magnet.depth * scale_val
                magnet.height = magnet.height * scale_val
            elif issubclass(magnet.__class__, Cylinder):
                magnet.radius = magnet.radius * scale_val
                magnet.length = magnet.length * scale_val


class Point_Array2(object):
    """2D point structure
    This is used to contain two position arrays (x, y), and the units
    ('mm', 'cm', etc)
    """

    def __init__(self, x, y, unit="mm"):
        """Init Method

        Args:
            x (ndarray): x coordinates
            y (ndarray): y coordinates
            unit (str, optional): Unit of length. Defaults to "mm".

        Raises:
            ValueError: Unit must an SI prefix, e.g. km, m, cm, mm
        """
        self.x = _np.asarray(x)
        self.y = _np.asarray(y)
        if get_unit_value_meter(unit) is not None:
            self.unit = unit
        else:
            raise ValueError("Error, not an SI prefix, e.g. km, m, cm, mm ")

    def __repr__(self) -> str:
        return f"[Unit: {self.unit}\nx: {self.x}\ny: {self.y}]"

    def __str__(self) -> str:
        return f"[Unit: {self.unit}\nx: {self.x}\ny: {self.y}]"

    def get_unit(self):
        return self.unit

    def change_unit(self, new_unit, get_unit_value=get_unit_value_meter):
        """Converts point array to a different unit. e.g from 'cm' to 'mm'

        Args:
            new_unit (str): unit to be converted to
            get_unit_value (function, optional): Function for checking unit type. Defaults to get_unit_value_meter.
        """
        from ..magnets import Magnet, Rectangle, Square, Circle, PolyMagnet

        current_unit_val = get_unit_value(self.get_unit())
        new_unit_val = get_unit_value(new_unit)
        scale_val = current_unit_val / new_unit_val

        self.x *= scale_val
        self.y *= scale_val
        self.unit = new_unit
        for magnet in Magnet.instances:
            magnet.center = scale_val * magnet.center
            if issubclass(magnet.__class__, Rectangle):
                magnet.a = magnet.a * scale_val
                magnet.b = magnet.b * scale_val
                magnet.width = magnet.width * scale_val
                magnet.height = magnet.height * scale_val
            elif issubclass(magnet.__class__, Square):
                magnet.a = magnet.a * scale_val
                magnet.width = magnet.width * scale_val
            elif issubclass(magnet.__class__, Circle):
                magnet.radius = magnet.radius * scale_val
            elif issubclass(magnet.__class__, PolyMagnet):
                magnet.polygon.vertices = (
                    scale_val * _np.array(magnet.polygon.vertices)
                ).tolist()


class Point_Array3(Point_Array2):
    """3D point structure
    This is used to contain three position arrays (x, y, z), and the units
    ('mm', 'cm', etc)
    """

    def __init__(self, x, y, z, unit="mm"):
        """Init Method

        Args:
            x (ndarray): x coordinates
            y (ndarray): y coordinates
            z (ndarray): z coordinates
            unit (str, optional): Unit of length. Defaults to "mm".

        Raises:
            ValueError: Unit must an SI prefix, e.g. km, m, cm, mm
        """
        super().__init__(x, y, unit=unit)
        self.z = _np.asarray(z)

    def __repr__(self) -> str:
        return f"[Unit: {self.unit}\nx: {self.x}\ny: {self.y}\nz: {self.z}]"

    def __str__(self) -> str:
        return f"[Unit: {self.unit}\nx: {self.x}\ny: {self.y}\nz: {self.z}]"

    def change_unit(self, new_unit, get_unit_value=get_unit_value_meter):
        """Converts point array to a different unit. e.g from 'cm' to 'mm'

        Args:
            new_unit (str): unit to be converted to
            get_unit_value (function, optional): Function for checking unit type. Defaults to get_unit_value_meter.
        """
        from ..magnets import Magnet, Prism, Cube, Cylinder, Sphere, Mesh

        current_unit_val = get_unit_value(self.get_unit())
        new_unit_val = get_unit_value(new_unit)
        scale_val = current_unit_val / new_unit_val
        self.x *= scale_val
        self.y *= scale_val
        self.z *= scale_val
        self.unit = new_unit

        for magnet in Magnet.instances:
            magnet.center = scale_val * magnet.center
            if issubclass(magnet.__class__, Prism):
                magnet.a = magnet.a * scale_val
                magnet.b = magnet.b * scale_val
                magnet.c = magnet.c * scale_val
                magnet.width = magnet.width * scale_val
                magnet.height = magnet.height * scale_val
                magnet.depth = magnet.depth * scale_val
            elif issubclass(magnet.__class__, Cube):
                magnet.a = magnet.a * scale_val
                magnet.width = magnet.width * scale_val
            elif issubclass(magnet.__class__, Cylinder):
                magnet.radius = magnet.radius * scale_val
                magnet.length = magnet.length * scale_val
            elif issubclass(magnet.__class__, Sphere):
                magnet.radius = magnet.radius * scale_val
            elif issubclass(magnet.__class__, Mesh):
                magnet.mesh_vectors = magnet.mesh_vectors * scale_val


class Field1(Point_Array1):
    """1D Field vector
    This is used to contain one component (Bz as z), and the units
    ('T', 'mT', etc)
    """

    def __init__(self, z, unit="T"):
        """Init method

        Args:
            z (ndarray): Magnetic field component
            unit (str, optional): Unit of field. Defaults to "T".

        Raises:
            ValueError: Unit must an SI prefix, e.g. T, mT, uT, nT
        """
        super().__init__(z)
        if get_unit_value_tesla(unit) is not None:
            self.unit = unit
        else:
            raise ValueError("Error, not an SI prefix, e.g. T, mT, uT, nT ")

    def __repr__(self) -> str:
        return f"[Unit: {self.unit}\nBz: {self.z}]"

    def __str__(self) -> str:
        return f"[Unit: {self.unit}\nBz: {self.z}]"

    def change_unit(self, new_unit, get_unit_value=get_unit_value_tesla):
        """Converts field array to a different unit. e.g from 'T' to 'mT'

        Args:
            new_unit (str): unit to be converted to
            get_unit_value (function, optional): Function for checking unit type. Defaults to get_unit_value_tesla.
        """
        super().change_unit(new_unit, get_unit_value)


class Field2(Point_Array2):
    """2D Field vector
    This is used to contain two components (x, y), and the units
    ('T', 'mT', etc)
    """

    def __init__(self, x, y, unit="T"):
        """Init method

        Args:
            x (ndarray): Magnetic field component Bx
            y (ndarray): Magnetic field component By
            unit (str, optional): Unit of field. Defaults to "T".

        Raises:
            ValueError: Unit must an SI prefix, e.g. T, mT, uT, nT
        """
        super().__init__(x, y)
        self.n = _np.zeros_like(x)
        if get_unit_value_tesla(unit) is not None:
            self.unit = unit
        else:
            raise ValueError("Error, not an SI prefix, e.g. T, mT, uT, nT ")

    def calc_norm(self):
        """Calculates the norm of the 2D vector"""
        self.n = _np.linalg.norm([self.x, self.y], axis=0)

    def __repr__(self) -> str:
        return f"[Unit: {self.unit}\nBx: {self.x}\nBy: {self.y}\nBn: {self.n}]"

    def __str__(self) -> str:
        return f"[Unit: {self.unit}\nBx: {self.x}\nBy: {self.y}\nBn: {self.n}]"

    def change_unit(self, new_unit, get_unit_value=get_unit_value_tesla):
        """Converts field array to a different unit. e.g from 'T' to 'mT'

        Args:
            new_unit (str): unit to be converted to
            get_unit_value (function, optional): Function for checking unit type. Defaults to get_unit_value_tesla.
        """
        super().change_unit(new_unit, get_unit_value)


class Field3(Point_Array3):
    """2D Field vector
    This is used to contain three components (x, y), and the units
    ('T', 'mT', etc)
    """

    def __init__(self, x, y, z, unit="T"):
        """Init method

        Args:
            x (ndarray): Magnetic field component Bx
            y (ndarray): Magnetic field component By
            z (ndarray): Magnetic field component Bz
            unit (str, optional): Unit of field. Defaults to "T".

        Raises:
            ValueError: Unit must an SI prefix, e.g. T, mT, uT, nT
        """
        super().__init__(x, y, z)
        self.n = _np.zeros_like(x)
        self.unit = unit

    def calc_norm(self):
        """Calculates the norm of the 3D vector"""
        self.n = _np.linalg.norm([self.x, self.y, self.z], axis=0)

    def change_unit(self, new_unit, get_unit_value=get_unit_value_tesla):
        """Converts field array to a different unit. e.g from 'T' to 'mT'

        Args:
            new_unit (str): unit to be converted to
            get_unit_value (function, optional): Function for checking unit type. Defaults to get_unit_value_tesla.
        """
        super().change_unit(new_unit, get_unit_value)

    def __repr__(self) -> str:
        return f"[Unit: {self.unit}\nBx: {self.x}\nBy: {self.y}\nBz: {self.z}\nBn: {self.n}]"

    def __str__(self) -> str:
        return f"[Unit: {self.unit}\nBx: {self.x}\nBy: {self.y}\nBz: {self.z}\nBn: {self.n}]"
