# This Source Code Form is subject to the terms of the Mozilla Public
# License, v.2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnets.utils._vector_structs

Private module consiting of vector and point array classes and their methods.
"""

__all__ = ["Field1", "Field2", "Field3", "Point_Array2", "Point_Array3"]

import numpy as _np
from ._conversions import get_unit_value_meter, get_unit_value_tesla


class Point_Array1(object):
    """2D vector class consisting of numpy arrays of x and y coordinates

    Args:
        x:
        y:

    Methods:
        calc_norm: calculates the magnitude of the fields at every point and
                    stores in self.n
    """

    def __init__(self, z, unit="m"):
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
        return self.unit

    def change_unit(self, new_unit, get_unit_value=get_unit_value_meter):
        """Converts vector from one unit scale to another

        Args:
            new_unit (string): SI prefix, eg cm, mm
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
    """2D vector class consisting of numpy arrays of x and y coordinates

    Args:
        x:
        y:

    Methods:
        calc_norm: calculates the magnitude of the fields at every point and
                    stores in self.n
    """

    def __init__(self, x, y, unit="m"):
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
        """Converts vector from one unit scale to another

        Args:
            new_unit (string): SI prefix, eg cm, mm
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
    def __init__(self, x, y, z, unit="m"):
        super().__init__(x, y, unit=unit)
        self.z = _np.asarray(z)

    def __repr__(self) -> str:
        return f"[Unit: {self.unit}\nx: {self.x}\ny: {self.y}\nz: {self.z}]"

    def __str__(self) -> str:
        return f"[Unit: {self.unit}\nx: {self.x}\ny: {self.y}\nz: {self.z}]"

    def change_unit(self, new_unit, get_unit_value=get_unit_value_meter):
        """Converts vector from one unit scale to another

        Args:
            new_unit (string): SI prefix, eg cm, mm
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
    """1D Field vector class consisting of numpy arrays of z values

    Args:
        x:
        y:

    Methods:
        calc_norm: calculates the magnitude of the fields at every point and
                    stores in self.n
    """

    def __init__(self, z, unit="T"):
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
        super().change_unit(new_unit, get_unit_value)


class Field2(Point_Array2):
    """2D Field vector class consisting of numpy arrays of x and y values

    Args:
        x:
        y:

    Methods:
        calc_norm: calculates the magnitude of the fields at every point and
                    stores in self.n
    """

    def __init__(self, x, y, unit="T"):
        super().__init__(x, y)
        self.n = _np.zeros_like(x)
        if get_unit_value_tesla(unit) is not None:
            self.unit = unit
        else:
            raise ValueError("Error, not an SI prefix, e.g. T, mT, uT, nT ")

    def calc_norm(self):
        self.n = _np.linalg.norm([self.x, self.y], axis=0)

    def __repr__(self) -> str:
        return f"[Unit: {self.unit}\nBx: {self.x}\nBy: {self.y}\nBn: {self.n}]"

    def __str__(self) -> str:
        return f"[Unit: {self.unit}\nBx: {self.x}\nBy: {self.y}\nBn: {self.n}]"

    def change_unit(self, new_unit, get_unit_value=get_unit_value_tesla):
        super().change_unit(new_unit, get_unit_value)


class Field3(Point_Array3):
    """3D Field vector class consisting of numpy arrays of x and y values

    Args:
        x:
        y:
        z:

    Methods:
        calc_norm: calculates the magnitude of the fields at every point and
                    stores in self.n
    """

    def __init__(self, x, y, z, unit="T"):
        super().__init__(x, y, z)
        self.n = _np.zeros_like(x)

    def calc_norm(self):
        self.n = _np.linalg.norm([self.x, self.y, self.z], axis=0)

    def change_unit(self, new_unit, get_unit_value=get_unit_value_tesla):
        super().change_unit(new_unit, get_unit_value)

    def __repr__(self) -> str:
        return f"[Unit: {self.unit}\nBx: {self.x}\nBy: {self.y}\nBz: {self.z}\nBn: {self.n}]"

    def __str__(self) -> str:
        return f"[Unit: {self.unit}\nBx: {self.x}\nBy: {self.y}\nBz: {self.z}\nBn: {self.n}]"
