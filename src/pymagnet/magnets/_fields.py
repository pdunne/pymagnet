# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne

import numpy as _np

__all__ = ['Point2', 'Point3', 'Vector2', 'Vector3']


class Point2(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self) -> str:
        return f"({self.x}, {self.y})"

    def __str__(self) -> str:
        return f"({self.x}, {self.y})"

    def __add__(self, other):
        return Point2(self.x + other.x, self.y + other.y)
    
    def __sub__(self, other):
        return Point2(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return Point2(self.x * other.x, self.y * other.y)
    
    def __div__(self, other):
        x = (self.x / other.x) if other.x != 0 else 0
        y = (self.y / other.y) if other.y != 0 else 0
        return Point2(x, y)

    def __lt__(self, other):
        self_mag = (self.x ** 2) + (self.y ** 2)
        other_mag = (other.x ** 2) + (other.y ** 2)
        return self_mag < other_mag

    def __le__(self, other):
        self_mag = (self.x ** 2) + (self.y ** 2)
        other_mag = (other.x ** 2) + (other.y ** 2)
        return self_mag <= other_mag
    
    def __gt__(self, other):
        self_mag = (self.x ** 2) + (self.y ** 2)
        other_mag = (other.x ** 2) + (other.y ** 2)
        return self_mag > other_mag

    def __ge__(self, other):
        self_mag = (self.x ** 2) + (self.y ** 2)
        other_mag = (other.x ** 2) + (other.y ** 2)
        return self_mag >= other_mag

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return self.x != other.x or self.y != other.y

    def distance_to(self, point):
        return _np.sqrt(_np.power(point.x - self.x, 2) +
                        _np.power(point.y - self.y, 2))

    def distance_to_origin(self):
        return self._norm()
    
    def _norm(self):
        return _np.linalg.norm([self.x, self.y], axis=0)

    def move(self, point):
        try:
            if type(point) is tuple:
                self.x += point[0]
                self.y += point[1]
            else:
                self + point 
        except ValueError as e:
            print("Invalid input, should be a 2-element tuple or Point2 object")




class Point3(Point2):
    def __init__(self, x, y, z):
        super().__init__(x, y)
        self.z = z

    def __repr__(self) -> str:
        return f"({self.x}, {self.y}, {self.z})"

    def __str__(self) -> str:
        return f"({self.x}, {self.y}, {self.z})"
    
    def __add__(self, other):
        return Point3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Point3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        return Point3(self.x * other.x, self.y * other.y, self.z * other.z)

    def __div__(self, other):
        x = (self.x / other.x) if other.x != 0 else 0
        y = (self.y / other.y) if other.y != 0 else 0
        z = (self.z / other.z) if other.z != 0 else 0
        return Point3(x, y, z)

    def __lt__(self, other):
        self_mag = (self.x ** 2) + (self.y ** 2)
        other_mag = (other.x ** 2) + (other.y ** 2)
        return self_mag < other_mag

    def __le__(self, other):
        self_mag = (self.x ** 2) + (self.y ** 2)
        other_mag = (other.x ** 2) + (other.y ** 2)
        return self_mag <= other_mag

    def __gt__(self, other):
        self_mag = (self.x ** 2) + (self.y ** 2)
        other_mag = (other.x ** 2) + (other.y ** 2)
        return self_mag > other_mag

    def __ge__(self, other):
        self_mag = (self.x ** 2) + (self.y ** 2)
        other_mag = (other.x ** 2) + (other.y ** 2)
        return self_mag >= other_mag

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return self.x != other.x or self.y != other.y

    def distance_to(self, point):
        return _np.sqrt(_np.power(point.x - self.x, 2) +
                        _np.power(point.y - self.y, 2) +
                        _np.power(point.z - self.z, 2))
    
    def _norm(self):
        return _np.linalg.norm([self.x, self.y, self.z], axis=0)

    def distance_to_origin(self):
        return self._norm()


class Vector2(object):
    """2D Field vector class consisting of numpy arrays of x and y values

    Args:
        x:
        y:
    
    Methods:
        calc_norm: calculates the magnitude of the fields at every point and
                    stores in self.n
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.n = _np.zeros_like(x)

    def calc_norm(self):
        self.n = _np.linalg.norm([self.x, self.y], axis=0)
    
    # def __repr__(self) -> str:
    #     return f"({self.x}, {self.y}, {self.n})"
        
    def __str__(self) -> str:
        return f"({self.x}, {self.y}, {self.n})"


class Vector3(Vector2):
    """3D Field vector class consisting of numpy arrays of x and y values

    Args:
        x:
        y:
        z:
    
    Methods:
        calc_norm: calculates the magnitude of the fields at every point and
                    stores in self.n
    """

    def __init__(self, x, y, z):
        super().__init__(x, y)
        self.z = z

    def calc_norm(self):
        self.n = _np.linalg.norm([self.x, self.y, self.z], axis=0)
    
    # def __repr__(self) -> str:
    #     return f"({self.x}, {self.y}, {self.z}, {self.n})"

    def __str__(self) -> str:
        return f"({self.x}, {self.y}, {self.z}, {self.n})"











       





