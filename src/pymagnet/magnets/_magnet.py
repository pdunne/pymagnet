# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""Magnet Base class

This private module implements the registry and base magnet classes
"""
# __all__ = ['Magnet', 'reset_magnets', 'list_magnets']
__all__ = ["Magnet"]

from typing import List, Any, Union, Type, Optional
from weakref import WeakSet


class Registry:
    """Registry class for tracking instances

    Instances are tracked in `class.instances` using weak references.
    This also includes any instances that are deleted manually or go out of
    scope.

    Class methods:

        `print_instances()`

        `get_instances()`

        `get_num_instances()`

        `reset()
    """

    instances = WeakSet()
    _class_instances = []

    def __new__(cls, *args, **kwargs):
        o = object.__new__(cls)
        cls._register_instance(o)
        return o

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.__class__._class_instances.append(self)

    def __del__(self):
        if self in self.__class__._class_instances:
            self.__class__._class_instances.remove(self)

    @classmethod
    def print_instances(cls):
        """Prints class instantiations"""
        if len(cls.instances) < 1:
            print("No Instances")
        else:
            for instance in cls.instances:
                print(instance)

    @classmethod
    def get_instances(cls):
        return cls.instances

    @classmethod
    def get_num_instances(cls, Print_Val=False):
        """Return number of instances of class

        Args:
            Print_Val (bool, optional): [Print to screen]. Defaults to False.

        Returns:
            num_instances [int]:
        """
        if Print_Val:
            print(len(cls.instances))
        return len(cls.instances)

    @classmethod
    def _register_instance(cls, instance):
        """[summary]

        Args:
            instance ([type]): [description]
        """
        cls.instances.add(instance)
        for b in cls.__bases__:
            if issubclass(b, Registry):
                b._register_instance(instance)

    @classmethod
    def reset(cls) -> None:
        """Removes all instances from registry."""
        for magnet in cls._class_instances:
            del magnet
        cls.instances = WeakSet()
        cls._class_instances = []

    def __init_subclass__(cls):
        cls.instances = WeakSet()
        cls._class_instances = []


class Magnet(Registry):
    """Magnet base class


    Returns:
        Magnet: magnet base class
    """

    tol = 1e-4  # tolerance for rotations, sufficient for 0.01 degree accuracy
    mag_type = "Magnet"

    def __init__(self, *args, **kwargs) -> None:
        super().__init__()


def reset_magnets():
    """Returns a list of all instantiated magnets."""
    from ._magnet2 import Magnet_2D, Rectangle, Square, Circle
    from ._magnet3 import Magnet_3D, Prism, Cube, Cylinder, Sphere

    magnet_classes = [
        # Registry,
        Magnet,
        Magnet_2D,
        Rectangle,
        Square,
        Circle,
        Magnet_3D,
        Prism,
        Cube,
        Cylinder,
        Sphere,
    ]
    for cls in magnet_classes:
        cls.reset()


def list_magnets():
    """Returns a list of all instantiated magnets.

    Assumes that the child class registries have not been modified outside of
    using `pymagnet.reset_magnets()`.
    """
    return Magnet.print_instances()
