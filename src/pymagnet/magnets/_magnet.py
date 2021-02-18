"""Magnet base class

This module demonstrates documentation as specified by the `Google
Python Style Guide`_. Docstrings may extend over multiple lines.
Sections are created with a section header and a colon followed by a
block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension
    * Type hinting

.. _Google Python Style Guide:   
http://google.github.io/styleguide/pyguide.html

"""
# __all__ = ['Magnet', 'reset_magnets', 'list_magnets']
__all__ = ['Magnet']

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
    """    
    instances = WeakSet()
    _hard_instances = []

    def __new__(cls, *args, **kwargs):
        o = object.__new__(cls)
        cls._register_instance(o)
        return o

    def __init__(self, *args, **kwargs):
        super().__init__()
        
    @classmethod
    def print_instances(cls):
        """Prints class instantiations
        """        
        if len(cls.instances) < 1:
            print("No Instances")
        else:    
            for instance in cls.instances:
                print(instance)

    @classmethod
    def get_instances(cls):
        return(cls.instances)
    
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
        cls._hard_instances.append(instance)
        for b in cls.__bases__:
            if issubclass(b, Registry):
                b._register_instance(instance)
                
    @classmethod
    def reset(cls) -> None:
        """Removes all magnet instances from registry.
        """
        for magnet in cls._hard_instances:
            del magnet
        cls.instances = WeakSet()
        cls._hard_instances = []

    def __init_subclass__(cls):
        cls.instances = WeakSet()
        cls._hard_instances = []


class Magnet(Registry):
    """Magnet base class


    Returns:
        Magnet: magnet base class
    """    
    counter = 0
    mag_type = 'Magnet'

    def __init__(self, *args, **kwargs) -> None:
        super().__init__()


def reset_magnets():
    """Returns a list of all instantiated magnets. 
    """
    from ._magnet2 import Magnet_2D, Rectangle, Square
    from ._magnet3 import Magnet_3D, Prism, Cube, Cylinder
    magnet_classes = [Registry, Magnet, 
                        Magnet_2D, Rectangle, Square,
                        Magnet_3D, Prism, Cube, Cylinder]
    for cls in magnet_classes:
        cls.reset()

def list_magnets():
    """Returns a list of all instantiated magnets. 

    Assumes that the child class registries have not been modified outside of 
    using `pymagnet.reset_magnets()`.
    """    
    return Magnet.print_instances()
