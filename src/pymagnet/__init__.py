# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https: // mozilla.org / MPL / 2.0 / .
# Copyright 2021 Peter Dunne
"""pymagnet

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension

.. _Google Python Style Guide:
   https://google.github.io/styleguide/pyguide.html

"""
from . import plot
from . import magnets
from .magnets._magnet import reset_magnets, list_magnets
from .magnets._routines2 import grid2D, B_calc_2D
from .magnets._routines3 import grid3D, B_calc_3D
from math import pi as PI
"""float: Module level PI.
"""

PI_2 = PI / 2.0
"""float: Module level PI/2.
"""
PI_4 = PI / 4.0
"""float: Module level PI/4.
"""

u0 = 4E-7 * PI
"""float: Module level u0, permittivity of free space.
"""

__all__ = ['magnets', 'plot']
