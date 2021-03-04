# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnets.magnets

This module imports the classes and functions in the private modules to create a public API.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

TODO:
    * Quaternions and arbitrary rotation
    * Conversion between cylindrical and cartesian coordinates


.. _Google Python Style Guide:
   https://google.github.io/styleguide/pyguide.html
"""
from ._magnet import *
from ._magnet1 import *
from ._magnet2 import *
from ._magnet3 import *
from ._fields import *
from ._routines import *
from ._routines2 import *
from ._routines3 import *
from ._assembly2 import *
