# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnets.magnets

This module imports the classes and functions in the private modules to create a public API.

TODO:
    * Quaternions and arbitrary rotation
    * Conversion between cylindrical and cartesian coordinates

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
from ._quaternion import *
