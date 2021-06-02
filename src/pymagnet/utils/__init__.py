# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnets.utils

This module imports the classes and functions in the private modules to create a public API, including:

- Quaternion()
- Global Constants
- Point structures
- Vector structures

"""

from ._quaternion import *
from .global_const import *
from ._point_structs import *
from ._routines2D import *
from ._routines3D import *
from ._trigonometry3D import *

# Not implemented yet:
# from ._assemblies2D import *
# from ._fit import *