# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnets.magnets

This module imports the classes and functions in the private modules to create a public API. 

"""
from ._magnet_base import *
from ._magnet1D import *
from ._magnet2D import *
from ._magnet3D import *
from ._polygon2D import *
from ._polygon3D import *