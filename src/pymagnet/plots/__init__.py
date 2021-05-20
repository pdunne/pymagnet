# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""pymagnets.plots

This module imports the classes and functions in the private modules to create a public API.

- Lines and contour plots are drawn using matplotlib
- 3D surface and volume plots are rendered using plotly.

"""
from ._plot1D import *
from ._plot2D import *
from ._plotly3D import *

# from ._plotly2D import *
