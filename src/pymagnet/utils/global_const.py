# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Global Constants
PI, PI/2, PI/4, µ0, and three internally used constants:
FP_CUTOFF = 1e-8, ALIGN_CUTOFF = 1e-5, MAG_TOL = 1e-4
"""

__all__ = ["PI", "PI_2", "PI_4", "MU0", "FP_CUTOFF", "ALIGN_CUTOFF", "MAG_TOL"]

# float: Module level PI.
from math import pi as PI

# float: Module level PI/2.
PI_2 = PI / 2.0

# float: Module level PI/4.
PI_4 = PI / 4.0

# float: Module level MU0, permittivity of free space.
MU0 = 4e-7 * PI

# float: Floating point cut off used for vector norm, and logs in mesh calculations
FP_CUTOFF = 1e-8

# Alignment precision for rotation of triangles in mesh magnets
ALIGN_CUTOFF = 1e-5

# Tolerance for computation of fields due to each magnetisation component Jx, Jy, Jz
# sufficient for 0.01 degree accuracy
MAG_TOL = 1e-4