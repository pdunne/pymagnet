__all__ = ["PI", "PI_2", "PI_4", "u0"]

"""float: Module level PI.
"""
from math import pi as PI


"""float: Module level PI/2.
"""
PI_2 = PI / 2.0

# float: Module level PI/4.
PI_4 = PI / 4.0

# float: Module level u0, permittivity of free space.
u0 = 4e-7 * PI


FP_CUTOFF = 1e-8

ALIGN_CUTOFF = 1e-5

MAG_TOL = 1e-4