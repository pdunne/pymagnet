from typing import Union

import numpy as _np
import numpy.typing as npt

Numeric64 = Union[float, npt.NDArray[_np.float64]]

Tuple2_64 = tuple[Numeric64, Numeric64]

Tuple3_64 = tuple[Numeric64, Numeric64, Numeric64]
