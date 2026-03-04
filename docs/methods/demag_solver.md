# Demagnetization Solver

Self-consistent solver for calculating the magnetization **M** of a soft magnetic element, accounting for the demagnetizing field. Given an applied field **H_ext** and a measured M(H) curve, the solver finds the operating point where the material response and the demagnetizing load line intersect.

## Physics

The internal field of a magnetized body is reduced from the applied field by its own demagnetizing contribution:

$$
H_\text{int} = H_\text{ext} - N \, M
$$

where **N** is the demagnetizing factor (geometry-dependent). The magnetization is determined self-consistently by the material's constitutive relation:

$$
M = f_{MH}(H_\text{int})
$$

The solver reformulates this as a 1D root-finding problem:

$$
g(M) = M - f_{MH}\!\left(H_\text{ext} - N \, M\right) = 0
$$

and solves over the bracket $[0,\; M_\text{sat}]$.

### Demagnetizing factors for common geometries

| Geometry | N |
|---|---|
| 2D circle (infinite cylinder) | 1/2 |
| Sphere | 1/3 |
| Thin film (normal to surface) | 1 |

## API

### `build_MH_interpolator(H_data, M_data)`

Constructs a monotonicity-preserving `PchipInterpolator` from measured M(H) data.

| Parameter | Type | Description |
|---|---|---|
| `H_data` | `np.ndarray` (1D) | Applied field values, strictly increasing |
| `M_data` | `np.ndarray` (1D) | Magnetization values, non-decreasing |

**Returns:** `PchipInterpolator` — callable mapping H → M.

**Raises:** `ValueError` if arrays differ in shape, are not 1D, or violate monotonicity.

### `solve_demagnetization(H_ext, MH_interp, M_sat, N=0.5)`

Solves the self-consistent demagnetization equation for a single H_ext value.

| Parameter | Type | Description |
|---|---|---|
| `H_ext` | `float` | Applied external field magnitude |
| `MH_interp` | `PchipInterpolator` | Prebuilt interpolator from `build_MH_interpolator` |
| `M_sat` | `float` | Saturation magnetization (upper bound for root search) |
| `N` | `float` | Demagnetizing factor (default 0.5) |

**Returns:** `DemagResult` (see below).

### `DemagResult`

Dataclass returned by `solve_demagnetization`.

| Field | Type | Description |
|---|---|---|
| `M_solution` | `float` | Self-consistent magnetization |
| `H_int` | `float` | Internal field: $H_\text{ext} - N \cdot M$ |
| `converged` | `bool` | Whether the solver converged |
| `solver_used` | `str` | `"brentq"`, `"fsolve"`, or `"exact"` |

## Usage

```python
import numpy as np
from demag_solver import build_MH_interpolator, solve_demagnetization

# Define a linear M(H) curve: M = chi * H
chi = 10.0
H_data = np.linspace(0, 1e6, 500)
M_data = chi * H_data

# Build interpolator (once)
interp = build_MH_interpolator(H_data, M_data)

# Solve for a given applied field
result = solve_demagnetization(
    H_ext=1e5,
    MH_interp=interp,
    M_sat=M_data[-1],
    N=0.5,  # 2D circle
)

print(f"M = {result.M_solution:.2f}")   # M = 166666.67
print(f"H_int = {result.H_int:.2f}")    # H_int = 16666.67

# Analytical check: M = chi * H_ext / (1 + N * chi) = 166666.67
```

The interpolator can be reused across many H_ext values without rebuilding:

```python
for H in [1e4, 5e4, 1e5, 5e5]:
    r = solve_demagnetization(H, interp, M_data[-1], N=0.5)
    print(f"H_ext={H:.0e}  →  M={r.M_solution:.2f}, H_int={r.H_int:.2f}")
```

## Solver strategy

1. **Exact endpoints** — if `g(0) ≈ 0` or `g(M_sat) ≈ 0`, the solution is returned immediately without iteration.
2. **Brent's method** (`brentq`) — used when the bracket `[0, M_sat]` is valid (opposite signs at endpoints). Guaranteed convergence, tolerance `1e-6`.
3. **Fallback** (`fsolve`) — if the bracket is invalid (e.g. when extrapolating beyond the M(H) data), a Newton-type solver is used with initial guess `M_sat / 2`. A warning is issued.

The solver is stateless and side-effect free — safe for use in loops or with multiprocessing.
