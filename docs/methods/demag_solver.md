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

## M(H) models

Two ways to provide the constitutive relation $M = f(H)$:

### `build_MH_interpolator(H_data, M_data)`

Constructs a monotonicity-preserving `PchipInterpolator` from measured M(H) data.

| Parameter | Type | Description |
|---|---|---|
| `H_data` | `np.ndarray` (1D) | Applied field values, strictly increasing |
| `M_data` | `np.ndarray` (1D) | Magnetization values, non-decreasing |

**Returns:** `PchipInterpolator` — callable mapping H -> M.

**Raises:** `ValueError` if arrays differ in shape, are not 1D, or violate monotonicity.

### `tanh_MH_model(Ms, chi)`

Returns an analytical M(H) function with built-in saturation:

$$
M = M_s \tanh\!\left(\frac{\chi \, H}{M_s}\right)
$$

At low fields this gives $M \approx \chi H$ (linear regime), and saturates smoothly to $M_s$ at high fields.

| Parameter | Type | Description |
|---|---|---|
| `Ms` | `float` | Saturation magnetization |
| `chi` | `float` | Dimensionless initial susceptibility (slope $dM/dH$ at $H=0$) |

**Returns:** Callable mapping H (float or array) -> M.

## Scipy solver

### `solve_demagnetization(H_ext, MH_interp, M_sat, N=0.5)`

General-purpose solver that accepts any callable M(H) model — either from `build_MH_interpolator` or `tanh_MH_model`.

| Parameter | Type | Description |
|---|---|---|
| `H_ext` | `float` | Applied external field magnitude |
| `MH_interp` | `Callable` | Any callable mapping H -> M |
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

## Numba-accelerated solver

When using the `tanh_MH_model`, the entire solve can run inside numba-compiled code, bypassing scipy. This gives ~80x speedup and supports parallel batch solving.

### `solve_demag_tanh(H_ext, Ms, chi, N=0.5)`

Solves for a single scalar H_ext value. Compiled with `@njit(cache=True)`.

| Parameter | Type | Description |
|---|---|---|
| `H_ext` | `float` | Applied external field magnitude |
| `Ms` | `float` | Saturation magnetization |
| `chi` | `float` | Dimensionless initial susceptibility |
| `N` | `float` | Demagnetizing factor (default 0.5) |

**Returns:** `(M_solution, H_int, converged)` tuple.

### `solve_demag_tanh_batch(H_ext_array, Ms, chi, N=0.5)`

Solves for an array of H_ext values in parallel. Compiled with `@njit(parallel=True, cache=True)` using `prange`.

| Parameter | Type | Description |
|---|---|---|
| `H_ext_array` | `np.ndarray` (1D) | Array of applied field values |
| `Ms` | `float` | Saturation magnetization |
| `chi` | `float` | Dimensionless initial susceptibility |
| `N` | `float` | Demagnetizing factor (default 0.5) |

**Returns:** `(M_solutions, H_int_solutions)` arrays with same shape as input.

## Usage

### Interpolated M(H) data

```python
import numpy as np
from demag_solver import build_MH_interpolator, solve_demagnetization

chi = 10.0
H_data = np.linspace(0, 1e6, 500)
M_data = chi * H_data

interp = build_MH_interpolator(H_data, M_data)

result = solve_demagnetization(
    H_ext=1e5, MH_interp=interp, M_sat=M_data[-1], N=0.5,
)
print(f"M = {result.M_solution:.2f}")   # M = 166666.67
print(f"H_int = {result.H_int:.2f}")    # H_int = 16666.67
```

### Analytical tanh model (scipy)

```python
from demag_solver import tanh_MH_model, solve_demagnetization

f_mh = tanh_MH_model(Ms=1e6, chi=10.0)

result = solve_demagnetization(
    H_ext=1e5, MH_interp=f_mh, M_sat=1e6, N=0.5,
)
print(f"M = {result.M_solution:.2f}")
```

### Numba-accelerated (single value)

```python
from demag_solver import solve_demag_tanh

M, H_int, converged = solve_demag_tanh(H_ext=1e5, Ms=1e6, chi=10.0, N=0.5)
```

### Numba-accelerated (batch, parallel)

```python
import numpy as np
from demag_solver import solve_demag_tanh_batch

H_ext = np.linspace(0, 5e5, 10_000)
M, H_int = solve_demag_tanh_batch(H_ext, Ms=1e6, chi=10.0, N=0.5)
```

## Solver strategy

### Scipy path (`solve_demagnetization`)

1. **Exact endpoints** — if `g(0) ≈ 0` or `g(M_sat) ≈ 0`, the solution is returned immediately without iteration.
2. **Brent's method** (`brentq`) — used when the bracket `[0, M_sat]` is valid (opposite signs at endpoints). Guaranteed convergence, tolerance `1e-6`.
3. **Fallback** (`fsolve`) — if the bracket is invalid (e.g. when extrapolating beyond the M(H) data), a Newton-type solver is used with initial guess `M_sat / 2`. A warning is issued.

### Numba path (`solve_demag_tanh` / `solve_demag_tanh_batch`)

Uses a fully compiled Brent's method implementation (`@njit`). Since the tanh residual is pure scalar math, the entire solve runs inside compiled code with no Python overhead. The batch solver parallelises across H_ext values using `prange`.

Both paths are stateless and side-effect free — safe for use in loops or with multiprocessing.
