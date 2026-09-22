# What's New

## Unreleased

### Mesh Force and Torque Rewrite

`calc_force_mesh()` (and `Mesh.get_force_torque()`) no longer loop over mesh
triangles in Python. All sub-triangle centroids are now packed up front, each
source magnet's field is evaluated in a single batched call, and the
force/torque accumulation runs in a parallel Numba kernel. Results are
unchanged.

---

### Distance Cutoff for Mesh Fields: `r_cut`

`Mesh.get_field()` accepts an optional `r_cut` argument. Triangles whose
centroid is farther than `r_cut` from an evaluation point are skipped:

```python
Bx, By, Bz = magnet.get_field(X, Y, Z, r_cut=40.0)
```

The default, `np.inf`, performs no culling and is exact.

!!! danger "`r_cut` is an approximation, not a free speedup"
    Because surface charges largely cancel at a distance, discarding them costs
    far more accuracy than the raw triangle count suggests. On a 10 000-triangle
    bunny evaluated over ±70 mm, `r_cut=40 mm` is 8× faster but carries a 37%
    nRMSE, while an `r_cut` accurate to ~2% is only 1.1× faster. See
    [3D Magnets → Distance cutoff](magnets/magnets_3d.md#distance-cutoff-r_cut)
    for the full measured trade-off table, and validate against an
    `r_cut=np.inf` reference before relying on it.

---

### Multi-Magnet Fused Mesh Field

`pymagnet.magnets.get_total_field_mesh()` evaluates several `Mesh` magnets in a
single parallel pass over the evaluation points instead of one pass per magnet:

```python
from pymagnet import get_total_field_mesh

B = get_total_field_mesh([m1, m2], X, Y, Z)
```

It returns a `Field3` object and is numerically identical to summing the
individual `get_field()` calls. Note that on current benchmarks it is not yet
measurably faster than the sequential sum for large meshes.

---

### Parallel Mesh Kernel Race Condition Fixed

A race condition in the parallel mesh field kernel has been fixed, and the
kernel now parallelises over evaluation points rather than triangles.

---

### Faster Cylinder Field

`Cylinder` field evaluation now converts cylindrical to Cartesian components
directly from `x/rho` and `y/rho`, avoiding `arctan2`, `cos` and `sin` calls per
point. Behaviour on the symmetry axis (`rho = 0`) is unchanged.

---

## v0.5.1

### TOML Configuration System

Pymagnet now supports a declarative TOML-based configuration system for defining and running simulations without writing Python code. Define your magnets, grids, plots, and force calculations in a `.toml` file and run them from the command line or Python API.

```bash
pymagnet my_simulation.toml
```

```python
from pymagnet.config import run, load, validate

result = run("my_simulation.toml")
```

See [Configuration](configuration.md) for full documentation and examples.

---

### Numba Parallel Mesh Acceleration

Field calculations for STL mesh magnets are now parallelised using Numba's `prange`, achieving up to **247x speedup** on multi-core systems. This makes complex non-convex geometries practical for interactive use.

---

### `slice3D()` Utility

A new `slice3D()` function generates planar evaluation grids in 3D space, supporting `xy`, `xz`, and `yz` planes with configurable bounds and offset values. This simplifies creating cross-sectional field visualisations.

```python
import pymagnet as pm

points = pm.slice3D(plane="xz", max1=30, max2=30, slice_value=5.0, num_points=100)
```

---

### Mesh Force and Torque Calculations

Force and torque calculations now support STL mesh magnets via `calc_force_mesh()`. The mesh surface is subdivided into smaller triangles for improved numerical integration accuracy.

---

### NaN Handling in 3D Field Summation

Fixed a critical bug where NaN values from field singularities (e.g. inside magnets) could propagate and corrupt the entire field array in multi-magnet systems. NaN values are now zeroed before accumulation and the interior mask is reapplied after summation.

---

### Python 3.13+ and NumPy 2.0+

The minimum Python version has been bumped to **3.13**. NumPy 2.0.1+ and Numba 0.60-0.64 are now required.

---

### Optional Plotting Dependencies

Matplotlib and plotly are now optional dependencies. Install with plotting support using:

```bash
pip install pymagnet[plots]
```

The library will raise a clear error if plotting functions are called without the required packages installed.

---

### CI/CD with GitHub Actions

Automated testing and release workflows have been added via GitHub Actions, with Ruff linting and pre-commit hooks for code quality.
