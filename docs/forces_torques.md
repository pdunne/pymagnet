# Forces and Torques

The total force acting on a magnetic body with magnetisation $\mathbf{M}$ due to an external field $\mathbf{B'}$ can be divided into the sum of a volume and surface integral:

$$
\mathbf{F} = -\int_V (\nabla \cdot \mathbf{M}) \mathbf{B'} \: dv + \int_S (\mathbf{M} \cdot \mathbf{\hat{n}}) \mathbf{B'} \: ds
$$

which can be rewritten as

$$
\mathbf{F} = \int_V \rho_m \mathbf{B'} \: dv + \int_S \sigma_m \mathbf{B'} \: ds
$$

where $\rho_m = -\nabla \cdot \mathbf{M}$ and $\sigma_m =  \mathbf{M} \cdot \mathbf{\hat{n}}$ are the volume and surface charge densities.

For a uniformly magnetised body, $\rho_m = -\nabla \cdot \mathbf{M} = 0$, thus the total magnetic force acting on a magnet due to an external field $\mathbf{B'}$ can be deduced using the surface charge alone:

$$
\mathbf{F} =  \int_S \sigma_m \mathbf{B'} \: ds
$$

Similarly, the torque on a uniformly magnetised body can be determined using

$$
\mathbf{T} = \int_S \sigma_m (\mathbf{r} \times \mathbf{B'}) \: ds
$$

where $\mathbf{r}$ is the vector from the point about which the torque is computed, usually the centre of mass.

---

## API Reference

The `pymagnet.forces` module provides functions to calculate the force and torque on a magnet due to all other instantiated magnets. Each magnet type has a dedicated calculation function.

```python
from pymagnet.forces import calc_force_prism, calc_force_cylinder, calc_force_sphere
```

### Available Functions

| Function | Magnet Type | Description |
|----------|-------------|-------------|
| `calc_force_prism()` | Prism (cuboid) | Force/torque on a cuboidal magnet |
| `calc_force_cylinder()` | Cylinder | Force/torque on a cylindrical magnet |
| `calc_force_sphere()` | Sphere | Force/torque on a spherical magnet |
| `calc_force_mesh()` | Mesh (STL) | Force/torque on an STL mesh magnet |

---

### `calc_force_prism`

Calculates the total force and torque on a cuboidal (prism) magnet due to all other instantiated magnets.

```python
from pymagnet.forces import calc_force_prism
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `active_magnet` | Prism | required | The target magnet to calculate forces on |
| `num_samples` | int | `20` | Number of grid points per face edge (total points per face = num_samples²) |
| `unit` | str | `"mm"` | Length scale for the calculation |

**Returns:**

| Type | Description |
|------|-------------|
| `tuple[ndarray, ndarray]` | `(force, torque)` - Force vector (N) and torque vector (N·m) as 3-element arrays |

**Example:**

```python
import pymagnet as pm
from pymagnet.forces import calc_force_prism

pm.reset_magnets()

# Create two magnets - force will be calculated on magnet1
magnet1 = pm.magnets.Prism(Jr=1.0, center=(0, 0, 0), size=(10, 10, 5))
magnet2 = pm.magnets.Prism(Jr=1.0, center=(0, 0, 15), size=(10, 10, 5))

force, torque = calc_force_prism(magnet1, num_samples=30)
print(f"Force: {force} N")
print(f"Torque: {torque} N·m")
```

---

### `calc_force_cylinder`

Calculates the total force and torque on a cylindrical magnet due to all other instantiated magnets.

```python
from pymagnet.forces import calc_force_cylinder
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `active_magnet` | Cylinder | required | The target magnet to calculate forces on |
| `num_segments` | int | `20` | Number of segments for discretizing the circular faces |
| `unit` | str | `"mm"` | Length scale for the calculation |

**Returns:**

| Type | Description |
|------|-------------|
| `tuple[ndarray, ndarray]` | `(force, torque)` - Force vector (N) and torque vector (N·m) as 3-element arrays |

**Example:**

```python
import pymagnet as pm
from pymagnet.forces import calc_force_cylinder

pm.reset_magnets()

# Create two cylindrical magnets
magnet1 = pm.magnets.Cylinder(Jr=1.0, center=(0, 0, 0), radius=5, length=10)
magnet2 = pm.magnets.Cylinder(Jr=1.0, center=(0, 0, 20), radius=5, length=10)

force, torque = calc_force_cylinder(magnet1, num_segments=25)
print(f"Force: {force} N")
print(f"Torque: {torque} N·m")
```

---

### `calc_force_sphere`

Calculates the total force and torque on a spherical magnet due to all other instantiated magnets.

```python
from pymagnet.forces import calc_force_sphere
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `active_magnet` | Sphere | required | The target magnet to calculate forces on |
| `num_samples` | int | `200` | Number of sample points on the sphere surface |
| `unit` | str | `"mm"` | Length scale for the calculation |

**Returns:**

| Type | Description |
|------|-------------|
| `tuple[ndarray, ndarray]` | `(force, torque)` - Force vector (N) and torque vector (N·m) as 3-element arrays |

**Example:**

```python
import pymagnet as pm
from pymagnet.forces import calc_force_sphere

pm.reset_magnets()

# Create a sphere near a prism magnet
sphere = pm.magnets.Sphere(Jr=1.0, center=(0, 0, 0), radius=5)
prism = pm.magnets.Prism(Jr=1.0, center=(0, 0, 20), size=(10, 10, 5))

force, torque = calc_force_sphere(sphere, num_samples=300)
print(f"Force: {force} N")
print(f"Torque: {torque} N·m")
```

---

### `calc_force_mesh`

Calculates the total force and torque on an STL mesh magnet due to all other instantiated magnets. The mesh surface is subdivided into smaller triangles using midpoint subdivision for improved numerical integration accuracy.

```python
from pymagnet.forces._mesh_force import calc_force_mesh
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `active_magnet` | Mesh | required | The target mesh magnet to calculate forces on |
| `depth` | int | `3` | Subdivision depth (each triangle becomes 4^depth sub-triangles) |
| `unit` | str | `"mm"` | Length scale for the calculation |

**Returns:**

| Type | Description |
|------|-------------|
| `tuple[ndarray, ndarray]` | `(force, torque)` - Force vector (N) and torque vector (N·m) as 3-element arrays |

**Example:**

```python
import pymagnet as pm
from pymagnet.forces._mesh_force import calc_force_mesh

pm.reset()

magnet1 = pm.magnets.Mesh("cube1.stl", Jr=1.0, center=(0, 0, 0))
magnet2 = pm.magnets.Mesh("cube2.stl", Jr=1.0, center=(0, 0, 15))

force, torque = calc_force_mesh(magnet1, depth=3)
print(f"Force: {force} N")
print(f"Torque: {torque} N·m")
```

!!! note "Triangle Subdivision"
    The `depth` parameter controls integration accuracy. At `depth=3`, each mesh triangle is split into 64 sub-triangles. Higher values improve accuracy but increase computation time.

---

### Mesh Utility Functions

The `pymagnet.forces` module also exports utility functions for working with triangular meshes:

| Function | Description |
|----------|-------------|
| `triangle_area(triangle)` | Compute triangle area (njit-compiled) |
| `get_centroid(triangle)` | Get triangle centroid (njit-compiled) |
| `get_midpoints(triangle)` | Get edge midpoints |
| `divide_triangle_centroid(triangle, depth)` | Subdivide into 3^depth sub-triangles (centroid method) |
| `divide_triangle_regular(triangle, depth)` | Subdivide into 4^depth sub-triangles (midpoint method) |
| `get_area_triangles(triangles, area)` | Compute areas for an array of triangles |

---

## Notes

!!! note "Numerical Integration"
    Force and torque calculations use numerical surface integration. Increase the `num_samples` or `num_segments` parameter for higher accuracy at the cost of computation time.

!!! note "Multi-magnet Systems"
    The force functions calculate the effect of **all other instantiated magnets** on the target magnet. To calculate forces between specific pairs, reset magnets and instantiate only the magnets of interest.

!!! note "Units"
    - Force is returned in Newtons (N)
    - Torque is returned in Newton-meters (N·m)
    - The `unit` parameter specifies the length scale of the magnet coordinates (default: "mm")

---

## Complete Example

```python
import pymagnet as pm
from pymagnet.forces import calc_force_prism, calc_force_cylinder

pm.reset_magnets()

# Create a magnet array
base = pm.magnets.Prism(Jr=1.2, center=(0, 0, 0), size=(20, 20, 5))
cylinder = pm.magnets.Cylinder(Jr=1.0, center=(0, 0, 10), radius=3, length=6)

# Calculate force on the cylinder due to the base prism
force, torque = calc_force_cylinder(cylinder, num_segments=30)

print(f"Force on cylinder: Fx={force[0]:.4f}, Fy={force[1]:.4f}, Fz={force[2]:.4f} N")
print(f"Torque on cylinder: Tx={torque[0]:.6f}, Ty={torque[1]:.6f}, Tz={torque[2]:.6f} N·m")
```
