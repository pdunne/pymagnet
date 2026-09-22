# Configuration

Pymagnet supports a declarative TOML-based configuration system for defining simulations without writing Python code. A single `.toml` file specifies magnets, evaluation grids, plot settings, and force calculations.

## Quick Start

Create a file called `single_cube.toml`:

```toml
[meta]
name = "Single Cube"
description = "A single 10mm cube magnet with slice views"
dimension = "3D"

[[magnets]]
type = "Cube"
Jr = 1.0
width = 10.0
center = [0.0, 0.0, 0.0]
mask_magnet = true

[grid]
type = "slice"
max1 = 30.0
max2 = 30.0
num_points = 100
unit = "mm"

[plot]
enabled = true
type = "slice"
cmin = 0.0
cmax = 0.5
planes = ["xz", "xy"]
show_magnets = true
```

Run it from the command line:

```bash
pymagnet single_cube.toml
```

Or from Python:

```python
from pymagnet.config import run

result = run("single_cube.toml")
# result["magnets"] — list of magnet instances
# result["results"] — list of dicts with points, field, figure
# result["force"]   — force/torque dict (if enabled)
```

---

## TOML File Structure

A configuration file has four main sections:

### `[meta]` — Simulation Metadata

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `name` | str | `""` | Simulation name (used in output filenames) |
| `description` | str | `""` | Description for documentation |
| `dimension` | str | `"3D"` | `"2D"` or `"3D"` |

### `[[magnets]]` — Magnet Definitions

Each `[[magnets]]` entry defines one magnet. The `type` and `Jr` fields are required; all others have defaults.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `type` | str | required | Magnet type (see below) |
| `Jr` | float | required | Remnant magnetisation (T) |
| `center` | list | `[0, 0, 0]` | Centre position |
| `phi` | float | `90.0` | Magnetisation azimuthal angle (degrees) |
| `theta` | float | `0.0` | Magnetisation polar angle (degrees, 3D only) |
| `alpha` | float | `0.0` | Rotation about z-axis (degrees) |
| `beta` | float | `0.0` | Rotation about y-axis (degrees, 3D only) |
| `gamma` | float | `0.0` | Rotation about x-axis (degrees, 3D only) |
| `mask_magnet` | bool | `false` | Set field inside magnet to NaN |

**Dimension-specific parameters:**

| Type | Required Parameters |
|------|-------------------|
| `Rectangle` | `width`, `height` |
| `Square` | `width` |
| `Circle` | `radius` |
| `Prism` | `width`, `depth`, `height` |
| `Cube` | `width` |
| `Cylinder` | `radius`, `length` |
| `Sphere` | `radius` |
| `Mesh` | `filename` (path to STL file, relative to TOML file) |

Mesh magnets also accept `mesh_scale` (float, default `1.0`) for unit conversion.

### `[grid]` / `[[grid]]` — Evaluation Grid

Use `[grid]` for a single grid or `[[grid]]` for multiple grids (paired with `[[plot]]` entries).

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `type` | str | `"grid2D"` | Grid type: `"grid2D"`, `"slice"`, `"grid3D"`, `"volume"` |
| `num_points` | int | `100` | Points per axis |
| `unit` | str | `"mm"` | Length scale |

**For `grid2D` (2D simulations):**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `xmax` | float | `50.0` | Maximum x value |
| `ymax` | float | `50.0` | Maximum y value |
| `xmin` | float | `-xmax` | Minimum x value |
| `ymin` | float | `-ymax` | Minimum y value |

**For `slice` (3D planar cross-section):**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `plane` | str | `"xy"` | Plane: `"xy"`, `"xz"`, or `"yz"` |
| `max1` | float | `30.0` | Maximum along first axis |
| `max2` | float | `30.0` | Maximum along second axis |
| `min1` | float | `-max1` | Minimum along first axis |
| `min2` | float | `-max2` | Minimum along second axis |
| `slice_value` | float | `0.0` | Position on the third axis |

**For `grid3D` / `volume` (3D volume):**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `xmax` | float | `30.0` | Maximum x value |
| `ymax` | float | `30.0` | Maximum y value |
| `zmax` | float | `30.0` | Maximum z value |
| `xmin` | float | `-xmax` | Minimum x value |
| `ymin` | float | `-ymax` | Minimum y value |
| `zmin` | float | `-zmax` | Minimum z value |

### `[plot]` / `[[plot]]` — Visualisation Settings

Each `[[plot]]` entry is paired with the corresponding `[[grid]]` entry by index.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `enabled` | bool | `true` | Enable/disable plotting |
| `type` | str | `"contour"` | Plot type (see below) |
| `cmap` | str | `"viridis"` | Colourmap name |
| `cmin` | float | `None` | Minimum colour scale value |
| `cmax` | float | `None` | Maximum colour scale value |
| `num_levels` | int | `11` | Number of contour/isosurface levels |
| `num_arrows` | int | `None` | Number of vector arrows |
| `show_magnets` | bool | `true` | Render magnets on plot |
| `figure_group` | str | `""` | Group multiple plots onto one figure |
| `save_fig` | bool | `false` | Save figure to disk |
| `filename` | str | `""` | Output filename (auto-generated if empty) |

**2D plot types:** `"contour"`, `"streamplot"`

**3D plot types:** `"slice"`, `"volume"`

**3D-specific parameters:**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `opacity` | float | `0.8` | Field surface opacity |
| `magnet_opacity` | float | `1.0` | Magnet rendering opacity |
| `cone_opacity` | float | `1.0` | Vector arrow opacity |
| `planes` | list | `["xy", "xz", "yz"]` | Planes to render for slice plots |

### `[force]` — Force and Torque Calculation

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `enabled` | bool | `false` | Enable force calculation |
| `target_magnet` | int | `0` | Index of the magnet to calculate forces on (0-based) |
| `num_samples` | int | `20` | Number of integration samples per face |
| `unit` | str | `"mm"` | Length scale |

---

## Python API

### `run(toml_path, *, output_dir=None)`

Load config, create magnets, calculate fields, generate plots, and optionally compute forces.

```python
from pymagnet.config import run

result = run("my_config.toml", output_dir="output/")
```

**Returns** a dict with:

- `config` — the `SimulationConfig` dataclass
- `magnets` — list of magnet instances
- `results` — list of dicts (one per grid/plot pair), each containing `points`, `field`, and optionally `figure` and `saved_to`
- `force` — dict with `force` and `torque` arrays (if force enabled)

### `load(toml_path)`

Load and validate a TOML file without running the simulation.

```python
from pymagnet.config import load

config = load("my_config.toml")
print(config.name, len(config.magnets), "magnets")
```

### `validate(toml_path)`

Validate a TOML file and return a list of error messages. An empty list means the config is valid.

```python
from pymagnet.config import validate

errors = validate("my_config.toml")
if errors:
    for e in errors:
        print(f"  - {e}")
```

---

## CLI Usage

```bash
pymagnet CONFIG_FILE [-o OUTPUT_DIR]
```

| Argument | Description |
|----------|-------------|
| `CONFIG_FILE` | Path to TOML configuration file |
| `-o`, `--output` | Output directory for saved figures |

**Example:**

```bash
pymagnet examples/configs/single_cube.toml -o output/
```

---

## Examples

### 2D Contour Plot

```toml
[meta]
name = "2D Rectangles"
dimension = "2D"

[[magnets]]
type = "Rectangle"
Jr = 1.0
width = 20.0
height = 10.0
center = [-15.0, 0.0]

[[magnets]]
type = "Rectangle"
Jr = -1.0
width = 20.0
height = 10.0
center = [15.0, 0.0]

[grid]
type = "grid2D"
xmax = 40.0
ymax = 30.0
num_points = 100

[plot]
type = "contour"
cmin = 0.0
cmax = 0.5
num_arrows = 11
save_fig = true
filename = "2d_rectangles.png"
```

### Force Calculation Between Two Cubes

```toml
[meta]
name = "Two Cubes - Force"
dimension = "3D"

[[magnets]]
type = "Cube"
Jr = 1.0
width = 10.0
center = [0.0, 0.0, 0.0]

[[magnets]]
type = "Cube"
Jr = 1.0
width = 10.0
center = [0.0, 0.0, 15.0]

[grid]
type = "slice"
plane = "xz"
max1 = 30.0
max2 = 30.0
num_points = 100

[plot]
type = "slice"
cmin = 0.0
cmax = 0.8
planes = ["xz"]

[force]
enabled = true
target_magnet = 1
num_samples = 20
```

### Grouped Slice Plots (Multiple Planes on One Figure)

Use `figure_group` to render multiple grid/plot pairs onto a single 3D figure. This is useful for showing different cross-sections of an assembly simultaneously.

```toml
[meta]
name = "Multi-Slice Assembly"
dimension = "3D"

[[magnets]]
type = "Prism"
Jr = 1.26
width = 30.0
depth = 60.0
height = 30.0
center = [0.0, 0.0, 30.0]
mask_magnet = true

# XY slice at z=15
[[grid]]
type = "slice"
max1 = 30.0
max2 = 60.0
min2 = -60.0
slice_value = 15.0
num_points = 100

[[plot]]
type = "slice"
figure_group = "main"
cmin = 0.0
cmax = 1.0
planes = ["xy"]
num_arrows = 12

# XZ slice at y=0
[[grid]]
type = "slice"
max1 = 30.0
max2 = 60.0
min1 = -30.0
min2 = -60.0
num_points = 100

[[plot]]
type = "slice"
figure_group = "main"
cmin = 0.0
cmax = 1.0
planes = ["xz"]
```

Plots sharing the same `figure_group` string are combined into one plotly figure.

---

## More Examples

See the [examples/configs/](https://github.com/pdunne/pymagnet/tree/main/examples/configs) directory for 14 ready-to-use configuration files covering single magnets, Halbach arrays, STL meshes, force calculations, and grouped plots.
