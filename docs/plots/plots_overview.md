# Plots Overview

The `pymagnet.plots` module provides visualization functions for magnetic field data calculated
by the library. Two plotting backends are supported:

- **Matplotlib** - For 1D line plots and 2D contour/streamplots
- **Plotly** - For interactive 3D slice and volume plots

## Available Functions

### Matplotlib-based (2D)

| Function | Description |
|----------|-------------|
| `plot_1D_field()` | Line plot along central axis of cylinder/prism |
| `plot_2D_line()` | Line plot of field components for 2D magnets |
| `plot_2D_contour()` | Contour or streamplot for 2D magnetic sources |
| `plot_3D_contour()` | 2D slice of a 3D field calculation |
| `plot_sub_contour_3D()` | Single field component contour plot |

### Plotly-based (3D)

| Function | Description |
|----------|-------------|
| `plot_magnet()` | Render magnet geometry in 3D |
| `slice_plot()` | Plot pre-computed field slices |
| `slice_quickplot()` | Compute and plot field slices |
| `volume_plot()` | Plot pre-computed volume data |
| `volume_quickplot()` | Compute and plot volume data |

## Quick Example

```python
import pymagnet as pm

# Create a magnet
magnet = pm.magnets.Cylinder(radius=5, length=10, Jr=1.0)

# 1D plot along symmetry axis
pm.plots.plot_1D_field(magnet, unit="mm")

# 3D interactive slice plot
fig, cache, data = pm.plots.slice_quickplot(
    planes=["xy", "xz"],
    max1=20,
    max2=20,
    cmax=0.5
)
```

## Dependencies

- **matplotlib** - Required for all 2D plotting functions
- **plotly** - Required for all 3D plotting functions

Both are optional dependencies. If a backend is not installed, the corresponding
functions will raise an `ImportError` with a descriptive message.

## See Also

- [2D Plotting](plots_2d.md) - Detailed documentation for matplotlib-based plots
- [3D Plotting](plots_3d.md) - Detailed documentation for plotly-based plots
