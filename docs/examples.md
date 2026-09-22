# Examples

Examples can be found [in the repository](https://github.com/pdunne/pymagnet/tree/main/examples).

**Or run on Colab:**

## 2D Examples

Getting Started [![First Steps](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/First%20Steps.ipynb)

1D Simple Plots [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Fields/1D%20Examples.ipynb)

2D Standard Shapes [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Fields/2D%20Examples.ipynb)

2D Arbitrary Polygons [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Fields/2D%20PolyMagnet.ipynb)

Misc. Examples [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Fields/Misc%20Examples.ipynb)

## 3D Examples

Analytical Shapes [![3D Examples](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Fields/3D%20Examples%20-%20Assemblies.ipynb)

Spheres [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Fields/3D_Examples%20Spheres.ipynb)

3D STL Meshes [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/STL%20Magnets/STL%20Examples.ipynb)

3D Stanford Bunny [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/STL%20Magnets/Stanford%20Bunny.ipynb)

Upload Your Own Custom Magnet [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/STL%20Magnets/Colab%20Custom%20Magnet.ipynb)

## Forces and Torques

Cubes [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Forces%20Torques/Cubes.ipynb)

Cylinders [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Forces%20Torques/Cylinders.ipynb)

Spheres [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/Forces%20Torques/Spheres.ipynb)

STL Cubes [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/STL%20Magnets/STL%20Examples.ipynb)

STL Pentagonal Prisms [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pdunne/pymagnet/blob/main/examples/notebooks/STL%20Magnets/STL%20Forces%20Pentagon.ipynb)

## Configuration Examples

Pymagnet includes ready-to-use TOML configuration files in [`examples/configs/`](https://github.com/pdunne/pymagnet/tree/main/examples/configs). Run any of them with:

```bash
pymagnet examples/configs/single_cube.toml
```

| Config File | Description |
|-------------|-------------|
| `single_cube.toml` | Single cube magnet with slice views |
| `single_cylinder.toml` | Single cylinder magnet |
| `2d_rectangles.toml` | Two 2D rectangular magnets with contour plot |
| `cube_halbach.toml` | 8-prism Halbach array |
| `pseudo_Halbach.toml` | Pseudo-Halbach field assembly |
| `pseudo_Halbach_plot_group.toml` | Grouped multi-plane slice plots |
| `quadrupolar.toml` | Quadrupolar magnet arrangement |
| `two_poles.toml` | Two-pole configuration |
| `cube_volume.toml` | Volume rendering of a cube magnet |
| `two_cubes_force.toml` | Force calculation between two cubes |
| `mesh_cube.toml` | STL mesh cube magnet |
| `mesh_star.toml` | STL mesh star shape |
| `mesh_bunny.toml` | Stanford bunny STL mesh |
| `mesh_two_cubes_force.toml` | Force calculation with STL mesh magnets |

See the [Configuration](configuration.md) page for full documentation.

## Binder

The example notebooks can be run as an instance using Binder:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pdunne/pymagnet/main?filepath=examples%2Fnotebooks)
