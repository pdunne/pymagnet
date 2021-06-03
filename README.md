# Pymagnet

User friendly magnetic field calculations in Python

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-blue.svg)](https://opensource.org/licenses/MPL-2.0)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)
[![DOI](https://zenodo.org/badge/339667292.svg)](https://zenodo.org/badge/latestdoi/339667292)
[![Anaconda-Server Badge](https://anaconda.org/pdunne/pymagnet/badges/version.svg)](https://anaconda.org/pdunne/pymagnet)

## Getting Started

Installing `pymagnet` can be done using

```bash
python -m pip install pymagnet 
```

or

```bash
conda install -c pdunne pymagnet
```

Pymagnet is a collection of routines to calculate and plot the magnetic field due to arbitrary 2D
and 3D objects, like cubes or cylinders, as well as complex non-convex structures stored in STL
files.

The approach assumes the magnets are uniformly magnetised, and fully transparent to magnetic fields.
There are some drawbacks to this compared to Finite Element Methods (FEM), but with the advantage of
significantly faster calculations.

The current version is written in Python with some speed up using [Numpy](https://numpy.org/) and
[Numba](https://numba.pydata.org/), but the backend is being ported to
[Rust](https://github.com/pdunne/magnet_rs) for improved performance.

## Documentation

Full documentation can be found here: [https://pdunne.github.io/pymagnet/](https://pdunne.github.io/pymagnet/)

### Examples

Examples can be found [in the repository](https://github.com/pdunne/pymagnet/tree/main/examples).

## Features

This code uses analytical expressions to calculate the magnetic field due to various magnets,
and forces and torques on a magnet due to all other magnets.

These include simple objects such as:

* 3D: cubes, prisms (cuboids), cylinders, spheres

<figure>
  <img src="docs/img/3d_example_slice_1.png" width=400/>
  <img src="docs/img/3d_example_volume_2.png" width=400/>
  <figcaption>Cylinder Plots</figcaption>
</figure>

* 2D: rectangles, squares, circles

<figure>
  <img src="docs/img/2d_circle_contour.png" width=400/>
  <img src="docs/img/2d_circle_stream.png" width=400/>
  <figcaption>2D contour plot and streamplot of a long bipolar rod</figcaption>
</figure>

and complex compound objects:

* 3D: Polyhedra stored as STL meshes
* 2D: Polygons constructed from collections of line elements

There are helper functions to plot the data as line, contour, slice, and volume plots,
but the underlying data is also accessible.

## Prerequisites

Ensure you have [Python](https://www.anaconda.com/) version >= 3.6
 (to use f-strings), and the following packages:

* numpy
* numpy-stl
* numba
* matplotlib
* plotly

## Usage

Forms of this library have been used in a number of projects including [Liquid flow and control without solid walls, Nature 2020](https://www.nature.com/articles/s41586-020-2254-4).

## Licensing

Source code licensed under the [Mozilla Public License Version 2.0](https://www.mozilla.org/en-US/MPL/2.0/)

Documentation is licensed under a Creative Commons Attribution-ShareAlike 4.0 International [(CC BY-SA 4.0)](https://creativecommons.org/licenses/by-sa/4.0/) license.

This is a human-readable summary of (and not a substitute for) the license, adapted from [CS50x](https://cs50.harvard.edu/x/2021/license/). Official translations of this license are available in other languages.

**You are free to:**

* **Share** — copy and redistribute the material in any medium or format.
* **Adapt** — remix, transform, and build upon the material.

**Under the following terms:**

* **Attribution** — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
* **ShareAlike** — If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original
* No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you shall be licensed as above, without any
additional terms or conditions.
