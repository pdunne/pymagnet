# Pymagnet

User friendly magnetic field calculations in Python

## Calculating Magnetic Fields and Forces

Forms of this library have been used in a number of projects including [Liquid flow and control without solid walls, Nature 2020](https://www.nature.com/articles/s41586-020-2254-4).

## Features

This code uses analytical expressions to calculate the magnetic field due to
simple magnets. These include:

* 3D objects: cubes, cuboids, cylinders
* 2D: rectangles, squares

There are helper functions to plot the data as either line or countour plots,
but the underlying data is also accessible.

## Prerequisites

Ensure you have [Python](https://www.anaconda.com/) version >= 3.6
 (to use f-strings), and the following packages:

* numpy
* matplotlib

## Partial TODOs

* Implement arbitrary orientation in public code
* 3D rendering of slices

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
