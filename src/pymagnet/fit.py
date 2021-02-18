"""Fit data for magnetic field models

This module demonstrates documentation as specified by the `Google
Python Style Guide`_. Docstrings may extend over multiple lines.
Sections are created with a section header and a colon followed by a
block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * Reimplement all 
    * You have to also use ``sphinx.ext.todo`` extension
    * 

.. _Google Python Style Guide:   
http://google.github.io/styleguide/pyguide.html

"""
pass

# def fit_cube(z, Jr, dz):
#     xDim = 10/2000
#     yDim = 10/2000
#     zDim = 10/2000
#     return(mag.Bzz1d(xDim, yDim, zDim, Jr, z-dz))
#     # return(mag.solBz(xDim,yDim,yDim + y-dy,n,I))


# def fit_cyl(y, Jr, dy):
#     R = 14/2000
#     HH = 10/2000
#     n = 1
#     I = Jr/(u0*n)
#     # return(mag.Byy2d(xDim, yDim, Jr, 0, yDim + y-dy))
#     return(mag.solBz(R, HH, HH+y-dy, n, I))
