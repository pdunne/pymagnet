# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""3D Plotting routines

This module contains all functions needed to plot 3D contours for 3D magnetic sources.
Unlike the plot2D module, here plotly is used as the backend.

TODO:
    * Update __str__ and __repr__ for polyhedra
"""
from ..magnets._magnet_base import Registry
import numpy as _np
from ..magnets import *
from ..utils import grid2D, grid3D, slice3D, B_calc_3D
from ..utils.global_const import PI, MAG_TOL
import plotly.graph_objects as _go


class Polyhedron(Registry):
    """Encodes magnet dimensions for drawing a polyhedon on 3D plots

    Polyhedra:
        Cuboid
        Cylinder
        Sphere
        Mesh
    """

    # Tolerance for minimum angle needed for rotation of object
    tol = MAG_TOL

    def __init__(self, center, size, **kwargs):
        """Initialises a cuboid

        Args:
            center (tuple): x,y,z
            size (tuple): x,y,z
        """
        super().__init__()

        self.center = _np.asarray(center)
        self.size = _np.asarray(size)

        self.alpha = kwargs.pop("alpha", 0.0)
        self.alpha_rad = _np.deg2rad(self.alpha)
        self.beta = kwargs.pop("beta", 0.0)
        self.beta_rad = _np.deg2rad(self.beta)
        self.gamma = kwargs.pop("gamma", 0.0)
        self.gamma_rad = _np.deg2rad(self.gamma)
        self.color = kwargs.pop("color", "white")

    def __repr__(self) -> str:
        return f"(center: {self.center}, size: {self.size} )"

    def __str__(self) -> str:
        return f"(center: {self.center}, size: {self.size} )"

    #     # FIXME: Move below function into quaternion class, do the same for Magnet_3D class
    #     def _generate_rotation_quaternions(self):
    #         """Generates single rotation quaternion for all non-zero rotation angles,
    #         which are:
    #
    #             alpha (float): angle in degrees around z-axis
    #             beta (float): angle in degrees around y-axis
    #             gamma (float): angle in degrees around x-axis
    #
    #         Returns:
    #             Quaternion: total rotation quaternion
    #         """
    #         from ..utils._quaternion import Quaternion
    #
    #         rotate_about_x = Quaternion()
    #         rotate_about_y = Quaternion()
    #         rotate_about_z = Quaternion()
    #
    #         forward_rotation, reverse_rotation = Quaternion(), Quaternion()
    #
    #         if _np.fabs(self.alpha_rad) > MAG_TOL:
    #             rotate_about_z = Quaternion.q_angle_from_axis(self.alpha_rad, (0, 0, 1))
    #
    #         if _np.fabs(self.beta_rad) > MAG_TOL:
    #             rotate_about_y = Quaternion.q_angle_from_axis(self.beta_rad, (0, 1, 0))
    #
    #         if _np.fabs(self.gamma_rad) > MAG_TOL:
    #             rotate_about_x = Quaternion.q_angle_from_axis(self.gamma_rad, (1, 0, 0))
    #
    #         forward_rotation = rotate_about_x * rotate_about_z * rotate_about_y
    #         reverse_rotation = forward_rotation.get_conjugate()
    #
    #         return forward_rotation, reverse_rotation

    def generate_vertices(self):
        """Generates vertices of a polyhedron

        This should be implemented for each Polyhedron child class
        """
        pass


class Graphic_Cuboid(Polyhedron):
    """Generates

    Args:
        center (tuple, optional): Cuboid center. Defaults to (0, 0, 0).
        size (tuple, optional): Size of cuboid. Defaults to (1, 1, 1).

    Kwargs:
        alpha (float):
        beta (float):
        gamma (float):
        color (float):
    """

    def __init__(self, center=(0, 0, 0), size=(1, 1, 1), **kwargs):
        """Init method

        Args:
            center (tuple, optional): [description]. Defaults to (0, 0, 0).
            size (tuple, optional): [description]. Defaults to (1, 1, 1).
        """
        super().__init__(center, size, **kwargs)

        self.vertices = self.generate_vertices()

    def generate_vertices(self):
        """Generates and rotates vertices of a cuboid based on orientation angles

        Returns:
            ndarray: 3xN array of vertex coordinates (columns are x, y, z)
        """
        # Generate and rotate the vertices
        if _np.any(
            _np.fabs(
                _np.array(
                    [
                        self.alpha_rad,
                        self.beta_rad,
                        self.gamma_rad,
                    ]
                )
            )
            > Polyhedron.tol
        ):

            # _, reverse_rotation = self._generate_rotation_quaternions()
            forward_rotation = Quaternion.gen_rotation_quaternion(
                self.alpha_rad, self.beta_rad, self.gamma_rad
            )
            reverse_rotation = forward_rotation.get_conjugate()
            # Generate 3xN array for quaternion rotation

            vertex_coords = self._gen_vertices(center=(0, 0, 0), size=self.size)

            # Rotate points
            x, y, z = reverse_rotation * vertex_coords

            # Reconstruct 3xN array and add center offset
            vertex_coords = _np.vstack([x, y, z])
            vertex_coords += _np.array(self.center).reshape(-1, 1)

            # finally return the coordinates
            return vertex_coords

        # Only generate
        else:
            vertex_coords = self._gen_vertices(self.center, self.size)
            return vertex_coords

    @staticmethod
    def _gen_vertices(center=(0, 0, 0), size=(1, 1, 1)):
        """Generates coordinates for all cuboid vertices

        Args:
            center (tuple, optional): x,y,z coordinates. Defaults to (0, 0, 0)
            size (tuple, optional): scale the cuboid in x, y, z directons. Defaults to (1, 1, 1).

        Returns:
            ndarray: numpy array of shape (3, 8)
        """
        # Center of this cube is (0.5, 0.5, 0.5)
        x = _np.array([0, 0, 1, 1, 0, 0, 1, 1]) - 0.5
        y = _np.array([0, 1, 1, 0, 0, 1, 1, 0]) - 0.5
        z = _np.array([0, 0, 0, 0, 1, 1, 1, 1]) - 0.5

        vertex_coords = _np.vstack([x, y, z])

        # Scale x, y, z by the size array
        vertex_coords = _np.multiply(vertex_coords.T, size).T

        # Offset by externally set center
        vertex_coords += _np.array(center).reshape(-1, 1)

        return vertex_coords


class Graphic_Sphere(Polyhedron):
    """Generates vertices for a sphere

    Args:
        center (tuple, optional): [description]. Defaults to (0, 0, 0).
        radius (float, optional): [description]. Defaults to 1.

    Kwargs:
        alpha (float):
        beta (float):
        gamma (float):
        color (float):
    """

    def __init__(self, center=(0, 0, 0), radius=1, **kwargs):
        """Init method

        Args:
            center (tuple, optional): [description]. Defaults to (0, 0, 0).
            size (float, optional): [description]. Defaults to (1).
        """
        super().__init__(center, size=(radius, radius, radius), **kwargs)
        self.radius = radius
        self.vertices = self.generate_vertices()

    def generate_vertices(self):
        """Generates and rotates vertices of a cuboid based on orientation angles

        Returns:
            ndarray: 3xN array of vertex coordinates (columns are x, y, z)
        """
        # Generate and rotate the vertices
        if _np.any(
            _np.fabs(
                _np.array(
                    [
                        self.alpha_rad,
                        self.beta_rad,
                        self.gamma_rad,
                    ]
                )
            )
            > Polyhedron.tol
        ):

            # _, reverse_rotation = self._generate_rotation_quaternions()
            forward_rotation = Quaternion.gen_rotation_quaternion(
                self.alpha_rad, self.beta_rad, self.gamma_rad
            )
            reverse_rotation = forward_rotation.get_conjugate()

            # Generate 3xN array for quaternion rotation

            vertex_coords = self._gen_vertices(center=(0, 0, 0), radius=self.radius)

            # Rotate points
            x, y, z = reverse_rotation * vertex_coords

            # Reconstruct 3xN array and add center offset
            vertex_coords = _np.vstack([x, y, z])
            vertex_coords += _np.array(self.center).reshape(-1, 1)

            # finally return the coordinates
            return vertex_coords

        # Only generate
        else:
            vertex_coords = self._gen_vertices(self.center, self.radius)
            return vertex_coords

    @staticmethod
    def _gen_vertices(center=(0, 0, 0), radius=1):
        """Generates coordinates for approximate sphere vertices

        Args:
            center (tuple, optional): x,y,z coordinates. Defaults to (0, 0, 0)
            radius (float, optional): radius of the sphere. Defaults to 1.

        Returns:
            ndarray: numpy array of shape (3, 8)
        """
        u, v = _np.mgrid[0 : 2 * PI : 20j, 0:PI:10j]
        x = radius * _np.cos(u) * _np.sin(v)
        y = radius * _np.sin(u) * _np.sin(v)
        z = radius * _np.cos(v)

        vertex_coords = _np.vstack([x.ravel(), y.ravel(), z.ravel()])

        # Offset by externally set center
        vertex_coords += _np.array(center).reshape(-1, 1)

        return vertex_coords


class Graphic_Cylinder(Polyhedron):
    """Generates vertices for a Cylinder

    Args:
        center (tuple, optional): [description]. Defaults to (0, 0, 0).
        radius (float, optional): [description]. Defaults to 1.
        length (float, optional): [description]. Defaults to 1.

    Kwargs:
        alpha (float):
        beta (float):
        gamma (float):
        color (float):
    """

    def __init__(self, center=(0, 0, 0), radius=1, length=1, **kwargs):
        """Init method

        Args:
            center (tuple, optional): [description]. Defaults to (0, 0, 0).
            size (float, optional): [description]. Defaults to (1).
        """
        super().__init__(center, size=(radius, radius, length), **kwargs)
        self.radius = radius
        self.length = length
        self.vertices = self.generate_vertices()

    def generate_vertices(self):
        """Generates and rotates vertices of a cuboid based on orientation angles

        Returns:
            ndarray: 3xN array of vertex coordinates (columns are x, y, z)
        """
        # Generate and rotate the vertices
        if _np.any(
            _np.fabs(
                _np.array(
                    [
                        self.alpha_rad,
                        self.beta_rad,
                        self.gamma_rad,
                    ]
                )
            )
            > Polyhedron.tol
        ):

            # _, reverse_rotation = self._generate_rotation_quaternions()
            forward_rotation = Quaternion.gen_rotation_quaternion(
                self.alpha_rad, self.beta_rad, self.gamma_rad
            )
            reverse_rotation = forward_rotation.get_conjugate()

            # Generate 3xN array for quaternion rotation

            vertex_coords = self._gen_vertices(
                center=(0, 0, 0), radius=self.radius, length=self.length
            )

            # Rotate points
            x, y, z = reverse_rotation * vertex_coords

            # Reconstruct 3xN array and add center offset
            vertex_coords = _np.vstack([x, y, z])
            vertex_coords += _np.array(self.center).reshape(-1, 1)

            # finally return the coordinates
            return vertex_coords

        # Only generate
        else:
            vertex_coords = self._gen_vertices(self.center, self.radius, self.length)
            return vertex_coords

    @staticmethod
    def _gen_vertices(center=(0, 0, 0), radius=1, length=1):
        """Generates coordinates for approximate cylinder vertices

        Args:
            center (tuple, optional): x,y,z coordinates. Defaults to (0, 0, 0)
            radius (float, optional): radius of the cylinder. Defaults to 1.
            length (float, optional): length of the cylinder. Defaults to 1.


        Returns:
            ndarray: numpy array of shape (3, 8)
        """
        rho, z = _np.mgrid[0 : 2 * PI : 20j, -length / 2 : length / 2 : 2j]
        x = radius * _np.cos(rho)
        y = radius * _np.sin(rho)

        vertex_coords = _np.vstack([x.ravel(), y.ravel(), z.ravel()])

        # Offset by externally set center
        vertex_coords += _np.array(center).reshape(-1, 1)
        return vertex_coords


class Graphic_Mesh(Polyhedron):
    """Generates

    Args:
        center (tuple, optional): Cuboid center. Defaults to (0, 0, 0).
        size (tuple, optional): Size of cuboid. Defaults to (1, 1, 1).

    Kwargs:
        color (float):
    """

    def __init__(self, mesh_vectors, start, stop, **kwargs):
        """Init method

        Args:
            center (tuple, optional): [description]. Defaults to (0, 0, 0).
            size (tuple, optional): [description]. Defaults to (1, 1, 1).
        """
        self.color = kwargs.pop("color", "white")
        self.mesh_vectors = mesh_vectors[start:stop]

    def generate_vertices(self):
        # return super().generate_vertices()
        p, q, r = self.mesh_vectors.shape  # (p, 3, 3)

        # the array stl_mesh.vectors.reshape(p*q, r) can contain multiple copies of the same vertex;
        # extract unique vertices from all mesh triangles
        vertices, ixr = _np.unique(
            self.mesh_vectors.reshape(p * q, r), return_inverse=True, axis=0
        )
        I = _np.take(ixr, [3 * k for k in range(p)])
        J = _np.take(ixr, [3 * k + 1 for k in range(p)])
        K = _np.take(ixr, [3 * k + 2 for k in range(p)])
        x, y, z = vertices.T
        trace = _go.Mesh3d(x=x, y=y, z=z, i=I, j=J, k=K, color=self.color)

        # optional parameters to make it look nicer
        trace.update(
            flatshading=True, lighting_facenormalsepsilon=0, lighting_ambient=0.7
        )
        return trace


def reset_polyhedra():
    """Returns a list of all instantiated polyhedra."""

    polyhedra_classes = [
        Polyhedron,
        Graphic_Cuboid,
        Graphic_Sphere,
        Graphic_Cylinder,
    ]
    for cls in polyhedra_classes:
        cls.reset()


def list_polyhedra():
    """Returns a list of all instantiated polyhedra.

    Assumes that the child class registries have not been modified outside of
    using `pymagnet.reset_magnets()`.
    """
    return Polyhedron.print_instances()


def _draw_surface_slice(
    points,
    field,
    colorscale="viridis",
    opacity=1.0,
    cmin=0.0,
    cmax=1.0,
    showscale=True,
):
    """Generates plotly surface object for plotting.

    Args:
        x (ndarray): x-coordinates
        y (ndarray): y-coordinates
        z (ndarray): z-coordinates
        B (Vector3): Magnetic field vector
        colorscale (str, optional): colormap. Defaults to "viridis".
        opacity (float, optional): surface opacity. Defaults to 1.0.
        cmin (float, optional): color scale lower limit. Defaults to 0.0.
        cmax (float, optional): color scale upper limit. Defaults to 1.0.
        showscale (bool, optional): show colorbar. Defaults to True.

    Returns:
        Surface object: plotly graphics object (Surface)
    """
    return _go.Surface(
        x=points.x,
        y=points.y,
        z=points.z,
        surfacecolor=field.n,
        colorscale=colorscale,
        opacity=opacity,
        showscale=showscale,
        cmin=cmin,
        cmax=cmax,
        colorbar=dict(title="|B| (" + field.unit + ")"),
    )


def _draw_mesh(vertices, magnet_opacity, magnet_color):
    """Generates mesh of magnet for drawing in plotly plot

    Args:
        vertices (ndarray): array of vertices
        magnet_opacity (float): opacity of object

    Returns:
        Mesh3d object: plotly graphics object (Mesh3d)
    """
    return _go.Mesh3d(
        x=vertices[0],
        y=vertices[1],
        z=vertices[2],
        alphahull=0,
        opacity=magnet_opacity,
        color=magnet_color,
        # vertexcolor=(0, 0, 0),
        flatshading=True,
        hoverinfo="skip",
        name="y",
        showscale=False,
    )


def _draw_cones(points, field, NA=10, cone_opacity=1.0):
    """Generates cones object for drawing vector field in 3D in plotly.
    Arrows are white and normalised.

    Args:
        x (ndarray): x-coordinates
        y (ndarray): y-coordinates
        z (ndarray): z-coordinates
        B (Vector3): Magnetic field vectors for plotting
        NA (int, optional): step size, number of elements to skip for drawing cones. Defaults to 10.
        cone_opacity (float, optional): Cone opacity. Defaults to 1.0.

    Returns:
        Cone object: plotly graphics object (Cone)
    """
    x = points.x.reshape(field.x.shape)
    y = points.y.reshape(field.y.shape)
    z = points.z.reshape(field.z.shape)
    maxB = field.n[::NA, ::NA].ravel()
    cscale = [[0, "rgb(255,255,255)"], [1, "rgb(255,255,255)"]]
    data_object = _go.Cone(
        x=x[::NA, ::NA].ravel(),
        y=y[::NA, ::NA].ravel(),
        z=z[::NA, ::NA].ravel(),
        u=field.x[::NA, ::NA].ravel() / maxB,
        v=field.y[::NA, ::NA].ravel() / maxB,
        w=field.z[::NA, ::NA].ravel() / maxB,
        sizemode="scaled",
        colorscale=cscale,
        showscale=False,
        opacity=cone_opacity,
    )
    return data_object


def _generate_all_meshes(magnet_opacity=1.0):

    data_objects = []
    for magnet in Magnet_3D.instances:
        if issubclass(magnet.__class__, Mesh):
            mesh_data = Graphic_Mesh(magnet.mesh_vectors, magnet.start, magnet.stop)
            data_objects.append(mesh_data.generate_vertices())

        else:
            if issubclass(magnet.__class__, Prism):
                polyhed = Graphic_Cuboid(
                    center=magnet.center,
                    size=magnet.get_size(),
                    alpha=magnet.alpha,
                    beta=magnet.beta,
                    gamma=magnet.gamma,
                )
            elif issubclass(magnet.__class__, Cylinder):
                polyhed = Graphic_Cylinder(
                    center=magnet.center,
                    radius=magnet.radius,
                    length=magnet.length,
                    alpha=magnet.alpha,
                    beta=magnet.beta,
                    gamma=magnet.gamma,
                )

            elif issubclass(magnet.__class__, Sphere):
                polyhed = Graphic_Sphere(
                    center=magnet.center,
                    radius=magnet.radius,
                    alpha=magnet.alpha,
                    beta=magnet.beta,
                    gamma=magnet.gamma,
                )

            data_objects.append(
                _draw_mesh(polyhed.vertices, magnet_opacity, polyhed.color)
            )

    return data_objects


def _generate_volume_data(points, field, **kwargs):

    cmin = kwargs.pop("cmin", 0.0)
    cmax = kwargs.pop("cmax", 0.5)
    colorscale = kwargs.pop("colorscale", "viridis")

    opacityscale = kwargs.pop("opacityscale", None)

    caps = kwargs.pop("no_caps", False)
    if caps:
        caps = dict(x_show=False, y_show=False, z_show=False)
    else:
        caps = dict(x_show=True, y_show=True, z_show=True)

    if type(opacityscale) is str:
        if opacityscale.lower() == "normal":
            opacityscale = [
                [cmin, 0],
                [(cmax - cmin) / 4, 0.5],
                [0.2, 0],
                [cmax, 1],
            ]
        elif opacityscale.lower() == "invert":
            opacityscale = [
                [cmin, 1],
                [(cmax - cmin) * 3 / 4, 0.5],
                [0.2, 0],
                [cmax, 0],
            ]

    return _go.Volume(
        x=points.x.flatten(),
        y=points.y.flatten(),
        z=points.z.flatten(),
        value=field.n.flatten(),
        colorscale=colorscale,
        cmin=cmin,
        cmax=cmax,
        isomin=cmin,
        isomax=cmax,
        # opacity needs to be small to see through all surfaces
        opacity=kwargs.pop("opacity", 0.1),
        opacityscale=opacityscale,
        surface_count=kwargs.pop("num_levels", 10),
        showscale=True,
        caps=caps,
        colorbar=dict(title="|B| (" + field.unit + ")"),
    )


def surface_slice3(**kwargs):
    """Calculates and plots magnetic field slices a

    Returns:
        dictionary: cached data for each plane with potential keys: 'xy', 'xz', 'yz'
        containing subdictionaries, whose keys are 'x','y', 'z', 'B'.
    """

    reset_polyhedra()

    max1 = kwargs.pop("max1", 30)
    max2 = kwargs.pop("max", 30)
    slice_value = kwargs.pop("slice_value", 0.0)
    unit = kwargs.pop("unit", "mm")

    opacity = kwargs.pop("opacity", 0.1)
    magnet_opacity = kwargs.pop("magnet_opacity", 1.0)
    cone_opacity = kwargs.pop("cone_opacity", 1.0)
    planes = kwargs.pop("planes", ["xy", "xz", "yz"])

    num_arrows = kwargs.pop("num_arrows", 20)
    num_points = kwargs.pop("num_arrows", 100)

    NA = num_points // num_arrows

    if NA < 1:
        NA = 1
    cmin = kwargs.pop("cmin", 0)
    cmax = kwargs.pop("cmax", 0.5)
    colorscale = kwargs.pop("colorscale", "viridis")

    data_objects = []
    cache = {}

    show_magnets = kwargs.pop("show_magnets", True)

    if show_magnets:
        data_objects.extend(_generate_all_meshes(magnet_opacity=magnet_opacity))

    for plane in planes:
        points = slice3D(
            plane=plane,
            max1=max1,
            max2=max2,
            slice_value=slice_value,
            unit=unit,
            num_points=num_points,
        )
        field = B_calc_3D(points)

        cache[plane] = {"points": points, "field": field}

        data_objects.append(
            _draw_surface_slice(
                points,
                field,
                colorscale,
                opacity=opacity,
                cmin=cmin,
                cmax=cmax,
                showscale=True,
            )
        )
        data_objects.append(
            _draw_cones(points, field, NA=NA, cone_opacity=cone_opacity)
        )

    fig = _go.Figure(data=data_objects)

    fig.update_layout(
        scene=dict(
            xaxis_title="x (" + points.unit + ")",
            yaxis_title="y (" + points.unit + ")",
            zaxis_title="z (" + points.unit + ")",
        ),
        width=700,
        margin=dict(r=20, b=10, l=10, t=10),
    )
    fig.show()
    return fig, cache, data_objects


def volume_plot(points, field, **kwargs):
    """Calculates and plots magnetic field vol

    Returns:
        dictionary: cached data for each plane with potential keys: 'xy', 'xz', 'yz'
        containing subdictionaries, whose keys are 'x','y', 'z', 'B'.
    """

    reset_polyhedra()

    opacity = kwargs.pop("opacity", 1)
    opacityscale = kwargs.pop("opacityscale", None)
    magnet_opacity = kwargs.pop("magnet_opacity", 1.0)
    cone_opacity = kwargs.pop("cone_opacity", 1.0)

    num_points = kwargs.pop("num_points", None)

    num_arrows = kwargs.pop("num_arrows", None)

    cmin = kwargs.pop("cmin", 0)
    cmax = kwargs.pop("cmax", 0.5)
    num_levels = kwargs.pop("num_levels", 5)

    show_magnets = kwargs.pop("show_magnets", True)

    data_objects = []

    if show_magnets:
        data_objects.extend(_generate_all_meshes(magnet_opacity=magnet_opacity))

    colorscale = kwargs.pop("colorscale", "viridis")

    #     kernel_size = 1
    #     kernel = np.ones([kernel_size, kernel_size, kernel_size]) / kernel_size
    #     B.n = ndimage.convolve(B.n, kernel)

    data_objects.append(
        _generate_volume_data(
            points,
            field,
            cmim=cmin,
            cmax=cmax,
            opacity=opacity,
            colorscale=colorscale,
            num_levels=num_levels,
            opacityscale=opacityscale,
        )
    )

    if num_arrows is not None:
        NA = num_points // num_arrows
        if NA < 1:
            NA = 1
        data_objects.append(
            _draw_cones(points, field, NA=NA, cone_opacity=cone_opacity)
        )

    fig = _go.Figure(data=data_objects)

    fig.update_layout(
        scene=dict(
            xaxis_title="x (" + points.unit + ")",
            yaxis_title="y (" + points.unit + ")",
            zaxis_title="z (" + points.unit + ")",
        ),
        width=700,
        margin=dict(r=20, b=10, l=10, t=10),
    )
    fig.show()

    return data_objects, fig


def volume_calculate_plot(**kwargs):

    """Calculates and plots magnetic field vol

    Returns:
        dictionary: cached data for each plane with potential keys: 'xy', 'xz', 'yz'
        containing subdictionaries, whose keys are 'x','y', 'z', 'B'.
    """

    num_points = kwargs.pop("num_points", 30)

    unit = kwargs.pop("unit", "mm")

    xmax = kwargs.pop("xmax", 30)
    ymax = kwargs.pop("ymax", 30)
    zmax = kwargs.pop("zmax", 30)

    xmin = kwargs.pop("xmin", -1 * xmax)
    ymin = kwargs.pop("ymin", -1 * ymax)
    zmin = kwargs.pop("zmin", -1 * zmax)

    points = grid3D(
        xmax,
        ymax,
        zmax,
        num_points=num_points,
        xmin=xmin,
        ymin=ymin,
        zmin=zmin,
        unit=unit,
    )
    field = B_calc_3D(points)

    data_objects, fig = volume_plot(points, field, num_points=num_points, **kwargs)

    return fig, data_objects, points, field