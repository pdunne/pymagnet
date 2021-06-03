"""2D Polygon Magnet class


"""
import numpy as _np
from ._magnet2D import Magnet2D
from ..utils.global_const import PI, MAG_TOL


def _sheet_field(x, y, h, Kr=1):
    """Calculates the magnetic field due to an infinite current sheet

    Args:
        x (ndarray): x-coordinates
        y (ndarray): y-coordinates
        h (float): sheet length
        Kr (float, optional): Sheet current Density as T. Defaults to 1.

    Returns:
        tuple: Bx (ndarray), By (ndarray) magnetic field vector
    """
    x = _np.asarray(x)
    y = _np.asarray(y)
    prefac = Kr / 4 / PI
    Bx = prefac * _np.log((x ** 2 + (y - h) ** 2) / (x ** 2 + (y + h) ** 2))
    By = 2 * prefac * _np.arctan2(2 * h * x, x ** 2 + y ** 2 - h ** 2)
    return Bx, By


class Polygon(object):
    """Polygon class for generating list of vertices"""

    def __init__(self, **kwargs):
        vertices = kwargs.pop("vertices", None)
        center = kwargs.pop("center", None)
        if vertices is not None:
            if type(vertices) is _np.ndarray:
                if center is not None:
                    vertices += center
                self.vertices = vertices.tolist()
            else:
                self.vertices = vertices

            if center is not None:
                self.center = center
            else:
                self.set_center()
        else:
            self.vertices = []
            self.center = _np.NaN

    def append(self, vertex):
        """Appends vertex to list of vertices

        Args:
            vertex (list): list of vertices
        """
        if len(vertex) != 2:
            print("Error")
        if type(vertex) == tuple:
            self.vertices.append(vertex)
        elif len(vertex) == 2:
            self.vertices.append(tuple(vertex))
        self.set_center()

    def num_vertices(self):
        """Gets number of vertices

        Returns:
            int: number of vertices
        """
        return len(self.vertices)

    def set_center(self):
        """Sets center of polygon to be centroid"""
        self.center = _np.mean(_np.asarray(self.vertices), axis=0)

    @staticmethod
    def gen_polygon(N=6, center=(0.0, 0.0), alpha=0.0, **kwargs):
        """Generates regular polygon. One of apothem, side length or radius must be defined.

        Args:
            N (int, optional): Number of sides. Defaults to 6.
            center (tuple, optional): Polygon center. Defaults to (0.0, 0.0).
            alpha (float, optional): Orientration with respect to x-axis. Defaults to 0.0.

        Raises:
            Exception: N must be > 2

        Returns:
            ndarray: polygon vertices
        """
        N = int(N)

        if N < 3:
            raise Exception("Error, N must be > 2.")

        apothem = kwargs.pop("apothem", None)
        length = kwargs.pop("length", None)
        radius = kwargs.pop("radius", None)

        radius = Polygon.check_radius(N, apothem, length, radius)

        k = _np.arange(0, N, 1)
        xc = center[0]
        yc = center[1]

        def f(N):
            if N % 2 == 0:
                return PI / N + _np.deg2rad(alpha)
            else:
                return PI / N + PI + _np.deg2rad(alpha)

        xv = xc + radius * _np.sin(2 * PI * k / N + f(N))
        yv = yc + radius * _np.cos(2 * PI * k / N + f(N))
        poly_verts = _np.vstack((xv, yv)).T.tolist()

        return poly_verts

    @staticmethod
    def check_radius(N, apothem, length, radius):
        """Checks which of apothem, side length, or radius has been passed as kwargs
        to `gen_polygon()`. Order of precendence is apothem, length, radius.

        Args:
            N (int): Number of sides
            apothem (float): polygon apothem
            length (float): side length
            radius (float): outcircle radius

        Raises:
            Exception: One of apothem, length, or raduis must be defined

        Returns:
            float: returns radius
        """
        if apothem is not None:
            return apothem / _np.around(_np.cos(PI / N), 4)
        elif length is not None:
            return length / _np.around(2 * _np.sin(PI / N), 4)
        elif radius is not None:
            return radius
        else:
            raise Exception("Error, one of apothem, length, or radius must be defined.")


class LineUtils(object):
    """Utility class consisting of rountines for 2D line elements"""

    @staticmethod
    def unit_norm(vertex_1, vertex_2, clockwise=True):
        """Get unit normal to vertex

        Args:
            vertex_1 (ndarray): vertex 1
            vertex_2 (ndarray): vertex 2
            clockwise (bool, optional): Clockwise orientation of points. Defaults to True.

        Returns:
            tuple: normal vector (ndarray), length i.e. distance between vertices (float)
        """

        dx = vertex_1[0] - vertex_2[0]
        dy = vertex_1[1] - vertex_2[1]

        # Clockwise winding of points:
        if clockwise:
            norm = _np.array([dy, -dx])
        else:
            norm = _np.array([-dy, dx])
        length = _np.linalg.norm(norm)
        norm = norm / length
        return norm, length

    @staticmethod
    def line_center(vertex_1, vertex_2):
        """Gets midpoint of two vertices

        Args:
            vertex_1 (ndarray): vertex 1
            vertex_2 (ndarray): vertex 2

        Returns:
            ndarray: midpoint
        """
        xc = (vertex_1[0] + vertex_2[0]) / 2
        yc = (vertex_1[1] + vertex_2[1]) / 2

        return _np.array([xc, yc])

    @staticmethod
    def signed_area(polygon):
        """Calculates signed area of a polygon

        Args:
            polygon (Polygon): Polygon instance

        Returns:
            float: signed area
        """
        j = 1
        NP = polygon.num_vertices()
        area = 0
        norm = _np.empty([NP, 2])
        center = _np.empty([NP, 2])
        beta = _np.empty(NP)  # angle w.r.t. y axis
        length = _np.empty(NP)

        for i in range(NP):
            j = j % NP
            area += (polygon.vertices[j][0] - polygon.vertices[i][0]) * (
                polygon.vertices[j][1] + polygon.vertices[i][1]
            )
            norm[i, :], length[i] = LineUtils.unit_norm(
                polygon.vertices[i], polygon.vertices[j]
            )
            center[i, :] = LineUtils.line_center(
                polygon.vertices[i], polygon.vertices[j]
            )
            j += 1

        # check winding order of polygon, area < 0 for clockwise ordering of points
        if area < 0:
            norm *= -1
        beta[:] = _np.rad2deg(_np.arctan2(norm[:, 1], norm[:, 0]))

        return area / 2.0, norm, beta, length, center


class Line(object):
    """Line Class for storing properties of a sheet manget"""

    def __init__(self, length, center, beta, K):
        """Init Method

        Args:
            length (float): side length
            center (ndarray): magnet center (x, y)
            beta (float): Orientation w.r.t. z-axis in degrees
            K (float): Sheet current density in tesla
        """
        self.length = length
        self.center = center
        self.beta = beta
        self.beta_rad = _np.deg2rad(beta)
        self.xc = center[0]
        self.yc = center[1]
        self.K = K
        self.tol = MAG_TOL

    def __str__(self):
        str = (
            f"K: {self.K} (T)\n"
            + f"Length: {self.length} (m)\n"
            + f"Center {self.center} (m)\n"
            + f"Orientation: {self.beta}\n"
        )
        return str

    def __repr__(self):
        str = (
            f"K: {self.K} (T)\n"
            + f"Length: {self.length}\n"
            + f"Center {self.center}\n"
            + f"Orientation: {self.beta}\n"
        )
        return str

    def get_center(self):
        """Returns line center

        Returns:
            ndarray: center (x,y)
        """

        return self.center

    def get_field(self, x, y):
        """Calculates the magnetic field due to a sheet magnet
        First a transformation into the local coordinates is made, the field calculated
        and then the magnetic field it rotated out to the global coordinates

        Args:
            x (ndarray): x-coordinates
            y (ndarray): y-coordinates

        Returns:
            tuple: Bx (ndarray), By (ndarray) magnetic field vector
        """
        from ..utils._routines2D import rotate_points_2D, _get_field_array_shape2

        array_shape = _get_field_array_shape2(x, y)
        Bx, By = _np.zeros(array_shape), _np.zeros(array_shape)
        if _np.fabs(self.beta_rad) > self.tol:
            xt, yt = rotate_points_2D(x - self.xc, y - self.yc, 2 * PI - self.beta_rad)
            Btx, Bty = _sheet_field(xt, yt, self.length / 2, self.K)
            Bx, By = rotate_points_2D(Btx, Bty, self.beta_rad)

        else:
            Bx, By = _sheet_field(x - self.xc, y - self.yc, self.length / 2, self.K)
        return Bx, By


class PolyMagnet(Magnet2D):
    """2D Magnet Polygon class."""

    mag_type = "PolyMagnet"

    def __init__(self, Jr, **kwargs) -> None:
        """Init method

        NOTE:
            * When creating a regular polygon, one of apothem, radius, or length must be defined as a kwarg or an exception will be raised.
            * When creating a regular polygon, the number of sides `num_sides` must be at least 3 or an exception will be raised.
            * When creating a custom polygon at least one vertex pair must be defined with `vertices` or an exception will be raised.

        Args:
            Jr (float): signed magnitude of remnant magnetisation

        Kwargs:
            alpha (float): Not used
            theta (float): Orientation of magnet w.r.t x-axis of magnet
            phi (float): Orientation of magnetisation w.r.t x-axis of magnet in degrees. Defaults to 90.0.
            center (ndarray): magnet center (x, y). Defaults to (0.0, 0.0)
            length (float): side length if creating a regular polygon
            apothem (float): apothem (incircle radius) if creating a regular polygon
            radius (float): radius (circumcircle radius) if  creating a regular polygon
            num_sides (int): number of sides of a regular polygon. Defaults to 6.
            custom_polygon (bool): Flag to define a custom polygon. Defaults to False.
            vertices (ndarray, list): List of custom vertices. Defaults to None.

        Raises:
            Exception: If creating a custom polygon, `vertices` must not be None.
        """
        from ..utils._routines2D import rotate_points_2D

        super().__init__(Jr, **kwargs)

        # Magnet rotation w.r.t. x-axis
        self.alpha = kwargs.pop("alpha", 0.0)
        self.alpha_radians = _np.deg2rad(self.alpha)

        self.theta = kwargs.pop("theta", 0.0)
        self.theta_radians = _np.deg2rad(self.theta)

        self.phi = kwargs.pop("phi", 90.0)
        self.phi_rad = _np.deg2rad(self.phi)

        self.Jx = _np.around(Jr * _np.cos(self.phi_rad), decimals=6)
        self.Jy = _np.around(Jr * _np.sin(self.phi_rad), decimals=6)
        self.tol = MAG_TOL
        self.area = None

        self.custom_polygon = kwargs.pop("custom_polygon", False)

        self.center = kwargs.pop("center", _np.array([0.0, 0.0, 0.0]))
        self.center = _np.asarray(self.center)

        if self.custom_polygon:
            vertices = kwargs.pop("vertices", None)
            if vertices is None:
                raise Exception("Error, no vertices were defined.")

            vertices = _np.atleast_2d(vertices)

            x_rot, y_rot = rotate_points_2D(
                vertices[:, 0],
                vertices[:, 1],
                self.theta_radians,  # + self.alpha_radians,
            )
            vertices = _np.stack([x_rot, y_rot]).T + self.center
            self.polygon = Polygon(vertices=vertices.tolist())
        else:
            self.length = kwargs.pop("length", None)
            self.apothem = kwargs.pop("apothem", None)
            self.radius = kwargs.pop("radius", None)
            self.num_sides = kwargs.pop("num_sides", 6)

            self.radius = Polygon.check_radius(
                self.num_sides,
                self.apothem,
                self.length,
                self.radius,
            )
            # Generate Polygon
            self.polygon = Polygon(
                vertices=Polygon.gen_polygon(
                    self.num_sides,
                    self.center,
                    self.theta,  # + self.alpha,
                    length=self.length,
                    apothem=self.apothem,
                    radius=self.radius,
                ),
                center=self.center,
            )

    def get_center(self):
        """Returns magnet centre

        Returns:
            center (ndarray): numpy array [xc, yc]
        """
        return self.center

    def get_orientation(self):
        """Returns magnet orientation, `alpha` in degrees

        Returns:
            float: alpha, rotation angle w.r.t x-axis.
        """

        return self.alpha

    def _gen_sheet_magnets(self):
        """Generates orientation, size, and centre of sheet magnets for a given
        polygon

        Returns:
            tuple: beta (ndarray), length (ndarray), centre (ndarray), K (ndarray) - sheet current density in tesla.
        """
        area, norms, beta, length, center = LineUtils.signed_area(self.polygon)
        K = self.Jx * norms[:, 1] - self.Jy * norms[:, 0]
        self.area = area
        return beta, length, center, K

    def get_field(self, x, y):
        """Calculates the magnetic field of a polygon due to each face

        Args:
            x (ndarray): x-coordinates
            y (ndarray): y-coordinates

        Returns:
            tuple: Bx (ndarray), By (ndarray) magnetic field vector
        """
        from ..utils._routines2D import rotate_points_2D, _get_field_array_shape2

        array_shape = _get_field_array_shape2(x, y)
        Bx, By = _np.zeros(array_shape), _np.zeros(array_shape)
        beta, length, center, K = self._gen_sheet_magnets()

        if _np.fabs(self.alpha_radians) > self.tol:
            pass
            print("Arbitrary rotation with alpha not yet implemented!!")

            # FIXME: rotate centres
            # xt, yt = rotate_points_2D(x - self.xc, y - self.yc, self.alpha_radians)
            # beta += self.alpha
            # xc_rot, yc_rot = rotate_points_2D(
            #     center[:, 0] - self.xc,
            #     center[:, 1] - self.yc,
            #     self.alpha_radians,
            # )
            # center[:, 0] = xc_rot
            # center[:, 1] = yc_rot
            #
            #
            # for i in range(len(K)):
            #     sheet = Line(length[i], center[i], beta[i], K[i])
            #     Btx, Bty = sheet.get_field(xt, yt)
            #     Btx, Bty = rotate_points_2D(Btx, Bty, 2 * PI - self.alpha_radians)
            #     Bx += Btx
            #     By += Bty

        for i in range(len(K)):
            sheet = Line(length[i], center[i], beta[i], K[i])
            Btx, Bty = sheet.get_field(x, y)
            Bx += Btx
            By += Bty
        return Bx, By
