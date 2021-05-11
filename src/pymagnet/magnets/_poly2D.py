"""2D Polygon Magnet class


TODO:
    * Add __del__ method for removing strong ref in class instance list

"""
import numpy as _np

# from ._magnet import Magnet
from ._magnet2 import Magnet_2D


PI = _np.pi
u0 = PI * 4e-7


def _furlani(x, y, h, Kr=1):
    """[summary]

    Args:
        x ([type]): [description]
        y ([type]): [description]
        h ([type]): [description]
        Kr (int, optional): [description]. Defaults to 1.

    Returns:
        [type]: [description]
    """
    x = _np.asarray(x)
    y = _np.asarray(y)
    prefac = Kr / 4 / PI
    Bx = prefac * _np.log((x ** 2 + (y - h) ** 2) / (x ** 2 + (y + h) ** 2))
    By = 2 * prefac * _np.arctan2(2 * h * x, x ** 2 + y ** 2 - h ** 2)
    return Bx, By


class Polygon(object):
    """[summary]

    Args:
        object ([type]): [description]
    """

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
        """[summary]

        Args:
            vertex ([type]): [description]
        """
        if len(vertex) != 2:
            print("Error")
        if type(vertex) == tuple:
            self.vertices.append(vertex)
        elif len(vertex) == 2:
            self.vertices.append(tuple(vertex))
        self.set_center()

    def num_vertices(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return len(self.vertices)

    def set_center(self):
        """[summary]"""
        self.center = _np.mean(_np.asarray(self.vertices), axis=0)

    @staticmethod
    def gen_polygon(N=6, center=(0, 0), alpha=0, **kwargs):
        """[summary]

        Args:
            N (int, optional): [description]. Defaults to 6.
            apo (int, optional): [description]. Defaults to 1.
            center (tuple, optional): [description]. Defaults to (0,0).
            offset (int, optional): [description]. Defaults to 0.

        Raises:
            Exception: [description]

        Returns:
            [type]: [description]
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
        if apothem is not None:
            return apothem / _np.around(_np.cos(PI / N), 4)
        elif length is not None:
            return length / _np.around(2 * _np.sin(PI / N), 4)
        elif radius is not None:
            return radius
        else:
            raise Exception("Error, one of apothem, length, or radius must be defined.")


class LineUtils(object):
    """[summary]

    Args:
        object ([type]): [description]

    Returns:
        [type]: [description]
    """

    @staticmethod
    def unit_norm(ver1, ver2, clockwise=True):
        """[summary]

        Args:
            ver1 ([type]): [description]
            ver2 ([type]): [description]
            clockwise (bool, optional): [description]. Defaults to True.

        Returns:
            [type]: [description]
        """

        dx = ver1[0] - ver2[0]
        dy = ver1[1] - ver2[1]

        # Clockwise winding of points:
        if clockwise:
            norm = _np.array([dy, -dx])
        else:
            norm = _np.array([-dy, dx])
        length = _np.linalg.norm(norm)
        norm = norm / length
        return norm, length

    @staticmethod
    def line_center(ver1, ver2):
        """[summary]

        Args:
            ver1 ([type]): [description]
            ver2 ([type]): [description]

        Returns:
            [type]: [description]
        """

        xc = (ver1[0] + ver2[0]) / 2
        yc = (ver1[1] + ver2[1]) / 2

        return _np.array([xc, yc])

    @staticmethod
    def signed_area(polygon):
        """Calculates signed area,

        Args:
            polygon ([type]): [description]

        Returns:
            [type]: [description]
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
    """[summary]

    Args:
        object ([type]): [description]
    """

    def __init__(self, length, center, beta, K):
        self.length = length
        self.center = center
        self.beta = beta
        self.beta_rad = _np.deg2rad(beta)
        self.xc = center[0]
        self.yc = center[1]
        self.K = K

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
            + f"Length: {self.length} (m)\n"
            + f"Center {self.center} (m)\n"
            + f"Orientation: {self.beta}\n"
        )
        return str

    def calcB(self, x, y):
        """[summary]

        Args:
            x ([type]): [description]
            y ([type]): [description]

        Returns:
            [type]: [description]
        """
        from ._routines2 import rotate_points_2D, _get_field_array_shape2

        array_shape = _get_field_array_shape2(x, y)
        Bx, By = _np.zeros(array_shape), _np.zeros(array_shape)
        if _np.fabs(self.beta_rad) > 1e-4:
            xt, yt = rotate_points_2D(x - self.xc, y - self.yc, 2 * PI - self.beta_rad)
            Btx, Bty = _furlani(xt, yt, self.length / 2, self.K)
            Bx, By = rotate_points_2D(Btx, Bty, self.beta_rad)

        else:
            Bx, By = _furlani(x - self.xc, y - self.yc, self.length / 2, self.K)
        return Bx, By


class PolyMagnet(Magnet_2D):
    """2D Magnet Polygon class.

    Args:
        width [float]: magnet width
        height [float]: magnet height
        Jr [float]: Remnant magnetisation
        **kwargs: Arbitrary keyword arguments.

    kwargs:
        center [Tuple(float, float), or Point2]: center of magnet, defaults to
        Point2(0.0, 0.0)
    """

    mag_type = "PolyMagnet"

    def __init__(self, Jr, **kwargs) -> None:
        from ._routines2 import rotate_points_2D

        super().__init__(Jr, **kwargs)

        # Magnet rotation w.r.t. x-axis
        self.alpha = kwargs.pop("alpha", 0.0)
        self.alpha_radians = _np.deg2rad(self.alpha)

        self.theta = kwargs.pop("theta", 0.0)
        self.theta_radians = _np.deg2rad(self.theta)

        self.phi = kwargs.pop("phi", 90)
        self.phi_rad = _np.deg2rad(self.phi)

        self.Jx = _np.around(Jr * _np.cos(self.phi_rad), decimals=6)
        self.Jy = _np.around(Jr * _np.sin(self.phi_rad), decimals=6)
        self.tol = 1e-4  # sufficient for 0.01 degree accuracy
        self.area = None

        # self.length = kwargs.pop("length", 10e-3)
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

        self.custom_polygon = kwargs.pop("custom_polygon", False)
        self.center = kwargs.pop("center", (0, 0))
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
            vertices = _np.stack([x_rot, y_rot]).T + _np.array(self.center)
            self.polygon = Polygon(vertices=vertices.tolist())
        else:
            self.center = kwargs.pop("center", (0, 0))
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

        self.xc = self.center[0]
        self.yc = self.center[1]

    def center(self):
        """Returns magnet centre

        Returns:
            center (ndarray): numpy array [xc, yc]
        """
        return _np.array([self.xc, self.yc])

    def get_orientation(self):
        """Returns magnet orientation, `alpha` in degrees

        Returns:
            float: alpha, rotation angle w.r.t x-axis.
        """

        return self.alpha

    def _gen_sheet_magnets(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        area, norms, beta, length, center = LineUtils.signed_area(self.polygon)
        K = self.Jx * norms[:, 1] - self.Jy * norms[:, 0]
        self.area = area
        return beta, length, center, K

    def calcB(self, x, y):
        from ._routines2 import rotate_points_2D, _get_field_array_shape2

        array_shape = _get_field_array_shape2(x, y)
        Bx, By = _np.zeros(array_shape), _np.zeros(array_shape)
        beta, length, center, K = self._gen_sheet_magnets()

        if _np.fabs(self.alpha_radians) > self.tol:
            pass
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
            #     Btx, Bty = sheet.calcB(xt, yt)
            #     Btx, Bty = rotate_points_2D(Btx, Bty, 2 * PI - self.alpha_radians)
            #     Bx += Btx
            #     By += Bty

        for i in range(len(K)):
            sheet = Line(length[i], center[i], beta[i], K[i])
            Btx, Bty = sheet.calcB(x, y)
            Bx += Btx
            By += Bty
        return Bx, By
