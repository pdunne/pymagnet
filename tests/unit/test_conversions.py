# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for coordinate conversion functions."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.utils._conversions import (
    cart2pol,
    cart2sph,
    get_unit_value_meter,
    get_unit_value_tesla,
    pol2cart,
    sph2cart,
    sphere_sph2cart,
    vector_pol2cart,
    vector_sph2cart,
)


class TestCart2Pol:
    """Tests for cartesian to polar conversion."""

    def test_origin(self):
        """Origin: (0,0) -> rho=0, phi=0."""
        rho, phi = cart2pol(0.0, 0.0)
        assert rho == 0.0
        assert phi == 0.0

    def test_positive_x_axis(self):
        """(1,0) -> rho=1, phi=0."""
        rho, phi = cart2pol(1.0, 0.0)
        npt.assert_allclose(rho, 1.0)
        npt.assert_allclose(phi, 0.0)

    def test_positive_y_axis(self):
        """(0,1) -> rho=1, phi=pi/2."""
        rho, phi = cart2pol(0.0, 1.0)
        npt.assert_allclose(rho, 1.0)
        npt.assert_allclose(phi, np.pi / 2)

    def test_negative_x_axis(self):
        """(-1,0) -> rho=1, phi=pi."""
        rho, phi = cart2pol(-1.0, 0.0)
        npt.assert_allclose(rho, 1.0)
        npt.assert_allclose(phi, np.pi)

    def test_negative_y_axis(self):
        """(0,-1) -> rho=1, phi=-pi/2."""
        rho, phi = cart2pol(0.0, -1.0)
        npt.assert_allclose(rho, 1.0)
        npt.assert_allclose(phi, -np.pi / 2)

    @pytest.mark.parametrize(
        "x,y,expected_phi",
        [
            (1, 1, np.pi / 4),  # First quadrant
            (-1, 1, 3 * np.pi / 4),  # Second quadrant
            (-1, -1, -3 * np.pi / 4),  # Third quadrant
            (1, -1, -np.pi / 4),  # Fourth quadrant
        ],
    )
    def test_quadrants(self, x, y, expected_phi):
        """Test all four quadrants."""
        rho, phi = cart2pol(float(x), float(y))
        npt.assert_allclose(rho, np.sqrt(2), rtol=1e-10)
        npt.assert_allclose(phi, expected_phi, rtol=1e-10)

    def test_array_input(self):
        """Test with numpy array inputs."""
        x = np.array([1.0, 0.0, -1.0, 0.0])
        y = np.array([0.0, 1.0, 0.0, -1.0])
        rho, _phi = cart2pol(x, y)
        npt.assert_allclose(rho, [1.0, 1.0, 1.0, 1.0])


class TestPol2Cart:
    """Tests for polar to cartesian conversion."""

    def test_zero_radius(self):
        """rho=0 gives origin regardless of phi."""
        x, y = pol2cart(0.0, np.pi / 4)
        npt.assert_allclose(x, 0.0, atol=1e-15)
        npt.assert_allclose(y, 0.0, atol=1e-15)

    def test_unit_circle_0(self):
        """Point on unit circle at phi=0."""
        x, y = pol2cart(1.0, 0.0)
        npt.assert_allclose(x, 1.0)
        npt.assert_allclose(y, 0.0, atol=1e-15)

    def test_unit_circle_pi_2(self):
        """Point on unit circle at phi=pi/2."""
        x, y = pol2cart(1.0, np.pi / 2)
        npt.assert_allclose(x, 0.0, atol=1e-15)
        npt.assert_allclose(y, 1.0)

    def test_unit_circle_pi(self):
        """Point on unit circle at phi=pi."""
        x, y = pol2cart(1.0, np.pi)
        npt.assert_allclose(x, -1.0)
        npt.assert_allclose(y, 0.0, atol=1e-15)

    def test_roundtrip(self):
        """cart2pol then pol2cart should return original."""
        x_orig, y_orig = 3.0, 4.0
        rho, phi = cart2pol(x_orig, y_orig)
        x_back, y_back = pol2cart(rho, phi)
        npt.assert_allclose(x_back, x_orig)
        npt.assert_allclose(y_back, y_orig)

    def test_array_roundtrip(self):
        """Roundtrip with arrays."""
        x_orig = np.array([1.0, 2.0, 3.0])
        y_orig = np.array([4.0, 5.0, 6.0])
        rho, phi = cart2pol(x_orig, y_orig)
        x_back, y_back = pol2cart(rho, phi)
        npt.assert_allclose(x_back, x_orig)
        npt.assert_allclose(y_back, y_orig)


class TestVectorPol2Cart:
    """Tests for vector field polar to cartesian conversion."""

    def test_radial_field_at_phi_0(self):
        """Radial field at phi=0: Br -> Bx."""
        Bx, By = vector_pol2cart(1.0, 0.0, 0.0)
        npt.assert_allclose(Bx, 1.0)
        npt.assert_allclose(By, 0.0, atol=1e-15)

    def test_azimuthal_field_at_phi_0(self):
        """Azimuthal field at phi=0: Bphi -> By."""
        Bx, By = vector_pol2cart(0.0, 1.0, 0.0)
        npt.assert_allclose(Bx, 0.0, atol=1e-15)
        npt.assert_allclose(By, 1.0)

    def test_radial_field_at_phi_90(self):
        """Radial field at phi=pi/2: Br -> By."""
        Bx, By = vector_pol2cart(1.0, 0.0, np.pi / 2)
        npt.assert_allclose(Bx, 0.0, atol=1e-15)
        npt.assert_allclose(By, 1.0)


class TestCart2Sph:
    """Tests for cartesian to spherical conversion."""

    def test_origin(self):
        """Origin handling (r=0)."""
        r, _theta, _phi = cart2sph(0.0, 0.0, 0.0)
        assert r == 0.0
        # theta is NaN due to arccos(0/0), this is expected

    def test_positive_z_axis(self):
        """(0,0,1) -> r=1, theta=0, phi=0."""
        r, theta, _phi = cart2sph(0.0, 0.0, 1.0)
        npt.assert_allclose(r, 1.0)
        npt.assert_allclose(theta, 0.0)
        # phi is undefined on z-axis

    def test_negative_z_axis(self):
        """(0,0,-1) -> r=1, theta=pi."""
        r, theta, _phi = cart2sph(0.0, 0.0, -1.0)
        npt.assert_allclose(r, 1.0)
        npt.assert_allclose(theta, np.pi)

    def test_positive_x_axis(self):
        """(1,0,0) -> r=1, theta=pi/2, phi=0."""
        r, theta, phi = cart2sph(1.0, 0.0, 0.0)
        npt.assert_allclose(r, 1.0)
        npt.assert_allclose(theta, np.pi / 2)
        npt.assert_allclose(phi, 0.0)

    def test_positive_y_axis(self):
        """(0,1,0) -> r=1, theta=pi/2, phi=pi/2."""
        r, theta, phi = cart2sph(0.0, 1.0, 0.0)
        npt.assert_allclose(r, 1.0)
        npt.assert_allclose(theta, np.pi / 2)
        npt.assert_allclose(phi, np.pi / 2)

    def test_diagonal_point(self):
        """Diagonal point in positive octant."""
        x = y = z = 1.0 / np.sqrt(3)
        r, theta, phi = cart2sph(x, y, z)
        npt.assert_allclose(r, 1.0, rtol=1e-10)
        npt.assert_allclose(theta, np.arccos(1 / np.sqrt(3)), rtol=1e-10)
        npt.assert_allclose(phi, np.pi / 4, rtol=1e-10)


class TestSph2Cart:
    """Tests for spherical to cartesian conversion."""

    def test_north_pole(self):
        """theta=0 should give point on positive z-axis."""
        x, y, z = sph2cart(1.0, 0.0, 0.0)
        npt.assert_allclose(x, 0.0, atol=1e-15)
        npt.assert_allclose(y, 0.0, atol=1e-15)
        npt.assert_allclose(z, 1.0)

    def test_south_pole(self):
        """theta=pi should give point on negative z-axis."""
        x, y, z = sph2cart(1.0, np.pi, 0.0)
        npt.assert_allclose(x, 0.0, atol=1e-15)
        npt.assert_allclose(y, 0.0, atol=1e-15)
        npt.assert_allclose(z, -1.0)

    def test_equator_phi_0(self):
        """theta=pi/2, phi=0: positive x-axis."""
        x, y, z = sph2cart(1.0, np.pi / 2, 0.0)
        npt.assert_allclose(x, 1.0)
        npt.assert_allclose(y, 0.0, atol=1e-15)
        npt.assert_allclose(z, 0.0, atol=1e-15)

    def test_equator_phi_90(self):
        """theta=pi/2, phi=pi/2: positive y-axis."""
        x, y, z = sph2cart(1.0, np.pi / 2, np.pi / 2)
        npt.assert_allclose(x, 0.0, atol=1e-15)
        npt.assert_allclose(y, 1.0)
        npt.assert_allclose(z, 0.0, atol=1e-15)

    def test_roundtrip(self):
        """cart2sph then sph2cart returns original."""
        x_orig, y_orig, z_orig = 1.0, 2.0, 3.0
        r, theta, phi = cart2sph(x_orig, y_orig, z_orig)
        x_back, y_back, z_back = sph2cart(r, theta, phi)
        npt.assert_allclose(x_back, x_orig, rtol=1e-10)
        npt.assert_allclose(y_back, y_orig, rtol=1e-10)
        npt.assert_allclose(z_back, z_orig, rtol=1e-10)

    def test_zero_radius(self):
        """r=0 gives origin."""
        x, y, z = sph2cart(0.0, np.pi / 4, np.pi / 4)
        npt.assert_allclose(x, 0.0, atol=1e-15)
        npt.assert_allclose(y, 0.0, atol=1e-15)
        npt.assert_allclose(z, 0.0, atol=1e-15)


class TestVectorSph2Cart:
    """Tests for vector field spherical to cartesian conversion."""

    def test_radial_at_north_pole(self):
        """Radial field at north pole: Br -> Bz."""
        Bx, By, Bz = vector_sph2cart(1.0, 0.0, 0.0, 0.0, 0.0)
        npt.assert_allclose(Bx, 0.0, atol=1e-15)
        npt.assert_allclose(By, 0.0, atol=1e-15)
        npt.assert_allclose(Bz, 1.0)

    def test_radial_on_equator(self):
        """Radial field on equator at phi=0: Br -> Bx."""
        Bx, By, Bz = vector_sph2cart(1.0, 0.0, 0.0, np.pi / 2, 0.0)
        npt.assert_allclose(Bx, 1.0)
        npt.assert_allclose(By, 0.0, atol=1e-15)
        npt.assert_allclose(Bz, 0.0, atol=1e-15)


class TestSphereSph2Cart:
    """Tests for sphere_sph2cart (no Bphi component)."""

    def test_radial_at_north_pole(self):
        """Radial field at north pole."""
        Bx, By, Bz = sphere_sph2cart(1.0, 0.0, 0.0, 0.0)
        npt.assert_allclose(Bx, 0.0, atol=1e-15)
        npt.assert_allclose(By, 0.0, atol=1e-15)
        npt.assert_allclose(Bz, 1.0)

    def test_theta_component_on_equator(self):
        """Theta field on equator at phi=0: Btheta -> -Bz."""
        Bx, By, Bz = sphere_sph2cart(0.0, 1.0, np.pi / 2, 0.0)
        npt.assert_allclose(Bx, 0.0, atol=1e-15)
        npt.assert_allclose(By, 0.0, atol=1e-15)
        npt.assert_allclose(Bz, -1.0)


class TestUnitValueMeter:
    """Tests for get_unit_value_meter function."""

    @pytest.mark.parametrize(
        "unit,expected",
        [
            ("m", 1.0),
            ("cm", 1e-2),
            ("mm", 1e-3),
            ("um", 1e-6),
            ("µm", 1e-6),
            ("nm", 1e-9),
            ("km", 1e3),
            ("Mm", 1e6),
        ],
    )
    def test_length_units(self, unit, expected):
        """Test get_unit_value_meter for various units."""
        result = get_unit_value_meter(unit)
        assert result == expected

    def test_invalid_unit_returns_none(self):
        """Invalid unit should return None."""
        result = get_unit_value_meter("invalid")
        assert result is None

    def test_empty_string_returns_none(self):
        """Empty string should return None."""
        result = get_unit_value_meter("")
        assert result is None


class TestUnitValueTesla:
    """Tests for get_unit_value_tesla function."""

    @pytest.mark.parametrize(
        "unit,expected",
        [
            ("T", 1.0),
            ("mT", 1e-3),
            ("uT", 1e-6),
            ("µT", 1e-6),
            ("nT", 1e-9),
            ("kT", 1e3),
            ("GT", 1e9),
        ],
    )
    def test_field_units(self, unit, expected):
        """Test get_unit_value_tesla for various units."""
        result = get_unit_value_tesla(unit)
        assert result == expected

    def test_invalid_unit_returns_none(self):
        """Invalid unit should return None."""
        result = get_unit_value_tesla("invalid")
        assert result is None

    def test_length_unit_returns_none(self):
        """Length unit should return None for Tesla function."""
        result = get_unit_value_tesla("mm")
        assert result is None
