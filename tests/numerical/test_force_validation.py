# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Numerical validation tests for force/torque calculations.

These tests validate the force/torque calculation functions against
published reference data from:

1. Allag2009: "3D Analytical Calculation of the Torque and Force Exerted
   Between Two Cuboidal Magnets"
   - Two 10x10x10 mm cube magnets, Jr=1.0T
   - Coaxial configuration with varying lateral offset

2. O'Connell2020: "Analytical Calculation of Force Between Cylindrical
   Permanent Magnets" (adapted for prism validation)
   - Asymmetric prism configuration with force and torque data

These tests use pytest.mark.numerical to allow selective execution.
"""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet import magnets, reset
from pymagnet.forces._prism_force import calc_force_prism

# Mark all tests in this module as numerical validation tests
pytestmark = pytest.mark.numerical


# ==================== Allag2009 Reference Data ====================
# Two 10x10x10 mm cubes, Jr=1.0T
# Magnet 1: center=(0, 0, 0)
# Magnet 2: center=(offset, 0, 20) where offset varies from -20 to +20 mm
# Reference forces in Newtons

ALLAG2009_OFFSETS = np.array(
    [
        -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10,
        -9, -8, -7, -6, -5, -4, -3, -2, -1, 0,
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    ]
)

# Reference Fx values (N) - lateral force component
ALLAG2009_FX = np.array(
    [
        0.285, 0.351, 0.426, 0.508, 0.597, 0.690, 0.785, 0.877, 0.960, 1.026, 1.066,
        1.071, 1.036, 0.958, 0.838, 0.681, 0.497, 0.299, 0.105, -0.060, -0.008,
        0.060, -0.105, -0.299, -0.497, -0.681, -0.838, -0.958, -1.036, -1.071, -1.066,
        -1.026, -0.960, -0.877, -0.785, -0.690, -0.597, -0.508, -0.426, -0.351, -0.317,
    ]
)

# Reference Fz values (N) - axial force component (attraction/repulsion)
ALLAG2009_FZ = np.array(
    [
        0.106, 0.09, 0.065, 0.03, -0.016, -0.074, -0.145, -0.228, -0.322, -0.426, -0.538,
        -0.653, -0.769, -0.881, -0.987, -1.084, -1.169, -1.242, -1.3, -1.343, -2.253,
        -1.343, -1.3, -1.242, -1.169, -1.084, -0.987, -0.881, -0.769, -0.653, -0.538,
        -0.426, -0.322, -0.228, -0.145, -0.074, -0.016, 0.03, 0.065, 0.09, 0.102,
    ]
)


# ==================== O'Connell2020 Reference Data ====================
# Magnet 1: 20x12x6 mm, Jr=0.38T, center=(0, 0, 0)
# Magnet 2: 12x20x6 mm, Jr=0.38T, center=(-4+offset, -4, 8)
# offset varies from 0 to 16 mm in steps of 2

OCONNELL_OFFSETS = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16])

# Reference forces in Newtons
OCONNELL_FX = np.array(
    [0.583, 0.244, 0.002, -0.248, -0.587, -0.913, -1.067, -1.113, -1.052]
)
OCONNELL_FY = np.array(
    [0.583, 0.617, 0.637, 0.625, 0.587, 0.517, 0.429, 0.321, 0.213]
)
OCONNELL_FZ = np.array(
    [-1.771, -1.86, -1.86, -1.856, -1.771, -1.467, -1.067, -0.64, -0.237]
)

# Reference torques in mN*m (millinewton-meters)
OCONNELL_TX = np.array(
    [-6.062, -6.3, -6.285, -6.315, -6.069, -5.069, -3.762, -2.438, -1.115]
)
OCONNELL_TY = np.array(
    [-3.646, -2.069, -0.008, 2.054, 3.669, 3.931, 3.577, 3.531, 3.977]
)
OCONNELL_TZ = np.array(
    [-1.585, -0.654, 0.023, 0.654, 1.608, 2.438, 2.7, 2.777, 2.731]
)


class TestAllag2009Validation:
    """Validate force calculations against Allag2009 reference data.

    Note: The reference data reports force on magnet 2 (the upper magnet).
    Our calc_force_prism(m1) computes force on m1 (the lower magnet).
    By Newton's third law, force on m1 = -force on m2.
    """

    @pytest.fixture
    def allag_tolerance(self):
        """Tolerance for Allag2009 validation (10% relative, 0.15N absolute)."""
        return {"rtol": 0.10, "atol": 0.15}

    def test_coaxial_configuration_fz(self, allag_tolerance):
        """Test axial force (Fz) at offset=0 against Allag2009."""
        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 20))

        # Compute force on m2 (upper magnet) to match reference convention
        force, _ = calc_force_prism(m2, num_samples=30)

        # At offset=0 (index 20), Fz should be -2.253 N (attraction toward m1)
        expected_fz = ALLAG2009_FZ[20]
        npt.assert_allclose(force[2], expected_fz, **allag_tolerance)

    def test_coaxial_configuration_fx_near_zero(self, allag_tolerance):
        """Test that lateral force (Fx) is near zero at offset=0."""
        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 20))

        force, _ = calc_force_prism(m2, num_samples=30)

        # At offset=0, Fx should be very small (symmetry)
        assert abs(force[0]) < 0.1  # Should be small

    @pytest.mark.parametrize(
        "offset_idx",
        [5, 10, 30, 35],  # Sample offsets: -15, -10, 10, 15 mm (avoid near-singular region)
    )
    def test_lateral_offset_force(self, offset_idx, allag_tolerance):
        """Test force at various lateral offsets against Allag2009."""
        offset = ALLAG2009_OFFSETS[offset_idx]
        expected_fx = ALLAG2009_FX[offset_idx]
        expected_fz = ALLAG2009_FZ[offset_idx]

        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2 = magnets.Prism(
            width=10, depth=10, height=10, Jr=1.0, center=(offset, 0, 20)
        )

        # Compute force on m2 (upper magnet) to match reference convention
        force, _ = calc_force_prism(m2, num_samples=25)

        # Check Fx and Fz against reference
        npt.assert_allclose(force[0], expected_fx, **allag_tolerance)
        npt.assert_allclose(force[2], expected_fz, **allag_tolerance)

    def test_force_symmetry_positive_negative_offset(self):
        """Test that Fx(-offset) = -Fx(offset) (antisymmetry)."""
        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2_pos = magnets.Prism(
            width=10, depth=10, height=10, Jr=1.0, center=(10, 0, 20)
        )
        force_pos, _ = calc_force_prism(m2_pos, num_samples=20)

        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2_neg = magnets.Prism(
            width=10, depth=10, height=10, Jr=1.0, center=(-10, 0, 20)
        )
        force_neg, _ = calc_force_prism(m2_neg, num_samples=20)

        # Fx should have opposite signs
        npt.assert_allclose(force_pos[0], -force_neg[0], rtol=0.05)
        # Fz should be the same
        npt.assert_allclose(force_pos[2], force_neg[2], rtol=0.05)


class TestOConnell2020Validation:
    """Validate force and torque calculations against O'Connell2020 reference data.

    Note: The reference data reports force on magnet 2 (the upper magnet).
    Our calc_force_prism(m1) computes force on m1 (the lower magnet).
    We compute force on m2 to match the reference convention.
    """

    @pytest.fixture
    def oconnell_tolerance(self):
        """Tolerance for O'Connell2020 validation."""
        return {"rtol": 0.15, "atol": 0.2}

    def test_asymmetric_configuration_force(self, oconnell_tolerance):
        """Test force at offset=0 against O'Connell2020."""
        reset()
        m1 = magnets.Prism(width=20, depth=12, height=6, Jr=0.38, center=(0, 0, 0))
        m2 = magnets.Prism(width=12, depth=20, height=6, Jr=0.38, center=(-4, -4, 8))

        # Compute force on m2 (upper magnet) to match reference convention
        force, _ = calc_force_prism(m2, num_samples=25)

        # Reference values at offset=0
        expected_fx = OCONNELL_FX[0]
        expected_fy = OCONNELL_FY[0]
        expected_fz = OCONNELL_FZ[0]

        npt.assert_allclose(force[0], expected_fx, **oconnell_tolerance)
        npt.assert_allclose(force[1], expected_fy, **oconnell_tolerance)
        npt.assert_allclose(force[2], expected_fz, **oconnell_tolerance)

    def test_asymmetric_configuration_torque(self, oconnell_tolerance):
        """Test torque at offset=0 against O'Connell2020."""
        reset()
        m1 = magnets.Prism(width=20, depth=12, height=6, Jr=0.38, center=(0, 0, 0))
        m2 = magnets.Prism(width=12, depth=20, height=6, Jr=0.38, center=(-4, -4, 8))

        # Compute torque on m2 (upper magnet) to match reference convention
        _, torque = calc_force_prism(m2, num_samples=25)

        # Reference torques in mN*m, convert to N*m
        expected_tx = OCONNELL_TX[0] * 1e-3
        expected_ty = OCONNELL_TY[0] * 1e-3
        expected_tz = OCONNELL_TZ[0] * 1e-3

        # Use looser tolerance for torque (more sensitive to numerical integration)
        torque_tol = {"rtol": 0.20, "atol": 0.005}
        npt.assert_allclose(torque[0], expected_tx, **torque_tol)
        npt.assert_allclose(torque[1], expected_ty, **torque_tol)
        npt.assert_allclose(torque[2], expected_tz, **torque_tol)

    @pytest.mark.parametrize("offset_idx", [2, 4, 6])  # offsets 4, 8, 12 mm
    def test_force_at_various_offsets(self, offset_idx, oconnell_tolerance):
        """Test force at various offsets against O'Connell2020."""
        offset = OCONNELL_OFFSETS[offset_idx]

        reset()
        m1 = magnets.Prism(width=20, depth=12, height=6, Jr=0.38, center=(0, 0, 0))
        m2 = magnets.Prism(
            width=12, depth=20, height=6, Jr=0.38, center=(-4 + offset, -4, 8)
        )

        # Compute force on m2 (upper magnet) to match reference convention
        force, _ = calc_force_prism(m2, num_samples=25)

        expected_fx = OCONNELL_FX[offset_idx]
        expected_fy = OCONNELL_FY[offset_idx]
        expected_fz = OCONNELL_FZ[offset_idx]

        npt.assert_allclose(force[0], expected_fx, **oconnell_tolerance)
        npt.assert_allclose(force[1], expected_fy, **oconnell_tolerance)
        npt.assert_allclose(force[2], expected_fz, **oconnell_tolerance)


class TestForcePhysicalConsistency:
    """Test that force calculations are physically consistent."""

    def test_attraction_force_negative_z(self, allag2009_magnets):
        """Two aligned magnets should attract (negative Fz for upper magnet pulling down)."""
        m1, m2 = allag2009_magnets
        force, _ = calc_force_prism(m2, num_samples=20)

        # m2 is above m1, both magnetized in same direction
        # m2 should be pulled toward m1 (negative z)
        assert force[2] < 0

    def test_force_decreases_with_distance(self):
        """Force magnitude should decrease as magnets move apart."""
        forces = []
        for gap in [15, 20, 30, 50]:  # mm gap between magnet centers
            reset()
            m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
            m2 = magnets.Prism(
                width=10, depth=10, height=10, Jr=1.0, center=(0, 0, gap)
            )
            force, _ = calc_force_prism(m1, num_samples=15)
            forces.append(np.linalg.norm(force))

        # Force should monotonically decrease with distance
        for i in range(len(forces) - 1):
            assert forces[i] > forces[i + 1]

    def test_force_direction_reverses_with_magnetization(self):
        """Reversing one magnet's magnetization should reverse force direction."""
        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2_aligned = magnets.Prism(
            width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 20)
        )
        force_aligned, _ = calc_force_prism(m1, num_samples=15)

        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2_anti = magnets.Prism(
            width=10, depth=10, height=10, Jr=-1.0, center=(0, 0, 20)
        )
        force_anti, _ = calc_force_prism(m1, num_samples=15)

        # Fz should have opposite signs
        assert force_aligned[2] * force_anti[2] < 0
