# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Unit tests for force/torque calculation functions.

Tests the force calculation functions in src/pymagnet/forces/:
- calc_force_prism (prism/cube magnets)
- calc_force_cylinder (cylindrical magnets)
- calc_force_sphere (spherical magnets)
- calc_force_mesh (STL/mesh magnets)
"""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet import magnets, reset
from pymagnet.forces._cylinder_force import calc_force_cylinder
from pymagnet.forces._mesh_force import calc_force_mesh
from pymagnet.forces._prism_force import calc_force_prism
from pymagnet.forces._sphere_force import calc_force_sphere

# ==================== Basic Return Type Tests ====================


class TestCalcForcePrism:
    """Test basic behavior of calc_force_prism."""

    def test_returns_tuple(self, prism_pair_aligned):
        """calc_force_prism returns a tuple of (force, torque)."""
        m1, m2 = prism_pair_aligned
        result = calc_force_prism(m1, num_samples=5)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_force_shape_is_3d(self, prism_pair_aligned):
        """Force array has shape (3,)."""
        m1, m2 = prism_pair_aligned
        force, torque = calc_force_prism(m1, num_samples=5)
        assert force.shape == (3,)

    def test_torque_shape_is_3d(self, prism_pair_aligned):
        """Torque array has shape (3,)."""
        m1, m2 = prism_pair_aligned
        force, torque = calc_force_prism(m1, num_samples=5)
        assert torque.shape == (3,)

    def test_force_values_are_finite(self, prism_pair_aligned):
        """Force values should be finite (no NaN or inf)."""
        m1, m2 = prism_pair_aligned
        force, torque = calc_force_prism(m1, num_samples=5)
        assert np.all(np.isfinite(force))

    def test_torque_values_are_finite(self, prism_pair_aligned):
        """Torque values should be finite (no NaN or inf)."""
        m1, m2 = prism_pair_aligned
        force, torque = calc_force_prism(m1, num_samples=5)
        assert np.all(np.isfinite(torque))


class TestCalcForceCylinder:
    """Test basic behavior of calc_force_cylinder."""

    def test_returns_tuple(self, cylinder_pair_aligned):
        """calc_force_cylinder returns a tuple of (force, torque)."""
        m1, m2 = cylinder_pair_aligned
        result = calc_force_cylinder(m1, num_segments=5)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_force_shape_is_3d(self, cylinder_pair_aligned):
        """Force array has shape (3,)."""
        m1, m2 = cylinder_pair_aligned
        force, torque = calc_force_cylinder(m1, num_segments=5)
        assert force.shape == (3,)

    def test_torque_shape_is_3d(self, cylinder_pair_aligned):
        """Torque array has shape (3,)."""
        m1, m2 = cylinder_pair_aligned
        force, torque = calc_force_cylinder(m1, num_segments=5)
        assert torque.shape == (3,)

    def test_force_values_are_finite(self, cylinder_pair_aligned):
        """Force values should be finite."""
        m1, m2 = cylinder_pair_aligned
        force, torque = calc_force_cylinder(m1, num_segments=5)
        assert np.all(np.isfinite(force))

    def test_torque_values_are_finite(self, cylinder_pair_aligned):
        """Torque values should be finite."""
        m1, m2 = cylinder_pair_aligned
        force, torque = calc_force_cylinder(m1, num_segments=5)
        assert np.all(np.isfinite(torque))


class TestCalcForceSphere:
    """Test basic behavior of calc_force_sphere."""

    def test_returns_tuple(self, sphere_pair_aligned):
        """calc_force_sphere returns a tuple of (force, torque)."""
        m1, m2 = sphere_pair_aligned
        result = calc_force_sphere(m1, num_samples=10)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_force_shape_is_3d(self, sphere_pair_aligned):
        """Force array has shape (3,)."""
        m1, m2 = sphere_pair_aligned
        force, torque = calc_force_sphere(m1, num_samples=10)
        assert force.shape == (3,)

    def test_torque_shape_is_3d(self, sphere_pair_aligned):
        """Torque array has shape (3,)."""
        m1, m2 = sphere_pair_aligned
        force, torque = calc_force_sphere(m1, num_samples=10)
        assert torque.shape == (3,)

    def test_force_values_are_finite(self, sphere_pair_aligned):
        """Force values should be finite."""
        m1, m2 = sphere_pair_aligned
        force, torque = calc_force_sphere(m1, num_samples=10)
        assert np.all(np.isfinite(force))

    def test_torque_values_are_finite(self, sphere_pair_aligned):
        """Torque values should be finite."""
        m1, m2 = sphere_pair_aligned
        force, torque = calc_force_sphere(m1, num_samples=10)
        assert np.all(np.isfinite(torque))


# ==================== Newton's Third Law Tests ====================


class TestNewtonThirdLaw:
    """Test that F12 = -F21 (Newton's third law)."""

    def test_prism_pair_force_reciprocity(self, prism_pair_aligned, force_tolerance):
        """Force on prism 1 from prism 2 equals -force on prism 2 from prism 1."""
        m1, m2 = prism_pair_aligned

        force1, _ = calc_force_prism(m1, num_samples=10)
        force2, _ = calc_force_prism(m2, num_samples=10)

        # F12 = -F21
        npt.assert_allclose(force1, -force2, **force_tolerance)

    def test_cylinder_pair_force_reciprocity(
        self, cylinder_pair_aligned, force_tolerance
    ):
        """Force on cylinder 1 from cylinder 2 equals -force on cylinder 2 from cylinder 1."""
        m1, m2 = cylinder_pair_aligned

        force1, _ = calc_force_cylinder(m1, num_segments=10)
        force2, _ = calc_force_cylinder(m2, num_segments=10)

        npt.assert_allclose(force1, -force2, **force_tolerance)

    def test_sphere_pair_force_reciprocity(self, sphere_pair_aligned):
        """Force on sphere 1 from sphere 2 equals -force on sphere 2 from sphere 1."""
        m1, m2 = sphere_pair_aligned

        # Sphere force calculation requires more samples for accuracy
        force1, _ = calc_force_sphere(m1, num_samples=50)
        force2, _ = calc_force_sphere(m2, num_samples=50)

        # Use looser tolerance for sphere (numerical integration on curved surface)
        npt.assert_allclose(force1, -force2, rtol=0.15, atol=0.01)


# ==================== Force Direction Tests ====================


class TestForceDirection:
    """Test that force direction is correct based on magnetization alignment."""

    def test_aligned_magnets_attract_along_z(self, prism_pair_aligned):
        """Two magnets with same magnetization direction attract (negative Fz on upper)."""
        m1, m2 = prism_pair_aligned
        # m1 at z=0, m2 at z=25, both Jr=+1.0 (magnetized in +z)
        # m2 should be pulled toward m1 (negative z direction)
        force2, _ = calc_force_prism(m2, num_samples=10)
        assert force2[2] < 0  # Attracted toward m1

    def test_anti_aligned_magnets_repel(self, prism_pair_anti_aligned):
        """Two magnets with opposite magnetization repel."""
        m1, m2 = prism_pair_anti_aligned
        # m1 at z=0 with Jr=+1.0, m2 at z=25 with Jr=-1.0
        # Opposite poles face each other, should repel
        force2, _ = calc_force_prism(m2, num_samples=10)
        assert force2[2] > 0  # Repelled away from m1


# ==================== Symmetry Tests ====================


class TestSymmetryConditions:
    """Test symmetric configurations produce expected results."""

    def test_single_magnet_zero_force(self, prism_default):
        """A single magnet with no others should have zero force."""
        # Only one magnet in the registry
        force, torque = calc_force_prism(prism_default, num_samples=5)
        npt.assert_allclose(force, [0, 0, 0], atol=1e-10)
        npt.assert_allclose(torque, [0, 0, 0], atol=1e-10)

    def test_symmetric_pair_zero_lateral_force(
        self, prism_pair_aligned, force_tolerance
    ):
        """On-axis pair should have zero lateral force (Fx=Fy=0)."""
        m1, m2 = prism_pair_aligned
        # Both magnets on z-axis, symmetric about it
        force1, _ = calc_force_prism(m1, num_samples=10)

        # Lateral forces should be zero
        assert abs(force1[0]) < force_tolerance["atol"] * 10  # Fx ~ 0
        assert abs(force1[1]) < force_tolerance["atol"] * 10  # Fy ~ 0

    def test_aligned_magnetization_zero_torque(
        self, prism_pair_aligned, force_tolerance
    ):
        """Aligned magnets on z-axis should have minimal torque."""
        m1, m2 = prism_pair_aligned
        _, torque1 = calc_force_prism(m1, num_samples=10)

        # Torque should be very small due to symmetry
        assert np.linalg.norm(torque1) < 0.1  # Small torque


# ==================== Force Scaling Tests ====================


class TestForceScaling:
    """Test that force scales correctly with parameters."""

    @pytest.mark.parametrize("Jr", [0.5, 1.0, 2.0])
    def test_force_scales_with_Jr_squared(self, Jr):
        """Force should scale with Jr^2 (both magnets have same Jr)."""
        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 25))
        force_ref, _ = calc_force_prism(m1, num_samples=8)

        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=Jr, center=(0, 0, 0))
        m2 = magnets.Prism(width=10, depth=10, height=10, Jr=Jr, center=(0, 0, 25))
        force_Jr, _ = calc_force_prism(m1, num_samples=8)

        # Force scales as Jr^2 since both magnets scale
        expected_scale = Jr * Jr
        npt.assert_allclose(
            np.linalg.norm(force_Jr),
            np.linalg.norm(force_ref) * expected_scale,
            rtol=0.05,
        )


# ==================== Unit Scaling Tests ====================


class TestUnitScaling:
    """Test that unit conversion works correctly."""

    def test_force_mm_vs_m_gives_same_result(self):
        """Force calculation should give consistent results regardless of unit."""
        # Using mm
        reset()
        m1_mm = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2_mm = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 25))
        force_mm, torque_mm = calc_force_prism(m1_mm, num_samples=8, unit="mm")

        # Using m (convert dimensions)
        reset()
        m1_m = magnets.Prism(
            width=0.01, depth=0.01, height=0.01, Jr=1.0, center=(0, 0, 0)
        )
        m2_m = magnets.Prism(
            width=0.01, depth=0.01, height=0.01, Jr=1.0, center=(0, 0, 0.025)
        )
        force_m, torque_m = calc_force_prism(m1_m, num_samples=8, unit="m")

        # Forces should be equal (same physical configuration)
        # Use atol for near-zero components (x, y) and rtol for significant component (z)
        npt.assert_allclose(force_mm, force_m, rtol=0.05, atol=1e-10)


# ==================== Registry Interaction Tests ====================


class TestRegistryInteraction:
    """Test interaction with Magnet3D.instances registry."""

    def test_force_excludes_self(self, prism_pair_aligned):
        """Force calculation should not include self-interaction."""
        m1, m2 = prism_pair_aligned
        # Should not crash or give NaN from self-interaction
        force, torque = calc_force_prism(m1, num_samples=5)
        assert np.all(np.isfinite(force))
        assert np.all(np.isfinite(torque))

    def test_force_includes_all_magnets(self):
        """Force on magnet should include contributions from all other magnets."""
        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 25))

        force_2mag, _ = calc_force_prism(m1, num_samples=8)

        # Add a third magnet
        m3 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, -25))
        force_3mag, _ = calc_force_prism(m1, num_samples=8)

        # With symmetric third magnet, z-force should decrease (forces partially cancel)
        # The force magnitude should be different
        assert not np.allclose(force_2mag, force_3mag)

    def test_reset_clears_force_contributions(self):
        """After reset(), force on single magnet should be zero."""
        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 25))

        force_before, _ = calc_force_prism(m1, num_samples=5)
        assert np.linalg.norm(force_before) > 0  # Non-zero force

        reset()
        m1_new = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        force_after, _ = calc_force_prism(m1_new, num_samples=5)

        # After reset, only one magnet, so force should be zero
        npt.assert_allclose(force_after, [0, 0, 0], atol=1e-10)


# ==================== Numerical Stability Tests ====================


class TestNumericalStability:
    """Test numerical stability in various configurations."""

    def test_distant_magnets_small_force(self):
        """Distant magnets should have small force."""
        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2_close = magnets.Prism(
            width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 25)
        )
        force_close, _ = calc_force_prism(m1, num_samples=8)

        reset()
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2_far = magnets.Prism(
            width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 100)
        )
        force_far, _ = calc_force_prism(m1, num_samples=8)

        # Far magnet should produce smaller force
        assert np.linalg.norm(force_far) < np.linalg.norm(force_close)

    def test_close_magnets_finite_force(self):
        """Close (but not touching) magnets should give finite force."""
        reset()
        # Gap of 1mm between 10mm cubes
        m1 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 0))
        m2 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 11))

        force, torque = calc_force_prism(m1, num_samples=10)
        assert np.all(np.isfinite(force))
        assert np.all(np.isfinite(torque))

    def test_rotated_magnets_force_finite(self, prism_rotated):
        """Rotated magnets should give finite force."""
        reset()
        # Create a second magnet to interact with
        m2 = magnets.Prism(width=10, depth=10, height=10, Jr=1.0, center=(0, 0, 50))

        force, torque = calc_force_prism(prism_rotated, num_samples=8)
        assert np.all(np.isfinite(force))
        assert np.all(np.isfinite(torque))


# ==================== Convergence Tests ====================


@pytest.mark.slow
class TestConvergence:
    """Test that increasing sampling improves accuracy."""

    def test_prism_force_converges_with_num_samples(self, prism_pair_aligned):
        """Force should converge as num_samples increases."""
        m1, m2 = prism_pair_aligned

        forces = []
        for num_samples in [5, 10, 20, 40]:
            force, _ = calc_force_prism(m1, num_samples=num_samples)
            forces.append(force.copy())

        # Calculate differences between successive refinements
        diff1 = np.linalg.norm(forces[1] - forces[0])
        diff2 = np.linalg.norm(forces[2] - forces[1])
        diff3 = np.linalg.norm(forces[3] - forces[2])

        # Differences should decrease (convergence)
        assert diff2 < diff1
        assert diff3 < diff2
