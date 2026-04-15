# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2026 Peter Dunne
"""Tests for elliptic integral implementations (Bulirsch cel and Carlson)."""

import math

import numpy as np
import pytest

from pymagnet.utils._elliptic import _elliprd, _elliprf, _elliprj, cel, cel_carlson


# ---------------------------------------------------------------------------
# Tests for Carlson primitives
# ---------------------------------------------------------------------------
class TestCarlsonRF:
    def test_rf_equal_args(self):
        """RF(a, a, a) = 1/sqrt(a)."""
        assert _elliprf(1.0, 1.0, 1.0) == pytest.approx(1.0, abs=1e-10)
        assert _elliprf(4.0, 4.0, 4.0) == pytest.approx(0.5, abs=1e-10)

    def test_rf_symmetric(self):
        """RF is symmetric in all three arguments."""
        val1 = _elliprf(0.5, 1.0, 2.0)
        val2 = _elliprf(2.0, 0.5, 1.0)
        val3 = _elliprf(1.0, 2.0, 0.5)
        assert val1 == pytest.approx(val2, abs=1e-10)
        assert val1 == pytest.approx(val3, abs=1e-10)

    def test_rf_vs_scipy(self):
        sp = pytest.importorskip("scipy")
        from scipy.special import elliprf as sp_rf

        cases = [(0.0, 1.0, 1.0), (0.5, 1.0, 2.0), (1.0, 2.0, 3.0)]
        for x, y, z in cases:
            assert _elliprf(x, y, z) == pytest.approx(sp_rf(x, y, z), abs=1e-10)


class TestCarlsonRD:
    def test_rd_vs_scipy(self):
        sp = pytest.importorskip("scipy")
        from scipy.special import elliprd as sp_rd

        cases = [(0.0, 2.0, 1.0), (1.0, 2.0, 3.0), (0.5, 1.0, 1.5)]
        for x, y, z in cases:
            assert _elliprd(x, y, z) == pytest.approx(sp_rd(x, y, z), abs=1e-10)


class TestCarlsonRJ:
    def test_rj_reduces_to_rd(self):
        """RJ(x, y, z, z) = RD(x, y, z)."""
        cases = [(0.0, 2.0, 1.0), (1.0, 2.0, 3.0)]
        for x, y, z in cases:
            assert _elliprj(x, y, z, z) == pytest.approx(
                _elliprd(x, y, z), rel=1e-6
            )

    def test_rj_vs_scipy(self):
        sp = pytest.importorskip("scipy")
        from scipy.special import elliprj as sp_rj

        cases = [
            (0.0, 1.0, 1.0, 2.0),
            (1.0, 2.0, 3.0, 0.5),
            (0.5, 1.0, 1.5, 2.0),
        ]
        for x, y, z, p in cases:
            assert _elliprj(x, y, z, p) == pytest.approx(
                sp_rj(x, y, z, p), rel=1e-6
            )


# ---------------------------------------------------------------------------
# Tests for cel vs standard elliptic integrals
# ---------------------------------------------------------------------------
class TestCelIdentities:
    def test_cel_1_1_1_1_is_pi_half(self):
        """cel(1, 1, 1, 1) = K(0) = π/2."""
        assert float(cel(1.0, 1.0, 1.0, 1.0)) == pytest.approx(
            math.pi / 2, abs=1e-10
        )

    def test_cel_vs_ellipk(self):
        """cel(kc, 1, 1, 1) = K(k) where k² = 1 - kc²."""
        sp = pytest.importorskip("scipy")
        from scipy.special import ellipk

        for kc in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
            k_sq = 1.0 - kc * kc
            expected = ellipk(k_sq)
            result = float(cel(kc, 1.0, 1.0, 1.0))
            assert result == pytest.approx(expected, rel=1e-8), (
                f"kc={kc}: cel={result}, K={expected}"
            )

    def test_cel_vs_ellipe(self):
        """cel(kc, 1, 1, kc²) = E(k)."""
        sp = pytest.importorskip("scipy")
        from scipy.special import ellipe

        for kc in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
            k_sq = 1.0 - kc * kc
            kc2 = kc * kc
            expected = ellipe(k_sq)
            result = float(cel(kc, 1.0, 1.0, kc2))
            assert result == pytest.approx(expected, rel=1e-8), (
                f"kc={kc}: cel={result}, E={expected}"
            )


# ---------------------------------------------------------------------------
# Three-way comparison: scipy vs Bulirsch vs Carlson
# ---------------------------------------------------------------------------
class TestThreeWayComparison:
    """Compare scipy, Bulirsch cel, and Carlson cel_carlson against each other."""

    KCS = np.array([0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99])

    def test_K_three_way(self):
        """K(k): scipy ellipk vs Bulirsch vs Carlson."""
        pytest.importorskip("scipy")
        from scipy.special import ellipk

        p = np.ones_like(self.KCS)
        c = np.ones_like(self.KCS)
        s = np.ones_like(self.KCS)
        k_sq = 1.0 - self.KCS**2

        ref = np.array([ellipk(k) for k in k_sq])
        res_bulirsch = cel(self.KCS, p, c, s)
        res_carlson = cel_carlson(self.KCS, p, c, s)

        np.testing.assert_allclose(res_bulirsch, ref, rtol=1e-10, err_msg="Bulirsch vs scipy K")
        np.testing.assert_allclose(res_carlson, ref, rtol=1e-8, err_msg="Carlson vs scipy K")

    def test_E_three_way(self):
        """E(k): scipy ellipe vs Bulirsch vs Carlson."""
        pytest.importorskip("scipy")
        from scipy.special import ellipe

        p = np.ones_like(self.KCS)
        c = np.ones_like(self.KCS)
        s = self.KCS**2
        k_sq = 1.0 - self.KCS**2

        ref = np.array([ellipe(k) for k in k_sq])
        res_bulirsch = cel(self.KCS, p, c, s)
        res_carlson = cel_carlson(self.KCS, p, c, s)

        np.testing.assert_allclose(res_bulirsch, ref, rtol=1e-10, err_msg="Bulirsch vs scipy E")
        np.testing.assert_allclose(res_carlson, ref, rtol=1e-8, err_msg="Carlson vs scipy E")


# ---------------------------------------------------------------------------
# Bulirsch vs Carlson direct comparison
# ---------------------------------------------------------------------------
class TestBulirschVsCarlson:
    """Compare Bulirsch cel against Carlson cel_carlson on representative inputs."""

    KCS = np.array([0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99])

    def _compare(self, kc_arr, p, c, s, atol=1e-10):
        p_arr = np.full_like(kc_arr, p)
        c_arr = np.full_like(kc_arr, c)
        s_arr = np.full_like(kc_arr, s)
        result_bulirsch = cel(kc_arr, p_arr, c_arr, s_arr)
        result_carlson = cel_carlson(kc_arr, p_arr, c_arr, s_arr)
        np.testing.assert_allclose(
            result_bulirsch,
            result_carlson,
            atol=atol,
            err_msg=f"Mismatch for p={p}, c={c}, s={s}",
        )

    def test_brho_case(self):
        """p=1, c=1, s=-1 — the Brho call pattern."""
        self._compare(self.KCS, p=1.0, c=1.0, s=-1.0)

    def test_bz_case(self):
        """p=γ², c=1, s=γ — the Bz call pattern (γ=0.5)."""
        self._compare(self.KCS, p=0.25, c=1.0, s=0.5)

    def test_complete_K(self):
        """p=1, c=1, s=1 — K(k)."""
        self._compare(self.KCS, p=1.0, c=1.0, s=1.0)

    def test_general_case(self):
        """p=0.5, c=2, s=3 — general parameters."""
        self._compare(self.KCS, p=0.5, c=2.0, s=3.0)


# ---------------------------------------------------------------------------
# Array broadcasting
# ---------------------------------------------------------------------------
class TestArrayBroadcasting:
    def test_array_input(self):
        """cel accepts ndarray arguments and returns element-wise results."""
        kc = np.array([0.1, 0.5, 0.9])
        p = np.ones(3)
        c = np.ones(3)
        s = -np.ones(3)
        result = cel(kc, p, c, s)
        assert result.shape == (3,)
        for i in range(3):
            scalar = float(cel(kc[i], p[i], c[i], s[i]))
            assert result[i] == pytest.approx(scalar, abs=1e-12)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------
class TestEdgeCases:
    def test_kc_zero_returns_nan(self):
        result = float(cel(0.0, 1.0, 1.0, 1.0))
        assert math.isnan(result)

    def test_kc_near_one(self):
        """kc ≈ 1 should give cel ≈ π/2 for (kc, 1, 1, 1)."""
        result = float(cel(0.9999, 1.0, 1.0, 1.0))
        assert result == pytest.approx(math.pi / 2, rel=1e-3)
