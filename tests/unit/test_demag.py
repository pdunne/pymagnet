# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for the demagnetization solver."""

import numpy as np
import numpy.testing as npt
import pytest

from pymagnet.utils._demag import (
    DemagResult,
    build_MH_interpolator,
    solve_demagnetization,
    solve_demag_tanh,
    solve_demag_tanh_batch,
    tanh_MH_model,
)


class TestBuildMHInterpolator:
    """Tests for build_MH_interpolator validation."""

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="same shape"):
            build_MH_interpolator(np.array([1, 2, 3]), np.array([1, 2]))

    def test_non_1d_raises(self):
        a = np.array([[1, 2], [3, 4]])
        with pytest.raises(ValueError, match="1D"):
            build_MH_interpolator(a, a)

    def test_non_increasing_H_raises(self):
        with pytest.raises(ValueError, match="strictly monotonically increasing"):
            build_MH_interpolator(np.array([1, 3, 2]), np.array([1, 2, 3]))

    def test_decreasing_M_raises(self):
        with pytest.raises(ValueError, match="non-decreasing"):
            build_MH_interpolator(np.array([1, 2, 3]), np.array([1, 3, 2]))

    def test_valid_input_returns_callable(self):
        H = np.linspace(0, 100, 50)
        M = 5.0 * H
        interp = build_MH_interpolator(H, M)
        assert callable(interp)


class TestTanhMHModel:
    """Tests for tanh_MH_model."""

    def test_low_field_linear(self):
        """At low fields M ≈ chi * H."""
        Ms = 1e6
        chi = 10.0
        f = tanh_MH_model(Ms, chi)
        H_low = 1.0  # very small relative to Ms/chi
        npt.assert_allclose(f(H_low), chi * H_low, rtol=1e-6)

    def test_high_field_saturation(self):
        """At very high fields M → Ms."""
        Ms = 1e6
        chi = 10.0
        f = tanh_MH_model(Ms, chi)
        npt.assert_allclose(f(1e12), Ms, rtol=1e-6)

    def test_accepts_array(self):
        f = tanh_MH_model(1e6, 10.0)
        result = f(np.array([0.0, 1e3, 1e6]))
        assert result.shape == (3,)


class TestSolveDemagnetization:
    """Tests for the scipy-based solver."""

    @pytest.fixture()
    def linear_interp(self):
        """Linear M(H) = chi * H with chi=10."""
        chi = 10.0
        H = np.linspace(0, 1e6, 500)
        M = chi * H
        return build_MH_interpolator(H, M), chi, M[-1]

    def test_linear_material(self, linear_interp):
        """Solver matches analytical solution for linear material."""
        interp, chi, M_sat = linear_interp
        H_ext = 1e5
        N = 0.5
        result = solve_demagnetization(H_ext, interp, M_sat, N=N)
        M_expected = chi * H_ext / (1 + N * chi)
        npt.assert_allclose(result.M_solution, M_expected, rtol=1e-4)
        assert result.converged

    def test_saturated_input(self):
        """Near-saturation with nonlinear tanh curve."""
        M_s = 1e6
        H = np.linspace(0, 1e6, 500)
        M = M_s * np.tanh(H / 1e4)
        interp = build_MH_interpolator(H, M)
        result = solve_demagnetization(1e6, interp, M_s, N=0.5)
        npt.assert_allclose(result.M_solution, M_s, rtol=0.01)

    def test_zero_field(self, linear_interp):
        """Zero applied field gives zero magnetization."""
        interp, _, M_sat = linear_interp
        result = solve_demagnetization(0.0, interp, M_sat, N=0.5)
        npt.assert_allclose(result.M_solution, 0.0, atol=1e-6)

    @pytest.mark.parametrize(
        "N_val,label",
        [(1 / 3, "sphere"), (1.0, "thin_film")],
    )
    def test_demag_factor_generalisation(self, linear_interp, N_val, label):
        """Solver works for different demagnetizing factors."""
        interp, chi, M_sat = linear_interp
        H_ext = 1e5
        result = solve_demagnetization(H_ext, interp, M_sat, N=N_val)
        M_expected = chi * H_ext / (1 + N_val * chi)
        H_int_expected = H_ext - N_val * M_expected
        npt.assert_allclose(result.M_solution, M_expected, rtol=1e-4)
        npt.assert_allclose(result.H_int, H_int_expected, rtol=1e-4)

    def test_tanh_model(self):
        """Solver works with tanh_MH_model callable."""
        Ms = 1e6
        chi = 10.0
        f_mh = tanh_MH_model(Ms, chi)
        result = solve_demagnetization(1e5, f_mh, Ms, N=0.5)
        M_approx = chi * 1e5 / (1 + 0.5 * chi)
        assert result.converged
        npt.assert_allclose(result.M_solution, M_approx, rtol=0.05)

    def test_returns_demag_result(self, linear_interp):
        interp, _, M_sat = linear_interp
        result = solve_demagnetization(1e5, interp, M_sat)
        assert isinstance(result, DemagResult)


class TestSolveDemagtanhNumba:
    """Tests for the numba-accelerated scalar solver."""

    def test_matches_scipy(self):
        """Numba scalar solver matches scipy result."""
        Ms, chi, H_ext, N = 1e6, 10.0, 1e5, 0.5
        f_mh = tanh_MH_model(Ms, chi)
        result_scipy = solve_demagnetization(H_ext, f_mh, Ms, N=N)
        M_nb, _, conv_nb = solve_demag_tanh(H_ext, Ms, chi, N=N)
        assert conv_nb
        npt.assert_allclose(M_nb, result_scipy.M_solution, rtol=1e-6)

    def test_zero_field(self):
        M, H_int, conv = solve_demag_tanh(0.0, 1e6, 10.0, N=0.5)
        assert conv
        npt.assert_allclose(M, 0.0, atol=1e-6)


class TestSolveDemagtanhBatch:
    """Tests for the numba-accelerated batch solver."""

    def test_monotonic_output(self):
        """M increases monotonically with increasing H_ext."""
        H = np.linspace(0, 5e5, 1000)
        M, _ = solve_demag_tanh_batch(H, 1e6, 10.0, N=0.5)
        assert np.all(np.diff(M) >= 0)

    def test_zero_at_zero_field(self):
        H = np.linspace(0, 5e5, 100)
        M, _ = solve_demag_tanh_batch(H, 1e6, 10.0, N=0.5)
        npt.assert_allclose(M[0], 0.0, atol=1e-6)

    def test_consistent_with_scalar(self):
        """Batch results match scalar solver at sampled points."""
        Ms, chi, N = 1e6, 10.0, 0.5
        H = np.linspace(0, 5e5, 1000)
        M_batch, _ = solve_demag_tanh_batch(H, Ms, chi, N=N)
        for idx in [0, 250, 500, 750, 999]:
            M_ref, _, _ = solve_demag_tanh(H[idx], Ms, chi, N=N)
            npt.assert_allclose(M_batch[idx], M_ref, atol=1e-4)

    def test_output_shapes(self):
        H = np.linspace(0, 1e5, 50)
        M, H_int = solve_demag_tanh_batch(H, 1e6, 10.0)
        assert M.shape == (50,)
        assert H_int.shape == (50,)
