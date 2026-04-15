# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Smoke tests for 2D gradient/force plot functions."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

import pymagnet as pm
import pymagnet.plots as mplt


@pytest.fixture(autouse=True)
def _reset_magnets():
    """Reset magnet registry before each test."""
    pm.reset()
    yield
    pm.reset()
    plt.close("all")


@pytest.fixture()
def field_and_points():
    """Create a simple magnet, grid, and field for testing."""
    _ = pm.magnets.Rectangle(width=5.0, height=10.0, Jr=1.0)
    points = pm.grid2D(20.0, 20.0, num_points=20)
    field = pm.get_field_2D(points)
    return points, field


class TestPlot2DContourGradient:
    """Tests for plot_2D_contour_gradient."""

    def test_returns_fig_ax(self, field_and_points):
        points, field = field_and_points
        fig, ax = mplt.plot_2D_contour_gradient(points, field, show_magnets=False)
        assert fig is not None
        assert ax is not None

    def test_custom_kwargs(self, field_and_points):
        points, field = field_and_points
        fig, ax = mplt.plot_2D_contour_gradient(
            points, field, cmin=0.0, cmax=1.0, cmap="plasma", show_magnets=False
        )
        assert fig is not None


class TestPlot2DContourForce:
    """Tests for plot_2D_contour_force."""

    def test_returns_fig_ax(self, field_and_points):
        points, field = field_and_points
        fig, ax = mplt.plot_2D_contour_force(
            points, field, chi_m=1e-3, c=1.0, show_magnets=False
        )
        assert fig is not None
        assert ax is not None


class TestPlot2DContourBdotgradB:
    """Tests for plot_2D_contour_BdotgradB."""

    def test_returns_fig_ax(self, field_and_points):
        points, field = field_and_points
        fig, ax = mplt.plot_2D_contour_BdotgradB(points, field, show_magnets=False)
        assert fig is not None
        assert ax is not None


class TestPlot2DContourGradB2:
    """Tests for plot_2D_contour_gradB2."""

    def test_returns_fig_ax(self, field_and_points):
        points, field = field_and_points
        fig, ax = mplt.plot_2D_contour_gradB2(points, field, show_magnets=False)
        assert fig is not None
        assert ax is not None

    def test_has_vector_arrows_by_default(self, field_and_points):
        """Default num_arrows=10 produces quiver overlay."""
        points, field = field_and_points
        fig, ax = mplt.plot_2D_contour_gradB2(points, field, show_magnets=False)
        from matplotlib.quiver import Quiver

        quivers = [c for c in ax.get_children() if isinstance(c, Quiver)]
        assert len(quivers) > 0


class TestPlot2DContourGradB:
    """Tests for plot_2D_contour_gradB."""

    def test_returns_fig_ax(self, field_and_points):
        points, field = field_and_points
        fig, ax = mplt.plot_2D_contour_gradB(points, field, show_magnets=False)
        assert fig is not None
        assert ax is not None

    def test_has_vector_arrows_by_default(self, field_and_points):
        """Default num_arrows=10 produces quiver overlay."""
        points, field = field_and_points
        fig, ax = mplt.plot_2D_contour_gradB(points, field, show_magnets=False)
        # Quiver collection should be present on the axes
        from matplotlib.quiver import Quiver

        quivers = [c for c in ax.get_children() if isinstance(c, Quiver)]
        assert len(quivers) > 0


class TestPlot2DContourJacobian:
    """Tests for plot_2D_contour_jacobian."""

    def test_returns_fig_and_2x2_axes(self, field_and_points):
        points, field = field_and_points
        fig, axes = mplt.plot_2D_contour_jacobian(points, field, show_magnets=False)
        assert fig is not None
        assert axes.shape == (2, 2)

    def test_custom_cmap(self, field_and_points):
        points, field = field_and_points
        fig, axes = mplt.plot_2D_contour_jacobian(
            points, field, cmap="seismic", show_magnets=False
        )
        assert fig is not None
