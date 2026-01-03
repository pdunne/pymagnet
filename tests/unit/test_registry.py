# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Tests for Registry class and magnet instance tracking."""

import pytest

from pymagnet import magnets, reset
from pymagnet.magnets._magnet_base import Magnet, Registry


class TestRegistry:
    """Tests for Registry base class."""

    def test_instances_is_weakset(self):
        """Instances should be stored in WeakSet."""
        from weakref import WeakSet

        assert isinstance(Registry.instances, WeakSet)

    def test_new_registers_instance(self):
        """Creating a magnet registers it."""
        initial_count = Magnet.get_num_instances()
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        assert Magnet.get_num_instances() == initial_count + 1

    def test_get_instances_returns_set(self):
        """get_instances returns set containing instances."""
        m = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        instances = Magnet.get_instances()
        assert m in instances

    def test_get_num_instances(self):
        """get_num_instances returns count."""
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        _ = magnets.Square(width=10.0, Jr=1.0)
        count = Magnet.get_num_instances()
        assert count == 2

    def test_get_num_instances_print(self, capsys):
        """get_num_instances with Print_Val=True prints count."""
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        magnets.Rectangle.get_num_instances(Print_Val=True)
        captured = capsys.readouterr()
        assert "1" in captured.out

    def test_reset_clears_all(self):
        """reset() clears all instances."""
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        _ = magnets.Prism(width=10.0, depth=10.0, height=10.0, Jr=1.0)
        reset()
        assert Magnet.get_num_instances() == 0

    def test_subclass_has_separate_instances(self):
        """Each subclass has its own instance set."""
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        _ = magnets.Prism(width=10.0, depth=10.0, height=10.0, Jr=1.0)

        rect_count = magnets.Rectangle.get_num_instances()
        prism_count = magnets.Prism.get_num_instances()

        assert rect_count == 1
        assert prism_count == 1

    def test_print_instances_no_magnets(self, capsys):
        """print_instances with no magnets prints 'No Instances'."""
        Magnet.print_instances()
        captured = capsys.readouterr()
        assert "No Instances" in captured.out

    def test_print_instances_with_magnets(self, capsys):
        """print_instances with magnets prints them."""
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        Magnet.print_instances()
        captured = capsys.readouterr()
        assert "Rectangle" in captured.out


class TestMagnetBase:
    """Tests for Magnet base class."""

    def test_mag_type(self):
        """Magnet has mag_type attribute."""
        assert Magnet.mag_type == "Magnet"

    def test_tolerance_attribute(self):
        """Magnet has tolerance attribute."""
        assert hasattr(Magnet, "tol")
        assert Magnet.tol > 0


class TestResetFunction:
    """Tests for module-level reset function."""

    def test_reset_clears_2d_magnets(self):
        """reset clears 2D magnets."""
        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        _ = magnets.Circle(radius=5.0, Jr=1.0)
        reset()
        assert magnets.Rectangle.get_num_instances() == 0
        assert magnets.Circle.get_num_instances() == 0

    def test_reset_clears_3d_magnets(self):
        """reset clears 3D magnets."""
        _ = magnets.Prism(width=10.0, depth=10.0, height=10.0, Jr=1.0)
        _ = magnets.Cylinder(radius=5.0, length=10.0, Jr=1.0)
        reset()
        assert magnets.Prism.get_num_instances() == 0
        assert magnets.Cylinder.get_num_instances() == 0

    def test_reset_is_idempotent(self):
        """Calling reset multiple times is safe."""
        reset()
        reset()
        reset()
        assert Magnet.get_num_instances() == 0


class TestListFunction:
    """Tests for module-level list function."""

    def test_list_prints_magnets(self, capsys):
        """list() prints all magnets via pymagnet.list()."""
        import pymagnet

        _ = magnets.Rectangle(width=10.0, height=20.0, Jr=1.0)
        pymagnet.list()  # Use pymagnet.list(), not magnets.list()
        # list() calls print_instances which should output the magnet
