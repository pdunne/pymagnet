# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""TOML-based configuration module for pymagnet.

Provides a declarative interface to define magnet configurations,
run field calculations, and generate plots from TOML files.

Usage:
    from pymagnet.config import run, load, validate

    # Run a full simulation from a TOML file
    result = run("my_config.toml")
    # result["results"] is a list of dicts, one per [[grid]]/[[plot]] pair
    # Each contains: points, field, figure (if plot enabled)

    # Load and validate only (no calculation)
    config = load("my_config.toml")

    # Validate and get error list
    errors = validate("my_config.toml")
"""

from pathlib import Path

from ._loader import load
from ._runner import run
from ._schema import SimulationConfig, validate_config


def validate(toml_path: str | Path) -> list[str]:
    """Validate a TOML configuration file without running it.

    Args:
        toml_path: Path to the TOML file.

    Returns:
        List of error messages. Empty list means valid.
    """
    config = load(toml_path)
    return validate_config(config)


__all__ = ["load", "run", "validate", "SimulationConfig"]
