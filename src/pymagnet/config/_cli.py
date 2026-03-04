# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Copyright 2021 Peter Dunne
"""Command-line interface for pymagnet TOML configuration runner."""

import argparse
import sys


def main(argv: list[str] | None = None) -> None:
    """Entry point for the pymagnet CLI."""
    parser = argparse.ArgumentParser(
        prog="pymagnet",
        description="Run pymagnet simulations from TOML configuration files.",
    )
    parser.add_argument(
        "config_file",
        help="Path to TOML configuration file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output directory for saved figures (default: current directory)",
    )

    args = parser.parse_args(argv)

    from ._runner import run

    try:
        result = run(args.config_file, output_dir=args.output)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    config = result["config"]
    print(f"Simulation: {config.name or config.dimension}")
    print(f"  Magnets: {len(result['magnets'])}")
    print(f"  Plots: {len(result['results'])}")

    for i, entry in enumerate(result["results"]):
        has_fig = entry.get("figure") is not None
        saved = entry.get("saved_to")
        status = "generated"
        if saved:
            status = f"saved to {saved}"
        print(f"    [{i}] figure {'yes' if has_fig else 'no'} — {status}")

    if "force" in result:
        f = result["force"]
        print(f"  Force:  {f['force']}")
        print(f"  Torque: {f['torque']}")


if __name__ == "__main__":
    main()
