"""Command-line interface: `fimtyping <genome.fna>`."""
import sys

from pgcore.report import run_cli

from . import analyze


def main(argv=None):
    return run_cli(
        "fimtyping",
        "Assign the FimA fimbrial type (I to V) of a P. gingivalis genome.",
        analyze, min_identity=60.0, min_coverage=50.0, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
