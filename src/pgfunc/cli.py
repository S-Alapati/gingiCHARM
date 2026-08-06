"""Command-line interface: `pgfunc <genome.fna>`."""
import sys

from pgcore.report import run_cli

from . import analyze


def main(argv=None):
    return run_cli(
        "pgfunc",
        "Assign functional categories to the genes of a P. gingivalis genome.",
        analyze, min_identity=30.0, min_coverage=50.0, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
