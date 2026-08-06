"""Command-line interface: `kantigen <genome.fna>`."""
import sys

from pgcore.report import run_cli

from . import analyze


def main(argv=None):
    return run_cli(
        "kantigen",
        "Assign the capsular K-antigen type of a P. gingivalis genome.",
        analyze, min_identity=80.0, min_coverage=40.0, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
