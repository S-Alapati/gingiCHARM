"""Command-line interface: `pgmlst <genome.fna>`."""
import sys

from pgcore.report import run_cli

from . import analyze


def main(argv=None):
    return run_cli(
        "pgmlst",
        "Call the PubMLST seven-locus sequence type of a P. gingivalis genome.",
        analyze, min_identity=90.0, min_coverage=95.0, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
