"""Command-line interface: `pgvfdb <genome.fna>`."""
import sys

from pgcore.report import run_cli

from . import analyze


def main(argv=None):
    return run_cli(
        "pgvfdb",
        "Screen a Porphyromonas gingivalis genome for virulence factors.",
        analyze, min_identity=80.0, min_coverage=70.0, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
