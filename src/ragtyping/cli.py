"""Command-line interface: `ragtyping <genome.fna>`."""
import sys

from pgcore.report import run_cli

from . import analyze


def main(argv=None):
    return run_cli(
        "ragtyping",
        "Assign the rag pathogenicity-island type (rag-1 through rag-4) of a P. gingivalis genome.",
        analyze, min_identity=85.0, min_coverage=70.0, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
