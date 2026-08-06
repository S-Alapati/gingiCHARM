"""Command-line interface: `mfatyping <genome.fna>`."""
import sys

from pgcore.report import run_cli

from . import analyze


def main(argv=None):
    return run_cli(
        "mfatyping",
        "Type the Mfa1 minor-fimbriae operon (mfa1/mfa2/mfa3/mfa4/mfa5) of a P. gingivalis genome.",
        analyze, min_identity=80.0, min_coverage=70.0, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
