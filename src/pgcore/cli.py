"""`pgtoolkit-info` - show package, database and BLAST+ status."""
import sys

from . import __version__, blast, db


def main(argv=None):
    print("pgtoolkit %s" % __version__)
    print("data directory: %s" % db.data_dir())
    try:
        db.check_ready()
        print("database      : OK")
        for k, v in sorted(db.counts().items()):
            print("    %-16s %d" % (k, v))
        print("    %-16s %d" % ("TOTAL", sum(db.counts().values())))
    except RuntimeError as e:
        print("database      : NOT READY - %s" % e)
    if blast.blast_available():
        print("BLAST+        : found (%s)" % blast.find_blast("blastn"))
    else:
        print("BLAST+        : NOT FOUND - run install.sh to download it")
    print("\nmodules: pgvfdb  fimtyping  mfatyping  ragtyping  kantigen  pgfunc")
    return 0


if __name__ == "__main__":
    sys.exit(main())
