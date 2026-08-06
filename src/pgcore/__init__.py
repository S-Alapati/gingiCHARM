"""
pgcore - shared internals for the pgtoolkit Porphyromonas gingivalis
analysis modules (pgvfdb, fimtyping, mfatyping, ragtyping, kantigen, pgfunc).
"""
from . import fasta, blast, db
from .report import Result
from .homology import screen, resolve_loci

__version__ = "1.0.0"
__all__ = ["fasta", "blast", "db", "Result", "screen", "resolve_loci",
           "__version__"]
