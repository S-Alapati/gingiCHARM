"""Access to the bundled Pg unified database."""
import csv
import os
from functools import lru_cache

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def data_dir():
    return _DATA


def nuc_db():
    """Prefix of the combined nucleotide BLAST database."""
    return os.path.join(_DATA, "blast_db", "Pg_unified_nuc")


def prot_db():
    """Prefix of the combined protein BLAST database."""
    return os.path.join(_DATA, "blast_db", "Pg_unified_prot")


def annotation_path():
    return os.path.join(_DATA, "Pg_unified_annotation.tsv")


@lru_cache(maxsize=1)
def annotation():
    """Return {Record_ID: row-dict} for the whole unified database."""
    out = {}
    with open(annotation_path()) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["Record_ID"]] = row
    return out


def records_of_type(*types):
    """Return the set of Record_IDs whose Record_type is in `types`."""
    want = set(types)
    return {rid for rid, r in annotation().items() if r["Record_type"] in want}


def counts():
    """Return {Record_type: count} for the bundled database."""
    out = {}
    for r in annotation().values():
        out[r["Record_type"]] = out.get(r["Record_type"], 0) + 1
    return out


def check_ready():
    """Raise RuntimeError if the bundled database files are missing."""
    missing = []
    if not os.path.exists(annotation_path()):
        missing.append("annotation table")
    if not (os.path.exists(nuc_db() + ".ndb") or os.path.exists(nuc_db() + ".nin")):
        missing.append("nucleotide BLAST database")
    if not (os.path.exists(prot_db() + ".pdb") or os.path.exists(prot_db() + ".pin")):
        missing.append("protein BLAST database")
    if missing:
        raise RuntimeError("pgtoolkit database files missing: " + ", ".join(missing))
