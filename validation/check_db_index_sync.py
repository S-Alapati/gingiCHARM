"""Fail if the shipped BLAST index and the sequence files have drifted apart.

The K-antigen panel was rebuilt on 3 October 2026 and every call came back
untypeable, because the package ships a pre-built BLAST index that the edit had
not regenerated. Nothing detected it: the index was internally valid, merely
describing a database that no longer existed. This check exists so that cannot
recur silently.
"""
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "src", "pgcore", "data")


def ids_in_fasta(path):
    return {l[1:].split()[0] for l in open(path) if l.startswith(">")}


def ids_in_blastdb(db):
    r = subprocess.run(["blastdbcmd", "-db", db, "-entry", "all", "-outfmt", "%a"],
                       capture_output=True, text=True)
    if r.returncode:
        return None
    return {x.strip() for x in r.stdout.split("\n") if x.strip()}


def main():
    fails = []
    for fa, db in [("sequences/Pg_unified_nuc.fasta", "blast_db/Pg_unified_nuc"),
                   ("sequences/Pg_unified_prot.fasta", "blast_db/Pg_unified_prot")]:
        f = os.path.join(DATA, fa)
        d = os.path.join(DATA, db)
        if not os.path.exists(f):
            continue
        want = ids_in_fasta(f)
        got = ids_in_blastdb(d)
        if got is None:
            fails.append("%s: BLAST index missing or unreadable" % db)
            continue
        missing, extra = sorted(want - got), sorted(got - want)
        if missing or extra:
            fails.append("%s: %d record(s) in the FASTA but not the index (%s), "
                         "%d in the index but not the FASTA (%s)"
                         % (db, len(missing), ", ".join(missing[:4]) or "-",
                            len(extra), ", ".join(extra[:4]) or "-"))
        else:
            print("  PASS  %s: %d records, index matches the FASTA" % (db, len(want)))
    # every annotated record that claims a sequence must have one
    ann = os.path.join(DATA, "Pg_unified_annotation.tsv")
    rows = [l.rstrip("\n").split("\t") for l in open(ann)]
    h = {x: i for i, x in enumerate(rows[0])}
    nuc = ids_in_fasta(os.path.join(DATA, "sequences/Pg_unified_nuc.fasta"))
    orphan = [r[h["Record_ID"]] for r in rows[1:]
              if r[h["Has_nucleotide"]] == "Y" and r[h["Record_ID"]] not in nuc]
    if orphan:
        fails.append("%d annotated record(s) claim a nucleotide sequence that is "
                     "absent: %s" % (len(orphan), ", ".join(orphan[:5])))
    else:
        print("  PASS  every record annotated Has_nucleotide=Y has a sequence")
    for f in fails:
        print("  FAIL  " + f)
    print("\n" + ("database and index are in sync" if not fails
                  else "%d problem(s); run makeblastdb to regenerate" % len(fails)))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
