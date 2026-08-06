"""Generic homology screen shared by every pgtoolkit analysis module."""
import os

from . import blast, db, fasta


def screen(query_path, record_types, *, min_identity=80.0, min_coverage=70.0,
           evalue=1e-5, threads=4, max_target_seqs=8000,
           require_nucleotide=False):
    """BLAST a query against the unified database and return matching hits.

    A nucleotide query (genome assembly, contigs or CDS) is searched with
    blastn against the nucleotide database; a protein query is searched with
    blastp against the protein database. Input type is detected automatically.

    Parameters
    ----------
    query_path : str
        FASTA file - genome assembly, contigs, CDS, or proteins.
    record_types : iterable of str
        Only hits to records of these Record_types are kept
        (e.g. ("VF_screen", "VF_allele"), ("typing_fimA",)).
    require_nucleotide : bool
        If True, reject a protein query (used by the typing modules, which
        need DNA-level resolution).

    Returns
    -------
    list of dict
        One dict per best-scoring HSP per (query, subject) pair, with all
        BLAST fields plus 'scov' and a merged 'annotation' row.
    """
    if not os.path.isfile(query_path):
        raise FileNotFoundError("query file not found: %s" % query_path)
    db.check_ready()
    seqs = fasta.read_fasta(query_path)
    if not seqs:
        raise ValueError("no FASTA sequences found in %s" % query_path)
    qtype = fasta.guess_type(seqs)

    if qtype == "nt":
        program, target = "blastn", db.nuc_db()
        extra = ["-task", "dc-megablast"]      # tolerant of divergent alleles
    else:
        if require_nucleotide:
            raise ValueError("this analysis needs a nucleotide query "
                             "(genome assembly, contigs or CDS); "
                             "a protein FASTA was supplied")
        program, target = "blastp", db.prot_db()
        extra = None

    hits = blast.run_blast(program, query_path, target, evalue=evalue,
                           threads=threads, max_target_seqs=max_target_seqs,
                           extra=extra)
    wanted = db.records_of_type(*record_types)
    ann = db.annotation()
    best = {}
    for h in hits:
        if h["sseqid"] not in wanted:
            continue
        if h["pident"] < min_identity or h["scov"] < min_coverage:
            continue
        key = (h["qseqid"], h["sseqid"])
        if key not in best or h["bitscore"] > best[key]["bitscore"]:
            best[key] = h
    out = []
    for h in best.values():
        h = dict(h)
        h["annotation"] = ann.get(h["sseqid"], {})
        out.append(h)
    out.sort(key=lambda x: x["bitscore"], reverse=True)
    return out


def resolve_loci(hits):
    """Reduce hits to one per non-overlapping query region.

    Greedy by bitscore: a hit is kept only if its query interval does not
    overlap a higher-scoring hit already kept. This turns raw HSPs against a
    genome assembly into one call per gene, and picks the single best hit per
    gene when the query is already a set of predicted proteins/CDS.
    """
    kept, occupied = [], {}
    for h in sorted(hits, key=lambda x: x["bitscore"], reverse=True):
        s, e = sorted((h["qstart"], h["qend"]))
        ivs = occupied.setdefault(h["qseqid"], [])
        if any(not (e < a or s > b) for a, b in ivs):
            continue
        ivs.append((s, e))
        kept.append(h)
    return kept
