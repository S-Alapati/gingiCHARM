"""
kantigen - capsular K-antigen typing for Porphyromonas gingivalis.

Reference K-types in the database: K1, K3, K4, K6 and K- (non-encapsulated).
Strains of other serotypes are present but not yet K-typed (K?).

    import kantigen
    result = kantigen.analyze("my_assembly.fna")
    print(result.k_type)
"""
from collections import Counter

import pgcore

__version__ = "1.0.0"
__all__ = ["analyze", "__version__"]

_REAL = {"K1", "K2", "K3", "K4", "K5", "K6", "K7", "K-"}


def analyze(query, *, min_identity=80.0, min_coverage=40.0, evalue=1e-20,
            threads=4):
    """Assign the capsular K-antigen type of a P. gingivalis genome.

    Uses the cps capsular locus: the best-matching reference cps region gives
    the primary call, refined by per-gene cps concordance. `query` must be
    nucleotide. Returns a pgcore.Result; `result.k_type` holds the call.
    """
    # The cps locus is 12-18 kb and serotype-variable, so a query aligns to a
    # reference in many fragments rather than one block. Thresholding each
    # fragment separately would reject every genuine match, so hits are
    # collected at a permissive per-fragment setting and aggregated per
    # reference before the identity and coverage thresholds are applied.
    hits = pgcore.screen(query, ("kantigen_region", "kantigen_gene"),
                         require_nucleotide=True, min_identity=min_identity,
                         min_coverage=1.0, evalue=evalue, threads=threads)
    agg = {}
    for h in hits:
        a = h["annotation"]
        rid = h["sseqid"]
        d = agg.setdefault(rid, dict(
            record=rid, K_type=a.get("Typing_type", "K?"),
            source_strain=a.get("Source_strain", ""),
            record_type=a.get("Record_type", ""),
            aligned=0.0, idsum=0.0, reflen=float(a.get("Seq_length_nt") or 0) or None))
        # scov is the fraction of the reference this fragment covers; recover
        # the fragment length from it when the record length is unavailable
        frag = h.get("length") or (h["scov"] / 100.0 * (d["reflen"] or 0))
        d["aligned"] += frag
        d["idsum"] += h["pident"] * frag
    regions, genes = [], []
    for d in agg.values():
        if not d["aligned"]:
            continue
        ident = d["idsum"] / d["aligned"]
        cov = 100.0 * d["aligned"] / d["reflen"] if d["reflen"] else 0.0
        if ident < min_identity or cov < min_coverage:
            continue
        rec = dict(record=d["record"], K_type=d["K_type"],
                   source_strain=d["source_strain"],
                   identity=round(ident, 1), coverage=round(min(cov, 100.0), 1),
                   aligned_bp=int(d["aligned"]),
                   match_score=round(ident * min(cov, 100.0) / 100.0, 1))
        (regions if d["record_type"] == "kantigen_region" else genes).append(rec)
    # rank on total aligned length: with a variable locus, how much of the
    # reference is recovered separates the serotypes, not identity alone
    regions.sort(key=lambda r: (r["aligned_bp"], r["match_score"]), reverse=True)

    # gene-level concordance among real K-types
    gene_votes = Counter(g["K_type"] for g in genes if g["K_type"] in _REAL)

    call, headline = None, "K-antigen type: untypeable - no cps locus hit"
    if regions:
        top = regions[0]                       # sorted by bitscore
        if top["K_type"] in _REAL:
            call = top["K_type"]
            headline = ("K-antigen type: %s  (best cps-region match: %s, "
                        "%.1f%% coverage)" % (call, top["source_strain"],
                                              top["coverage"]))
        else:
            call = "untyped-reference"
            headline = ("K-antigen type: closest cps region is an un-serotyped "
                        "strain (%s); not assignable to K1-K7" % top["source_strain"])
            if gene_votes:
                best = gene_votes.most_common(1)[0]
                headline += " | cps genes lean %s (%d)" % best
    res = pgcore.Result(
        "Capsular K-antigen typing",
        ["record", "K_type", "source_strain", "identity", "coverage",
         "match_score"],
        regions,
        summary={"k_type": call or "untypeable",
                 "cps-region references hit": len(regions),
                 "cps-gene concordance": dict(gene_votes) or "-"},
        headline=headline)
    res.k_type = call
    return res
