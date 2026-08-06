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
    hits = pgcore.screen(query, ("kantigen_region", "kantigen_gene"),
                         require_nucleotide=True, min_identity=min_identity,
                         min_coverage=min_coverage, evalue=evalue,
                         threads=threads)
    regions, genes = [], []
    for h in hits:
        a = h["annotation"]
        # match score balances identity and coverage so a perfect full-length
        # match outranks a longer reference matched at lower identity
        score = round(h["pident"] * min(h["scov"], 100.0) / 100.0, 1)
        rec = dict(record=h["sseqid"], K_type=a.get("Typing_type", "K?"),
                   source_strain=a.get("Source_strain", ""),
                   identity=round(h["pident"], 1), coverage=round(h["scov"], 1),
                   match_score=score)
        if a.get("Record_type") == "kantigen_region":
            regions.append(rec)
        else:
            genes.append(rec)
    regions.sort(key=lambda r: r["match_score"], reverse=True)

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
