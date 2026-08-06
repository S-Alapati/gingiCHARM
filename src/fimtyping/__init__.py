"""
fimtyping - FimA major-fimbriae typing for Porphyromonas gingivalis,
covering types I to V.

    import fimtyping
    result = fimtyping.analyze("my_assembly.fna")
    print(result.fimA_type)
"""
import pgcore

__version__ = "1.0.0"
__all__ = ["analyze", "__version__"]


def analyze(query, *, min_identity=60.0, min_coverage=50.0, evalue=1e-10,
            threads=4):
    """Assign the FimA fimbrial type of a P. gingivalis genome assembly.

    `query` must be nucleotide (a genome assembly, contigs or CDS). The call is
    the best-matching reference type; identity and coverage are reported so the
    user can judge how strong the match is. Returns a pgcore.Result, with the
    type also on result.fimA_type.
    """
    hits = pgcore.screen(query, ("typing_fimA",), require_nucleotide=True,
                         min_identity=min_identity, min_coverage=min_coverage,
                         evalue=evalue, threads=threads)
    records = []
    for h in hits:
        a = h["annotation"]
        records.append(dict(fimA_type=a.get("Typing_type", "?"),
                            reference=h["sseqid"],
                            identity=round(h["pident"], 1),
                            coverage=round(h["scov"], 1),
                            match_score=round(h["pident"] * min(h["scov"], 100.0)
                                              / 100.0, 1),
                            query_contig=h["qseqid"]))
    # rank by identity x coverage, not bitscore: fimA type references differ
    # widely in length and a longer reference must not outscore a better match
    records.sort(key=lambda r: r["match_score"], reverse=True)
    if records:
        call = records[0]["fimA_type"]
        headline = ("fimA type: %s  (%.1f%% identity, %.1f%% coverage of reference)"
                    % (call, records[0]["identity"], records[0]["coverage"]))
    else:
        call = None
        headline = "fimA type: untypeable - no reference hit above thresholds"
    res = pgcore.Result(
        "FimA fimbrial typing",
        ["fimA_type", "reference", "identity", "coverage", "match_score",
         "query_contig"],
        records,
        summary={"fimA_type": call or "untypeable"},
        headline=headline)
    res.fimA_type = call
    return res
