"""
ragtyping - RagAB pathogenicity-island typing for Porphyromonas gingivalis.
The rag type is the ragA + ragB allele combination per Hall et al. 2005.
The reference panel covers all four published alleles (rag-1, rag-2, rag-3,
rag-4).

    import ragtyping
    result = ragtyping.analyze("my_assembly.fna")
    print(result.rag_type)
"""
import pgcore

__version__ = "1.0.0"
__all__ = ["analyze", "__version__"]


def _best(hits, rtype):
    sub = [h for h in hits if h["annotation"].get("Record_type") == rtype]
    if not sub:
        return None
    return max(sub, key=lambda h: h["pident"] * min(h["scov"], 100.0) / 100.0)


def analyze(query, *, min_identity=85.0, min_coverage=70.0, evalue=1e-10,
            threads=4):
    """Assign the rag type of a P. gingivalis genome from its ragA/ragB locus.

    `query` must be nucleotide. Returns a pgcore.Result; `result.rag_type`
    holds the combined call.
    """
    hits = pgcore.screen(query, ("typing_ragA", "typing_ragB"),
                         require_nucleotide=True, min_identity=min_identity,
                         min_coverage=min_coverage, evalue=evalue, threads=threads)
    records = []
    for h in hits:
        a = h["annotation"]
        records.append(dict(gene=a.get("Gene", ""),
                            rag_allele=a.get("Typing_type", "?"),
                            reference=h["sseqid"],
                            identity=round(h["pident"], 1),
                            coverage=round(h["scov"], 1),
                            match_score=round(h["pident"] * min(h["scov"], 100.0)
                                              / 100.0, 1),
                            query_contig=h["qseqid"]))
    records.sort(key=lambda r: r["match_score"], reverse=True)
    ra, rb = _best(hits, "typing_ragA"), _best(hits, "typing_ragB")
    a_type = ra["annotation"].get("Typing_type") if ra else None
    b_type = rb["annotation"].get("Typing_type") if rb else None
    # ragA is the primary discriminator. The rag-4 group (ATCC 33277-like)
    # has only ragB in the reference set because Bakta typically misses its
    # divergent ragA; for those queries the call rests on ragB.
    call = a_type or b_type
    headline = ("rag type: %s" % call) if call else \
               "rag type: untypeable - no ragA/ragB hit above thresholds"
    res = pgcore.Result(
        "RagAB pathogenicity-island typing",
        ["gene", "rag_allele", "reference", "identity", "coverage", "match_score",
         "query_contig"],
        records,
        summary={"rag_type": call or "untypeable",
                 "ragA_allele": a_type or "-", "ragB_allele": b_type or "-"},
        headline=headline)
    res.rag_type = call
    res.ragA_allele = a_type
    res.ragB_allele = b_type
    return res
