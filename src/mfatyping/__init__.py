"""
mfatyping - Mfa1 minor-fimbriae operon typing for Porphyromonas gingivalis.

Reports the full operon profile (mfa1 + mfa2 + mfa3 + mfa4 + mfa5) following
Watanabe T et al. (PLOS One 2021) and Nagano K et al. (J Oral Microbiol 2023).
mfa1, mfa2, mfa3 and mfa4 split into 53 and 70 (with 70A/70B subtypes for
mfa1); mfa5 has its own A/B/C/D/E classification independent of the rest.
A gene whose locus cannot be confidently typed - either genuinely absent or
present in a form too divergent to match the reference panel - is reported
as 'X'.
The reference panel currently covers mfa1 (53, 70A, 70B), mfa2/3/4 (53, 70)
and mfa5 (A1, plus the divergent Ando allele).

    import mfatyping
    result = mfatyping.analyze("my_assembly.fna")
    print(result.mfa1_genotype)        # '70A'
    print(result.operon_profile)       # ('70A', '70', '70', '70', 'A1')
"""
import pgcore

__version__ = "1.1.0"
__all__ = ["analyze", "__version__"]

# all five operon components, in operon order; the dict value is the attribute
# name attached to the Result for convenience
_OPERON = [
    ("typing_mfa1", "mfa1_genotype"),
    ("typing_mfa2", "mfa2_genotype"),
    ("typing_mfa3", "mfa3_genotype"),
    ("typing_mfa4", "mfa4_genotype"),
    ("typing_mfa5", "mfa5_genotype"),
]


def _best(hits, rtype):
    """Best hit of a given Record_type, ranked by identity x coverage."""
    sub = [h for h in hits if h["annotation"].get("Record_type") == rtype]
    if not sub:
        return None
    return max(sub, key=lambda h: h["pident"] * min(h["scov"], 100.0) / 100.0)


def analyze(query, *, min_identity=80.0, min_coverage=70.0, evalue=1e-10,
            threads=4):
    """Type the mfa operon of a P. gingivalis genome.

    `query` must be nucleotide. Returns a pgcore.Result with one row per
    operon gene. Genes with no qualifying hit are reported as 'X' (unknown
    or missing). X is not a confirmed deletion: it can mean the gene is
    truly absent, present but too divergent to match the reference panel,
    or in an incomplete part of the assembly. Manual follow-up is needed
    to distinguish those cases.
    """
    record_types = tuple(rt for rt, _ in _OPERON)
    hits = pgcore.screen(query, record_types, require_nucleotide=True,
                         min_identity=min_identity, min_coverage=min_coverage,
                         evalue=evalue, threads=threads)

    records = []
    calls = {}
    for rtype, attr in _OPERON:
        gene = rtype.replace("typing_", "")
        best = _best(hits, rtype)
        if best is None:
            calls[gene] = "X"
            records.append(dict(gene=gene, genotype="X", reference="-",
                                identity=0.0, coverage=0.0, match_score=0.0,
                                query_contig="-"))
            continue
        a = best["annotation"]
        gtype = a.get("Typing_type", "?")
        calls[gene] = gtype
        records.append(dict(gene=gene, genotype=gtype,
                            reference=best["sseqid"],
                            identity=round(best["pident"], 1),
                            coverage=round(best["scov"], 1),
                            match_score=round(best["pident"]
                                              * min(best["scov"], 100.0)
                                              / 100.0, 1),
                            query_contig=best["qseqid"]))

    profile = tuple(calls[g] for g in ("mfa1", "mfa2", "mfa3", "mfa4", "mfa5"))
    profile_str = " / ".join(profile)

    # cross-check: mfa1/2/3/4 should co-segregate (all 53 or all 70-family).
    # mfa5 is independent.
    body = [calls[g] for g in ("mfa1", "mfa2", "mfa3", "mfa4")]
    family = lambda g: ("53" if g == "53"
                        else "70" if g and g.startswith("70")
                        else "X" if g == "X"
                        else "?")
    body_families = {family(g) for g in body if g != "X"}
    discordant = len(body_families - {"X"}) > 1

    if all(c == "X" for c in profile):
        headline = "mfa operon: untypeable - no mfa1-5 hits above thresholds"
    else:
        headline = "mfa operon: " + profile_str
        if discordant:
            headline += "  [WARNING: mfa1/2/3/4 disagree on 53 vs 70 family]"

    res = pgcore.Result(
        "Mfa1 minor-fimbriae operon typing",
        ["gene", "genotype", "reference", "identity", "coverage",
         "match_score", "query_contig"],
        records,
        summary={"mfa1_genotype": calls["mfa1"],
                 "mfa2_genotype": calls["mfa2"],
                 "mfa3_genotype": calls["mfa3"],
                 "mfa4_genotype": calls["mfa4"],
                 "mfa5_genotype": calls["mfa5"],
                 "operon_profile": profile_str,
                 "co_segregation": "OK" if not discordant else "discordant"},
        headline=headline)

    # convenience attrs - keep mfa1_genotype for backward compatibility
    res.mfa1_genotype = calls["mfa1"]
    res.mfa2_genotype = calls["mfa2"]
    res.mfa3_genotype = calls["mfa3"]
    res.mfa4_genotype = calls["mfa4"]
    res.mfa5_genotype = calls["mfa5"]
    res.operon_profile = profile
    return res
