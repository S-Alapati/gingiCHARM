"""
pgvfdb - virulence-factor screening for Porphyromonas gingivalis.

    import pgvfdb
    result = pgvfdb.analyze("my_assembly.fna")
    print(result)                 # summary
    result.to_tsv("vf_report.tsv")
"""
from collections import Counter

import pgcore

__version__ = "1.0.0"
__all__ = ["analyze", "__version__"]


def analyze(query, *, min_identity=80.0, min_coverage=70.0, evalue=1e-5,
            threads=4):
    """Screen a P. gingivalis genome / proteome for virulence factors.

    `query` may be a genome assembly, contigs, CDS (nucleotide) or predicted
    proteins (amino acid) in FASTA format. Returns a pgcore.Result with one
    row per virulence factor detected.
    """
    hits = pgcore.screen(query, ("VF_screen", "VF_allele"),
                         min_identity=min_identity, min_coverage=min_coverage,
                         evalue=evalue, threads=threads)
    vf = {}
    for h in hits:
        a = h["annotation"]
        name = a.get("VF_name") or a.get("Gene") or h["sseqid"]
        cur = vf.get(name)
        if cur is None or h["pident"] > cur["identity"]:
            lo, hi = sorted((h["qstart"], h["qend"]))
            vf[name] = dict(
                VF_name=name,
                VF_category=a.get("VF_category", ""),
                Biological_Category=a.get("Biological_Category", ""),
                Evidence_tier=a.get("VF_evidence_tier", ""),
                identity=round(h["pident"], 1),
                coverage=round(h["scov"], 1),
                query_locus="%s:%d-%d" % (h["qseqid"], lo, hi),
                hit_record=h["sseqid"])
    records = sorted(vf.values(), key=lambda r: (r["VF_category"], r["VF_name"]))
    bycat = Counter(r["VF_category"] for r in records)
    res = pgcore.Result(
        "Pg virulence-factor screen",
        ["VF_name", "VF_category", "Biological_Category", "Evidence_tier",
         "identity", "coverage", "query_locus", "hit_record"],
        records,
        summary={"virulence factors detected": len(records),
                 "VF categories represented": len(bycat)},
        headline="%d virulence factors detected" % len(records))
    res.by_category = dict(bycat)
    return res
