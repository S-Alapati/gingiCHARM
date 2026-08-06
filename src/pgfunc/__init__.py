"""
pgfunc - functional-category assignment for Porphyromonas gingivalis genes.

Transfers the 10-category biological scheme (+ COG) from the core/pan-genome
backbone to the genes of a query genome by homology.

    import pgfunc
    result = pgfunc.analyze("my_assembly.fna")
    result.to_tsv("functional_categories.tsv")
    print(result.by_category)
"""
from collections import Counter

import pgcore

__version__ = "1.0.0"
__all__ = ["analyze", "__version__"]


def analyze(query, *, min_identity=30.0, min_coverage=50.0, evalue=1e-5,
            threads=4):
    """Assign a functional category to each gene of a P. gingivalis genome.

    `query` may be nucleotide (genome / contigs / CDS) or proteins. Returns a
    pgcore.Result with one row per query gene; `result.by_category` holds the
    category breakdown.
    """
    hits = pgcore.screen(query, ("core_gene",),
                         min_identity=min_identity, min_coverage=min_coverage,
                         evalue=evalue, threads=threads)
    genes = pgcore.resolve_loci(hits)        # one call per gene / query region
    records = []
    for h in sorted(genes, key=lambda x: (x["qseqid"],
                                          min(x["qstart"], x["qend"]))):
        a = h["annotation"]
        lo, hi = sorted((h["qstart"], h["qend"]))
        records.append(dict(
            query_gene="%s:%d-%d" % (h["qseqid"], lo, hi),
            Biological_Category=a.get("Biological_Category", ""),
            COG_category=a.get("COG_category", ""),
            Gene=a.get("Gene", ""),
            Product=a.get("Product", ""),
            Core_status=a.get("Core_status", ""),
            core_cluster=h["sseqid"],
            identity=round(h["pident"], 1),
            coverage=round(h["scov"], 1)))
    bycat = Counter(r["Biological_Category"] for r in records)
    res = pgcore.Result(
        "Functional-category assignment",
        ["query_gene", "Biological_Category", "COG_category", "Gene", "Product",
         "Core_status", "core_cluster", "identity", "coverage"],
        records,
        summary={"genes categorized": len(records),
                 "categories represented": len(bycat)},
        headline="%d genes assigned to %d functional categories"
                 % (len(records), len(bycat)))
    res.by_category = dict(bycat.most_common())
    return res
