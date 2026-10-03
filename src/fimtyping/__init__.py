"""
fimtyping - FimA major-fimbriae typing for Porphyromonas gingivalis,
covering all six published types: I, Ib, II, III, IV and V.

    import fimtyping
    result = fimtyping.analyze("my_assembly.fna")
    print(result.fimA_type)

Types II to V are called by best match against the reference panel, ranked by
identity x coverage. Type Ib is not, and cannot be: Nakagawa et al. (2002)
state plainly that type Ib cannot be distinguished from type I by
type-specific PCR because the two are 97.1 % identical, and they define Ib
instead by a diagnostic RsaI (AfaI) site. Their assay amplifies types I and Ib
alike with the Ib-F/Ib-R pair and then digests: only type Ib is cut.

That distinction is reproduced here in silico rather than approximated by
sequence identity, because identity gets it wrong. Across 104 public genomes
the published assay calls fifteen as Ib, of which whole-gene best match
recovers only six at a one-point margin, leaves four inside the margin, and
assigns the remaining five to type I outright, since they sit marginally
closer to the type I reference overall while carrying the diagnostic site.
A type defined by one base is not recoverable from a whole-gene average.
"""
import re

import pgcore

__version__ = "1.1.0"
__all__ = ["analyze", "__version__"]

# Nakagawa et al. 2002, J Periodont Res 37:425-432. The primer pair amplifies
# types I and Ib and nothing else; RsaI cuts the type Ib product only.
_IB_F = "CAGCAGAGCCAAAAACAATCG"
_IB_R = "TGTCAGATAATTAGCGTCTGC"
_RSAI = "GTAC"


def _rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _ib_assay(query, max_product=600):
    """Nakagawa's type Ib assay in silico.

    Returns (call, product_length, fragments) where call is "Ib", "I" or None.
    None means no product, which is the expected result for types II to V and
    carries no information about them.
    """
    rev = _rc(_IB_R)
    seqs = []
    with open(query) as fh:
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
    for raw in seqs:
        up = raw.upper()
        for strand in (up, _rc(up)):
            start = 0
            while True:
                i = strand.find(_IB_F, start)
                if i < 0:
                    break
                j = strand.find(rev, i)
                if 0 <= j - i <= max_product:
                    product = strand[i:j + len(rev)]
                    cuts = [m.start() + 2 for m in re.finditer(_RSAI, product)]
                    if cuts:
                        frags, prev = [], 0
                        for c in cuts:
                            frags.append(c - prev)
                            prev = c
                        frags.append(len(product) - prev)
                        return "Ib", len(product), sorted(frags, reverse=True)
                    return "I", len(product), [len(product)]
                start = i + 1
    return None, 0, []


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

    # I and Ib are separated by the published restriction assay, not by
    # identity. The assay is authoritative where it gives a product; where it
    # does not, the genome is not type I or Ib and the best match stands.
    ib_call, ib_len, ib_frags = _ib_assay(query)
    assay_note = None
    if ib_call:
        if call in ("I", "Ib") and call != ib_call:
            assay_note = ("best match gave type %s, but the Nakagawa RsaI assay "
                          "gives type %s and is authoritative for this pair"
                          % (call, ib_call))
        call = ib_call
        headline = ("fimA type: %s  (Nakagawa Ib assay: %d bp product, %s)"
                    % (call, ib_len,
                       "RsaI cut to " + " + ".join(str(f) for f in ib_frags)
                       if ib_call == "Ib" else "uncut by RsaI"))
    res = pgcore.Result(
        "FimA fimbrial typing",
        ["fimA_type", "reference", "identity", "coverage", "match_score",
         "query_contig"],
        records,
        summary={"fimA_type": call or "untypeable",
                 "ib_assay": ib_call or "no product (not type I or Ib)",
                 "ib_product_bp": str(ib_len) if ib_len else "-",
                 "ib_rsai_fragments": " + ".join(str(f) for f in ib_frags) or "-"},
        headline=headline)
    res.fimA_type = call
    res.ib_assay = ib_call
    res.ib_fragments = ib_frags
    res.notes = ([assay_note] if assay_note else []) + [
        "Types I and Ib are separated by the RsaI (AfaI) assay of Nakagawa "
        "et al. 2002, not by sequence identity, because the two differ by "
        "roughly 2 % over the gene and the type is defined by the restriction "
        "site. Types II to V give no product in that assay and are called by "
        "best match against the reference panel."]
    return res
