"""
pgmlst - PubMLST seven-locus sequence typing for Porphyromonas gingivalis.

    import pgmlst
    result = pgmlst.analyze("my_assembly.fna")
    print(result.ST, result.profile)

The scheme (ftsQ, gpdxJ, hagB, mcmA, pepO, pga, recA) was defined by Enersen
et al. 2008 and is curated at PubMLST. Allele sequences and the profile table
are redistributed with this package; the versions bundled here were retrieved
in October 2026 and hold 210 of the 220 defined sequence types, together with
the complete allele set: 42/42/41/31/40/29/16 alleles for
ftsQ/gpdxJ/hagB/mcmA/pepO/pga/recA (241 in all). The ten sequence types not
bundled here, ST201 to ST207 and ST214 to ST216, are profiles recovered from
the PRJDB2925 assemblies whose numbers were assigned after this release was
cut; they will be added once the designations are confirmed.

Part of the scheme came from this work. ST178 to ST199 were defined from
twenty-two profiles recovered from public NCBI RefSeq assemblies, and seven
alleles recovered from the same assemblies were defined as ftsQ 41, gpdxJ 38
and 39, hagB 38, mcmA 31, pga 28 and recA 15. A further eleven alleles were
defined as gpdxJ 40 to 42, hagB 39 to 41, pepO 38 to 40, pga 29 and recA 16,
and ten of the profiles they completed were defined as ST208 to ST213 and
ST217 to ST220.

A sequence type places a strain in the nomenclature the rest of bacteriology
uses, which is why gingiCHARM reports one. It should be read alongside the
fimA, mfa, rag and K-antigen calls rather than in place of them. The seven
loci are housekeeping genes chosen for neutrality, so an ST carries no
information about the fimbrial, capsular or rag loci; and with 220 STs defined
over 221 deposited isolates and no clonal complexes assigned, the scheme
identifies individual strains far better than it groups them into lineages.
"""
import csv
import os
import subprocess

import pgcore

__version__ = "1.0.0"
__all__ = ["analyze", "LOCI", "__version__"]

LOCI = ["ftsQ", "gpdxJ", "hagB", "mcmA", "pepO", "pga", "recA"]
N_ST = 220   # sequence types defined in the scheme, October 2026
N_ST_BUNDLED = 210   # of those, the profiles held in data/profiles.tsv
N_ALLELES = 241
SOURCE = "PubMLST (https://pubmlst.org/organisms/porphyromonas-gingivalis)"

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _profiles():
    with open(os.path.join(_DATA, "profiles.tsv")) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


_SITE_BIN = 3000     # nt window used to group hits into one genomic locus


def _call_locus(query, locus, min_coverage=95.0, min_identity=90.0, threads=4):
    """Call one locus, returning (candidates, identity, coverage, status).

    Four outcomes are kept apart because they mean different things:

      known      full-length, 100% identical to exactly one deposited allele
      ambiguous  the genome carries more than one copy of the locus and the
                 copies match different deposited alleles
      novel      full-length but not identical - a real new allele, usually one
                 or two SNPs from its nearest deposited neighbour
      partial    the locus is not covered end to end, which almost always means
                 a contig boundary runs through it rather than that the gene is
                 missing, so no allele is claimed

    The ambiguous case is the one that bites. hagB lies in a tandem pair of
    near-identical haemagglutinin genes about 3.1 kb apart, and 48 of 109
    public assemblies carry both copies. Usually both give the same allele, but
    when they differ, taking whichever hit BLAST reported first silently yields
    the wrong allele and therefore the wrong ST. Every candidate is returned so
    the caller can resolve it against the profile table rather than guess.
    """
    blastn = pgcore.blast.find_blast("blastn")
    fasta = os.path.join(_DATA, locus + ".fas")
    out = subprocess.run(
        [blastn, "-query", fasta, "-subject", query, "-dust", "no",
         "-num_threads", str(threads), "-perc_identity", str(min_identity),
         "-max_target_seqs", "500",
         "-outfmt", "6 qseqid pident length qlen sseqid sstart send"],
        capture_output=True, text=True)

    hits = []
    for line in out.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 7:
            continue
        qid, pid, alen, qlen = p[0], float(p[1]), int(p[2]), int(p[3])
        contig, s1, e1 = p[4], int(p[5]), int(p[6])
        hits.append({"allele": qid.split("_")[-1], "pid": pid,
                     "cov": 100.0 * alen / qlen,
                     "site": (contig, min(s1, e1) // _SITE_BIN)})
    if not hits:
        return [], 0.0, 0.0, "absent"

    full = [h for h in hits if h["cov"] >= min_coverage]
    if not full:
        b = max(hits, key=lambda h: (h["pid"], h["cov"]))
        return [b["allele"]], b["pid"], b["cov"], "partial"

    top = max(h["pid"] for h in full)
    best = [h for h in full if h["pid"] == top]
    alleles = sorted({h["allele"] for h in best},
                     key=lambda x: int(x) if x.isdigit() else x)
    cov = max(h["cov"] for h in best)
    if top < 100.0:
        return alleles, top, cov, "novel"
    return alleles, top, cov, ("known" if len(alleles) == 1 else "ambiguous")


def analyze(query, *, min_identity=90.0, min_coverage=95.0, threads=4):
    """Call the PubMLST sequence type of a P. gingivalis genome assembly.

    `query` must be nucleotide (a genome assembly or contigs). `min_coverage`
    is the fraction of the allele that must align before an allele is called
    at all; `min_identity` only bounds the BLAST search and does not relax the
    exact-match requirement for assigning a deposited allele number. Returns a
    pgcore.Result, with the call also on result.ST and result.profile.
    """
    if not os.path.exists(query):
        raise FileNotFoundError(query)

    records, cands = [], {}
    novel, partial, ambig = [], [], []
    for locus in LOCI:
        alleles, pid, cov, status = _call_locus(query, locus,
                                                min_coverage=min_coverage,
                                                min_identity=min_identity,
                                                threads=threads)
        cands[locus] = alleles if status in ("known", "ambiguous") else []
        shown = "/".join(alleles) if alleles else "-"
        if status == "novel":
            novel.append(locus)
            shown = "~" + shown
        elif status in ("partial", "absent"):
            partial.append(locus)
            shown = ("?" + shown) if alleles else "-"
        elif status == "ambiguous":
            ambig.append(locus)
        records.append(dict(locus=locus, allele=shown,
                            identity=round(pid, 1), coverage=round(cov, 1),
                            status=status, candidates=",".join(alleles)))

    st, resolved, note = None, None, ""
    if not novel and not partial:
        matches = [pr for pr in _profiles()
                   if all(pr[l] in cands[l] for l in LOCI)]
        if len(matches) == 1:
            st = matches[0]["ST"]
            resolved = {l: matches[0][l] for l in LOCI}
        elif len(matches) > 1:
            note = " (ambiguous: matches ST%s)" % ", ST".join(m["ST"] for m in matches)
        elif not matches:
            note = " (all alleles known, combination not in the bundled table)"
        if st and ambig:
            note = " (resolved; %s is multi-copy in this assembly)" % ", ".join(ambig)
        call = "ST" + st if st else "novel ST"
    elif partial:
        call = "incomplete"
        note = " (%s not covered end to end - likely a contig break)" % ", ".join(partial)
    else:
        call = "novel ST"
        note = " (novel allele%s at %s)" % ("s" if len(novel) > 1 else "",
                                            ", ".join(novel))

    if resolved:
        for r in records:
            if r["status"] == "ambiguous":
                r["allele"] = resolved[r["locus"]]
    profile = "-".join(r["allele"] for r in records)

    notes = ["Scheme and allele definitions from %s." % SOURCE,
             "The scheme defines %d sequence types across 221 "
             "deposited isolates and assigns no clonal complexes, so an ST identifies a "
             "strain but says little about lineage. Read it alongside the "
             "fimA, mfa, rag and K-antigen calls." % N_ST]
    if call == "novel ST" and not novel and not partial:
        notes.append(
                "The profile is complete and every allele is deposited, but it "
                "matches no profile in the bundled table, which holds %d of the "
                "%d defined sequence types. Check the profile against PubMLST "
                "before treating it as new." % (N_ST_BUNDLED, N_ST))
    if novel and not partial:
        notes.append("Novel allele%s at %s; the nearest deposited allele is "
                     "shown after the tilde. Worth submitting to PubMLST."
                     % ("s" if len(novel) > 1 else "", ", ".join(novel)))
    if partial:
        notes.append("No sequence type is reported because %s could not be "
                     "recovered at full length. This reflects assembly "
                     "contiguity, not gene absence." % ", ".join(partial))
    if ambig:
        notes.append("%s is present in more than one copy and the copies carry "
                     "different alleles (%s). %s"
                     % (", ".join(ambig),
                        "; ".join("%s=%s" % (l, "/".join(cands[l])) for l in ambig),
                        "Resolved against the profile table." if resolved else
                        "Not resolvable against the profile table; all "
                        "candidates are reported."))

    res = pgcore.Result(
        "MLST (PubMLST seven-locus scheme)",
        ["locus", "allele", "identity", "coverage", "status", "candidates"],
        records,
        summary={"ST": call, "profile": profile,
                 "novel_loci": ", ".join(novel) or "none",
                 "partial_loci": ", ".join(partial) or "none",
                 "multicopy_loci": ", ".join(ambig) or "none"},
        headline="MLST: %s%s  [%s]" % (call, note, profile))
    res.ST = ("ST" + st) if st else None
    res.profile = profile
    res.novel_loci = novel
    res.partial_loci = partial
    res.notes = notes
    return res
