"""
pgmlst - PubMLST seven-locus sequence typing for Porphyromonas gingivalis.

    import pgmlst
    result = pgmlst.analyze("my_assembly.fna")
    print(result.ST, result.profile)

The scheme (ftsQ, gpdxJ, hagB, mcmA, pepO, pga, recA) was defined by Enersen
et al. 2008 and is curated at PubMLST. Allele sequences and the profile table
are redistributed with this package; the versions bundled here were retrieved
on 14 August 2026. PubMLST serves unauthenticated clients only data submitted
on or before 31 December 2024, so the bundled copy holds 172 of the 177
sequence types the scheme contained on 5 February 2026. A profile called novel
here may therefore be a recently deposited type; an authenticated user can
refresh the bundled files.

A sequence type places a strain in the nomenclature the rest of bacteriology
uses, which is why gingiCHARM reports one. It should be read alongside the
fimA, mfa, rag and K-antigen calls rather than in place of them. The seven
loci are housekeeping genes chosen for neutrality, so an ST carries no
information about the fimbrial, capsular or rag loci; and with 177 STs defined
over roughly 190 deposited isolates and no clonal complexes assigned, the scheme
identifies individual strains far better than it groups them into lineages.
"""
import csv
import os
import subprocess

import pgcore

__version__ = "1.0.0"
__all__ = ["analyze", "LOCI", "__version__"]

LOCI = ["ftsQ", "gpdxJ", "hagB", "mcmA", "pepO", "pga", "recA"]
N_ST = 177   # scheme total at 2026-02-05; 172 in the public release
SOURCE = "PubMLST (https://pubmlst.org/organisms/porphyromonas-gingivalis)"

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _profiles():
    with open(os.path.join(_DATA, "profiles.tsv")) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _call_locus(query, locus, min_coverage=95.0, min_identity=90.0, threads=4):
    """Call one locus, returning (allele, identity, coverage, status).

    Three outcomes are kept apart because they mean different things:

      known    full-length and 100% identical to a deposited allele
      novel    full-length but not identical - a real new allele, usually one
               or two SNPs from its nearest deposited neighbour, and worth
               submitting to PubMLST
      partial  the locus is not covered end to end, which almost always means
               a contig boundary runs through it rather than that the gene is
               missing, so no allele is claimed

    An exact-match test that reported only "matched / did not match" would put
    the last two in the same bucket, and a submittable result and an assembly
    limitation are not the same thing.
    """
    blastn = pgcore.blast.find_blast("blastn")
    fasta = os.path.join(_DATA, locus + ".fas")
    out = subprocess.run(
        [blastn, "-query", fasta, "-subject", query, "-dust", "no",
         "-num_threads", str(threads), "-perc_identity", str(min_identity),
         "-max_target_seqs", "100",
         "-outfmt", "6 qseqid pident length qlen"],
        capture_output=True, text=True)
    best = None
    for line in out.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 4:
            continue
        qid, pid, alen, qlen = p[0], float(p[1]), int(p[2]), int(p[3])
        cov = 100.0 * alen / qlen
        full = cov >= min_coverage
        score = (1 if (full and pid == 100.0) else 0, pid, alen)
        if best is None or score > best[0]:
            best = (score, qid, pid, cov, full)
    if best is None:
        return None, 0.0, 0.0, "absent"
    _, qid, pid, cov, full = best
    allele = qid.split("_")[-1]
    if not full:
        return allele, pid, cov, "partial"
    return allele, pid, cov, ("known" if pid == 100.0 else "novel")


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

    records, alleles, novel, partial = [], {}, [], []
    for locus in LOCI:
        allele, pid, cov, status = _call_locus(query, locus,
                                               min_coverage=min_coverage,
                                               min_identity=min_identity,
                                               threads=threads)
        alleles[locus] = allele if status == "known" else None
        if status == "novel":
            novel.append(locus)
            shown = "~" + allele
        elif status == "partial":
            partial.append(locus)
            shown = "?" + allele
        elif status == "absent":
            partial.append(locus)
            shown = "-"
        else:
            shown = allele
        records.append(dict(locus=locus, allele=shown,
                            identity=round(pid, 1), coverage=round(cov, 1),
                            status=status))

    profile = "-".join(r["allele"] for r in records)
    st = None
    if not novel and not partial:
        for prof in _profiles():
            if all(prof[l] == alleles[l] for l in LOCI):
                st = prof["ST"]
                break

    if st:
        call, note = "ST" + st, ""
    elif partial:
        call = "incomplete"
        note = " (%s not covered end to end - likely a contig break)" % ", ".join(partial)
    elif novel:
        call = "novel ST"
        note = " (novel allele%s at %s)" % ("s" if len(novel) > 1 else "",
                                            ", ".join(novel))
    else:
        call = "novel ST"
        note = " (all alleles known, combination not in PubMLST)"

    notes = ["Scheme and allele definitions from %s." % SOURCE,
             "The scheme defines %d sequence types across roughly 190 "
             "deposited isolates and assigns no clonal complexes, so an ST identifies a "
             "strain but says little about lineage. Read it alongside the "
             "fimA, mfa, rag and K-antigen calls." % N_ST]
    if novel and not partial:
        notes.append("Novel allele%s at %s; the nearest deposited allele is "
                     "shown after the tilde. Worth submitting to PubMLST."
                     % ("s" if len(novel) > 1 else "", ", ".join(novel)))
    if partial:
        notes.append("No sequence type is reported because %s could not be "
                     "recovered at full length. This reflects assembly "
                     "contiguity, not gene absence." % ", ".join(partial))

    res = pgcore.Result(
        "MLST (PubMLST seven-locus scheme)",
        ["locus", "allele", "identity", "coverage", "status"],
        records,
        summary={"ST": call, "profile": profile,
                 "novel_loci": ", ".join(novel) or "none",
                 "partial_loci": ", ".join(partial) or "none"},
        headline="MLST: %s%s  [%s]" % (call, note, profile))
    res.ST = ("ST" + st) if st else None
    res.profile = profile
    res.novel_loci = novel
    res.partial_loci = partial
    res.notes = notes
    return res
