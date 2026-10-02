# MLST module — validation record

*P. gingivalis*, PubMLST seven-locus scheme. Data retrieved 14 August 2026.

## Background: an error in the earlier drafts

Both manuscript drafts previously stated that none of the three red-complex
organisms has a published MLST scheme, and that this had been confirmed against
the PubMLST registry. That claim was wrong for *P. gingivalis*. The absence was
verified for *Prevotella*, *Treponema* and *Tannerella* but never separately for
*Porphyromonas*, and the negative result was generalised across the genus
without being checked. The scheme exists, has existed since 2008, and the paper
that defines it (Enersen *et al.* 2008) was already cited in the gingiCHARM
manuscript for its *fimA* content.

Verified directly against the REST API:

| Item | Value |
|---|---|
| Databases | `pubmlst_pgingivalis_seqdef`, `pubmlst_pgingivalis_isolates` |
| Scheme 1 | MLST |
| Loci | ftsQ, gpdxJ, hagB, mcmA, pepO, pga, recA |
| Sequence types | 177 total (172 in the public, unauthenticated release) |
| Isolates | ~190 |
| Genomes in collection | 0 |
| Clonal complexes | none assigned (empty field on all 172 STs) |
| Alleles per locus (public release) | ftsQ 40, gpdxJ 37, hagB 37, mcmA 30, pepO 37, pga 27, recA 14 |
| Last ST added | 2026-02-05 |

**Access caveat.** The REST API serves unauthenticated clients only records
submitted on or before 31 December 2024. The bundled allele and profile files
are therefore the public subset: 172 of 177 sequence types. Some fraction of
the profiles called novel in the sweep below may correspond to types deposited
during 2025–2026 that are not visible without an account. This does not affect
correctness — a novel call is still a correct statement that the profile is not
in the bundled release — but it means the 71 % novel rate is an upper bound.

Note the genome count. The scheme holds 192 isolate records but no assembled
genomes, so a tree "of the PubMLST genomes" cannot be built from that source.
The phylogeny below therefore uses the 109 RefSeq assemblies instead, which is
both possible and more useful.

## Method

All alleles within a locus are the same length (ftsQ 420, gpdxJ 380, hagB 380,
mcmA 420, pepO 400, pga 400, recA 400 bp) with no indels, so no alignment step
is needed anywhere in the module or in the phylogeny.

Each locus is searched against the assembly with `blastn` (`-dust no`,
`-perc_identity 90`). Three outcomes are kept apart:

| Status | Criterion | Reported as |
|---|---|---|
| known | ≥95 % coverage, 100 % identity | allele number |
| novel | ≥95 % coverage, <100 % identity | `~N` against nearest allele |
| partial | <95 % coverage | `?N`, no ST assigned |

An ST is assigned only when all seven loci are `known` and the combination
appears in `profiles.tsv`.

The known/novel/partial split matters. An earlier version of this module used a
plain exact-match test, which put novel alleles and contig-broken loci in the
same "not exact" bucket. Those mean opposite things: a novel allele is a
finished result worth depositing, a partial locus is an assembly limitation
where no allele should be claimed at all.

## Validation 1 — profile reproduction

Every genome returning a defined ST was checked against the deposited PubMLST
profile row, locus by locus:

| Strain | Called | PubMLST profile for that ST | Match |
|---|---|---|---|
| ATCC 33277 | ST96, 5-9-1-1-1-5-5 | 5-9-1-1-1-5-5 | yes |
| W83 | ST24, 16-6-1-4-1-3-5 | 16-6-1-4-1-3-5 | yes |
| A7A1-28 | ST1, 1-12-1-1-1-1-1 | 1-12-1-1-1-1-1 | yes |

No mismatches across the full sweep.

## Validation 2 — strain identity recovered without metadata

The stronger check. No strain name, accession or metadata enters the
calculation, so if the loci are being recovered correctly then independent
deposits of the same strain must return the same ST. They do:

| ST | Genomes |
|---|---|
| ST96 | ATCC 33277 (×2), DSM 20709, NBRC 115150, 381 |
| ST24 | W83 (×2), W50 (×2), W50 BE1, W50 BR1 |
| ST25 | A7436, A7436-C, 5A85Pg |
| ST1 | A7A1-28 (×2) |
| ST30 | W4087, WW2952, WW5019 |

ATCC 33277, DSM 20709 and NBRC 115150 are the same strain under three culture
collection numbers, and 381 is its close relative. The W50/W83 group and the
A7436 derivatives behave the same way. ST25 differs from ST24 at *hagB* alone,
consistent with the known relatedness of A7436 and W83.

## Sweep results — 109 RefSeq assemblies

| Outcome | n | % |
|---|---|---|
| Defined ST | 24 | 22 |
| Novel allele combination | 77 | 71 |
| Incomplete (locus not recovered) | 8 | 7 |

Ten distinct STs among the 24: ST24 (6), ST96 (5), ST25 (3), ST30 (3), ST1 (2),
ST14, ST29, ST67, ST69, ST160.

Novel alleles per locus: ftsQ 21, gpdxJ 20, pga 17, pepO 12, hagB 8, mcmA 8,
recA 6.

All eight incomplete calls failed at *hagB* and nothing else. *hagB* sits in a
family of near-identical haemagglutinin genes that short-read assemblers
routinely collapse, so this is an expected and locus-specific artefact rather
than a threshold problem. Affected: 13_1, AFR5B1, F0570, H3, KCOM 2797,
NBRC 115149, SJD2, SJD5.

The 71 % novel rate is the headline finding. It reflects how thinly the scheme
has been sampled — roughly 190 isolates deposited in eighteen years — rather than any
deficiency in the loci. The user's prior intuition that *P. gingivalis* MLST is
"fragmented and does not give conclusive relationships" is borne out by the
numbers (177 STs over ~190 isolates, no clonal complexes), but the cause is
under-deposition, and running gingiCHARM across a genome collection is a cheap
way to correct it.

## Validation 3 — phylogeny and trait association

Concatenated seven-locus alignment: 109 taxa × 2,800 bp, 139 variable sites
(5.0 %), 84 unique sequences. Eight taxa carry 380 Ns at the *hagB* block.
Neighbour-joining on p-distances with missing positions excluded pairwise;
midpoint-rooted.

Phylogeny–trait association by Fitch parsimony, null from 1,000 tip-label
permutations. Untypeable strains excluded per trait.

| Trait | n typed | states | observed changes | null mean | P | CI |
|---|---|---|---|---|---|---|
| fimA | 109 | 5 | 34 | 49.0 | 0.001 | 0.12 |
| mfa1–mfa4 | 109 | 4 | 26 | 34.3 | 0.001 | 0.12 |
| mfa5 | 109 | 8 | 42 | 54.3 | 0.001 | 0.17 |
| rag | 89 | 4 | 22 | 42.4 | 0.001 | 0.14 |
| K-antigen | 105 | 5 | 36 | 54.3 | 0.001 | 0.11 |

CI = consistency index = (states − 1) / observed changes.

Both halves of this result matter and they pull in opposite directions. Every
trait is significantly clustered on the tree (P = 0.001 is the resolution floor
of 1,000 permutations, so the true P is lower), which supports the reviewer's
point that virulence alleles map onto lineages. But consistency indices of
0.11–0.14 mean each type has been gained, lost or exchanged ten or more times
across the tree, which supports the user's point that the relationships are not
clean. *rag* shows the strongest signal relative to its null, as expected for a
horizontally acquired island that is then inherited vertically.

The defensible conclusion is that an ST constrains but does not determine the
virulence-locus profile, and neither typing framework substitutes for the other.

## Files

- `pg_mlst_calls.tsv` — ST and seven-locus profile for all 109 assemblies
- `pg_typing_calls.tsv` — fimA / mfa / rag / K-antigen for the same 109
- `pg_mlst_concat.fasta` — 2,800 bp concatenated alignment
- `pg_mlst_nj.newick` — neighbour-joining tree
- `phylo_trait_association.json` — parsimony test output
- `Fig_MLST_tree_heatmap.svg` / `.png` — tree with virulence-type heatmap

## Known limitations

- No *T. denticola* or *T. forsythia* scheme exists; the module says so
  explicitly rather than returning an empty result.
- *hagB* recovery depends on assembly contiguity. Long-read assemblies resolve
  it; short-read drafts may not.
- The scheme is nucleotide-only and cannot be run from a protein FASTA.
- The bundled allele set is the public PubMLST release (data submitted on or
  before 31 December 2024). Five sequence types added since are not included.
- Novel alleles are reported but not assigned numbers. Numbering requires
  submission to PubMLST, which the module does not automate.

---

# October 2026 release: database update and re-validation

Eleven alleles accepted on 30 September 2026 and ten sequence types assigned
after them were added to the bundled database on 2 October 2026. The scripts
that produced every number below are in `validation/`.

## What changed

| | Before | After |
|---|---|---|
| Alleles, ftsQ/gpdxJ/hagB/mcmA/pepO/pga/recA | 42/39/38/31/37/28/15 | 42/42/41/31/40/29/16 |
| Alleles, total | 230 | 241 |
| Sequence types defined in the scheme | 200 | 220 |
| Sequence types bundled | 200 | 210 |

The eleven alleles were assigned gpdxJ 40 to 42, hagB 39 to 41, pepO 38 to 40,
pga 29 and recA 16. The ten sequence types were assigned as below, and each is
the profile of one RefSeq assembly.

| ST | Strain | Profile |
|---|---|---|
| ST208 | 1439 | 14-38-1-21-8-3-2 |
| ST209 | KCOM_2799 | 4-39-1-3-6-3-2 |
| ST210 | Kyudai-3 | 41-38-1-6-8-7-2 |
| ST211 | Kyudai-4 | 4-1-38-31-5-7-2 |
| ST212 | SJD12 | 30-26-3-5-29-7-15 |
| ST213 | WW2866 | 21-8-3-3-21-28-14 |
| ST217 | JKG10 | 30-26-41-1-29-2-11 |
| ST218 | LyEC01 | 40-35-12-5-38-7-2 |
| ST219 | NBRC_115148 | 24-13-13-26-40-3-15 |
| ST220 | WW2096 | 1-8-40-3-8-17-4 |

Ten sequence types defined in the scheme are deliberately absent from the
bundled table: ST201 to ST207 and ST214 to ST216 are the PRJDB2925 profiles,
whose designations had not been confirmed when the release was cut. A genome
matching one of those profiles is reported as complete but undefined rather
than as a new type, and the module now says so in its notes.

## Re-validation

Fifty checks were run by `validation/verify_mlst_database.py` and all passed.

| Check | Result |
|---|---|
| Allele numbering per locus | 1..n contiguous, no gaps, at all seven loci |
| Allele length | every allele matches the locus length (420/380/380/420/400/400/400 bp) |
| Allele composition | ACGT only |
| Allele redundancy | no two alleles of a locus share a sequence |
| Allele total | 241 |
| Profile uniqueness | all 210 profiles distinct |
| Profile integrity | every allele referenced by a profile exists in the allele files |
| Profile gaps | exactly ST201 to ST207 and ST214 to ST216 absent |
| New profiles against the submitted file | all ten identical |
| PubMLST isolate validation | 221/221 re-typed to the deposited profile and ST |
| Cross-package sync | all eight data files byte-identical in gingiCHARM and perioCHARM |

The ten new profiles were then re-derived from the assemblies themselves by
exact allele matching on both strands (`validation/rederive_new_profiles.py`),
independently of the submission files, and all ten returned a single
unambiguous allele at every locus matching the assigned type.

**WW3039.** Each of the seven sequences PubMLST exports for isolate 221 carries
a trailing `=` character, which makes every locus one base longer than the
scheme length and matches no allele. The character is an artefact of the export
and is stripped before matching; with that done, WW3039 types to ST199 as
curated. The validation script strips it explicitly rather than silently, so the
behaviour is visible to anyone re-running it.

## Re-typing the RefSeq assemblies

`validation/retype_refseq_assemblies.py` re-typed all 109 RefSeq assemblies
against the updated scheme and rewrote Supplementary Table S3.

| Outcome | Assemblies | Before |
|---|---|---|
| Defined ST | **59 (54 %)** | 48 (44 %) |
| Undeposited allele at one or more loci | 41 | 44 |
| Locus not recovered end to end | 3 | 8 |
| All alleles deposited, type not yet defined | 6 | 7 |
| Unresolved multi-copy locus | 0 | 3 |

Distinct types rose from 33 to 43, of which 32 were defined by this work. The
multi-copy *hagB* pair was detected in 11 assemblies and resolved against the
profile table in all 11, where previously three were left unresolved; the three
extra deposited alleles at *hagB* are what closed them.

Six assemblies, representing five distinct strains, now carry seven deposited
alleles in a combination that is not yet a defined type: AJW4, KCOM_2796, TDC60
(deposited twice, as TDC60 and TDC_60), WW3040 and WW5127. These are candidates
for a further profile submission.
