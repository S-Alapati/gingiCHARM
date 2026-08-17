# MLST correction: multi-copy loci and tied alleles

Raised by review of the tree/heatmap figure: several isolates that carry a
sequence type in PubMLST were reported as novel. Investigated 14 August 2026.
One real bug found and fixed; two other classes of discrepancy characterised
and shown not to be calling errors.

## What was checked

All 192 *P. gingivalis* isolate records were pulled from the PubMLST isolates
database with their seven-locus designations, and every genome in the
109-assembly set whose name could be matched to an isolate record was compared
locus by locus. Thirteen genomes matched. (PubMLST isolate names are mostly
clinical collection identifiers — SA1–46, ME/*n*-*n*, L*n*A*nn*, 4###-C — that
do not correspond to NCBI assembly names, so the small overlap is real and not
a matching failure.)

An earlier attempt at this comparison reported that only 36 of 192 isolates had
public allele designations. That was wrong: the bulk fetch was silently
swallowing HTTP failures and returning empty designations. With retries, all
192 have full seven-locus designations. The conclusion drawn from the broken
fetch — that the Spanish isolates had no public STs — was therefore also wrong.

## Result

| Strain | PubMLST | Called | Per-locus |
|---|---|---|---|
| ATCC 33277 (×2) | 5-9-1-1-1-5-5 | 5-9-1-1-1-5-5 | ======= |
| 5A85Pg | 16-6-2-4-1-3-5 | 16-6-2-4-1-3-5 | ======= |
| A7A1-28 (×2) | 1-12-1-1-1-1-1 | 1-12-1-1-1-1-1 | ======= |
| B42 | 24-22-1-17-8-7-2 | 24-22-1-17-8-7-2 | ======= |
| Pg152 | 19-10-18-1-16-2-1 | 19-10-18-1-16-2-1 | ======= |
| 144_3 | 5-5-1-4-16-16-4 | 5-**~5**-1-4-16-16-4 | =~===== |
| 153Pg | 14-5-1-3-8-2-3 | 14-5-1-**~1/3**-8-2-3 | ===~=== |
| P4_Pig1 | 19-1-1-12-28-3-3 | **~19**-1-1-12-28-3-3 | ~====== |
| UB14Pg | 17-8-17-5-21-17-1 | 17-8-**12/17**-5-21-**~17**-1 | ==+==~= |
| AJW4 | 16-9-22-17-22-15-1 | 21-23-1/6-3-20-3-7 | XXXXXXX* |
| HG1691old | 3-6-35-28-12-3-3 | 14-1-1/6-6-8-~23-2 | XXXXXXX* |

*Compared against the isolate record. For HG1691old, a scan of the assembly itself agrees with gingiCHARM at all seven loci — see below.

`=` exact · `+` designated allele is among the candidates · `~` designated
allele is the nearest, one base away in the assembly · `X` disagree

Seven of thirteen agree exactly. The remaining six fall into three groups.

## Bug 1 — multi-copy loci resolved arbitrarily (fixed)

*hagB* lies in a tandem pair of near-identical haemagglutinin genes about
3.1 kb apart. **48 of 109 assemblies carry both copies, and 38 give a
full-length 100 % match at both.** In most strains both copies carry the same
allele and the distinction is invisible, but in 12 genomes they differ.

The old code scored hits as `(exact, identity, length)` and took the maximum,
which is undefined when two hits tie at 100 %. It returned whichever BLAST
happened to emit first. For UB14Pg that produced *hagB* 12 where PubMLST
designates 17 — and both are genuinely present in the assembly, at
2,126,683 and 2,123,533 respectively.

No positional rule fixes this. In UB14Pg the correct copy is the
lower-coordinate one; in Pg152 it is the higher. Assemblies also differ in
origin and orientation, so absolute coordinates are not comparable.

The fix returns **every** candidate allele rather than one, and resolves the
ambiguity where it can be resolved: the seven-locus combinations the candidates
allow are enumerated and looked up in the profile table, and if exactly one is
a defined ST that ST is reported with a note naming the multi-copy locus. If
none or several match, all candidates are reported and no ST is claimed.

Verified no regression: all 24 previously defined STs are unchanged, and
W83, W50, ATCC 33277, A7436, HG66, A7A1-28, B42 and Pg152 still reproduce
their deposited profiles exactly.

## Bug 2 — tied nearest alleles reported arbitrarily (fixed)

When a novel allele is equidistant from several deposited alleles, the old code
named one of them arbitrarily. 153Pg *mcmA* is 1 SNP from both allele 1 and
allele 3; the old output named 1, PubMLST designates 3. Both are now reported
(`~1/3`), which shows the designated allele is among the nearest.

Affected genomes per locus: *hagB* 12, ftsQ 6, mcmA 5, pepO 5, gpdxJ 4, pga 3,
recA 3.

## Not a bug — single-base assembly discrepancies

For 144_3, 153Pg, P4_Pig1 and UB14Pg, exactly one locus sits exactly one base
from the allele PubMLST designates for that same isolate. These are not calling
errors. Searching the deposited allele against the whole assembly with no
identity floor and word size 11 finds it exactly once, at 99.74–99.76 %, and
never at 100 %. The assembly genuinely carries the variant base.

All five differences are transitions, four of them C→T or G→A, which is the
signature of cytosine deamination and a common library-preparation artefact.
Two explanations remain and cannot be separated without read data or curator
input: the assembly has a base-calling error at that position, or the PubMLST
designation was assigned from separate template. Assembly level does not
predict novelty overall (51 % of complete genomes carry a novel allele versus
47–63 % of drafts), so this is not a general draft-quality effect.

**Consequence for the submission set.** Four candidate novel alleles are
withdrawn — 144_3 *gpdxJ*, 153Pg *mcmA*, P4_Pig1 *ftsQ*, UB14Pg *pga*. The
isolate is already typed in PubMLST and the deposited allele is one base away,
so a new allele should not be claimed.

## Independent validation against the PubMLST genome scanner

PubMLST's own scanner was run on the HG1691old assembly (NZ_CP116613.1) and its
output compared with gingiCHARM's, locus by locus and coordinate by coordinate:

| Locus | PubMLST scanner | gingiCHARM | Agree |
|---|---|---|---|
| ftsQ | 14, 1795047–1795466 | 14, same coordinates | yes |
| gpdxJ | 1, 1750464–1750843 | 1, same coordinates | yes |
| hagB | **6**, 2272663–2273042 **and 1**, 2269513–2269892 | 1/6, same coordinates | yes, both copies |
| mcmA | 6, 461492–461911 | 6, same coordinates | yes |
| pepO | 8, 246573–246972 | 8, same coordinates | yes |
| pga | no allele assigned | novel, 99.75 % to allele 23 | yes, neither assigns |
| recA | 2, 1289278–1289677 | 2, same coordinates | yes |

Agreement is exact at all seven loci, to the base. Three things follow.

The *hagB* duplication is confirmed by the curators' own tool, which likewise
reports two copies carrying different alleles at the two positions this analysis
identified. Reporting both candidates rather than silently choosing one is
therefore the behaviour that matches curator practice, not a workaround.

The *pga* call is corroborated. PubMLST's scanner assigns no allele at that
locus, which is what it does when no deposited allele matches exactly. That is
independent evidence the allele is genuinely new rather than an artefact, and it
promotes this candidate — *pga*, 1 SNP from allele 23, otherwise a single-strain
singleton in the "verify" tier — to a supportable submission.

## Correction to an earlier conclusion — HG1691

An earlier version of this document listed HG1691old as a case where
gingiCHARM disagreed with PubMLST at every locus. That was the wrong framing.
gingiCHARM does not disagree with PubMLST: run on the same sequence, the two
tools agree exactly, as above.

What differs is the *isolate record*. PubMLST isolate "HG1691" carries the
designation 3-6-35-28-12-3-3, while the genome deposited as HG1691old
(NZ_CP116613.1) types as 14-1-1/6-6-8-(novel *pga*)-2. The record and the
assembly describe different organisms. The evidence for that conclusion stands,
but it is evidence about the two public records, not about the typing method.

The same caution applies to **AJW4**, still listed above as a whole-profile
conflict. That comparison was also made against the isolate record rather than
against a scan of the assembly. The assembly (CP011996.1) differs from the
deposited AJW4 designation at all seven loci by 21 SNPs across 2,800 bp, which
is far more than duplicate deposits of one strain show, so a labelling conflict
between record and assembly remains the most likely reading — but running the
PubMLST scanner on CP011996.1 would settle it the same way it settled HG1691,
and that check has not yet been done.

## Files

- `pg_mlst_calls_v2.tsv` — corrected calls, with `multicopy_loci` column
- `pubmlst_comparison.tsv` — the thirteen-genome comparison above
- `ambiguity_audit.json` — per-genome, per-locus site and tie counts
- `withdrawn_from_submission.json` — the four withdrawn candidates
