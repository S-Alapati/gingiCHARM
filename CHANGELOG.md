# Changelog

All dates are the dates the work was done. Where a change depended on an
external record, such as a PubMLST curation decision, the date that record
carries is given too, so the two can be checked against each other.

---

## 2026-10-03 — Panel audit

Each typing panel was checked against the published reference sequences
rather than against the tool's own output. Three faults were found. None was
visible from the results, which is the reason each is recorded here rather
than quietly fixed.

### The K-antigen panel was built on the wrong locus

The thirteen region-level references described as the *cps* neighbourhood sat
200 kb to 2.15 Mb away from the capsular locus in every strain tested.
Locating AJ969093 (strain 381, K−) and AJ969094 (HG1703, K4), the only two
serotype-annotated capsular records in GenBank, places the locus at
241,416–254,299 in ATCC 33277 and 2,274,267–2,277,267 in ATCC 49417, matching
at 99.99 % and 99.93 % over their full lengths. The former references sat at
2,308,574 and 2,475,342 in those same genomes. A direct search of AJ969094
against the reference we called K4 returned nothing above 15 bp.

The gene content confirmed it. Among the 390 gene-level entries the commonest
products were TPR domain proteins, MFS transporters, excinuclease UvrA,
peptide deformylase, Holliday junction resolvase RuvX and the TrkA/TrkH
potassium uptake system. Only 41 of 390 carried any sugar or polysaccharide
term.

The panel nonetheless returned the published serotype for each of its own
reference strains, which is why the fault persisted. Each reference had been
cut from a strain of known serotype and the region is strain-variable, so the
module was matching strains rather than determining serotypes. The difference
is measurable: across every serotype pair the old references covered one
another at 92–105 %, leaving no sequence difference to discriminate on, while
the rebuilt references cover one another at 40–68 %.

Rebuilt as five references, one per serotype, located by homology to AJ969093
and extracted with 1.5 kb flanks: K− (ATCC 33277), K1 (W83), K3 (A7A1-28),
K4 (ATCC 49417) and K6 (GCF_028335085.1). All six reference strains now call
correctly, including A7436 as K1. **Re-typing changed 71 of 104 calls.**

Matching was changed with it. The locus is serotype-variable, so a query
aligns in fragments rather than one block and per-fragment thresholding
rejected every genuine match. Fragments are now aggregated per reference
before the identity and coverage thresholds apply.

### The rag panel had no rag-3 reference

RAGA_2 was labelled rag-3 but is 99.8 % identical to the deposited ragA-1
allele over its full length, and A7436's *ragB* is 99.9 % to ragB-1, so A7436
is rag-1. Relabelled, and a genuine rag-3 added from the allele-designated
GenBank records PQ658524 and PQ658651. Twenty genomes previously untypeable
are now rag-3; none remain untypeable.

### fimA type Ib cannot be called by sequence identity

Nakagawa *et al.* 2002 (*J Periodont Res* 37:425–432) state that type Ib
cannot be distinguished from type I by type-specific PCR, the two being
97.1 % identical, and define it instead by a diagnostic RsaI (AfaI) site.
Their Ib-F/Ib-R pair amplifies types I and Ib alike and only the Ib product
is cut. That assay is now reproduced in silico and is authoritative for the
I/Ib pair.

It gives a 269 bp product from 28 of 104 genomes and none from the other 76,
matching the published behaviour, with Ib products cutting to 160 and 109 bp
against the 162 and 109 reported. Fifteen genomes carry the site. Whole-gene
best match recovers only six of those at a one-point margin, leaves four
inside the margin and assigns the remaining five to type I outright, each
sitting about 0.4 points closer to the type I reference while carrying the
base that defines Ib.

An in-house Ib primer previously held in the digital panel was withdrawn: its
261 bp product in strain CP3 aligns to *fimA* type IV at 96.9 % identity.

### The shipped BLAST index can go stale silently

The package ships a prebuilt BLAST index that edits to the sequence files do
not regenerate. A stale index returns nothing while remaining internally
valid, so nothing detects it; FIMA_6 had never been indexed.
`validation/check_db_index_sync.py` now compares the index against the FASTA
and fails if they have drifted.

### Nearest-profile reporting for MLST

A genome with no exact sequence type was reported as "novel ST" or
"incomplete", which says nothing about where the strain sits. The nearest
defined profile is now reported with the loci matched, the mismatch count,
the mismatching loci and any equally close profiles. A trailing symbol marks
the result as nearest rather than assigned: `*` where the profile is complete
and every allele deposited, `¶` where a novel allele is carried, `§` where a
locus was not recovered end to end.

A PubMLST profile query defaults to exact-or-nearest with no threshold and
offers an optional "N or more matches" filter. The filtered form is used
here, at four of seven loci, because a nearest profile sharing three loci or
fewer invites being read as an assignment. `min_loci_matched=0` restores the
unfiltered behaviour.

Validated against PubMLST's own answer for strain WW2885: it returns ST-111
and two others at 4/7 with 3 mismatches; this returns ST111 with the same
4/7, the same three mismatches at *ftsQ*, *gpdxJ* and *recA*, and the same
two tied profiles, ST120 and ST148. All 210 bundled profiles return
themselves at 7/7.

Of the 109 RefSeq assemblies, 28 of the 50 without an exact type now reach
the threshold, ten of them matching a named type at six of seven loci.

### Database

3,399 → 3,001 records. fimA calls unchanged, 109/109.

---

## 2026-10-02 — Fourth PubMLST contribution, and fimA type Ib

Eleven alleles accepted by the curators on 30 September 2026 (*gpdxJ* 40–42,
*hagB* 39–41, *pepO* 38–40, *pga* 29, *recA* 16) and ten sequence types
assigned on 2 October (ST208–ST213, ST217–ST220) were added to the bundled
scheme. All ten profiles were re-derived from their assemblies by exact
matching on both strands, independently of the submission files, and each
returned a single unambiguous allele at every locus.

ST201–ST207 and ST214–ST216 are defined in the scheme but deliberately not
bundled: their designations went to profiles from a separate clinical
collection and were not reported back. A genome matching one is reported as
complete but undefined rather than as a new type.

*fimA* type Ib was added as FIMA_6 from AB058848.1, strain HG1691. The coding
sequence annotated on that record begins at an internal GTG and encodes the
mature protein; aligning the type I coding sequence onto it places the
initiator ATG at position 107, giving 1,152 bp that start ATG, end TAA and
encode a 383-residue protein with the expected lipoprotein leader. Positions
107–1258 are used.

Two records from ATCC 33277, PgCG_02243 and VFA_055, carried the product name
"FimA type Ib". UniProt B2RH54 names that protein FimA type-1 and the typing
panel calls the same sequence type I, so both were corrected.

Fifty validation checks on the bundled database, all passing, including the
221-isolate PubMLST validation and byte-identity with the perioCHARM copy.

---

## 2026-09-29 — Second and third PubMLST contributions bundled

Scheme bundled at 200 sequence types and 230 alleles. Seven alleles recovered
from the RefSeq assemblies had been accepted: *ftsQ* 41, *gpdxJ* 38 and 39,
*hagB* 38, *mcmA* 31, *pga* 28 and *recA* 15.

*recA* 15 is worth recording. COM30Pg1.1 had been deposited as *recA* 10,
whereas its assembly carries a variant one base away, and the same variant
occurs in two further RefSeq assemblies with distinct seven-locus profiles,
each sequenced to more than 200x. On being shown this the curators uploaded
the assembly, extracted the alleles directly and retyped the isolate as
ST200. Recurrence across independently sequenced strains distinguished a
genuine allele from a base-calling artefact, which one assembly could not
have done.

Alleles were submitted only when carried by at least two strains with
distinct seven-locus profiles and present in at least one assembly sequenced
to 100x or more. Thirteen single-strain candidates above 200x coverage were
held back on that rule.

---

## 2026-09-18 — September PubMLST release

Bundled scheme updated to the September 2026 release.

---

## 2026-08-17 — First PubMLST contribution

Twenty-two seven-locus profiles recovered from public RefSeq assemblies, each
consisting entirely of deposited alleles in a combination not yet defined,
were submitted and accepted as **ST178 to ST199**, taking the scheme from 177
types to 199. Bundled, and the isolate count in the module notes corrected to
199.

---

## 2026-08-14 — MLST module

`pgmlst` added: PubMLST seven-locus sequence typing (*ftsQ*, *gpdxJ*, *hagB*,
*mcmA*, *pepO*, *pga*, *recA*), with the allele files and profile table
redistributed separately from the unified record database.

Four locus-call outcomes are kept apart because they differ in what they
imply: a deposited allele, a new allele reported against its nearest
neighbour, a locus not recovered end to end, which is an assembly limitation
rather than a biological result, and a multi-copy locus with differing
copies.

The last case is not hypothetical. *hagB* lies in a tandem pair about 3.1 kb
apart in 48 of 109 assemblies, and the copies match different deposited
alleles in 11 of them. No positional rule selects the right copy, since the
designated copy is the lower-coordinate one in some strains and the higher in
others, so every candidate is carried forward and resolved against the
profile table. Fixing this also fixed tied-allele resolution.

---

## 2026-08-06 — v1.0.0

First release. One curated reference database behind six modules: virulence
factor screening, *fimA* typing, the mfa operon, *rag*, the capsular
K-antigen, and a ten-category functional classification.

Built from five complete *P. gingivalis* genomes (ATCC 33277, W83, TDC60,
A7436, AJW4), re-annotated with Bakta under a single parameter set and
clustered with MMseqs2 at 80 % identity over 80 % coverage into 2,616
non-redundant gene families, 1,476 of them present in all five strains.

Zenodo deposition metadata and citation file added.

---

## Known at the time of writing

- K2, K5 and K7 cannot be called. No capsular sequence has ever been
  deposited for their reference strains (HG184, HG1690, 34-4). Strains of
  those serotypes are reported untypeable rather than mis-assigned.
- The K6 reference comes from GCF_028335085.1, an assembly carrying no strain
  designation at NCBI and disagreeing with the deposited HG1691 isolate
  record at all seven MLST loci. It is the only K6 sequence in existence, so
  it is kept and cited by accession rather than by strain name.
- Ten sequence types defined at PubMLST are not in the bundled profile table
  (see 2026-10-02).
- Strains carrying a tandem *mfa5-1*/*mfa5-2* pair get the best-scoring hit,
  not both copies.
- *rag-4* calls rest on *ragB*, since standard annotation does not recognise
  *ragA* in that group, and are flagged as ragB-supported.
- MFA1_2 and MFA1_3 are byte-identical duplicates, and one mfa5 reference is
  labelled "X-Ando", which is not a published genotype name. Neither is wrong
  in the way the K-antigen panel was, but neither is clean.
