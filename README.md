<p align="center">
  <img src="assets/gingicharm_logo.svg" alt="gingiCHARM" width="200">
</p>

<h1 align="center">gingiCHARM</h1>

<p align="center">
  <em>P. <strong>gingi</strong>valis <strong>CH</strong>aracterisation, <strong>A</strong>llele typing and <strong>R</strong>eference <strong>M</strong>odule</em>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.8%2B-blue" alt="Python 3.8+">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="MIT licence">
  <img src="https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey" alt="Linux and macOS">
</p>

---

Characterise a *Porphyromonas gingivalis* genome in one pass — virulence
factors, functional categories, and full fimA / mfa-operon / rag / capsular
K-antigen typing — from a single curated reference database.

Sequencing has got cheap enough that oral anaerobes are whole-genome sequenced
routinely. The bioinformatics has not kept up. *E. coli*, the streptococci and
the staphylococci all have dedicated typing servers, and the general typing
platforms inherit their reference sets from those species. *P. gingivalis* is
in none of them. Typing an isolate today means reading four separate typing
literatures, pulling reference alleles out of papers published between 1998 and
2025, and writing a script per scheme — then repeating it genome by genome.
Most groups skip it.

gingiCHARM makes it one command.

```
fimtyping  my_isolate.fna
# fimA type: I  (100.0% identity, 100.0% coverage of reference)
```

## What it does

| Module | Question it answers |
|---|---|
| `pgvfdb` | which virulence factors does this genome carry? |
| `fimtyping` | what is the fimA fimbrial type (I, Ib, II–V)? |
| `mfatyping` | what is the minor-fimbriae operon profile (mfa1/2/3/4/5)? |
| `ragtyping` | what is the rag type (rag-1 to rag-4)? |
| `kantigen` | what is the capsular K-antigen serotype? |
| `pgfunc` | what functional category does each gene fall into? |
| `pgmlst` | what is the PubMLST seven-locus sequence type? |

Every module reads the same bundled database — 3,001 curated records, plus the
PubMLST allele set used by `pgmlst`. Nothing extra to download.

## Install

Python 3.8 or newer. The first install needs internet access to fetch NCBI
BLAST+, the one external dependency.

```bash
git clone https://github.com/S-Alapati/gingicharm.git
cd gingicharm
bash install.sh
```

This creates an isolated environment, installs the package, downloads BLAST+
and puts the commands on your PATH. Check it worked:

```bash
pgtoolkit-info
```

## Use it

Every command takes one FASTA — assembly, contigs, CDS or predicted proteins:

```bash
fimtyping  my_strain.fna
mfatyping  my_strain.fna
ragtyping  my_strain.fna
kantigen   my_strain.fna
pgvfdb     my_strain.fna -o virulence.tsv
pgfunc     my_strain.fna -o categories.tsv
```

Add `-o file.tsv` for the full table. `--min-identity`, `--min-coverage` and
`--threads` are there if you need them.

From Python:

```python
import mfatyping

result = mfatyping.analyze("my_strain.fna")
print(result.operon_profile)     # ('70A', '70', '70', '70', 'A1')
result.to_tsv("mfa_report.tsv")
```

### Browser interface

```bash
cd web
bash start.sh
```

Opens on <http://127.0.0.1:5057>. Upload one or many genomes, tick the modules
you want, browse the results and download them as TSV.

## What is in the database

Built from five complete *P. gingivalis* genomes — ATCC 33277, W83, TDC60,
A7436 and AJW4 — re-annotated with Bakta and clustered with MMseqs2 at 80 %
identity over 80 % coverage into 2,616 non-redundant gene families, 1,476 of
which are present in all five strains. To that core we added 154 curated
virulence-factor families and the reference alleles for all four typing
schemes.

| Scheme | Types covered | Not yet covered |
|---|---|---|
| fimA | I, Ib, II, III, IV, V | — |
| mfa operon | mfa1: 53 / 70A / 70B · mfa2–4: 53 / 70 · mfa5: A1 A2 B C D E | — |
| rag | rag-1, rag-2, rag-3, rag-4 | — |
| K-antigen | K1, K3, K4, K6, K− | K2, K5, K7 |
| MLST | 220 sequence types, 241 alleles | ST201–207, ST214–216 not yet bundled |

K2, K5 and K7 are missing because no capsular sequence has ever been deposited
for their reference strains (HG184, HG1690, 34-4). GenBank holds exactly two
serotype-annotated *P. gingivalis* capsular records, AJ969093 (strain 381, K−)
and AJ969094 (HG1703, K4), and the K-antigen references are built from the
locus those two anchor.

The MLST module carries the PubMLST seven-locus scheme. Thirty-two sequence
types and seven alleles recovered from public RefSeq assemblies were returned
to PubMLST and defined there, taking the scheme from 177 types to 220.

Every gene also carries a functional category from a ten-class scheme curated
for this organism. GO and KEGG annotate fewer than half the *P. gingivalis*
genome, so an enrichment analysis run against them discards most of the genome
before it starts; these categories cover all of it.

## How typing works

Sequence comparison, no PCR step. Each module BLASTs the query against the
reference alleles for its scheme and reports the closest match, ranked by
**identity × coverage** rather than raw bitscore — reference alleles within a
scheme differ a lot in length (the mfa5 genotypes span 3,681 to 6,615 bp) and
bitscore ranking would systematically favour the longest one.

An mfa component with no qualifying hit is reported as `X`, meaning unknown or
missing rather than confirmed deletion: the locus may be genuinely absent, too
divergent for the current panel, or in an incomplete part of the assembly.

## Can you trust it

Validated against twelve genomes whose types are documented in the literature —
every scheme returned the expected call. Functional categories were
cross-checked against curated Swiss-Prot homologs, which found and fixed 22
mislabels. The MLST caller reproduces the curated allele designations and
sequence type for all 221 isolate records held at PubMLST.

Checking the panels against the published reference sequences, rather than
against our own output, found three faults that are now fixed. They are
recorded here because each was invisible from the results alone.

The K-antigen panel was built on the wrong locus. Its references sat 200 kb to
2.15 Mb from the capsular locus in every strain tested, and the 390 gene-level
entries were largely housekeeping genes. It nonetheless returned the published
serotype for its own reference strains, because each had been cut from a strain
of known serotype and the region is strain-variable, so the module was matching
strains rather than determining serotypes. The old references covered one
another at 92–105 % across every serotype pair, leaving nothing to discriminate
on; the rebuilt ones cover one another at 40–68 %. Re-typing changed 71 of 104
calls.

The rag panel labelled A7436 as rag-3 when it is rag-1 (99.8 % to the deposited
ragA-1 allele at full length), so there was no rag-3 reference at all. A
genuine one was added from the allele-designated GenBank records.

fimA type Ib cannot be called by sequence identity. Nakagawa et al. 2002 define
it by a diagnostic RsaI site, not by similarity, and their assay is reproduced
in silico here. Across 104 genomes it finds 15 carrying the site; whole-gene
best match recovers only 6 of those and assigns 5 of the rest to type I.

The package also ships a prebuilt BLAST index that edits to the sequence files
do not regenerate, and a stale index returns nothing while remaining internally
valid. `validation/check_db_index_sync.py` compares the two and fails if they
have drifted.

## Limitations

- K2, K5 and K7 capsular references are unavailable (see above); strains of
  those serotypes are reported untypeable rather than mis-assigned.
- Strains carrying tandem mfa5-1 + mfa5-2 copies get the single best mfa5 hit,
  not both. An assembly-aware extension is planned.
- The rag locus is divergent enough that annotation misses *ragA* in the rag-4
  group; those calls rest on *ragB* and are flagged as ragB-supported.
- Ten sequence types defined at PubMLST (ST201–207, ST214–216) are not in the
  bundled profile table; a genome matching one is reported as complete but
  undefined rather than as a new type.
- The K6 capsular reference comes from an assembly whose strain attribution is
  unresolved. The sequence is the only K6 reference in existence, so it is kept
  and cited by accession rather than by strain name.
- Calibrated for *P. gingivalis*. Confirm species identity before typing —
  a *P. gulae* genome will produce confusing calls.

## Repository layout

```
gingicharm/
├── install.sh, install.py    one-step installer
├── pyproject.toml
├── src/                      the Python package
│   ├── pgcore/               shared BLAST layer + bundled database
│   ├── pgvfdb/  fimtyping/  mfatyping/
│   ├── ragtyping/ kantigen/ pgfunc/
│   └── pgmlst/               PubMLST scheme, bundled separately
├── validation/               database and typing validation scripts
├── web/                      Flask browser interface
└── assets/                   logo and figures
```

The package is distributed under the module name `pgtoolkit`; gingiCHARM is
the product name. The console commands and Python imports use the module
names shown above.

## Citing

If gingiCHARM is useful in your work, please cite it. A `CITATION.cff` is
included, and GitHub will render a "Cite this repository" button from it.
The accompanying paper is in preparation — this section will be updated with
the citation once it is published.

## Licence

MIT — see [LICENSE](LICENSE).
