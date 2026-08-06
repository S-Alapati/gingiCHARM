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
2023, and writing a script per scheme — then repeating it genome by genome.
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
| `fimtyping` | what is the fimA fimbrial type (I–V)? |
| `mfatyping` | what is the minor-fimbriae operon profile (mfa1/2/3/4/5)? |
| `ragtyping` | what is the rag type (rag-1 to rag-4)? |
| `kantigen` | what is the capsular K-antigen serotype? |
| `pgfunc` | what functional category does each gene fall into? |

Every module reads the same bundled database — 3,388 curated records. Nothing
extra to download.

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
| fimA | I, II, III, IV, V | — |
| mfa operon | mfa1: 53 / 70A / 70B · mfa2–4: 53 / 70 · mfa5: A1 A2 B C D E | — |
| rag | rag-1, rag-2, rag-3, rag-4 | — |
| K-antigen | K1, K3, K4, K6, K− | K2, K5, K7 |

K2, K5 and K7 are missing because no complete-genome assembly has ever been
deposited for their reference strains (HG184, HG1690, 34-4) — only short *cps*
fragments from the original typing papers.

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
every scheme returned the expected call. For fimA we also ran an independent
in-silico PCR with the published typing primers; those calls matched the
sequence-based calls for all twelve. Functional categories were cross-checked
against curated Swiss-Prot homologs, which found and fixed 22 mislabels. Every
NCBI accession in the reference set was re-verified through Entrez, which
caught two mislabels during construction (ATCC 49417 capsular serotype and rag
allele).

## Limitations

- K2, K5 and K7 capsular references are unavailable (see above); strains of
  those serotypes are reported untypeable rather than mis-assigned.
- Strains carrying tandem mfa5-1 + mfa5-2 copies get the single best mfa5 hit,
  not both. An assembly-aware extension is planned.
- The rag locus is divergent enough that annotation misses *ragA* in the rag-4
  group; those calls rest on *ragB* and are flagged as ragB-supported.
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
│   └── ragtyping/ kantigen/ pgfunc/
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
