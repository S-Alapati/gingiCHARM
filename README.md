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
in none of them. Typing an isolate today means reading several typing
literatures, pulling reference alleles out of papers published between 1998 and
2025, and writing a script per scheme, then repeating it genome by genome.
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
| MLST | 210 of the 220 defined profiles, 241 alleles |-|

K2, K5 and K7 are missing because no capsular sequence has ever been deposited
for their reference strains (HG184, HG1690, 34-4).

The MLST module carries the PubMLST seven-locus scheme. Thirty-two sequence
types and seven alleles recovered from public RefSeq assemblies were returned
to PubMLST and defined there, taking the scheme from 177 types to 220.

## Functional categories

Every gene family carries exactly one of ten categories, curated for this
organism rather than adapted from a general ontology.

![How every gene gets a functional category](assets/functional_categories.png)

### The ten categories

| Category | Scope | Families | With a COG letter |
|---|---|---|---|
| Metabolic Homeostasis | central and intermediary metabolism; energy production; carbohydrate, amino-acid, nucleotide, lipid and cofactor metabolism; metabolite transport; housekeeping enzymes | 755 | 78 % |
| Genetic Information Processing | replication, recombination and repair; transcription and its regulation; translation, including ribosomal proteins, tRNA and rRNA modification and translation factors; RNA processing; signal-transduction regulators acting through DNA binding | 311 | 84 % |
| Cell Envelope & Morphogenesis | biogenesis and maintenance of wall and membranes; peptidoglycan synthesis and remodelling; outer-membrane proteins and porins; lipoproteins; division and elongation machinery | 113 | 35 % |
| Surface Structures & Colonization | surface-exposed structures mediating host attachment: major and minor fimbriae, pili, adhesins, capsular polysaccharide, LPS and O-antigen biosynthesis | 88 | 45 % |
| Secretion & Molecular Export | protein secretion and solute export across the envelope: T9SS/Por, ABC transporters, efflux pumps, dedicated export permeases | 106 | 76 % |
| Proteolytic Virulence & Effectors | secreted and surface proteolytic virulence factors: gingipains, collagenases, haemagglutinins and other extracellular proteases implicated in tissue degradation and immune evasion | 80 | 74 % |
| Heme & Ion Homeostasis | acquisition, transport and storage of iron, haem and metal ions: TonB-dependent receptors, the Hmu system, ferrous and ferric transporters, metal-storage proteins | 55 | 78 % |
| Environmental Stress & Defense | protection against environmental and host-derived stress: oxidative-stress defence, chaperones and heat-shock proteins, restriction–modification and CRISPR–Cas anti-phage systems | 78 | 51 % |
| Mobile Genetic Elements & HGT | mobile and horizontally transferred DNA: transposases and IS elements, integrases, prophage genes, conjugative-transfer machinery | 152 | 10 % |
| Uncharacterized Proteins | no assignable function: hypothetical and DUF proteins, and proteins whose only homologues are themselves uncharacterized | 878 | 14 % |

Category names are stored verbatim in the `Biological_Category` field, so the
American spellings above are the literal data values and are kept as such.

### Why not COG, GO or KEGG

A COG letter is available for 1,292 of the 2,616 families, 49.4 %. The missing
half is the obvious problem; the structure of the gap is the real one. COG
reaches 84 % of genetic information processing and 78 % of metabolism but
10 % of mobile genetic elements, 35 % of the cell envelope and 45 % of surface
structures. That follows from how orthologous groups are built, since COGs are
defined by conservation across distant taxa and the categories that fall out
worst are the lineage-restricted and recently acquired ones. A COG-based
enrichment test asking about surface structures or horizontally acquired
material therefore runs on a tenth to a half of the relevant genes, and the
genes it drops are not a random sample with respect to the question. The
result is biased, not merely underpowered.

### Evidence gathered per family

Four independent lines, for every family: the COG functional-category letter
from the Bakta annotation; the best `blastp` hit against UniProtKB/Swiss-Prot
at e-value ≤ 1e-10, with the hit product name retained; the Bakta product name
and Pfam/InterPro domains; and a curated list of *P. gingivalis* gene-name and
product-name fragments that unambiguously flag a category, for instance
*rgpA*, *rgpB* and *kgp* for proteolytic virulence and *porK*, *porL*, *porM*
and *sov* for secretion.

### How assignment actually resolved

Each family runs through the routes in order and the first to fire decides.
The route that fired is recorded per family in `Category_evidence`, so every
assignment carries its provenance and can be filtered on.

| Route | `Category_evidence` | What it matches | Families |
|---|---|---|---|
| Curated product name | `exact-product` | the Bakta product name against a curated product-to-category table | 2,248 |
| Keyword and Pfam domain | `keyword` | product name and Pfam accession against a conservative rule table, where a term must be functionally definitive to fire: "ribosomal protein" to Genetic Information Processing, "TonB-dependent receptor" to Heme & Ion Homeostasis, `transposase\|IS_family` to Mobile Genetic Elements | 192 |
| Organism-specific marker | `gene-symbol` | gene symbol against the curated *P. gingivalis* marker list | 35 |
| COG letter | `COG-fallback` | the COG letter, where nothing above fired | 16 |
| Default | `default` | nothing fired; the family is placed in Uncharacterized Proteins | 103 |
| Swiss-Prot review | `swissprot-corrected` | reassigned by hand after the verification pass below | 22 |

Two things in that table are worth reading carefully, because they are not
what a four-tier description would lead you to expect.

The curated product-name table does **86 % of the work**. The organism-marker
list, which sounds like the primary mechanism, fires for 35 families, and the
COG letter is a genuine last resort at 16. Most of the curation effort is
therefore in the product-to-category mapping, not in the organism-specific
patterns.

And of the 878 families in Uncharacterized Proteins, **745 got there because
the product name said so** — hypothetical protein, DUF-something — and only
103 through the default route. The category is mostly a positive statement
that the annotation declares the protein unknown, rather than a residue of
classification failure.

### Verification against Swiss-Prot

Every representative protein was searched against UniProtKB/Swiss-Prot with
`blastp` in batches of about 550, and the best hit compared with the assigned
category. Swiss-Prot names beginning UPF, Uncharacterized, DUF or hypothetical,
and names that are only locus tags, were excluded by an informative-name
filter, since such a name can neither confirm nor refute a category.

Where an informative name disagreed, the family was reviewed by hand. Of the
disagreements, 97 proved to be systematic naming differences, for instance
"UvrABC system protein C" against "excinuclease UvrC", and were resolved by
synonym-aware comparison and a gene-symbol substring match. Twenty-two were
real mislabels and were corrected: seven named genes (*panC*, *ompR*, *xthA*,
*ompA*, *porG*, *rarA* and PgCG_01627) and fifteen families that had sat in
Uncharacterized Proteins despite unambiguous Swiss-Prot homology. All twenty-two
carry `Category_evidence = swissprot-corrected`.

### Do the categories mean anything

They can be checked against something they were not built from. The pangenome
partition is computed from presence and absence across the five source genomes
and uses no category information, and the categories line up with it: 595 of
755 metabolic families are core, as central metabolism should be; 143 of 152
mobile element families are accessory, as mobile DNA should be; and 360 of 878
uncharacterized families are strain-unique, which is what an unknown residue
should look like if it holds genuinely unknown material.

### Using them

```bash
pgfunc my_strain.fna -o categories.tsv
```

Each predicted coding sequence inherits the category of its best-scoring
database hit at 70 % identity over 70 % coverage or better, rather than being
re-derived. A query with no qualifying hit is reported explicitly rather than
dropped. A gene family absent from all five source genomes therefore comes
back as uncharacterized by default rather than by evidence, which is worth
remembering when the query is divergent.

## How typing works

Sequence comparison, no PCR step. Each module BLASTs the query against the
reference alleles for its scheme and reports the closest match, ranked by
**identity × coverage** rather than raw bitscore — reference alleles within a
scheme differ a lot in length (the mfa5 genotypes span 3,681 to 6,615 bp) and
bitscore ranking would systematically favour the longest one.

An mfa component with no qualifying hit is reported as `X`, meaning unknown or
missing rather than confirmed deletion: the locus may be genuinely absent, too
divergent for the current panel, or in an incomplete part of the assembly.

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

If gingiCHARM is useful in your work, please cite it.

> Alapati, S. gingiCHARM: *P. gingivalis* characterisation, allele typing and
> reference module. GitHub https://github.com/S-Alapati/gingiCHARM (2026)

The accompanying paper is in preparation, and this section will be updated with
the citation once it is published.

## Licence

MIT — see [LICENSE](LICENSE).
