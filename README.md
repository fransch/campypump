<img src="assets/campypump_logo.png" alt="campypump" width="150" align="right"/>

# campypump

Classify the cmeB gene of the *Campylobacter* **cmeABC** efflux operon as **RE-cmeB** — the
resistance-enhancing variant that raises fluoroquinolone and macrolide MICs — or as **wildtype**,
directly from genome assemblies.

```bash
campypump run assemblies/ -t 8 -o results.tsv
```

It replaces a manual BLAST-and-align workflow and reproduces its identity values exactly: across
1,513 *C. coli* genomes the calls agreed with the manual method for every genome, and 1,497 of 1,508
complete genes matched to within 0.01 percentage points.

| | |
|---|---|
| Input | genome assemblies (FASTA), one file per isolate |
| Output | one row per genome: cmeB type, identity to reference alleles, operon context, promoter, QC flags |
| Requires | Python ≥ 3.10, Biopython ≥ 1.80, BLAST+ ≥ 2.10 |
| Speed | ~1.7 CPU-seconds per genome (3,736 genomes in 13 minutes on 8 CPUs) |
| Tested on | Python 3.10–3.12 |

---

## Contents

- [Installation](#installation) · [Quick start](#quick-start) · [How it works](#how-it-works)
- [Commands](#commands) · [Output](#output) · [Interpreting the results](#interpreting-the-results)
- [Optional extras](#optional-extras) · [Reference panel](#reference-panel)
- [Troubleshooting](#troubleshooting) · [Documentation](#documentation) · [Citing](#citing)

---

## Installation

### With conda or micromamba (recommended)

```bash
micromamba create -n campypump -c conda-forge -c bioconda \
    python=3.11 biopython blast pip pytest -y
micromamba activate campypump

git clone https://github.com/fransch/campypump.git
cd campypump
pip install .
```

Substitute `conda` or `mamba` for `micromamba` if that is what your site provides. Once the bioconda
recipe is accepted this becomes `conda install -c conda-forge -c bioconda campypump`.

### Optional dependencies

| Package | Needed for |
|---|---|
| `pandas`, `openpyxl` | reading isolate id lists from Excel (`--ids file.xlsx`) |
| `bigsdb-downloader` | PubMLST records that are not publicly released |
| `skani` | whole-genome ANI species check; campypump runs it for you |
| `mafft`, `iqtree` | allele alignment and phylogeny |

None are required for typing.

### Check the installation

```bash
bash verify_install.sh
```

This prints software versions, runs the 60 tests, types the eight bundled reference genomes and
compares the output with expected values — three of them identity values from the manual workflow
campypump replaces. It ends with `ALL CHECKS PASSED` or names the mismatch. Run it before using
campypump on real data.

---

## Quick start

```bash
# type a folder of assemblies
campypump run assemblies/ -t 8 -o results.tsv \
    --alleles cmeB_alleles.fasta --segment-plots plots/

# how many of each type
cut -f5 results.tsv | sort | uniq -c

# what each column means
cat results_column_guide.txt

# check the classes separate in your data
campypump figure results.tsv -o figures/
```

Accepted extensions are `.fasta`, `.fas`, `.fa` and `.fna`; assembled contigs, not reads. Files
named `id-<number>.fasta` get that number carried into a `pubmlst_id` column.

On a cluster, see [`docs/hpc_quickstart.md`](docs/hpc_quickstart.md) and the templates in `hpc/`.

---

## How it works

First the operon is located: all four genes (cmeR, cmeA, cmeB, cmeC) are found with sensitive BLAST
settings — the default megablast misses divergent *C. coli* cmeR — and each is extracted as a
complete ORF, anchored on the reference start codon and read to the first in-frame stop.

Then cmeB is compared with a bundled panel of reference alleles using a global alignment scored like
Geneious Prime's default (match 5, mismatch −4, gap open 12, gap extend 3), with identity computed
over the alignment including gaps. That definition was chosen because it reproduces values from the
manual workflow exactly, so old and new results stay comparable.

The gene is then classified: RE-cmeB at ≥ 95 % identity to the RE reference, wildtype at ≤ 85 %,
ambiguous in between. A sliding window along the gene detects mosaic alleles, where some segments
are RE-like and others wildtype-like, and separate flags mark genes broken across contigs (partial,
with a provisional call) and genes unlike every reference (divergent).

Finally the rest of the operon is reported — per-gene identities, and the sequence of the
CmeR-binding inverted repeat in the cmeR–cmeA intergenic region — along with plain-language QC
warnings.

In the collection campypump was developed on, the two classes were separated by a wide empty gap:
RE alleles were 97.4–100 % identical to the RE reference and wildtype alleles ≤ 80.2 %, with nothing
in between. A maximum-likelihood tree of the alleles recovered the RE set as a single clade with
full support, so the cut-off corresponds to a real evolutionary group rather than an arbitrary
threshold. Full method in [`docs/technical_description.md`](docs/technical_description.md).

---

## Commands

### `campypump run` — type assemblies

```bash
campypump run assemblies/ -t 8 -o results.tsv
```

| Option | Default | What it does |
|---|---|---|
| `-o, --output` | `campypump_results.tsv` | results table |
| `-t, --threads` | 1 | genomes typed in parallel |
| `--species` | — | label written to a `species` column, e.g. `"C. coli"` |
| `--ani` | — | use this skani table instead of running skani |
| `--no-ani` | off | skip the whole-genome species call |
| `--alleles FILE` | — | write the extracted cmeB sequences (needed for diversity figures and trees) |
| `--segment-plots DIR` | — | SVG identity profiles for mosaic, divergent and mixed alleles |
| `--plot-all` | off | with `--segment-plots`, plot every genome |
| `--refs FILE` | bundled panel | use your own reference alleles |
| `--re-min` | 95 | min % identity to the RE reference to call RE |
| `--wt-max-re` | 85 | max % identity to the RE reference to call wildtype |
| `--window`, `--step` | 300, 50 | mosaic scan window and step, bp |
| `--mosaic-frac` | 0.10 | fraction of RE-like *and* wildtype-like windows needed to flag mosaic |

Accepts files or directories, and several at once. Every run also writes
`<output>_column_guide.txt`.

### `campypump figure` — classification figure

```bash
campypump figure results.tsv -o figures/
```

Two panels per species: identity to the RE reference against identity to the wildtype reference for
each genome, and the distribution of identity to the RE reference. This is how you check the classes
separate in *your* collection. `--exclude-duplicates` skips rows marked `use_for_prevalence = no`;
`--include-partial` adds genes broken across contigs, whose whole-gene identity is deflated;
`--combined` adds an all-species figure.

### `campypump diversity` — allele diversity figure

```bash
campypump diversity cmeB_alleles.fasta -o figures/diversity.svg --table figures/alleles.tsv
```

Three panels: alleles ranked by how many genomes carry them, a rarefaction curve showing whether
sampling is saturated, and distinct alleles per cmeB type. `--table` writes one row per distinct
allele. Pass several FASTA files to compare species; `--labels` names them.

### `campypump combine` — stack result tables

```bash
campypump combine ccoli.tsv cjejuni.tsv -o all.tsv
```

### `campypump annotate` — build the analysis file

```bash
campypump annotate results.tsv --metadata pubmlst_metadata.tsv \
    --ani species_ani.tsv --exclude exclude_ids.txt -o annotated.csv
```

Merges PubMLST metadata, a whole-genome ANI species call and an exclusion list, and flags duplicate
records by isolate name. `--columns core` (default) writes the 49 columns analyses need; `--columns
full` writes every column the run produced (95 in the Peru run, with PubMLST metadata and ANI). Every column is documented in
[`docs/codebook_annotated_csv.md`](docs/codebook_annotated_csv.md).

### `campypump alleles` — prepare tree input

```bash
campypump alleles cmeB_alleles.fasta --results annotated.csv -o tree/cmeB_alleles
```

Collapses genomes to distinct alleles with tree-safe names (`Cc_a003_RE_n63`: species, rank by
frequency, type, genomes carrying it) and writes an annotation table, iTOL colour strips and a
Microreact CSV. `--references` adds published alleles, `--outgroup` a root.

### `campypump treeplot` — annotated tree figure

```bash
campypump treeplot cmeB_tree.treefile --annotation tree/cmeB_alleles.tsv -o figures/tree.svg
```

Branches are coloured by cmeB type, or by species inside the RE clade where every tip is RE. With
`--annotation` it adds the dominant ST per allele, scales tips by genome count and draws a source
composition bar. `--subtree RE` draws only the RE clade; `--cladogram` and `--drop-outgroup` help
when one clade is far more divergent than the rest; `--min-support`, `--labels`, `--width`,
`--tip-height` and `--colour-by` adjust the look.

### `campypump fetch` and `campypump metadata` — PubMLST

```bash
campypump fetch --ids ids.txt --outdir assemblies/
campypump metadata --ids ids.txt --token-dir ~/.bigsdb_tokens -o pubmlst_metadata.tsv
```

Download assemblies and isolate metadata (ST, clonal complex, provenance) by isolate id. Both are
resumable and skip what is already present. `--token-dir` is only needed for records that are not
publicly released.

---

## Output

`campypump run` writes one row per genome, 52 columns. The ones that matter first:

| Column | Example | Meaning |
|---|---|---|
| `sample` | `id-152096` | genome file name without extension |
| `cmeB_type_final` | `RE-cmeB` | the call: RE-cmeB, wildtype-cmeB, mosaic, ambiguous, unresolved, not_found |
| `re_cmeB` | `1` | 1 / 0 for analysis; empty for mosaic and unresolved |
| `cmeB_identity_RE_PA001_pct` | `97.438` | % identity to the RE reference — the value the call rests on |
| `cmeB_identity_WT_11168_pct` | `79.776` | % identity to the wildtype reference |
| `cmeB_closest_reference` | `CAA03704-009522 (Cc RE)` | most similar allele in the panel |
| `cmeB_partial_gene` | `no` | gene broken across contigs; its call is provisional |
| `cmeB_divergent` | `no` | complete gene < 90 % identical to every reference |
| `cmeB_RE_like_segments_pct` | `94.7` | % of windows clearly closer to the RE reference |
| `cmeB_WT_like_segments_pct` | `5.3` | % clearly closer to a wildtype reference; both ≥ 10 % ⇒ mosaic |
| `warnings` | | plain-language QC notes; empty means no problems |

Also reported: identity to the published DH161 RE allele and whether the genome meets that
definition, identities for cmeR, cmeA and cmeC, the promoter length and CmeR-binding site sequence,
the locus coordinates, and identity to each panel allele.

The CmeR operator is described by its structure rather than by identity to one strain:

| Column | Example | Meaning |
|---|---|---|
| `CmeR_IR_halfsite_bp` | `7` | longest n for which the first and last n bases of the window are reverse complements; 0 when the half-sites are broken |
| `CmeR_IR_spacer_bp` | `2` | bases between the half-sites; blank when the half-sites are broken and the boundary is arbitrary |
| `CmeR_IR_length_bp` | `16` | total length of the window |
| `CmeR_IR_class` | `intact` | intact, half_site_changed, spacer_changed, half_site_and_spacer_changed, unclassified |

NCTC 11168 and *C. jejuni* 81-176 differ by one base here and are both intact — 11168 a 6-bp
repeat with a 4-bp spacer, 81-176 a 7-bp repeat with a 2-bp spacer — so classifying against
either strain alone makes the other look altered. `CmeR_binding_site_as_11168` is still written
for compatibility but is deprecated; use `CmeR_IR_class`.

Example, from the bundled reference genomes:

```
sample                        cmeB_type_final   re_cmeB   %RE      %WT      closest reference
id-105572_PA001               RE-cmeB           1         100.0    79.598   PA001 (Cj RE)
id-150796_EF0003911-007703    wildtype-cmeB     0         79.662   97.983   EF0003911-007703 (Cj WT)
id-152096_CAA03704-009522     RE-cmeB           1         97.438   79.776   CAA03704-009522 (Cc RE)
```

Optional outputs: the extracted alleles (`--alleles`), SVG segment plots (`--segment-plots`), and
the column guide written on every run.

---

## Interpreting the results

**Check the gap before trusting the calls.** Run `campypump figure` and look at panel A: RE and
wildtype should sit in opposite corners with nothing between them. Genomes in the middle mean the
cut-offs need examining for your collection, and `--re-min` / `--wt-max-re` can be moved.

**`ambiguous` means what it says.** The allele is between 85 % and 95 % identical to the RE
reference. Look at its segment plot: a mosaic shows blocks of each type, while a uniformly
intermediate allele is a divergent lineage.

**`incomplete` is about assembly, not biology.** cmeB is broken across contigs. `cmeB_type_final`
still gives a provisional call from the part present, flagged by `cmeB_partial_gene = yes`. A
premature stop codon also appears this way, so a real nonsense mutation looks like a partial gene.

**`divergent` often means the species is not what you think.** A complete gene under 90 % identical
to every reference is a long way outside *C. jejuni* and *C. coli* diversity. In the development
collection the one divergent allele turned out to be a *C. upsaliensis* genome filed under the wrong
species.

**No single cme gene identifies the species.** The operon crosses between *C. coli* and *C. jejuni*,
and *C. coli* DH161 carries a 99.4 % *C. jejuni*-like cmeR. Use whole-genome ANI for species; the
`operon_*` columns describe the operon only.

---

## Optional extras

### Species check by whole-genome ANI

```bash
# nothing to set up: campypump ships the four reference genomes and runs skani itself
campypump run assemblies/ -t 8 -o results.tsv

# skani not installed? install it and the species columns appear:
micromamba install -c conda-forge -c bioconda skani -y

# a big collection, typed in batches: make the table once, reuse it
skani dist --ql queries.txt --rl refs.txt -t 8 -o species_ani.tsv
campypump run assemblies/ --ani species_ani.tsv -o results.tsv
campypump annotate results.tsv --ani species_ani.tsv -o annotated.csv
```

Adds `species_by_ANI`, the ANI value, the margin to the next species, the aligned fraction and a
`CHECK` flag where ANI disagrees with the label. In the development collection this found 13
genomes whose submitted species was wrong, all confirmed by rMLST.

### Allele phylogeny

```bash
micromamba install -c conda-forge -c bioconda mafft iqtree -y
campypump alleles cmeB_alleles.fasta --results annotated.csv -o tree/cmeB_alleles
mafft --auto --thread 8 tree/cmeB_alleles.fasta > tree/cmeB_alleles.aln.fasta
iqtree -s tree/cmeB_alleles.aln.fasta -m MFP -B 1000 -alrt 1000 -T 8 --prefix tree/cmeB_tree
campypump treeplot tree/cmeB_tree.treefile --annotation tree/cmeB_alleles.tsv -o figures/tree.svg
```

campypump prepares the input and draws the result; MAFFT and IQ-TREE do the alignment and
inference. Full walkthrough in [`docs/tree_howto.md`](docs/tree_howto.md).

---

## Reference panel

`src/campypump/data/references.fasta` holds cmeR/A/B/C for *C. jejuni* NCTC 11168 (wildtype),
*C. jejuni* PA001 (RE) and *C. coli* DH161 (the published RE-cmeABC, GenBank KT778508.1), plus six
further cmeB alleles spanning both species and both types, and the 11168 promoter. Use your own with
`--refs`.

Both definitions of RE are reported: campypump's (≥ 95 % identity to PA001) and the published one
(> 90 % identity to DH161). They agreed on every genome in the validation set.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `blastn not found on PATH` | the environment is not active: `micromamba activate campypump` |
| `campypump: command not found` | same, or `pip install .` ran in a different environment |
| Every genome comes back `not_found` | the folder holds reads or annotations rather than assembled contigs |
| Many genomes flagged `incomplete` | fragmented assemblies; cmeB is broken across contigs and its call is provisional |
| Many flagged `divergent` | the genomes may not be *C. jejuni* or *C. coli*; check with ANI or rMLST |
| `verify_install.sh` reports MISMATCH | note the ids and values, and your Biopython and BLAST versions, and open an issue |

---

## Documentation

| File | Content |
|---|---|
| [`docs/hpc_quickstart.md`](docs/hpc_quickstart.md) | start here for a cluster: install, run, read the results, optional extras |
| [`docs/technical_description.md`](docs/technical_description.md) | the method in full, with code, parameters, validation and known limitations |
| [`docs/codebook_annotated_csv.md`](docs/codebook_annotated_csv.md) | every column of the analysis CSV |
| [`docs/tree_howto.md`](docs/tree_howto.md) | building and drawing a cmeB allele tree |
| `hpc/*.slurm` | Slurm job templates: typing, PubMLST download, metadata |
| `CHANGELOG.md` | what changed in each version |

Figure colours come from `src/campypump/palettes.py`; edit that one file to recolour every figure.

---

## Citing

See `CITATION.cff`. Please also cite BLAST+ and Biopython, and MAFFT, IQ-TREE or skani where used.

## Licence

MIT (`LICENSE`).
