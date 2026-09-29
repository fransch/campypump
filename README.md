<img src="assets/campypump_logo.png" alt="campypump" width="150" align="right"/>

# campypump

Classify the cmeB gene of the *Campylobacter* **cmeABC** efflux operon as **RE-cmeB** — the
resistance-enhancing variant that raises fluoroquinolone and macrolide MICs — or as **wildtype**,
directly from genome assemblies.

It replaces a manual BLAST-and-align workflow and reproduces its identity values exactly: on 1,513
*C. coli* genomes, calls agreed with the manual method for every genome, and 1,497 of 1,508 complete
genes matched to within 0.01 percentage points.

```bash
campypump run genomes/ -t 8 -o results.tsv
```

| | |
|---|---|
| Input | genome assemblies (FASTA), one file per isolate |
| Output | one row per genome: cmeB type, identity to reference alleles, operon context, promoter, QC flags |
| Dependencies | Python ≥ 3.9, Biopython, BLAST+ |
| Speed | ~1.7 CPU-seconds per genome (3,736 genomes in 13 minutes on 8 CPUs) |

## Install

```bash
conda install -c conda-forge -c bioconda campypump      # once the bioconda recipe is accepted
```

or from source:

```bash
git clone https://github.com/fransch/campypump && cd campypump
conda env create -f environment.yml && conda activate campypump
pip install .
bash verify_install.sh        # tests + the 8 bundled reference genomes, with expected values
```

On a cluster, follow `docs/hpc_quickstart.md`.

## What it does

1. **Finds** cmeR, cmeA, cmeB and cmeC in each assembly with sensitive BLAST settings, and extracts
   each as a complete ORF.
2. **Measures** identity of cmeB against a bundled panel of reference alleles, using a global
   alignment scored like Geneious Prime's default, so values are comparable with manual work.
3. **Classifies** each gene as RE-cmeB (≥ 95 % identity to the RE reference), wildtype (≤ 85 %),
   mosaic (segments of both), or flags it as partial, divergent or ambiguous.
4. **Reports** the rest of the operon, the CmeR-binding inverted repeat in the promoter, and
   plain-language QC warnings.

In a collection of 3,736 Peru genomes the two classes were separated by a wide empty gap: RE alleles
were 97.4–100 % identical to the RE reference and wildtype alleles ≤ 80.2 %, with nothing in
between. A maximum-likelihood tree of the alleles recovered the RE set as a single clade with full
support, so the cut-off corresponds to a real evolutionary group.

## Commands

| Command | Purpose |
|---|---|
| `campypump run` | type assemblies; writes the results table, extracted alleles, segment plots |
| `campypump combine` | stack result tables (e.g. one per species) |
| `campypump annotate` | merge metadata, a whole-genome ANI species call, and exclusions into one analysis CSV |
| `campypump figure` | classification figure (identity scatter + distribution), per species |
| `campypump diversity` | allele diversity figure: rank abundance, rarefaction, alleles per type |
| `campypump alleles` | deduplicate alleles into tree input, with annotation, iTOL and Microreact files |
| `campypump treeplot` | annotated tree figure from a Newick file |
| `campypump fetch`, `campypump metadata` | download assemblies and isolate metadata from PubMLST |

Run `campypump <command> --help` for options. Every run also writes a column guide explaining each
output column.

## Typical workflow with your own genomes

```bash
# 1. type them (add --species to label the output)
campypump run assemblies/ -t 8 --species "C. jejuni" -o results.tsv \
    --alleles cmeB_alleles.fasta --segment-plots plots/

# 2. optional: species check by whole-genome ANI (skani, run separately)
skani dist -q assemblies/*.fasta -r species_refs/*.fasta -t 8 -o species_ani.tsv
campypump annotate results.tsv --ani species_ani.tsv -o annotated.csv

# 3. figures
campypump figure annotated.csv -o figures/
campypump diversity cmeB_alleles.fasta -o figures/diversity.svg --table figures/alleles.tsv

# 4. optional: allele phylogeny (MAFFT + IQ-TREE, run separately)
campypump alleles cmeB_alleles.fasta --results annotated.csv -o tree/cmeB_alleles
mafft --auto tree/cmeB_alleles.fasta > tree/cmeB_alleles.aln.fasta
iqtree -s tree/cmeB_alleles.aln.fasta -m MFP -B 1000 -T 8 --prefix tree/cmeB_tree
campypump treeplot tree/cmeB_tree.treefile --annotation tree/cmeB_alleles.tsv -o figures/tree.svg
```

campypump does not infer phylogenies or compute ANI: those stay with MAFFT, IQ-TREE and skani, and
campypump prepares their input and annotates their output.

## Reference panel

`src/campypump/data/references.fasta` holds cmeR/A/B/C for *C. jejuni* NCTC 11168 (wildtype),
*C. jejuni* PA001 (RE) and *C. coli* DH161 (the published RE-cmeABC, GenBank KT778508.1), plus six
further cmeB alleles spanning *C. coli* and *C. jejuni*, wildtype and RE, and the 11168 promoter.
Supply your own with `--refs`.

Both the campypump definition of RE (≥ 95 % to PA001) and the published one (> 90 % to DH161) are
reported, and they agreed on every genome in the validation set.

## Documentation

| File | Content |
|---|---|
| `docs/technical_description.md` | the method in full, with code, parameters, validation and known limitations |
| `docs/codebook_annotated_csv.md` | every column of the analysis CSV |
| `docs/tree_howto.md` | building and drawing a cmeB allele tree |
| `docs/hpc_quickstart.md` | **start here for a cluster:** install, run, read the results, optional extras |
| `hpc/*.slurm` | Slurm job templates (typing, PubMLST download, metadata) |

## Colours

Every figure takes its colours from `src/campypump/palettes.py`. Edit that one file to recolour the
classification, diversity, segment and tree figures together.

## Citing

See `CITATION.cff`. Please also cite BLAST+, Biopython, and MAFFT/IQ-TREE/skani where used.

## Licence

MIT (`LICENSE`).
