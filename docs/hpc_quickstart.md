# Running campypump on an HPC cluster

From nothing to results in about 25 minutes of your own time, most of it installation. Written for a
Slurm cluster; adapt the job scripts if yours uses something else.

You need: an account on a cluster with conda or micromamba available, and genome assemblies in FASTA
format, one file per isolate.

---

## 1. Install (once, about 15 minutes)

Work on a compute node, not a login node. On a Slurm cluster that usually means something like
`srun --pty -A YOUR_ACCOUNT -n 4 -t 02:00:00 bash` or a site-specific wrapper such as `interactive`.

```bash
module load micromamba          # or: module load anaconda / mamba, whatever your site provides
micromamba create -n campypump -c conda-forge -c bioconda \
    python=3.11 biopython blast pip pytest -y
micromamba activate campypump

git clone https://github.com/fransch/campypump.git
cd campypump
pip install .
```

Only Python, Biopython and BLAST+ are required. `pandas` and `openpyxl` are optional, for reading
Excel id lists.

**Check it works** before running anything real:

```bash
bash verify_install.sh
```

This runs the test suite, types the eight bundled reference genomes, and compares the output with
expected values, three of which are identity values from the manual Geneious workflow campypump
replaces. It should end with `ALL CHECKS PASSED`.

## 2. Put your assemblies somewhere with space

```bash
mkdir -p /path/to/your/project/assemblies
# copy or link your FASTA files there: one per isolate
ls /path/to/your/project/assemblies | head
```

Accepted extensions are `.fasta`, `.fas`, `.fa` and `.fna`. If files are named `id-<number>.fasta`
(as `campypump fetch` writes them) the number is carried through to a `pubmlst_id` column.

Use scratch or project storage rather than your home directory; genome collections run to several
GB. Remember that scratch is usually not backed up, so copy results elsewhere.

## 3. Type them

For a few hundred genomes, an interactive session is enough:

```bash
cd /path/to/your/project
campypump run assemblies/ -t 8 -o results/campypump_results.tsv \
    --alleles results/cmeB_alleles.fasta --segment-plots results/segment_plots
```

Expect roughly 1.7 CPU-seconds per genome: about 15 minutes for 3,700 genomes on 8 CPUs. Add
`--species "C. jejuni"` to label the output when a run covers one species; run each species
separately and `campypump combine` the tables if you have both.

For larger sets, use the batch template:

```bash
cp campypump/hpc/campypump_job.slurm .
# edit the lines marked <-- EDIT: account, partition, project path, assemblies folder
sbatch campypump_job.slurm
```

## 4. Read the results

`results/campypump_results.tsv` has one row per genome. The columns that matter first:

| Column | Meaning |
|---|---|
| `cmeB_type_final` | RE-cmeB, wildtype-cmeB, mosaic, ambiguous, incomplete or not_found |
| `re_cmeB` | 1 / 0, for analysis |
| `cmeB_identity_RE_PA001_pct` | % identity to the RE reference — the value the call rests on |
| `cmeB_divergent`, `cmeB_partial_gene` | flags for genes that need a look |
| `warnings` | plain-language QC notes; empty means no problems |

Every run also writes `<output>_column_guide.txt` describing each column, and
`docs/codebook_annotated_csv.md` documents the fuller annotated file.

A quick summary:

```bash
cut -f4 results/campypump_results.tsv | sort | uniq -c
```

## 5. Figures (seconds)

```bash
campypump figure results/campypump_results.tsv -o results/figures
campypump diversity results/cmeB_alleles.fasta \
    -o results/figures/cmeB_allele_diversity.svg --table results/figures/alleles.tsv
```

The first plots each genome's identity to the RE and wildtype references, which is how you check the
classes separate in *your* collection rather than trusting a threshold. The second counts distinct
alleles: rank abundance, a rarefaction curve, and alleles per cmeB type.

SVG files open in any browser. Colours come from `src/campypump/palettes.py`; edit that one file to
recolour every figure.

---

## Optional extras

### Genomes from PubMLST

If your isolates are in PubMLST, campypump can fetch them by id rather than you downloading by hand.
Put one isolate id per line in `ids.txt`, then:

```bash
campypump fetch --ids ids.txt --outdir assemblies/
```

Records that are not publicly released return HTTP 404 and are listed in
`assemblies/failed_ids.txt`. To download those you need an account with access to them, plus
BIGSdb_downloader and a personal API key from your PubMLST profile:

```bash
pip install bigsdb-downloader
mkdir -p ~/.bigsdb_tokens && chmod 700 ~/.bigsdb_tokens
bigsdb-downloader --key_name PubMLST --site PubMLST \
    --db pubmlst_campylobacter_isolates --setup --token_dir ~/.bigsdb_tokens
```

`campypump metadata --ids ids.txt --token-dir ~/.bigsdb_tokens -o results/pubmlst_metadata.tsv`
then retrieves ST, clonal complex and provenance for each isolate, and
`campypump annotate` merges that into the results.

### Species check by whole-genome ANI

No single cme gene identifies the species: the operon crosses between *C. coli* and *C. jejuni*. If
species assignment matters for your analysis, compute whole-genome ANI separately and merge it in:

```bash
micromamba install -c conda-forge -c bioconda skani -y
mkdir species_refs      # one reference genome per species, named for the organism
skani dist -q assemblies/*.fasta -r species_refs/*.fasta -t 8 -o species_ani.tsv
campypump annotate results/campypump_results.tsv --ani species_ani.tsv \
    -o results/campypump_annotated.csv
```

In the collection campypump was developed on, this found thirteen genomes whose submitted species
was wrong, every one confirmed by rMLST.

### Allele phylogeny

```bash
micromamba install -c conda-forge -c bioconda mafft iqtree -y
campypump alleles results/cmeB_alleles.fasta --results results/campypump_annotated.csv \
    -o tree/cmeB_alleles
mafft --auto --thread 8 tree/cmeB_alleles.fasta > tree/cmeB_alleles.aln.fasta
iqtree -s tree/cmeB_alleles.aln.fasta -m MFP -B 1000 -alrt 1000 -T 8 --prefix tree/cmeB_tree
campypump treeplot tree/cmeB_tree.treefile --annotation tree/cmeB_alleles.tsv \
    -o results/figures/cmeB_tree.svg
```

`campypump alleles` collapses genomes to distinct alleles and writes tree-safe names, an annotation
table, iTOL colour strips and a Microreact CSV. Full detail in `docs/tree_howto.md`.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `blastn not found on PATH` | the environment is not active: `micromamba activate campypump` |
| `campypump: command not found` | same, or `pip install .` was run in a different environment |
| `verify_install.sh` reports MISMATCH | note the ids and values, plus your Biopython and BLAST versions, and open an issue |
| Every genome comes back `not_found` | the folder holds reads or annotations rather than assembled contigs |
| Many genomes flagged `incomplete` | the assemblies are fragmented; cmeB is broken across contigs and its call is provisional |
| A job stays pending | normal on a busy cluster; it starts on its own |

For what the method does and why, see `docs/technical_description.md`; it also lists the known
limitations.
