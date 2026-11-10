# Building a cmeB allele tree

campypump prepares the input and annotates the output; alignment and tree inference use MAFFT and
IQ-TREE. This keeps campypump's dependencies to Biopython and BLAST, and means the phylogeny is
produced by tools a reviewer already trusts.

Total time: about 20 minutes of setup, then 10–30 minutes of compute for ~380 alleles.

---

## 1. Install the tree software (once)

```bash
module load micromamba
micromamba activate campypump
micromamba install -c conda-forge -c bioconda mafft iqtree -y

mafft --version
ls $CONDA_PREFIX/bin | grep -i iqtree     # the binary may be iqtree, iqtree2 or iqtree3
iqtree --version
```

Both install in a couple of minutes. **The command name varies by version** — bioconda has shipped
`iqtree` (v1 and v3), `iqtree2` (v2) and `iqtree3`. Substitute whichever you have in the commands
below. Versions 2 and 3 take the same options. Version 1.6 differs: use `-bb` for `-B`, `-nt` for
`-T`, and `-pre` for `--prefix`.

Confirm the options exist in your build before the long run:

```bash
iqtree -h 2>&1 | grep -E "^\s+(-B|-alrt|--prefix|-o|-T)\b"
```

## 2. Download an outgroup

*C. lari* RM2100 (GenBank CP000932.1) is a closed, complete genome and the standard *C. lari*
reference. Download the genome and let campypump pull cmeB out of it:

```bash
cd "$PROJECT"        # the folder holding genomes/ and results/
mkdir -p tree
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000932.1&rettype=fasta&retmode=text" \
    -o tree/C_lari_RM2100_CP000932.fasta
head -c 80 tree/C_lari_RM2100_CP000932.fasta; echo

campypump run tree/C_lari_RM2100_CP000932.fasta -o tree/lari_check.tsv --alleles tree/lari_cmeB.fasta
cut -f1,4,5,6,9,10 tree/lari_check.tsv
```

**Check the identity before using it.** A usable outgroup should sit somewhere around 65–80 %
identity to the *C. jejuni* and *C. coli* alleles: close enough to align confidently, far enough to
fall outside the ingroup. If it comes back much below ~60 %, the alignment will be unreliable and a
closer outgroup (*C. upsaliensis*, *C. hepaticus*) is the better choice. If it comes back inside the
ingroup range (> 85 %), it is not an outgroup at all.

Rename the record so the tree tip is readable:

```bash
sed -i '1s/.*/>OUT_C_lari_RM2100/' tree/lari_cmeB.fasta
```

## 3. Prepare the alleles

```bash
campypump alleles results/cmeB_alleles_ccoli.fasta results/cmeB_alleles_cjejuni.fasta \
    --results results/campypump_all_peru_annotated.csv \
    --references campypump/src/campypump/data/references.fasta \
    --outgroup tree/lari_cmeB.fasta \
    -o tree/cmeB_alleles
```

This writes:

| File | Content |
|---|---|
| `cmeB_alleles.fasta` | one record per distinct allele, ids like `Cc_a003_RE_n63` (species, allele rank, cmeB type, genomes carrying it) |
| `cmeB_alleles.tsv` | one row per allele: species, type, genomes, top STs, top sources, mean wildtype-like segments, block status, example genome |
| `cmeB_alleles_itol_type.txt` | iTOL colour strip by cmeB type |
| `cmeB_alleles_itol_species.txt` | iTOL colour strip by species |
| `cmeB_alleles_itol_genomes.txt` | iTOL bar chart of genomes per allele (log scale) |

`--references` adds the bundled panel (NCTC 11168, PA001, DH161, CAA03704 and the rest) tagged
`REF_`, so published alleles appear in the tree. Ids are restricted to letters, digits and
underscores, which Newick requires.

## 4. Align

```bash
cd tree
mafft --auto --thread 8 cmeB_alleles.fasta > cmeB_alleles.aln.fasta
grep -c ">" cmeB_alleles.aln.fasta
python3 -c "
from Bio import AlignIO
a = AlignIO.read('cmeB_alleles.aln.fasta', 'fasta')
print(len(a), 'sequences,', a.get_alignment_length(), 'columns')"
```

For ~380 sequences of ~3.1 kb this takes a minute or two. The alignment length should be close to
3,120–3,150; much longer means the outgroup or a fragment is forcing large gaps, so check before
continuing.

## 5. Infer the tree

```bash
iqtree -s cmeB_alleles.aln.fasta -m MFP -B 1000 -alrt 1000 -T 8 --prefix cmeB_tree \
    -o OUT_C_lari_RM2100
```

- `-m MFP` chooses the substitution model automatically (ModelFinder).
- `-B 1000` ultrafast bootstrap, `-alrt 1000` SH-aLRT; a branch is well supported when both are
  high (roughly ≥ 95 and ≥ 80).
- `-o` sets the outgroup for rooting; drop it to leave the tree unrooted.

The tree is `cmeB_tree.treefile` (Newick), with the full log in `cmeB_tree.iqtree`, which records
the model chosen and is what you cite in the methods.

As a batch job for larger runs:

```bash
sbatch --account=YOUR_SLURM_ACCOUNT --partition=standard --nodes=1 --ntasks=8 --time=02:00:00 \
    --wrap "module load micromamba && micromamba activate campypump && cd $PWD && \
    iqtree -s cmeB_alleles.aln.fasta -m MFP -B 1000 -alrt 1000 -T 8 --prefix cmeB_tree -o OUT_C_lari_RM2100"
```

## 6. Recombination caveat

Some *C. coli* RE alleles carry a ~300 bp wildtype-like block, so different parts of cmeB have
different histories and a single tree averages over them. Two ways to handle it:

**Split the gene and compare trees.** Build one tree from the block and one from the flanks; if the
block-carrying alleles move between them, that is direct evidence of recombination.

```bash
python3 - <<'PY'
from Bio import AlignIO
a = AlignIO.read('cmeB_alleles.aln.fasta', 'fasta')
# adjust to the block coordinates seen in the segment plots (1-based, inclusive)
start, end = 1100, 1400
AlignIO.write(a[:, start-1:end], 'cmeB_block.aln.fasta', 'fasta')
flanks = a[:, :start-1] + a[:, end:]
AlignIO.write(flanks, 'cmeB_flanks.aln.fasta', 'fasta')
PY

for part in block flanks; do
  iqtree -s cmeB_$part.aln.fasta -m MFP -B 1000 -T 8 --prefix cmeB_$part -o OUT_C_lari_RM2100
done
```

**Or test formally** with RDP, GARD (HyPhy) or 3SEQ on the alignment. The campypump window scan is a
screen, not a test.

## 7. View and annotate

**campypump treeplot** (no upload, publication-ready SVG):

```bash
# the RE clade, coloured by species — the figure that shows the interspersion
campypump treeplot cmeB_tree.treefile --outgroup OUT_C_lari_RM2100 --subtree RE \
    --annotation cmeB_alleles.tsv \
    -o figures/cmeB_tree_RE.svg --title "cmeB: the RE-cmeB clade"

# whole tree; RE and wildtype are ~0.4 substitutions/site apart, so use a cladogram
campypump treeplot cmeB_tree.treefile --outgroup OUT_C_lari_RM2100 \
    --drop-outgroup --cladogram --tip-height 5 --min-support 95 \
    -o figures/cmeB_tree_all.svg
```

With `--annotation` the figure also scales each tip by how many genomes carry the allele, rings the
alleles carrying the wildtype-like block, and draws a source-composition bar per allele
(`--no-size-by-genomes`, `--no-source-bars` to turn those off). `--annotation` adds the dominant ST
to each labelled tip, as `ST 12298 (59/63)` — the ST and how
many of that allele's genomes carry it. `--annotate-cols top_sources` shows sources instead, or
both. `--colour-by type|species`, `--labels notable|all|none`, `--min-support`, `--width` and
`--tip-height` adjust the rest.

**iTOL** (best fit for an allele tree): upload `cmeB_tree.treefile` at https://itol.embl.de, then
drag the three `*_itol_*.txt` files onto the tree. Show support with Advanced → Bootstrap/metadata
→ Display; the labels are `SH-aLRT/UFboot`, so set a threshold (e.g. 80) to avoid clutter.

To change the colours, edit the generated files: each data line is `<tip id> <tab> <hex> <tab>
<label>`, and the `LEGEND_COLORS` line near the top sets the legend swatches. Change both so the
legend matches the tips. Colours can also be edited in iTOL's dataset panel after upload.

**Microreact** (https://microreact.org): upload `cmeB_tree.treefile` and
`cmeB_alleles_microreact.csv` together. The CSV has an `ID` column matching the tips and
`<column>__colour` columns that fix the colours; edit those hex values to recolour, or change them
in the colour palette box in Microreact. Note that the tips are alleles rather than isolates, so the
map and timeline panels have nothing to show — Microreact is the better choice once you have a
core-genome tree of isolates.

**ggtree** in R, if you want a figure you can version:

```r
library(ggtree); library(treeio); library(ggplot2)
tr   <- read.newick("cmeB_tree.treefile")
meta <- read.delim("cmeB_alleles.tsv")
ggtree(tr) %<+% meta +
  geom_tippoint(aes(colour = cmeB_type, size = genomes)) +
  scale_colour_manual(values = c("RE-cmeB" = "#B2182B", "wildtype-cmeB" = "#2166AC",
                                 "mosaic" = "#E08214", "REF" = "#111111", "OUT" = "#4D9221")) +
  scale_size_continuous(range = c(1, 6), trans = "log10")
```

## 8. What the tree should answer

1. **Do the *C. coli* RE alleles form one clade or two?** Two clades would support two separate
   arrivals of RE-cmeB in Peruvian *C. coli*; one clade with internal structure would favour a single
   arrival followed by diversification.
2. **Where does `Cc_a003_RE_n63` sit?** That is the block-carrying allele, almost entirely ST-12298.
   Whether it groups with the other *C. coli* RE alleles or with wildtype alleles should differ
   between the block tree and the flank tree.
3. **Do *C. coli* RE alleles fall inside the *C. jejuni* RE clade?** That is the direct test of the
   *C. jejuni* origin inferred by Cooper et al. (2025).
4. **Where do DH161, PA001 and CAA03704 fall** relative to the Peru alleles.
