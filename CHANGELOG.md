# Changelog

## 0.2.1

- **`docs/hpc_quickstart.md`:** install, run, read the results and the optional extras (PubMLST
  download, ANI species check, allele phylogeny) for someone bringing their own genomes to a
  cluster.
- The Slurm templates in `hpc/` are now genuinely generic: they type a folder of assemblies rather
  than referring to one dataset's isolate id lists, and no longer assume a particular site.

## 0.2.0

- **Repository layout for public release:** the tool is separated from the dataset it was developed
  on. `hpc/` holds generic cluster templates; dataset-specific results, isolate id lists and as-run
  scripts are kept outside the repository until the corresponding manuscript is published.
- Added `CITATION.cff`, `CHANGELOG.md`, `.gitignore`, and a GitHub Actions workflow running the
  tests on Python 3.9, 3.10, 3.11 and 3.12.
- `README.md` rewritten for someone bringing their own genomes.
- No change to any analysis logic: results from 0.1.19 are unaffected.

## 0.1.19

- **Fix:** 0.1.18 failed to import on Python 3.11 — a nested same-quote f-string, which is only
  legal from 3.12. A test now parses every module under the 3.9 and 3.11 grammars, so this cannot
  recur regardless of the interpreter used to develop.

## 0.1.18

- **One palette for every figure**, in `palettes.py`: crimson / thistle / royal blue / green / gold,
  with two darkened or lightened variants for line work and paired categories. Editing that file
  recolours the classification, diversity, segment, tree and iTOL/Microreact outputs together.
- **Fix:** `campypump figure` fell back to no species split when reading an annotated CSV, which
  carries `species_by_ANI` rather than `species`; it now reads either.

## 0.1.17

- **`treeplot` shows more per allele:** tip size scales with the number of genomes carrying the
  allele, an orange ring marks alleles carrying the wildtype-like block, and a column of stacked
  bars gives each allele's source composition. `--no-size-by-genomes` and `--no-source-bars` turn
  them off.

## 0.1.16

- **`treeplot --annotation`:** appends the dominant ST (and how many of the allele's genomes carry
  it) to each labelled tip, read from the `campypump alleles` table. `--annotate-cols` picks other
  columns, e.g. `top_sources`.

## 0.1.15

- **`campypump treeplot`:** annotated SVG tree from a Newick file — branches coloured by cmeB type
  (or by species inside the RE clade), species shading, support labels above a threshold, scale bar.
  `--subtree RE`, `--drop-outgroup` and `--cladogram` handle the very long RE/wildtype branch.

## 0.1.14

- **`campypump alleles`** also writes `<prefix>_microreact.csv`: an ID column matching the tree tips
  plus `__colour` columns, for https://microreact.org.

## 0.1.13

- The core annotation now also carries the cmeR, cmeA and cmeC identities (to both NCTC 11168 and
  PA001) and the promoter length, alongside the CmeR-binding-site sequence it already had — 39
  columns.

## 0.1.12

- **`annotate --columns`:** the annotated CSV now defaults to the ~39 columns analyses actually use
  (`core`); `--columns full` writes all 82. Everything dropped is still in the results table the
  annotation was built from.

## 0.1.11

- **Species call from whole-genome ANI:** `annotate --ani species_ani.tsv` merges a skani table and
  adds `species_by_ANI`, `ANI_pct`, `ANI_margin`, `ANI_align_fraction_pct`, `ANI_call_quality` and
  `species_check_ANI`. ANI is unaffected by a single transferred locus, so it is the species call.
- **The operon flag is demoted to what it measures:** `species_check_operon` becomes
  `operon_atypical_for_label` — a database error *or* complete introgression of the operon, not a
  species call.
- **`annotate --exclude`:** mark genomes as not-for-analysis (`use_for_prevalence = no` plus an
  `exclude_reason`) without deleting their results.

## 0.1.10

- **Species screen:** new `operon_species_signal`, `operon_Cj_like_genes` and `species_check_operon`
  columns. cmeR, cmeA and cmeC are each scored against NCTC 11168; the operon is called for a
  species only when all three agree, and `CHECK` marks a disagreement with the `--species` label.
  No single cme gene can do this — the operon crosses species, and *C. coli* DH161 carries a
  99.4 % *C. jejuni*-like cmeR. Confirm any CHECK with rMLST.

## 0.1.9

- **Fix:** `alleles --outgroup`/`--references` no longer doubles the tag when a record is already
  named `OUT_...` or `REF_...` (which produced `OUT_OUT_...` and broke IQ-TREE's `-o`).

## 0.1.8

- **`campypump alleles`:** deduplicates cmeB alleles into tree input — a FASTA with Newick-safe ids
  (`Cc_a003_RE_n63`), an annotation table (species, type, genomes, top STs, sources, block status),
  and iTOL colour strips. `--references` and `--outgroup` add published alleles and a root.
- **`docs/tree_howto.md`:** MAFFT + IQ-TREE walkthrough, including the recombination split.

## 0.1.7

- **Fix:** in the classification figure, a bin holding a single genome was drawn with zero height and
  vanished — this hid the mosaic isolate. Bars now have a minimum height.
- The figure subtitle states that partial genes are excluded, so the n is unambiguous.

## 0.1.6

- **`campypump figure`:** classification figure (panels A–B) from a results table — identity to the
  RE and wildtype references per genome, and the distribution of identity to the RE reference.
- **`campypump diversity`:** allele diversity figure (panels A–C) from the allele FASTAs — rank
  abundance, rarefaction, and distinct alleles per cmeB type; `--table` writes the allele table.

## 0.1.5

- **Fix:** `RE_by_DH161_definition` now uses identity over the part of the gene present for partial
  genes. Previously a fragment's deflated whole-gene identity could read as "not RE" even when the
  part present was 100 % RE.

## 0.1.4

- **Published RE reference added:** the *C. coli* DH161 operon (GenBank KT778508.1, Yao et al. 2016)
  is now in the panel, with `cmeB_identity_DH161_pct` and `RE_by_DH161_definition` (> 90 % identity,
  the definition used by Cooper et al. 2025) reported for every genome.

## 0.1.3

- **Partial genes:** provisional RE or wildtype call (`cmeB_type_final`), plus identity over the part
  present (`cmeB_aligned_identity_*`).
- **Divergent alleles:** complete genes below 90 % identity to every reference are flagged.
- **`--segment-plots`:** SVG and TSV plots of identity along cmeB for mosaic, divergent and mixed alleles.
- **`campypump metadata`:** MLST and provenance from PubMLST, including private records via
  BIGSdb_downloader.
- **`campypump annotate`:** merged analysis CSV with duplicate-record and species-mismatch checks.

## Validation status

- **All 1,513 Peru *C. coli*:** calls agree 100 % with Geneious; 1,497 of 1,508 complete genes are within
  0.01 %. One Geneious transcription error was identified (CAA04817-012617).
- **8 reference genomes:** all called correctly. The three *C. coli* genomes reproduce the Geneious
  values to within 0.005 %.
- **Synthetic tests:** a RE/WT chimera is called `mosaic`, a gene cut at a contig edge is called
  `incomplete`, and a genome without the gene is called `not_found`.
- **Peru *C. jejuni* (2,223):** classes are fully separated (RE >= 97.6 %, wildtype <= 79.7 % to PA001).
