## 0.5.0

The CmeR operator is now described by its structure instead of by identity to one
reference strain.

CmeR binds as a dimer, one subunit to each half-site of an inverted repeat, either side
of a spacer that must hold the two half-sites the right distance apart. NCTC 11168 and
C. jejuni 81-176 differ by a single base in this window: 11168 forms a perfect 6-bp
inverted repeat with a 4-bp spacer, 81-176 a perfect 7-bp repeat with a 2-bp spacer.
Both are intact operators, but classifying against either one alone makes the other look
altered, and 11168 carries the rarer of the two forms.

Added, all computed from the window campypump already reported, so no genome needs
retyping:

* `CmeR_IR_halfsite_bp` - longest n for which the first and last n bases of the window
  are reverse complements
* `CmeR_IR_spacer_bp` - bases between the half-sites; blank when the half-sites are
  broken and the boundary is arbitrary
* `CmeR_IR_length_bp` - total length of the window
* `CmeR_IR_class` - intact / half_site_changed / spacer_changed /
  half_site_and_spacer_changed / unclassified

`CmeR_binding_site_as_11168` is kept so nothing downstream breaks, but is deprecated in
the codebook. Thirteen regression tests cover the new classifier.

# Changelog

## 0.4.1

- **A directory argument to `run` is now searched at every depth.** Collections are
  routinely laid out one directory per species or per isolate
  (`genomes/ccoli/id-48.fasta`), and pointing `run` at the collection root found nothing
  and said only "No genome FASTA files given". It now walks the tree, reports how many
  files it found under each directory given, and names the directory and the extensions
  it searched when it finds none.
- A path that does not exist is reported as such instead of being passed to BLAST.

## 0.4.0

- **The four species reference genomes now ship with campypump**, gzipped in
  `data/species/` (1.9 MB): *C. jejuni* NCTC 11168 (AL111168.1), *C. coli* RM4661
  (CP007181.1), *C. lari* RM2100 (CP000932.1) and *C. upsaliensis* RM3940 (CP053849.1).
  They are unpacked once into a cache directory, so skani is always handed plain FASTA.
- **`campypump run` makes the species call itself.** It runs `skani dist` against those
  genomes and writes `species_by_ANI`, `ANI_pct`, `ANI_margin`,
  `ANI_align_fraction_pct`, `ANI_call_quality` and `species_check_ANI` straight into the
  results table. Until now the user had to find four genomes and run a second tool, which
  is exactly the user most likely to be working from a mislabelled record.
- Query and reference paths are passed to skani in list files (`--ql`, `--rl`), not on the
  command line: a few thousand genomes would otherwise overflow the shell's argument limit.
- `--ani TSV` on `run` reuses a table you made yourself; `--no-ani` skips the call. The
  existing `annotate --ani` is unchanged, so the HPC workflow still works as before.
- **Without skani the run still types everything.** It prints the install command, leaves
  the species columns empty and carries on. skani is a run dependency of the conda package,
  so this fallback is for pip installs.
- `CAMPYPUMP_SPECIES_REFS` points the call at your own directory of reference genomes; file
  names supply the species, so keep the naming convention.
- No change to the reference panel, the thresholds, the mosaic scan or any cmeB, cmeR, cmeA
  or cmeC value. 0.3.2 results remain reproducible.

## 0.3.2

- **DH161's cmeR, cmeA and cmeC are removed from the panel.** DH161 is a *C. coli* strain whose
  cmeABC operon is partly *C. jejuni*-derived. Its cmeR and cmeA are *C. jejuni*-like (75.2 % and
  89.2 % to the native *C. coli* reference, against 99.4 % and 98.6 % to NCTC 11168) while
  carrying the class `Cc_RE`. Those two records were inert before 0.3.0, but 0.3.0 made every
  reference eligible for `<gene>_closest_reference`, so a genome with a *C. jejuni*-type flanking
  gene could be reported as matching a *C. coli* reference. In the Peru collection that affected
  93 *C. coli* and 138 *C. jejuni* at cmeR, and 120 and 195 at cmeA.
- DH161's cmeC is **not** mislabelled: at 97.8 % to RM4661 against 85.1 % to NCTC 11168 it is
  genuinely native *C. coli*, which is consistent with the introgressed block in that strain
  spanning cmeR to cmeB and stopping before cmeC. It is removed only for consistency and because
  RM4661 and CAA03009-008347 already cover native *C. coli* cmeC.
- A test now fails if any reference whose class says *C. coli* is closer to the *C. jejuni*
  reference than to the native *C. coli* one, for any of the three flanking genes.
- DH161's **cmeB is unchanged** and still provides the published > 90 % RE definition.
- cmeR, cmeA and cmeC now have four references each: NCTC 11168, PA001, RM4661 and
  CAA03009-008347.

## 0.3.1

- **Fixes a crash in 0.3.0.** The new per-gene loop for cmeR, cmeA and cmeC bound a local
  variable `ids`, which shadowed the cmeB identities computed earlier in the same function. The
  `panel_<class>_<strain>_pct` columns are written from those cmeB identities, so every genome
  raised `KeyError` and 0.3.0 could not type anything. The flanking-gene loop now uses `gids`,
  the cmeB panel columns are written while the cmeB identities are still in scope rather than at
  the end of the function, and a test asserts `ids` is bound exactly once in `classify_genome`.
  0.3.0 should not be used.

## 0.3.0

- **cmeR, cmeA and cmeC are now scored against a native *C. coli* reference.** Until now these
  three genes were compared only with NCTC 11168 and PA001, both *C. jejuni*, so a *C. coli* gene
  could only be reported as its distance from the *C. jejuni* version and its absolute value was
  not interpretable. The panel gains *C. coli* RM4661 (CP007181.1) as the native wildtype
  yardstick, giving a new column per gene: `cmeR_identity_CcWT_RM4661_pct` and the same for cmeA
  and cmeC. RM4661 was chosen because it is a published complete *C. coli* genome and is already
  the *C. coli* reference used for the ANI species call.
- **Each of the three genes now reports its closest reference**, `<gene>_closest_reference` and
  `<gene>_closest_identity_pct`, mirroring what cmeB has always done. This answers "which allele
  is this gene actually like" rather than only "how far is it from a *C. jejuni* gene".
- **Fixed: a third primary reference silently overwrote the NCTC 11168 column.** The per-gene
  scoring selected its output column with "PA001, or else 11168", so any primary reference that
  was not PA001 was written into the 11168 column. Columns are now selected by reference name.
  This had no effect before 0.3.0, when only two primary references existed per gene.
- The flanking genes of *C. coli* CAA03009-008347 are included as `panel` references. They feed
  the closest-match lookup but produce no column, so there is never a second *C. coli* identity
  to choose between.
- **cmeB is unchanged**: same nine references, same thresholds, same RE/wildtype logic. Every
  cmeB value produced by earlier versions is reproduced exactly. RM4661's cmeB was evaluated and
  deliberately left out: it is the closest wildtype reference for 39 of 215 *C. coli* wildtype
  alleles and improves the genome-weighted best match by 0.4 points, but it would also shift the
  wildtype baseline of the mosaic window scan, which the ST-12298 block finding rests on. Adding
  it belongs in a change where that scan is re-run and the block re-checked.
- DH161's cmeR, cmeA and cmeC, carried unused since 0.1.4, are now scored. Note that DH161 is a
  *C. coli* strain whose operon is *C. jejuni*-derived, so its flanking genes are not a native
  *C. coli* comparator.

## 0.2.3

- **Species names are read from any reference file named for its organism**, not only the handful
  that were hardcoded: `Clanienae_NCTC13004.fasta` now gives "C. lanienae" rather than the file
  stem. Add a reference genome for any *Campylobacter* species and ANI will name it correctly.

## 0.2.2

- **README expanded** into a full front page: installation with optional dependencies, a quick
  start, how the method works, every command with its options, the output columns with a worked
  example, how to read the flags, the optional ANI and phylogeny workflows, and troubleshooting.

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
