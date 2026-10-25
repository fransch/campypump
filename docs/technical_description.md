# campypump: technical description for review

**Tool:** campypump (*Campylobacter* efflux PUMP typer), version 0.1.19
**Purpose of this document:** to let a genomic scientist or bioinformatician check, step by step,
that the analysis does what we intend, and to flag design choices and limitations for review.
**Companion documents:** `docs/codebook_annotated_csv.md` (every column of the analysis file),
the results summary (not yet public) (findings), `docs/tree_howto.md` (phylogeny), `INSTRUCTIONS.md` (how to run it).
**Project context:** RE-cmeABC typing of all Peru *C. coli* (n = 1,513) and *C. jejuni* (n = 2,223)
genomes in PubMLST, replacing a manual Geneious Prime workflow. Analysis run on the UArizona HPC
(Puma), September 2026. Biological findings are summarised separately in the results summary (not yet public).

All code excerpts below are copied from the package source (`src/campypump/`). Line-level details
can be checked against the repository.

---

## 1. Question and definitions

The tool answers one primary question per genome: **is the cmeB gene of the cmeABC efflux operon
the resistance-enhancing variant (RE-cmeB) or a wildtype allele?** Secondary outputs describe
the rest of the operon (cmeR, cmeA, cmeC), the CmeR-binding site in the cmeR–cmeA intergenic
region, and evidence of mosaic (recombinant) cmeB alleles.

| Term | Operational definition in campypump |
|---|---|
| RE-cmeB | cmeB with ≥ 95 % global nucleotide identity to the RE-cmeB of *C. jejuni* PA001. Reported alongside the published definition (> 90 % identity to *C. coli* DH161 cmeB, GenBank KT778508.1); the two agree on every Peru genome. |
| Wildtype cmeB | cmeB with ≤ 85 % identity to PA001 RE-cmeB. These alleles are typically 84–100 % identical to NCTC 11168 cmeB. |
| % identity | Geneious-equivalent global identity: identical columns / alignment length, gaps included (Section 5.3) |
| Mosaic | ≥ 10 % of 300-bp windows clearly RE-like **and** ≥ 10 % clearly wildtype-like (Section 5.5) |
| Partial gene | cmeB interrupted by a contig edge, or < 95 % of reference length |
| Divergent | complete cmeB < 90 % identical to every reference allele in the panel |
| Species | assigned by whole-genome ANI (skani) against one reference genome per species, ≥ 95 % (Section 5.10). **Not** assigned from the cme operon, which crosses species |

The operational definition is nucleotide identity to a reference RE allele. The tool does **not**
test individual amino-acid substitutions reported for RE-CmeB, and it does **not** measure
expression or phenotype.

---

## 2. Software and dependencies

| Component | Version tested | Role |
|---|---|---|
| Python | ≥ 3.9 (3.11 on HPC; 3.12 in development) | runtime |
| Biopython | ≥ 1.80 (1.88 in development) | FASTA parsing, translation, pairwise alignment (`Bio.Align.PairwiseAligner`) |
| BLAST+ | ≥ 2.10 (2.12.0 in development) | locating operon genes (`blastn`) |
| bigsdb-downloader | 1.0.6 | OAuth-authenticated PubMLST API access (private records) |
| skani | >=0.2 | whole-genome ANI for the species call; campypump runs it against the four bundled reference genomes |
| MAFFT / IQ-TREE | 7.526 / 3.1.3 | alignment and phylogeny (run separately; campypump prepares the input and draws the tree) |
| pandas + openpyxl | optional | reading Excel id lists |

The package installs with `pip install .` into a conda/micromamba environment that provides
BLAST+. A draft bioconda recipe is included (`recipe/meta.yaml`).

Package layout:

```
src/campypump/
  core.py        gene finding, ORF extraction, alignment, classification, promoter
  cli.py         command line: run, fetch, metadata, annotate, combine
  plots.py       per-allele segment plots (plain SVG)
  figures.py     classification and allele-diversity figures (plain SVG)
  species.py     whole-genome ANI species calls from a skani table
  alleles.py     deduplicated alleles + annotation for phylogenetics
  treeplot.py    annotated tree figures from a Newick file
  palettes.py    the one colour palette every figure uses
  pubmlst.py     PubMLST isolate metadata (MLST + provenance) via REST API
  annotate.py    merge results + metadata, duplicate and species checks
  data/references.fasta   bundled reference panel
tests/           regression tests on 8 reference genomes + unit tests (25 tests)
hpc/             Slurm job scripts and HPC guide
```

---

## 3. Inputs and outputs

**Input:** genome assemblies in FASTA format (contigs), one file per isolate. Files downloaded
from PubMLST are named `id-<PubMLST id>.fasta`, and the id is parsed from the file name.

**Main output:** a tab-delimited table with one row per genome. Key columns:

| Column | Content |
|---|---|
| `cmeB_type` | call from the full gene: RE-cmeB, wildtype-cmeB, ambiguous, mosaic, incomplete, not_found |
| `cmeB_type_final`, `re_cmeB` | type used in analyses (partial genes receive a provisional call); 1/0 coding |
| `cmeB_identity_RE_PA001_pct`, `cmeB_identity_WT_11168_pct` | Geneious-equivalent % identity to the two primary references |
| `cmeB_identity_DH161_pct`, `RE_by_DH161_definition` | % identity to the published DH161 RE-cmeB and the > 90 % published call |
| `cmeB_closest_reference`, `cmeB_closest_identity_pct`, `panel_*_pct` | identity to every panel allele |
| `cmeB_RE_like_segments_pct`, `cmeB_WT_like_segments_pct` | window scan summary |
| `cmeB_aligned_identity_RE_pct`, `..._WT_pct` | identity over the part of the gene present (partial genes) |
| `cmeR/cmeA/cmeC_identity_WT_11168_pct`, `..._RE_PA001_pct` | per-gene identities |
| `promoter_length_bp`, `CmeR_binding_site_seq`, `CmeR_binding_site_as_11168` | promoter analysis |
| `cmeB_contig`, `cmeB_start`, `cmeB_end`, `cmeB_strand`, `cmeB_length_bp` | locus coordinates |
| `warnings` | plain-language QC flags |

**Optional outputs:** a FASTA file of all extracted cmeB alleles (`--alleles`), SVG/TSV segment
plots (`--segment-plots`), and an annotated CSV merged with PubMLST metadata (`annotate`).

### 3.1 Every output file

**`campypump run`**

| File | Format | Content |
|---|---|---|
| `<output>.tsv` (`-o`) | tab-delimited, one row per genome | all results columns (Section 3 and the column guide) |
| `<output>_column_guide.txt` | text | plain-language description of every column, written on every run |
| `--alleles <file>.fasta` | FASTA | extracted cmeB ORF per genome; header `>sample cmeB_type_final`. Input for trees. |
| `--segment-plots <dir>/<sample>_cmeB_segments.svg` | SVG image | identity to PA001 (red) and best wildtype reference (blue) along cmeB |
| `--segment-plots <dir>/<sample>_cmeB_segments.tsv` | tab-delimited | the plotted values: `position_bp`, `identity_RE_PA001_pct`, `identity_best_wildtype_pct` (window midpoints) |

Segment plots are written for alleles that are mosaic, ambiguous, unresolved or divergent, or that
have both ≥ 5 % RE-like and ≥ 5 % wildtype-like windows. `--plot-all` plots every genome. A one-line
summary of counts per type is printed to the log.

**`campypump combine`**

| File | Content |
|---|---|
| `<output>.tsv` | result tables stacked (header once); columns are the union of the inputs |

**`campypump fetch`** (public genomes)

| File | Content |
|---|---|
| `<outdir>/id-<PubMLST id>.fasta` | contigs as served by PubMLST (headers are sequence-bin ids) |
| `<outdir>/failed_ids.txt` | `id <tab> reason` for every id not downloaded in the last run (e.g. `HTTP 404 (not public / not available)`) |
| `<outdir>/id-<id>.fasta.part` | temporary file during download; renamed on success, never left behind on success |

**`hpc/download_private.slurm`** (private genomes, BIGSdb_downloader)

| File | Content |
|---|---|
| `genomes/<species>/id-<id>.fasta` | contigs of private records |
| `genomes/<species>/private_failed_ids.txt` | ids that failed in the last run (empty when all succeeded) |
| `download_private_<jobid>.log` | progress every 100 genomes, one line per failure, `DONE` summary per species |

**`campypump metadata`**

| File | Content |
|---|---|
| `<output>.tsv` | one row per isolate: `pubmlst_id`, `isolate`, `pubmlst_species`, `ST`, `clonal_complex`, `year`, `source`, then all other provenance fields (Stata-safe names), `aliases`, `assembly_contigs`, `assembly_length_bp`, `mlst_loci_designated` |
| `<output>.cache.jsonl` | one JSON record per retrieved isolate. This makes the command resumable; delete it to force a full re-download. |
| `<output>.failed.txt` | `id <tab> reason` for records not retrieved (only written if any failed) |

**`campypump figure` / `campypump diversity`**

| File | Content |
|---|---|
| `figures/cmeB_classification_<species>.svg` | per genome: identity to the RE and wildtype references (A) and the distribution of identity to the RE reference (B) |
| `<output>.svg` (diversity) | allele rank abundance (A), rarefaction (B), distinct alleles per cmeB type (C) |
| `--table <file>.tsv` | one row per distinct allele: rank, genomes carrying it, length, cmeB type, example genome |

**`campypump alleles`** (input for phylogenetics)

| File | Content |
|---|---|
| `<prefix>.fasta` | one record per distinct allele, Newick-safe ids (`Cc_a003_RE_n63`), plus any `--references` and `--outgroup` records |
| `<prefix>.tsv` | per allele: species, type, genomes, top STs, top sources, mean wildtype-like segments, block status, example genome |
| `<prefix>_itol_type.txt`, `_itol_species.txt`, `_itol_genomes.txt` | iTOL colour strips and a genomes-per-allele bar chart |

**`campypump annotate`**

| File | Content |
|---|---|
| `<output>.csv` | results merged with metadata, ANI species calls and exclusion flags. `--columns core` (default) writes the 39 columns analyses use; `--columns full` writes all 82. Documented column by column in `docs/codebook_annotated_csv.md` |

**HPC job scripts** (`hpc/`, project folder `/xdisk/cooperma/fschiaffino/campypump_project`)

| File | Written by | Content |
|---|---|---|
| `campypump_<jobid>.log` | `campypump_job.slurm` | download summary, pending counts, type counts per species, timestamps |
| `results/campypump_ccoli.tsv`, `results/campypump_cjejuni.tsv` | `campypump_job.slurm` | per-species results (+ column guides) |
| `results/campypump_all_peru.tsv` | `campypump_job.slurm` | both species combined |
| `results/cmeB_alleles_ccoli.fasta`, `results/cmeB_alleles_cjejuni.fasta` | `campypump_job.slurm` | extracted cmeB alleles |
| `results/segment_plots_ccoli/`, `results/segment_plots_cjejuni/` | `campypump_job.slurm` | segment plots |
| `results/pending_ccoli_ids.txt`, `results/pending_cjejuni_ids.txt` | `campypump_job.slurm` | ids still not downloaded (empty once complete) |
| `metadata_<jobid>.log`, `results/pubmlst_metadata.tsv` | `metadata_job.slurm` | metadata retrieval |
| `results/campypump_all_peru_annotated.csv` | `metadata_job.slurm`, refreshed by `campypump_job.slurm` | analysis-ready file |

---

## 4. Reference panel

`data/references.fasta` holds 15 sequences. Headers follow `gene|class|strain|role`:

| Gene | Class | Strain | PubMLST id | Role | Length (bp) |
|---|---|---|---|---|---|
| cmeR, cmeA, cmeB, cmeC | *C. jejuni* wildtype | NCTC 11168 (AL111168.1) | 48 | primary WT | 633 / 1,104 / 3,123 / 1,479 |
| cmeR, cmeA, cmeB, cmeC | *C. jejuni* RE | PA001 (ST-10288) | 105572 | primary RE | 633 / 1,104 / 3,120 / 1,479 |
| cmeB | *C. jejuni* wildtype | EF0003911-007703 (ST-10315) | 150796 | panel | 3,123 |
| cmeB | *C. jejuni* wildtype | CAM013653-067612 (ST-11962) | 151401 | panel | 3,123 |
| cmeB | *C. jejuni* RE | CAM015192-073790 (ST-1036) | 151430 | panel | 3,120 |
| cmeB | *C. coli* wildtype | CAA03417-009068 | 152083 | panel | 3,123 |
| cmeB | *C. coli* wildtype | CAA03009-008347 | 151856 | panel | 3,123 |
| cmeB | *C. coli* RE | CAA03704-009522 | 152096 | panel | 3,120 |
| cmeR, cmeA, cmeB, cmeC | *C. coli* RE (published) | DH161 (GenBank KT778508.1) | — | panel | 633 / 1,104 / 3,120 / 1,479 |
| promoter | *C. jejuni* wildtype | NCTC 11168 | 48 | primary | 95 |

**How the panel was built.** The cmeR/A/B/C sequences in CARD (as distributed with abricate) were
used once as seed queries to locate the operon in the eight assemblies above, which the project
PI chose as representative isolates. Each gene was then extracted as a full ORF with the same
function used at runtime (Section 5.2). Every sequence starts with ATG, ends with a stop codon,
and has no internal stop. The CARD seed is not used at runtime.

The promoter reference is the 95-bp region upstream of cmeA in NCTC 11168, in cmeA orientation.
It contains the 16-bp inverted repeat `TGTAATAAAAATTACA`, which campypump treats as the
CmeR-binding site.

**Why PA001 is the RE reference.** Identity values computed against PA001 reproduce the
laboratory's Geneious values exactly (Section 9), so PA001's RE-cmeB is the allele the laboratory
used, or an identical one.

**Why a *C. jejuni* allele is the primary RE reference.** The *C. coli* RE allele CAA03704-009522
is in the panel, and every genome's identity to it is reported (`panel_CcRE_CAA03704_pct`). It is
not used for the call, for three reasons:

- PA001 reproduces the laboratory's Geneious values, which keeps results continuous with earlier
  work.
- CAA03704 itself carries a roughly 300-bp wildtype-like block (Section 5.5). As the yardstick, it
  would hide that block in related alleles.
- Cooper et al. (2025) inferred a *C. jejuni* origin for these alleles.

Using the best of the three RE alleles instead of PA001 would not change any call: every *C. coli*
RE allele is ≥ 98.6 % identical to at least one RE reference, and ≥ 97.4 % to PA001.

**The published RE-cmeABC sequence is in the panel (since 0.1.4).** RE-cmeABC was first detected
in ***C. coli* DH161*** (Yao et al. 2016), GenBank KT778508.1. Cooper et al. (2025) define RE-cmeB
as > 90 % nucleotide identity to DH161 cmeB. The complete DH161 operon was added to the panel with
the same ORF extraction used at runtime, and two columns report the published definition directly:

- `cmeB_identity_DH161_pct`
- `RE_by_DH161_definition` (yes when > 90 %). For partial genes this uses identity over the part
  present (Section 5.6), since a fragment's whole-gene identity is deflated by the missing section.

KT778507.1 (strain NT161, the transformant that received DH161 cmeABC) has an operon sequence
**identical** to KT778508.1, so only KT778508.1 was added.

**The two definitions agree on the Peru data (verified on the full run).** Across all 3,733
non-duplicate genomes, the campypump call and the published DH161 rule classify every RE and every
wildtype genome identically; the single mosaic allele is the only isolate they treat differently. DH161 cmeB is 99.04 % identical to PA001 and
99.17 % to CAM015192, so the campypump RE class and the published class coincide: every Peru RE
allele is ≥ 97.4 % identical to PA001, far above the 90 % DH161 threshold, and every wildtype
allele is ≤ 80.2 % to PA001, far below it. The only genome between the two is the mosaic
(91.7 % to PA001), which is > 90 % to DH161 and therefore RE under the published definition while
campypump reports it as mosaic. That difference is deliberate: the mosaic flag carries more
information than a single identity cut-off.

In detail:

**DH161 groups with the *C. jejuni* RE alleles.** At the protein level, DH161 CmeB differs from
PA001 CmeB at 11 of 1,039 residues and from CAM015192 at 9, but from the *C. coli* RE allele
CAA03704 at 13, and it is 97.76 % identical to CAA03704 at the nucleotide level versus 99.04–99.17 %
to the two *C. jejuni* RE alleles. All RE alleles, including DH161, are 1,039 aa, whereas wildtype
alleles are 1,040 aa. This is consistent with the *C. jejuni* origin inferred by Cooper et al.
(2025), even though DH161 itself is *C. coli*.

---

## 5. Per-genome algorithm

### 5.1 Locating the operon genes (BLAST)

All cme gene references (both primary and panel) are searched against the assembly in a single
`blastn` call. For each gene, the single highest-bitscore HSP across all reference alleles is
kept.

```python
out = subprocess.run(
    ["blastn", "-task", "blastn", "-query", qfa, "-subject", genome, "-evalue", "1e-20",
     "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send bitscore"],
    capture_output=True, text=True, check=True).stdout
best: dict[str, tuple] = {}
for line in out.splitlines():
    f = line.split("\t")
    ref = refs[int(f[0])]
    bits = float(f[9])
    if ref.gene not in best or bits > best[ref.gene][-1]:
        best[ref.gene] = (ref, f[1], int(f[5]), int(f[6]), int(f[7]), int(f[8]), bits)
```

Design notes:

- **`-task blastn`, not the default megablast.** During development, megablast failed to detect
  the divergent *C. coli* cmeR in CAA03009-008347. The more sensitive task finds it.
- **One copy assumed.** The single best HSP per gene means the operon is assumed to be
  single-copy. A gene split across two contigs contributes only its best fragment, and is then
  flagged as partial.
- **Paralogs.** The e-value threshold (1e-20) would also admit hits to other RND transporters,
  such as the CmeDEF system. The best-bitscore rule selects the true cmeB whenever it is present.
  If cmeB were absent, a paralog hit would be reported, but it would be caught as low identity
  (`divergent`). See the limitations section.

### 5.2 ORF extraction and QC

The gene is reconstructed from the BLAST hit by projecting the reference start codon onto the
contig, then reading in frame to the first stop codon.

```python
def extract_orf(contig_seq: Seq, ref: Ref, qs: int, ss: int, se: int) -> GeneHit:
    """Anchor on the reference start codon and extend to the first in-frame stop."""
    hit = GeneHit(gene=ref.gene)
    L = len(contig_seq)
    if ss < se:
        strand, s = "+", contig_seq
        start = (ss - 1) - (qs - 1)
    else:
        strand, s = "-", contig_seq.reverse_complement()
        start = (L - ss) - (qs - 1)
    if start < 0:
        hit.flags.append("contig_edge_5prime")
        start = start % 3  # keep frame, start at contig edge
    i, stop_found = start, False
    while i + 3 <= len(s):
        if str(s[i:i + 3]).upper() in STOPS and i > start:
            stop_found = True
            break
        i += 3
    if not stop_found:
        hit.flags.append("contig_edge_3prime")
    ...
    if orf[:3] not in STARTS:
        hit.flags.append("no_start_codon")
```

After extraction, a gene shorter than 95 % of its reference is flagged `truncated`:

```python
if len(h.seq) < p.min_len_frac * len(ref.seq):   # min_len_frac = 0.95
    h.flags.append("truncated")
```

Things to verify:

- **Start projection.** The start position is projected as `hit start − (query start − 1)`,
  which assumes no indel between the true start codon and the beginning of the HSP. If there is
  one, the ORF can start out of frame, which is caught by `no_start_codon` or `truncated`. In the
  Peru data, complete cmeB alleles were 3,111, 3,120, 3,123 or 3,126 bp. All are multiples of 3
  and consistent with the references.
- **Premature stop codons.** Reading stops at the first in-frame stop. A premature stop (or a
  frameshift from an assembly error) therefore shows up as a short ORF, flagged `truncated` and
  classified as `incomplete`. The separate `internal_stop` check is effectively never triggered,
  because the ORF always ends at its first stop.

### 5.3 Percent identity (Geneious-equivalent)

Each extracted cmeB is globally aligned to every cmeB reference with the scoring used by Geneious
Prime's default global alignment (65 % similarity matrix, match 5 / mismatch −4, gap open 12,
gap extend 3). Identity is identical columns divided by alignment length, with gap columns
counted.

```python
def geneious_aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score = 5
    al.mismatch_score = -4
    al.open_gap_score = -12
    al.extend_gap_score = -3
    return al

def pct_identity(ref: str, query: str) -> tuple[float, object]:
    """Geneious-style global % identity: identical sites / alignment length (gaps counted)."""
    aln = _ALIGNER.align(ref, query)[0]
    c = aln.counts()
    length = c.identities + c.mismatches + c.gaps
    return 100.0 * c.identities / length, aln
```

**How this definition was chosen.** During development, several definitions were compared
against the laboratory's Geneious values for three *C. coli* genomes. These were: BLAST pident,
identity/(alignment length − gaps), identity/reference length, and global identity under three
scoring schemes. Only the Geneious scoring with gap-inclusive alignment length reproduced all six
values (two references × three genomes) to the third decimal (Section 9).

**Tie-breaking.** Only the first optimal alignment is used. Equally scoring alternative gap
placements can shift the identity by a few hundredths of a percent, and by up to 0.33 points in
one allele with a 12-bp in-frame deletion.

### 5.4 Classification

```python
if any(f in b.flags for f in ("truncated", "contig_edge_5prime", "contig_edge_3prime", "internal_stop")):
    call = "incomplete"
elif "mosaic" in flags:
    call = "mosaic"
elif pid_re >= p.re_min:          # 95.0
    call = "RE-cmeB"
elif pid_re <= p.wt_max_re:       # 85.0
    call = "wildtype-cmeB"
else:
    call = "ambiguous"
```

The order of these checks matters. Incomplete genes are never classified from full-length
identity, because a missing section lowers global identity. Mosaic takes precedence over RE and
wildtype. The call depends only on identity to PA001; identity to 11168 is reported but not used
for the call.

**Why the call uses PA001, not 11168.** The laboratory's original rule was wildtype if > 96 % to
11168, RE if < 80 % to 11168. Applied to the 1,513 *C. coli*, it left 1,207 (80 %) unclassified,
because native *C. coli* cmeB alleles are 84–92 % identical to 11168. It also missed 42 RE
alleles that score 80–81 % to 11168. Identity to PA001, by contrast, separates the classes with a
wide empty gap:

| Species | RE-cmeB, % to PA001 | Wildtype, % to PA001 | Isolates between 80.2 % and 97.4 % |
|---|---|---|---|
| *C. coli* (n = 1,507 complete) | 97.4–99.9 | 76.6–80.2 (plus one divergent allele at 72.6) | 0 |
| *C. jejuni* (n = 2,205 complete) | 97.6–100 | 77.4–79.7 | 1 (the mosaic, 91.7 %) |

The 95 % and 85 % cut-offs sit inside this gap, and the 85–95 % zone is reported as `ambiguous`.
No complete non-mosaic gene in either species fell in that zone.

### 5.5 Window scan for mosaic alleles

The global alignments from 5.3 are converted to a per-position vector for the query: 1 where the
query base matches the aligned reference base, 0 otherwise.

```python
def query_match_vector(aln) -> list[int]:
    """1/0 per query position: identical to the aligned reference base or not."""
    ref_idx, q_idx = aln.indices
    ref, q = aln.sequences
    vec = [0] * len(q)
    for ri, qi in zip(ref_idx, q_idx):
        if qi >= 0 and ri >= 0 and ref[ri] == q[qi]:
            vec[qi] = 1
    return vec
```

Then, in 300-bp windows with a 50-bp step along the query, identity to PA001 is compared with the
**best** identity among all wildtype references in the panel (*C. jejuni* and *C. coli*):

```python
for w0 in range(0, max(n - p.window, 0) + 1, p.step):
    w1 = min(w0 + p.window, n)
    re_id = 100 * sum(re_vec[w0:w1]) / (w1 - w0)
    wt_id = max(100 * sum(v[w0:w1]) / (w1 - w0) for v in wt_vecs)
    total += 1
    if windows is not None:
        windows.append(((w0 + w1) // 2, re_id, wt_id))
    if re_id - wt_id >= p.window_diff:        # 5 percentage points
        re_like += 1
    elif wt_id - re_id >= p.window_diff:
        wt_like += 1
```

A gene is flagged `mosaic` when both `re_like/total` and `wt_like/total` are at least 0.10.

Review points:

- **Heuristic, not a formal test.** This is a screening heuristic, not a statistical
  recombination test; there is no breakpoint inference and no significance. For publication,
  candidate alleles should be confirmed with dedicated methods (for example RDP, GARD, or
  ClonalFrameML/Gubbins on the operon alignment).
- **RE side uses only PA001.** The RE-cmeB panel alleles (CAM015192, CAA03704) are not included
  in the RE side of the window comparison, while the wildtype side takes the maximum over all
  wildtype alleles. This makes the scan slightly conservative towards calling RE-like windows.
- **Heterogeneous wildtype side.** Because the wildtype side includes *C. coli* wildtype alleles,
  a window "closer to wildtype" may mean closer to a *C. coli*-type allele or to a *C. jejuni*
  wildtype allele. The segment plots show the best wildtype identity without saying which
  allele supplied it.
- **RE subgroups in *C. coli*.** Of the 211 *C. coli* RE alleles, 142 are closest to the
  *C. coli* RE reference CAA03704 and 69 to a *C. jejuni* RE reference (41 PA001, 28 CAM015192).
  All 63 *C. coli* RE alleles with ≥ 5 % wildtype-like windows are ≥ 99 % identical to CAA03704.
  This suggests the wildtype-like block is shared by one *C. coli* RE subgroup, consistent with a
  single recombination event followed by spread, rather than repeated independent events. This
  should be tested phylogenetically.
- **Example.** The *C. coli* RE reference CAA03704-009522 (97.4 % to PA001) shows a block at
  about 1,100–1,400 bp where it is about 96 % identical to a wildtype reference and about 87 % to
  PA001. It has 5.3 % wildtype-like windows, below the 10 % mosaic threshold, so it is not flagged
  mosaic but is plotted. In the full Peru data, 63 of 211 *C. coli* RE alleles (and 12 of 1,072
  *C. jejuni* RE alleles) had ≥ 5 % wildtype-like windows. All 63 *C. coli* alleles are ≥ 99 %
  identical to one another and 59 of them belong to a single sequence type (ST-12298), so the
  screen is picking out one lineage rather than scattered noise. This is the main evidence to
  review for partial recombination in *C. coli* RE-cmeB.

### 5.6 Partial genes: provisional call and aligned identity

Partial genes (`incomplete`) get a provisional `cmeB_type_final` when at least 90 % of the
windows over the part that is present agree:

```python
partial = call == "incomplete"
if partial:
    if re_frac >= p.partial_min_frac:          # 0.90
        final = "RE-cmeB"
    elif wt_frac >= p.partial_min_frac:
        final = "wildtype-cmeB"
    else:
        final = "unresolved"
```

A second identity is also reported, computed only over the reference region covered by the
fragment. It uses a semi-global alignment in which terminal gaps in the query are free:

```python
def aligned_identity(ref: str, query: str) -> float:
    aln = _SEMIGLOBAL.align(ref, query)[0]
    ref_idx, q_idx = aln.indices
    cols = [k for k, qi in enumerate(q_idx) if qi >= 0]
    ...
    for k in range(cols[0], cols[-1] + 1):
        ...
    return 100.0 * ident / length
```

For complete genes this equals the global identity. For fragments it gives a meaningful value:
a 780-bp 5′ fragment of 11168 cmeB scores 100 % aligned identity to 11168, versus 25 % global
identity. This appears to be how Geneious reports partial hits. For example, the laboratory's
Geneious value for the 780-bp fragment in PE1009159 was 95.58 % to 11168.

**Validation of the provisional calls.** In the Peru data, all 21 partial genes (5 *C. coli*,
16 *C. jejuni*) had ≥ 95 % of windows agreeing, and all 5 *C. coli* provisional calls matched the
laboratory's Geneious classification. The limitation is that a provisional call reflects only
the part of the gene that is present.

### 5.7 Divergent alleles

```python
divergent = (not partial) and ids[nearest.label] < p.divergent_max   # 90.0
```

A complete gene that is less than 90 % identical to every panel allele is flagged. In the Peru
data, one *C. coli* genome is flagged: CAM004289-088753 (PubMLST 164376), with a 3,126-bp cmeB
that is 72.6 % identical to PA001 and 73.7 % to 11168. Geneious gave the same values, so this is
a real, unusual allele rather than a tool artefact. It is counted as wildtype because it is
clearly not RE. We recommend checking its species assignment (rMLST) and checking whether it is
cmeB rather than a paralog, for example by translated identity to CmeB.

### 5.8 cmeR, cmeA and cmeC

These genes are extracted with the same functions as cmeB. They are globally aligned only to the
primary 11168 and PA001 versions:

```python
for gene in ("cmeR", "cmeA", "cmeC"):
    ...
    for r in refs:
        if r.gene == gene and r.role == "primary":
            row[col_re if r.name == PRIMARY_RE else col_wt] = round(pct_identity(r.seq, h.seq)[0], 3)
```

These values are descriptive and do not affect the cmeB call. There are no *C. coli* references
for cmeR/A/C, so in *C. coli* these identities measure distance from the *C. jejuni* versions.
For example, in *C. coli* RE isolates the median cmeA identity to 11168 is 94.1 %, versus 85.6 %
in *C. coli* wildtype isolates.

### 5.9 Promoter and CmeR-binding inverted repeat

The cmeR–cmeA intergenic region is extracted when cmeR and cmeA sit on the same contig and strand
in the expected order. Otherwise, the 150 bp upstream of cmeA is used. The region is put in cmeA
orientation and locally aligned to the 95-bp 11168 promoter reference. The query bases aligned to
the 16 IR positions, including any insertions within them, are reported:

```python
al = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3,
                           open_gap_score=-5, extend_gap_score=-2)
aln = al.align(promoter_ref, region)[0]
ir0 = promoter_ref.find(IR_11168)          # "TGTAATAAAAATTACA"
ir1 = ir0 + len(IR_11168)
ref_idx, q_idx = aln.indices
cols = [k for k, ri in enumerate(ref_idx) if ir0 <= ri < ir1]
c0, c1 = cols[0], cols[-1]
ir_seq = "".join(region[q_idx[k]] for k in range(c0, c1 + 1) if q_idx[k] >= 0)
```

Review points:

- **Which IR is used.** The IR is defined from the 11168 sequence (`TGTAATAAAAATTACA`, a
  6-bp arm, a 4-bp spacer, and the reverse-complement arm). Please confirm it matches the
  published CmeR-binding site (Lin et al. 2005).
- **Alignment artefacts.** Local alignment of a 16-bp motif can place gaps ambiguously, so
  reported IR sequences that include indels should be spot-checked. In *C. jejuni*, 17 distinct
  IR sequences were observed among 2,203 genomes with an IR call; the ten most frequent account
  for 98.3 %. No IR was reported for 20 *C. jejuni* and 54 *C. coli* genomes.
- **C. coli.** *C. coli* promoters diverge from 11168, so the mapped IR in *C. coli* is less
  reliable. `CmeR_binding_site_as_11168` only tests identity to the 11168 IR, not function.

### 5.10 Species assignment (whole-genome ANI)

**A single cme gene cannot identify the species.** The operon moves between *C. coli* and
*C. jejuni*, so identity at any one of its genes reflects the operon's history, not the genome's.
Concretely, in the Peru collection *C. coli* cmeR identity to NCTC 11168 spans 42–99.7 % and
*C. jejuni* spans 1.4–100 %, and the published *C. coli* strain DH161 carries a cmeR **99.4 %**
identical to 11168. Any threshold on a cme gene would call DH161 misfiled.

campypump therefore takes the species call from whole-genome ANI, which a single transferred locus
cannot move. ANI is computed outside campypump:

```bash
skani dist -q genomes/*/*.fasta -r species_refs/*.fasta -t 8 -o species_ani.tsv
```

with one reference genome per species — here *C. jejuni* NCTC 11168 (AL111168.1), *C. coli* RM4661
(CP007181.1), *C. lari* RM2100 (CP000932.1) and *C. upsaliensis* RM3940 (CP053849.1). The species is
read from the reference file name. `campypump annotate --ani` merges the table:

```python
ANI_SPECIES_MIN = 95.0        # conventional species boundary
ANI_MARGIN_MIN = 1.0          # best hit must beat the runner-up by this much
ANI_ALIGN_MIN = 60.0          # % of the query aligned; below this the call rests on little data
```

For each genome the best hit gives `species_by_ANI` and `ANI_pct`; `ANI_margin` is the gap to the
best hit of a *different* species; `ANI_call_quality` records why a call is weak (nothing above
95 %, low alignment fraction, or another species within 1 point); and `species_check_ANI` is `CHECK`
when the call disagrees with the submitted label.

**Separately, the operon signal** (`operon_species_signal`) scores cmeR, cmeA and cmeC against
11168, calling the operon for a species only when all three agree (each > 92 % ⇒ *C. jejuni*-like;
none ⇒ *C. coli*-like; otherwise mixed; any incomplete gene ⇒ not scored). `operon_atypical_for_label`
marks a disagreement with the label. **This is deliberately not a species call:** a disagreement
means either a database error or complete introgression of the operon, and the two cannot be
separated at this locus. Section 9.7 reports how each performed.

### 5.11 Exclusions and column sets

`annotate --exclude <file>` marks listed genomes `use_for_prevalence = no` with an `exclude_reason`;
nothing is deleted, so an exclusion is visible and reversible. `annotate --columns` writes either
the 39 columns analyses use (`core`, the default) or all 82 (`full`); every column is documented in
`docs/codebook_annotated_csv.md`.

### 5.12 Figures

All figures are written as plain SVG with no plotting dependency, from one palette in
`palettes.py`; editing that file recolours every output together. `campypump figure` draws the
classification panels, `campypump diversity` the allele-diversity panels, `campypump run
--segment-plots` the per-allele identity profiles, and `campypump treeplot` an annotated tree from a
Newick file.

`treeplot` colours each branch by what lies below it (cmeB type, or species inside the RE clade
where every tip is RE), shades single-species clades, and shows support only at or above a
threshold and only on clades of three or more tips, since cherries of near-identical alleles would
otherwise fill the figure with uninformative labels. With `--annotation` it adds the dominant ST per
allele, scales tips by the number of genomes carrying the allele, rings alleles carrying the
wildtype-like block, and draws a per-allele source-composition bar. `--cladogram` and
`--drop-outgroup` exist because RE and wildtype cmeB differ by roughly 0.4 substitutions per site,
which on a phylogram compresses everything else.

### 5.13 Phylogenetics input (`campypump alleles`)

campypump does not infer trees. It deduplicates genomes to distinct alleles, gives them Newick-safe
ids (`Cc_a003_RE_n63`: species, allele rank by frequency, cmeB type, genomes carrying it), and
writes an annotation table plus iTOL colour strips and a Microreact CSV. Alignment and inference are
done with MAFFT and IQ-TREE, so the phylogeny comes from established tools rather than a bespoke
implementation. The procedure used for the Peru collection is in `docs/tree_howto.md`.

Partial genes are excluded from the allele set: a gene interrupted by a contig edge is a fragment,
not a distinct allele, and counting it both inflates the allele count and distorts the alignment.

---

## 6. Data acquisition

### 6.1 Public assemblies (`campypump fetch`)

Contigs are downloaded anonymously from the PubMLST REST API:

```python
PUBMLST = "https://rest.pubmlst.org/db/pubmlst_campylobacter_isolates/isolates/{}/contigs_fasta"
```

- **Checks on each download.** A response is accepted only if it begins with `>`. Files are
  written to a temporary `.part` name and renamed when complete, so an interrupted download never
  looks finished.
- **Failures.** HTTP 401, 403 and 404 responses are treated as definitive (private or absent
  records) and are not retried. Other errors are retried three times with back-off.

In the Peru collection, all isolates with PubMLST id ≤ 151958 were public (1,600 genomes). The
remaining 2,136 are private and returned 404 to anonymous requests ("Isolate … does not exist").

### 6.2 Private assemblies (`hpc/download_private.slurm`)

Private records are downloaded with BIGSdb_downloader, using a personal OAuth token tied to the
PI's PubMLST account (created once with `--setup`). It uses the same `contigs_fasta` route and the
same check that each file starts with `>`. All 2,136 private genomes were downloaded with no
failures.

### 6.3 Metadata (`campypump metadata`)

For each isolate, the full isolate record (`GET /db/pubmlst_campylobacter_isolates/isolates/{id}`)
is retrieved, authenticated when `--token-dir` is given. From it the tool extracts:

- provenance fields
- assembly contig count and total length
- the MLST scheme fields (ST, clonal complex) taken from the scheme whose description is "MLST"

```python
for sch in rec.get("schemes", []):
    if str(sch.get("description", "")).strip().upper() == "MLST":
        for k, v in (sch.get("fields") or {}).items():
            row[safe_name(k)] = ...
        row["mlst_loci_designated"] = sch.get("loci_designated_count", "")
```

STs are PubMLST's own designations and are not recomputed by campypump. The retrieval can be
resumed: records are cached one per line in `<output>.cache.jsonl`.

---

## 7. Annotation (`campypump annotate`)

This step merges results and metadata on `pubmlst_id` and adds three columns:

- **`duplicate_of` / `use_for_prevalence`:** records whose isolate names are identical after
  upper-casing and removing leading zeros inside numbers are treated as duplicates. For example,
  `CAM000744-3795` and `CAM000744-003795` both normalise to `CAM744-3795`. The first record is
  kept; later ones get `use_for_prevalence = no`. In the Peru data this found two re-uploaded
  isolates, one per species, with identical cmeB results.
- **`species_check`:** compares the last word of the `--species` label with PubMLST's species
  field. Note that PubMLST's species field is submitter-provided, so agreement here does not mean
  the species is right — the ANI call (Section 5.10) is what does.
- **ANI and operon species columns** (Section 5.10), and the exclusion flags (Section 5.11).

```python
def norm_isolate(name: str) -> str:
    s = re.sub(r"\s+", "", str(name).upper())
    return re.sub(r"(?<![0-9])0+(?=[0-9])", "", s)
```

Duplicate detection is by name only. Genomes that are identical but carry different names are not
detected; a cgMLST distance check would catch those.

---

## 8. Parameters (defaults)

| Parameter | Default | Where | Meaning |
|---|---|---|---|
| BLAST task / e-value | blastn / 1e-20 | 5.1 | gene finding |
| `min_len_frac` | 0.95 | 5.2 | shorter ORF → truncated |
| Alignment scoring | 5 / −4 / −12 / −3, global | 5.3 | Geneious-equivalent identity |
| `re_min` | 95 % | 5.4 | RE call (identity to PA001) |
| `wt_max_re` | 85 % | 5.4 | wildtype call (identity to PA001) |
| `window` / `step` | 300 / 50 bp | 5.5 | window scan |
| `window_diff` | 5 points | 5.5 | window counted as RE- or wildtype-like |
| `mosaic_min_frac` | 0.10 | 5.5 | mosaic flag |
| `partial_min_frac` | 0.90 | 5.6 | provisional call for partial genes |
| `divergent_max` | 90 % | 5.7 | divergent flag |

`re_min`, `wt_max_re`, `window`, `step` and `mosaic_min_frac` can be changed on the command line
(`campypump run --help`).

---

## 9. Validation

### 9.1 Regression tests (25, all passing)

- **8 reference genomes** (the panel isolates). Expected calls are checked, and for the three
  *C. coli* genomes the exact Geneious identities are checked (tolerance 0.01):

| Genome | Geneious (RE / 11168) | campypump (RE / 11168) |
|---|---|---|
| CAA03009-008347 | 77.714 / 89.978 | 77.714 / 89.978 |
| CAA03417-009068 | 78.172 / 95.453 | 78.172 / 95.457 |
| CAA03704-009522 | 97.438 / 79.776 | 97.438 / 79.776 |

- **Unit tests:**
  - aligned identity of gene fragments
  - provisional call for a gene running off a contig end
  - isolate-name normalisation
  - PubMLST record parsing
  - the three-gene operon signal, including that one *C. jejuni*-like gene does not decide it
  - ANI parsing, margin, low-alignment-fraction quality flag, and the species check
  - exclusion handling, allele deduplication and Newick-safe names, figure structure, and that a
    single-genome histogram bin is still drawn
  - tree parsing, tip naming, ST annotation and the source/block layers of the tree figure
  - that no colour carries two meanings within a figure
  - that every module parses under the Python 3.9 and 3.11 grammars, not only the interpreter used
    to develop (a nested same-quote f-string, legal only from 3.12, broke 0.1.18 on 3.11)

### 9.2 Synthetic tests (development)

| Synthetic genome | Result |
|---|---|
| Chimera: RE 5′ half + wildtype 3′ half | `mosaic` (46 % RE-like, 49 % wildtype-like windows). Its global identity was 89.5 % to PA001 and 90.0 % to 11168, which is indistinguishable from a native *C. coli* allele by identity alone. |
| Gene cut at a contig end | `incomplete`, provisional wildtype, 100 % aligned identity |
| Genome with no cmeB | `not_found` |
| Simulated 404 from PubMLST | skipped without retry, listed in `failed_ids.txt` |

### 9.3 Full concordance with Geneious: 1,513 Peru *C. coli*

- **Classification:** 100 % agreement. There were 211 RE and 1,297 wildtype. The 5 partial genes
  were all non-RE in Geneious, and all were provisional wildtype in campypump.
- **Identity to PA001:** 1,497 of 1,508 complete genes within 0.01 points, 1,507 within 0.1, and
  a maximum difference of 0.33 (151308, a 3,111-bp allele).
- **Identity to 11168:** 1,498 within 0.01, and all but one within 0.5. The exception is
  **CAA04817-012617 (164310)**. Geneious lists it as 97.44 % to PA001 (reproduced exactly) but
  94.97 % to 11168, whereas campypump finds 79.36 %. An allele 97 % identical to PA001 cannot be
  95 % identical to 11168, since PA001 and 11168 cmeB are only 79.6 % identical, so the Geneious
  entry is most likely a transcription error. The isolate is RE by both methods.

### 9.4 *C. jejuni* (no Geneious comparison available)

- The four *C. jejuni* panel isolates in the Peru set returned their expected calls.
- The classes are fully separated (Section 5.4). Excluding duplicate records there were 1,072 RE,
  1,148 wildtype and 1 mosaic, including 16 partial genes all resolved provisionally.

### 9.5 Performance

Typing all 3,736 genomes took 13 min 9 s on 8 CPUs (Puma), about 1.7 CPU-seconds per genome.

### 9.6 Final run totals

Three duplicate records (the same isolate uploaded twice under slightly different names) are
excluded; PubMLST species matched the analysis label for all 3,733 remaining genomes, and cmeB was
found in every one.

Grouped by the ANI species call, excluding duplicate records and one mixed assembly:

| | Genomes | RE-cmeB | Wildtype | Mosaic | % RE |
|---|---|---|---|---|---|
| *C. coli* | 1,516 | 208 | 1,308 | 0 | 13.7 % |
| *C. jejuni* | 2,206 | 1,075 | 1,130 | 1 | 48.8 % |
| *C. upsaliensis* | 1 | 0 | 1 | 0 | — |

Twenty-one genes were interrupted by a contig edge and carry provisional calls. Nine genomes have
no ANI species call (94.3–94.9 % to *C. jejuni*, just below the 95 % boundary) and are unresolved.

### 9.7 Species assignment: how the two methods performed

The Peru collection contained 13 records whose submitted species was wrong, all confirmed
independently by rMLST through PubMLST: nine *C. coli* filed as *C. jejuni*, three *C. jejuni* filed
as *C. coli*, and one *C. upsaliensis* filed as *C. coli*.

| Method | Flagged | True errors found | False positives |
|---|---|---|---|
| Whole-genome ANI (`species_check_ANI`) | 13 | 13 / 13 | 0 |
| Operon signal (`operon_atypical_for_label`) | 13 | 10 / 13 | 3 |

All 13 ANI calls were confirmed independently by rMLST. A 14th genome (80361) is flagged by ANI with
only 43.5 % of the assembly aligning — a mixed assembly — and is excluded rather than assigned.

ANI recovered every error with no false positives, and the separation is wide: the median genome is
9.2 percentage points closer to its own species than to the next, with three genomes below a
2-point margin.

The operon signal missed one (a mixed 1-of-3 signal) and flagged three genomes that are genuine
*C. coli* by both ANI and rMLST but whose **entire operon is *C. jejuni*-like** (115954, 115976,
152066). Those three are a biological result rather than an error, and they are the reason the
operon columns are not used for the species call.

One further genome (80361) is flagged by ANI with a 43.5 % alignment fraction against 81–92 % for
the rest — a mixed or contaminated assembly. It is excluded via `annotate --exclude` rather than
assigned a species.

---

### 9.8 Phylogeny (Peru collection)

MAFFT aligned the 391 sequences into 3,130 columns, close to the longest input gene (3,126 bp), so
no sequence forces large gaps. IQ-TREE selected GTR+F+I+R4 by BIC; log-likelihood −49,755.9;
17.5 minutes on 8 CPUs. Ultrafast bootstrap did not converge (split-frequency correlation 0.980
against a 0.99 target), which is expected where many alleles differ by one or two bases.

The RE clade is recovered with 100/100 support and contains no wildtype allele, so the 95 % identity
cut-off used for classification corresponds to a clade rather than an arbitrary threshold. Within
it, the *C. coli* alleles are not monophyletic: they form 16 separate clusters, which is the minimum
number of independent acquisitions. Findings are in the results summary (not yet public).

## 10. Known limitations and suggested reviewer checks

1. **Species assignment depends on the reference set.** ANI can only call species for which a
   reference genome is supplied; nine genomes at 94.3–94.9 % to *C. jejuni* fall below the 95 %
   boundary and are left unresolved. Adding references (e.g. *C. jejuni* subsp. *doylei*) would
   resolve some. Suggested check: rMLST on the unresolved genomes.
2. **The operon signal is not a species test** (Section 5.10) and is reported only as an observation
   about the operon.
3. **The phylogeny is a single-gene tree** (Section 9.8). It resolves the RE clade and the position
   of *C. coli* alleles within it, but not the wildtype backbone, and it cannot date transfers.
   "At least 16 acquisitions" is a minimum, since a transferred allele that later became identical
   to an existing one would not be counted separately.
4. **Definition of RE.** The RE call rests on nucleotide identity to PA001, but identity to the
   published DH161 allele is reported and agrees on every Peru genome (Section 4). campypump does
   not test the individual amino-acid substitutions shown to be responsible for the functional
   gain (PNAS 2024). Suggested extension: score those specific residues in RE-CmeB.
5. **One *C. jejuni* RE reference for both species.** *C. coli* RE alleles are 97.4–99.9 %
   identical to PA001, so they are clearly the same variant. But a hypothetical *C. coli*-native
   RE-like lineage, divergent from PA001, would be called wildtype or ambiguous. The ambiguous
   zone (85–95 %) was empty in these data.
6. **Mosaic detection is heuristic** (Section 5.5). Confirm candidates with formal recombination
   analysis on an alignment of the operon.
7. **Single best hit per gene.** Assumes one operon copy. A second copy, or a paralog in a genome
   lacking cmeB, would not be reported separately.
8. **ORF reconstruction.** Assumes the reference start codon projects correctly. Premature stops
   and frameshifts appear as `truncated` → `incomplete`, not as a separate "disrupted gene" call,
   so a real nonsense mutation would be treated as a partial gene. Suggested check: inspect the
   `warnings` column of all incomplete genes whose `cmeB_contig` is not at a contig edge.
9. **Assembly quality.** Contig breaks inside the operon produce partial genes (21 of 3,736 here).
   Their provisional calls rely on the part present.
10. **IR analysis.** Defined from the 11168 sequence; mapping is less reliable in divergent
   (*C. coli*) promoters. The IR is descriptive only.
11. **Lineage structure.** RE-cmeB is strongly clustered by ST and clonal complex (for example
   99 % of ST-1036 versus 0 % of ST-179 complex). Associations with source, time or the IR
   sequence should be modelled with lineage as a covariate or with phylogenetic methods.
12. **Duplicates.** Detected by name only (Section 7).
13. **Metadata.** Taken as recorded in PubMLST; STs are not recomputed.

---

## 11. Reproducing the Peru analysis

```bash
# environment (once)
micromamba create -n campypump -c conda-forge -c bioconda python=3.11 biopython blast pip pytest pandas openpyxl
micromamba activate campypump
pip install ./campypump        # campypump 0.1.3
pip install bigsdb-downloader  # private PubMLST access
pytest -q campypump/tests      # 12 passed

# genomes: public, then private (ids in the case-study id lists (not yet public), )
campypump fetch --ids peru_ccoli_ids.txt --outdir genomes/ccoli/
campypump fetch --ids peru_cjejuni_ids.txt --outdir genomes/cjejuni/
sbatch hpc/download_private.slurm          # remaining private records via BIGSdb_downloader

# typing
campypump run genomes/ccoli/   -t 8 --species "C. coli"   -o campypump_ccoli.tsv \
    --alleles cmeB_alleles_ccoli.fasta   --segment-plots segment_plots_ccoli
campypump run genomes/cjejuni/ -t 8 --species "C. jejuni" -o campypump_cjejuni.tsv \
    --alleles cmeB_alleles_cjejuni.fasta --segment-plots segment_plots_cjejuni
campypump combine campypump_ccoli.tsv campypump_cjejuni.tsv -o campypump_all_peru.tsv

# metadata
campypump metadata --ids peru_ccoli_ids.txt peru_cjejuni_ids.txt --token-dir ~/.bigsdb_tokens \
    -o pubmlst_metadata.tsv

# species by whole-genome ANI (one reference genome per species in species_refs/)
skani dist -q genomes/ccoli/*.fasta genomes/cjejuni/*.fasta -r species_refs/*.fasta \
    -t 8 -o species_ani.tsv

# annotation
campypump annotate campypump_all_peru.tsv --metadata pubmlst_metadata.tsv \
    --ani species_ani.tsv --exclude exclude_ids.txt \
    -o campypump_all_peru_annotated.csv

# figures and tree input
campypump figure campypump_all_peru_annotated.csv -o figures/ --exclude-duplicates
campypump diversity cmeB_alleles_ccoli.fasta cmeB_alleles_cjejuni.fasta \
    -o figures/cmeB_allele_diversity.svg --table figures/alleles.tsv
campypump alleles cmeB_alleles_ccoli.fasta cmeB_alleles_cjejuni.fasta \
    --results campypump_all_peru_annotated.csv --references tree/refs_cmeB.fasta \
    --outgroup tree/lari_cmeB.fasta -o tree/cmeB_alleles

# alignment and phylogeny (outside campypump)
mafft --auto --thread 8 tree/cmeB_alleles.fasta > tree/cmeB_alleles.aln.fasta
iqtree -s tree/cmeB_alleles.aln.fasta -m MFP -B 1000 -alrt 1000 -T 8 \
    --prefix tree/cmeB_tree -o OUT_C_lari_RM2100

# tree figures
campypump treeplot tree/cmeB_tree.treefile --outgroup OUT_C_lari_RM2100 --subtree RE \
    --annotation tree/cmeB_alleles.tsv -o figures/cmeB_tree_RE.svg
campypump treeplot tree/cmeB_tree.treefile --outgroup OUT_C_lari_RM2100 \
    --drop-outgroup --cladogram -o figures/cmeB_tree_all.svg
```

The Slurm scripts in `hpc/` run these steps on UArizona Puma (`campypump_job.slurm`,
`download_private.slurm`, `metadata_job.slurm`). Results reported in Section 9 were generated with
campypump 0.1.19. The classification logic has been unchanged since 0.1.2; later versions added the
partial-gene, divergent, plotting, metadata, annotation, ANI and phylogenetics features. All of it
runs in one step with `the case-study run script (not yet public)`.

---

## 12. References

Verify the citation details below before using them in a manuscript.

- Yao H, Shen Z, Wang Y, et al. Emergence of a potent multidrug efflux pump variant that enhances
  *Campylobacter* resistance to multiple antibiotics. *mBio* 2016;7(5):e01543-16.
- Mutation-based mechanism and evolution of the potent multidrug efflux pump RE-CmeABC in
  *Campylobacter*. *PNAS* 2024 (doi:10.1073/pnas.2415823121). This study shows that RE-CmeB amino
  acid substitutions, together with a promoter mutation, explain the functional gain, and it uses
  DH161 (KT778508.1) as the RE-cmeABC source.
- Lin J, Akiba M, Sahin O, Zhang Q. CmeR functions as a transcriptional repressor for the
  multidrug efflux pump CmeABC in *Campylobacter jejuni*. *Antimicrob Agents Chemother*
  2005;49:1067–1075.
- Schiaffino F, et al. Genomic resistant determinants of multidrug-resistant *Campylobacter* spp.
  isolates in Peru. *J Glob Antimicrob Resist* 2024;36:309–318.
- Cooper KK, et al. Sharing of *cmeRABC* alleles between *C. coli* and *C. jejuni* associated with
  extensive drug resistance in *Campylobacter* isolates from infants and poultry in the Peruvian
  Amazon. *mBio* 2025;16(2):e0205424.
- Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the
  PubMLST.org website and their applications. *Wellcome Open Res* 2018;3:124.
- Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani.
  *Nat Methods* 2023;20:1661–1665.
- Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7. *Mol Biol Evol*
  2013;30:772–780.
- Minh BQ, et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the
  genomic era. *Mol Biol Evol* 2020;37:1530–1534.
- Jolley KA, et al. Ribosomal multilocus sequence typing: universal characterization of bacteria
  from domain to strain. *Microbiology* 2012;158:1005–1015.
- Camacho C, et al. BLAST+: architecture and applications. *BMC Bioinformatics* 2009;10:421.
- Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology
  and bioinformatics. *Bioinformatics* 2009;25:1422–1423.
