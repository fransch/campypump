# Codebook: `campypump_all_peru_annotated.csv`

**File:** the annotated results table written by `campypump annotate` (core column set, 49 columns).
**Rows:** one per genome, 3,736 in the Peru collection.
**Produced by:** campypump 0.5.0 and later. Method: `docs/technical_description.md`.
**Full column set:** `--columns full` adds the remaining columns the run produced (46 in the Peru
run: per-panel identities, locus coordinates, accessions, full provenance). Everything omitted here is also in
`results/campypump_all_peru.tsv` and `results/pubmlst_metadata.tsv`.

Counts and ranges below are from the September 2026 run of the Peru collection, as a sanity check
when you regenerate the file.

---

## Before you analyse

Three columns decide which rows and which grouping to use:

1. **`use_for_prevalence == "yes"`** — drops duplicate records and excluded genomes (4 rows).
2. **Group by `species_by_ANI`, not `species`** — whole-genome ANI is the species call; the
   PubMLST label is wrong for 13 records in this collection.
3. **`re_cmeB`** is the 1/0 outcome variable; it is empty for mosaic and unresolved alleles.

```stata
import delimited "campypump_all_peru_annotated.csv", clear
keep if use_for_prevalence == "yes"
tab species_by_ani cmeb_type_final
```

---

## 1. Identity

| Variable | Type | Description |
|---|---|---|
| `pubmlst_id` | integer | PubMLST isolate id; the key for merging with PubMLST exports |
| `isolate` | text | isolate name as recorded in PubMLST (e.g. `CAA03704-009522`) |
| `sample` | text | genome file name without extension (`id-152096`); the key for merging with the results table and allele FASTAs |

## 2. Species (from whole-genome ANI)

No single cme gene identifies the species — the operon crosses between *C. coli* and *C. jejuni*.
These columns come from skani against reference genomes for *C. jejuni* (NCTC 11168, AL111168.1),
*C. coli* (RM4661, CP007181.1), *C. lari* (RM2100, CP000932.1) and *C. upsaliensis* (RM3940,
CP053849.1).

| Variable | Type | Values / range | Description |
|---|---|---|---|
| `species_by_ANI` | text | C. jejuni 2,208; C. coli 1,518; C. upsaliensis 1; empty 9 | **The species call to use.** Best-matching reference above 95 % ANI; empty when nothing reaches 95 % |
| `species_check_ANI` | text | consistent 3,714; CHECK 13; empty 9 | CHECK = ANI disagrees with the submitted species label |
| `ANI_pct` | numeric | 94.34–100 | ANI to the best-matching reference |
| `ANI_call_quality` | text | ok 3,709; "no reference above 95% ANI" 9; "low alignment fraction (n%)"; "close to another species" | Why a call is weak. A low alignment fraction means a fragmented or mixed assembly |

**Reading it.** 13 records are CHECK: 9 *C. coli* filed as *C. jejuni*, 3 *C. jejuni* filed as
*C. coli*, 1 *C. upsaliensis* filed as *C. coli* — all confirmed by rMLST. The 9 empty calls sit at
94.3–94.9 % to *C. jejuni*, just below the species boundary, and are unresolved.

## 3. The cmeB call

| Variable | Type | Values / range | Description |
|---|---|---|---|
| `cmeB_type_final` | text | wildtype-cmeB 2,451; RE-cmeB 1,283; mosaic 2 | **The result.** RE-cmeB = the resistance-enhancing variant. Partial genes carry a provisional call (see `cmeB_partial_gene`) |
| `re_cmeB` | 1/0 | 1 = 1,283; 0 = 2,451; empty 2 | RE-cmeB as a binary outcome; empty for mosaic, ambiguous and unresolved |
| `cmeB_partial_gene` | yes/no | yes 21 | The gene is interrupted by a contig edge; its type is provisional, from the part present |
| `cmeB_divergent` | yes/no | yes 1 | A complete gene under 90 % identical to every reference allele. In this collection the one `yes` is the *C. upsaliensis* record |

**Classification rule:** RE-cmeB at ≥ 95 % nucleotide identity to the RE-cmeB of *C. jejuni* PA001,
wildtype at ≤ 85 %, mosaic when RE-like and wildtype-like segments are both present. The classes do
not overlap in these data (RE ≥ 97.4 %, wildtype ≤ 80.2 %).

## 4. Evidence behind the call

| Variable | Type | Range | Description |
|---|---|---|---|
| `cmeB_identity_RE_PA001_pct` | numeric | 19.98–100 | % identity to the RE-cmeB of PA001 — the value the call is based on. Equals the laboratory's Geneious `% identity-cjREcmeB` |
| `cmeB_identity_WT_11168_pct` | numeric | 24.94–100 | % identity to the cmeB of NCTC 11168. Equals Geneious `% identity-cjcmeB11168` |
| `cmeB_closest_reference` | text | 9 panel alleles | The most similar allele in the reference panel, with its class in brackets |
| `cmeB_closest_identity_pct` | numeric | — | % identity to that allele |
| `cmeB_RE_like_segments_pct` | numeric | 0–100 | % of 300 bp windows along cmeB clearly closer to RE-cmeB |
| `cmeB_WT_like_segments_pct` | numeric | 0–100 | % of windows clearly closer to a wildtype allele. Both ≥ 10 % ⇒ mosaic |
| `cmeB_length_bp` | integer | 780–3,126 | Length of the extracted gene. Complete alleles are 3,111–3,126; shorter means partial |

**Note on low identities.** Values in the 20–50 % range belong to partial genes, where the whole-gene
identity is deflated by the missing section. Use `cmeB_partial_gene` to exclude them, or the
aligned-identity columns in the full set.

**The wildtype-like block.** An RE allele with `cmeB_RE_like_segments_pct` high but
`cmeB_WT_like_segments_pct` ≥ 5 carries a recombinant block — the signal that defines the *C. coli*
ST-12298 RE alleles.

## 5. Operon context

Identities to the two primary references for the other three operon genes. Descriptive; they do not
affect the cmeB call.

| Variable | Type | Description |
|---|---|---|
| `cmeR_identity_WT_11168_pct`, `cmeR_identity_RE_PA001_pct` | numeric | cmeR, the repressor, immediately upstream |
| `cmeA_identity_WT_11168_pct`, `cmeA_identity_RE_PA001_pct` | numeric | cmeA, the periplasmic fusion protein |
| `cmeC_identity_WT_11168_pct`, `cmeC_identity_RE_PA001_pct` | numeric | cmeC, the outer membrane channel |

Empty when the gene was not found (9 genomes for cmeR). In *C. coli*, identity to 11168 measures
distance from the *C. jejuni* version: median cmeR 75.6 %, cmeA 85.6 %, cmeC 84.0 % in wildtype
isolates, versus cmeA 94.1 % in RE isolates — the introgression signal.

## 6. Operon species signal

| Variable | Type | Values | Description |
|---|---|---|---|
| `operon_species_signal` | text | C. jejuni-like 2,144; C. coli-like 1,097; mixed 393; incomplete 102 | Agreement of cmeR, cmeA and cmeC against NCTC 11168 (each >92 % = *C. jejuni*-like). All three ⇒ C. jejuni-like; none ⇒ C. coli-like; 1 or 2 ⇒ mixed; any gene incomplete ⇒ incomplete |
| `operon_atypical_for_label` | yes/no/empty | yes 13; no 3,228; empty 495 | The operon does not match the labelled species. Empty when the signal is mixed or incomplete |

**`operon_atypical_for_label` is not a species call.** A `yes` means either a database error or
complete introgression of the operon. In this collection, 10 of the 13 were misfiled records and 3
(115954, 115976, 152066) are genuine *C. coli* whose entire operon is *C. jejuni*-like. Use
`species_by_ANI` for the species, and read a `yes` on a genome whose ANI agrees with its label as
evidence of whole-operon introgression.

## 7. Promoter and CmeR-binding site

| Variable | Type | Values / range | Description |
|---|---|---|---|
| `promoter_length_bp` | integer | 0–663 | Length of the cmeR–cmeA intergenic region recovered. Around 95 bp is typical; 0 means it could not be extracted |
| `CmeR_binding_site_seq` | text | 31 distinct sequences | The sequence aligned to the 16 bp inverted repeat of NCTC 11168 (`TGTAATAAAAATTACA`) |
| `CmeR_binding_site_as_11168` | yes/no | yes 249; no 3,487 | **Deprecated.** Identical to the 11168 repeat. Kept so nothing downstream breaks; use `CmeR_IR_class` |
| `CmeR_IR_halfsite_bp` | integer | 0–8 | Longest n for which the first and last n bases of the window are reverse complements. 0 means the half-sites are broken |
| `CmeR_IR_spacer_bp` | integer | 0–6, blank | Bases between the half-sites. Blank when the half-sites are broken, because the boundary is then arbitrary |
| `CmeR_IR_length_bp` | integer | 12–19 | Total length of the recovered window; 16 in an unaltered operator |
| `CmeR_IR_class` | text | intact 1,750; half_site_changed 1,223; spacer_changed 448; half_site_and_spacer_changed 198 | Structural state of the operator |

CmeR binds as a dimer, one subunit to each half-site of an inverted repeat, either side of a
spacer that holds the two half-sites the right distance apart. The class is therefore read from
the structure of the window — the longest perfect inverted repeat and the gap between its arms —
rather than from identity to one reference strain. NCTC 11168 and *C. jejuni* 81-176 differ by a
single base here and are both intact, and the 81-176 form is the commoner of the two, so anchoring
on 11168 alone would call the majority wildtype form altered.

Descriptive only: campypump reports the sequence and does not test function. The mapping is less
reliable in *C. coli*, whose promoters diverge from the 11168 reference. In *C. jejuni* the repeat
variants track cmeB type, including within a single ST.

## 8. Usability flags

| Variable | Type | Values | Description |
|---|---|---|---|
| `use_for_prevalence` | yes/no | yes 3,732; no 4 | **Filter on this.** `no` = a duplicate record or an excluded genome |
| `duplicate_of` | integer/empty | 3 filled | The `pubmlst_id` of the first record with the same isolate name; this row is a re-upload |
| `exclude_reason` | text/empty | 1 filled | Why a genome was excluded (`annotate --exclude`). Here: 80361, a mixed assembly |

## 9. Metadata (from PubMLST)

| Variable | Type | Description |
|---|---|---|
| `ST` | integer/empty | MLST sequence type as assigned in PubMLST; empty when untypable (e.g. a different species) |
| `clonal_complex` | text/empty | MLST clonal complex; empty when the ST is unassigned |
| `year` | integer | Year of isolation |
| `source` | text | Isolation source (chicken, human stool, pig, cattle faeces, …) |
| `town_or_city` | text/empty | Iquitos, Lima and others; often empty |
| `country` | text | Peru for this collection |

These are recorded values, not recomputed: STs come from PubMLST and are not reassigned by
campypump.

## 10. QC

| Variable | Type | Description |
|---|---|---|
| `warnings` | text/empty | Plain-language notes per gene: contig edges, truncation, missing start codon, divergent allele, provisional partial call. Empty means no problems |

---

## Quick recipes

```stata
* analysis set
keep if use_for_prevalence == "yes"

* RE prevalence by species
tab species_by_ani re_cmeb, row

* by source, C. coli only
tab source re_cmeb if species_by_ani == "C. coli", row

* lineage-aware: the mixed STs are the informative ones
tab st re_cmeb if species_by_ani == "C. coli" & st == "825", row

* the recombinant C. coli RE alleles
list pubmlst_id st if cmeb_type_final == "RE-cmeB" & cmeb_wt_like_segments_pct >= 5

* records whose species label is wrong
list pubmlst_id isolate species_by_ani if species_check_ani == "CHECK"
```

Column names are at most 32 characters with no special symbols, so Stata's `import delimited`
reads them without renaming. Stata lower-cases names on import, as shown above.
