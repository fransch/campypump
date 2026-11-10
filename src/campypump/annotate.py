"""Merge campypump results with PubMLST metadata into one analysis-ready CSV.

Adds:
  duplicate_of        pubmlst_id of the first record with the same isolate name
                      (names compared ignoring leading zeros: CAM000744-3795 == CAM744-003795)
  use_for_prevalence  'no' for duplicate records, else 'yes'
  species_check       'match' / 'MISMATCH' between the --species label and PubMLST's species field
"""
from __future__ import annotations

import csv
import re
import sys

# A single cme gene cannot tell the species apart: the operon moves between C. coli and
# C. jejuni, and C. coli strains (DH161 included) can carry a fully C. jejuni-like cmeR.
# What does separate them is agreement across the three genes that are not the target of
# selection here — cmeR, cmeA and cmeC — scored against the C. jejuni NCTC 11168 versions.
CJ_LIKE_MIN = 92.0
SIGNAL_GENES = ("cmeR", "cmeA", "cmeC")

# What the analysis file needs: the species call, the cmeB call and its evidence, the flags that
# decide whether a row is usable, and the metadata most analyses group by. Everything else —
# per-panel identities, locus coordinates, accessions, the full provenance record — is kept in the
# "full" preset and in the results table the annotation was built from.
CORE = [
    "pubmlst_id", "isolate", "sample",
    "species_by_ANI", "species_check_ANI", "ANI_pct", "ANI_call_quality",
    "cmeB_type_final", "re_cmeB", "cmeB_partial_gene", "cmeB_divergent",
    "operon_species_signal", "operon_atypical_for_label",
    "duplicate_of", "use_for_prevalence", "exclude_reason",
    "ST", "clonal_complex", "year", "source", "town_or_city", "country",
    "cmeB_identity_RE_PA001_pct", "cmeB_identity_WT_11168_pct",
    "cmeB_closest_reference", "cmeB_closest_identity_pct",
    "cmeB_RE_like_segments_pct", "cmeB_WT_like_segments_pct",
    "cmeR_identity_WT_11168_pct", "cmeR_identity_RE_PA001_pct",
    "cmeR_identity_CcWT_RM4661_pct", "cmeR_closest_reference",
    "cmeA_identity_WT_11168_pct", "cmeA_identity_RE_PA001_pct",
    "cmeA_identity_CcWT_RM4661_pct", "cmeA_closest_reference",
    "cmeC_identity_WT_11168_pct", "cmeC_identity_RE_PA001_pct",
    "cmeC_identity_CcWT_RM4661_pct", "cmeC_closest_reference",
    "CmeR_binding_site_seq", "CmeR_binding_site_as_11168",
    "CmeR_IR_halfsite_bp", "CmeR_IR_spacer_bp", "CmeR_IR_length_bp",
    "CmeR_IR_class", "promoter_length_bp",
    "cmeB_length_bp", "warnings",
]

LEAD = ["pubmlst_id", "isolate", "species", "pubmlst_species", "species_check",
        "species_by_ANI", "species_check_ANI", "ANI_pct", "ANI_margin",
        "ANI_align_fraction_pct", "ANI_call_quality",
        "operon_species_signal", "operon_Cj_like_genes", "operon_atypical_for_label",
        "cmeB_type_final", "re_cmeB", "cmeB_partial_gene", "cmeB_divergent",
        "duplicate_of", "use_for_prevalence", "ST", "clonal_complex", "year", "source"]


def norm_isolate(name: str) -> str:
    """Upper-case and drop leading zeros inside number blocks."""
    s = re.sub(r"\s+", "", str(name).upper())
    return re.sub(r"(?<![0-9])0+(?=[0-9])", "", s)


def _read(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        first = fh.readline()
        fh.seek(0)
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        return list(csv.DictReader(fh, delimiter=delim))


def _species_word(s: str) -> str:
    words = re.findall(r"[a-z]+", str(s).lower())
    return words[-1] if words else ""


def species_signal(row: dict) -> tuple[str, str]:
    """Return (n_Cj_like_genes, signal) from cmeR/cmeA/cmeC identity to NCTC 11168.

    Only scored when all three genes are complete: a gene cut at a contig edge has a deflated
    identity and would look like the wrong species.
    """
    warn = row.get("warnings", "") or ""
    values = []
    for gene in SIGNAL_GENES:
        if f"{gene}:" in warn or f"{gene} not found" in warn:
            return "", "incomplete"
        try:
            values.append(float(row[f"{gene}_identity_WT_11168_pct"]))
        except (KeyError, ValueError, TypeError):
            return "", "incomplete"
    votes = sum(v > CJ_LIKE_MIN for v in values)
    if votes == 3:
        return "3", "C. jejuni-like"
    if votes == 0:
        return "0", "C. coli-like"
    return str(votes), "mixed"


def add_species_signal(rows: list[dict]) -> int:
    """Add the operon species-signal columns; return how many rows are flagged for review."""
    n_flagged = 0
    for r in rows:
        votes, signal = species_signal(r)
        r["operon_Cj_like_genes"] = votes
        r["operon_species_signal"] = signal
        label = _species_word(r.get("species", ""))
        if not label or signal in ("incomplete", "mixed"):
            r["operon_atypical_for_label"] = ""
            continue
        expected = "jejuni" if signal == "C. jejuni-like" else "coli"
        ok = label == expected
        r["operon_atypical_for_label"] = "no" if ok else "yes"
        n_flagged += not ok
    return n_flagged


def _apply_exclusions(rows: list[dict], path: str) -> int:
    """Mark genomes listed in a file as not for analysis (one id per line, optional <tab>reason).

    Nothing is deleted: the row keeps its results and gains use_for_prevalence = no plus a reason,
    so an exclusion is visible and reversible.
    """
    reasons = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            reasons[parts[0].strip()] = parts[1].strip() if len(parts) > 1 else "excluded"
    n = 0
    for r in rows:
        for key in (r.get("pubmlst_id", ""), r.get("sample", "")):
            if key and key in reasons:
                r["use_for_prevalence"] = "no"
                r["exclude_reason"] = reasons[key]
                n += 1
                break
    return n


def annotate(results: str, metadata: str | None, out_csv: str,
             ani: str | None = None, exclude: str | None = None,
             columns: str = "core") -> None:
    rows = _read(results)
    meta = {r["pubmlst_id"]: r for r in _read(metadata)} if metadata else {}
    meta_cols = []
    if meta:
        meta_cols = list(next(iter(meta.values())).keys())

    seen: dict[str, str] = {}
    n_dup = n_mis = n_nometa = 0
    for r in rows:
        m = meta.get(r.get("pubmlst_id", ""), {})
        if meta and not m:
            n_nometa += 1
        for k, v in m.items():
            if k not in r or r[k] == "":
                r[k] = v
        key = norm_isolate(r.get("isolate", "")) if r.get("isolate") else ""
        if key and key in seen:
            r["duplicate_of"] = seen[key]
            r["use_for_prevalence"] = "no"
            n_dup += 1
        else:
            if key:
                seen[key] = r.get("pubmlst_id", "")
            r["duplicate_of"] = ""
            r["use_for_prevalence"] = "yes"
        if r.get("species") and r.get("pubmlst_species"):
            ok = _species_word(r["species"]) == _species_word(r["pubmlst_species"])
            r["species_check"] = "match" if ok else "MISMATCH"
            n_mis += not ok
        else:
            r["species_check"] = ""

    n_flagged = add_species_signal(rows)
    n_ani_missing = n_ani_flagged = 0
    if ani:
        from .species import add_ani_columns, print_summary
        n_ani_missing, n_ani_flagged = add_ani_columns(rows, ani)
        print_summary(rows)
    n_excluded = _apply_exclusions(rows, exclude) if exclude else 0
    all_cols = list(dict.fromkeys(list(rows[0].keys()) + meta_cols +
                                  ["duplicate_of", "use_for_prevalence", "exclude_reason",
                                   "species_check", "operon_species_signal",
                                   "operon_Cj_like_genes", "operon_atypical_for_label"]))
    if columns == "core":
        cols = [c for c in CORE if c in all_cols]
        dropped = len(all_cols) - len(cols)
    else:
        cols = [c for c in LEAD if c in all_cols] + [c for c in all_cols if c not in LEAD]
        dropped = 0
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(rows)
    msg = f"{len(rows)} rows ({len(cols)} columns) -> {out_csv}  duplicates flagged: {n_dup}"
    if n_flagged:
        msg += f"  operon atypical for label: {n_flagged}"
    if ani:
        msg += f"  ANI/species mismatches: {n_ani_flagged}"
        if n_ani_missing:
            msg += f"  without ANI: {n_ani_missing}"
    if exclude:
        msg += f"  excluded: {n_excluded}"
    if columns == "core" and dropped:
        msg += f"\n{dropped} further columns are available with --columns full"
    if meta:
        msg += f"  species mismatches: {n_mis}  without metadata: {n_nometa}"
    print(msg, file=sys.stderr)
