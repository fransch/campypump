"""Species assignment from whole-genome ANI (skani output).

No single cme gene identifies the species: the operon crosses between *C. coli* and *C. jejuni*,
and *C. coli* DH161 carries a 99.4 % *C. jejuni*-like cmeR. Whole-genome ANI is unaffected by a
single transferred locus, so it is what campypump uses for the species call; the operon columns
describe the operon only.

Input is the table written by:

    skani dist -q genomes/*/*.fasta -r species_refs/*.fasta -t 8 -o species_ani.tsv

Reference file names supply the species, so name them for the organism
(`Ccoli_RM4661.fasta`, `C.jejuni_NCTC11168.fasta`, `Cupsaliensis_RM3940.fasta`, ...).
"""
from __future__ import annotations

import csv
import os
import re
import sys

ANI_SPECIES_MIN = 95.0        # conventional species boundary
ANI_MARGIN_MIN = 1.0          # best hit must beat the runner-up by this much
ANI_ALIGN_MIN = 60.0          # % of the query aligned; below this the call rests on little data

_SPECIES_PATTERNS = [
    (re.compile(r"jejuni", re.I), "C. jejuni"),
    (re.compile(r"coli", re.I), "C. coli"),
    (re.compile(r"lari", re.I), "C. lari"),
    (re.compile(r"upsaliensis", re.I), "C. upsaliensis"),
    (re.compile(r"helveticus", re.I), "C. helveticus"),
    (re.compile(r"fetus", re.I), "C. fetus"),
    (re.compile(r"hyointestinalis", re.I), "C. hyointestinalis"),
    (re.compile(r"concisus", re.I), "C. concisus"),
]


def species_from_reference(ref_file: str, ref_name: str = "") -> str:
    """Work out the species from the reference file name, falling back to its FASTA header."""
    for text in (os.path.basename(ref_file), ref_name):
        for pattern, species in _SPECIES_PATTERNS:
            if pattern.search(text or ""):
                return species
    return os.path.splitext(os.path.basename(ref_file))[0]


def sample_from_query(query_file: str) -> str:
    name = os.path.basename(query_file)
    for ext in (".fasta", ".fas", ".fa", ".fna"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def read_ani(path: str) -> dict[str, dict]:
    """Best hit per genome, with the margin to the runner-up."""
    hits: dict[str, list[tuple[float, str, float]]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                ani = float(row["ANI"])
                align = float(row["Align_fraction_query"])
            except (KeyError, ValueError):
                continue
            sample = sample_from_query(row["Query_file"])
            species = species_from_reference(row["Ref_file"], row.get("Ref_name", ""))
            hits.setdefault(sample, []).append((ani, species, align))

    out = {}
    for sample, rows in hits.items():
        rows.sort(reverse=True)
        ani, species, align = rows[0]
        # margin against the best hit of a *different* species
        other = next((a for a, s, _ in rows if s != species), None)
        out[sample] = {
            "species_by_ANI": species if ani >= ANI_SPECIES_MIN else "",
            "ANI_pct": f"{ani:.2f}",
            "ANI_margin": f"{ani - other:.2f}" if other is not None else "",
            "ANI_align_fraction_pct": f"{align:.1f}",
            "ANI_call_quality": _quality(ani, other, align),
        }
    return out


def _quality(ani: float, other: float | None, align: float) -> str:
    if ani < ANI_SPECIES_MIN:
        return "no reference above 95% ANI"
    if align < ANI_ALIGN_MIN:
        return f"low alignment fraction ({align:.0f}%)"
    if other is not None and (ani - other) < ANI_MARGIN_MIN:
        return "close to another species"
    return "ok"


def add_ani_columns(rows: list[dict], ani_path: str) -> tuple[int, int]:
    """Merge ANI calls into rows; return (n_without_ANI, n_disagreeing_with_label)."""
    from .annotate import _species_word

    ani = read_ani(ani_path)
    missing = flagged = 0
    for r in rows:
        key = r.get("sample") or ""
        rec = ani.get(key) or ani.get(f"id-{r.get('pubmlst_id', '')}")
        if not rec:
            missing += 1
            for c in ("species_by_ANI", "ANI_pct", "ANI_margin",
                      "ANI_align_fraction_pct", "ANI_call_quality", "species_check_ANI"):
                r.setdefault(c, "")
            continue
        r.update(rec)
        label, called = _species_word(r.get("species", "")), _species_word(rec["species_by_ANI"])
        if not label or not called:
            r["species_check_ANI"] = ""
        elif label == called:
            r["species_check_ANI"] = "consistent"
        else:
            r["species_check_ANI"] = "CHECK"
            flagged += 1
    return missing, flagged


def print_summary(rows: list[dict], out=sys.stderr) -> None:
    flagged = [r for r in rows if r.get("species_check_ANI") == "CHECK"]
    if not flagged:
        return
    print(f"{len(flagged)} genomes: whole-genome ANI disagrees with the species label", file=out)
    for r in sorted(flagged, key=lambda x: x.get("pubmlst_id", "")):
        print(f"    {r.get('sample', '')}  labelled {r.get('species', '?')}  "
              f"-> {r.get('species_by_ANI', '?')} (ANI {r.get('ANI_pct', '?')}%, "
              f"aligned {r.get('ANI_align_fraction_pct', '?')}%)", file=out)
