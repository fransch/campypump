"""Prepare cmeB alleles for phylogenetics: deduplicate, name, and annotate.

campypump does not infer trees — alignment and tree inference belong to established tools
(MAFFT, IQ-TREE). What it does is produce their input reliably: one record per distinct allele,
with tree-safe names, plus an annotation table keyed by those names for iTOL/ggtree/Bio.Phylo.

Outputs
  <out>.fasta   one record per distinct allele; id like  Cc_a001_RE_n74
  <out>.tsv     one row per allele: species, type, genomes carrying it, STs, block status, ...
  <out>_itol_*.txt  optional iTOL colour-strip files (type, species)
"""
from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter, defaultdict

# Newick chokes on these; keep ids to letters, digits and underscore
SAFE = re.compile(r"[^A-Za-z0-9_]+")

SPECIES_ABBREV = {"c. coli": "Cc", "c. jejuni": "Cj", "campylobacter coli": "Cc",
                  "campylobacter jejuni": "Cj"}
TYPE_ABBREV = {"RE-cmeB": "RE", "wildtype-cmeB": "WT", "mosaic": "MOS",
               "ambiguous": "AMB", "unresolved": "UNR", "": "NA"}

from .palettes import SPECIES as _SP, TYPE_SHORT as _TY

ITOL_TYPE_COLORS = dict(_TY)
ITOL_SPECIES_COLORS = {"Cc": _SP["Cc"], "Cj": _SP["Cj"], "ref": _SP["REF"], "out": _SP["OUT"]}


def abbrev_species(name: str) -> str:
    key = str(name).strip().lower()
    if key in SPECIES_ABBREV:
        return SPECIES_ABBREV[key]
    return SAFE.sub("", str(name))[:4] or "Sp"


def read_table(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        first = fh.readline()
        fh.seek(0)
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        return list(csv.DictReader(fh, delimiter=delim))


def collect_alleles(fasta_paths, labels, results_rows, min_frac=0.95, include_partial=False):
    """Group genomes by identical cmeB sequence.

    Returns (alleles, n_partial) where each allele is a dict with the sequence, species,
    cmeB type, the genomes carrying it, and whatever metadata the results table provides.
    """
    from Bio import SeqIO

    by_sample = {r.get("sample") or r.get("pubmlst_id", ""): r for r in results_rows}
    groups: dict[tuple, dict] = {}
    n_partial = 0
    for path, label in zip(fasta_paths, labels):
        recs = list(SeqIO.parse(path, "fasta"))
        if not recs:
            sys.exit(f"No sequences in {path}")
        full = max(len(r.seq) for r in recs)
        for r in recs:
            if len(r.seq) < min_frac * full and not include_partial:
                n_partial += 1
                continue
            seq = str(r.seq).upper()
            meta = by_sample.get(r.id, {})
            cls = (meta.get("cmeB_type_final") or meta.get("cmeB_type")
                   or (r.description.split()[-1] if " " in r.description else ""))
            key = (label, seq)
            g = groups.setdefault(key, {"species": label, "seq": seq, "samples": [],
                                        "types": Counter(), "sts": Counter(),
                                        "sources": Counter(), "wt_segments": []})
            g["samples"].append(r.id)
            g["types"][cls] += 1
            if meta.get("ST"):
                g["sts"][meta["ST"]] += 1
            if meta.get("source"):
                g["sources"][meta["source"]] += 1
            try:
                g["wt_segments"].append(float(meta["cmeB_WT_like_segments_pct"]))
            except (KeyError, ValueError, TypeError):
                pass
    alleles = sorted(groups.values(), key=lambda g: (g["species"], -len(g["samples"])))
    return alleles, n_partial


def name_alleles(alleles: list[dict]) -> None:
    """Assign tree-safe ids: <species>_a<rank>_<type>_n<genomes>, e.g. Cc_a001_RE_n74."""
    rank = defaultdict(int)
    for a in alleles:
        sp = abbrev_species(a["species"])
        rank[sp] += 1
        cls = a["types"].most_common(1)[0][0] if a["types"] else ""
        a["type"] = cls
        a["type_short"] = TYPE_ABBREV.get(cls, SAFE.sub("", cls)[:3] or "NA")
        a["species_short"] = sp
        a["rank"] = rank[sp]
        a["id"] = f"{sp}_a{rank[sp]:03d}_{a['type_short']}_n{len(a['samples'])}"


def add_extra(alleles: list[dict], path: str, kind: str) -> None:
    """Add reference or outgroup sequences as their own 'alleles' (kind: REF or OUT)."""
    from Bio import SeqIO

    for r in SeqIO.parse(path, "fasta"):
        sid = SAFE.sub("_", r.id)[:40]
        # don't double up the tag if the record is already named for it
        if sid.upper().startswith(kind.upper() + "_"):
            sid = sid[len(kind) + 1:]
        alleles.append({"species": kind.lower(), "species_short": kind.lower(), "seq": str(r.seq).upper(),
                        "samples": [r.id], "types": Counter(), "sts": Counter(), "sources": Counter(),
                        "wt_segments": [], "type": kind, "type_short": kind, "rank": 0,
                        "id": f"{kind}_{sid}"})


def write_outputs(alleles: list[dict], out_base: str, itol: bool = True) -> list[str]:
    made = []
    fa = out_base + ".fasta"
    with open(fa, "w") as fh:
        for a in alleles:
            fh.write(f">{a['id']}\n{a['seq']}\n")
    made.append(fa)

    tsv = out_base + ".tsv"
    cols = ["allele_id", "species", "cmeB_type", "genomes", "length_bp", "allele_rank",
            "top_STs", "n_STs", "top_sources", "mean_WT_like_segments_pct",
            "block_carrier", "example_sample", "all_samples"]
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for a in alleles:
            wt = a["wt_segments"]
            mean_wt = sum(wt) / len(wt) if wt else ""
            w.writerow({
                "allele_id": a["id"],
                "species": a["species"],
                "cmeB_type": a["type"],
                "genomes": len(a["samples"]),
                "length_bp": len(a["seq"]),
                "allele_rank": a["rank"] or "",
                "top_STs": "; ".join(f"{s}({n})" for s, n in a["sts"].most_common(3)),
                "n_STs": len(a["sts"]),
                "top_sources": "; ".join(f"{s}({n})" for s, n in a["sources"].most_common(3)),
                "mean_WT_like_segments_pct": f"{mean_wt:.1f}" if mean_wt != "" else "",
                # the ~300 bp wildtype-like block seen in some RE alleles; only meaningful for
                # RE alleles, since wildtype alleles are wildtype-like along their whole length
                "block_carrier": "" if (mean_wt == "" or a["type"] != "RE-cmeB")
                                 else ("yes" if mean_wt >= 5 else "no"),
                "example_sample": a["samples"][0],
                "all_samples": "; ".join(a["samples"][:50]) + ("; ..." if len(a["samples"]) > 50 else ""),
            })
    made.append(tsv)

    if itol:
        made.append(_itol_strip(alleles, out_base + "_itol_type.txt", "cmeB type",
                                lambda a: a["type_short"], ITOL_TYPE_COLORS))
        made.append(_itol_strip(alleles, out_base + "_itol_species.txt", "species",
                                lambda a: a["species_short"], ITOL_SPECIES_COLORS))
        made.append(_itol_sizes(alleles, out_base + "_itol_genomes.txt"))
    made.append(write_microreact(alleles, out_base + "_microreact.csv"))
    return made


def _itol_strip(alleles, path, label, key, colours) -> str:
    seen = []
    lines = ["DATASET_COLORSTRIP", "SEPARATOR TAB", f"DATASET_LABEL\t{label}",
             "COLOR\t#333333", "STRIP_WIDTH\t28", "MARGIN\t4", "BORDER_WIDTH\t0.5",
             "SHOW_STRIP_LABELS\t0"]
    body = []
    for a in alleles:
        k = key(a)
        colour = colours.get(k, "#999999")
        if k not in seen:
            seen.append(k)
        body.append(f"{a['id']}\t{colour}\t{k}")
    lines.append("LEGEND_TITLE\t" + label)
    lines.append("LEGEND_SHAPES\t" + "\t".join("1" for _ in seen))
    lines.append("LEGEND_COLORS\t" + "\t".join(colours.get(k, "#999999") for k in seen))
    lines.append("LEGEND_LABELS\t" + "\t".join(seen))
    lines += ["DATA"] + body
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def write_microreact(alleles: list[dict], path: str) -> str:
    """Metadata CSV for Microreact (https://microreact.org): an ID column matching the tree tips,
    plus <column>__colour columns so colours are fixed rather than assigned at random.

    Note that the tips here are alleles, not isolates, so Microreact's map and timeline panels have
    nothing to show; the value is the tree plus a sortable table. For a map and timeline, build a
    tree whose tips are isolates.
    """
    cols = ["ID", "species", "species__colour", "cmeB_type", "cmeB_type__colour",
            "genomes", "allele_rank", "top_STs", "top_sources",
            "block_carrier", "block_carrier__colour", "length_bp", "example_sample"]
    block_colour = {"yes": _TY["MOS"], "no": "#D9D9D9", "": "#FFFFFF"}
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for a in alleles:
            wt = a["wt_segments"]
            mean_wt = sum(wt) / len(wt) if wt else ""
            block = "" if (mean_wt == "" or a["type"] != "RE-cmeB") else ("yes" if mean_wt >= 5 else "no")
            w.writerow({
                "ID": a["id"],
                "species": a["species"],
                "species__colour": ITOL_SPECIES_COLORS.get(a["species_short"], "#999999"),
                "cmeB_type": a["type"] or "reference/outgroup",
                "cmeB_type__colour": ITOL_TYPE_COLORS.get(a["type_short"], "#999999"),
                "genomes": len(a["samples"]),
                "allele_rank": a["rank"] or "",
                "top_STs": "; ".join(f"{s}({n})" for s, n in a["sts"].most_common(3)),
                "top_sources": "; ".join(f"{s}({n})" for s, n in a["sources"].most_common(3)),
                "block_carrier": block,
                "block_carrier__colour": block_colour.get(block, "#FFFFFF"),
                "length_bp": len(a["seq"]),
                "example_sample": a["samples"][0],
            })
    return path


def _itol_sizes(alleles, path) -> str:
    lines = ["DATASET_SIMPLEBAR", "SEPARATOR TAB", "DATASET_LABEL\tgenomes carrying the allele",
             "COLOR\t#666666", "WIDTH\t120", "DATASET_SCALE\t1\t10\t100", "LOG_SCALE\t1", "DATA"]
    lines += [f"{a['id']}\t{len(a['samples'])}" for a in alleles]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path
