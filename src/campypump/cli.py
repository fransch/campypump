"""Command-line interface.

  campypump run   -o results.tsv genomes/*.fasta
  campypump fetch --ids ids.txt --outdir genomes/
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
import time
import urllib.error
import urllib.request
from multiprocessing import Pool

from . import __version__
from .alleles import add_extra, collect_alleles, name_alleles, read_table, write_outputs
from .annotate import add_species_signal, annotate
from .figures import classification_figure, diversity_figure
from .treeplot import tree_figure
from .plots import write_segment_plot
from .pubmlst import collect_metadata
from .core import FASTA_EXTS, Params, check_blast, classify_genome, load_refs, pubmlst_id_from_name

# column order and plain-language descriptions (also written to <output>_column_guide.txt)
COLUMN_GUIDE = [
    ("sample", "Genome file name without extension"),
    ("pubmlst_id", "PubMLST isolate id, taken from file names like id-152096.fasta"),
    ("species", "Species label given with --species (empty if not given)"),
    ("cmeB_type", "Call from the full gene: RE-cmeB | wildtype-cmeB | ambiguous | mosaic | incomplete | not_found"),
    ("cmeB_type_final", "Type to use in analyses: as cmeB_type, but partial (incomplete) genes get a provisional RE/wildtype call when >=90% of the part present agrees"),
    ("re_cmeB", "1 = RE-cmeB, 0 = wildtype-cmeB, empty = mosaic/ambiguous/unresolved (from cmeB_type_final)"),
    ("cmeB_partial_gene", "yes = cmeB broken across contigs or truncated; its type is provisional"),
    ("cmeB_divergent", "yes = complete gene <90% identical to every reference (unusual allele; check species)"),
    ("cmeB_identity_RE_PA001_pct", "% identity of cmeB to RE-cmeB of PA001 (same as Geneious column '% identity-cjREcmeB')"),
    ("cmeB_identity_WT_11168_pct", "% identity of cmeB to wildtype cmeB of NCTC 11168 (same as Geneious '% identity-cjcmeB11168')"),
    ("cmeB_identity_DH161_pct", "% identity of cmeB to the published RE-cmeB of C. coli DH161 (GenBank KT778508.1)"),
    ("RE_by_DH161_definition", "yes = >90% identity to DH161 cmeB, the published RE-cmeB definition (Cooper et al. 2025)"),
    ("cmeB_closest_reference", "Reference allele in the panel most similar to this cmeB"),
    ("cmeB_closest_identity_pct", "% identity to that closest reference"),
    ("cmeB_RE_like_segments_pct", "% of 300-bp segments along cmeB clearly closer to RE-cmeB"),
    ("cmeB_WT_like_segments_pct", "% of 300-bp segments along cmeB clearly closer to a wildtype cmeB (mosaic if both >= 10%)"),
    ("cmeB_aligned_identity_RE_pct", "% identity to PA001 over the part of the gene present (differs from the full-gene value only for partial genes; Geneious-style)"),
    ("cmeB_aligned_identity_WT_pct", "% identity to 11168 over the part of the gene present"),
    ("species_by_ANI", "Species from whole-genome ANI against the reference genomes (skani); the species call to use"),
    ("species_check_ANI", "CHECK = ANI disagrees with the --species label; consistent = agrees"),
    ("ANI_pct", "ANI to the best-matching reference genome"),
    ("ANI_margin", "How far the best hit beats the best hit of a different species (percentage points)"),
    ("ANI_align_fraction_pct", "% of the genome that aligned to that reference; a low value means the call rests on little data"),
    ("ANI_call_quality", "ok, or why the ANI call is weak (low alignment fraction, close to another species, nothing above 95%)"),
    ("exclude_reason", "Why a genome was excluded from analysis (from annotate --exclude); use_for_prevalence is set to no"),
    ("operon_species_signal", "Whether cmeR, cmeA and cmeC look C. jejuni-like, C. coli-like, mixed, or could not be scored (incomplete). No single cme gene identifies the species; this is the agreement of the three"),
    ("operon_Cj_like_genes", "How many of cmeR, cmeA, cmeC are >92% identical to NCTC 11168 (0-3)"),
    ("operon_atypical_for_label", "yes = the cme operon does not match the labelled species. This is NOT a species call: it means either a database error or complete introgression of the operon. Use species_by_ANI for the species"),
    ("cmeR_identity_WT_11168_pct", "% identity of cmeR to NCTC 11168"),
    ("cmeR_identity_RE_PA001_pct", "% identity of cmeR to PA001"),
    ("cmeR_identity_CcWT_RM4661_pct", "% identity of cmeR to C. coli RM4661 (CP007181.1), the native C. coli wildtype reference"),
    ("cmeR_closest_reference", "Which reference allele of cmeR this gene matches best, with its class"),
    ("cmeR_closest_identity_pct", "% identity of cmeR to that closest reference"),
    ("cmeA_identity_WT_11168_pct", "% identity of cmeA to NCTC 11168"),
    ("cmeA_identity_RE_PA001_pct", "% identity of cmeA to PA001"),
    ("cmeA_identity_CcWT_RM4661_pct", "% identity of cmeA to C. coli RM4661 (CP007181.1), the native C. coli wildtype reference"),
    ("cmeA_closest_reference", "Which reference allele of cmeA this gene matches best, with its class"),
    ("cmeA_closest_identity_pct", "% identity of cmeA to that closest reference"),
    ("cmeC_identity_WT_11168_pct", "% identity of cmeC to NCTC 11168"),
    ("cmeC_identity_RE_PA001_pct", "% identity of cmeC to PA001"),
    ("cmeC_identity_CcWT_RM4661_pct", "% identity of cmeC to C. coli RM4661 (CP007181.1), the native C. coli wildtype reference"),
    ("cmeC_closest_reference", "Which reference allele of cmeC this gene matches best, with its class"),
    ("cmeC_closest_identity_pct", "% identity of cmeC to that closest reference"),
    ("promoter_length_bp", "Length of the cmeR-cmeA intergenic (promoter) region"),
    ("CmeR_binding_site_seq", "Sequence aligned to the 16-bp CmeR-binding inverted repeat of 11168 (TGTAATAAAAATTACA)"),
    ("CmeR_binding_site_as_11168", "DEPRECATED: yes = identical to the 11168 inverted repeat. NCTC 11168 carries a minority form of the operator, so this is a poor baseline; use CmeR_IR_class"),
    ("CmeR_IR_halfsite_bp", "Length of the inverted-repeat half-sites: the longest n for which the first and last n bases are reverse complements"),
    ("CmeR_IR_spacer_bp", "Bases between the two half-sites; blank when the half-sites are broken and the boundary is arbitrary"),
    ("CmeR_IR_length_bp", "Total length of the operator window"),
    ("CmeR_IR_class", "intact / half_site_changed / spacer_changed / half_site_and_spacer_changed"),
    ("cmeB_length_bp", "Length of the extracted cmeB gene (bp)"),
    ("cmeB_contig", "Contig carrying cmeB"),
    ("cmeB_start", "cmeB start position on the contig"),
    ("cmeB_end", "cmeB end position on the contig"),
    ("cmeB_strand", "Strand of cmeB on the contig (+/-)"),
    ("panel_<class>_<strain>_pct", "% identity of cmeB to each extra reference in the panel (CjWT, CjRE, CcWT, CcRE)"),
    ("warnings", "Quality or interpretation notes; empty = no problems"),
]
FIRST_COLS = [c for c, _ in COLUMN_GUIDE if not c.startswith("panel_") and c != "warnings"]


def cmd_combine(a):
    """Stack several campypump result tables into one (header written once)."""
    header, out_rows = None, []
    for path in a.tables:
        with open(path) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if header is None:
                header = rd.fieldnames
            else:
                header += [c for c in rd.fieldnames if c not in header]
            out_rows += list(rd)
    with open(a.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(out_rows)
    print(f"{len(out_rows)} rows from {len(a.tables)} tables -> {a.output}", file=sys.stderr)

PUBMLST = "https://rest.pubmlst.org/db/pubmlst_campylobacter_isolates/isolates/{}/contigs_fasta"


def _worker(args):
    genome, refs_path, params = args
    try:
        return classify_genome(genome, load_refs(refs_path), params)
    except Exception as e:  # keep batch running; report the failure in the table
        return {"sample": os.path.basename(genome), "cmeB_type": "error", "warnings": repr(e)}


def cmd_run(a):
    check_blast()
    genomes = _collect_genomes(a.genomes)
    p = Params(re_min=a.re_min, wt_max_re=a.wt_max_re, window=a.window, step=a.step,
               mosaic_min_frac=a.mosaic_frac)
    jobs = [(g, a.refs, p) for g in genomes]
    if a.threads > 1:
        with Pool(a.threads) as pool:
            rows = pool.map(_worker, jobs, chunksize=4)
    else:
        rows = [_worker(j) for j in jobs]

    if a.species:
        for r in rows:
            r["species"] = a.species
    n_flagged = add_species_signal(rows)
    n_ani_flagged = _add_species_call(a, genomes, rows)
    extra = sorted({k for r in rows for k in r if not k.startswith("_")} - set(FIRST_COLS) - {"warnings"})
    cols = FIRST_COLS + extra + ["warnings"]
    with open(a.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    if a.alleles:
        with open(a.alleles, "w") as fh:
            for r in rows:
                if r.get("_cmeB_seq"):
                    fh.write(f">{r['sample']} {r.get('cmeB_type_final', r['cmeB_type'])}\n{r['_cmeB_seq']}\n")
    guide = os.path.splitext(a.output)[0] + "_column_guide.txt"
    with open(guide, "w") as fh:
        fh.write("campypump output columns\n\n")
        for c, d in COLUMN_GUIDE:
            fh.write(f"{c}\n    {d}\n")
    if a.segment_plots:
        n_plots = 0
        for r in rows:
            if not r.get("_windows"):
                continue
            flagged = (r.get("cmeB_type_final") in ("mosaic", "ambiguous", "unresolved")
                       or r.get("cmeB_divergent") == "yes"
                       or (float(r.get("cmeB_RE_like_segments_pct") or 0) >= 5
                           and float(r.get("cmeB_WT_like_segments_pct") or 0) >= 5))
            if flagged or a.plot_all:
                write_segment_plot(a.segment_plots, r)
                n_plots += 1
        print(f"{n_plots} segment plots -> {a.segment_plots}/", file=sys.stderr)
    if n_flagged:
        print(f"{n_flagged} genomes: the cmeR/cmeA/cmeC signal disagrees with --species "
              f"(species_check_operon = CHECK) — confirm with rMLST", file=sys.stderr)
    if n_ani_flagged:
        print(f"{n_ani_flagged} genomes: the whole-genome species call disagrees with the "
              f"label (species_check_ANI = CHECK) — use species_by_ANI", file=sys.stderr)
    counts = {}
    for r in rows:
        counts[r.get("cmeB_type_final", r["cmeB_type"])] = counts.get(r.get("cmeB_type_final", r["cmeB_type"]), 0) + 1
    print(f"{len(rows)} genomes -> {a.output}  " + "  ".join(f"{k}: {v}" for k, v in sorted(counts.items())),
          file=sys.stderr)




def _collect_genomes(paths):
    """FASTA files from the paths given; a directory is searched at every depth.

    Collections are routinely laid out one directory per species or per isolate
    (`genomes/ccoli/id-48.fasta`). Searching only the top level meant pointing at the
    collection root silently found nothing, so the search walks the tree.
    """
    genomes, searched = [], []
    for g in paths:
        if os.path.isdir(g):
            found = sorted(os.path.join(root, f)
                           for root, _dirs, files in os.walk(g)
                           for f in files if f.endswith(FASTA_EXTS))
            searched.append((g, len(found)))
            genomes += found
        elif os.path.exists(g):
            genomes.append(g)
        else:
            sys.exit(f"No such file or directory: {g}")
    for d, n in searched:
        print(f"{n} FASTA files under {d}", file=sys.stderr)
    if not genomes:
        where = ", ".join(d for d, _ in searched) or ", ".join(paths)
        sys.exit(f"No genome FASTA files found in {where} "
                 f"(looked at every depth for {', '.join(FASTA_EXTS)})")
    return genomes


def _add_species_call(a, genomes, rows):
    """Whole-genome species call for `run`. Returns the number of label disagreements.

    The operon crosses between species, so it cannot define the species; ANI against four
    reference genomes can. skani is run here rather than left to the user, because a user
    who skips it is the user most likely to be working from a mislabelled record.
    """
    if getattr(a, "no_ani", False):
        return 0
    from .species import add_ani_columns, run_skani, skani_available, SKANI_MISSING_HINT
    table = getattr(a, "ani", None)
    if not table:
        if not skani_available():
            print(SKANI_MISSING_HINT, file=sys.stderr)
            return 0
        table = os.path.splitext(a.output)[0] + "_species_ani.tsv"
        try:
            run_skani(genomes, table, threads=a.threads)
        except (OSError, RuntimeError) as exc:
            # a species call is worth having, but not at the cost of losing the typing
            print(f"skani did not run, so no species call was made: {exc}", file=sys.stderr)
            return 0
    try:
        missing, flagged = add_ani_columns(rows, table)
    except (OSError, KeyError) as exc:
        print(f"could not read the ANI table {table}: {exc}", file=sys.stderr)
        return 0
    if missing:
        print(f"{missing} genomes had no ANI result; their species columns are empty",
              file=sys.stderr)
    print(f"species call from whole-genome ANI -> {table}", file=sys.stderr)
    return flagged


def _read_ids(path, column):
    if path.endswith((".xlsx", ".xls")):
        try:
            import pandas as pd
        except ImportError:
            sys.exit("Reading Excel files needs pandas + openpyxl; or give a plain text file of IDs.")
        return [str(int(x)) for x in pd.read_excel(path)[column].dropna()]
    with open(path) as fh:
        return [l.strip().split()[0] for l in fh if l.strip() and not l.startswith("#")]


def _figure_rows(path, exclude_dups):
    import csv as _csv
    with open(path, newline="") as fh:
        first = fh.readline()
        fh.seek(0)
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        rows = list(_csv.DictReader(fh, delimiter=delim))
    out = []
    for r in rows:
        if exclude_dups and r.get("use_for_prevalence") == "no":
            continue
        try:
            re_, wt = float(r["cmeB_identity_RE_PA001_pct"]), float(r["cmeB_identity_WT_11168_pct"])
        except (KeyError, ValueError):
            continue
        out.append({"re": re_, "wt": wt,
                    "cls": r.get("cmeB_type_final") or r.get("cmeB_type", ""),
                    # the annotated file carries species_by_ANI; the raw results table carries species
                    "species": r.get("species_by_ANI") or r.get("species", ""),
                    "partial": r.get("cmeB_partial_gene") == "yes"})
    return out


def cmd_figure(a):
    rows = _figure_rows(a.results, a.exclude_duplicates)
    # partial genes have deflated whole-gene identities, so they plot far from their class
    n_partial = sum(r["partial"] for r in rows)
    if not a.include_partial:
        rows = [r for r in rows if not r["partial"]]
    if not rows:
        sys.exit("No usable rows found (need cmeB_identity_RE_PA001_pct and cmeB_identity_WT_11168_pct).")
    os.makedirs(a.outdir, exist_ok=True)
    made = []
    species = sorted({r["species"] for r in rows if r["species"]})
    unit = "genomes" if a.include_partial else "genomes, complete cmeB only"
    for sp in species:
        sub_rows = [r for r in rows if r["species"] == sp]
        slug = sp.lower().replace(".", "").replace(" ", "_")
        made.append(classification_figure(sub_rows, os.path.join(a.outdir, f"cmeB_classification_{slug}.svg"),
                                          subtitle=f"{sp}  (n = {len(sub_rows):,} {unit})"))
    if a.combined or not species:
        made.append(classification_figure(rows, os.path.join(a.outdir, "cmeB_classification_all.svg"),
                                          subtitle=f"all species (n = {len(rows):,} {unit})"))
    for f in made:
        print(f, file=sys.stderr)
    if n_partial and not a.include_partial:
        print(f"({n_partial} partial genes omitted; --include-partial to show them)", file=sys.stderr)


def _label_from_name(path):
    n = os.path.basename(path).lower()
    if "ccoli" in n or "c_coli" in n:
        return "C. coli"
    if "cjejuni" in n or "c_jejuni" in n:
        return "C. jejuni"
    return os.path.splitext(os.path.basename(path))[0]


def cmd_diversity(a):
    import random
    from collections import Counter
    from Bio import SeqIO

    labels = a.labels or [_label_from_name(f) for f in a.alleles]
    if len(labels) != len(a.alleles):
        sys.exit("--labels must give one label per allele file")
    per_species, table_rows, n_partial_all = [], [], 0
    for path, label in zip(a.alleles, labels):
        recs = list(SeqIO.parse(path, "fasta"))
        if not recs:
            sys.exit(f"No sequences in {path}")
        full = max(len(r.seq) for r in recs)
        keep, n_partial = [], 0
        for r in recs:
            # a gene broken across contigs is a fragment, not a distinct allele
            if len(r.seq) < 0.95 * full and not a.include_partial:
                n_partial += 1
                continue
            keep.append(r)
        n_partial_all += n_partial
        seqs = [str(r.seq) for r in keep]
        cls_of = {r.id: (r.description.split()[-1] if " " in r.description else "") for r in keep}
        counts = Counter(seqs)
        ranked = sorted(counts.values(), reverse=True)

        # rarefaction: mean distinct alleles seen in a random subsample of n genomes
        random.seed(0)
        acc, n = [], len(seqs)
        steps = sorted({int(round(x)) for x in
                        [n * i / 40 for i in range(1, 41)]} | {1})
        for k in steps:
            reps = 30 if k < n else 1
            acc.append((k, sum(len(set(random.sample(seqs, k))) for _ in range(reps)) / reps))

        by_cls = {}
        for r in keep:
            c = cls_of[r.id]
            by_cls.setdefault(c, {"genomes": 0, "alleles": set()})
            by_cls[c]["genomes"] += 1
            by_cls[c]["alleles"].add(str(r.seq))
        byclass = [(c, v["genomes"], len(v["alleles"])) for c, v in by_cls.items()
                   if v["genomes"] >= 0.01 * len(keep)]
        byclass.sort(key=lambda t: -t[1])
        per_species.append((label, ranked, acc, byclass))

        print(f"{label}: {len(keep):,} genomes, {len(counts)} distinct alleles, "
              f"{sum(1 for v in counts.values() if v == 1)} singletons"
              + (f" ({n_partial} partial genes omitted)" if n_partial else ""), file=sys.stderr)
        for c, v in by_cls.items():
            print(f"    {c or 'unlabelled'}: {v['genomes']:,} genomes, {len(v['alleles'])} alleles",
                  file=sys.stderr)

        if a.table:
            order = {s: i for i, (s, _) in enumerate(counts.most_common(), 1)}
            for seq, n_gen in counts.most_common():
                ids = [r.id for r in keep if str(r.seq) == seq]
                table_rows.append({"species": label, "allele_rank": order[seq], "genomes": n_gen,
                                   "length_bp": len(seq), "cmeB_type": cls_of[ids[0]],
                                   "example_sample": ids[0]})
    if a.table:
        with open(a.table, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["species", "allele_rank", "genomes", "length_bp",
                                               "cmeB_type", "example_sample"], delimiter="\t")
            w.writeheader()
            w.writerows(table_rows)
        print(f"allele table -> {a.table}", file=sys.stderr)
    note = "complete genes only" if n_partial_all and not a.include_partial else ""
    print(diversity_figure(per_species, a.output, note=note), file=sys.stderr)


def cmd_alleles(a):
    labels = a.labels or [_label_from_name(f) for f in a.fasta]
    if len(labels) != len(a.fasta):
        sys.exit("--labels must give one label per FASTA")
    rows = read_table(a.results) if a.results else []
    alleles, n_partial = collect_alleles(a.fasta, labels, rows, include_partial=a.include_partial)
    name_alleles(alleles)
    if a.references:
        add_extra(alleles, a.references, "REF")
    if a.outgroup:
        add_extra(alleles, a.outgroup, "OUT")
    outdir = os.path.dirname(os.path.abspath(a.output))
    os.makedirs(outdir, exist_ok=True)
    made = write_outputs(alleles, a.output, itol=not a.no_itol)

    from collections import Counter
    per_species = Counter(x["species"] for x in alleles)
    for sp, n in per_species.items():
        print(f"{sp}: {n} distinct alleles", file=sys.stderr)
    if n_partial:
        print(f"({n_partial} partial genes omitted; --include-partial to keep them)", file=sys.stderr)
    if not rows:
        print("(no --results given: ST, source and block columns will be empty)", file=sys.stderr)
    for f in made:
        print(f, file=sys.stderr)


def cmd_metadata(a):
    ids = []
    for f in a.ids:
        ids += [i for i in _read_ids(f, a.id_column) if i not in ids]
    token_dir = os.path.expanduser(a.token_dir) if a.token_dir else None
    collect_metadata(ids, a.output, token_dir, a.delay)


def cmd_fetch(a):
    ids = _read_ids(a.ids, a.id_column)
    os.makedirs(a.outdir, exist_ok=True)
    failed_path = os.path.join(a.outdir, "failed_ids.txt")
    # genomes already present, whatever their naming (id-152010.fasta, 152010.fas, ...)
    present = {pubmlst_id_from_name(f) for f in os.listdir(a.outdir)
               if f.endswith(FASTA_EXTS) and os.path.getsize(os.path.join(a.outdir, f)) > 0}
    ok = skip = 0
    failed = []  # (id, reason)
    for n, iid in enumerate(ids, 1):
        out = os.path.join(a.outdir, f"id-{iid}.fasta")
        if iid in present:
            skip += 1
            continue
        for attempt in range(3):
            try:
                with urllib.request.urlopen(PUBMLST.format(iid), timeout=120) as resp:
                    data = resp.read()
                if not data.startswith(b">"):
                    raise ValueError("response is not FASTA")
                tmp = out + ".part"          # write to a temporary name, then rename,
                with open(tmp, "wb") as fh:  # so an interrupted download never looks complete
                    fh.write(data)
                os.replace(tmp, out)
                ok += 1
                break
            except urllib.error.HTTPError as e:
                if e.code in (401, 403, 404):
                    # definitive answer from PubMLST (private or unavailable record): do not retry
                    failed.append((iid, f"HTTP {e.code} (not public / not available)"))
                    print(f"[{n}/{len(ids)}] id {iid}: HTTP {e.code}, skipped", file=sys.stderr)
                    break
                if attempt == 2:
                    failed.append((iid, f"HTTP {e.code}"))
                    print(f"[{n}/{len(ids)}] id {iid}: failed (HTTP {e.code})", file=sys.stderr)
                time.sleep(5 * (attempt + 1))
            except (urllib.error.URLError, ValueError, TimeoutError, OSError) as e:
                if attempt == 2:
                    failed.append((iid, str(e)))
                    print(f"[{n}/{len(ids)}] id {iid}: failed ({e})", file=sys.stderr)
                time.sleep(5 * (attempt + 1))
        time.sleep(a.delay)
    with open(failed_path, "w") as fh:
        for iid, why in failed:
            fh.write(f"{iid}\t{why}\n")
    print(f"downloaded {ok}, skipped (already present) {skip}, failed {len(failed)}"
          + (f"  -> failed ids listed in {failed_path}" if failed else ""), file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="campypump", description="Classify RE-cmeB vs wildtype cmeB in Campylobacter genomes.")
    ap.add_argument("--version", action="version", version=f"campypump {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="type genome assemblies")
    r.add_argument("genomes", nargs="+", help="FASTA files or directories of FASTA files")
    r.add_argument("-o", "--output", default="campypump_results.tsv")
    r.add_argument("--alleles", help="write extracted cmeB alleles to this FASTA")
    r.add_argument("-t", "--threads", type=int, default=1)
    r.add_argument("--species", default="", help="label written to a 'species' column, e.g. 'C. coli'")
    r.add_argument("--ani", metavar="TSV",
                   help="use this skani table for the species call instead of running skani")
    r.add_argument("--no-ani", action="store_true",
                   help="skip the whole-genome species call entirely")
    r.add_argument("--segment-plots", metavar="DIR",
                   help="write SVG segment plots for mosaic/ambiguous/divergent/mixed alleles to DIR")
    r.add_argument("--plot-all", action="store_true", help="with --segment-plots: plot every genome")
    r.add_argument("--refs", help="custom reference panel FASTA (default: bundled panel)")
    r.add_argument("--re-min", type=float, default=95.0, help="min %% identity to PA001 RE-cmeB to call RE [95]")
    r.add_argument("--wt-max-re", type=float, default=85.0, help="max %% identity to PA001 to call wildtype [85]")
    r.add_argument("--window", type=int, default=300, help="mosaic scan window, bp [300]")
    r.add_argument("--step", type=int, default=50)
    r.add_argument("--mosaic-frac", type=float, default=0.10,
                   help="min fraction of RE-like AND WT-like windows to flag mosaic [0.10]")
    r.set_defaults(func=cmd_run)

    m = sub.add_parser("metadata", help="get MLST (ST, clonal complex) and provenance from PubMLST")
    m.add_argument("--ids", required=True, nargs="+", help="id files (one PubMLST id per line), or .xlsx")
    m.add_argument("--id-column", default="Isolate id", help="column name when --ids is .xlsx")
    m.add_argument("-o", "--output", default="pubmlst_metadata.tsv")
    m.add_argument("--token-dir", help="BIGSdb_downloader token folder (e.g. ~/.bigsdb_tokens) "
                                       "to include private records; omit for public records only")
    m.add_argument("--delay", type=float, default=0.5, help="seconds between requests [0.5]")
    m.set_defaults(func=cmd_metadata)

    n = sub.add_parser("annotate", help="merge results with PubMLST metadata; flag duplicates")
    n.add_argument("results", help="campypump results table (.tsv or .csv)")
    n.add_argument("--metadata", help="table from 'campypump metadata'")
    n.add_argument("--ani", help="skani table (species_ani.tsv) for the whole-genome species call")
    n.add_argument("--exclude", help="file of ids to mark use_for_prevalence=no "
                                     "(one id per line, optional <tab>reason)")
    n.add_argument("--columns", choices=["core", "full"], default="core",
                   help="core: the ~39 columns analyses need (default); full: every column")
    n.add_argument("-o", "--output", required=True, help="output CSV")
    n.set_defaults(func=lambda a: annotate(a.results, a.metadata, a.output, a.ani, a.exclude,
                                          a.columns))

    fig = sub.add_parser("figure", help="classification figure (panels A and B) from a results table")
    fig.add_argument("results", help="campypump results table (.tsv or .csv)")
    fig.add_argument("-o", "--outdir", default="figures", help="output folder [figures]")
    fig.add_argument("--by-species", action="store_true", default=True,
                     help="one figure per species (default)")
    fig.add_argument("--combined", action="store_true", help="also write one figure for all genomes")
    fig.add_argument("--include-partial", action="store_true",
                     help="also plot genes broken across contigs (their whole-gene identity is deflated)")
    fig.add_argument("--exclude-duplicates", action="store_true",
                     help="skip rows where use_for_prevalence is 'no'")
    fig.set_defaults(func=cmd_figure)

    dv = sub.add_parser("diversity", help="allele diversity figure (panels A-C) from cmeB allele FASTAs")
    dv.add_argument("alleles", nargs="+", help="FASTA files written by 'campypump run --alleles'")
    dv.add_argument("-o", "--output", default="figures/cmeB_allele_diversity.svg")
    dv.add_argument("--labels", nargs="+", help="species label per file (default: from the file name)")
    dv.add_argument("--include-partial", action="store_true",
                    help="count genes broken across contigs as alleles (they inflate the count)")
    dv.add_argument("--table", help="also write the allele table (one row per distinct allele) here")
    dv.set_defaults(func=cmd_diversity)

    al = sub.add_parser("alleles", help="deduplicate cmeB alleles and build tree inputs")
    al.add_argument("fasta", nargs="+", help="allele FASTAs from 'campypump run --alleles'")
    al.add_argument("-o", "--output", default="tree/cmeB_alleles",
                    help="output prefix; writes .fasta, .tsv and iTOL files [tree/cmeB_alleles]")
    al.add_argument("--labels", nargs="+", help="species label per FASTA (default: from the file name)")
    al.add_argument("--results", help="annotated results table, to carry ST/source/block status into the "
                                      "annotation file")
    al.add_argument("--references", help="extra FASTA of reference alleles to include (tagged REF)")
    al.add_argument("--outgroup", help="FASTA of outgroup sequence(s) for rooting (tagged OUT)")
    al.add_argument("--include-partial", action="store_true",
                    help="keep genes broken across contigs (they are fragments, not alleles)")
    al.add_argument("--no-itol", action="store_true", help="skip the iTOL annotation files")
    al.set_defaults(func=cmd_alleles)

    tp = sub.add_parser("treeplot", help="draw an annotated tree (SVG) from a Newick file")
    tp.add_argument("treefile", help="Newick tree, e.g. cmeB_tree.treefile from IQ-TREE")
    tp.add_argument("-o", "--output", default="figures/cmeB_tree.svg")
    tp.add_argument("--outgroup", help="tip name to root on")
    tp.add_argument("--subtree", choices=["RE"], help="draw only the clade holding the RE alleles")
    tp.add_argument("--labels", choices=["notable", "all", "none"], default="notable",
                    help="which tips to label [notable: RE, mosaic, references, outgroup]")
    tp.add_argument("--min-support", type=float, default=90.0,
                    help="show node support at or above this bootstrap value [90]")
    tp.add_argument("--width", type=int, default=1000)
    tp.add_argument("--tip-height", type=float, default=11.0)
    tp.add_argument("--title", default="")
    tp.add_argument("--drop-outgroup", action="store_true",
                    help="remove the outgroup tip after rooting so the ingroup fills the figure")
    tp.add_argument("--annotation", help="table from 'campypump alleles' (<prefix>.tsv); adds the "
                                        "dominant ST to each labelled tip")
    tp.add_argument("--annotate-cols", nargs="+", default=["top_STs"],
                    help="columns from --annotation to append to labels [top_STs]")
    tp.add_argument("--no-size-by-genomes", action="store_true",
                    help="draw every tip the same size instead of scaling by genomes carrying it")
    tp.add_argument("--no-source-bars", action="store_true",
                    help="omit the per-allele source composition bars")
    tp.add_argument("--cladogram", action="store_true",
                    help="ignore branch lengths; useful when one clade is far more divergent")
    tp.add_argument("--colour-by", choices=["auto", "type", "species"], default="auto",
                    help="colour branches by cmeB type or species [auto: species inside the RE clade]")
    tp.set_defaults(func=lambda a: print(tree_figure(
        a.treefile, a.output, a.outgroup, a.labels, a.min_support, a.width, a.tip_height,
        a.subtree, a.title, a.colour_by, a.drop_outgroup, a.cladogram,
        a.annotation, tuple(a.annotate_cols), not a.no_size_by_genomes,
        not a.no_source_bars), file=sys.stderr))

    c = sub.add_parser("combine", help="combine several result tables into one")
    c.add_argument("tables", nargs="+")
    c.add_argument("-o", "--output", required=True)
    c.set_defaults(func=cmd_combine)

    f = sub.add_parser("fetch", help="download assemblies from PubMLST by isolate id")
    f.add_argument("--ids", required=True, help="text file with one PubMLST isolate id per line, or .xlsx")
    f.add_argument("--id-column", default="Isolate id", help="column name when --ids is .xlsx")
    f.add_argument("--outdir", default="genomes")
    f.add_argument("--delay", type=float, default=1.0, help="seconds between requests [1]")
    f.set_defaults(func=cmd_fetch)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
