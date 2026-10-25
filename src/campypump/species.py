"""Species assignment from whole-genome ANI (skani output).

No single cme gene identifies the species: the operon crosses between *C. coli* and *C. jejuni*,
and *C. coli* DH161 carries a 99.4 % *C. jejuni*-like cmeR. Whole-genome ANI is unaffected by a
single transferred locus, so it is what campypump uses for the species call; the operon columns
describe the operon only.

Input is the table written by:

    skani dist -q genomes/*/*.fasta -r species_refs/*.fasta -t 8 -o species_ani.tsv

Reference file names supply the species, so name them for the organism
(`Ccoli_RM4661.fasta`, `C.jejuni_NCTC11168.fasta`, `Cupsaliensis_RM3940.fasta`, ...).

Since 0.4.0 campypump ships those four genomes and runs skani itself, so `campypump run`
gives a species call without a separate step. `run_skani` below is that call; everything
after it is unchanged, so a table produced by hand still works with `annotate --ani`.
"""
from __future__ import annotations

import csv
import glob
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from importlib import resources

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


# a reference named for its organism but not in the list above, e.g. Clanienae_NCTC13004.fasta
# or C.hepaticus_HV10.fna -> "C. lanienae", "C. hepaticus"
_GENERIC = re.compile(r"^C\.?_?([a-z]{4,})(?:[_.-]|$)")


def species_from_reference(ref_file: str, ref_name: str = "") -> str:
    """Work out the species from the reference file name, falling back to its FASTA header."""
    base = os.path.basename(ref_file)
    for text in (base, ref_name):
        for pattern, species in _SPECIES_PATTERNS:
            if pattern.search(text or ""):
                return species
    m = _GENERIC.match(base)
    if m:
        return f"C. {m.group(1)}"
    return os.path.splitext(base)[0]


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


# --------------------------------------------------------------------------- #
# running skani
# --------------------------------------------------------------------------- #

SPECIES_REF_ENV = "CAMPYPUMP_SPECIES_REFS"
SKANI_MISSING_HINT = (
    "skani was not found on PATH, so no whole-genome species call was made.\n"
    "    install it:   conda install -c bioconda skani\n"
    "    or skip it:   --no-ani\n"
    "    or supply a table you made yourself:   --ani species_ani.tsv"
)


def skani_available() -> bool:
    return shutil.which("skani") is not None


def bundled_species_refs() -> str:
    """Directory holding the four species reference genomes, as plain FASTA.

    They ship gzipped inside the package and are unpacked once into a cache directory,
    so skani is always handed plain files and nothing depends on which input formats a
    given skani build accepts. Set CAMPYPUMP_SPECIES_REFS to use your own directory
    instead; file names supply the species, so keep the naming convention.
    """
    override = os.environ.get(SPECIES_REF_ENV)
    if override:
        return override
    src = resources.files("campypump").joinpath("data/species")
    cache = os.path.join(os.environ.get("XDG_CACHE_HOME")
                         or os.path.join(os.path.expanduser("~"), ".cache"),
                         "campypump", "species_refs")
    os.makedirs(cache, exist_ok=True)
    for entry in src.iterdir():
        if not entry.name.endswith(".fasta.gz"):
            continue
        out = os.path.join(cache, entry.name[:-3])
        if os.path.exists(out) and os.path.getsize(out) > 0:
            continue
        # write to a temporary name first: two campypump runs may start at once, and a
        # half-written reference is worse than none
        fd, tmp = tempfile.mkstemp(dir=cache, suffix=".part")
        os.close(fd)
        with entry.open("rb") as raw, gzip.open(raw, "rb") as fin, open(tmp, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        os.replace(tmp, out)
    return cache


def species_reference_files(refs_dir: str | None = None) -> list[str]:
    d = refs_dir or bundled_species_refs()
    return sorted(f for ext in ("*.fasta", "*.fa", "*.fna")
                  for f in glob.glob(os.path.join(d, ext)))


def run_skani(genomes: list[str], out_path: str, threads: int = 1,
              refs_dir: str | None = None) -> str:
    """Run `skani dist` of these genomes against the species references.

    Query and reference paths go in list files rather than on the command line: a run of
    several thousand genomes would otherwise overflow the shell's argument limit.
    """
    refs = species_reference_files(refs_dir)
    if not refs:
        raise FileNotFoundError(
            f"no species reference genomes in {refs_dir or bundled_species_refs()}")
    work = tempfile.mkdtemp(prefix="campypump_skani_")
    try:
        ql, rl = os.path.join(work, "query.txt"), os.path.join(work, "ref.txt")
        for path, items in ((ql, genomes), (rl, refs)):
            with open(path, "w") as fh:
                fh.write("\n".join(os.path.abspath(i) for i in items) + "\n")
        cmd = ["skani", "dist", "--ql", ql, "--rl", rl,
               "-t", str(max(1, threads)), "-o", out_path]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError("skani failed (" + " ".join(cmd) + "):\n"
                               + (res.stderr or res.stdout).strip())
        return out_path
    finally:
        shutil.rmtree(work, ignore_errors=True)
