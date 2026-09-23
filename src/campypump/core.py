"""Core functions for campypump.

Pipeline per genome:
  1. BLAST (blastn -task blastn) the bundled reference panel against the assembly
     to locate cmeR, cmeA, cmeB and cmeC.
  2. Extract each gene as a full ORF (anchored on the reference start, extended to
     the first in-frame stop codon) and run QC.
  3. Globally align the extracted cmeB to every cmeB reference using the same
     scoring as Geneious Prime's default global alignment (5/-4, gap open 12,
     gap extend 3). % identity = identical sites / alignment length (gaps
     included), which reproduces the Geneious values in the Peru spreadsheet.
  4. Classify cmeB as RE-cmeB / WT-cmeB / ambiguous, and scan along the gene in
     windows to flag mosaic (recombinant) alleles.
  5. Report cmeR, cmeA and cmeC identity to the 11168 and PA001 references, and
     the sequence of the CmeR-binding inverted repeat in the cmeR-cmeA intergenic
     region.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from importlib import resources

from Bio import Align, SeqIO
from Bio.Seq import Seq

GENES = ("cmeR", "cmeA", "cmeB", "cmeC")
STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}

# 16-bp inverted repeat in the NCTC 11168 cmeR-cmeA intergenic region
IR_11168 = "TGTAATAAAAATTACA"

PRIMARY_WT = "NCTC11168"
PRIMARY_RE = "PA001"
# published RE-cmeABC source: C. coli DH161, GenBank KT778508.1 (Yao et al. 2016).
# Cooper et al. 2025 define RE-cmeB as >90% nucleotide identity to DH161 cmeB.
PUBLISHED_RE = "DH161-KT778508"
PUBLISHED_RE_MIN = 90.0


@dataclass
class Ref:
    gene: str
    cls: str       # Cj_WT, Cj_RE, Cc_WT, Cc_RE
    name: str
    role: str      # primary or panel
    seq: str

    @property
    def label(self) -> str:
        return f"{self.cls}_{self.name}"


@dataclass
class Params:
    re_min: float = 95.0          # >= this % identity to PA001 RE-cmeB -> RE
    wt_max_re: float = 85.0       # <= this % identity to PA001 -> wildtype
    window: int = 300             # mosaic scan window (bp)
    step: int = 50
    window_diff: float = 5.0      # window counts as RE-like / WT-like if it wins by this margin
    mosaic_min_frac: float = 0.10  # both window types must reach this fraction to flag mosaic
    min_len_frac: float = 0.95    # extracted ORF must be >= this fraction of reference length
    partial_min_frac: float = 0.90  # partial gene: provisional call if >= this share of segments agree
    divergent_max: float = 90.0     # complete gene < this % identity to every reference -> divergent


@dataclass
class GeneHit:
    gene: str
    seq: str = ""
    contig: str = ""
    start: int = 0   # 1-based, genomic coordinates
    end: int = 0
    strand: str = ""
    flags: list = field(default_factory=list)


def load_refs(path: str | None = None) -> list[Ref]:
    if path is None:
        path = str(resources.files("campypump").joinpath("data/references.fasta"))
    refs = []
    for r in SeqIO.parse(path, "fasta"):
        gene, cls, name, role = r.id.split("|")
        refs.append(Ref(gene, cls, name, role, str(r.seq).upper()))
    return refs


def check_blast() -> None:
    if shutil.which("blastn") is None:
        raise SystemExit("blastn not found on PATH. Install BLAST+ (e.g. conda install -c bioconda blast).")


# ---------------------------------------------------------------- alignment
def geneious_aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score = 5
    al.mismatch_score = -4
    al.open_gap_score = -12
    al.extend_gap_score = -3
    return al


_ALIGNER = geneious_aligner()


def pct_identity(ref: str, query: str) -> tuple[float, object]:
    """Geneious-style global % identity: identical sites / alignment length (gaps counted)."""
    aln = _ALIGNER.align(ref, query)[0]
    c = aln.counts()
    length = c.identities + c.mismatches + c.gaps
    return 100.0 * c.identities / length, aln


def _semiglobal_aligner() -> Align.PairwiseAligner:
    """Geneious scoring, but gaps at the ends of the query are free (for partial genes)."""
    al = geneious_aligner()
    try:  # Biopython >= 1.86 names (gaps in the query = deletions)
        al.open_left_deletion_score = al.extend_left_deletion_score = 0
        al.open_right_deletion_score = al.extend_right_deletion_score = 0
    except AttributeError:  # older Biopython
        al.query_left_open_gap_score = al.query_left_extend_gap_score = 0
        al.query_right_open_gap_score = al.query_right_extend_gap_score = 0
    return al


_SEMIGLOBAL = _semiglobal_aligner()


def aligned_identity(ref: str, query: str) -> float:
    """% identity over the part of the reference covered by the query (terminal gaps ignored).

    Equals the global identity for complete genes; for a gene fragment it gives the identity
    of the piece that is present, which is how Geneious reports partial hits."""
    aln = _SEMIGLOBAL.align(ref, query)[0]
    ref_idx, q_idx = aln.indices
    cols = [k for k, qi in enumerate(q_idx) if qi >= 0]
    if not cols:
        return 0.0
    ident = length = 0
    for k in range(cols[0], cols[-1] + 1):
        ri, qi = ref_idx[k], q_idx[k]
        length += 1
        if ri >= 0 and qi >= 0 and ref[ri] == query[qi]:
            ident += 1
    return 100.0 * ident / length


def query_match_vector(aln) -> list[int]:
    """1/0 per query position: identical to the aligned reference base or not."""
    ref_idx, q_idx = aln.indices
    ref, q = aln.sequences
    vec = [0] * len(q)
    for ri, qi in zip(ref_idx, q_idx):
        if qi >= 0 and ri >= 0 and ref[ri] == q[qi]:
            vec[qi] = 1
    return vec


# ---------------------------------------------------------------- locating genes
def blast_hits(refs: list[Ref], genome: str, tmpdir: str) -> dict[str, tuple]:
    """Best hit per gene (highest bitscore across all references for that gene)."""
    qfa = os.path.join(tmpdir, "q.fa")
    with open(qfa, "w") as fh:
        for i, r in enumerate(refs):
            if r.gene in GENES:
                fh.write(f">{i}\n{r.seq}\n")
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
    return best


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
    orf = str(s[start:i + 3 if stop_found else i]).upper()
    hit.seq, hit.strand = orf, strand
    if strand == "+":
        hit.start, hit.end = start + 1, start + len(orf)
    else:
        hit.start, hit.end = L - (start + len(orf)) + 1, L - start
    if orf[:3] not in STARTS:
        hit.flags.append("no_start_codon")
    prot = str(Seq(orf[: len(orf) - len(orf) % 3]).translate())
    if "*" in prot[:-1]:
        hit.flags.append("internal_stop")
    return hit


# ---------------------------------------------------------------- analysis
def mosaic_scan(query: str, re_aln, wt_alns: list, p: Params, windows: list | None = None) -> tuple[float, float]:
    """Fraction of windows where query is clearly closer to RE vs to the best WT reference.

    If a list is passed as `windows`, (window midpoint, RE identity, best WT identity) is
    appended for every window (used for segment plots)."""
    re_vec = query_match_vector(re_aln)
    wt_vecs = [query_match_vector(a) for a in wt_alns]
    n = len(query)
    re_like = wt_like = total = 0
    for w0 in range(0, max(n - p.window, 0) + 1, p.step):
        w1 = min(w0 + p.window, n)
        re_id = 100 * sum(re_vec[w0:w1]) / (w1 - w0)
        wt_id = max(100 * sum(v[w0:w1]) / (w1 - w0) for v in wt_vecs)
        total += 1
        if windows is not None:
            windows.append(((w0 + w1) // 2, re_id, wt_id))
        if re_id - wt_id >= p.window_diff:
            re_like += 1
        elif wt_id - re_id >= p.window_diff:
            wt_like += 1
    return (re_like / total, wt_like / total) if total else (0.0, 0.0)


def promoter_ir(contig_seq: Seq, cmeR: GeneHit | None, cmeA: GeneHit, promoter_ref: str) -> tuple[str, str]:
    """Return (intergenic sequence, IR sequence aligned to the 11168 IR)."""
    if cmeA.strand == "+":
        a_start = cmeA.start - 1
        if cmeR and cmeR.contig == cmeA.contig and cmeR.strand == "+" and cmeR.end < cmeA.start:
            region = contig_seq[cmeR.end:a_start]
        else:
            region = contig_seq[max(a_start - 150, 0):a_start]
    else:
        a_start = cmeA.end
        if cmeR and cmeR.contig == cmeA.contig and cmeR.strand == "-" and cmeR.start > cmeA.end:
            region = contig_seq[a_start:cmeR.start - 1].reverse_complement()
        else:
            region = contig_seq[a_start:a_start + 150].reverse_complement()
    region = str(region).upper()
    if not region:
        return "", ""
    al = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3,
                               open_gap_score=-5, extend_gap_score=-2)
    aln = al.align(promoter_ref, region)[0]
    ir0 = promoter_ref.find(IR_11168)
    ir1 = ir0 + len(IR_11168)
    ref_idx, q_idx = aln.indices
    # alignment columns spanning the IR (query insertions inside it are kept)
    cols = [k for k, ri in enumerate(ref_idx) if ir0 <= ri < ir1]
    if not cols:
        return region, ""
    c0, c1 = cols[0], cols[-1]
    ir_seq = "".join(region[q_idx[k]] for k in range(c0, c1 + 1) if q_idx[k] >= 0)
    return region, ir_seq


# plain-language descriptions for QC / analysis flags
WARNING_TEXT = {
    "contig_edge_5prime": "gene starts at contig edge (5' end missing)",
    "contig_edge_3prime": "gene runs off contig end (3' end missing)",
    "truncated": "gene shorter than 95% of reference",
    "internal_stop": "premature stop codon",
    "no_start_codon": "no start codon",
    "mosaic": "mixed RE-like and wildtype-like segments",
    "divergent": "divergent allele: <90% identity to every reference (check species / novel allele)",
    "partial_provisional": "partial gene: type is provisional, based on the part present",
    "partial_unresolved": "partial gene: part present is not clearly RE or wildtype",
}


def _warn(flag: str) -> str:
    for gene in ("cmeR", "cmeA", "cmeC"):
        if flag == f"{gene}_not_found":
            return f"{gene} not found"
        if flag.startswith(gene + "_"):
            return f"{gene}: {WARNING_TEXT.get(flag[len(gene) + 1:], flag)}"
    return WARNING_TEXT.get(flag, flag)


FASTA_EXTS = (".fasta", ".fas", ".fa", ".fna", ".fsa", ".fasta.txt")


def pubmlst_id_from_name(name: str) -> str:
    """id-152010.fasta, 152010.fas or 152010_CAA03704.fasta -> '152010'."""
    m = re.match(r"^(?:id-)?(\d+)(?:[._-]|$)", os.path.basename(name))
    return m.group(1) if m else ""


def panel_column(ref: Ref) -> str:
    """e.g. Cc_WT + CAA03009-008347 -> panel_CcWT_CAA03009_pct (Stata-safe, <=32 chars)."""
    strain = ref.name.split("-")[0].replace(".", "")
    return f"panel_{ref.cls.replace('_', '')}_{strain}_pct"[:32]


def classify_genome(genome: str, refs: list[Ref], p: Params) -> dict:
    name = os.path.basename(genome)
    for ext in FASTA_EXTS:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    row: dict = {"sample": name,
                 "pubmlst_id": pubmlst_id_from_name(name)}
    contigs = {r.id: r.seq for r in SeqIO.parse(genome, "fasta")}
    with tempfile.TemporaryDirectory() as td:
        best = blast_hits(refs, genome, td)

    hits: dict[str, GeneHit] = {}
    for gene in GENES:
        if gene not in best:
            continue
        ref, contig, qs, qe, ss, se, _ = best[gene]
        h = extract_orf(contigs[contig], ref, qs, ss, se)
        h.contig = contig
        if len(h.seq) < p.min_len_frac * len(ref.seq):
            h.flags.append("truncated")
        hits[gene] = h

    cmeB_refs = [r for r in refs if r.gene == "cmeB"]
    b = hits.get("cmeB")
    if b is None:
        row.update(cmeB_type="not_found", cmeB_type_final="not_found", re_cmeB="",
                   warnings="cmeB not found in assembly")
        return row

    ids, alns = {}, {}
    for r in cmeB_refs:
        ids[r.label], alns[r.label] = pct_identity(r.seq, b.seq)
    re_ref = next(r for r in cmeB_refs if r.name == PRIMARY_RE)
    wt_ref = next(r for r in cmeB_refs if r.name == PRIMARY_WT)
    dh_ref = next((r for r in cmeB_refs if r.name == PUBLISHED_RE), None)
    pid_re, pid_wt = ids[re_ref.label], ids[wt_ref.label]
    nearest = max(cmeB_refs, key=lambda r: ids[r.label])

    wt_labels = [r.label for r in cmeB_refs if r.cls.endswith("WT")]
    windows: list = []
    re_frac, wt_frac = mosaic_scan(b.seq, alns[re_ref.label], [alns[l] for l in wt_labels], p, windows)

    flags = list(b.flags)
    if re_frac >= p.mosaic_min_frac and wt_frac >= p.mosaic_min_frac:
        flags.append("mosaic")
    if any(f in b.flags for f in ("truncated", "contig_edge_5prime", "contig_edge_3prime", "internal_stop")):
        call = "incomplete"
    elif "mosaic" in flags:
        call = "mosaic"
    elif pid_re >= p.re_min:
        call = "RE-cmeB"
    elif pid_re <= p.wt_max_re:
        call = "wildtype-cmeB"
    else:
        call = "ambiguous"

    # final type: partial genes get a provisional call when the part present is clearly one type
    partial = call == "incomplete"
    divergent = (not partial) and ids[nearest.label] < p.divergent_max
    if partial:
        if re_frac >= p.partial_min_frac:
            final = "RE-cmeB"
            flags.append("partial_provisional")
        elif wt_frac >= p.partial_min_frac:
            final = "wildtype-cmeB"
            flags.append("partial_provisional")
        else:
            final = "unresolved"
            flags.append("partial_unresolved")
    else:
        final = call
    if divergent:
        flags.append("divergent")

    row.update(
        cmeB_type=call,
        cmeB_type_final=final,
        re_cmeB={"RE-cmeB": 1, "wildtype-cmeB": 0}.get(final, ""),
        cmeB_partial_gene="yes" if partial else "no",
        cmeB_divergent="yes" if divergent else "no",
        cmeB_identity_RE_PA001_pct=round(pid_re, 3),
        cmeB_identity_WT_11168_pct=round(pid_wt, 3),
        cmeB_closest_reference=f"{nearest.name} ({nearest.cls.replace('_', ' ')})",
        cmeB_closest_identity_pct=round(ids[nearest.label], 3),
        cmeB_RE_like_segments_pct=round(100 * re_frac, 1),
        cmeB_WT_like_segments_pct=round(100 * wt_frac, 1),
        cmeB_identity_DH161_pct=round(ids[dh_ref.label], 3) if dh_ref else "",
        # the published >90% rule needs identity over the part present, or a fragment's
        # deflated whole-gene identity would read as "not RE"
        RE_by_DH161_definition=("" if dh_ref is None else
                                "yes" if (aligned_identity(dh_ref.seq, b.seq) if partial
                                          else ids[dh_ref.label]) > PUBLISHED_RE_MIN else "no"),
        cmeB_aligned_identity_RE_pct=round(aligned_identity(re_ref.seq, b.seq), 3),
        cmeB_aligned_identity_WT_pct=round(aligned_identity(wt_ref.seq, b.seq), 3),
    )

    # other operon genes: identity to the 11168 and PA001 versions
    for gene in ("cmeR", "cmeA", "cmeC"):
        h = hits.get(gene)
        col_wt, col_re = f"{gene}_identity_WT_11168_pct", f"{gene}_identity_RE_PA001_pct"
        if h is None:
            row[col_wt] = row[col_re] = ""
            flags.append(f"{gene}_not_found")
            continue
        for r in refs:
            if r.gene == gene and r.role == "primary":
                row[col_re if r.name == PRIMARY_RE else col_wt] = round(pct_identity(r.seq, h.seq)[0], 3)
        flags += [f"{gene}_{f}" for f in h.flags]

    # promoter / CmeR-binding inverted repeat
    a = hits.get("cmeA")
    prom_ref = next((r.seq for r in refs if r.gene == "promoter"), None)
    if a is not None and prom_ref:
        region, ir = promoter_ir(contigs[a.contig], hits.get("cmeR"), a, prom_ref)
        row.update(promoter_length_bp=len(region), CmeR_binding_site_seq=ir,
                   CmeR_binding_site_as_11168="yes" if ir == IR_11168 else "no")

    row.update(cmeB_length_bp=len(b.seq), cmeB_contig=b.contig, cmeB_start=b.start,
               cmeB_end=b.end, cmeB_strand=b.strand)
    for r in cmeB_refs:
        row[panel_column(r)] = round(ids[r.label], 3)
    row["warnings"] = "; ".join(_warn(f) for f in flags)
    row["_cmeB_seq"] = b.seq
    row["_windows"] = windows
    return row
