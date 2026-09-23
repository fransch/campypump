"""Figures summarising cmeB typing across a collection (plain SVG, no plotting dependencies).

Panel A  identity to RE-cmeB (PA001) vs identity to wildtype cmeB (NCTC 11168), one point per
         genome, coloured by cmeB type. Shows the separation the classification relies on.
Panel B  distribution of identity to PA001, in 0.5 % bins, split by cmeB type.

Panels carry letters only; describing what they show is the figure legend's job.
"""
from __future__ import annotations

import os
from collections import Counter

from .palettes import SPECIES_LINE, TYPE as COLORS
ORDER = ["RE-cmeB", "wildtype-cmeB", "mosaic", "ambiguous", "unresolved", "not_found"]

W, H = 900, 430
AX = dict(left=62, right=18, top=34, bottom=52)       # panel A margins
PANEL_W = (W - 40) // 2


def _fmt(v: float) -> str:
    return f"{v:g}"


def _esc(s: str) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class _Axes:
    """Linear axes for one panel, with its own origin inside the figure."""

    def __init__(self, x0, y0, w, h, xlim, ylim):
        self.x0, self.y0, self.w, self.h = x0, y0, w, h
        self.xlim, self.ylim = xlim, ylim

    def x(self, v):
        lo, hi = self.xlim
        return self.x0 + self.w * (min(max(v, lo), hi) - lo) / (hi - lo)

    def y(self, v):
        lo, hi = self.ylim
        return self.y0 + self.h - self.h * (min(max(v, lo), hi) - lo) / (hi - lo)

    def frame(self, xticks, yticks, xlabel, ylabel, xfmt=_fmt, yfmt=_fmt):
        out = [f'<rect x="{self.x0}" y="{self.y0}" width="{self.w}" height="{self.h}" '
               f'fill="none" stroke="#CCCCCC" stroke-width="1"/>']
        for t in xticks:
            x = self.x(t)
            out.append(f'<line x1="{x:.1f}" y1="{self.y0 + self.h}" x2="{x:.1f}" '
                       f'y2="{self.y0 + self.h + 4}" stroke="#666"/>')
            out.append(f'<text x="{x:.1f}" y="{self.y0 + self.h + 17}" text-anchor="middle" '
                       f'font-size="11" fill="#333">{xfmt(t)}</text>')
        for t in yticks:
            y = self.y(t)
            out.append(f'<line x1="{self.x0 - 4}" y1="{y:.1f}" x2="{self.x0}" y2="{y:.1f}" stroke="#666"/>')
            out.append(f'<text x="{self.x0 - 7}" y="{y + 4:.1f}" text-anchor="end" '
                       f'font-size="11" fill="#333">{yfmt(t)}</text>')
        out.append(f'<text x="{self.x0 + self.w / 2:.0f}" y="{self.y0 + self.h + 38}" '
                   f'text-anchor="middle" font-size="12" fill="#222">{_esc(xlabel)}</text>')
        out.append(f'<text transform="translate({self.x0 - 44},{self.y0 + self.h / 2:.0f}) rotate(-90)" '
                   f'text-anchor="middle" font-size="12" fill="#222">{_esc(ylabel)}</text>')
        return "".join(out)


def _panel_letter(x, y, letter):
    return (f'<text x="{x}" y="{y}" font-size="15" font-weight="bold" fill="#111">{letter}</text>')


def scatter_panel(rows, x0, y0, w, h):
    ax = _Axes(x0, y0, w, h, (70, 101), (70, 101))
    parts = [ax.frame(range(70, 101, 5), range(70, 101, 5),
                      "% identity to RE-cmeB (PA001)", "% identity to cmeB (NCTC 11168)",
                      lambda t: f"{t}", lambda t: f"{t}")]
    for cut in (85, 95):
        x = ax.x(cut)
        parts.append(f'<line x1="{x:.1f}" y1="{y0}" x2="{x:.1f}" y2="{y0 + h}" stroke="#999" '
                     f'stroke-width="1" stroke-dasharray="4,3"/>')
        parts.append(f'<text x="{x + 3:.1f}" y="{y0 + 12}" font-size="10" fill="#777">{cut}%</text>')
    # draw wildtype first so rarer classes stay visible on top
    for cls in ["wildtype-cmeB", "RE-cmeB", "ambiguous", "unresolved", "mosaic"]:
        pts = [r for r in rows if r["cls"] == cls]
        if not pts:
            continue
        colour = COLORS.get(cls, "#7F7F7F")
        # rarer classes get a larger, outlined marker
        big = cls in ("mosaic", "ambiguous", "unresolved") or len(pts) <= 5
        r = 4.2 if big else 2.4
        op = 1.0 if big else 0.45
        stroke = ' stroke="#222" stroke-width="0.9"' if big else ""
        parts.append(f'<g fill="{colour}" fill-opacity="{op}"{stroke}>')
        parts += [f'<circle cx="{ax.x(p["re"]):.1f}" cy="{ax.y(p["wt"]):.1f}" r="{r}"/>' for p in pts]
        parts.append("</g>")
    return "".join(parts)


def histogram_panel(rows, x0, y0, w, h, binw=0.5):
    lo, hi = 70.0, 101.0
    nbins = int((hi - lo) / binw)
    counts = {cls: [0] * nbins for cls in ORDER}
    for r in rows:
        b = min(int((r["re"] - lo) / binw), nbins - 1)
        if 0 <= b:
            counts.setdefault(r["cls"], [0] * nbins)[b] += 1
    top = max((max(v) for v in counts.values() if v), default=1)
    # log-ish scale keeps single genomes visible next to bins of hundreds
    import math
    ymax = math.ceil(math.log10(max(top, 10)) * 10) / 10
    floor = math.log10(0.55)      # keep single-genome bins visible on the log axis
    ax = _Axes(x0, y0, w, h, (lo, hi), (floor, math.log10(top) + 0.25))
    ticks = [t for t in (1, 10, 100, 1000) if t <= top * 1.5]
    parts = [ax.frame(range(70, 101, 5), [math.log10(t) for t in ticks],
                      "% identity to RE-cmeB (PA001)", "genomes",
                      lambda t: f"{t}", lambda v: f"{int(round(10 ** v))}")]
    for cut in (85, 95):
        x = ax.x(cut)
        parts.append(f'<line x1="{x:.1f}" y1="{y0}" x2="{x:.1f}" y2="{y0 + h}" stroke="#999" '
                     f'stroke-width="1" stroke-dasharray="4,3"/>')
        parts.append(f'<text x="{x + 3:.1f}" y="{y0 + 12}" font-size="10" fill="#777">{cut}%</text>')
    bw = w / nbins
    for cls in ["wildtype-cmeB", "RE-cmeB", "mosaic", "ambiguous", "unresolved"]:
        v = counts.get(cls)
        if not v or not max(v):
            continue
        parts.append(f'<g fill="{COLORS.get(cls, "#7F7F7F")}" fill-opacity="0.85">')
        for i, c in enumerate(v):
            if not c:
                continue
            ytop = ax.y(math.log10(c))
            height = max(y0 + h - ytop, 2.0)
            parts.append(f'<rect x="{x0 + i * bw:.2f}" y="{y0 + h - height:.1f}" '
                         f'width="{max(bw - 0.4, 0.8):.2f}" height="{height:.1f}"/>')
        parts.append("</g>")
    return "".join(parts)


def legend(x, y, classes, counts):
    parts = []
    for cls in ORDER:
        if cls not in classes:
            continue
        parts.append(f'<circle cx="{x + 5}" cy="{y - 4}" r="4.2" fill="{COLORS.get(cls, "#7F7F7F")}"/>')
        label = f"{cls} (n = {counts[cls]:,})"
        parts.append(f'<text x="{x + 15}" y="{y}" font-size="11.5" fill="#222">{_esc(label)}</text>')
        x += 16 + 7.0 * len(label)
    return "".join(parts)


SPECIES_COLORS = SPECIES_LINE


def _log_axes(x0, y0, w, h, xlim, ylim_log):
    """Axes with a log10 y scale; ylim_log given as (min_count, max_count)."""
    import math
    lo, hi = ylim_log
    return _Axes(x0, y0, w, h, xlim, (math.log10(max(lo, 0.9)), math.log10(hi) + 0.05))


def rank_abundance_panel(series, x0, y0, w, h):
    """series: [(label, colour, [counts sorted desc])]"""
    import math
    nmax = max(len(c) for _, _, c in series)
    cmax = max(max(c) for _, _, c in series)
    ax = _log_axes(x0, y0, w, h, (0, nmax + 2), (1, cmax))
    ticks = [t for t in (1, 10, 100, 1000) if t <= cmax * 1.5]
    xt = [t for t in range(0, nmax + 1, 50)]
    parts = [ax.frame(xt, [math.log10(t) for t in ticks],
                      "allele rank", "genomes carrying the allele",
                      lambda t: f"{t}", lambda v: f"{int(round(10 ** v))}")]
    for label, colour, counts in series:
        pts = " ".join(f"{ax.x(i + 1):.1f},{ax.y(math.log10(c)):.1f}" for i, c in enumerate(counts))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{colour}" stroke-width="2"/>')
    return "".join(parts)


def accumulation_panel(series, x0, y0, w, h):
    """series: [(label, colour, [(n_genomes, mean_unique_alleles), ...])]"""
    nmax = max(p[0] for _, _, pts in series for p in pts)
    amax = max(p[1] for _, _, pts in series for p in pts)
    ax = _Axes(x0, y0, w, h, (0, nmax * 1.02), (0, amax * 1.08))
    xstep = 500 if nmax > 1200 else 250
    parts = [ax.frame(range(0, int(nmax) + 1, xstep), range(0, int(amax) + 1, 50),
                      "genomes sampled", "distinct alleles",
                      lambda t: f"{t:,}", lambda t: f"{t}")]
    for label, colour, pts in series:
        line = " ".join(f"{ax.x(n):.1f},{ax.y(a):.1f}" for n, a in pts)
        parts.append(f'<polyline points="{line}" fill="none" stroke="{colour}" stroke-width="2"/>')
    return "".join(parts)


def alleles_per_class_panel(groups, x0, y0, w, h):
    """groups: [(species, colour, [(cls, n_genomes, n_alleles), ...])]"""
    bars = [(sp, colour, cls, g, a) for sp, colour, rows in groups for cls, g, a in rows]
    amax = max(b[4] for b in bars)
    ax = _Axes(x0, y0, w, h, (0, len(bars)), (0, amax * 1.18))
    parts = [ax.frame([], range(0, int(amax) + 1, 50), "", "distinct alleles",
                      lambda t: "", lambda t: f"{t}")]
    bw = w / len(bars) * 0.62
    for i, (sp, colour, cls, ngen, nall) in enumerate(bars):
        cx = ax.x(i + 0.5)
        ytop = ax.y(nall)
        parts.append(f'<rect x="{cx - bw / 2:.1f}" y="{ytop:.1f}" width="{bw:.1f}" '
                     f'height="{y0 + h - ytop:.1f}" fill="{COLORS.get(cls, "#7F7F7F")}" fill-opacity="0.9"/>')
        parts.append(f'<text x="{cx:.1f}" y="{ytop - 14:.1f}" text-anchor="middle" font-size="10.5" '
                     f'fill="#222">{nall}</text>')
        parts.append(f'<text x="{cx:.1f}" y="{ytop - 4:.1f}" text-anchor="middle" font-size="9" '
                     f'fill="#666">{nall / ngen * 100:.0f} per 100 genomes</text>')
        short = cls.replace("-cmeB", "").replace("wildtype", "WT").replace("RE", "RE")
        parts.append(f'<text x="{cx:.1f}" y="{y0 + h + 15}" text-anchor="middle" font-size="10.5" '
                     f'fill="#222">{_esc(short)}</text>')
        parts.append(f'<text x="{cx:.1f}" y="{y0 + h + 27}" text-anchor="middle" font-size="9.5" '
                     f'fill="{colour}">{_esc(sp)}</text>')
    return "".join(parts)


def diversity_figure(per_species, outfile, note=""):
    """per_species: [(species, counts_desc, accumulation_pts, [(cls, n_genomes, n_alleles)])]"""
    w, h = 1180, 380
    left, right, top, bottom = 58, 16, 40, 62
    pw = (w - left - right - 2 * 46) // 3
    ph = h - top - bottom
    series_rank, series_acc, groups = [], [], []
    for i, (sp, counts, acc, byclass) in enumerate(per_species):
        colour = SPECIES_COLORS[i % len(SPECIES_COLORS)]
        series_rank.append((sp, colour, counts))
        series_acc.append((sp, colour, acc))
        groups.append((sp, colour, byclass))
    xs = [left, left + pw + 46, left + 2 * (pw + 46)]
    body = [
        _panel_letter(14, 20, "A"), rank_abundance_panel(series_rank, xs[0], top, pw, ph),
        _panel_letter(xs[1] - 44, 20, "B"), accumulation_panel(series_acc, xs[1], top, pw, ph),
        _panel_letter(xs[2] - 44, 20, "C"), alleles_per_class_panel(groups, xs[2], top, pw, ph),
    ]
    lx = left
    for sp, colour, _ in groups:
        body.append(f'<line x1="{lx}" y1="{h - 12}" x2="{lx + 22}" y2="{h - 12}" stroke="{colour}" stroke-width="2.5"/>')
        body.append(f'<text x="{lx + 28}" y="{h - 8}" font-size="11.5" font-style="italic" fill="#222">{_esc(sp)}</text>')
        lx += 40 + 7.2 * len(sp)
    if note:
        body.append(f'<text x="{w - right}" y="19" text-anchor="end" font-size="12" fill="#333">{_esc(note)}</text>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" '
           f'font-family="Helvetica, Arial, sans-serif">'
           f'<rect width="100%" height="100%" fill="white"/>{"".join(body)}</svg>')
    os.makedirs(os.path.dirname(os.path.abspath(outfile)) or ".", exist_ok=True)
    with open(outfile, "w") as fh:
        fh.write(svg)
    return outfile


def classification_figure(rows: list[dict], outfile: str, subtitle: str = "") -> str:
    """rows: dicts with 're' (% id to PA001), 'wt' (% id to 11168), 'cls' (cmeB type)."""
    counts = Counter(r["cls"] for r in rows)
    pw = PANEL_W - AX["left"] - AX["right"]
    ph = H - AX["top"] - AX["bottom"] - 26
    body = [
        _panel_letter(14, 20, "A"),
        scatter_panel(rows, AX["left"], AX["top"], pw, ph),
        _panel_letter(PANEL_W + 14, 20, "B"),
        histogram_panel(rows, PANEL_W + AX["left"], AX["top"], pw, ph),
        legend(AX["left"], H - 8, set(counts), counts),
    ]
    head = ""
    if subtitle:
        head = (f'<text x="{W - AX["right"]}" y="19" text-anchor="end" font-size="12.5" '
                f'fill="#333">{_esc(subtitle)}</text>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
           f'font-family="Helvetica, Arial, sans-serif">'
           f'<rect width="100%" height="100%" fill="white"/>{head}{"".join(body)}</svg>')
    os.makedirs(os.path.dirname(os.path.abspath(outfile)) or ".", exist_ok=True)
    with open(outfile, "w") as fh:
        fh.write(svg)
    return outfile
