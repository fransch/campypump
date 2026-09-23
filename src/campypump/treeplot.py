"""Draw an annotated cmeB allele tree as SVG (parsing by Biopython, drawing by hand).

Branches are coloured by what lies below them: a branch whose descendants are all one cmeB type
takes that colour, mixed branches stay grey. Clades that are all one species are shaded, so
interspersed alleles are visible at a glance. Support is drawn only where it is high enough to
mean something.
"""
from __future__ import annotations

import csv
import io
import math as math0
import os
import re

from .figures import COLORS, _esc
from .palettes import GREY_PALE, SOURCES as SOURCE_COLORS, SOURCE_OTHER, SPECIES as SPECIES_COLORS, SPECIES_FILL as _FILL

from .palettes import TYPE_SHORT as TYPE_COLORS
# fill + fill-opacity, not 8-digit hex: some SVG renderers ignore the alpha channel
SPECIES_FILL = _FILL
LEGEND_NAMES = {"RE": "RE-cmeB", "WT": "wildtype cmeB", "MOS": "mosaic", "REF": "reference",
                "OUT": "outgroup", "Cc": "C. coli", "Cj": "C. jejuni"}


def tip_class(name: str) -> str:
    if name.startswith("OUT"):
        return "OUT"
    if name.startswith("REF"):
        return "REF"
    for tag in ("RE", "WT", "MOS"):
        if f"_{tag}_" in name:
            return tag
    return "mixed"


def tip_species(name: str) -> str:
    if name.startswith("OUT"):
        return "OUT"
    if name.startswith("REF"):
        return "Cc" if "_Cc_" in name else "Cj"
    return name.split("_")[0]


def pretty_tip(name: str) -> str:
    """Cc_a003_RE_n63 -> 'Cc a003 RE (n=63)'; REF_cmeB_Cj_RE_PA001_primary -> 'PA001 (Cj RE)'."""
    if name.startswith("OUT"):
        return name.replace("OUT_", "").replace("_", " ")
    if name.startswith("REF"):
        parts = name.replace("REF_cmeB_", "").replace("_panel", "").replace("_primary", "").split("_")
        sp, typ, strain = parts[0], parts[1], "_".join(parts[2:])
        return f"{strain.replace('_', '-')} ({sp} {typ})"
    m = re.match(r"^(C[cj])_a(\d+)_([A-Z]+)_n(\d+)$", name)
    if m:
        return f"{m[1]} a{m[2]} {m[3]} (n={m[4]})"
    return name


def _support(clade) -> tuple[float | None, float | None]:
    """IQ-TREE writes 'SH-aLRT/UFboot' as the node label."""
    raw = clade.confidence if clade.confidence is not None else clade.name
    if raw is None:
        return None, None
    m = re.match(r"^([\d.]+)/([\d.]+)$", str(raw).strip())
    if m:
        return float(m[1]), float(m[2])
    try:
        return None, float(raw)
    except (TypeError, ValueError):
        return None, None


def _n_from_name(name: str) -> int:
    m = re.search(r"_n(\d+)$", name)
    return int(m[1]) if m else 0


def tree_figure(newick_path: str, outfile: str, outgroup: str | None = None,
                label: str = "notable", min_support: float = 90.0,
                width: int = 1000, tip_height: float = 11.0, subtree: str | None = None,
                title: str = "", colour_by: str = "auto", drop_outgroup: bool = False,
                cladogram: bool = False, annotation: str | None = None,
                annotate_cols: tuple[str, ...] = ("top_STs",), size_by_genomes: bool = True,
                source_bars: bool = True) -> str:
    """Render a Newick tree.

    label: 'notable' (references, outgroup, mosaic and every non-wildtype allele), 'all', or 'none'.
    drop_outgroup: remove the outgroup tip after rooting, so the ingroup fills the width.
    annotation: the table from `campypump alleles` (`<prefix>.tsv`), keyed by allele id; the columns
      named in annotate_cols are appended to each tip label (default: the dominant STs).
    subtree: 'RE' draws only the clade containing every RE allele.
    """
    from Bio import Phylo

    extra: dict[str, str] = {}
    n_genomes: dict[str, int] = {}
    sources: dict[str, list] = {}
    blocks: dict[str, str] = {}
    if annotation:
        with open(annotation, newline="") as fh:
            first = fh.readline()
            fh.seek(0)
            delim = "\t" if first.count("\t") >= first.count(",") else ","
            for row in csv.DictReader(fh, delimiter=delim):
                key = row.get("allele_id") or row.get("ID") or ""
                bits = []
                for col in annotate_cols:
                    v = (row.get(col) or "").strip()
                    if not v:
                        continue
                    if col == "top_STs":
                        # "12298(59); 13856(3); 14512(1)" -> "ST 12298 (59 of 63)"
                        top = v.split(";")[0].strip()
                        m = re.match(r"^(\S+)\((\d+)\)$", top)
                        total = row.get("genomes", "")
                        bits.append(f"ST {m[1]} ({m[2]}/{total})" if m and total else f"ST {top}")
                    else:
                        bits.append(v)
                if key and bits:
                    extra[key] = "  ·  " + "  ·  ".join(bits)
                if key:
                    try:
                        n_genomes[key] = int(row.get("genomes") or 0)
                    except ValueError:
                        pass
                    blocks[key] = (row.get("block_carrier") or "").strip()
                    parsed = []
                    for item in (row.get("top_sources") or "").split(";"):
                        m2 = re.match(r"^\s*(.+?)\((\d+)\)\s*$", item)
                        if m2:
                            parsed.append((m2[1].strip(), int(m2[2])))
                    sources[key] = parsed

    tree = Phylo.read(io.StringIO(open(newick_path).read()), "newick")
    if outgroup:
        try:
            tree.root_with_outgroup(outgroup)
        except ValueError:
            pass
    if outgroup and drop_outgroup:
        # rooting done: the outgroup's long branch would otherwise compress the whole ingroup
        node = next((t for t in tree.get_terminals() if t.name == outgroup), None)
        if node is not None:
            tree.prune(node)
    tree.ladderize()
    root = tree.root
    if subtree == "RE":
        # the mosaic sits on a long branch at the base and would compress the rest
        re_tips = [t.name for t in tree.get_terminals() if tip_class(t.name) == "RE"]
        if re_tips:
            root = tree.common_ancestor(re_tips)

    # inside the RE clade every tip is RE, so colour by species there instead
    if colour_by == "auto":
        colour_by = "species" if subtree == "RE" else "type"
    palette = SPECIES_COLORS if colour_by == "species" else TYPE_COLORS
    key_of = tip_species if colour_by == "species" else tip_class

    tips = [t for t in root.get_terminals()]
    n = len(tips)
    y_of, x_of = {}, {}
    for i, t in enumerate(tips):
        y_of[id(t)] = i

    def depth(clade, acc=0.0):
        x_of[id(clade)] = acc
        for ch in clade.clades:
            depth(ch, acc + (1.0 if cladogram else (ch.branch_length or 0.0)))
    depth(root)
    if cladogram:      # push every tip to the right edge
        deepest = max(x_of[id(t)] for t in root.get_terminals())
        for t in root.get_terminals():
            x_of[id(t)] = deepest

    def y_pos(clade):
        if clade.is_terminal():
            return y_of[id(clade)]
        ys = [y_pos(c) for c in clade.clades]
        y_of[id(clade)] = (min(ys) + max(ys)) / 2
        return y_of[id(clade)]
    y_pos(root)

    def classify(clade):
        kinds = {key_of(t.name) for t in clade.get_terminals()}
        kinds.discard("REF")
        kinds.discard("OUT")
        return kinds.pop() if len(kinds) == 1 else "mixed"

    def species_of(clade):
        sp = {tip_species(t.name) for t in clade.get_terminals()}
        sp.discard("OUT")
        return sp.pop() if len(sp) == 1 else None

    max_x = max(x_of.values()) or 1.0
    left, right, top, bottom = 18, (400 if annotation else 250), 46 if title else 26, \
        (104 if (annotation and source_bars) else 72)
    plot_w = width - left - right
    height = int(top + bottom + n * tip_height)
    sx = lambda v: left + plot_w * v / max_x
    sy = lambda v: top + (v + 0.5) * tip_height

    block_ring = TYPE_COLORS["MOS"]      # named here: nested quotes inside an f-string need 3.12
    parts = [f'<rect width="100%" height="100%" fill="white"/>']
    if title:
        parts.append(f'<text x="{width/2}" y="22" text-anchor="middle" font-size="13" '
                     f'font-weight="bold" fill="#111">{_esc(title)}</text>')

    # species shading behind clades that are all one species and worth marking
    shaded = []
    def shade(clade):
        sp = species_of(clade)
        tps = clade.get_terminals()
        if sp in SPECIES_FILL and len(tps) >= 6:
            ys = [y_of[id(t)] for t in tps]
            shaded.append((sp, min(ys), max(ys)))
            return
        for ch in clade.clades:
            shade(ch)
    shade(root)
    for sp, y0, y1 in shaded:
        parts.append(f'<rect x="{left}" y="{sy(y0)-tip_height/2:.1f}" width="{plot_w}" '
                     f'height="{(y1-y0+1)*tip_height:.1f}" fill="{SPECIES_FILL[sp]}" '
                     f'fill-opacity="0.07"/>')

    # branches
    def draw(clade):
        colour = palette.get(classify(clade), palette["mixed"])
        x, y = sx(x_of[id(clade)]), sy(y_of[id(clade)])
        if clade.clades:
            ys = [sy(y_of[id(c)]) for c in clade.clades]
            parts.append(f'<line x1="{x:.1f}" y1="{min(ys):.1f}" x2="{x:.1f}" y2="{max(ys):.1f}" '
                         f'stroke="{colour}" stroke-width="1.4"/>')
            for c in clade.clades:
                cx = sx(x_of[id(c)])
                cc = palette.get(classify(c), palette["mixed"])
                parts.append(f'<line x1="{x:.1f}" y1="{sy(y_of[id(c)]):.1f}" x2="{cx:.1f}" '
                             f'y2="{sy(y_of[id(c)]):.1f}" stroke="{cc}" stroke-width="1.4"/>')
                draw(c)
            alrt, boot = _support(clade)
            if boot is not None and boot >= min_support and len(clade.get_terminals()) >= 3:
                txt = f"{alrt:g}/{boot:g}" if alrt is not None else f"{boot:g}"
                parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.1" fill="#444"/>')
                parts.append(f'<text x="{x-4:.1f}" y="{y-3:.1f}" text-anchor="end" font-size="7.5" '
                             f'fill="#666">{txt}</text>')
    draw(root)

    # tips
    present = []
    for t in tips:
        k, x, y = key_of(t.name), sx(x_of[id(t)]), sy(y_of[id(t)])
        if k not in present:
            present.append(k)
        colour = palette.get(k, palette["mixed"])
        n_gen = n_genomes.get(t.name) or _n_from_name(t.name)
        r = 2.4
        if size_by_genomes and n_gen:
            r = 2.0 + 1.5 * math0.log10(n_gen)          # 1 genome -> 2.0, 300 -> ~5.7
        cls = tip_class(t.name)
        show = (label == "all") or (label == "notable" and cls in ("RE", "MOS", "REF", "OUT"))
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="{colour}"/>')
        if blocks.get(t.name) == "yes":
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r+2.4:.1f}" fill="none" '
                         f'stroke="{block_ring}" stroke-width="1.5"/>')
        if show:
            weight = ' font-weight="bold"' if cls in ("REF", "OUT") else ""
            txt = pretty_tip(t.name)
            parts.append(f'<text x="{x+5:.1f}" y="{y+3:.1f}" font-size="8.5" fill="#222"{weight}>'
                         f'{_esc(txt)}<tspan fill="#777" font-weight="normal">'
                         f'{_esc(extra.get(t.name, ""))}</tspan></text>')

    import math
    if cladogram:
        y = height - bottom + 26
        parts.append(f'<text x="{left}" y="{y+2}" font-size="10" fill="#333">'
                     f'cladogram: branch lengths not to scale</text>')
    # source composition: a small stacked bar per allele, in its own column
    bar_x = width - right + 196
    bar_w = 46
    drew_sources = False
    if source_bars and sources:
        for t in tips:
            comp = sources.get(t.name) or []
            if not comp:
                continue
            drew_sources = True
            total = sum(c for _, c in comp) or 1
            y0 = sy(y_of[id(t)]) - 3.0
            xx = bar_x
            for src, c in comp:
                col = next((v for k, v in SOURCE_COLORS if src.lower().startswith(k)), SOURCE_OTHER)
                w = bar_w * c / total
                parts.append(f'<rect x="{xx:.1f}" y="{y0:.1f}" width="{max(w,0.6):.1f}" height="6" '
                             f'fill="{col}"/>')
                xx += w
    if drew_sources:
        parts.append(f'<text x="{bar_x:.1f}" y="{top-6}" font-size="8.5" fill="#555">'
                     f'top sources</text>')

    # scale bar: a round number about a fifth of the tree width
    raw = max_x / 5 if not cladogram else None
    if not cladogram:
        step = 10 ** math.floor(math.log10(raw))
        unit = min([m * step for m in (1, 2, 5, 10)], key=lambda v: abs(v - raw))
        y = height - bottom + 26
        parts.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_w*unit/max_x:.1f}" y2="{y}" '
                     f'stroke="#333" stroke-width="1.4"/>')
        parts.append(f'<text x="{left + plot_w*unit/max_x/2:.1f}" y="{y+13}" text-anchor="middle" '
                     f'font-size="10" fill="#333">{unit:g} substitutions per site</text>')
        lx = left + plot_w * unit / max_x + 40
    else:
        lx = left + 210
    order = ["Cc", "Cj", "RE", "WT", "MOS", "REF", "OUT"]
    for key in [k for k in order if k in present]:
        name = LEGEND_NAMES.get(key, key)
        parts.append(f'<line x1="{lx}" y1="{y-3}" x2="{lx+16}" y2="{y-3}" '
                     f'stroke="{palette.get(key, palette["mixed"])}" stroke-width="3"/>')
        parts.append(f'<text x="{lx+21}" y="{y}" font-size="10" fill="#222">{_esc(name)}</text>')
        lx += 30 + 6.2 * len(name)
    if drew_sources:
        lx2, y2 = left, y + 32
        parts.append(f'<text x="{lx2}" y="{y2}" font-size="9" fill="#555">sources:</text>')
        lx2 += 52
        for name, col in SOURCE_COLORS[:6] + [("other", SOURCE_OTHER)]:
            parts.append(f'<rect x="{lx2}" y="{y2-7}" width="10" height="7" fill="{col}"/>')
            parts.append(f'<text x="{lx2+13}" y="{y2}" font-size="9" fill="#333">{_esc(name)}</text>')
            lx2 += 18 + 5.4 * len(name)
    note = (f'branch colour: {"species" if colour_by == "species" else "cmeB type"} of everything '
            f'below it, grey where mixed; node labels: SH-aLRT/bootstrap, shown at or above '
            f'{min_support:g}')
    note += ('; tip size: genomes carrying the allele (log); gold ring: allele carries the '
             'wildtype-like block') if size_by_genomes else ''
    parts.append(f'<text x="{left}" y="{y + (50 if drew_sources else 30)}" font-size="9" '
                 f'fill="#777">{_esc(note)}</text>')

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}" font-family="Helvetica, Arial, sans-serif">'
           + "".join(parts) + "</svg>")
    os.makedirs(os.path.dirname(os.path.abspath(outfile)) or ".", exist_ok=True)
    with open(outfile, "w") as fh:
        fh.write(svg)
    return outfile
