"""Segment plots: identity to RE-cmeB vs wildtype cmeB along the gene (plain SVG, no dependencies)."""
from __future__ import annotations

import os

W, H = 820, 340
LEFT, RIGHT, TOP, BOTTOM = 70, 20, 70, 60
Y_MIN, Y_MAX = 60.0, 100.0
from .palettes import SEGMENT_RE as RE_COLOR, SEGMENT_WT as WT_COLOR


def _x(pos: float, length: int) -> float:
    return LEFT + (W - LEFT - RIGHT) * pos / max(length, 1)


def _y(pct: float) -> float:
    pct = min(max(pct, Y_MIN), Y_MAX)
    return TOP + (H - TOP - BOTTOM) * (Y_MAX - pct) / (Y_MAX - Y_MIN)


def write_segment_plot(outdir: str, row: dict) -> str:
    """Write <sample>_cmeB_segments.svg and .tsv; return the SVG path."""
    os.makedirs(outdir, exist_ok=True)
    sample = row["sample"]
    windows = row.get("_windows") or []
    length = int(row.get("cmeB_length_bp") or 3123)
    base = os.path.join(outdir, f"{sample}_cmeB_segments")

    with open(base + ".tsv", "w") as fh:
        fh.write("position_bp\tidentity_RE_PA001_pct\tidentity_best_wildtype_pct\n")
        for pos, re_id, wt_id in windows:
            fh.write(f"{pos}\t{re_id:.1f}\t{wt_id:.1f}\n")

    def line(idx, color):
        pts = " ".join(f"{_x(w[0], length):.1f},{_y(w[idx]):.1f}" for w in windows)
        return f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2.2"/>'

    grid = []
    for pct in range(int(Y_MIN), int(Y_MAX) + 1, 10):
        y = _y(pct)
        grid.append(f'<line x1="{LEFT}" y1="{y:.1f}" x2="{W - RIGHT}" y2="{y:.1f}" stroke="#ddd"/>'
                    f'<text x="{LEFT - 8}" y="{y + 4:.1f}" text-anchor="end" font-size="12">{pct}%</text>')
    for bp in range(0, length + 1, 500):
        x = _x(bp, length)
        grid.append(f'<line x1="{x:.1f}" y1="{H - BOTTOM}" x2="{x:.1f}" y2="{H - BOTTOM + 5}" stroke="#333"/>'
                    f'<text x="{x:.1f}" y="{H - BOTTOM + 20}" text-anchor="middle" font-size="12">{bp}</text>')

    title = (f"{sample}: cmeB {row.get('cmeB_type_final', row.get('cmeB_type', ''))} "
             f"({row.get('cmeB_identity_RE_PA001_pct', '')}% to PA001 RE, "
             f"{row.get('cmeB_identity_WT_11168_pct', '')}% to 11168)")
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="Arial, sans-serif">
<rect width="100%" height="100%" fill="white"/>
<text x="{W / 2}" y="24" text-anchor="middle" font-size="14" font-weight="bold">{title}</text>
{"".join(grid)}
<line x1="{LEFT}" y1="{H - BOTTOM}" x2="{W - RIGHT}" y2="{H - BOTTOM}" stroke="#333"/>
<line x1="{LEFT}" y1="{TOP}" x2="{LEFT}" y2="{H - BOTTOM}" stroke="#333"/>
{line(1, RE_COLOR)}
{line(2, WT_COLOR)}
<text x="{W / 2}" y="{H - 15}" text-anchor="middle" font-size="12">Position along cmeB (bp); 300-bp windows, step 50 bp</text>
<text transform="translate(18,{(TOP + H - BOTTOM) / 2}) rotate(-90)" text-anchor="middle" font-size="12">% identity</text>
<line x1="{LEFT + 150}" y1="46" x2="{LEFT + 180}" y2="46" stroke="{RE_COLOR}" stroke-width="2.5"/>
<text x="{LEFT + 187}" y="50" font-size="12">RE-cmeB (PA001)</text>
<line x1="{LEFT + 330}" y1="46" x2="{LEFT + 360}" y2="46" stroke="{WT_COLOR}" stroke-width="2.5"/>
<text x="{LEFT + 367}" y="50" font-size="12">best wildtype reference (11168 or panel)</text>
</svg>
'''
    with open(base + ".svg", "w") as fh:
        fh.write(svg)
    return base + ".svg"
