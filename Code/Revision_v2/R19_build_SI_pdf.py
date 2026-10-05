#!/usr/bin/env python
"""
R19_build_SI_pdf.py -- build the Supporting Information PDF for the v2 submission.

Replaces the old Matplotlib PdfPages builder, which was never committed to the
repository, so the v1 supplementary PDF could not be regenerated from source.
That is the root cause of the v1 defect where the shipped supplement came from a
superseded analysis and contradicted the main text.

Layout rules
------------
*   Captions are parsed out of Manuscript_PLOSONE_v2.md itself, so the supplement
    and the manuscript cannot drift apart: a caption that is edited in the
    manuscript is edited in the supplement on the next build.
*   Figures render on portrait A4; tables render on landscape A4.
*   Tables wider than MAX_COLS are split across pages by column blocks, with the
    first column (the row key) repeated on every block so each page stands alone.
*   Tables longer than MAX_FULL_ROWS are printed as a bounded preview with an
    explicit notice, because two of them are machine-scale tables (S7 is 20,426
    rows, S8 is 2,967). Silent truncation would misrepresent coverage, so the
    notice states the true row count and names the CSV that carries the rest.
*   A missing input renders a placeholder page rather than aborting the build,
    so the SI can be built while R14 is still running.

Run: python R19_build_SI_pdf.py
"""
import glob
import os
import re
import sys
import textwrap

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.stdout.reconfigure(encoding="utf-8")

ROOT = "D:/R_projects/revision_analysis"
TAB = f"{ROOT}/v2_output/tables"
FIG = f"{ROOT}/v2_output/figures"
MS = f"{ROOT}/Manuscript_PLOSONE_v2.md"
OUT = f"{ROOT}/v2_output/Supporting_Information_v2.pdf"

MAX_COLS = 12        # columns per page block before splitting
MAX_FULL_ROWS = 200  # tables longer than this are previewed, with a notice
PREVIEW_ROWS = 60

A4_PORTRAIT = (8.27, 11.69)
A4_LANDSCAPE = (11.69, 8.27)

# Row geometry. SCALE_Y multiplies every row's height, so the pagination
# calculation must include it: an earlier version omitted it, drew 1.6x taller
# rows than it had budgeted for, and silently dropped the rows that overflowed
# the bottom of the page.
SCALE_Y = 1.6
ROW_FACTOR = 1.7     # row height in font-size units before SCALE_Y
TABLE_TOP = 0.80     # axes top, figure fraction
TABLE_BOTTOM = 0.03  # axes bottom, figure fraction
CAPTION_RESERVE = 1.35  # inches held for the caption banner

# ----------------------------------------------------------------- caption parse

def load_captions():
    """Pull '**S3 Table.** ...' / '**S2 Fig.** ...' captions out of the manuscript."""
    txt = open(MS, encoding="utf-8").read()
    fig, tab = {}, {}
    for m in re.finditer(r"^\*\*(S\d+[a-z]?) (Table|Fig)\.\*\*\s*(.+?)(?=\n\n|\Z)",
                         txt, re.M | re.S):
        label, kind, body = m.group(1), m.group(2), m.group(3)
        body = " ".join(body.split())
        (tab if kind == "Table" else fig)[label] = clean(body)
    return fig, tab

def _script_run(s):
    """Rewrite superscript/subscript runs as ^n / _n for the PDF's Helvetica.

    Enumerating codepoints one at a time does not work here. The manuscript
    mixes blocks -- "I^2" uses Latin-1 U+00B2 while the exponents in "10^-57"
    use U+2075/U+2077, and "10^-15" combines U+207B with the Latin-1 U+00B9 --
    so a per-codepoint table silently converts only the digits it happens to
    list, leaving "10^-" followed by a stray raised glyph. Handling both blocks
    as runs fixes the whole family at once.
    """
    sup = {0x2070 + i: str(i) for i in range(10)}
    sup.update({0x00B9: "1", 0x00B2: "2", 0x00B3: "3"})
    sup[0x207B] = "-"          # without this, "10^-57" prints as "10⁻^57"
    sup[0x207A] = "+"
    sub = {0x2080 + i: str(i) for i in range(10)}
    sub.update({0x1D62: "i", 0x2C7C: "j"})   # Latin subscript i / j
    out, i, n = [], 0, len(s)
    while i < n:
        table, mark = (sup, "^") if ord(s[i]) in sup else \
                      (sub, "_") if ord(s[i]) in sub else (None, None)
        if table is None:
            out.append(s[i]); i += 1; continue
        run = []
        while i < n and ord(s[i]) in table:
            run.append(table[ord(s[i])]); i += 1
        out.append(mark + "".join(run))
    return "".join(out)

def clean(s):
    """Strip markdown emphasis and collapse <sub>/<sup> markup for print."""
    s = re.sub(r"<sub>(.*?)</sub>", r"\1", s)
    s = re.sub(r"<sup>(.*?)</sup>", r"^\1", s)
    s = s.replace("**", "").replace("*", "").replace("`", "")
    s = s.replace("ρ", "rho").replace("×", "x").replace("−", "-")
    return _script_run(s)

# ------------------------------------------------------------------- rendering

def draw_caption(fig, text, title=None):
    fig.text(0.06, 0.955, text, ha="left", va="top", fontsize=9.5,
             wrap=True, linespacing=1.45,
             bbox=dict(boxstyle="round,pad=0.55", fc="#f4f6f8", ec="#b8c2cc", lw=0.7))
    if title:
        fig.text(0.06, 0.985, title, ha="left", va="top", fontsize=12, weight="bold")

def fmt_val(v):
    """Format one cell for print. Numbers are rounded for legibility only --
    the CSVs keep full precision. Unrounded floats rendered at 15 significant
    digits, which forced every numeric column absurdly wide and squeezed the
    header row into overlapping ellipses."""
    if v is None:
        return ""
    if isinstance(v, float):
        if v != v:            # NaN
            return ""
        if v == int(v) and abs(v) < 1e15:
            return str(int(v))
        return f"{v:.4g}"
    s = str(v)
    return "" if s.lower() == "nan" else s

def wrap_header(name, width=13):
    """Wrap a column header onto multiple lines instead of truncating it.
    Truncation was silently colliding adjacent headers into unreadable runs."""
    parts = re.split(r"[._]", str(name))
    lines, cur = [], ""
    for p in parts:
        cand = f"{cur}.{p}" if cur else p
        if len(cand) <= width or not cur:
            cur = cand
        else:
            lines.append(cur); cur = p
    if cur:
        lines.append(cur)
    return "\n".join(lines[:4])

def wrap_cell(v, width=34):
    s = fmt_val(v)
    return textwrap.shorten(s, width=width, placeholder="...") if len(s) > width else s

def col_widths(sub, header_lines, total=1.0):
    """Proportional widths, so a 4-char column does not get the same space as a
    40-character cohort list."""
    w = []
    for j, c in enumerate(sub.columns):
        hdr = max((len(x) for x in header_lines[j].split("\n")), default=4)
        cell = max((len(fmt_val(v)) for v in sub[c].head(80)), default=4)
        w.append(max(hdr, min(cell, 30), 4))
    tot = float(sum(w))
    return [x / tot * total for x in w]

def _page_frame(label, caption, suf, note):
    fig = plt.figure(figsize=A4_LANDSCAPE)
    draw_caption(fig, caption + (("\n\n" + note) if note else ""),
                 title=f"{label} Table{suf}")
    ax = fig.add_axes([0.035, TABLE_BOTTOM, 0.93, TABLE_TOP - TABLE_BOTTOM])
    ax.axis("off")
    return fig, ax

def _place_table(ax, sub, hdr, fs, numcol):
    cells = [[wrap_cell(v) for v in row] for row in sub.itertuples(index=False)]
    tbl = ax.table(cellText=cells or [["(no rows)"]],
                   colLabels=hdr,
                   colWidths=col_widths(sub, hdr),
                   loc="upper left", cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fs)
    tbl.scale(1, SCALE_Y)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.35)
        if r == 0:
            cell.set_facecolor("#e8edf2")
            cell.set_text_props(weight="bold", va="center")
        else:
            if numcol[c]:
                cell.get_text().set_ha("right")
            if r % 2 == 0:
                cell.set_facecolor("#fafbfc")
    return tbl

def _fits(fig, ax, tbl):
    """Does the drawn table stay inside the axes? Row height depends on font
    metrics, wrapped header lines and the scale factor, so it is measured on a
    real draw rather than predicted -- the prediction is what clipped rows."""
    fig.canvas.draw()
    tb = tbl.get_window_extent(fig.canvas.get_renderer())
    ab = ax.get_window_extent()
    return tb.y0 >= ab.y0 - 0.5

def _layout_pages(sub, hdr, fs, numcol, caption, label, suf, note):
    """Greedily fit as many rows per page as actually render. Returns a list of
    (start, stop) row index pairs.

    'Rows fit' is monotone in the row count, so the largest fitting count is
    found by binary search (about 6 draws) rather than by counting down from 60.
    """
    spans, i, n = [], 0, len(sub)
    while i < n:
        lo, hi, best = 1, min(n - i, 60), 1
        while lo <= hi:
            mid = (lo + hi) // 2
            fig, ax = _page_frame(label, caption, suf, note)
            tbl = _place_table(ax, sub.iloc[i:i + mid], hdr, fs, numcol)
            ok = _fits(fig, ax, tbl)
            plt.close(fig)
            if ok:
                best = mid
                lo = mid + 1
            else:
                hi = mid - 1
        spans.append((i, i + best))
        i += best
    return spans

def render_table_pages(pdf, label, caption, df, src_note):
    """Landscape A4; split by columns, then paginate rows by measured fit."""
    cols = list(df.columns)
    blocks = []
    if len(cols) <= MAX_COLS:
        blocks.append(cols)
    else:
        key = cols[:1]
        rest = cols[1:]
        step = MAX_COLS - 1
        for i in range(0, len(rest), step):
            blocks.append(key + rest[i:i + step])

    truncated = len(df) > MAX_FULL_ROWS
    body = df.head(PREVIEW_ROWS) if truncated else df
    note = ""
    if truncated:
        note = (f"Showing the first {PREVIEW_ROWS} of {len(df):,} rows. "
                f"The complete table is provided as the machine-readable file "
                f"{os.path.basename(src_note)}.")

    for bi, block in enumerate(blocks):
        sub = body[block]
        hdr = [wrap_header(c) for c in block]
        fs = 7.0 if len(block) <= 9 else 6.2 if len(block) <= 16 else 5.6
        numcol = [pd.api.types.is_numeric_dtype(sub[c]) for c in block]

        col_suf = f"  (columns {bi+1} of {len(blocks)})" if len(blocks) > 1 else ""
        spans = _layout_pages(sub, hdr, fs, numcol, caption, label, col_suf, note)
        total = len(spans)
        for p, (a, b) in enumerate(spans, 1):
            suf = col_suf + (f"  (page {p} of {total})" if total > 1 else "")
            # the truncation notice belongs on the final page of the final block
            pg_note = note if (p == total and bi == len(blocks) - 1) else ""
            fig, ax = _page_frame(label, caption, suf, pg_note)
            _place_table(ax, sub.iloc[a:b], hdr, fs, numcol)
            pdf.savefig(fig)
            plt.close(fig)

def render_figure_page(pdf, label, caption, png):
    fig = plt.figure(figsize=A4_PORTRAIT)
    draw_caption(fig, caption, title=f"{label} Fig.")
    ax = fig.add_axes([0.05, 0.04, 0.90, 0.76])
    ax.axis("off")
    ax.imshow(plt.imread(png))
    pdf.savefig(fig)
    plt.close(fig)

def render_placeholder(pdf, label, kind, reason, caption=""):
    fig = plt.figure(figsize=A4_PORTRAIT)
    draw_caption(fig, caption, title=f"{label} {kind}")
    fig.text(0.5, 0.5, reason, ha="center", va="center", fontsize=11,
             color="#a03030", wrap=True,
             bbox=dict(boxstyle="round,pad=0.8", fc="#fdf2f2", ec="#d09090"))
    pdf.savefig(fig)
    plt.close(fig)

def manuscript_title():
    """The title exactly as the manuscript states it.

    Read rather than typed. The cover page used to carry its own copy of the v1
    title, and it was still printing that copy after the title was changed on
    2026-10-05 -- the supplement would have gone out under a title the
    manuscript no longer used. A title that is written down twice is a title
    that will disagree with itself, which is the same defect class as the v1
    supplement having no script behind it.
    """
    with open(MS, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("# "):
                return " ".join(re.sub(r"[*_`]", "", line[2:]).split())
    raise SystemExit(f"no H1 title line found in {MS}")


def render_cover(pdf, n_fig, n_tab):
    fig = plt.figure(figsize=A4_PORTRAIT)
    fig.text(0.5, 0.72, "Supporting Information", ha="center", fontsize=24, weight="bold")
    title_lines = textwrap.wrap(manuscript_title(), width=62)
    # Block-anchor so the title sits the same distance below the heading whatever
    # number of lines it wraps to; the old fixed y-positions only worked for the
    # three lines of the v1 title.
    y = 0.665
    for ln in title_lines:
        fig.text(0.5, y, ln, ha="center", fontsize=11)
        y -= 0.025
    fig.text(0.5, y - 0.020,
             f"{n_fig} supplementary figures  ·  {n_tab} supplementary tables",
             ha="center", fontsize=11, color="#40506a")
    fig.text(0.5, 0.44,
             "Every figure and table in this supplement is generated by the frozen\n"
             "pipeline under Code/Revision_v2/ from public data. The configuration\n"
             "hash recorded in S1 Table identifies the exact settings used.",
             ha="center", fontsize=9.5, color="#40506a", linespacing=1.6)
    fig.text(0.5, 0.30,
             "Tables larger than 200 rows are printed here as a bounded preview;\n"
             "the complete table is supplied as a machine-readable CSV alongside\n"
             "this document. Each such preview states its true row count.",
             ha="center", fontsize=9, color="#7a3020", linespacing=1.6,
             bbox=dict(boxstyle="round,pad=0.6", fc="#fdf6f2", ec="#e0c0b0"))
    fig.text(0.5, 0.10, "Generated by Code/Revision_v2/R19_build_SI_pdf.py",
             ha="center", fontsize=8, color="#8899aa")
    pdf.savefig(fig)
    plt.close(fig)

# ------------------------------------------------------------------------ main

def main():
    caps_fig, caps_tab = load_captions()
    fig_labels = sorted(caps_fig, key=lambda s: int(re.sub(r"\D", "", s)))
    tab_labels = sorted(caps_tab, key=lambda s: (int(re.sub(r"\D", "", s)), s))

    print(f"captions parsed: {len(fig_labels)} figures, {len(tab_labels)} tables")

    missing = []
    with PdfPages(OUT) as pdf:
        render_cover(pdf, len(fig_labels), len(tab_labels))

        for lab in fig_labels:
            png = f"{FIG}/Fig{lab}_*.png"
            hits = [p for p in glob.glob(f"{FIG}/Fig{lab}_*.png")]
            if hits:
                render_figure_page(pdf, lab, caps_fig[lab], hits[0])
                print(f"  fig {lab}: {os.path.basename(hits[0])}")
            else:
                missing.append(f"{lab} Fig")
                render_placeholder(pdf, lab, "Fig.",
                                   f"Figure file not yet generated.\nExpected {png}",
                                   caps_fig[lab])
                print(f"  fig {lab}: MISSING -> placeholder")

        for lab in tab_labels:
            hits = glob.glob(f"{TAB}/Table_{lab}_*.csv")
            if hits:
                src = hits[0]
                df = pd.read_csv(src)
                render_table_pages(pdf, lab, caps_tab[lab], df, src)
                print(f"  tab {lab}: {os.path.basename(src)} ({len(df)}r x {len(df.columns)}c)")
            else:
                missing.append(f"{lab} Table")
                render_placeholder(pdf, lab, "Table",
                                   f"Table file not yet generated.\nExpected Table_S{lab}_*.csv",
                                   caps_tab[lab])
                print(f"  tab {lab}: MISSING -> placeholder")

    size_mb = os.path.getsize(OUT) / 1e6
    print(f"\nwrote {OUT}  ({size_mb:.1f} MB)")
    if missing:
        print(f"placeholders rendered for: {', '.join(missing)}")
        return 0
    print("complete: no placeholder pages.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
