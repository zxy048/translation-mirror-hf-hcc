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

# ------------------------------------------------------------------ page budget
#
# The supplement grows with its inputs unless something stops it, and the
# something used to be a row count: print in full up to 200 rows, preview at 60
# past that. Two tables showed why a row count is the wrong unit.
#
#   * S37 holds 163 rows, so it sat *under* the 200-row threshold, printed every
#     row, and cost 24 pages.
#   * S19 holds 45 rows -- one page of content -- and cost 12, because its 37
#     columns split into four blocks and each block reprinted the same 45 rows.
#
# The driver is rows x column blocks, so the budget is stated in pages and the
# row count is derived from it: each column block is given MAX_PAGES_PER_BLOCK
# pages and carries however many rows actually fit on them. A table small enough
# to fit is still printed whole and carries no notice, so the rule degrades to
# the previous behaviour for the tables that were never the problem.
# The budget is per TABLE, not per column block. Counting per block made a
# two-block table four pages before anything was looked at, because each block
# reprints the same rows through different columns and each was given the full
# budget of its own.
MAX_PAGES_PER_TABLE = 2

def max_cols_for(n_rows):
    """Columns per block, by row count.

    Splitting on a flat count wastes pages on short tables: S20 panel a is two
    rows across 37 numeric columns and cost four pages, because 12 columns per
    block split it four ways. A block is limited by width, and a table with two
    rows has vertical space to spend on narrow columns.
    """
    if n_rows > 20:
        return 12
    if n_rows > 6:
        return 20
    return 40

# Tables the manuscript cites by value are printed as a structural preview when
# they do not fit, and the preview narrows the column set as well as the row set.
# S37 carries 25 columns and S19 37, and printing the other eighteen is most of
# what the pages were spent on; the columns kept are the ones the cited claim is
# made of. Narrowing columns is only done for tables listed here, so the
# reduction is a decision recorded per table rather than a rule applied blind.
PREVIEW_COLS = {
    "S37": ["cohort", "set", "tier", "n", "g_unadjusted", "g_adjusted",
            "attenuation", "survives"],
    "S19": ["reference", "test", "module", "moduleSize", "Zsummary.pres",
            "medianRank.pres", "preservation_class"],
    "S17": ["set", "g_HCC", "g_HF", "mirror", "tier", "cluster",
            "hcc_cohort", "hf_cohort"],
}

# Rows the manuscript cites, ordered to the front of a preview so that a reader
# checking a number in the text finds the row it came from on the printed page.
# Each predicate is applied to the table itself rather than to a hard-coded row
# list, so it cannot go stale when the table is regenerated:
#   S37, S17, S21 -- the manuscript reports counts over the Tier A sets
#                    ("GSE57338 gives 10 of 15 Tier A sets", "10 of 15 Tier A
#                    sets are mirror-perturbed", median I2 across the Tier A
#                    sets), which is a statement about every Tier A row.
#   S19           -- the two designated modules, which are the only rows whose
#                    Zsummary the manuscript and the response letter quote.
#   S30           -- the seven transcription factors the text names by value in
#                    TCGA-LIHC (ATF4, MYC, DDIT3, EIF4EBP1, EIF2AK3, ERN1,
#                    MTOR), on the TCS score. Narrowed to that score because the
#                    quoted coefficients are all partial correlations against
#                    TCS; widening it to every score triples the pinned set and
#                    pushes the table back over its budget for no extra citation.
#   S39           -- the continuous per-IQR Cox rows, which is the model every
#                    hazard ratio in the text and in the response letter is
#                    quoted from.
#
# Tables cited only through a *derived* statistic get no pin and are not meant
# to: S14's "mean |g| is 1.44" and S25's "55 of 103 Tier A set-cohort pairs" are
# computed over every row, so no bounded preview can make them checkable by eye.
# They are checkable against the CSV, which the cover names as the authoritative
# copy, and that is the distinction the preview notice exists to draw.
# Written as column -> allowed values rather than as a predicate, so that a pin
# can be checked against the frame it is applied to. A merged item publishes
# several CSVs under one number, and a pin written for the parent table does not
# describe all of them: S39's pin names a `score` column that S39 panel b -- the
# annotation-coverage table -- does not have. As a predicate that was a KeyError
# that stopped the build; as data it is a pin that simply does not apply.
PIN_ROWS = {
    "S37": {"tier": ["A"]},
    "S17": {"tier": ["A"]},
    "S21": {"tier": ["A"]},
    "S19": {"module": ["magenta", "purple"]},
    "S30": {"tf": ["ATF4", "MYC", "DDIT3", "EIF4EBP1", "EIF2AK3", "ERN1", "MTOR"],
            "cohort": ["TCGA_LIHC"], "score": ["TCS"]},
    "S39": {"score": ["TCS"], "model": ["continuous_per_IQR"]},
}


def pin_mask(df, label):
    """Rows the pin for `label` selects, or None where no pin applies.

    None rather than an all-False mask, so a caller cannot mistake "this table
    has no cited rows to protect" for "this table's cited rows are empty".
    """
    spec = PIN_ROWS.get(label)
    if not spec:
        return None
    if any(col not in df.columns for col in spec):
        return None
    mask = pd.Series(True, index=df.index)
    for col, allowed in spec.items():
        mask &= df[col].isin(allowed)
    return mask

# PLOS requires supporting-information item numbers to be whole numbers
# ("Combine separate parts ... into one file, or rename with whole numbers"), so
# the five tables that carried a letter suffix in the v2 draft are published as
# panels of their parent table rather than as separate items. The CSVs keep the
# filenames they were built under -- they are repository artefacts, not
# submitted item names -- so the label-to-source mapping has to be recorded
# here. A label absent from this map is its own single source and renders
# exactly as before.
PANELS = {
    "S20": ["S20", "S20b"],
    "S22": ["S22", "S22b", "S22c", "S22d"],
    "S39": ["S39", "S39b"],
}
PANEL_LETTER = "abcdefgh"

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

MAX_ROWS_PER_PAGE = 60   # ceiling for the per-page fit search, both here and in
                         # rows_that_fit: if the two disagreed, the budget would
                         # be computed against one row height and paginated at
                         # another, and the pages would exceed the budget.


def _layout_pages(sub, hdr, fs, numcol, caption, label, suf, note):
    """Greedily fit as many rows per page as actually render. Returns a list of
    (start, stop) row index pairs.

    'Rows fit' is monotone in the row count, so the largest fitting count is
    found by binary search (about 6 draws) rather than by counting down from 60.
    """
    spans, i, n = [], 0, len(sub)
    while i < n:
        lo, hi, best = 1, min(n - i, MAX_ROWS_PER_PAGE), 1
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

def table_panels(label):
    """Every CSV published under `label`, as (panel_letter, csv_path) pairs.

    Returns [(None, path)] for an ordinary single-source table, and
    [('a', path), ('b', path), ...] for a merged one. A panel whose CSV is
    missing comes back as (letter, None) so the caller can render a placeholder
    in its place rather than dropping the panel silently -- the failure this
    whole file exists to avoid.
    """
    srcs = PANELS.get(label, [label])
    merged = len(srcs) > 1
    out = []
    for i, s in enumerate(srcs):
        letter = PANEL_LETTER[i] if merged else None
        hits = glob.glob(f"{TAB}/Table_{s}_*.csv")
        out.append((letter, hits[0] if hits else None))
    return out


def _block_meta(body, block):
    """Header lines, font size and numeric-column flags for one column block."""
    hdr = [wrap_header(c) for c in block]
    fs = 7.0 if len(block) <= 9 else 6.2 if len(block) <= 16 else 5.6
    numcol = [pd.api.types.is_numeric_dtype(body[c]) for c in block]
    return hdr, fs, numcol


def rows_that_fit(sub, hdr, fs, numcol, caption, label, suf, note):
    """The largest number of leading rows of `sub` that renders on ONE page.

    Binary search on the drawn extent rather than a predicted count: row height
    depends on font metrics, wrapped headers and SCALE_Y, and predicting it is
    what clipped rows when this file was first written.
    """
    if len(sub) == 0:
        return 0
    lo, hi, best = 1, min(len(sub), MAX_ROWS_PER_PAGE), 0
    while lo <= hi:
        mid = (lo + hi) // 2
        fig, ax = _page_frame(label, caption, suf, note)
        tbl = _place_table(ax, sub.iloc[:mid], hdr, fs, numcol)
        ok = _fits(fig, ax, tbl)
        plt.close(fig)
        if ok:
            best, lo = mid, mid + 1
        else:
            hi = mid - 1
    return best


def order_rows(df, label):
    """Rows the manuscript cites first, everything else in file order.

    A preview is a selection, so which rows it takes is a real choice and not a
    detail. The first N would cut whichever rows happen to sort first, and the
    cited rows are not those: S37's are the Tier A rows of three cohorts, which
    sit behind five rows of another tier within every cohort block. Putting them
    first is what makes the truncation safe to make at all.
    """
    mask = pin_mask(df, label)
    if mask is None or not mask.any():
        return df
    return pd.concat([df[mask], df[~mask]])


def table_notice(csv_name, n_rows, n_cols, shown_rows, shown_cols, ordered):
    """One sentence naming what was reduced and where the whole table is."""
    if shown_rows == n_rows and shown_cols == n_cols:
        return f"Complete table: {csv_name} ({n_rows} rows x {n_cols} columns)."
    parts = []
    if shown_cols < n_cols:
        parts.append(f"{shown_cols} of {n_cols} columns")
    if shown_rows < n_rows:
        parts.append(f"{shown_rows} of {n_rows} rows")
    tail = " Rows cited in the manuscript are printed first." if ordered else ""
    return (f"Printed from {csv_name} ({n_rows} rows x {n_cols} columns): "
            + " and ".join(parts)
            + f".{tail} The CSV is the complete and authoritative version.")


def table_plan(label, caption, df, csv_name, panel=None):
    """Decide what is printed for one table, and say so.

    Returns (blocks, ordered_rows_shown, note, panel_suffix). Kept in one place
    and called from both the builder and R20: a verifier that recomputed the row
    count for itself would be checking its own arithmetic rather than the drawn
    page, and the two would agree exactly until the day one of them was edited.
    """
    n_rows, n_cols = len(df), len(df.columns)

    cols = PREVIEW_COLS.get(label, list(df.columns))
    absent = [c for c in cols if c not in df.columns]
    if absent:
        # A preview column set naming a column the table no longer has would
        # drop it without a word, which is the failure class this budget exists
        # to prevent -- so it aborts the build instead.
        raise SystemExit(f"{label}: PREVIEW_COLS names columns absent from "
                         f"{csv_name}: {absent}")

    ordered = order_rows(df, label)
    pinned = pin_mask(ordered, label)
    was_ordered = pinned is not None and bool(pinned.any())
    body = ordered[cols]

    per_block = max_cols_for(n_rows)
    if len(cols) <= per_block:
        blocks = [cols]
    else:
        key, rest, step = cols[:1], cols[1:], per_block - 1
        blocks = [key + rest[i:i + step] for i in range(0, len(rest), step)]

    panel_suf = f"  (panel {panel})" if panel else ""
    n_pin = int(pinned.sum()) if pinned is not None else 0

    def plan(note):
        """How many rows this table prints, given the caption it will carry."""
        cap = [rows_that_fit(body[b], *_block_meta(body, b),
                             caption, label, panel_suf, note) for b in blocks]
        if not cap or min(cap) == 0:
            return 0
        # Pages are summed over blocks: the blocks reprint the same rows through
        # different columns, so a two-block table printed one page each is two
        # pages of the document, and that is what the budget has to count.
        n = min(cap) * MAX_PAGES_PER_TABLE
        while n > 0 and sum(-(-n // c) for c in cap) > MAX_PAGES_PER_TABLE:
            n -= 1
        # The cited rows are a floor and not an option. A preview that hides the
        # row a manuscript number came from is worse than a long supplement, so
        # where the floor needs more pages than the budget, the floor wins and
        # the notice reports the larger count.
        return min(n_rows, max(n, n_pin))

    shown = plan("")
    note = table_notice(csv_name, n_rows, n_cols, shown, len(cols), was_ordered)
    if shown < n_rows or len(cols) < n_cols:
        # The notice is part of the caption and takes vertical space, so re-fit
        # with it present. It can only reduce the row count, so one pass is
        # enough; the second notice then reports the reduced count.
        shown = plan(note)
        note = table_notice(csv_name, n_rows, n_cols, shown, len(cols), was_ordered)

    return blocks, body.iloc[:shown], note, panel_suf


def render_table_pages(pdf, label, caption, df, src_note, panel=None):
    """Landscape A4; narrow the columns where the label has a preview set, split
    into column blocks, and give each block at most MAX_PAGES_PER_BLOCK pages."""
    blocks, shown_body, note, panel_suf = table_plan(
        label, caption, df, os.path.basename(src_note), panel)

    for bi, block in enumerate(blocks):
        sub = shown_body[block]
        hdr, fs, numcol = _block_meta(shown_body, block)
        col_suf = f"  (columns {bi+1} of {len(blocks)})" if len(blocks) > 1 else ""
        base_suf = panel_suf + col_suf
        spans = _layout_pages(sub, hdr, fs, numcol, caption, label, base_suf, note)
        total = len(spans)
        for p, (a, b) in enumerate(spans, 1):
            suf = base_suf + (f"  (page {p} of {total})" if total > 1 else "")
            # the notice closes each block, so a block read on its own still
            # says which CSV it came from and what it left out
            pg_note = note if p == total else ""
            fig, ax = _page_frame(label, caption, suf, pg_note)
            _place_table(ax, sub.iloc[a:b], hdr, fs, numcol)
            pdf.savefig(fig)
            plt.close(fig)

def render_figure_page(pdf, label, caption, png):
    """Portrait A4, with the axes box cut to the figure's own aspect ratio.

    The box used to be a fixed 0.90 x 0.76 of the page and `imshow` letterboxed
    inside it, so a wide, short figure -- S2, S4 and S5 -- sat in the top third
    of an otherwise blank page. Sizing the box to the image removes the gap
    without touching the image itself.
    """
    img = plt.imread(png)
    h_px, w_px = img.shape[0], img.shape[1]
    pw, ph = A4_PORTRAIT
    avail_w, avail_h = 0.90 * pw, 0.76 * ph      # the box the caption leaves free
    scale = min(avail_w / w_px, avail_h / h_px)
    w_in, h_in = w_px * scale, h_px * scale
    x0 = (0.05 * pw + (avail_w - w_in) / 2) / pw
    y0 = (0.04 * ph + (avail_h - h_in) / 2) / ph
    fig = plt.figure(figsize=A4_PORTRAIT)
    draw_caption(fig, caption, title=f"{label} Fig.")
    ax = fig.add_axes([x0, y0, w_in / pw, h_in / ph])
    ax.axis("off")
    ax.imshow(img)
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
    # The CSV is the authoritative full text and the PDF is a rendering of it.
    # Saying so on the title page is what makes the reductions below legitimate
    # rather than lossy: a set name shortened to "REACTOME_SRP_...", a table
    # narrowed to the columns its claim is made of, or a table shown as the rows
    # that fit are each a display choice, and the reader has to be told where the
    # untruncated version is before the choice is theirs to make.
    #
    # This paragraph states the policy in terms of a page budget rather than a
    # row count. The row count it used to name was also written into the
    # manuscript body, so the two could disagree the moment either was edited;
    # a policy that names no number cannot go stale that way.
    fig.text(0.5, 0.30,
             "Every table is also supplied complete as a machine-readable CSV\n"
             "alongside this document, and the CSV, not this PDF, is the\n"
             "authoritative full text. A table too large to print within its\n"
             "page budget is shown as a bounded preview that states its full row\n"
             "and column counts and names the CSV carrying the rest, and the rows\n"
             "the text cites are printed first. Long cell values are shortened\n"
             "with \"...\" in the printed preview only; the CSV holds them whole.",
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
            for panel, src in table_panels(lab):
                tag = f"{lab}{' panel ' + panel if panel else ''}"
                if src:
                    df = pd.read_csv(src)
                    render_table_pages(pdf, lab, caps_tab[lab], df, src, panel)
                    print(f"  tab {tag}: {os.path.basename(src)} ({len(df)}r x {len(df.columns)}c)")
                else:
                    missing.append(f"{lab} Table"
                                   + (f" panel {panel}" if panel else ""))
                    render_placeholder(pdf, lab, "Table",
                                       f"Table file not yet generated.\nExpected "
                                       f"Table_S{lab}{panel or ''}_*.csv",
                                       caps_tab[lab])
                    print(f"  tab {tag}: MISSING -> placeholder")

    size_mb = os.path.getsize(OUT) / 1e6
    print(f"\nwrote {OUT}  ({size_mb:.1f} MB)")
    if missing:
        print(f"placeholders rendered for: {', '.join(missing)}")
        return 0
    print("complete: no placeholder pages.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
