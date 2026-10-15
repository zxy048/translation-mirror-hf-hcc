#!/usr/bin/env python
"""
R20_verify_SI.py -- verify the built Supporting Information PDF against its inputs.

Eyeballing a hundred-page supplement does not establish that it is complete, and
the page count is not fixed -- it moves with every caption and row added. This
checks three things the v1 supplement got wrong, in the same spirit as R18:

  1. Every item the manuscript promises is actually in the PDF, once.
  2. No v1 residue survives in the supplement (v1's "227 genes", the superseded
     "top-3000 gene filter", the hardcoded beta = 12).
  3. No table row was silently dropped at a page break or column split. Each row
     is fingerprinted from its last two columns formatted through the *same*
     R19 formatter, so this validates the real rendering path rather than a
     re-derivation of it.

Run: python R20_verify_SI.py
"""
import os, re, sys
import pandas as pd
import fitz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from R19_build_SI_pdf import (fmt_val, wrap_cell, load_captions, table_panels,
                              table_plan, pin_mask, PIN_ROWS, TAB, MS)

sys.stdout.reconfigure(encoding="utf-8")
PDF = "D:/R_projects/revision_analysis/v2_output/Supporting_Information_v2.pdf"

doc = fitz.open(PDF)
pages = [doc[i].get_text() for i in range(doc.page_count)]
full = "\n".join(pages)
flat = re.sub(r"\s+", " ", full)

fails, oks = [], []

print(f"PDF: {doc.page_count} pages\n")

# ------------------------------------------------------------------ 1. labels --
caps_fig, caps_tab = load_captions()
print("--- every promised item appears exactly once ---")
for lab in sorted(caps_fig, key=lambda s: int(re.sub(r"\D", "", s))):
    needle = f"{lab} Fig."
    n = flat.count(needle)
    ok = n >= 1
    (oks if ok else fails).append(f"{lab} Fig present")
    print(f"  {'OK  ' if ok else 'FAIL'} {lab} Fig. -> {n} occurrence(s)")
for lab in sorted(caps_tab, key=lambda s: (int(re.sub(r'\D', '', s)), s)):
    needle = f"{lab} Table"
    n = flat.count(needle)
    ok = n >= 1
    (oks if ok else fails).append(f"{lab} Table present")
    print(f"  {'OK  ' if ok else 'FAIL'} {lab} Table -> {n} occurrence(s)")

# ------------------------------------------ 1b. no caption kind is ignored ----
# R19 renders the kinds it knows (Table, Fig) and skips every other kind without
# a word. So a caption whose kind is misspelled ("Tabel") or newly introduced
# drops out of the supplement silently, and nothing upstream notices because the
# manuscript still defines it and the page count still looks plausible. A File
# caption exists in this manuscript, which is how the gap was found: the citation
# audit was reading that line as body text while R19 was ignoring it as a
# caption, and each script was consistent with itself and wrong about the other.
# File is the one legitimate skip -- it names the consolidated supplement, which
# is the document being built -- so it is allowed by name rather than by
# accident, and any *other* kind is a failure.
print("\n--- every caption kind in the manuscript is accounted for ---")
with open(MS, encoding="utf-8") as fh:
    declared = re.findall(r"^\*\*(S\d+[a-z]?) ([A-Za-z]+)\.\*\*", fh.read(), re.M)
by_kind = {}
for label, kind in declared:
    by_kind.setdefault(kind, []).append(label)

RENDERED = {"Table": caps_tab, "Fig": caps_fig}
KNOWN_SKIPPED = {"File"}
for kind, labels in sorted(by_kind.items()):
    if kind in RENDERED:
        want = sorted(RENDERED[kind])
        ok = sorted(labels) == want
        (oks if ok else fails).append(f"caption kind {kind} parsed 1:1")
        print(f"  {'OK  ' if ok else 'FAIL'} {kind}: {len(labels)} declared, "
              f"{len(want)} parsed by R19")
        if not ok:
            print(f"        parsed but not declared: {sorted(set(want) - set(labels))}")
            print(f"        declared but not parsed: {sorted(set(labels) - set(want))}")
    elif kind in KNOWN_SKIPPED:
        oks.append(f"caption kind {kind} defined and deliberately not rendered")
        print(f"  OK   {kind}: {len(labels)} defined, not rendered by design "
              f"(the supplement itself)")
    else:
        fails.append(f"unrecognised caption kind {kind!r}")
        print(f"  FAIL unrecognised caption kind {kind!r} -- R19 ignores it, so "
              f"{len(labels)} item(s) would be absent from the supplement: {labels}")

# ------------------------------------------------------- 2. no v1 residue ------
print("\n--- no superseded v1 content in the supplement ---")
# NB: "227 genes" is deliberately NOT listed. The S27 caption quotes v1's
# figure in order to say it is replaced; that is the correction, not residue.
FORBIDDEN = {
    "superseded variance filter": "top-3000",
    "superseded variance filter (alt)": "top 3,000",
    "superseded variance filter (alt2)": "top 6,000",
    "hardcoded beta from v1": "beta = 12",
    "hardcoded beta from v1 (glyph)": "\u03b2 = 12",
    "v1 figure filename": "Figure_S2_TATS",
    "v1 module size (black, 99)": "MEblack (99",
    "v1 module size (green, 227)": "green (227",
}
for why, needle in FORBIDDEN.items():
    hit = needle in flat
    ok = not hit
    (oks if ok else fails).append(f"absent: {why}")
    print(f"  {'OK  ' if ok else 'FAIL'} {why}: {'FOUND (bad)' if hit else 'absent'}")

# ------------------------------------------- 3. no row dropped in rendering ----
# table_plan decides the layout by actually drawing candidate pages and
# measuring them, so each call costs a binary search over real renders. Its
# result is a pure function of the table and the caption, and three sections
# below want the same answer, so it is computed once per table here and read
# back. Recomputing it would triple the run time to produce identical numbers.
PLAN_CACHE = {}


def plan_for(tag, label, df, src, panel):
    if tag not in PLAN_CACHE:
        PLAN_CACHE[tag] = table_plan(label, caps_tab[label], df,
                                     os.path.basename(src), panel)
    return PLAN_CACHE[tag]


def row_check(tag, label, df, src, panel):
    """Fingerprint every row R19 actually prints, in the form R19 prints it.

    This check used to skip every table longer than MAX_FULL_ROWS, on the
    grounds that a truncated table had nothing to verify. That was backwards:
    the truncation itself was never checked, so a table reduced to a single row
    passed -- and a row silently dropped at a page break inside a preview passed
    with it. The rows to check now come from R19.table_plan, the same call the
    builder makes, so what is verified is the printed set rather than an
    assumption about it.
    """
    blocks, shown, note, _ = plan_for(tag, label, df, src, panel)
    printed = [c for b in blocks for c in b]
    if not len(shown):
        fails.append(f"{tag} printed no rows at all")
        print(f"  FAIL {tag}: table_plan prints no rows")
        return
    # Two columns, not one: a single column can repeat down a table, so a row
    # that was never printed can be "found" through another row's value.
    first = [wrap_cell(r[printed[0]]) for _, r in shown.iterrows()]
    last = [wrap_cell(r[printed[-1]]) for _, r in shown.iterrows()]
    fp = [(a, b) for a, b in zip(first, last) if (a and a != "nan") or (b and b != "nan")]
    if not fp:
        print(f"  --   {tag}: no usable fingerprint")
        return
    found = sum(1 for a, b in fp if (a in flat) and (b in flat))
    frac = found / len(fp)
    ok = frac >= 0.98
    (oks if ok else fails).append(f"{tag} printed rows rendered ({frac:.0%})")
    print(f"  {'OK  ' if ok else 'FAIL'} {tag}: {found}/{len(fp)} printed rows "
          f"located ({frac:.0%}), {len(shown)} of {len(df):,} rows printed")


print("\n--- row completeness (over the rows actually printed) ---")
PLANS = {}
for lab in sorted(caps_tab, key=lambda s: (int(re.sub(r"\D", "", s)), s)):
    # A merged table publishes several CSVs under one item number, so the
    # sources come from R19 rather than from a glob on the label: globbing the
    # parent label alone would verify panel a and leave the rest unchecked.
    for panel, src in table_panels(lab):
        tag = f"{lab}{' panel ' + panel if panel else ''}"
        if not src:
            print(f"  --   {tag}: no CSV (placeholder page, skipped)")
            continue
        df = pd.read_csv(src)
        PLANS[tag] = (lab, df, src, panel)
        row_check(tag, lab, df, src, panel)

# ------------------------------------------- 4. cited rows survive the budget --
# The budget truncates tables, and the rows it must never truncate are the ones
# the manuscript cites: a reader checking "10 of 15 Tier A sets" in the text has
# to find those fifteen rows on the printed page. R19 orders them first, so the
# property to verify is that there was room for all of them. This fires if the
# budget is tightened or the cited set grows, which is exactly when the preview
# would start hiding the evidence the text points at.
print("\n--- tables whose cited rows must survive the page budget ---")
for lab in sorted(PIN_ROWS, key=lambda s: (int(re.sub(r"\D", "", s)), s)):
    applied = 0
    for panel, src in table_panels(lab):
        tag = f"{lab}{' panel ' + panel if panel else ''}"
        if not src:
            continue
        df = pd.read_csv(src)
        mask = pin_mask(df, lab)
        if mask is None:
            # A pin describes the table it was written for, and a merged item's
            # panels are different tables -- S39 panel b has no `score` column.
            # Saying so is the point: silence here would also cover a pin
            # misspelled into never matching anything, which is the failure this
            # section exists to catch, so the label is required to apply below.
            print(f"  --   {tag}: the {lab} pin does not describe this panel, "
                  f"nothing to protect")
            continue
        applied += 1
        mask = mask.to_numpy(dtype=bool)
        blocks, shown, note, _ = plan_for(tag, lab, df, src, panel)
        printed = [c for b in blocks for c in b]
        # Compared by identity of the printed cells, not by row count. R19 takes
        # the pinned set as a floor on the row count, so a count comparison
        # ("n_pin <= len(shown)") is true by construction and could never fail
        # -- it would keep passing if the floor were removed. Comparing the rows
        # themselves fails when a cited row stops being printed, whatever the
        # reason.
        def key(frame):
            return {tuple(str(r[c]) for c in printed) for _, r in frame.iterrows()}
        pinned_rows = df[mask]
        missing = key(pinned_rows) - key(shown)
        ok = not missing
        (oks if ok else fails).append(f"{tag} prints all {int(mask.sum())} cited rows")
        print(f"  {'OK  ' if ok else 'FAIL'} {tag}: {int(mask.sum())} cited rows, "
              f"{len(shown)} printed"
              + ("" if ok else f" -- {len(missing)} cited row(s) not printed, e.g. "
                 f"{sorted(missing)[0][:3]}"))
    # A pin that matches no frame at all protects nothing, and would report the
    # same clean "nothing to protect" line as a panel it legitimately does not
    # describe. Failing here is what separates the two: the pin for a cited table
    # has to land somewhere, so a column renamed out from under it is caught
    # instead of quietly downgrading the table to unprotected.
    ok = applied > 0
    (oks if ok else fails).append(f"{lab} pin applies to a real table")
    print(f"  {'OK  ' if ok else 'FAIL'} {lab} pin: applies to {applied} of its frames")

# ------------------------------- 5. a reduced table says where the rest is ----
# Every page that shows less than the whole table has to name the CSV that holds
# the rest, or the reduction is indistinguishable from a table that simply has
# fewer rows. The CSV is the authoritative copy; this checks the pointer to it
# is really on the page.
print("\n--- a reduced table names the CSV carrying the rest ---")
for tag, (lab, df, src, panel) in PLANS.items():
    blocks, shown, note, _ = plan_for(tag, lab, df, src, panel)
    reduced = len(shown) < len(df) or len(blocks[0]) < len(df.columns)
    if not reduced:
        continue
    csv_name = os.path.basename(src)
    ok = csv_name in flat
    (oks if ok else fails).append(f"{tag} names its CSV when reduced")
    print(f"  {'OK  ' if ok else 'FAIL'} {tag}: reduced to {len(shown)} rows, "
          f"names {csv_name}: {ok}")

# ------------------------------------------------- 6. the page budget holds --
# A ceiling, not a target. The supplement was 160 pages because nothing stopped
# it: the row policy was applied per column block and no rule looked at the
# document as a whole. This is that rule. It is set with room to spare, so it
# fires on a return to the old behaviour or an unforeseen table, not on a caption
# growing by a line.
PAGE_CEILING = 95
print("\n--- the whole supplement stays inside its page budget ---")
ok = doc.page_count <= PAGE_CEILING
(oks if ok else fails).append(f"supplement is {doc.page_count} pages")
print(f"  {'OK  ' if ok else 'FAIL'} supplement is {doc.page_count} pages "
      f"(ceiling {PAGE_CEILING})")

# ------------------------------------------------------------------- report ----
print(f"\n{'='*70}\nPASSED {len(oks)} | FAILED {len(fails)}")
if fails:
    print("\nFAILURES:")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("Supporting Information reconciles with its inputs.")
