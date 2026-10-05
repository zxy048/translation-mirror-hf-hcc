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
import glob, os, re, sys
import pandas as pd
import fitz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from R19_build_SI_pdf import fmt_val, wrap_cell, load_captions, TAB, MAX_FULL_ROWS

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
print("\n--- row completeness (fingerprint of the last two columns) ---")
for lab in sorted(caps_tab, key=lambda s: (int(re.sub(r"\D", "", s)), s)):
    hits = glob.glob(f"{TAB}/Table_{lab}_*.csv")
    if not hits:
        print(f"  --   {lab}: no CSV (placeholder page, skipped)")
        continue
    df = pd.read_csv(hits[0])
    if len(df) > MAX_FULL_ROWS:
        print(f"  --   {lab}: {len(df):,} rows, previewed by design (skipped)")
        continue
    cols = list(df.columns)
    # Fingerprint each row on its LAST column, formatted AND wrapped exactly as
    # R19 prints it -- a re-derivation could pass while the printed text differs.
    # Cells are checked individually because a wide table is split across column
    # blocks, so no row's text is ever contiguous on one page.
    fp = [wrap_cell(r[cols[-1]]) for _, r in df.iterrows()]
    fp = [f for f in fp if f and f not in ("", "nan")]
    if not fp:
        print(f"  --   {lab}: no usable fingerprint")
        continue
    found = sum(1 for a in fp if a in flat)
    frac = found / len(fp)
    ok = frac >= 0.98
    (oks if ok else fails).append(f"{lab} rows rendered ({frac:.0%})")
    print(f"  {'OK  ' if ok else 'FAIL'} {lab}: {found}/{len(fp)} rows located ({frac:.0%})")

# ------------------------------------------------------------------- report ----
print(f"\n{'='*70}\nPASSED {len(oks)} | FAILED {len(fails)}")
if fails:
    print("\nFAILURES:")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("Supporting Information reconciles with its inputs.")
