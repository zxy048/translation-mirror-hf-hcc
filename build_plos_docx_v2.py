# -*- coding: utf-8 -*-
"""
Build the PLOS ONE revision docx from Manuscript_PLOSONE_v2.md:
  1. Read the v2 markdown
  2. Convert HTML <sup>...</sup> -> ^...^ and <sub>...</sub> -> ~...~
  3. pandoc md -> docx (--from markdown+superscript+subscript)
  4. format_plos.format_plos (TNR 12pt, double spacing, Letter, 1-in margins,
     continuous line numbers, footer page numbers, page break before Abstract)

This mirrors build_plos_docx.py, which stays v1-hardcoded on purpose: the Word
"Compare Documents" step needs a regenerable v1 docx to diff the revision
against, so the two versions must not share a mutable source of truth.

Differences from v1 that the downstream QA (qa_plos_v2.py) asserts, because a
silent carry-over of v1's counts would pass a copy-pasted check while being
wrong: v2 cites 25 references (v1: 26), has 3 main figures (v1: 5), and adds
Data Availability, Author Contributions, Funding, Competing Interests and
Acknowledgments sections (v1: 12 top-level sections, v2: 14).
"""
import os
import re
import subprocess
import sys

ROOT = r"D:\R_projects\revision_analysis"
SRC_MD = os.path.join(ROOT, "Manuscript_PLOSONE_v2.md")
TMP_MD = os.path.join(ROOT, "Manuscript_PLOSONE_v2_build.md")
DST_DOCX = os.path.join(ROOT, "Manuscript_PLOSONE_v2.docx")

with open(SRC_MD, encoding="utf-8") as f:
    md = f.read()

# pandoc's superscript/subscript extensions take ^...^ / ~...~ with no spaces
# inside the delimiters; the manuscript's tags never span a space, so a direct
# substitution is safe here (checked below rather than assumed).
bad = re.findall(r"<(?:sup|sub)>[^<]*\s[^<]*</(?:sup|sub)>", md)
if bad:
    print(f"WARNING: {len(bad)} sup/sub tag(s) contain whitespace and will not "
          f"convert cleanly, e.g. {bad[0]!r}")

md = re.sub(r"<sup>(.*?)</sup>", r"^\1^", md)
md = re.sub(r"<sub>(.*?)</sub>", r"~\1~", md)

if "<sup>" in md or "<sub>" in md:
    sys.exit("ERROR: unconverted sup/sub tags remain after substitution")

with open(TMP_MD, "w", encoding="utf-8") as f:
    f.write(md)

cmd = ["pandoc", TMP_MD, "-o", DST_DOCX,
       "--from", "markdown+superscript+subscript"]
subprocess.run(cmd, check=True)
print(f"pandoc OK: {DST_DOCX}")

sys.path.insert(0, ROOT)
from format_plos import format_plos  # noqa: E402

format_plos(DST_DOCX)
print(f"\nbuilt {DST_DOCX}")
