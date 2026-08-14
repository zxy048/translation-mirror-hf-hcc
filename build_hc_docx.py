# -*- coding: utf-8 -*-
"""
Build Hepatology Communications manuscript docx:
  1. Read Manuscript_HC_v1.md
  2. Convert HTML <sup>...</sup> -> ^...^ and <sub>...</sub> -> ~...~
     (pandoc renders RawInline HTML without vertAlign; native syntax does apply
      superscript/subscript in the docx writer)
  3. pandoc md -> docx (--from markdown+superscript+subscript)
  4. format_hc.py (Times New Roman 12pt, double-spaced, Letter, line numbers)
"""
import re
import subprocess
import sys

SRC_MD = r"D:\R_projects\revision_analysis\Manuscript_HC_v1.md"
TMP_MD = r"D:\R_projects\revision_analysis\Manuscript_HC_v1_build.md"
DST_DOCX = r"D:\R_projects\revision_analysis\Manuscript_HC_v1.docx"

with open(SRC_MD, encoding="utf-8") as f:
    md = f.read()

md = re.sub(r"<sup>(.*?)</sup>", r"^\1^", md)
md = re.sub(r"<sub>(.*?)</sub>", r"~\1~", md)

with open(TMP_MD, "w", encoding="utf-8") as f:
    f.write(md)

cmd = [
    "pandoc", TMP_MD, "-o", DST_DOCX,
    "--from", "markdown+superscript+subscript",
]
subprocess.run(cmd, check=True)
print(f"pandoc OK: {DST_DOCX}")

from format_hc import format_hc  # noqa: E402
format_hc(DST_DOCX)
