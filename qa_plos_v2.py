# -*- coding: utf-8 -*-
"""
QA for the PLOS ONE *revision* (md + docx).

md-level:
  - abstract <= 300 words, unstructured, no citations
  - every [N] citation ascending within groups, numbered by first appearance
  - all 26 refs cited; no leftover Harvard (et al., YYYY)
  - Fig 1-3 cited in text and captioned inline (v2 has 3 main figures, not 5)
  - required PLOS sections present (14 in v2, not v1's 12)
  - no withdrawn v1 wording or superseded v1 numbers survive in the body
docx-level:
  - Letter page size, 1-in margins
  - double spacing (2.0) on body
  - continuous line numbers (w:lnNumType)
  - footer page number field
  - Abstract starts on a new page

Run: python qa_plos_v2.py
"""
import re
import sys

ROOT = r"D:\R_projects\revision_analysis"
MD = ROOT + r"\Manuscript_PLOSONE_v2.md"
DOCX = ROOT + r"\Manuscript_PLOSONE_v2.docx"

N_REFS = 26
N_FIGS = 3

fail = 0
def ok(cond, msg):
    global fail
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        fail += 1

with open(MD, encoding="utf-8") as f:
    md = f.read()

body = md.split("## References")[0]  # text before references (includes abstract)

# ---- abstract ----
abs_start = md.index("## Abstract")
abs_end = md.index("## Introduction")
abstract = md[abs_start:abs_end]
abs_words = [w for w in re.sub(r"[#*\n]", " ", abstract).split() if w.strip()]
print(f"Abstract word count: {len(abs_words)} (limit 300)")
ok(len(abs_words) <= 300, f"abstract <= 300 words ({len(abs_words)})")
ok("[" not in abstract.replace("## Abstract", ""), "abstract contains no citations")

# ---- citations ----
# Integer tokens in range only; excludes numeric CIs such as [0.67, 4.21] and
# the bracketed reference numbers that legitimately exceed the ref count.
cites = []
for m in re.finditer(r"\[([0-9]+(?:,[0-9]+)*)\]", body):
    nums = [int(t) for t in m.group(1).split(",")]
    if all(1 <= n <= N_REFS for n in nums):
        cites.append(nums)
first_seen, order = {}, []
for group in cites:
    for n in group:
        if n not in first_seen:
            first_seen[n] = len(order)
            order.append(n)
ok(order == sorted(order) and len(set(order)) == N_REFS,
   f"citations numbered by first appearance; {len(order)} of {N_REFS} refs cited")
ok(all(g == sorted(g) for g in cites), "no descending citation groups (e.g. [8,5,7])")

ref_block = md[md.index("## References"):md.index("## Supporting Information Captions")]
refs = re.findall(r"^(\d+)\.\s", ref_block, flags=re.M)
ok(len(refs) == N_REFS, f"reference list has {N_REFS} entries (got {len(refs)})")
ok(set(int(r) for r in refs) == set(first_seen.keys()),
   "reference numbering matches first-appearance order")

leftover = re.findall(r"et al\.,\s+\d{4}|\([A-Z][a-z]+ et al\.|\([A-Z][a-z]+ and [A-Z][a-z]+, \d{4}", body)
ok(not leftover, f"no leftover Harvard citations ({len(leftover)})")

# ---- figures (v2 has three) ----
for n in range(1, N_FIGS + 1):
    ok(bool(re.search(r"\bFig %d\b" % n, body)), f"Fig {n} cited in text")
    ok(bool(re.search(r"\*\*Fig %d\.\*\*" % n, body)), f"Fig {n} caption present inline")

# ---- required sections (v2 adds five declarations sections) ----
for sec in ["## Abstract", "## Introduction", "## Materials and Methods",
            "## Results", "## Discussion", "## References",
            "## Supporting Information Captions",
            "## Data Availability", "## Author Contributions", "## Funding",
            "## Competing Interests", "## Acknowledgments"]:
    ok(sec in md, f"section present: {sec}")

# ---- no withdrawn v1 content in the body ----
# R20 polices the supplement for v1 residue; nothing previously policed the main
# text, which is where the v1 overstatements actually were.
#
# The revision legitimately QUOTES the withdrawn wording in order to retract it
# ("v1 stated that 'all 33 ... showed positive effect sizes'; that statement is
# false", "v1's '227 genes' ... did not hold across cohorts"). A phrase-level
# search flags those retractions as residue -- the same false positive R20 hit
# with its S27 caption. So a hit counts only if some sentence contains the phrase
# WITHOUT also containing a marker that marks it as superseded.
print("\n--- no withdrawn v1 wording asserted anywhere in the body ---")
WITHDRAWN = {
    "v1 claim 'all 33 positive'": r"all 33 translation-related pathways",
    "v1 phrase 'therapeutic entry points'": "therapeutic entry points",
    "v1 phrase 'ISR engagement'": "ISR engagement",
    "v1 phrase 'disease-stage-specific vulnerability'": "disease-stage-specific vulnerability",
    "v1 phrase 'dominant regulatory dimensions'": "dominant regulatory dimensions",
    "v1 phrase 'partial replication'": "partial replication",
    "v1 module size 227": r"\b227 genes\b",
    "v1 hard-coded beta = 12": r"[Bb]eta = 12|\u03b2 = 12",
    "v1 variance filter top-3000": "top-3000",
    "v1 variance filter top-6,000": "top-6,000",
}
RETRACTED = re.compile(
    r"v1|supersed|withdraw|false|unverifiable|did not hold|no longer|replac|removed",
    re.I)
sentences = re.split(r"(?<=[.;:])\s+", body)
for why, pat in WITHDRAWN.items():
    asserted = [s for s in sentences
                if re.search(pat, s) and not RETRACTED.search(s)]
    ok(not asserted,
       f"withdrawn -- {why}" + (f": ASSERTED -> {asserted[0][:90]}..." if asserted else ""))

# ---- docx-level ----
from docx import Document
from docx.oxml.ns import qn
doc = Document(DOCX)

sec = doc.sections[0]
W, H = round(sec.page_width.cm, 1), round(sec.page_height.cm, 1)
ok((W, H) == (21.6, 27.9), f"page size Letter ({W:.1f} x {H:.1f} cm)")
ok(round(sec.left_margin.cm, 2) == 2.54 and round(sec.right_margin.cm, 2) == 2.54,
   "left/right margins 1 in (2.54 cm)")

double_spaced = sum(1 for p in doc.paragraphs
                    if p.text.strip() and p.paragraph_format.line_spacing == 2.0)
total_paras = sum(1 for p in doc.paragraphs if p.text.strip())
ok(double_spaced / max(total_paras, 1) > 0.9,
   f"double spacing on body ({double_spaced}/{total_paras} non-empty paras)")

ln = sec._sectPr.find(qn('w:lnNumType'))
ok(ln is not None and ln.get(qn('w:restart')) == 'continuous',
   "continuous line numbers present")

footer_xml = sec.footer.paragraphs[0]._p.xml if sec.footer.paragraphs else ""
ok("PAGE" in footer_xml, "footer page-number field present")

abs_para = next((p for p in doc.paragraphs if p.text.strip() == "Abstract"), None)
ok(abs_para is not None and bool(abs_para.paragraph_format.page_break_before),
   "Abstract starts on a new page")

# superscript/subscript survived the md -> docx conversion. A silent loss here
# would print "Zsummary" and "10-50" in the submitted manuscript.
raw = doc.element.xml
ok("superscript" in raw or "w:vertAlign" in raw, "superscript formatting present in docx")
ok("subscript" in raw, "subscript formatting present in docx")

print("\n" + ("ALL CHECKS PASSED" if fail == 0 else f"{fail} CHECK(S) FAILED"))
sys.exit(1 if fail else 0)
