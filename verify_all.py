# -*- coding: utf-8 -*-
"""
Single entry point for every cross-artifact check in this project.

WHY THIS FILE EXISTS
--------------------
The v1 submission was checked by `qa_plos.py`, which read exactly two files --
the manuscript .md and its .docx -- and asserted twelve properties of them:
abstract length, citation numbering, reference count, section presence, page
size, margins, line spacing, line numbers, footer page number. Every one of
those is a property of the manuscript text alone.

The defects that mattered were not, and none of them was caught:

  * "all 33 translation-related pathways showed positive effect sizes" is
    refuted by the submission's own S2 table, which contains six negatives.
    Catching it requires opening the table, not the manuscript.
  * The reported binomial test P < 0.0001 had no implementation anywhere in the
    codebase. Catching it requires grepping the code, not reading the text.
  * The HF soft threshold was hard-coded in a script that is public.
  * The submitted supplementary PDF was an abandoned analysis that contradicted
    the manuscript. qa_plos.py never opened it.

The fix is not "be more careful". Carefulness is not checkable. The fix is to
make verification cross-artifact by construction -- every check below reads at
least one file that is NOT the manuscript -- and to make running the checks a
single command that cannot be partially remembered or partially skipped.

Run this before assembling any submission package. It exits non-zero if any
check fails, so it can gate a release rather than merely inform one.
"""
import os
import re
import subprocess
import sys
import zlib

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

# (label, script path relative to ROOT, what it reads that the manuscript does not)
CHECKS = [
    ("manuscript numbers vs regenerated tables",
     "Code/Revision_v2/R18_verify_manuscript.py",
     "opens the regenerated S-table CSVs and re-derives every number quoted in the text"),
    ("supporting information",
     "Code/Revision_v2/R20_verify_SI.py",
     "opens the built SI PDF and checks every table and figure it holds against its source CSV"),
    ("manuscript structure and PLOS formatting",
     "qa_plos_v2.py",
     "checks abstract length, citation order, section presence, page setup, line numbers"),
    ("tracked-changes docx round-trip",
     "build_tracked_docx.py",
     "rebuilds the tracked docx and confirms accept-all reproduces v2 and reject-all reproduces v1"),
]

# v1 figures were rendered at 13 x 11.5 in and downscaled by the journal to
# 7.5 x 8.75 in, which shrank every label with them. The ceiling is asserted
# here rather than trusted, because "we rendered it at the right size" is
# exactly the kind of claim that silently goes stale after a re-render.
PLOS_MAX_W_IN = 7.5
PLOS_MAX_H_IN = 8.75
FIGURES = ["PLOS_ONE_Submission_v2/Fig1.pdf",
           "PLOS_ONE_Submission_v2/Fig2.pdf",
           "PLOS_ONE_Submission_v2/Fig3.pdf"]


def media_box_inches(path):
    """Width and height in inches from a PDF's first /MediaBox.

    Read raw bytes rather than with a PDF library, so this needs no dependency
    the checking environment might lack. The first version of this only scanned
    the raw file, which works for R's own pdf() output but silently reports
    "no /MediaBox found" for the PDF 1.7 files written here, where the page
    dictionary lives inside a FlateDecode object stream. A check that cannot
    find its input must not read as a pass, so the streams are decompressed
    before giving up.
    """
    with open(path, "rb") as fh:
        blob = fh.read()

    pattern = re.compile(rb"/MediaBox\s*\[\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)")

    def first_box(data):
        m = pattern.search(data)
        return m.groups() if m else None

    groups = first_box(blob)
    if groups is None:
        for m in re.finditer(rb"stream\r?\n", blob):
            start = m.end()
            end = blob.find(b"endstream", start)
            if end < 0:
                continue
            try:
                data = zlib.decompress(blob[start:end])
            except zlib.error:
                continue
            groups = first_box(data)
            if groups:
                break
    if groups is None:
        return None
    x0, y0, x1, y1 = (float(v) for v in groups)
    return abs(x1 - x0) / 72.0, abs(y1 - y0) / 72.0


def check_figure_specs():
    """Every main figure must sit at or under the PLOS ceiling."""
    fails = []
    for rel in FIGURES:
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            fails.append("%s is missing" % rel)
            continue
        box = media_box_inches(path)
        if box is None:
            fails.append("%s: no /MediaBox found" % rel)
            continue
        w, h = box
        ok = w <= PLOS_MAX_W_IN + 0.01 and h <= PLOS_MAX_H_IN + 0.01
        print("  %s %s  %.2f x %.2f in (%.0f x %.0f px @300dpi)"
              % ("PASS" if ok else "FAIL", rel, w, h, w * 300, h * 300))
        if not ok:
            fails.append("%s is %.2f x %.2f in, over the %.2f x %.2f ceiling"
                         % (rel, w, h, PLOS_MAX_W_IN, PLOS_MAX_H_IN))
    return fails


def run(label, script, reads):
    path = os.path.join(ROOT, script)
    print("\n" + "=" * 74)
    print("[%s]" % label)
    print("  reads: %s" % reads)
    print("=" * 74)
    if not os.path.exists(path):
        print("  FAIL  script not found: %s" % script)
        return ["%s: %s not found" % (label, script)]
    # Console is GBK on this machine; the scripts print non-ASCII and die on
    # encode errors unless the child is told the stream is UTF-8.
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    p = subprocess.run([PY, path], cwd=ROOT, env=env,
                       capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    out = (p.stdout or "") + (p.stderr or "")
    for line in out.rstrip().splitlines()[-14:]:
        print("  | " + line)
    if p.returncode != 0:
        return ["%s failed (exit %d)" % (label, p.returncode)]
    return []


def main():
    print("Cross-artifact verification.")
    print("Every check below reads at least one file that is not the manuscript.")
    print("A check that only reads the manuscript cannot see the failures that")
    print("matter, which is how the v1 submission shipped six of them.")
    failures = []
    for label, script, reads in CHECKS:
        failures += run(label, script, reads)

    print("\n" + "=" * 74)
    print("[figure size against the PLOS ceiling]")
    print("=" * 74)
    failures += check_figure_specs()

    print("\n" + "=" * 74)
    if failures:
        print("VERIFICATION FAILED -- %d problem(s):" % len(failures))
        for f in failures:
            print("  * " + f)
        return 1
    print("ALL CROSS-ARTIFACT CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
