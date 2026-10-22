# -*- coding: utf-8 -*-
"""Apply the 2026-10-22 cross-document micro-patch to the v2 manuscript.

Five text-level replacements, each asserted to occur exactly once before it is
applied. No analysis output, table, figure, response letter or cover letter is
touched. The manuscript remains the single source for the SI captions; the SI
is rebuilt separately by R19_build_SI_pdf.py.

Edits:
  N-1  Results: network-construction rationale no longer implies a causal
       guarantee; wording aligned with the cover letter ("held constant").
  N-2  Results: drop "activity-based" from the secondary designation criterion.
  N-4  S8 Table caption: same removal of "activity-based".
  N-7  S3 Table caption: the table lists all Tier A members with retention
       status, not only genes lost from the common universe.
  B    Discussion: "mutually consistent" -> "compatible", matching the
       observational strength of overlapping confidence intervals.
"""
import hashlib
import os
import re
import sys

ROOT = r"D:/R_projects/revision_analysis"
MS = os.path.join(ROOT, "Manuscript_PLOSONE_v2.md")

EDITS = [
    (
        "N-1",
        "so that disease state could not drive either network's module structure",
        "so that disease status is held constant within each network",
    ),
    (
        "N-2",
        "As an activity-based secondary criterion",
        "As a secondary criterion",
    ),
    (
        "N-4",
        "this was the secondary (activity-based) designation criterion",
        "this was the secondary designation criterion",
    ),
    (
        "N-7",
        "**S3 Table.** Tier A gene loss. Members of the 15 Tier A sets that are absent from the common universe *U*, listed per set with the reason for exclusion, so that the effect of the universe definition on the translation sets is auditable.",
        "**S3 Table.** Tier A gene retention and loss. All members of the 15 Tier A sets, with retention status and reason for loss, so that the effect of the universe definition on the translation sets is auditable.",
    ),
    (
        "B",
        "All four confidence intervals overlap, so the cohorts are mutually consistent rather than contradictory.",
        "All four confidence intervals overlap, so the cohorts are compatible rather than contradictory.",
    ),
]


def line_of(text, pos):
    return text.count("\n", 0, pos) + 1


def main():
    with open(MS, "r", encoding="utf-8", newline="") as fh:
        text = fh.read()
    before_sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    original = text

    for tag, old, new in EDITS:
        n = text.count(old)
        if n != 1:
            print(f"ABORT {tag}: expected 1 occurrence, found {n}")
            sys.exit(1)
        pos = text.index(old)
        print(f"{tag}: line {line_of(text, pos)}: {old[:70]}...")
        text = text.replace(old, new, 1)

    # The two "activity-based" occurrences must both be gone.
    leftover = [i + 1 for i, line in enumerate(text.splitlines()) if "activity-based" in line]
    if leftover:
        print(f"ABORT: 'activity-based' remains at lines {leftover}")
        sys.exit(1)

    # Guard rails: no change to the abstract block or to any numeric result.
    abs_before = re.search(r"(?ms)^## Abstract\s*\n(.*?)(?=^## )", original)
    abs_after = re.search(r"(?ms)^## Abstract\s*\n(.*?)(?=^## )", text)
    if abs_before is None or abs_after is None or abs_before.group(1) != abs_after.group(1):
        print("ABORT: abstract changed")
        sys.exit(1)

    with open(MS, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
    after_sha = hashlib.sha256(text.encode("utf-8")).hexdigest()

    print(f"before {before_sha}")
    print(f"after  {after_sha}")
    print(f"bytes  {len(original.encode('utf-8'))} -> {len(text.encode('utf-8'))}")
    print("OK: 5 replacements applied")


if __name__ == "__main__":
    main()
