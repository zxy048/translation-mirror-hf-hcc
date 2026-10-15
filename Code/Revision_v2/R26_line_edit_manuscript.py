#!/usr/bin/env python
"""
R26_line_edit_manuscript.py -- apply the cross-check line-editing pass to the manuscript.

The edits are applied as exact-string replacements rather than by hand, for the
same reason R25 gives for the response letter: a small language pass over a long
document is exactly the kind of work where one silent miss is invisible -- the
sentence still reads, it just keeps the defect. Each pair below asserts it
matches exactly once, so a pair that has drifted fails loudly instead of
no-op'ing. Nothing is written unless every pair matched.

Provenance of the items (all were verified against the sources before the
strings were written, not taken on report):

  M1  triple "against" in one sentence -- the same preposition doing three jobs
      ("correlations run ... against X ... against Y, against 9 of 14"). Only
      the third is the operative comparison; the first two compare one
      correlation against a cohort, which "with" states without the collision.
  M2  the mirror count is given twice, for two cohorts, in identical words
      ("GSE141910 1 of 15 both times, GSE116250 1 of 15 both times"). Both facts
      are correct and both are kept; only the duplicated clause is merged.
      R17_log.txt:26-27 is the source of the underlying ranks.
  M3/M4  -ise house style. The manuscript is British throughout ("harmonised",
      "labelled", "modelled"); these two were the only -ize stragglers.
  M5  "four cohorts" is wrong arithmetic for the row it describes. The four
      units listed at L153 are TCGA-LIHC, GSE14520_GPL3921, GSE14520_GPL571 and
      GSE76427 -- three studies, four dataset-arms. The cover letter was already
      corrected to "four dataset-arms"; this syncs the manuscript to it. No
      statistic changes: the count of four, and the rho range, are untouched.
  M6  "scored sets" is used in three places (Fig 2a, Fig 3a legends and the
      Results) and defined in none of them, so a reader cannot recover the 41.
      The pool is built in R06_ssgsea_cohorts.R:24 as the union of the three
      tiers, then filtered at L40-46 to sets with >=5 measured genes
      (minSize = 5). Table_S13_scored_set_counts.csv shows the resulting count
      is cohort-specific -- 41 for GSE57338/TCGA_LIHC, but 40, 36, 35 and 34 for
      the others -- so the caption says which two cohorts the 41 belongs to.
  M7  the y axis of Fig 1d is median rank, where smaller means stronger. That is
      conventional in WGCNA and opaque to everyone else. One parenthetical.

Run: python R26_line_edit_manuscript.py [--dry-run]
"""
import sys

MD = "D:/R_projects/revision_analysis/Manuscript_PLOSONE_v2.md"

PAIRS = [
    # -----------------------------------------------------------------------
    # M1  triple "against"
    # -----------------------------------------------------------------------
    ("and the effect-size correlations run from \u22120.393 to +0.211 against "
     "GSE116250 and from \u22120.457 to \u22120.059 against GSE141910, against "
     "9 of 14 to 10 of 15",
     "and the effect-size correlations run from \u22120.393 to +0.211 with "
     "GSE116250 and from \u22120.457 to \u22120.059 with GSE141910, against "
     "9 of 14 to 10 of 15"),

    # -----------------------------------------------------------------------
    # M2  duplicated cohort clause
    # -----------------------------------------------------------------------
    ("GSE57338 gives 10 of 15 Tier A sets before and after adjustment, "
     "GSE141910 1 of 15 both times, GSE116250 1 of 15 both times \u2014",
     "GSE57338 gives 10 of 15 Tier A sets before and after adjustment, "
     "GSE116250 and GSE141910 1 of 15 both times \u2014"),

    # -----------------------------------------------------------------------
    # M3/M4  -ise house style
    # -----------------------------------------------------------------------
    ("is increasingly being recognized as a targetable vulnerability",
     "is increasingly being recognised as a targetable vulnerability"),

    ("is characterized by suppression of anabolic processes",
     "is characterised by suppression of anabolic processes"),

    # -----------------------------------------------------------------------
    # M5  three studies, four dataset-arms
    # -----------------------------------------------------------------------
    ("(\u03c1 = \u22120.55 to \u22120.76 across four cohorts and three platforms)",
     "(\u03c1 = \u22120.55 to \u22120.76 across four cohort arms and three platforms)"),

    # -----------------------------------------------------------------------
    # M6  define "scored sets" where the term is first used
    # -----------------------------------------------------------------------
    ("for the 15 Tier A translation sets and the other scored sets (grey).",
     "for the 15 Tier A translation sets and the other scored sets (grey; "
     "scored sets are the union of Tier A, Tier B and Tier C restricted to sets "
     "with at least five measured genes in that cohort \u2014 41 here)."),

    # -----------------------------------------------------------------------
    # M7  which way the median-rank axis runs
    # -----------------------------------------------------------------------
    ("*Z*<sub>summary</sub> against median rank for every named module, "
     "200 permutations; dotted lines mark *Z* = 2 and *Z* = 10.",
     "*Z*<sub>summary</sub> against median rank for every named module, "
     "200 permutations (lower median rank indicates stronger preservation); "
     "dotted lines mark *Z* = 2 and *Z* = 10."),
]


def main():
    dry = "--dry-run" in sys.argv
    with open(MD, encoding="utf-8") as f:
        body = f.read()

    ok = True
    for i, (old, new) in enumerate(PAIRS, 1):
        n = body.count(old)
        if n != 1:
            print(f"  FAIL M{i}: {n} occurrence(s), expected exactly 1")
            print(f"        {old[:90]}...")
            ok = False
        else:
            print(f"  ok   M{i}: 1 occurrence")

    if not ok:
        print("\nNothing written: a pair did not match exactly once.")
        return 1

    if dry:
        print("\n--dry-run: all pairs matched; nothing written.")
        return 0

    for old, new in PAIRS:
        body = body.replace(old, new, 1)
    with open(MD, "w", encoding="utf-8", newline="") as f:
        f.write(body)
    print(f"\nWrote {len(PAIRS)} replacements to {MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
