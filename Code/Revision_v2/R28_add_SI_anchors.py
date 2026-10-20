#!/usr/bin/env python
"""
R28_add_SI_anchors.py -- add the 13 approved in-text citations to the SI items.

These are the 13 Supporting Information items that had a caption but no body
citation. The permanent gate in verify_all.py (check_si_citations) fails on
them, so this script is what turns the gate green. Each pair asserts exactly one
match; nothing is written unless every pair matches.

All edits are pure pointer insertions: no number changes, no new claim, no
caption change. Two items were revised at approval:

  * S5 Fig is a heterogeneity figure (I-squared, Cochran's Q), so it is anchored
    to the Methods sentence that defines that quantity, beside S21 -- not to
    S22, which is the aetiology table.
  * S1 Fig and S31 Table are two halves of one analysis (threshold figure and
    underlying counts). They are cited together, in the sentence supplied at
    approval, which also puts S31 Table into the body for the first time.

Run: python R28_add_SI_anchors.py [--dry-run]
"""
import sys

MD = "D:/R_projects/revision_analysis/Manuscript_PLOSONE_v2.md"

PAIRS = [
    # --- S3 Table: common universe, after the S1 Table citation ------------
    ("The complete loss chain from platform feature to `U` is reported per "
     "platform (S1 Table).",
     "The complete loss chain from platform feature to `U` is reported per "
     "platform (S1 Table). Members of the Tier A sets that are absent from `U` "
     "are listed per set in S3 Table."),

    # --- S3 Fig: soft-threshold section, beside S5 Table -------------------
    ("at which the HCC fit is R\u00b2 = 0.778 (S5 Table).",
     "at which the HCC fit is R\u00b2 = 0.778 (S5 Table; fit curves in S3 Fig)."),

    # --- S6 Table: Samples used, the reason the table exists --------------
    ("Sensitivity networks built on all 313 GSE57338 samples and on TCGA-LIHC "
     "are reported separately.",
     "Sensitivity networks built on all 313 GSE57338 samples and on TCGA-LIHC "
     "are reported separately. The module\u2013clinical-trait relationships that "
     "remain definable in the HF network are reported in S6 Table."),

    # --- S7 Table and S1 Fig/S31 Table: module designation -----------------
    # S7 is gene-level membership, S8 is the module x set matrix; the sentence
    # now names both rather than citing S8 for content it does not hold.
    ("The full module \u00d7 set score matrix is published (S8 Table) so that "
     "the designation can be recomputed from the deposited data.",
     "The full module \u00d7 set score matrix (S8 Table) and the gene-level "
     "module membership (S7 Table) are published so that the designation can be "
     "recomputed from the deposited data. Hub-gene counts under v1's fixed "
     "|GS| > 0.20 threshold, the stricter |GS| > 0.50, and the FDR < 0.05 "
     "criterion adopted here are given in S1 Fig (counts in S31 Table)."),

    # --- S10 Table: designation sentence ----------------------------------
    ("The module with the highest score was designated the translation module.",
     "The module with the highest score was designated the translation module "
     "(member genes in S10 Table)."),

    # --- S4 Fig: the null-report sentence (not the 20-bin sentence) -------
    ("Both the unmatched (uniform) and matched nulls are reported (S32 Table) "
     "for transparency.",
     "Both the unmatched (uniform) and matched nulls are reported (S32 Table; "
     "S4 Fig) for transparency."),

    # --- S12 Table: ssGSEA parameters -------------------------------------
    ("using an identical parameter object in every cohort.",
     "using an identical parameter object in every cohort (set coverage per "
     "cohort in S12 Table)."),

    # --- S5 Fig: heterogeneity definition in Methods, beside S21 ----------
    ("Between-cohort heterogeneity was quantified by Cochran's *Q* and *I*\u00b2 "
     "[24] within each disease.",
     "Between-cohort heterogeneity was quantified by Cochran's *Q* and *I*\u00b2 "
     "[24] within each disease (S21 Table; S5 Fig)."),

    # --- S33 Table: cross-disease coherence paragraph ---------------------
    ("Both are far above chance; neither is symmetric.",
     "Both are far above chance; neither is symmetric (S33 Table)."),

    # --- S13 Table: the 41 scored sets ------------------------------------
    ("Across all 41 scored sets the cross-disease correlation of effect sizes "
     "is \u03c1 = \u22120.455;",
     "Across all 41 scored sets (S13 Table) the cross-disease correlation of "
     "effect sizes is \u03c1 = \u22120.455;"),

    # --- S34 Table: stage spectrum, Results (not the figure caption) ------
    ("The hepatic translation program is highest in non-malignant chronic liver "
     "disease and falls as malignancy emerges; it is not switched on at the "
     "malignant transition.",
     "The hepatic translation program is highest in non-malignant chronic liver "
     "disease and falls as malignancy emerges; it is not switched on at the "
     "malignant transition (S34 Table)."),

    # --- S2 Fig: score families -------------------------------------------
    ("module eigengene scores are reported alongside it in both directions "
     "(S28, S29 Tables).",
     "module eigengene scores are reported alongside it in both directions "
     "(S28, S29 Tables; S2 Fig)."),
]


def main():
    dry = "--dry-run" in sys.argv
    with open(MD, encoding="utf-8") as fh:
        text = fh.read()

    bad = []
    for old, new in PAIRS:
        n = text.count(old)
        if n != 1:
            bad.append((n, old[:90]))
    if bad:
        print("REFUSING TO WRITE -- %d pair(s) do not match exactly once:" % len(bad))
        for n, head in bad:
            print("  %d match(es): %s..." % (n, head))
        return 1

    for old, new in PAIRS:
        text = text.replace(old, new, 1)

    if dry:
        print("Dry run: all %d pair(s) match exactly once; nothing written." % len(PAIRS))
        return 0

    with open(MD, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("Applied %d pair(s) to %s" % (len(PAIRS), MD))
    return 0


if __name__ == "__main__":
    sys.exit(main())
