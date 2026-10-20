#!/usr/bin/env python
"""
R27_apply_v3_corrections.py -- apply the approved v3 corrections to the manuscript.

Every pair below asserts it matches exactly once. Nothing is written unless all
pairs match, so a sentence that has drifted since the strings were collected
fails loudly instead of silently skipping. This is the same contract as R25/R26.

Provenance of the items (verified against the sources before the strings were
written):

  F1   Fig 1e's panel title asserted "most modules are not shared" on the
       strength of a count that aggregates module-direction tests. 17/43 tests
       show no evidence; 26/43 are at or above the weak-moderate threshold, so
       the direction of the summary was wrong even though the arithmetic was
       right. The caption now states the counts and makes the narrower claim
       about the designated modules only.
  TF   "two-fifths of modules" divides 17 by a denominator of 20 modules
       (2 of 20 show no evidence in both directions), not 43 tests. Replaced
       everywhere with the count over the unit actually computed.
  U1   "the three ... and ... the five" asserts that the two sets are the same
       five; they overlap but are not identical.
  U2   *p*_Q -> *P*_Q, matching the manuscript's P-value notation.
  U3   three coordinated predicates sharing one auxiliary; the "was ... was" form
       in the GSE89377 sentence (which R26 already fixed) is this one's model.
  U4   "soft-threshold" is the WGCNA term for the power parameter; the
       per-dataset choice is of the power, not of the threshold.
  U5   "this data" -> "these data" (two captions).
  N3   the third "against" collides with the two that precede it; the sentence
       compares cells, so "whereas ... gives" states the comparison.
  N4   the exclusion list omitted the post-translational-modification,
       single-gene-miRNA and proliferation (E2F/G2M) families that
       R00_config.R SETS_EXCLUDED actually excludes. The list is completed from
       the config, which is the arbiter.
  U6   Oxford comma at the three named places only; the rest of the class is
       deferred to production copyediting and recorded in the checklist.

Run: python R27_apply_v3_corrections.py [--dry-run]
"""
import sys

MD = "D:/R_projects/revision_analysis/Manuscript_PLOSONE_v2.md"

PAIRS = [
    # --- F1 -----------------------------------------------------------------
    ("showing that most modules are not shared while the translation modules are.",
     "showing that 17 of 43 module\u2013direction tests show no evidence of "
     "preservation, 21 are weak-to-moderate and 5 are strong, and the designated "
     "translation modules are preserved in the direction in which each was tested."),

    # --- TF, abstract -------------------------------------------------------
    ("; about two-fifths of modules show no preservation in either direction.",
     "; about two-fifths of module\u2013direction tests show no evidence of "
     "preservation (17 of 43)."),

    # --- TF, Results --------------------------------------------------------
    ("About two-fifths of the transcriptome's modular structure is not shared "
     "between these diseases; the translation modules are among the minority that is.",
     "About two-fifths of the module\u2013direction preservation tests show no "
     "evidence (17 of 43), and the designated translation modules are preserved "
     "in the direction in which each was tested."),

    # --- TF, Discussion -----------------------------------------------------
    ("about two-fifths of all named modules show no evidence of preservation in "
     "either direction, so this is not a statement that the two transcriptomes "
     "are organised alike.",
     "about two-fifths of the module\u2013direction preservation tests show no "
     "evidence (17 of 43), so this is not a statement that the two "
     "transcriptomes are organised alike."),

    # --- N3 -----------------------------------------------------------------
    ("with GSE116250 and from \u22120.457 to \u22120.059 with GSE141910, against "
     "9 of 14 to 10 of 15 and \u03c1 = \u22120.55 to \u22120.76 in every cell "
     "with GSE57338 (Fig 2d).",
     "with GSE116250 and from \u22120.457 to \u22120.059 with GSE141910, whereas "
     "every cell with GSE57338 gives 9 of 14 to 10 of 15 mirror sets and "
     "\u03c1 = \u22120.55 to \u22120.76 (Fig 2d)."),

    # --- U1 -----------------------------------------------------------------
    ("which are the three Tier A sets that discord between cirrhosis and HCC "
     "and are also the five that fail to mirror in the primary comparison.",
     "which include the three Tier A sets that are discordant between cirrhosis "
     "and HCC and are among the five that fail to mirror in the primary comparison."),

    # --- U2, Results --------------------------------------------------------
    ("median *p*<sub>Q</sub> = 2.1 \u00d7 10\u207b\u2074 and 2.2 \u00d7 "
     "10\u207b\u2076 respectively",
     "median *P*<sub>Q</sub> = 2.1 \u00d7 10\u207b\u2074 and 2.2 \u00d7 "
     "10\u207b\u2076 respectively"),

    # --- U2, S21 caption (same quantity, kept in the same notation) ---------
    ("Cochran's *Q*, *p*<sub>Q</sub>, *I*\u00b2, and pooled effect size",
     "Cochran's *Q*, *P*<sub>Q</sub>, *I*\u00b2, and pooled effect size"),

    # --- U3 -----------------------------------------------------------------
    ("it was essentially unchanged after adjustment for tumour purity (HR 1.52), "
     "only modestly attenuated by proliferation (HR 1.47), and remained "
     "significant with purity, proliferation and immune content together",
     "it remained essentially unchanged after adjustment for tumour purity "
     "(HR 1.52), was only modestly attenuated by proliferation (HR 1.47), and "
     "remained significant with purity, proliferation, and immune content together"),

    # --- U4 -----------------------------------------------------------------
    ("different soft-thresholds chosen by per-dataset optimisation",
     "different soft-thresholding powers chosen per dataset"),

    # --- U5, two captions ---------------------------------------------------
    ("The criteria are nested on this data \u2014 every gene passing |GS| > 0.50 "
     "also passes FDR < 0.05 \u2014",
     "The criteria are nested on these data \u2014 every gene passing |GS| > 0.50 "
     "also passes FDR < 0.05 \u2014"),

    ("The three criteria are nested on this data: every gene passing |GS| > 0.50 "
     "also passes FDR < 0.05.",
     "The three criteria are nested on these data: every gene passing |GS| > 0.50 "
     "also passes FDR < 0.05."),

    # --- N4 -----------------------------------------------------------------
    ("MYC targets, mTORC1, UPR, and all viral and miRNA sets are excluded from "
     "Tiers A and B.",
     "Post-translational modification sets, all viral sets, single-gene "
     "miRNA-regulation sets, and the proliferation-family sets (MYC targets, "
     "mTORC1, E2F targets, G2M checkpoint, and UPR) are excluded from Tiers A "
     "and B."),

    # --- U6, three named places only ---------------------------------------
    ("the union of Tier A, Tier B and Tier C",
     "the union of Tier A, Tier B, and Tier C"),

    ("0.74 (TCGA-LIHC) and 0.54 (GSE89377)",
     "0.74 (TCGA-LIHC), and 0.54 (GSE89377)"),

    ("holding cohort, platform and pipeline constant",
     "holding cohort, platform, and pipeline constant"),
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
