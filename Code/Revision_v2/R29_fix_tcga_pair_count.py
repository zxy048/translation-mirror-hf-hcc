#!/usr/bin/env python
"""
R29_fix_tcga_pair_count.py -- correct the TCGA-LIHC paired-patient count.

V01 audit (2026-10-19), computed from colData of TCGA_LIHC_se.rds:

  424 sample barcodes (371 Primary Tumor + 50 Solid Tissue Normal + 3 Recurrent)
  371 unique patients
   50 unique patients with a Solid Tissue Normal
   50 patients with both a Primary Tumor and a Solid Tissue Normal

The manuscript and the response letter said 53. Fifty normal samples can pair
with at most fifty patients, so 53 was arithmetically impossible; the barcode
count gives 50. (The three recurrent samples come from two patients, which is
the likely origin of the extra three.)

Each pair asserts exactly one match; nothing is written unless all match.

Run: python R29_fix_tcga_pair_count.py [--dry-run]
"""
import sys

FILES = {
    "manuscript": "D:/R_projects/revision_analysis/Manuscript_PLOSONE_v2.md",
    "response": "D:/R_projects/revision_analysis/v2_output/Response_to_Reviewers.md",
}

PAIRS = [
    ("manuscript",
     "50 normal); 53 patients contribute both a tumour and an adjacent-normal sample.",
     "50 normal); 50 patients contribute both a tumour and an adjacent-normal sample."),
    ("response",
     "53 patients in the analysed set contribute both a primary tumour and a solid-tissue-normal sample.",
     "50 patients in the analysed set contribute both a primary tumour and a solid-tissue-normal sample."),
]


def main():
    dry = "--dry-run" in sys.argv
    texts, bad = {}, []
    for name, path in FILES.items():
        with open(path, encoding="utf-8") as fh:
            texts[name] = fh.read()

    for name, old, new in PAIRS:
        n = texts[name].count(old)
        if n != 1:
            bad.append((name, n, old[:80]))
    if bad:
        print("REFUSING TO WRITE -- pair(s) do not match exactly once:")
        for name, n, head in bad:
            print("  %s: %d match(es): %s..." % (name, n, head))
        return 1

    for name, old, new in PAIRS:
        texts[name] = texts[name].replace(old, new, 1)

    if dry:
        print("Dry run: all %d pair(s) match exactly once; nothing written." % len(PAIRS))
        return 0

    for name, path in FILES.items():
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(texts[name])
    print("Applied %d pair(s): %s" % (len(PAIRS), ", ".join(sorted(FILES))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
