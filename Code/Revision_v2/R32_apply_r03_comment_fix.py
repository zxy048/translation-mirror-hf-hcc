# -*- coding: utf-8 -*-
"""Correct the stale header comment in R03_translation_module.R.

The executed designation rule is hypergeometric fold enrichment among modules
at FDR < alpha (primary), with eigengene tracking as the secondary check. The
header block still described the old order (tracking primary, hypergeometric
secondary). This script replaces that block only; the R code is not touched.

Guards:
  * the old block must occur exactly once;
  * after removing every full-line comment, the remaining text must be
    identical before and after, so no executable statement, literal, or
    parameter can change;
  * R00_config.R and the S1 table are hashed before and after and must not
    change; R03 does not write either, but the check makes that explicit.
"""
import hashlib
import os
import sys

ROOT = r"D:/R_projects/revision_analysis"
R03 = os.path.join(ROOT, "Code", "Revision_v2", "R03_translation_module.R")
R00 = os.path.join(ROOT, "Code", "Revision_v2", "R00_config.R")
S1 = os.path.join(ROOT, "v2_output", "tables", "Table_S1_gene_universe_loss_chain.csv")

OLD = """#   1. ssGSEA score per sample for every canonical translation set
#   2. correlate each set score with each module eigengene
#   3. translation module = argmax over modules of mean |r| across Tier A sets,
#      with a sample-level permutation p-value (max-statistic null, so the
#      p-value is corrected across the whole sets x modules matrix)
#   4. independent secondary check: hypergeometric over-representation of
#      Tier A genes inside each module
"""

NEW = """#   1. ssGSEA score per sample for every canonical translation set
#   2. hypergeometric over-representation of Tier A genes inside each module
#      (one-sided; fold enrichment)
#   3. translation module = the module with the largest fold enrichment among
#      modules significant at FDR < alpha; the identical rule is applied to
#      both diseases
#   4. secondary check: eigengene tracking, i.e. mean |r| between each module
#      eigengene and the Tier A scores, with a sample-level permutation
#      p-value (max-statistic null, corrected across the sets x modules matrix)
"""


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    with open(R03, "r", encoding="utf-8", newline="") as fh:
        before = fh.read()

    n_old = before.count(OLD)
    n_new = before.count(NEW)
    if n_old == 0 and n_new == 1:
        print("already applied: new comment block present once, old block absent")
        print("no change made")
        return
    if n_old != 1:
        print(f"ABORT: old comment block occurs {n_old} times, expected 1")
        sys.exit(1)

    after = before.replace(OLD, NEW, 1)

    # R is not Python: prove comment-only by removing every full-line comment
    # and requiring the remaining executable text to be identical.
    code_before = "\n".join(l for l in before.splitlines() if not l.lstrip().startswith("#"))
    code_after = "\n".join(l for l in after.splitlines() if not l.lstrip().startswith("#"))
    if code_before != code_after:
        print("ABORT: executable code changed; refusing to write")
        sys.exit(1)

    hashes_before = {"R00": sha(R00), "S1": sha(S1)}
    hashes_after = {"R00": sha(R00), "S1": sha(S1)}
    if hashes_before != hashes_after:
        print("ABORT: R00_config.R or S1 table changed unexpectedly")
        sys.exit(1)

    # The only diff may be one replacement inside the header block; line 8
    # (the unchanged first item) acts as the anchor, so the replacement spans
    # 0-based lines 8-14 -> 8-16. Every inserted line must itself be a comment.
    import difflib
    old_lines = before.splitlines(keepends=True)
    new_lines = after.splitlines(keepends=True)
    ops = [op for op in difflib.SequenceMatcher(None, old_lines, new_lines).get_opcodes()
           if op[0] != "equal"]
    if ops != [("replace", 8, 14, 8, 16)]:
        print(f"ABORT: unexpected line diff {ops}")
        sys.exit(1)
    if any(not l.lstrip().startswith("#") for l in new_lines[8:16]):
        print("ABORT: a non-comment line entered the replacement block")
        sys.exit(1)

    with open(R03, "w", encoding="utf-8", newline="") as fh:
        fh.write(after)

    print(f"R03 sha256 {hashlib.sha256(before.encode()).hexdigest()} -> "
          f"{hashlib.sha256(after.encode()).hexdigest()}")
    print(f"diff op: {ops}")
    print(f"R00_config.R sha256 {hashes_after['R00']} (unchanged)")
    print(f"S1 table sha256 {hashes_after['S1']} (unchanged)")
    print("OK: comment-only change; executable text identical")


if __name__ == "__main__":
    main()
