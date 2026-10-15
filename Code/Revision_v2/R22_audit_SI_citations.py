# -*- coding: utf-8 -*-
"""
Audit every Supporting Information item against the manuscript body text.

WHY THIS EXISTS
---------------
A reviewer-side audit asserted that some SI tables -- S20b named specifically --
are never called out in the text, and that others are cited inconsistently.
That is a format claim, so it is checkable, and a claim about citations checked
by eye is exactly the kind of claim this repository has been burned by before.

The manuscript cites SI items with three shorthands, all of which have to be
expanded before a "never cited" verdict means anything:

    S19, S20 Tables        an explicit list
    S23-S26 Tables         a numeric range
    S22b-d Tables          a letter range that repeats no number

The caption definitions themselves are excluded from the body, so a table that
only mentions another table in its own caption does not count as cited. That
exclusion covers the S1 File caption, which names every table and figure in the
supplement and would otherwise vouch for all of them (see CAPTION below).

S1 File is the supplement as a whole rather than an item inside it, so it is
defined but never tallied: a body sentence citing the container says nothing
about whether a reader can find any particular table in the text.

Read-only. Writes nothing; prints a report and exits non-zero if any item is
never cited, so it can be wired into verify_all.py later if wanted.

Run:  python Code/Revision_v2/R22_audit_SI_citations.py
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MS = os.path.join(ROOT, "Manuscript_PLOSONE_v2.md")

# "File" is matched as well as "Table" and "Fig", and it has to be. The S1 File
# caption reads "A single PDF containing Tables S1-S39 and Figures S1-S5", so
# while it went unrecognised it was left in the body -- where the numeric-range
# expander below read "S1-S39" and marked *every table* as cited. The audit then
# reported 39/39 on any manuscript at all, which is the one failure mode an
# audit must not have. Treating the line as the caption definition it is takes
# it out of the body and restores the count.
CAPTION = re.compile(r"^\*\*(S\d+[a-z]?) (Table|Fig|File)\.\*\*")
LETTERS = "abcdefgh"


def key(token):
    """'S22b' -> (22, 'b'); 'S7' -> (7, '')."""
    m = re.fullmatch(r"S(\d+)([a-z]?)", token)
    return int(m.group(1)), m.group(2)


def expand(text):
    """Every SI id the body text cites, with the three shorthands expanded."""
    cited = set()
    for token in re.findall(r"\bS\d+[a-z]?\b", text):
        cited.add(token)

    # S23-S26, S22b-S22d, and S22b-d, which repeats neither the number nor the S
    for a, ad, b, bd in re.findall(
            r"\bS(\d+)([a-z]?)\s*[–—-]\s*S?(\d*)([a-z]?)", text):
        if not a:
            continue
        na, la = int(a), ad
        if b:                             # S23-S26: a numeric range
            for n in range(na, int(b) + 1):
                cited.add("S%d" % n)
        elif bd:                          # S22b-d: a letter range, same number
            for c in LETTERS[LETTERS.index(la):LETTERS.index(bd) + 1]:
                cited.add("S%d%s" % (na, c))
    return cited


def main():
    with open(MS, encoding="utf-8") as fh:
        lines = fh.read().split("\n")

    defined, body = [], []
    for ln in lines:
        m = CAPTION.match(ln)
        if m:
            defined.append((m.group(1), m.group(2)))
        else:
            body.append(ln)

    cited = expand("\n".join(body))

    print("manuscript: %s" % MS)
    # Counted per kind and summed, rather than printing len(defined) beside the
    # table and figure counts: once "File" was recognised the three no longer
    # added up (45 items printed as "39 tables, 5 figures") and a breakdown that
    # does not total is worse than no breakdown.
    print("SI items defined: %s" % ", ".join(
        "%d %s" % (sum(1 for _, k in defined if k == kind),
                   kind if sum(1 for _, k in defined if k == kind) == 1 else kind + "s")
        for kind in ("Table", "Fig", "File")
        if any(k == kind for _, k in defined)))
    print()

    # The other direction, and the one that actually breaks a submission: a body
    # citation with no matching caption definition is a dangling reference. The
    # reader is sent to an item that does not exist in the supplement, and PLOS
    # renders it as a broken link. An uncited item is merely unhelpful; a
    # dangling one is wrong, so this direction fails hard while the other stays
    # advisory.
    defined_ids = set(d for d, _ in defined)
    dangling = sorted(cited - defined_ids, key=key)
    if dangling:
        print("DANGLING: cited in the body but never defined: %s"
              % ", ".join(dangling))
        print()

    uncited = 0
    for kind in ("Table", "Fig"):
        ids = [d for d, k in defined if k == kind]
        missing = [d for d in ids if d not in cited]
        uncited += len(missing)
        print("%ss cited in the body: %d/%d" % (kind, len(ids) - len(missing), len(ids)))
        if missing:
            print("  never cited: %s" % ", ".join(missing))
        print()

    if dangling:
        print("FAILED: %d dangling reference(s). A body citation with no caption"
              % len(dangling))
        print("definition points the reader at an item the supplement does not")
        print("contain, so this must be fixed before submission.")
        return 1

    if uncited:
        # Advisory, not a compliance failure. PLOS: "We recommend that you cite
        # supporting information in the manuscript text, but this is not a
        # requirement." The list is what a human needs in order to decide.
        print("%d SI item(s) are not cited in the manuscript body." % uncited)
        print("PLOS recommends citation but does not require it, so this is not a")
        print("compliance failure. It is listed because v1's body did cite its SI")
        print("figures and the v2 rewrite did not carry those citations over.")
        print("No dangling references.")
        return 1
    print("PASSED: every SI item is cited in the manuscript body, and no citation")
    print("in the body lacks a caption definition.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
