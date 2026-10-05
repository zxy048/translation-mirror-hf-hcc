# -*- coding: utf-8 -*-
"""
Refuse to publish peer-review material to the public repository.

WHY THIS FILE EXISTS
--------------------
This repository is public (github.com/zxy048/translation-mirror-hf-hcc). It is
also the working directory, so the four journal submission packages, the
manuscript drafts and the cover letters all sit in the same tree as the
analysis code. `git add -A` therefore sweeps in whatever happens to be lying
around -- and what lies around after a revision round includes the editor's
decision letter and all four reviewers' comments.

Reviewer comments are confidential. Publishing them is a breach of the review
process, of the journal's terms, and of the reviewers' expectation -- and it is
irreversible once pushed, because the content is indexed and forked within
minutes. This is not a hypothetical: `PLOS_ONE_Submission/00_Decision_Letter_
Reviewer_Comments.txt` is currently untracked and NOT ignored, so a single
`git add -A && git push` would publish it.

Two layers, because one is not enough. `.gitignore` stops the sweep, but it
only covers paths someone thought to list, and it silently stops working the
moment a file is renamed into a directory that is not covered. This script
reads what git would actually stage and inspects the content, so it catches the
file that .gitignore does not know about yet.

Run before every push. Exits non-zero if anything confidential would be
published.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))

# Filename patterns that must never be published. Peer-review correspondence
# and reviewer comments are the hard blocks; manuscript drafts are flagged for
# a human decision rather than blocked, since the repo already carries some.
BLOCK_NAME = [
    (r"decision[_ -]?letter", "editor's decision letter"),
    (r"reviewer[_ -]?comment", "reviewers' comments"),
    (r"\breviewer\b.*\.(txt|md|pdf|docx)$", "reviewer correspondence"),
    (r"response[_ -]?to[_ -]?reviewers", "response letter (quotes reviewer comments)"),
    (r"rebuttal", "rebuttal (quotes reviewer comments)"),
]

# Content markers, searched in text files that would be staged. A file whose
# name gives nothing away can still contain the review -- which is how this
# check earns its place: the first run caught four session-transcript files
# under v2_output/logs/ that quote reviewers' comments in full and are named
# in a way that reveals nothing.
BLOCK_CONTENT = [
    (r"Reviewer\s*#?\s*\d\s*[:.]", "reviewer comment block"),
    (r"Decision Letter", "decision letter heading"),
    (r"confidential", "explicitly marked confidential"),
]

# The checker's own source names these patterns, so it always matches itself.
SELF = "pre_publish_check.py"

# A file may discuss this policy without containing the material it protects --
# the README's disclosure section quotes the phrase "decision letter" while
# explaining why the letter itself is excluded. Such a file opts out by carrying
# this sentinel, so the exemption is explicit and visible in the diff rather
# than an implicit hole in the pattern. This guards against accident, not
# against an author who deliberately labels a leak as policy.
POLICY_SENTINEL = "pre-publish-check: policy-documentation"

TEXT_EXT = {".txt", ".md", ".py", ".R", ".csv", ".json", ".tex", ".html"}
MAX_SCAN = 4 * 1024 * 1024  # skip anything larger; review material is not


def git(*args):
    p = subprocess.run(["git"] + list(args), cwd=ROOT,
                       capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    return p.stdout


def would_publish():
    """Every path git would stage, tracked or not, deleted files excluded."""
    out = git("status", "--porcelain", "--untracked-files=all")
    paths = []
    for line in out.splitlines():
        if len(line) < 4:
            continue
        status, path = line[:2], line[3:].strip()
        if status.strip() == "D":
            continue
        # git quotes paths containing non-ASCII or spaces
        if path.startswith('"') and path.endswith('"'):
            path = path[1:-1].encode("utf-8").decode("unicode_escape")
        paths.append((status.strip(), path))
    return paths


def scan_content(path):
    ext = os.path.splitext(path)[1].lower()
    if ext not in TEXT_EXT:
        return []
    full = os.path.join(ROOT, path)
    try:
        if os.path.getsize(full) > MAX_SCAN:
            return []
        with open(full, encoding="utf-8", errors="ignore") as fh:
            text = fh.read()
    except OSError:
        return []
    if POLICY_SENTINEL in text:
        return []
    hits = []
    for pattern, what in BLOCK_CONTENT:
        m = re.search(pattern, text, re.IGNORECASE)
        if m:
            hits.append("%s (matched %r)" % (what, m.group(0)[:40]))
    return hits


def main():
    entries = would_publish()
    if not entries:
        print("Nothing staged or untracked. Nothing to publish.")
        return 0

    blocked, flagged = [], []
    for status, path in entries:
        if os.path.basename(path) == SELF:
            continue
        for pattern, what in BLOCK_NAME:
            if re.search(pattern, path, re.IGNORECASE):
                blocked.append((path, what))
                break
        else:
            for hit in scan_content(path):
                blocked.append((path, hit))
                break
        # Flag rather than block: drafts and cover letters are already part of
        # this repository's history, so the decision is the author's, but it
        # should be a decision and not an accident.
        low = path.lower()
        if any(k in low for k in ("manuscript_", "cover_letter", "conflict_of_interest",
                                  "declarations", "submission/", "_submission")):
            flagged.append(path)

    print("=" * 74)
    print("What a push would publish: %d path(s)" % len(entries))
    print("=" * 74)
    for status, path in entries:
        print("  %-2s %s" % (status, path))

    if flagged and not blocked:
        print("\n" + "-" * 74)
        print("Manuscript and submission material also in scope (%d). These are" % len(flagged))
        print("not blocked -- the repository already carries such files -- but they")
        print("will be public. Confirm each is intended:")
        for path in flagged:
            print("    ? %s" % path)

    print("\n" + "=" * 74)
    if blocked:
        print("REFUSING TO PUBLISH -- %d confidential item(s) in scope:" % len(blocked))
        for path, why in blocked:
            print("  * %s" % path)
            print("      %s" % why)
        print("\nPeer-review material must not be pushed to a public repository.")
        print("Add these paths to .gitignore, or move them outside the repository,")
        print("then re-run this check.")
        return 1
    print("No confidential material in scope.")
    print("Review the list above before pushing; a clean check is not the same as")
    print("a reviewed one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
