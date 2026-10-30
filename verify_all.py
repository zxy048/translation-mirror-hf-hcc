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
import tempfile

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

# (label, script path relative to ROOT, what it reads that the manuscript does
# not, file the check cannot run without)
#
# The fourth field exists for exactly one check. The response letter is
# peer-review correspondence and is gitignored, so it is in the working tree
# and in no public repository. A check that needs it must be *skipped* there,
# not failed: a red "response letter numbers" line reads as "the manuscript's
# numbers do not reconcile", which is the opposite of what has happened. The
# skip prints why, so the missing check cannot pass for a passing one.
CHECKS = [
    ("manuscript numbers vs regenerated tables",
     "Code/Revision_v2/R18_verify_manuscript.py",
     "opens the regenerated S-table CSVs and re-derives every number quoted in the text",
     None),
    ("supporting information",
     "Code/Revision_v2/R20_verify_SI.py",
     "opens the built SI PDF and checks every table and figure it holds against its source CSV",
     None),
    ("response letter numbers",
     "Code/Revision_v2/R24_audit_response_letter.py",
     "opens the response letter and checks its numbers against the manuscript, the "
     "regenerated tables and the built figure files",
     os.path.join("v2_output", "Response_to_Reviewers_condensed_draft.md")),
    ("manuscript structure and PLOS formatting",
     "qa_plos_v2.py",
     "checks abstract length, citation order, section presence, page setup, line numbers",
     None),
    ("tracked-changes docx round-trip",
     "build_tracked_docx.py",
     "rebuilds the tracked docx and confirms accept-all reproduces v2 and reject-all reproduces v1",
     None),
]

# v1 figures were rendered at 13 x 11.5 in and downscaled by the journal to
# 7.5 x 8.75 in, which shrank every label with them. The ceiling is asserted
# here rather than trusted, because "we rendered it at the right size" is
# exactly the kind of claim that silently goes stale after a re-render.
#
# The figures are TIFF, not PDF: PLOS takes TIFF or EPS, and rasterising the
# vector page is also what removes Fig 2's 6 dpi colour-bar bitmap. The checked
# file is therefore the uploaded one. Reading the TIFF's own tag instead of
# recomputing the render is deliberate -- a check that re-derives the number it
# is checking cannot notice that the render and the check disagree.
PLOS_MAX_W_IN = 7.5
PLOS_MAX_H_IN = 8.75
PLOS_MIN_DPI = 300
PLOS_MAX_MB = 10.0
FIGURES = ["PLOS_ONE_Submission_v2/Fig1.tiff",
           "PLOS_ONE_Submission_v2/Fig2.tiff",
           "PLOS_ONE_Submission_v2/Fig3.tiff"]


def tiff_specs(path):
    """(width_px, height_px, dpi_x, dpi_y) as the TIFF file itself records them."""
    from PIL import Image
    with Image.open(path) as im:
        w, h = im.size
        dpi = im.info.get("dpi") or (0, 0)
    return w, h, dpi[0], dpi[1]


def check_figure_specs():
    """Every main figure must be a TIFF at or above PLOS's resolution floor and
    at or under its print-size ceiling."""
    fails = []
    for rel in FIGURES:
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            fails.append("%s is missing" % rel)
            continue
        try:
            w, h, dx, dy = tiff_specs(path)
        except Exception as exc:                      # unreadable or not a TIFF
            fails.append("%s: not a readable TIFF (%s)" % (rel, exc))
            continue
        if dx < PLOS_MIN_DPI or dy < PLOS_MIN_DPI:
            fails.append("%s is %g dpi, under the %d dpi floor"
                         % (rel, min(dx, dy), PLOS_MIN_DPI))
            continue
        w_in, h_in = w / dx, h / dy
        mb = os.path.getsize(path) / 1e6
        ok = (w_in <= PLOS_MAX_W_IN + 0.01 and h_in <= PLOS_MAX_H_IN + 0.01
              and mb <= PLOS_MAX_MB)
        print("  %s %s  %.2f x %.2f in (%d x %d px @ %g dpi, %.2f MB)"
              % ("PASS" if ok else "FAIL", rel, w_in, h_in, w, h, dx, mb))
        if not ok:
            fails.append("%s is %.2f x %.2f in / %.2f MB, over the "
                         "%.2f x %.2f in / %.0f MB ceiling"
                         % (rel, w_in, h_in, mb, PLOS_MAX_W_IN, PLOS_MAX_H_IN,
                            PLOS_MAX_MB))
    return fails


# ---------------------------------------------------------------------------
# Every Supporting Information item the manuscript defines must also be cited.
#
# PLOS treats an uncited supporting file as an incomplete submission, and a
# caption plus no citation is invisible in every other check: the supplement
# builds, the caption table is complete, the page count is plausible. Three
# rounds of manual audit left ten or so items uncited because nothing made the
# omission fail. This makes it a build error.
#
# The citation forms are parsed rather than searched for literally, because the
# manuscript writes "S4, S5 Tables" and "S23–S26 Tables" as often as "S4 Table".
# A literal search reported 17 uncited items; the real number was 12.
#
# A reference counts anywhere in the manuscript except the item's own caption
# and the consolidated S1 File index, which names every item and would otherwise
# mark the whole supplement as cited. Caption-to-caption references do count:
# a reader following either one reaches the same table.
# ---------------------------------------------------------------------------
SI_CAPTION = re.compile(r"^\*\*(S\d+[a-z]?) (Table|Fig)\.\*\*")
SI_INDEX = re.compile(r"^\*\*S1 File\.\*\*")


def check_si_citations():
    path = os.path.join(ROOT, "Manuscript_PLOSONE_v2.md")
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()

    caps = {}
    for i, ln in enumerate(lines):
        m = SI_CAPTION.match(ln)
        if m:
            caps.setdefault((m.group(1), m.group(2)), i)

    # The index names every item, so it is not a citation. Every other line --
    # including another item's caption -- is: S31 is reachable only from the
    # S1 Fig caption and S33 only from S36's, and a reader following either
    # reaches the table. Only the item's *own* caption is excluded.
    index = {i for i, ln in enumerate(lines) if SI_INDEX.match(ln)}

    refs = []                        # (line index, kind, number, letter suffix)
    for i, ln in enumerate(lines):
        if i in index:
            continue
        for m in re.finditer(
                r"((?:S\d+[a-z]?)(?:\s*[,–-]\s*S\d+[a-z]?)*)\s+Tables?", ln):
            for part in re.split(r"\s*,\s*", m.group(1)):
                r = re.match(r"S(\d+)([a-z]?)(?:[–-]S(\d+)([a-z]?))?", part)
                if r:
                    a = int(r.group(1))
                    b = int(r.group(3)) if r.group(3) else a
                    if a == b:
                        # S18b is its own item, not a synonym for S18: carry
                        # the letter through so a plain "S18 Table" cannot
                        # stand in for a caption that is never cited.
                        refs.append((i, "Table", a, r.group(2) or r.group(4) or ""))
                    else:
                        refs.extend((i, "Table", n, "") for n in range(a, b + 1))
        for m in re.finditer(r"S(\d+)([a-z]?)\s*Fig", ln):
            refs.append((i, "Fig", int(m.group(1)), m.group(2)))
        for m in re.finditer(r"Fig(?:ure)?s?\.?\s*[–-]?\s*S(\d+)([a-z]?)", ln):
            refs.append((i, "Fig", int(m.group(1)), m.group(2)))

    def label_key(lab):
        m = re.match(r"S(\d+)([a-z]?)", lab)
        return int(m.group(1)), m.group(2)

    fails = []
    for (label, kind), own in sorted(caps.items(),
                                     key=lambda kv: (kv[0][1],) + label_key(kv[0][0])):
        n, suf = label_key(label)
        if any(k == kind and num == n and s == suf and i != own
               for i, k, num, s in refs):
            continue
        fails.append("%s %s (caption L%d) is defined but never cited"
                     % (label, kind, own + 1))

    print("  %d caption(s), %d cited, %d uncited"
          % (len(caps), len(caps) - len(fails), len(fails)))
    for f in fails:
        print("  FAIL " + f)
    return fails


# ---------------------------------------------------------------------------
# Cross-document consistency.
#
# Each validator above checks one document against the pipeline. A quantity
# that is right in one document and stale in another therefore passes every
# one of them: R24 could not see the response letter's 2.9% because 2.9 is
# arithmetically derivable from the Tier A subset, while the manuscript
# reported the all-fits 2.8%; and the manuscript, the cover letter and the
# response letter each carried their own R00-Rnn range. Those are the seams
# between documents, and no per-document validator can see across them.
#
# Three shared quantities are checked here: the median signed attenuation after
# age/sex adjustment (S37), the median signed attenuation after tissue-
# composition adjustment (S25), and the pipeline range cited to the public
# repository. All three are values a reader will compare across documents.
# ---------------------------------------------------------------------------
# Two *different* adjustments are both called "attenuation", and both are
# reported in the manuscript and the response letter:
#
#   age and sex (S37)      -- 163 set-cohort fits, median signed attenuation 2.8%
#   tissue composition (S25) -- 103 Tier A pairs,  median signed attenuation 34.1%
#
# They are separate quantities, so each is compared only with itself: a check
# that pooled them would fail on a correct document (the value would be "not
# unique") and a check that silently took the first would pass one of them
# unexamined. The scope is read from the sentence the value sits in, because the
# phrase alone does not distinguish them and the two documents order the two
# clues differently -- the manuscript writes "2.8% across the 163 fits" and
# "over the 103 pairs ... 34.1%".
#
# The wording variants are accepted rather than normalised: the manuscript has
# "Median signed attenuation is", the condensed response letter has "median
# attenuation is" and "median signed attenuation is", and the long version had
# "Median attenuation of |g| is".
#
# The two letters are the files assemble_plos_v2.py actually builds 01 and 02
# from. They were Cover_Letter_PLOSONE_v2.md and Response_to_Reviewers.md until
# 2026-10-27; both were superseded that day, and a cross-document check still
# reading the superseded pair would report agreement between the manuscript and
# documents nobody is submitting.
ATTEN = re.compile(
    r"[Mm]edian (?:signed )?attenuation (?:of \|\*g\*\| )?(?:is|of) ([\d.]+)%")
ATTEN_SCOPE = (
    ("age/sex (S37, 163 fits)", re.compile(r"\b163\b")),
    ("composition (S25, 103 pairs)", re.compile(r"\b103\b")),
)
# Documents that must carry the figures if they are present at all.
ATTEN_REQUIRED = ("manuscript", "response letter")
PIPE = re.compile(r"R00\s*[–-]\s*R(\d+)")
DOCS = [("manuscript", "Manuscript_PLOSONE_v2.md"),
        ("cover letter", "Cover_Letter_PLOSONE_R34.md"),
        ("response letter", os.path.join("v2_output",
                                         "Response_to_Reviewers_condensed_draft.md"))]


def attenuation_values(text):
    """(scope, value) for every attenuation figure in `text`.

    The scope is whichever of the two fit counts sits nearest the figure in its
    own sentence, measured from the end of the matched phrase (just past the
    percentage). An earlier version took the first scope marker found in the
    sentence, which assigned both figures to age/sex whenever a document
    reported them in one sentence; measuring from the start of the phrase was
    no better, because the phrase begins with "median", which pulls the match
    back towards whatever precedes it.
    """
    out = []
    for sent in re.split(r"(?<=[.)])\s+", text):
        for m in ATTEN.finditer(sent):
            best, best_d = None, None
            for name, rx in ATTEN_SCOPE:
                for a in rx.finditer(sent):
                    d = abs(a.start() - m.end())
                    if best_d is None or d < best_d:
                        best, best_d = name, d
            out.append((best or "unspecified", m.group(1)))
    return out


def check_cross_document():
    fails, texts = [], {}
    for name, rel in DOCS:
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            # The letters are gitignored, so a clone of the public repository
            # has none of them. That is the expected state there, not a defect:
            # the checks still run on whatever is present, and the skip says so
            # rather than passing silently.
            print("  SKIP  %s (%s) is not present." % (name, rel))
            print("        Peer-review correspondence is never published, so a")
            print("        clone of the repository legitimately lacks it.")
            continue
        with open(path, encoding="utf-8") as fh:
            texts[name] = fh.read()

    # 1. each attenuation quantity, one value per document and the same value
    #    in every document that reports it.
    atten = {}
    for name, text in texts.items():
        by_scope = {}
        for scope, v in attenuation_values(text):
            by_scope.setdefault(scope, set()).add(v)
        for scope in sorted(by_scope):
            vals = by_scope[scope]
            if len(vals) != 1:
                fails.append("%s: %s attenuation is not unique (%s)"
                             % (name, scope, ", ".join(sorted(vals))))
                continue
            atten.setdefault(scope, {})[name] = next(iter(vals))
        if name in ATTEN_REQUIRED:
            for scope, _ in ATTEN_SCOPE:
                if scope not in by_scope:
                    fails.append("%s: no %s attenuation figure found"
                                 % (name, scope))
    for scope in sorted(atten):
        got = atten[scope]
        print("  attenuation, %s: %s" % (scope, ", ".join(
            "%s=%s%%" % (n, got[n]) for n in sorted(got))))
        if len(set(got.values())) > 1:
            fails.append("%s attenuation disagrees across documents: %s"
                         % (scope, ", ".join("%s=%s%%" % (n, v)
                                             for n, v in sorted(got.items()))))
    if not texts:
        print("  nothing to compare: no cross-document text is present")

    # 2. one pipeline range, not one per document.
    #
    # The manuscript states it (under Data Availability). The cover letter does
    # not, and no longer needs to: it cites the repository rather than a script
    # range. So a document that states no range is silent, not wrong, and only
    # the manuscript is required to state one. A document that states two
    # different ranges is wrong.
    pipe = {}
    for name, text in texts.items():
        vals = PIPE.findall(text)
        if not vals:
            if name == "manuscript":
                fails.append("%s: pipeline range not found" % name)
            continue
        if len(set(vals)) != 1:
            fails.append("%s: pipeline range is not unique (%s)"
                         % (name, ", ".join(sorted(set(vals)))))
            continue
        pipe[name] = "R00–R" + vals[0]
    if pipe:
        print("  pipeline range: " + ", ".join(
            "%s=%s" % (n, pipe[n]) for n in sorted(pipe)))
    if len(set(pipe.values())) > 1:
        fails.append("pipeline range disagrees across documents: %s"
                     % ", ".join("%s=%s" % (n, v) for n, v in sorted(pipe.items())))

    for f in fails:
        print("  FAIL " + f)
    return fails


# ---------------------------------------------------------------------------
# Scoped cohort/module counts that the text quotes and the tables define.
#
# Each number below was audited by hand and then written into the manuscript,
# the response letter or the cover letter. Nothing else in this file re-derives
# them, so a regenerated table could change a count while every other check
# stayed green -- exactly the failure mode that produced the v1 submission.
# These assertions pin the counts to the table that defines them, with the
# scope stated explicitly:
#
#   S14   268 rows over the 7 cohorts that have a control arm. GSE141198 is
#         tumour-only (n_control = 0) and is not a disease-control cohort.
#   S40a  1,965 sample rows over all 8 harmonised cohorts. This is the
#         all-cohort scope, not the paired scope, and GSE89377 is not in it.
#   S40b  334 complete patient pairs = 668 samples. Paired analysis is a
#         sensitivity analysis; these are not the primary independent samples.
#   S4    23 named HF modules and 20 named HCC modules; 24 and 21 including
#         the unassigned grey bin.
#   Fig 1 the schematic's scoring box describes the disease-control cohorts,
#         not a hard-coded cohort count, and the caption cites Figs 2-3.
# ---------------------------------------------------------------------------
def check_scoped_counts():
    import csv

    def read_csv(rel):
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            return None
        with open(path, encoding="utf-8-sig", newline="") as fh:
            return list(csv.DictReader(fh))

    fails = []

    # 1. S14: disease-control contrasts only.
    s14 = read_csv(os.path.join("v2_output", "tables",
                                "Table_S14_effect_sizes_within_disease.csv"))
    if s14 is None:
        fails.append("S14 table is missing")
    else:
        cohorts = {r["cohort"] for r in s14}
        no_control = sorted({r["cohort"] for r in s14
                             if float(r["n_control"]) == 0})
        print("  S14: %d rows, %d cohorts" % (len(s14), len(cohorts)))
        if len(s14) != 268:
            fails.append("S14 has %d rows, expected 268" % len(s14))
        if len(cohorts) != 7:
            fails.append("S14 has %d cohorts, expected 7" % len(cohorts))
        if no_control:
            fails.append("S14 contains cohort(s) without controls: %s"
                         % ", ".join(no_control))
        if "GSE141198" in cohorts:
            fails.append("S14 contains GSE141198 (tumour-only, no control arm)")

    # 2. S40a: all-cohort paired manifest, not the paired-only subset.
    s40a = read_csv(os.path.join("v2_output", "tables",
                                 "Table_S40a_paired_sample_manifest.csv"))
    if s40a is None:
        fails.append("S40a table is missing")
    else:
        cohorts = {r["cohort"] for r in s40a}
        print("  S40a: %d sample rows, %d cohorts" % (len(s40a), len(cohorts)))
        if len(s40a) != 1965:
            fails.append("S40a has %d sample rows, expected 1,965" % len(s40a))
        if len(cohorts) != 8:
            fails.append("S40a has %d cohorts, expected 8" % len(cohorts))
        if "GSE141198" not in cohorts:
            fails.append("S40a is missing GSE141198 (tumour-only cohort)")
        if "GSE89377" in cohorts:
            fails.append("S40a contains GSE89377, which is outside its scope")

    # 3. S40b: complete pairs and their samples.
    s40b = read_csv(os.path.join("v2_output", "tables",
                                 "Table_S40b_pairing_summary.csv"))
    if s40b is None:
        fails.append("S40b table is missing")
    else:
        pairs = sum(int(r["n_complete_pairs"]) for r in s40b)
        samples = sum(int(r["n_samples_in_complete_pairs"]) for r in s40b)
        print("  S40b: %d complete pairs, %d paired samples" % (pairs, samples))
        if pairs != 334:
            fails.append("S40b has %d complete pairs, expected 334" % pairs)
        if samples != 668:
            fails.append("S40b has %d paired samples, expected 668" % samples)
        if samples != 2 * pairs:
            fails.append("S40b paired samples (%d) are not 2 x pairs (%d)"
                         % (samples, pairs))

    # 4. S4: named modules exclude the grey bin.
    s4 = read_csv(os.path.join("v2_output", "tables",
                               "Table_S4_module_inventory.csv"))
    if s4 is None:
        fails.append("S4 table is missing")
    else:
        for disease, named, total in (("HF", 23, 24), ("HCC", 20, 21)):
            mods = [r["module"] for r in s4 if r["disease"] == disease]
            n_named = sum(1 for m in mods if m.lower() != "grey")
            print("  S4 %s: %d modules total, %d named" % (disease, len(mods), n_named))
            if len(mods) != total:
                fails.append("S4 %s has %d modules, expected %d"
                             % (disease, len(mods), total))
            if n_named != named:
                fails.append("S4 %s has %d named modules, expected %d"
                             % (disease, n_named, named))

    # 5. Fig 1: caption and schematic agree with the disease-control scope.
    ms_path = os.path.join(ROOT, "Manuscript_PLOSONE_v2.md")
    r17_path = os.path.join(ROOT, "Code", "Revision_v2", "R17_figures.R")
    caption = ""
    if not os.path.exists(ms_path):
        fails.append("manuscript is missing")
    else:
        with open(ms_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("**Fig 1.**"):
                    caption = line
                    break
        if not caption:
            fails.append("Fig 1 caption not found")
        else:
            for need in ("each disease–control cohort", "Figs 2–3"):
                if need not in caption:
                    fails.append("Fig 1 caption lacks %r" % need)
            for bad in ("six cohorts", "Figs 2–4"):
                if bad in caption:
                    fails.append("Fig 1 caption still contains %r" % bad)
            print("  Fig 1 caption: disease-control scope, Figs 2-3")
    if not os.path.exists(r17_path):
        fails.append("R17_figures.R is missing")
    else:
        with open(r17_path, encoding="utf-8") as fh:
            r17 = fh.read()
        if "per two-arm cohort" not in r17:
            fails.append("R17 lacks 'per two-arm cohort'")
        if "six cohorts" in r17:
            fails.append("R17 still contains 'six cohorts'")
        # R17 stores the condition as the escape ≠ inside an R string.
        if "CI\\u2260 0" not in r17 and "CI≠ 0" not in r17:
            fails.append("R17 mirror box lacks the CI-excludes-zero condition")

    for f in fails:
        print("  FAIL " + f)
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
    # Output goes to a temporary file, not to a pipe. `capture_output=True`
    # deadlocks on the tracked-docx check: that script drives Word through COM,
    # WINWORD.EXE inherits the pipe handles, outlives the direct child, and the
    # parent then waits for an EOF that arrives only when Word exits -- which it
    # does not. The symptom was a run that printed nothing and never finished.
    # A file handle has the same inheritance problem but not the same
    # consequence: subprocess.run waits for the child it started, and the file
    # is read after that, so a lingering Word process costs nothing.
    with tempfile.TemporaryFile(mode="w+", encoding="utf-8",
                                errors="replace") as fh:
        p = subprocess.run([PY, path], cwd=ROOT, env=env,
                           stdout=fh, stderr=subprocess.STDOUT)
        fh.seek(0)
        out = fh.read()
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
    for check in CHECKS:
        label, script, reads = check[:3]
        requires = check[3] if len(check) > 3 else None
        if requires and not os.path.exists(os.path.join(ROOT, requires)):
            print("\n" + "=" * 74)
            print("[%s]" % label)
            print("=" * 74)
            print("  SKIP  %s is not present." % requires)
            print("        It is peer-review correspondence and is never published,")
            print("        so this check runs only in the working tree that holds it.")
            print("        Nothing was verified here; nothing failed.")
            continue
        failures += run(label, script, reads)

    print("\n" + "=" * 74)
    print("[figure size against the PLOS ceiling]")
    print("=" * 74)
    failures += check_figure_specs()

    print("\n" + "=" * 74)
    print("[every Supporting Information item is cited]")
    print("=" * 74)
    failures += check_si_citations()

    print("\n" + "=" * 74)
    print("[cross-document consistency]")
    print("=" * 74)
    failures += check_cross_document()

    print("\n" + "=" * 74)
    print("[scoped cohort, module and Fig 1 counts]")
    print("=" * 74)
    failures += check_scoped_counts()

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
