#!/usr/bin/env python
"""
R24_audit_response_letter.py -- check the numbers in the response letter.

Why this exists
---------------
R18 verifies the manuscript against the generated tables, and R22 verifies the
supplement. Nothing verified the response letter, which quotes more than two
hundred distinct numbers back at four reviewers. A letter is the one document a
reviewer reads closely line by line, so a number in it that appears nowhere in
the manuscript or the tables is both a claim of ours that nothing supports and
the kind of thing an editor checks first.

The pool a letter number may come from is the manuscript plus every cell of
every generated table, in the rounded forms prose actually uses. The manuscript
is in the pool rather than the tables alone because R18 has already established
that every manuscript number traces to a table; going through the manuscript
lets this file stay a check on *agreement between the two documents* instead of
re-deriving R18's work.

Three parts:

  A. Sweep. Every numeric token in the letter must be in the pool. Anything
     absent is either a defect or a value the letter derives itself, and the
     second kind has to be declared in DERIVED -- where it is re-*computed*
     rather than exempted, so the declaration cannot rot into a silent
     allow-list. (An allow-list rots in the silent direction: an entry added
     while its number was wrong stays behind after the number is fixed, and
     from then on vouches for it.)

     The sweep's power is bounded, and it is worth being plain about the bound.
     The pool holds every table value in every rounding prose might use, some
     65,000 of them, so a wrong number that happens to collide with an unrelated
     cell passes. Mutation testing showed this directly: changing 2.9% to 2.94%
     in a prose sentence is not caught, because something in some table rounds
     to 2.94. Part A therefore catches numbers that fit the pipeline nowhere --
     which is what it found on first run, v1's 1,315 and 5,003 and the
     undeclared 53.4 -- and part C is where the load-bearing figures are pinned
     to the rows they come from.

  B. Headline claims. The mirror counts, the preservation Zsummary values, the
     gene-significance thresholds and the covariate correlations are recomputed
     from the tables they should come from and required to appear in the letter
     in the form the letter uses. Part A only proves a number exists somewhere;
     this proves it is attached to the right claim.

  C. No superseded values. Numbers that were corrected during the revision must
     not survive anywhere in the letter.

Reading the letter is not publishing it: this script opens the file and writes
nothing. The letter is gitignored and stays out of the repository.

Run:  python Code/Revision_v2/R24_audit_response_letter.py
Exit code 0 = every letter number is corroborated, 1 = at least one is not.
"""
import csv
import glob
import os
import re
import statistics
import sys

sys.stdout.reconfigure(encoding="utf-8")

ROOT = "D:/R_projects/revision_analysis"
TAB = os.path.join(ROOT, "v2_output", "tables")
MS = os.path.join(ROOT, "Manuscript_PLOSONE_v2.md")
# Overridable so the checks can be pointed at a deliberately broken copy. Every
# part of this file has been mutation-tested this way: a check that has never
# been seen to fail is not known to be a check, and R22 in this repository spent
# a full build reporting a false pass because exactly that was never done.
LETTER = os.environ.get(
    "R24_LETTER", os.path.join(ROOT, "v2_output", "Response_to_Reviewers.md"))


def norm(t):
    """Fold the unicode dashes to ASCII so -0.550 and U+2212 0.550 compare equal."""
    return t.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")


def tokens(t):
    """Numeric tokens, with thousands separators removed first. Without that
    step "5,003" splits into "5" and "003" and the fragment "003" is reported as
    an uncorroborated number, which is a false alarm about our own tokenizer
    rather than a finding about the letter.

    The sign is only read as a sign when nothing word-like precedes it. "top-
    3,000" is one number with a hyphenated prefix, not negative 3,000; a range
    written "0.90-1.42" is two positive bounds; and the "-9" in the pooled-SD
    formula "4(n1+n2)-9" is subtraction, not a negative constant. Each would
    otherwise enter as a negative value matching nothing and be reported as a
    defect rather than as the tokenizer's own reading error.
    """
    t = re.sub(r"(?<=\d),(?=\d{3}(?!\d))", "", norm(t))
    # A trailing sentence period is not part of the number: "set to 200." would
    # otherwise yield the token "200." and miss the table cell "200".
    return re.findall(r"(?<![\w.)])-?\d+(?:\.\d+)?", t)


def load_rows(name):
    with open(os.path.join(TAB, name), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def fmt(v, dec=4):
    """Signed display form with the unicode minus the documents use."""
    return ("\u2212" if v < 0 else "") + f"{abs(v):.{dec}f}"


def spearman(x, y):
    """Rank correlation, computed here rather than imported: scipy is not a
    dependency of this repository's verification chain, and the only property
    needed is that ties break consistently (average ranks are not required for
    the tables checked here, none of which contain tied effect sizes)."""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return cov / den


# --------------------------------------------------------------------- pool ---
pool = set()
with open(MS, encoding="utf-8") as fh:
    pool |= set(tokens(fh.read()))

TABLE_FILES = sorted(glob.glob(os.path.join(TAB, "*.csv")))
for path in TABLE_FILES:
    with open(path, encoding="utf-8-sig") as fh:
        for row in csv.reader(fh):
            for cell in row:
                for tok in tokens(cell):
                    pool.add(tok)
                    v = num(tok)
                    if v != v:
                        continue
                    # Prose rounds table cells, so the rounding is part of the
                    # pool. Without it every rounded value reads as uncorroborated.
                    for d in (1, 2, 3, 4):
                        pool.add(f"{round(v, d):.{d}f}")

with open(LETTER, encoding="utf-8") as fh:
    letter = fh.read()
letter_n = norm(letter)

# ------------------------------------------------------------------ part A ---
print("--- A. every number in the letter is corroborated ---")

# Values the letter computes for itself. Re-verified below in part B rather than
# exempted here, so this dict cannot silently vouch for a wrong number.
DERIVED = {"53.4": "55 of 103 Tier A pairs survive composition adjustment"}

# Values the letter quotes *because* they are superseded: v1's module sizes and
# gene filters, reproduced from the reviewer's own comment or from the withdrawn
# analysis in order to record that they no longer hold. These are legitimately
# absent from the manuscript, and that absence is the check -- the letter's
# framing of each as superseded is true exactly as long as the manuscript does
# not report it. Should one ever reappear upstream, this fires.
QUOTED = {
    "1315": "v1 HCC module size, quoted in the corrections table as not reproducible",
    "5003": "v1 HCC pre-filter gene count, quoted from the reviewer's comment",
    "6000": "v1 HF variance filter, quoted from the reviewer's comment",
    "3000": "v1 HF variance filter, quoted as the mechanism that was removed",
}

# Facts about the submission package rather than results, so no table holds
# them. Each names the part of this file that verifies it against the artefact
# itself; exempting one without such a part would be an allow-list entry again.
OPERATIONAL = {
    "600": "figure export resolution, checked against the built TIFFs in part E",
}

toks = tokens(letter)
distinct = sorted(set(toks), key=lambda s: (len(s), s))
uncorroborated = [t for t in distinct
                  if t not in pool and t not in DERIVED
                  and t not in QUOTED and t not in OPERATIONAL]

print(f"  {len(toks)} numeric tokens, {len(distinct)} distinct")
print(f"  {len(TABLE_FILES)} table file(s), {len(pool):,} values in the pool")
if uncorroborated:
    print(f"  {len(uncorroborated)} not corroborated:")
    for t in uncorroborated:
        ctx = ""
        for ln in letter_n.split("\n"):
            if re.search(r"(?<![\w.])" + re.escape(t) + r"(?![\w])", ln):
                ctx = ln.strip()[:90]
                break
        print(f"    {t:>10}  {ctx}")
else:
    print("  OK   none")

# ------------------------------------------------------------------ part B ---
print("\n--- B. declared values are re-checked, not exempted ---")
FAILS = []
for shown, why in sorted(DERIVED.items()):
    s25 = [r for r in load_rows("Table_S25_composition_adjustment.csv")
           if r["tier"] == "A"]
    surv = sum(1 for r in s25 if r["survives"].strip().upper() == "TRUE")
    got = f"{100 * surv / len(s25):.1f}"
    # Recomputed *and* required to be the string the letter prints. Checking
    # only the recomputation leaves the declaration free to agree with itself
    # while the letter says something else: mutation testing changed the letter
    # to 53.5% and this part still passed, because 53.4 == 53.4.
    printed = shown in letter_n
    ok = got == shown and printed
    if not ok:
        why_bad = ("recomputed " + got) if got != shown else "not printed"
        FAILS.append(f"DERIVED {shown}: {why_bad}")
    print(f"  {'OK  ' if ok else 'FAIL'} derived {shown} = {got}, printed  ({why})")

for shown, why in sorted(OPERATIONAL.items()):
    print(f"  --   operational {shown}  ({why})")

for shown, why in sorted(QUOTED.items()):
    # Must be genuinely absent from the manuscript and every table, or the
    # letter is calling a live value superseded.
    ok = shown not in pool
    if not ok:
        FAILS.append(f"QUOTED {shown}: now present upstream, so it is not superseded")
    print(f"  {'OK  ' if ok else 'FAIL'} quoted {shown} still absent upstream  ({why})")

# ------------------------------------------------------------------ part C ---
print("\n--- C. headline claims match the tables they come from ---")

s18 = load_rows("Table_S18_mirror_summary.csv")
s20 = load_rows("Table_S20_designated_module_preservation.csv")
s31 = load_rows("Table_S31_gene_significance_threshold.csv")
s38 = load_rows("Table_S38_age_sex_mirror_recheck.csv")


def expect(claim, text):
    """The recomputed string has to be in the letter, as written.

    Both sides are dash-folded. fmt() emits U+2212 because that is what the
    documents use, while letter_n is ASCII-normalised, so comparing them raw
    fails on strings that are in fact present -- a false alarm, which in an
    audit costs more than a miss because it teaches the reader to ignore it.
    """
    ok = norm(text) in letter_n
    if not ok:
        FAILS.append(f"{claim}: {text!r} not in the letter")
    print(f"  {'OK  ' if ok else 'FAIL'} {claim}: {text!r}")


def expect_all(claim, pattern, allowed):
    """Every value the letter gives for a statistic must be one the table holds.

    `expect` only requires the right value to appear *somewhere*, which a
    corrupted duplicate evades: the letter quotes Zsummary = 11.03 in two
    places, so overwriting one of them leaves the other to satisfy the check and
    the letter then contradicts itself in print. Mutation testing found exactly
    that, so any statistic the letter quotes more than once is checked this way
    instead.

    `allowed` is the set of values the source table actually contains, so this
    fails on a value from nowhere as well as on a mistyped one.
    """
    found = re.findall(pattern, letter_n)
    bad = [f for f in found if norm(f) not in {norm(a) for a in allowed}]
    ok = bool(found) and not bad
    if not ok:
        FAILS.append(f"{claim}: letter gives {bad or 'no value at all'}, "
                     f"table allows {sorted(allowed)}")
    detail = (f"{len(found)} occurrence(s) {sorted(set(found))}, all in table"
              if ok else f"{bad} not in table")
    print(f"  {'OK  ' if ok else 'FAIL'} {claim}: {detail}")


def s18row(pair, tier):
    return [r for r in s18 if r["pair"] == pair and r["tier"] == tier][0]


# The letter reports the two primary rows as table cells (Reviewer 2, Major 7)
# and the same-technology rows as "9/15" in prose (Reviewer 3, Major 1). Each
# row is checked in the form the letter actually uses for it: a slash string
# would never appear for the primary rows, and asserting it would report a
# mismatch on a number that is present.
for pair, tier, label, form in [
    ("primary (v1 cohorts)", "A", "primary Tier A", "cells"),
    ("primary (v1 cohorts)", "C", "primary Tier C", "cells"),
    ("HCC array vs HF array-adjacent", "A", "same-technology Tier A", "slash"),
    ("HCC array vs HF RNA-seq 2", "C", "GSE141910 Tier C", "slash"),
]:
    r = s18row(pair, tier)
    text = (f"| {r['n_sets']} | **{r['n_mirror_sets']}** |" if form == "cells"
            else f"{r['n_mirror_sets']}/{r['n_sets']}")
    expect(f"mirror count {label} ({form} form)", text)

r = s18row("primary (v1 cohorts)", "A")
expect("primary Tier A mirror (prose form)",
       f"{r['n_mirror_sets']} of {r['n_sets']}")
expect("primary Tier A cluster count", f"| {r['n_mirror_clusters']} |")

# The two cohort pairings the letter reports as null, which is the result the
# whole response turns on -- these are the ones a reader would check.
for pair, label in [("HCC array vs HF RNA-seq", "GSE116250"),
                    ("HCC array vs HF RNA-seq 2", "GSE141910")]:
    r = s18row(pair, "A")
    expect(f"{label} Tier A mirror count (null result)",
           f"{r['n_mirror_sets']}/{r['n_sets']}")

# Preservation: the values that replaced the Fisher's exact test. The column is
# WGCNA's own dotted name, not the "Zsummary" the letter prints. Both are pinned
# across every occurrence, because the letter quotes them in two places each (a
# response and the corrections table) and a drifted duplicate would otherwise
# leave the correct value standing in the other.
zs = {f"{num(x['Zsummary.pres']):.2f}" for x in s20}
for ref, label in [("HF", "HF purple to HCC"), ("HCC", "HCC magenta to HF")]:
    r = [x for x in s20 if x["reference"] == ref][0]
    expect(f"preservation Zsummary {label}",
           f"Zsummary = {num(r['Zsummary.pres']):.2f}")
expect_all("every Zsummary the letter prints is a designated-module value",
           r"Zsummary = (\d+\.\d+)", zs)

# Gene-significance thresholds, the reviewer's own requested cutoffs.
for cohort, label in [("GSE57338", "HF module"), ("TCGA_LIHC", "HCC module")]:
    r = [x for x in s31 if x["cohort"] == cohort][0]
    expect(f"|GS| > 0.50 count, {label}", f"{r['n_absGS_gt_050']}")
    expect(f"FDR < 0.05 count, {label}", f"{r['n_FDR05_alone']}")

# Age/sex adjustment: the article's claim is that it changes nothing, which is
# only meaningful if the before and after values are the ones in the table.
r = [x for x in s38 if x["hf_cohort"] == "GSE57338" and x["tier"] == "A"][0]
expect("GSE57338 age/sex rho before", f"\u03c1 = {fmt(num(r['rho_unadjusted']), 3)}")
expect("GSE57338 age/sex rho after", f"\u03c1 = {fmt(num(r['rho_adjusted']), 3)}")

# Paired-sample sensitivity (S40 Table), added during the revision.
s40c = load_rows("Table_S40c_paired_effect_sizes.csv")
s40d = load_rows("Table_S40d_paired_summary.csv")
total_pairs = sum(int(num(r["n_pairs_analysed"])) for r in s40d)
n_ci = sum(1 for r in s40c
           if str(r["paired_ci_excludes_0"]).lower() == "true")
expect("paired total", f"{total_pairs} matched pairs")
expect("paired estimate count", f"All {len(s40c)} cohort–set estimates")
expect("paired CI support", f"{n_ci} had confidence intervals excluding zero")
for r in s40d:
    expect(f"paired count {r['cohort']}",
           f"{int(num(r['n_pairs_analysed']))} pairs")
expect("paired TCGA subset", "50 of 371 tumours")

# Heterogeneity, which the letter reports instead of a list of possible causes.
s21 = load_rows("Table_S21_cohort_heterogeneity.csv")
for disease in ("HF", "HCC"):
    v = [num(x["I2_pct"]) for x in s21
         if x["disease"] == disease and x["tier"] == "A"]
    expect(f"median I2 across Tier A sets, {disease}",
           f"{statistics.median(v):.1f}%")

# The aetiology stratification: the letter's claim is that the two arms agree
# with each other, which is a claim about the mean effect in each arm.
s22 = [x for x in load_rows("Table_S22_HF_aetiology_stratified.csv")
       if x["tier"] == "A"]
mean_dcm = statistics.mean(num(x["g_DCM_vs_NF"]) for x in s22)
mean_icm = statistics.mean(num(x["g_ICM_vs_NF"]) for x in s22)
expect("aetiology arm means agree in sign and size",
       f"{fmt(mean_dcm, 3)} and {fmt(mean_icm, 3)}")

# Same aetiology across the two cohorts. The letter reports this as a Spearman
# rho, and it has to be checked as one: the Pearson r on the same 15 cells is
# +0.667, so checking the wrong statistic would report the letter as wrong when
# it is right. What makes the result striking is that only 4 of 15 sets point
# the same way while the ranks still order consistently -- so both are checked.
s22d = load_rows("Table_S22d_DCM_across_cohorts.csv")
xa = [num(x["g_GSE57338_DCM"]) for x in s22d]
xb = [num(x["g_GSE141910_DCM"]) for x in s22d]
expect("same aetiology across cohorts (Spearman)",
       f"{'+' if spearman(xa, xb) >= 0 else ''}{spearman(xa, xb):.3f}")
expect("same aetiology across cohorts (concordant sets)",
       f"{sum(1 for u, v in zip(xa, xb) if (u > 0) == (v > 0))} of {len(xa)}")

# Survival. S39 carries a row per aetiology stratum as well as the whole-cohort
# row, and the whole-cohort rows are the ones marked stratum "NA" -- an earlier
# version of this check took the first row matching the cohort and silently read
# the HBV stratum of GSE141198 (HR 1.24) instead of the whole cohort (1.04). It
# passed, because both are above 1, which is exactly the kind of near-miss a
# check has to be written to exclude. So the stratum is named, and each cohort's
# hazard ratio is matched together with its confidence interval rather than
# merely checked for being on the right side of 1.
s39 = load_rows("Table_S39_translation_score_survival.csv")


def s39row(cohort, cov, stratum="NA"):
    return [x for x in s39 if x["cohort"] == cohort and x["score"] == "TCS"
            and x["model"] == "continuous_per_IQR" and x["stratum"] == stratum
            and x["covariates"] == cov][0]


for cohort, form in [("TCGA_LIHC", "{hr:.2f} per IQR [{lo:.2f}-{hi:.2f}]"),
                     ("GSE141198", "{hr:.2f} [{lo:.2f}-{hi:.2f}]"),
                     ("GSE76427", "{hr:.2f} [{lo:.2f}-{hi:.2f}]"),
                     ("GSE14520_GPL3921", "{hr:.2f} [{lo:.2f}-{hi:.2f}]")]:
    x = s39row(cohort, "none")
    expect(f"survival HR {cohort} (whole cohort)",
           form.format(hr=num(x["hr"]), lo=num(x["ci_lo"]), hi=num(x["ci_hi"])))

x = s39row("GSE14520_GPL3921", "none")
# Emphasis markers included: the letter writes these in italics, and a check
# written against the plain text fails on a string that is present.
expect("survival sample and events", f"*n* = {x['n']}, {x['n_events']} deaths")
x = s39row("GSE14520_GPL3921", "purity+proliferation+immune")
expect("survival HR fully adjusted",
       f"HR {num(x['hr']):.2f} [{num(x['ci_lo']):.2f}-{num(x['ci_hi']):.2f}]")

# ------------------------------------------------------------------ part D ---
print("\n--- D. superseded values appear only where they are being corrected ---")
# A plain "must not appear" test is wrong here, and testing it is how that was
# found: the letter *must* contain v1's "P < 0.0001" and the equal-weight pooled
# SD, because the corrections table lists them in order to say what replaced
# them. Absence would mean the correction had been dropped, not that the letter
# was clean.
#
# The invariant that does hold is that every occurrence sits in a correction
# context. A superseded value standing on its own is the defect: the letter
# would be reporting it as current. So each occurrence is required to fall
# within a window of a marker saying it is withdrawn.
MARKERS = ("withdrawn", "removed", "no longer", "false", "superseded",
           "not reproducible", "replaced", "original", "v1", "had no",
           "before it was corrected", "invalid")
WINDOW = 320
SUPERSEDED = {
    "v1 binomial P": "P < 0.0001",
    "v1 equal-weight pooled SD": "sqrt((s\u2081\u00b2+s\u2082\u00b2)/2)",
    "v1 pathway-label permutation": "permuted *pathway labels*",
    "v1 top-6,000 variance filter": "top-6,000",
    "v1 top-3,000 variance filter": "top-3,000",
}
for why, needle in SUPERSEDED.items():
    hay = norm(needle)
    starts = [m.start() for m in re.finditer(re.escape(hay), letter_n)]
    if not starts:
        # Also fine: a superseded value absent entirely needs no context.
        print(f"  OK   {why}: absent from the letter")
        continue
    bad = [s for s in starts
           if not any(k in letter_n[max(0, s - WINDOW): s + WINDOW].lower()
                      for k in MARKERS)]
    if bad:
        FAILS.append(f"superseded value stands uncorrected: {why}")
    print(f"  {'FAIL' if bad else 'OK  '} {why}: {len(starts)} occurrence(s), "
          f"{len(bad)} with no correction marker nearby")

# ------------------------------------------------------------------ part E ---
print("\n--- E. the letter's claim about the figure files is true of them ---")
# The letter states the figures are supplied as 600 dpi TIFF. That is a claim
# about the artefacts in the package, so it is checked against the artefacts.
try:
    from PIL import Image
except ImportError:
    Image = None
PKG = os.path.join(ROOT, "PLOS_ONE_Submission_v2")
for name in ("Fig1.tiff", "Fig2.tiff", "Fig3.tiff"):
    path = os.path.join(PKG, name)
    if not os.path.exists(path):
        FAILS.append(f"{name} missing from the package")
        print(f"  FAIL {name}: not in the package")
        continue
    if Image is None:
        print(f"  --   {name}: present (Pillow unavailable, dpi unchecked)")
        continue
    with Image.open(path) as im:
        # Pillow hands back a rational, not a float, so it has to be cast before
        # it will format or subtract.
        dpi = float(im.info.get("dpi", (0, 0))[0])
        fmt = im.format
    ok = abs(dpi - 600) < 1 and fmt == "TIFF"
    if not ok:
        FAILS.append(f"{name}: {fmt} at {dpi:.0f} dpi, letter says 600 dpi TIFF")
    print(f"  {'OK  ' if ok else 'FAIL'} {name}: {fmt}, {dpi:.0f} dpi")

# ------------------------------------------------------------------ report ---
n_bad = len(uncorroborated) + len(FAILS)
print(f"\n{'=' * 72}")
if n_bad:
    print(f"FAILED: {len(uncorroborated)} uncorroborated number(s), "
          f"{len(FAILS)} claim mismatch(es)")
    for f in FAILS:
        print("  -", f)
    sys.exit(1)
print("Every number in the response letter reconciles with the manuscript "
      "and the generated tables.")
