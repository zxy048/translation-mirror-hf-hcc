#!/usr/bin/env python
"""
R18_verify_manuscript.py -- cross-check every numeric claim in
Manuscript_PLOSONE_v2.md against the tables the pipeline actually produced.

Why this exists: v1 shipped a sentence ("all 33 translation-related pathways
showed positive effect sizes in HCC") that its own supplementary table
contradicted, and a binomial P-value with no implementation anywhere in the
repository. Neither was caught by reading the manuscript.

Design note: an earlier draft of this script only asserted that a *prose*
literal appeared somewhere in the manuscript. That version passed its TF check
while actually reading the wrong row of S30 (+0.583 instead of the quoted
+0.237), because the quoted string happened to exist elsewhere. So each claim
below now asserts two independent things:
  (1) the display string computed FROM THE TABLE appears in the manuscript, and
  (2) the surrounding prose literal appears in the manuscript.
A claim is a hard failure unless both hold.

Run: python R18_verify_manuscript.py
Exit code 0 = all claims verified, 1 = at least one mismatch.
"""
import csv, math, os, re, sys

sys.stdout.reconfigure(encoding="utf-8")

TAB = "D:/R_projects/revision_analysis/v2_output/tables"
MS = "D:/R_projects/revision_analysis/Manuscript_PLOSONE_v2.md"

with open(MS, encoding="utf-8") as fh:
    ms = fh.read()

def rows(name):
    with open(os.path.join(TAB, name), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")

def fnum(v, dec=3):
    """Signed display form with U+2212 minus, matching the manuscript."""
    return ("+" if v >= 0 else "\u2212") + f"{abs(v):.{dec}f}"

PASS, FAIL = [], []
def check(claim, computed, literal=None):
    """computed: display string derived from the table; must appear in the ms.
       literal:  surrounding prose; must also appear. None = table-only claim."""
    ok_v = computed in ms
    ok_p = True if literal is None else literal in ms
    ok = ok_v and ok_p
    (PASS if ok else FAIL).append((claim, computed, literal, ok_v, ok_p))
    tag = "OK  " if ok else "FAIL"
    why = "" if ok else ("  <-- value not in ms" if not ok_v else "  <-- prose not in ms")
    print(f"  {tag} {claim}: computed={computed!r} prose={(literal or '-')!r}{why}")

def forbid(claim, text):
    """A claim the manuscript must NOT make."""
    ok = text not in ms
    (PASS if ok else FAIL).append((claim, f"absent:{text}", None, ok, True))
    print(f"  {'OK  ' if ok else 'FAIL'} {claim}: must not contain {text!r}")

def spearman(xs, ys):
    def rank(z):
        order = sorted(range(len(z)), key=lambda i: z[i])
        rk = [0.0] * len(z)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and z[order[j + 1]] == z[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                rk[order[k]] = avg
            i = j + 1
        return rk
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    nu = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    de = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return nu / de

print("\n--- universe and network construction ---")
s1 = rows("Table_S1_gene_universe_loss_chain.csv")
u = [r for r in s1 if "universe" in r["stage"].lower()]
check("common universe U", f"{int(num(u[0]['n'])):,}", "10,213 genes")

s4 = rows("Table_S4_module_inventory.csv")
nHF = sum(1 for r in s4 if r["disease"] == "HF" and r["module"].lower() != "grey")
nHCC = sum(1 for r in s4 if r["disease"] == "HCC" and r["module"].lower() != "grey")
check("HF named modules", f"{nHF} named co-expression modules in HF")
check("HCC named modules", f"{nHCC} in HCC")

s5 = rows("Table_S5_soft_threshold_fit.csv")
comm = [r for r in s5 if r["chosen_common"] in ("TRUE", "True", "1")][0]
check("common beta", "β = 14")
for d, lit in (("HF", "0.861"), ("HCC", "0.778")):
    row = [r for r in s5 if r["disease"] == d and r["power"] == comm["power"]][0]
    check(f"{d} fit at common beta", f"{num(row['fit']):.3f}", lit)

print("\n--- translation module designation ---")
s9 = rows("Table_S9_translation_overrepresentation.csv")
for d, m in (("HF", "purple"), ("HCC", "magenta")):
    r = [x for x in s9 if x["disease"] == d and x["module"].lower() == m][0]
    check(f"{d} {m} module size", f"{d} {m} ({int(num(r['n_module']))} genes)")
    check(f"{d} {m} Tier A overlap", str(int(num(r["n_overlap"]))))
    check(f"{d} {m} fold", f"{num(r['fold']):.1f}-fold")

print("\n--- module preservation ---")
s19 = rows("Table_S19_module_preservation.csv")
for ref, mod, prose in (("HF", "purple", "rank 2 of 23"), ("HCC", "magenta", "rank 6 of 20")):
    sub = sorted([r for r in s19 if r["reference"] == ref and r["module"].lower() != "grey"],
                 key=lambda r: -num(r["Zsummary.pres"]))
    i = [k for k, r in enumerate(sub) if r["module"].lower() == mod][0] + 1
    check(f"{ref} {mod} Zsummary", f"{num(sub[i-1]['Zsummary.pres']):.2f}")
    check(f"{ref} {mod} rank", f"rank {i} of {len(sub)}", prose)
for ref in ("HF", "HCC"):
    sub = [r for r in s19 if r["reference"] == ref and r["module"].lower() != "grey"]
    strong = sum(1 for r in sub if num(r["Zsummary.pres"]) > 10)
    none = sum(1 for r in sub if num(r["Zsummary.pres"]) < 2)
    n = len(sub)
    # the HF sentence is spelled out; the HCC one is elliptical, so build each
    # in the form the manuscript actually uses.
    display = (f"{strong} of {n} are strongly preserved and {none} of {n} show no evidence"
               if ref == "HF" else f"{strong} of {n} and {none} of {n} respectively")
    check(f"{ref} preservation classes", display)

print("\n--- mirror perturbation ---")
s18 = rows("Table_S18_mirror_summary.csv")
prim = max([r for r in s18 if r["tier"] == "A"], key=lambda r: int(r["n_mirror_sets"]))
check("primary mirror (Tier A)", f"{prim['n_mirror_sets']}/{prim['n_sets']}",
      f"**{prim['n_mirror_sets']} of {prim['n_sets']} Tier A sets**")
check("primary mirror P", "permutation *P* = 1 × 10⁻⁴")
check("cluster mirror", f"{prim['n_mirror_clusters']} of {prim['n_clusters']} clusters")
check("cluster mirror P", f"*P* = {num(prim['p_perm_clusters']):.4f}")
t = max([r for r in s18 if r["tier"] != "A"], key=lambda r: int(r["n_mirror_sets"]))
check("Tier C mirror", f"{t['n_mirror_sets']} of {t['n_sets']} sets")

print("\n--- effect sizes and cross-cohort correlation ---")
s14 = rows("Table_S14_effect_sizes_within_disease.csv")
g = {}
for r in s14:
    g.setdefault((r["set"], r["tier"]), {})[r["cohort"]] = num(r["hedges_g"])

def rho(a, b, tier=None):
    xs, ys = [], []
    for (s, t), v in g.items():
        if tier and t != tier:
            continue
        if a in v and b in v and not (math.isnan(v[a]) or math.isnan(v[b])):
            xs.append(v[a]); ys.append(v[b])
    return spearman(xs, ys)

check("rho Tier A (TCGA vs GSE57338)", f"ρ = {fnum(rho('TCGA_LIHC','GSE57338','A'))}")
check("rho all sets", f"ρ = {fnum(rho('TCGA_LIHC','GSE57338'))}")
for hcc, prose in (("TCGA_LIHC", "−0.550 (TCGA-LIHC, 10/15"),
                   ("GSE14520_GPL3921", "−0.702 (GSE14520 GPL3921, 9/14)"),
                   ("GSE14520_GPL571", "−0.764 (GSE14520 GPL571, 9/15)"),
                   ("GSE76427", "−0.587 (GSE76427, 9/14)")):
    check(f"matrix rho {hcc} vs GSE57338", fnum(rho(hcc, "GSE57338", "A")), prose)
for hf, prose in (("GSE116250", "−0.393 to +0.211"), ("GSE141910", "−0.457 to −0.059")):
    v = [rho(h, hf, "A") for h in ("TCGA_LIHC", "GSE14520_GPL3921",
                                   "GSE14520_GPL571", "GSE76427")]
    check(f"{hf} rho range", f"{fnum(min(v))} to {fnum(max(v))}", prose)

print("\n--- cohort heterogeneity and aetiology ---")
s21 = rows("Table_S21_cohort_heterogeneity.csv")
for d in ("HF", "HCC"):
    # Tier A only, which is the family both the sentence and S5 Fig are about.
    # This used to average all 41 scored sets while quoting a figure that plots
    # Tier A alone; the two medians differed by about 1 percentage point, and
    # the fix was to scope the text to Tier A rather than to change the figure.
    # I2 is undefined for sets present in fewer than two cohorts, so those rows
    # are dropped.
    v = sorted(num(r["I2_pct"]) for r in s21
               if r["disease"] == d and r.get("tier") == "A"
               and not math.isnan(num(r["I2_pct"])))
    med = v[len(v)//2] if len(v) % 2 else (v[len(v)//2-1] + v[len(v)//2]) / 2
    # the manuscript gives the HF value in full and the HCC one as "... and X% for HCC"
    display = (f"median *I*² = {med:.1f}% for {d}" if d == "HF" else f"{med:.1f}% for {d}")
    check(f"median I2 {d}", display)
    pcol = [num(r["p_Q"]) for r in s21
            if r["disease"] == d and r.get("tier") == "A" and not math.isnan(num(r["p_Q"]))]
    pcol.sort()
    pmed = pcol[len(pcol)//2] if len(pcol) % 2 else (pcol[len(pcol)//2-1] + pcol[len(pcol)//2]) / 2
    # the manuscript writes these as "2.1 × 10⁻⁴", not as "2.1e-04"
    mant, expo = f"{pmed:.1e}".split("e")
    sup = "".join(chr(0x2070 + int(c)) if c.isdigit() else chr(0x207B) for c in f"{int(expo):d}")
    check(f"median p_Q {d}", f"{mant} × 10{sup}")

s22d = rows("Table_S22d_DCM_across_cohorts.csv")
xs = [num(r["g_GSE57338_DCM"]) for r in s22d]
ys = [num(r["g_GSE141910_DCM"]) for r in s22d]
same = sum(1 for a, b in zip(xs, ys) if (a > 0) == (b > 0))
check("DCM cross-cohort rho", f"ρ = {fnum(spearman(xs, ys))}")
check("DCM cross-cohort concordance", f"{same}/{len(xs)} sets concordant")

print("\n--- R15 GSE89377 cirrhosis control ---")
s35 = rows("Table_S35_GSE89377_HCC_vs_cirrhosis_specificity.csv")
ta = [r for r in s35 if r["tier"] == "A"]
same = sum(1 for r in ta if r["cirrhosis_same_direction_as_HCC"] in ("TRUE", "True", "1"))
check("R15 Tier A concordant", f"{same} of {len(ta)} Tier A sets")
check("R15 Tier A rho", f"ρ = {fnum(spearman([num(r['Cirrhosis']) for r in ta], [num(r['HCC']) for r in ta]))}")
up = sum(1 for r in s35 if r["cirrhosis_same_direction_as_HCC"] in ("TRUE", "True", "1")
         and num(r["HCC"]) > 0)
dn = sum(1 for r in s35 if r["cirrhosis_same_direction_as_HCC"] in ("TRUE", "True", "1")
         and num(r["HCC"]) < 0)
pct = round(100 * (up + dn) / len(s35))
rho35 = fnum(spearman([num(r["Cirrhosis"]) for r in s35], [num(r["HCC"]) for r in s35]))
check("R15 all-set concordance", f"{pct}% ({up} up, {dn} down; ρ = {rho35})")
rib = [r for r in s35 if r["set"] == "KEGG_RIBOSOME"][0]
stages = (("CH_low", "chronic hepatitis, low grade"), ("CH_high", "chronic hepatitis, high grade"),
          ("Cirrhosis", "cirrhosis"), ("DN_low", "dysplastic nodule, low grade"),
          ("DN_high", "dysplastic nodule, high grade"), ("eHCC", "early HCC"), ("HCC", "HCC"))
check("KEGG_RIBOSOME trajectory",
      " → ".join(f"{num(rib[c]):.2f} ({lab})" for c, lab in stages))

print("\n--- TF association (TCGA-LIHC, all samples, canonical TCS) ---")
s30 = rows("Table_S30_TF_translation_association.csv")
for tf, prose in (("ATF4", "1.6 × 10⁻⁶"), ("MYC", "2.2 × 10⁻⁶"), ("DDIT3", "8.2 × 10⁻⁵"),
                  ("EIF4EBP1", "9.3 × 10⁻⁴"), ("EIF2AK3", None), ("ERN1", None), ("MTOR", None)):
    r = [x for x in s30 if x["tf"] == tf and x["cohort"] == "TCGA_LIHC"
         and x["arm"] == "all" and x["score"] == "TCS"][0]
    v = fnum(num(r["r_adjusted"]))
    check(f"TF {tf} adjusted r (TCS, all)", v,
          f"{v}, FDR = {prose}" if prose else None)

print("\n--- composition adjustment ---")
s25 = rows("Table_S25_composition_adjustment.csv")
ta25 = [r for r in s25 if r["tier"] == "A"]
surv = sum(1 for r in ta25 if r["survives"] in ("TRUE", "True", "1"))
att = sorted(num(r["attenuation"]) for r in ta25)
med = att[len(att)//2] if len(att) % 2 else (att[len(att)//2-1] + att[len(att)//2]) / 2
vifs = [num(r["max_vif"]) for r in ta25]
check("composition survivorship", f"{surv} of {len(ta25)} Tier A set–cohort pairs")
check("median attenuation", f"median attenuation of {100*med:.1f}%")
check("max VIF", f"maximum VIF {max(vifs):.2f}")

print("\n--- age and sex adjustment (R2 Major 3) ---")
s37 = rows("Table_S37_clinical_covariate_adjustment.csv")
s38 = rows("Table_S38_age_sex_mirror_recheck.csv")
ta37 = [r for r in s37 if r["tier"] == "A"]
att37 = sorted(num(r["attenuation"]) for r in ta37)
med37 = att37[len(att37)//2] if len(att37) % 2 else \
        (att37[len(att37)//2-1] + att37[len(att37)//2]) / 2
# NB: `check` matches the manuscript as raw markdown, so these strings carry the
# emphasis markers the manuscript actually contains. Written without them the
# check fails on text that is present, which is a false alarm rather than a
# finding -- the opposite of the failure mode this file exists to catch.
check("age/sex median attenuation",
      f"Median attenuation of |*g*| is {100*med37:.1f}%")
check("age/sex max VIF", f"maximum VIF {max(num(r['max_vif']) for r in s37):.2f} "
                         f"across {len(s37)} set–cohort fits")
check("age/sex sex confound count",
      f"{sum(1 for r in ta37 if num(r['p_sex']) < 0.05)} of {len(ta37)} "
      f"Tier A set–cohort fits at *P* < 0.05")
check("age/sex age confound count",
      f"{sum(1 for r in ta37 if num(r['p_age']) < 0.05)} of {len(ta37)} for age")
# the claim the paragraph actually rests on: the mirror count does not move
for r in s38:
    if r["tier"] != "A":
        continue
    u, a = int(num(r["mirror_unadjusted"])), int(num(r["mirror_adjusted"]))
    if r["hf_cohort"] == "GSE57338":
        check("mirror count, discovery HF, before/after",
              f"GSE57338 gives {u} of {int(num(r['n_sets']))} Tier A sets "
              f"before and after adjustment")
        check("discovery HF rho, adjusted",
              f"from ρ = {num(r['rho_unadjusted']):.3f} to "
              f"{num(r['rho_adjusted']):.3f}".replace("-", "−"))
    else:
        check(f"mirror count, {r['hf_cohort']}, before/after",
              f"{r['hf_cohort']} {u} of {int(num(r['n_sets']))} both times")

print("\n--- hub-gene threshold (S31) ---")
s31 = rows("Table_S31_gene_significance_threshold.csv")
for r in s31:
    lab = f"{r['cohort']} {r['module']}"
    a, b, c = (int(num(r["n_absGS_gt_020"])), int(num(r["n_absGS_gt_050"])),
               int(num(r["n_FDR05_alone"])))
    nested = int(num(r["n_absGS_gt_050_and_FDR05"])) == b
    # the criteria are nested and |GS| > 0.50 is the strictest: it must retain
    # fewer genes than v1's cut and fewer than the FDR criterion.
    ok = nested and b < a and b < c
    (PASS if ok else FAIL).append((f"{lab} nesting", f"{a}/{b}/{c}", None, ok, True))
    print(f"  {'OK  ' if ok else 'FAIL'} {lab} nesting: |GS|>0.20={a} |GS|>0.50={b} FDR<0.05={c} "
          f"(all {b} top genes pass FDR: {nested})")
def trio(r):
    return (int(num(r["n_absGS_gt_020"])), int(num(r["n_absGS_gt_050"])),
            int(num(r["n_FDR05_alone"])))

hf31 = [r for r in s31 if r["cohort"] == "GSE57338"][0]
hm31 = [r for r in s31 if r["cohort"] == "TCGA_LIHC"][0]
check("S31 counts quoted for HF purple",
      "%d, %d and %d for HF purple" % trio(hf31))
check("S31 counts quoted for HCC magenta",
      "%d, %d and %d for HCC magenta" % trio(hm31))
forbid("stringency claim (false: FDR admits the most genes)",
       "the more stringent of the three on this data")
forbid("stringency claim, table caption form",
       "showing that the FDR criterion is the more stringent")

# ---------------------------------------------------------------------------
# orphan sweep
# ---------------------------------------------------------------------------
# The checks above are a hand-written list, so they only cover numbers someone
# thought to enumerate. This sweep is the complement: every numeric literal in
# the body must match at least one cell in at least one generated table, modulo
# a short, explicit list of values that are legitimately *derived* rather than
# stored (medians and percentages computed from a table's rows).
#
# It exists because a hand-written list is exactly how the cohesion figure
# "mean pairwise r = 0.3916 against a null of 0.0239 (Z = 26.0)" survived
# unchallenged: it was correct, but it came from an ad-hoc diagnostic log and
# no table carried it. The sweep asks the general question that spot-checking
# cannot -- is there any *other* number in this manuscript with no artefact
# behind it?
#
# Matching is numeric, not textual: 11.03 in the prose is the correct rounding
# of the stored 11.0259251028192, and a substring test calls that a failure.
print("\n--- orphan sweep (every body number must match a table cell) ---")

# Written as escapes rather than literals: the superscript block is easy to
# mistype, and a wrong codepoint here fails silently (the sweep simply matches
# nothing) rather than raising.
#
# BOTH blocks are needed. The manuscript is internally inconsistent: "I^2" uses
# the Latin-1 U+00B2 (16 times) while the exponents in "10^-57" use the
# superscript block U+2075/U+2077, and one exponent -- "10^-15" -- is built from
# U+207B followed by the Latin-1 U+00B9. That renders identically and is not
# worth churning the text over, but tooling that knows only one block breaks on
# it: a class covering U+2070-U+2079 stops at the U+00B9 and yields an exponent
# of "-", which is a ValueError rather than a wrong answer.
SUP = {chr(0x2070 + i): str(i) for i in range(10)}
SUP[chr(0x207B)] = "-"
SUP[chr(0x207A)] = "+"
SUP[chr(0x00B9)] = "1"
SUP[chr(0x00B2)] = "2"
SUP[chr(0x00B3)] = "3"
SUP_EXP = "[" + "".join(SUP) + "]"

cells = set()
for fn in os.listdir(TAB):
    if not fn.startswith("Table_") or not fn.endswith(".csv"):
        continue            # SMOKE_ tables are not artefacts of record
    with open(os.path.join(TAB, fn), encoding="utf-8-sig") as fh:
        for row in csv.reader(fh):
            for c in row:
                try:
                    cells.add(float(c))
                except (TypeError, ValueError):
                    pass

body = ms.split("## References")[0]          # DOIs are not results
body = body[:body.find("**S1 Fig.**")]      # captions describe, they do not assert

def matches(v, tol):
    return any(abs(c - v) <= tol for c in cells)

orphans = set()
# Spans already consumed by a scientific-notation expression. The plain-decimal
# pass must skip them, or the mantissa of "5.7 x 10^-57" is judged a second time
# on its own and fails against a table that stores 5.68e-57 rather than 5.7.
spent = []
# The exponent class must include U+207B (superscript minus), not stop at U+207A
# (superscript plus): every negative exponent in the manuscript uses ⁻, so a
# class ending at ⁺ silently matches no exponent and this whole pass becomes a
# no-op -- which is how "5.7 x 10^-57" slipped through as a bare "5.7".
for m in re.finditer(r"(\d+(?:\.\d+)?)\s*[×x]\s*10(" + SUP_EXP + r"+)", body):
    mant, sup = m.group(1), m.group(2)
    spent.append(m.span())
    digits = "".join(SUP.get(ch, "") for ch in sup)
    if not any(ch.isdigit() for ch in digits):
        # An exponent with no digits is a source defect, not a number to check.
        # Report it rather than crashing on int("").
        orphans.add(m.group(0))
        continue
    exp = int(digits)
    dp = len(mant.split(".")[1]) if "." in mant else 0
    if not matches(float(mant) * 10 ** exp, 0.5 * 10 ** (exp - dp)):
        orphans.add(m.group(0))
for m in re.finditer(r"(?<![\w.])(\d+\.\d+)(?![\w])", body):
    if any(a <= m.start() < b for a, b in spent):
        continue
    s = m.group(1)
    if not matches(float(s), 0.5 * 10 ** -len(s.split(".")[1])):
        orphans.add(s)

# Values the manuscript states as an aggregate of a table's rows rather than as
# a stored cell. Both are recomputed and asserted earlier in this script, so
# they are checked -- just not by this sweep.
ALLOWED = {
    "34.1": "median attenuation over Tier A pairs (recomputed above from S25)",
}

unexplained = sorted(o for o in orphans if o not in ALLOWED)
for o in unexplained:
    FAIL.append((f"orphan number {o}", o, None, False, True))
    print(f"  FAIL {o}: no table cell rounds to this value")
for o, why in sorted(ALLOWED.items()):
    print(f"  ok   {o}: {why}")

# An allow-list rots in the silent direction: an entry added while its number was
# unsourced stays behind after that number acquires a table, and from then on
# quietly whitelists it. So each entry is re-checked in BOTH directions -- the
# number must still be in the manuscript, and must still be unsourced. An earlier
# version of this check tested `o in orphans and matches(o)`, which is a
# contradiction and could never fire; the allow-list it guarded had already gone
# stale, and the sweep reported a clean pass anyway.
for o in sorted(ALLOWED):
    if o in orphans:
        continue
    present = re.search(r"(?<![\w.])" + re.escape(o) + r"(?![\w])", body)
    why = ("now matches a table cell" if present
           else "no longer appears in the manuscript")
    FAIL.append((f"stale allow-list entry {o}", o, None, False, True))
    print(f"  FAIL {o}: {why} -- delete it from ALLOWED")

if not unexplained:
    print(f"  no unexplained numbers ({len(orphans)} scanned, "
          f"{len(ALLOWED)} allow-listed aggregates)")

print(f"\n{'='*72}\nVERIFIED {len(PASS)} claims | FAILED {len(FAIL)}")
if FAIL:
    print("\nFAILURES:")
    for c, v, l, okv, okp in FAIL:
        why = "value missing" if not okv else "prose missing"
        print(f"  - {c}: {v!r} ({why})")
    sys.exit(1)
print("All manuscript numbers reconcile with the generated tables.")
