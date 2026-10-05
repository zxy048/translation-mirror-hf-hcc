# Module-internal robustness on the harmonised modules (R14)

Date: 2026-10-01. Script: `Code/Revision_v2/R14_tom_robustness.R`.
Diagnostic: `Code/Revision_v2/R14_diag.R` → `R14_diag_log.txt`.
Tables: `Table_S32_module_internal_robustness.csv`, `Table_S33_cross_disease_module_coherence.csv`.
Raw: `intermediate/R14_tom_robustness.rds`.

Full run: 10,000 permutations for **both** nulls, both diseases (4 × 10,000).

## What v1 claimed

`Manuscript_PLOSONE_v1.md:51` (Methods) and `:105` (Results):
"mean intra-modular TOM for the green module (227 genes) was compared against a
null distribution from 10,000 random gene sets of equal size" → **Z = 62.8,
permutation P < 0.0001**, jackknife CV = 1.6%; cross-disease coherence Z = 2.9.
Used at `:105` to argue the module is "a genuine co-expression structure rather
than an artifact of the β = 17 soft-threshold parameter".

## The verification chain (two guessed formulas, both caught)

The metric had to be reimplemented because R00–R13 compute nothing like it. Two
independent mistakes were caught by the script's own SELFTEST rather than by
inspection, and **both would have produced plausible-looking numbers**:

1. **Numerator.** Guessed `(A Aᵀ)_ij + a_ij`; the true form is `(A Aᵀ)_ij − a_ij`
   (since `A Aᵀ` sums over all `u`, `Σ_{u≠i,j} a_iu a_uj = (A Aᵀ)_ij − 2 a_ij`).
   Error ≈ 0.83 in absolute TOM.
2. **Denominator.** Guessed `min(k_i,k_j) + 1 − a_ij`; the form that matches is
   `min(k_i,k_j) − a_ij`, i.e. the documented formula with a connectivity that
   **excludes self**. My "simplification" had double-counted the `+1`. Error
   exactly 1.

Verified form, off-diagonal, with `A = (0.5(1+cor))^power`, `a_ii = 1`,
`k = rowSums(A)`:

    TOM_ij = ( (A Aᵀ)_ij − a_ij ) / ( min(k_i, k_j) − a_ij )

SELFTEST now pins **both** the adjacency (vs `WGCNA::adjacency`) and the TOM (vs
`WGCNA::TOMsimilarity`) and `stop()`s on mismatch. Currently 3.5e-16 (HF) and
7.2e-16 (HCC) on TOM, 2e-14/6e-14 on adjacency.

**Why this matters for the letter:** a wrong numerator and a wrong null shift the
observed statistic and the null in the *same* direction, so Z can look entirely
reasonable while being wrong. This is exactly the failure mode that produced the
v1–v2 inconsistencies the reviewers detected. It is worth one sentence in the
response letter as evidence the pipeline is now self-checking.

## Finding: v1's null is not exchangeable, and it changes the answer

`TOM` is connectivity-normalised — its denominator is `min(k_i,k_j) − a_ij`.
A uniformly random gene set is therefore **not** a valid null unless the module's
genes share the universe's `k` distribution. They do not
(`R14_diag_log.txt`, full universe `n = 10,213`):

| disease | module | module mean `k` | universe mean `k` | ratio |
|---|---|---|---|---|
| HF | purple (153) | 6.6 (median 6.0) | 15.2 (median 10.0) | 0.43 |
| HCC | magenta (146) | 21.1 (median 12.9) | 74.1 (median 30.5) | 0.28 |

Both designated modules sit at **1/4 to 1/2 the universe's connectivity**. Low-`k`
genes have a small numerator *and* a small denominator, so their TOM ratio is
dragged back toward the random ratio. The consequence under v1's own null:

| disease | module | obs mean intra-TOM | unmatched null | Z | P |
|---|---|---|---|---|---|
| HF | purple | 0.00787 | 0.00277 (sd 3.5e-4) | 14.70 | < 0.003 |
| HCC | magenta | 0.01798 | 0.01843 (sd 2.6e-3) | **−0.17** | 0.58 |

Taken at face value this says the HCC translation module is **not** internally
robust — which is false. A connectivity-free measure of the same property says
the opposite, and says HCC magenta is the *tighter* of the two:

| disease | module | mean pairwise `r` | random | Z |
|---|---|---|---|---|
| HF | purple | 0.3403 | 0.0355 (sd 0.0074) | 40.96 |
| HCC | magenta | 0.3916 | 0.0239 (sd 0.0141) | 26.00 |

So the module is real; the null was wrong. **Fix: a connectivity-matched null** —
20 quantile bins of `k` over the universe; each module gene is replaced by a
uniformly drawn member of its own bin (bins excluded of module members). This is
a null-specification fix, not a new biological analysis.

Under the matched null: HF Z = 404, HCC Z = 89.6, both P < 0.003.

### Do NOT compare Z across the two nulls

Matching collapses the null's variance (sd 1.8e-5 vs 3.5e-4 for HF), which
inflates Z **mechanically**. Z = 404 does not mean "more significant" than v1's
62.8, and it is not comparable to it. The interpretable effect size is the
observed/null mean ratio, and that is the column to quote:

| disease | module | obs / unmatched null | obs / matched null |
|---|---|---|---|
| HF | purple | 2.84 | 11.03 |
| HCC | magenta | 0.98 | 4.69 |

The 0.98 → 4.69 change for HCC is the whole story: under a valid null the module
shows ~4.7× enrichment of intra-modular TOM over connectivity-matched random
sets, versus no enrichment at all under v1's null.

## Jackknife: report the right CV

Two different quantities, and v1 did not say which:

- **CV of the mean intra-modular connectivity across iterations** — stable,
  averaging over 153 genes cancels per-gene swings: **2.83% (HF), 2.99% (HCC)**.
  This is the quantity comparable to v1's "CV = 1.6%".
- **Median per-gene CV** — ~50% in both. This is *not* instability: at power 14
  the adjacency is dominated by a few strong partners (r = 0.8 → a = 0.229 vs
  r = 0.5 → a = 0.0178), so dropping 10% of genes frequently removes a gene's
  dominant partner. No gene has kIM < 0.05, so the tail is not degenerate.

Report the first; mention the second only if asked. In the letter, state that v1's
1.6% is taken to be the first quantity — the Methods text was ambiguous.

## Cross-disease coherence

Mean pairwise `r` is not connectivity-normalised, so the uniform null **is**
exchangeable here and no matching is needed. Unchanged from the smoke run:

| module | tested in | genes present | obs mean `r` | null | Z | P |
|---|---|---|---|---|---|---|
| HF purple | HCC | 153/153 | 0.3256 | 0.0264 | 21.14 | < 0.01 |
| HCC magenta | HF | 146/146 | 0.1400 | 0.0358 | 12.16 | < 0.01 |

Both are far above v1's coherence Z = 2.9 — v1 measured the green module's genes
in an *external* HCC matrix (GSE141198), a harder test than the harmonised
common-universe version here. Note the asymmetry: HF purple holds up better in
HCC (r = 0.33) than HCC magenta does in HF (r = 0.14), consistent with the R04
preservation asymmetry (11.03 vs 7.09).

## Consequences for the manuscript and the letter

1. v1's S4 Fig claim **survives**, but only under the corrected null, and the
   number changes (Z = 62.8 for HF green → ratio-based reporting for purple and
   magenta). It must be reported for **both** modules, not just the favourable one.
2. The Methods must state the null used. Writing "10,000 random gene sets of
   equal size" without saying they are connectivity-matched would be the same
   omission v1 made, and R2/R3 are already pressing on exactly this.
3. **The connectivity confound is worth disclosing to the reviewers as a caught
   error.** It is the second instance (after the "all 33 positive" claim) of a
   statistic that was wrong in a direction that flattered the conclusion, and
   saying so explicitly is the strongest available answer to R2 Major 4 and the
   "has the statistical analysis been performed rigorously?" question.
4. R03's module designation is unaffected — it used canonical ssGSEA enrichment,
   not TOM, so low connectivity did not bias which module was chosen.

## Housekeeping

`SMOKE_Table_S32/S33` and `SMOKE_R14_tom_robustness.rds` are smoke-run leftovers
(300 permutations) and must be deleted before packaging. `R14_probe.R` and
`R14_diag.R` are diagnostics, not part of the pipeline; keep them (they document
the verification) but they produce no shipping artefact.
