# The mirror test under the corrected null (R09)

Source: `intermediate/R09_mirror.rds` -> `summary` (8 rows).
Config: sample-level permutation, N_PERM_NULL = 10 000, seed 42.
Mirror requires |g| >= 0.20 in **both** diseases **and** opposite signs.

| comparison | tier | sets | clusters | mirror sets | perm P (sets) | mirror clusters | perm P (clusters) | binom (v1 style) | n_eff Nyholt | n_eff Li&Ji |
|---|---|---|---|---|---|---|---|---|---|---|
| primary (v1 cohorts) | A | 15 | 5 | **10** | **1e-4** | 3 | **0.0031** | 1.0000 | 9.44 | 5.0 |
| primary (v1 cohorts) | C | 33 | 19 | **21** | **1e-4** | 14 | **1e-4** | 0.0636 | 32.37 | 15.0 |
| HCC array vs HF RNA-seq | A | 14 | 5 | 2 | 0.174 | 1 | 0.270 | 0.375 | 9.18 | 5.5 |
| HCC array vs HF RNA-seq | C | 28 | 19 | 5 | 0.093 | 4 | 0.087 | 0.0192 | 30.98 | 16.5 |
| HCC array vs HF RNA-seq 2 | A | 14 | 5 | 0 | 1.000 | 0 | 1.000 | 0.0625 | 8.68 | 5.0 |
| HCC array vs HF RNA-seq 2 | C | 27 | 18 | 11 | **2e-4** | 10 | **1e-4** | 0.815 | 27.03 | 14.0 |
| HCC array vs HF array-adjacent | A | 15 | 5 | **9** | **0.0014** | 3 | **0.0091** | 1.0000 | 9.39 | 5.5 |
| HCC array vs HF array-adjacent | C | 29 | 19 | **19** | **2e-4** | 14 | **2e-4** | 0.0636 | 28.80 | 13.0 |

## 1. The mirror survives the corrected null

This was the plan's flagged risk (mirror may shrink under honest thresholds). It
did not collapse:

- Primary, Tier A: **10/15** mirror sets against a permutation null mean of
  **0.11** sets, perm **P = 1e-4**; at cluster level 3 mirror clusters vs null
  0.06, **P = 0.0031**.
- Primary, Tier C (v1's 33 pathways): 21/33, perm P = 1e-4 sets and clusters.

So the primary finding holds under sample-level permutation, under Jaccard
de-redundancy, and with the effect-size gate applied in both diseases.

## 2. v1's headline statistic was wrong, and differently wrong than assumed

v1 asserted a "binomial test P < 0.0001" for 24/33. The exact two-sided binomial
for 24/33 at p = 0.5 is **P = 0.0135**, and for the v2 counts it is **P = 1.0000**
(Tier A, 10/15) and **P = 0.0636** (Tier C, 21/33). The v1 number is not
reproducible under any reading.

Why the binomial and the permutation disagree so sharply, and why the
permutation is the right test: the mirror event requires two independent
conditions to coincide — a large effect in **both** diseases **and** opposite
signs. Under sample permutation that coincidence occurs in only ~0.7 % of
set-slots (null mean 0.11 of 15), not 50 %. The binomial's p = 0.5 null assumes
each set independently has a 50 % chance of mirroring, which is wrong on the
data. The permutation null is the honest one, and it is **more** stringent, not
less — and the observed count still clears it by a wide margin.

## 3. R3 #1 (same-technology validation) — answered, and the answer is not the hoped one

Within matched technology the mirror is present when the HF arm is the
**array** cohort and absent when the HF arm is either **RNA-seq** cohort:

- HCC array vs HF **array**-adjacent: Tier A 9/15, P = 0.0014 (present)
- HCC array vs HF RNA-seq (GSE116250): Tier A 2/14, P = 0.174 (absent)
- HCC array vs HF RNA-seq 2 (GSE141910): Tier A 0/14, P = 1.000 (absent)

But this is **not** a clean platform effect, because the **primary** comparison
is itself cross-technology (GSE57338 array vs TCGA-LIHC RNA-seq) and it is the
strongest result in the table. The discordance therefore attaches to the two HF
RNA-seq **cohorts**, not to the technology as such — consistent with R10/R10b,
where every aetiology subset of GSE141910 fails equally.

One Tier C row is inconsistent with this reading (HCC array vs HF RNA-seq 2:
11/27, P = 2e-4, significant) while its Tier A row is 0/14. That has to be
reported rather than smoothed over: at Tier C the HF RNA-seq 2 comparison is
significant, at Tier A it is null. The discrepancy is likely driven by the
low-coverage cohorts in Tier C; it must be stated explicitly in the manuscript
and the response letter.

## 4. Consequence for the manuscript

The plan's contingency ("mirror 降为探索性, module architecture 升为主发现") is
**not** triggered — the mirror is significant under the corrected null. But its
scope must be stated precisely: the mirror is a property of the GSE57338-anchored
comparisons and is not recovered in either independent HF RNA-seq cohort. The
abstract, title and conclusions must not claim cross-disease generality beyond
that, which is what R2 #3 and R4 #1 demanded.

## 5. What R04 must still settle

Whether the translation modules are preserved across diseases decides whether
"disease-context-dependent organisation" can be stated as the primary claim.
R04 is running; if the modules are **not** preserved, the module-architecture
finding becomes the headline and the mirror becomes secondary supporting
evidence rather than the title claim.
