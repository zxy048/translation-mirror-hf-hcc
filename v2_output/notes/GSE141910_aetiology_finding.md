# GSE141910: what the cohort actually is, and what the aetiology stratification shows

Date: 2026-10-01.  Script: `Code/Revision_v2/R10b_hf_aetiology.R`.
Tables: `Table_S22b/c/d`.  Raw: `intermediate/R10b_hf_aetiology.rds`.

## The annotation

`GSE141910_series_matrix.txt.gz`, 366 samples, per-sample CSV expression:

| aetiology | n | group assigned in R05 |
|---|---|---|
| Non-Failing Donor | 166 | control |
| Dilated cardiomyopathy (DCM) | 166 | disease |
| Hypertrophic cardiomyopathy (HCM) | 28 | disease |
| Peripartum cardiomyopathy (PPCM) | 6 | disease |

## v1 vs v2 use of this dataset — two different uses of one series

**What v1 actually wrote** (`Manuscript_PLOSONE_v1.md:79`): *"hypertrophic
cardiomyopathy (HCM) vs. non-failing myocardium (GSE141910; n = 28 HCM,
n = 166 NF)"*, concluding the mirror was absent (rho = -0.036).

So v1 **did** state the sample sizes, and 28 + 166 does disclose that only the
HCM subset was analysed. An earlier reading of mine — that v1 "described the
dataset as though it were an HCM cohort" — overstated it, and is withdrawn.
**The actual, narrower defect is two omissions:**

1. v1 never disclosed that GSE141910 is a **366-sample cardiomyopathy
   collection** whose HCM arm is 28 of 200 diseased samples; the reader is given
   no way to know the series is 83 % DCM.
2. v1 never justified **selecting the HCM arm** for the same-organ control, when
   the same series contains a larger DCM arm (and the HF discovery cohort,
   GSE57338, is itself DCM + ICM). Selecting the one aetiology that is *not* a
   failing phenotype, to serve as the cardiac control for a heart-failure study,
   needs an explicit rationale. R3 minor #8 asks for exactly this.

**v2 (R05)** instead pools **all 200 diseased samples** as one
`disease = "HF"` group. That is a *change* from v1, and it imports its own
imprecision: the pooled group is cardiomyopathy (83 % DCM / 14 % HCM / 3 %
PPCM), and HCM is a hypertrophic rather than failing phenotype. v2 must
therefore either relabel the group as "cardiomyopathy vs non-failing" or report
the DCM arm as primary — see "Decision taken" below.

## Stratified effect sizes (Hedges g on Tier A sets, vs the 166 donors)

| cohort | subgroup | n | Spearman rho vs TCGA-LIHC | mirror sets (gated) | n gated | mean g |
|---|---|---|---|---|---|---|
| GSE57338 | DCM | 82 | **-0.554** | 11 | 13 | -0.291 |
| GSE57338 | ICM | 95 | **-0.596** | 11 | 13 | -0.431 |
| GSE141910 | DCM | 166 | -0.193 | 1 | 13 | +0.132 |
| GSE141910 | HCM | 28 | -0.100 | 2 | 13 | -0.019 |
| GSE141910 | PPCM | 6 | +0.114 | 3 | 13 | -0.126 |

## What follows

1. **v1's HCM control claim is substantively correct.** It really did test the 28
   HCM samples against the donors, and the effect really is ~0 (rho -0.100 here
   vs -0.036 in v1; the difference is the set definition — 15 Tier A sets vs
   v1's 33 keyword pathways — not the samples). The defect is **transparency**:
   v1 described the dataset as though it were an HCM cohort, when it used a
   28-of-366 subset of a DCM-dominated cardiomyopathy collection. That must be
   stated plainly in v2.

2. **The failure to replicate in GSE141910 is not attributable to its aetiology
   mix.** GSE141910's DCM arm alone — 166 vs 166, the aetiology that dominates,
   and a precisely balanced design — still gives rho -0.193 with 1/13 mirror
   sets. Excluding HCM and PPCM does not rescue the replication.

3. **The discovery result is not driven by aetiology composition.** DCM and ICM
   arms of GSE57338 agree with each other (both 11/13 mirror sets, mean g -0.291
   and -0.431), so the GSE57338 mirror does not depend on which HF aetiology is
   sampled. Confirms R10.

4. **Same-aetiology, same-organ cohorts disagree.** DCM in GSE57338 vs DCM in
   GSE141910: Spearman rho = **+0.350**, and only **4 / 15** Tier A sets point
   the same way. This is the strongest available statement of the cohort
   dependence, and it is the honest answer to R2 #3 and R3 #7 — v1 offered
   "severity, etiology, platform, sample size" as untested explanations;
   etiology is now tested and excluded as *the* explanation.

5. **Caveat that must accompany (4).** GSE57338 is a microarray (Affymetrix
   HG-U133 Plus 2.0) and GSE141910 is RNA-seq, so the DCM-vs-DCM discordance
   still conflates cohort with platform. It excludes aetiology, not platform.
   The same-technology pairing required by R3 #1 is a separate analysis.

## Decision taken

Primary grouping in R05/R06–R12 is left as the pooled 200 ("cardiomyopathy vs
non-failing"). The discordance holds in **every** subset, so no conclusion
depends on the choice; re-keying R05 would invalidate R06–R12 (including the
permutation and bootstrap stages) for no change in any result. The
stratification is reported as a sensitivity analysis (Table S22b/c/d), and the
cohort is relabelled honestly in the v2 manuscript and the response letter.

## Marker-coverage caveat carried alongside this

The GSE57338 discovery array retains only 10/24 `immune_pan` markers — losing
the entire lymphoid panel (CD3D/E/G, CD2, CD8A, CD19, MS4A1, GNLY, KLRD1,
FCGR3A, CXCL9, CXCL10, FOXP3, CTLA4) and both the `lymphocyte` and
`proliferation` panels outright. Its "immune-adjusted" model is therefore
myeloid-only (R11, Table_S23). Must be disclosed wherever that adjustment is
cited.
