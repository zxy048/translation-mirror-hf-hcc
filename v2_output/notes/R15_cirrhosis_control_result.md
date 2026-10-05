# GSE89377 cirrhosis control: is the hepatic translation signal cancer-specific? (R15)

Date: 2026-10-01. Script: `Code/Revision_v2/R15_cirrhosis_control.R`.
Probes: `R15_probe_magnitude.R`, `R15_probe_compression.R`. Log: `R15_log.txt`.
Tables: `Table_S34_GSE89377_stage_spectrum_effect_sizes.csv`,
`Table_S35_GSE89377_HCC_vs_cirrhosis_specificity.csv`.
Raw: `intermediate/R15_cirrhosis_control.rds`.

## Why this cohort

R2 Major 3 and R4 #1 press the same point: every HCC cohort v1 used contrasts
tumour against **adjacent non-tumour liver from the same organ**, so "translation
programs are up in HCC" is equally consistent with a liver-disease effect
(inflammation, fibrosis, regeneration). v1 could not separate the two, because it
had no non-malignant liver arm.

GSE89377 (GPL6947, Illumina HumanHT-12 V3) has one. 107 samples across the
multistep hepatocarcinogenesis spectrum:

| stage | prefix | n |
|---|---|---|
| Normal | N | 13 |
| Chronic hepatitis, low grade | FL | 8 |
| Chronic hepatitis, high grade | FH | 12 |
| **Cirrhosis** | **CS** | **12** |
| Dysplastic nodule, low grade | DL | 11 |
| Dysplastic nodule, high grade | DH | 11 |
| Early HCC | eHCC | 5 |
| HCC (TG1/TG2/TG3 pooled) | TG | 35 |

Cirrhosis is the decisive arm: diseased but not malignant.

## Pipeline and its guards

Raw intensities with per-probe detection p-values → keep probes detected
(p < 0.05) in ≥ 50% of samples (17,793 / 48,803) → log2 → quantile normalise →
GPL6947 annot → symbol collapse (11,756 genes). ssGSEA with R06's exact
parameter object, so scores are comparable.

**Guard.** The whole script rests on `SAMPLE k` in the non-normalized file
mapping to column *k* of the series matrix. Rather than assume it, the title
prefixes (`N-`, `FL-`, `CS-`, `TG1-` …) are independently decoded to a stage and
required to agree with the phenotype characteristic sample for sample.
**Agreement was 107/107**, and the script `stop()`s otherwise. Getting this wrong
would relabel every sample with someone else's disease stage with no visible
symptom.

## The HCC arm of this cohort is sound

Before interpreting anything, GSE89377's HCC arm was checked against the four HCC
cohorts already in the pipeline (Tier A sets):

| cohort | sign agreement | Spearman ρ |
|---|---|---|
| TCGA_LIHC | 100% | 0.607 |
| GSE14520_GPL3921 | 100% | 0.790 |
| GSE14520_GPL571 | 93% | 0.596 |
| GSE76427 | 100% | 0.833 |

The cohort reproduces the established HCC pattern.

## Finding: the direction is largely shared with non-malignant cirrhosis

Cirrhosis vs HCC, within GSE89377 (Hedges g vs the same Normal arm):

| scope | both up | both down | discordant | concordance | Spearman ρ |
|---|---|---|---|---|---|
| Tier A (n = 15) | 11 | 1 | 3 | **80%** | **+0.668** |
| all tiers (n = 40) | 21 | 3 | 16 | 60% | +0.532 |

80% of Tier A sets are concordant — 12 of 15 (11 rise in both, 1 falls in both).
Per-set values (Tier A, `Table_S35`):

| set | Cirrhosis | HCC | same direction |
|---|---|---|---|
| REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING | 2.14 | 0.81 | yes |
| KEGG_RIBOSOME | 2.03 | 0.47 | yes |
| REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION | 1.98 | 0.63 | yes |
| REACTOME_EUKARYOTIC_TRANSLATION_INITIATION | 1.87 | 0.58 | yes |
| REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL | 1.76 | 0.93 | yes |
| REACTOME_TRNA_AMINOACYLATION | −0.60 | 0.47 | **no** |
| REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION | −0.85 | 0.35 | **no** |
| KEGG_AMINOACYL_TRNA_BIOSYNTHESIS | −0.75 | 0.22 | **no** |

The three discordant sets are the aminoacyl-tRNA and mitochondrial ones — the
same family the manuscript already flags as behaving differently. The bulk
cytosolic ribosome/translation program, by contrast, is **liver-disease-shared**.

**Consequence for the manuscript:** the hepatic arm of the mirror may not be
described as cancer-specific reprogramming. It is a translation-associated
program shared with non-malignant chronic liver disease, with a smaller
cancer-associated component concentrated in the tRNA-aminoacylation and
mitochondrial sets.

## The trajectory is non-monotone

Across the ordered spectrum, the ribosomal sets peak in chronic hepatitis and
cirrhosis and fall at dysplasia and early HCC (Tier A, KEGG_RIBOSOME):
`CH_low 0.61 → CH_high 1.40 → Cirrhosis 2.03 → DN_low 1.07 → DN_high 0.15 →
eHCC 0.20 → HCC 0.47`.

So the hepatic translation program is **not** switched on at the malignant
transition; it is highest in non-malignant chronic liver disease and drops as
malignancy emerges. A dose-response claim in either direction would be wrong.

## The magnitude comparison survives, but only after a scope check

Within GSE89377, |g| is larger for cirrhosis than for HCC in **14/15** Tier A
sets (mean |g| 1.13 vs 0.51), and for KEGG_RIBOSOME cirrhosis (2.10) exceeds
every external HCC cohort (TCGA 0.64, GPL3921 1.30, GPL571 0.63, GSE76427 1.01).

This invited a confound: GSE89377's HCC arm is compressed relative to the other
HCC cohorts —

| HCC cohort | mean \|g\|, all shared sets |
|---|---|
| GSE14520_GPL3921 | 1.44 |
| GSE76427 | 1.10 |
| GSE14520_GPL571 | 1.08 |
| TCGA_LIHC | 0.74 |
| **GSE89377** | **0.54** |

1.4× to 2.7× smaller. `R15_probe_compression.R` tests whether that deficit is
cohort-wide or specific to the translation sets: the GSE89377/TCGA |g| ratio is
**0.68 for Tier A** versus **0.76 for non-Tier-A** sets — essentially uniform
across tiers. It is a platform/pipeline scaling effect, so the *internal*
cirrhosis-vs-HCC comparison is interpretable, because both arms carry the same
compression.

The reverse check is also the reason the **direction** result is the headline and
the **magnitude** result is not: compression preserves signs and attenuates
magnitudes, so the concordance claim is robust to it and the "cirrhosis moves
more" claim is the fragile one. Report the direction; report the magnitude only
with the uniform-scaling check attached.

## Limitations to disclose

1. **No composition correction.** Unlike R11 for the main cohorts, R15 does not
   adjust for tissue composition. Chronically diseased liver and tumour both
   carry shifted immune/fibrotic/stromal content relative to normal liver, so
   part of the *shared* signal could be a shared composition shift rather than a
   shared cell-intrinsic program. This is a real alternative explanation and
   should be stated rather than buried.
2. **eHCC n = 5.** Every eHCC estimate, and the eHCC-vs-HCC agreement
   (13/15 up-in-both, ρ = +0.343), rests on 5 samples.
3. **TG1/TG2/TG3 pooled.** The 35-sample HCC arm mixes three groups; pooling
   heterogeneous grades inflates within-arm variance and attenuates |g|, which
   may explain part of the compression in the row above.
4. GSE89377 is an array cohort with an aggressive detection filter (36% of probes
   retained); the gene universe is 11,756 and not the discovery universe `U`.

## What this changes

- The mirror's **HF arm is untouched** — this cohort says nothing about heart
  failure.
- The mirror's **HCC arm must be reframed** from "cancer-specific" to
  "shared with non-malignant chronic liver disease, with a cancer-associated
  residual in the aminoacyl-tRNA/mitochondrial sets".
- The Title's "disease-context-dependent" framing survives; a claim that the
  hepatic signal is specific to malignancy does not, and v1 never made one
  explicitly — this prevents v2 from accidentally making it in the Discussion.
- This is a **narrowing** of claims, which is what R2/R3/R4 asked for, and it
  gives a concrete, data-backed answer to R2 Major 3 instead of an argument.
