# Translation-related transcriptional programs in heart failure and hepatocellular carcinoma

Analysis code for a cross-disease comparison of translation-related transcriptional
programs between failing myocardium and hepatocellular carcinoma.

---

## Which pipeline is authoritative

**`Code/Revision_v2/` is the authoritative pipeline.** Everything else in this
repository is retained for provenance only.

The scripts under `Code/WGCNA/`, `Code/ssGSEA/`, `Code/TF_analysis/` and
`Code/Figure_generation/` produced the first submission. They are superseded, and
the results they produce should not be cited or reused. Read
[What the first analysis got wrong](#what-the-first-analysis-got-wrong) before
using anything in this repository.

> **If you are reading this because you found a number in one of these scripts
> that does not appear in the paper: that is expected.** The v1 numbers are not
> the paper's numbers. The paper's numbers come from `Code/Revision_v2/`.

---

## What the first analysis got wrong

These are disclosed rather than quietly corrected, because the errors are
instructive and because some of them were in the public record.

**The HF soft threshold was hard-coded.** `Code/WGCNA/02_HF_WGCNA.R` sets
β = 12 by hand. Its own comment says the value was chosen "to ensure Figure 2
legend matches the text", and it records a computed fit of R² = 0.788 while the
manuscript reported 0.865 from a different script. The manuscript and the
deposited code therefore described two different networks.

**The two networks were built incomparably.** HF used a variance filter
(top 6,000 genes), a manual `cor → adjacency → TOM → cutreeDynamic` procedure,
and one soft threshold; HCC used all expressed genes, `blockwiseModules`, and a
different threshold. A difference in module composition between two networks
built by different methods cannot distinguish a biological difference from a
methodological one — so the "disease-specific network organisation" conclusion
the first submission drew from it was not supported by its own construction.
The revision rebuilds both networks under one frozen configuration on a common
gene universe.

**The variance filter removed the genes under study.** Selecting the top 6,000
genes by variance systematically drops ribosomal protein genes, which are
exactly the translation machinery the analysis is about. The revised pipeline
applies no variance pre-filtering.

**A reported test had no implementation.** The manuscript reported a binomial
test with *P* < 0.0001 for its mirror-pathway count. No such test exists
anywhere in this codebase, and the value the data actually give is *P* = 0.0135.
The revision replaces it with a sample-level permutation test.

**A factual claim contradicted the deposited table.** The manuscript stated that
"all 33 translation-related pathways showed positive effect sizes in HCC". The
submission's own `Table_S2_ssGSEA_Effect_Sizes.csv` contains six negative values
among those pathways.

**The supplementary materials were from an abandoned analysis.** The submitted
PDF was byte-identical to four other journals' submissions and internally
inconsistent — its Table S1 reported a 99-gene black module at β = 12 while the
manuscript reported a 227-gene green module at β = 17. The revision regenerates
every supplementary table and figure from the pipeline.

### Why none of these was caught before submission

The pre-submission check, `qa_plos.py`, read exactly two files — the manuscript
`.md` and its `.docx` — and asserted twelve properties of them: abstract length,
citation numbering, reference count, section presence, page size, margins, line
spacing, line numbers, footer page number. Every one of those is a property of
the manuscript text alone.

**Every defect listed above is invisible from inside the manuscript.** "All 33
pathways are positive" is refuted by a CSV. The hard-coded threshold is in a
script. The contradictory supplementary is a PDF the checker never opened. The
unimplemented test can only be found by grepping the code. A check confined to
the document cannot see any of them, however carefully it is run.

The response is not to be more careful — carefulness is not checkable. It is to
make verification cross-artifact by construction and to make running it a single
command. That is `verify_all.py`.

---

## Verifying the submission

```bash
python verify_all.py
```

Runs every check and exits non-zero if any fails. Each one reads at least one
file that is **not** the manuscript:

| Check | Reads |
|---|---|
| `Code/Revision_v2/R18_verify_manuscript.py` | the regenerated S-table CSVs; re-derives every number quoted in the text |
| `Code/Revision_v2/R20_verify_SI.py` | the built SI PDF; all 40 tables and 5 figures |
| `qa_plos_v2.py` | manuscript structure, citation order, page setup, and that withdrawn v1 wording has not crept back |
| `build_tracked_docx.py` | rebuilds the tracked docx and confirms accept-all reproduces v2, reject-all reproduces v1 |
| *(inline)* | main-figure PDF page size against the PLOS ceiling |

Before pushing, also run:

```bash
python pre_publish_check.py
```

This refuses to publish peer-review material. **This repository is public and is
also the working directory**, so submission packages sit in the same tree as the
code; a plain `git add -A` would publish the editor's decision letter and the
reviewers' comments.

<!-- pre-publish-check: policy-documentation -->
<!-- This file describes the exclusion policy. It quotes the phrase "decision
     letter" in prose and therefore trips the content scan; the sentinel above
     exempts it. It contains no peer-review material. -->

---

## Pipeline: `Code/Revision_v2/`

Runs in order. `R00_config.R` freezes every parameter and is `source()`d by
every script; the configuration hash is recorded alongside each output table.

| Script | Purpose |
|---|---|
| `R00_config.R`, `R00_io.R` | frozen configuration and I/O helpers |
| `R01_build_universe.R` | per-platform filtering, symbol mapping, common gene universe `U` |
| `R02_wgcna_coordinated.R` | both networks, one configuration, disease samples only |
| `R03_translation_module.R` | quantitative translation-module designation by ssGSEA |
| `R04_module_preservation.R` | `modulePreservation` in both directions |
| `R05_load_cohorts.R` | cohort loading and harmonisation |
| `R06_ssgsea_cohorts.R` | ssGSEA scoring across cohorts |
| `R07_effect_sizes.R` | Hedges *g* with *n*-weighted pooled SD and bootstrap CIs |
| `R08_redundancy.R` | Jaccard clustering, effective number of tests |
| `R09_mirror_test.R` | sample-level permutation mirror test |
| `R10_heterogeneity.R`, `R10b_hf_aetiology.R` | Cochran's *Q*, *I*², aetiology stratification |
| `R11_composition.R` | marker-based composition scores and adjusted correlations |
| `R12_tats_and_tf.R` | translation scores, survival, transcription-factor association |
| `R13_headline_numbers.R` | the numbers quoted in the manuscript |
| `R14_tom_robustness.R`, `R14b_cohesion.R` | connectivity-matched null, jackknife stability |
| `R15_cirrhosis_control.R` | same-organ non-malignant liver disease control |
| `R17_figures.R` | all main and supplementary figures, at PLOS page size |
| `R18_verify_manuscript.py` | manuscript numbers against regenerated tables |
| `R19_build_SI_pdf.py`, `R20_verify_SI.py` | build and verify the supporting information |

The `00_probe*.R` and `00_diag_*.R` scripts are exploratory and are not part of
the pipeline; they are kept because they record how the parameters were settled.

## Superseded pipeline

`Code/WGCNA/`, `Code/ssGSEA/`, `Code/TF_analysis/`, `Code/Figure_generation/`,
and the numbered scripts at the repository root (`00_MASTER_execution_guide.R`
through `19b_*.R`) — the first submission's analysis, retained for provenance.

**[`archive/README.md`](../archive/README.md) indexes all of it**, including why
the files stay at their current paths: they bind to those paths by absolute
reference, so relocating them would break the provenance the archive exists to
preserve.

## Data sources

All datasets are public; no new data were generated.

| Dataset | Description |
|---|---|
| GSE57338 | left ventricular myocardium (HF discovery, network construction) |
| GSE141198 | HCC tumours (network construction) |
| TCGA-LIHC | HCC, GDC portal (primary inference) |
| GSE14520, GSE76427 | independent HCC cohorts |
| GSE116250 | independent HF RNA-seq cohort |
| GSE141910 | cardiomyopathy collection |
| GSE89377 | multistep hepatocarcinogenesis, same-organ non-malignant control |

## Environment

R 4.6.0. Packages: WGCNA 1.74, DESeq2 1.52.0, GSVA 2.6.2, msigdbr 26.1.0,
clusterProfiler 4.20.0, survival 3.8-6, survminer 0.5.2, limma 3.68.4,
org.Hs.eg.db 3.23.1, TCGAbiolinks 2.40.0, dynamicTreeCut 1.63-1.

Python ≥ 3.9 with `python-docx` (submission assembly) and `pypdf`-free standard
library only for the verification scripts.

Seeds are fixed throughout (`set.seed(42)`): 10,000 iterations for the
sample-level permutation null, 2,000 for bootstrap resampling, 200 for module
preservation, 10,000 for the module-robustness nulls.
