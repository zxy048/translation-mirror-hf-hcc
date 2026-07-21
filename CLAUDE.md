# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R bioinformatics project for the manuscript: **"Translational Co-expression Programs Are Conserved but Exhibit Opposite Perturbation Patterns in Heart Failure and Hepatocellular Carcinoma"** (target: BBA – Molecular Basis of Disease, SCI Q1, IF 6–10).

Core finding: Translation co-expression architecture is conserved across HF and HCC but perturbed in opposite directions — upregulated in HCC, downregulated in HF. This is termed **"mirror regulation."**

GitHub: `https://github.com/zxy048/translation-mirror-hf-hcc`

## R Environment

- R 4.6.0 (2026-04-24 ucrt), Windows 11
- Project root: `D:/R_projects/revision_analysis`
- Key packages: WGCNA 1.74, DESeq2 1.52.0, GSVA 2.6.2, msigdbr 26.1.0, clusterProfiler 4.20.0, limma 3.68.4
- Working directory is set via `PROJ_DIR <- "D:/R_projects/revision_analysis"` at the top of each script

## Analysis Pipeline (numbered scripts)

Scripts are numbered by execution order. Run in RStudio with `source()`, not from CLI:

| Script | Purpose | Time |
|--------|---------|------|
| `01_download_GSE141198.R` | Download GSE141198 from GEO, parse ExpressionSet | ~5 min |
| `01b_parse_GSE141198.R` | Parse clinical variables from downloaded data | ~1 min |
| `01c_diagnose_GSE141198.R` | Diagnose data structure (read counts vs normalized) | ~1 min |
| `01d_load_readcounts.R` | Load raw read counts from supplementary files | ~1 min |
| `01e_ensembl_to_symbol.R` | Ensembl ID → gene symbol conversion | ~2 min |
| `01f_fix_vst.R` | VST normalization with correct parameters | ~2 min |
| `02_TGS_RPL39_validation.R` | TGS + RPL39 survival validation across 4 cohorts | 10–20 min |
| `02b_TGS_GSE141198.R` | TGS survival analysis in GSE141198 | ~5 min |
| `03_ssGSEA_cross_disease_pathway.R` | Cross-disease ssGSEA pathway comparison (original) | 30–60 min |
| `03b_ssGSEA_final.R` | ssGSEA using actual data with updated msigdbr API | 30–60 min |
| `03c_check_msigdbr.R` | Verify msigdbr API behavior | Quick |
| `04_GSE141198_WGCNA.R` | GSE141198 independent WGCNA (original) | 1–3 hr |
| `04b_GSE141198_WGCNA_final.R` | GSE141198 WGCNA final version (5,003 genes, β=4) | 1–3 hr |
| `05_TF_upstream_prediction.R` | Upstream TF correlation analysis (original) | 5–15 min |
| `05b_TF_upstream_prediction_final.R` | TF prediction final version (19 TFs) | 5–15 min |
| `06_direction_consistency_test.R` | Hub gene cross-disease direction consistency | <1 min |
| `07_final_figures_and_tables.R` | Assemble final figures and tables | Variable |
| `08_generate_remaining_figures.R` | Generate supplementary materials | Variable |
| `09_extract_clinical_characteristics.R` | Extract clinical characteristics for Table 1 | ~5 min |
| `10_generate_table_S1_GO_comparison.R` | GO enrichment comparison table | ~5 min |
| `11_generate_table_S1_final.R` | Final Table S1 generation | ~5 min |
| `00_MASTER_execution_guide.R` | Execution guide (readable, not executable) | — |

**Convention**: Scripts suffixed `b` (e.g., `03b`, `04b`, `05b`) are the final/production versions. Scripts suffixed `c`–`g` are diagnostics run once to debug data issues — not part of the main pipeline.

## Figure Generation Scripts

- `fig1_clean.R` / `fig6_clean.R` — Final versions of study design and mechanistic model figures
- `fix_soft_threshold_dpi.R` — Regenerates Figure S4A at 300 DPI (GSE141198, β=4, R²=0.84)
- `quick_HF_soft_threshold.R` — Regenerates Figure S4B at 300 DPI (GSE57338, 3,000 genes, β=12, R²=0.92)
- `regenerate_figures.R` — Batch figure regeneration
- `get_HF_WGCNA_params.R` — GSE57338 WGCNA parameter extraction (manual stepwise, 3,000 most variable genes)
- Output: `figures/` directory

## Key Data Files (RDS intermediates)

| File | Content |
|------|---------|
| `GSE141198_raw.rds` | Raw ExpressionSet from GEO |
| `GSE141198_pdata.rds` | Phenotype/clinical data |
| `GSE141198_counts_raw.rds` | Raw read counts |
| `GSE141198_counts_filt.rds` | Filtered counts (≥10 in ≥20% samples) |
| `GSE141198_vst.rds` | VST-normalized expression (DESeq2, blind=TRUE, nsub=1000) |
| `GSE141198_wgcna_input.rds` | 5,003 genes × 148 samples for WGCNA |
| `GSE141198_WGCNA_result.rds` | Complete WGCNA output |
| `GSE141198_clinical.rds` | Clinical annotations |
| `tgs_GSE141198_result.rds` | TGS survival analysis results |
| `ssgsea_cross_disease_result.rds` | ssGSEA pathway activity scores |
| `TF_upstream_result.rds` | TF-TGS correlation results |

## Manuscript Files

- **`Manuscript_EN_v1_draft.md`** — Primary English draft (the authoritative version). Contains the full manuscript with embedded figure/table placeholders. All Steps 8–13 revision edits were applied here.
- `Manuscript_CN_v8_draft.md` / `Manuscript_CN_v9_draft.md` — Earlier Chinese drafts (v8 is pre-revision, v9 is post-revision).
- `Manuscript_Revision_Plan.md` — The revision framework: narrative shift from "shared signals" to "mirror regulation," evidence matrix, figure/table plan, journal targeting, and anticipated reviewer responses.
- `Cover_Letter_BBA_MBD.md` — Cover letter for BBA-MBD submission.
- `Highlights.md` — 3–5 bullet points for journal submission.

## Manuscript Revision Audit Notes (Steps 8–13 applied to v1)

When editing the manuscript, follow these conventions:

### Terminology (three-tier system)
- **Tier 1** (WGCNA technical): "translation co-expression module"
- **Tier 2** (Results): "translation co-expression program"
- **Tier 3** (Abstract/Discussion/Conclusions): "translation co-expression architecture"

### Banned words
- **AI trace words** (never use): "Interestingly", "Importantly", "Notably", "markedly", "Collectively", "Together" (as discourse marker), "Remarkably", "Strikingly", "Profoundly"
- **Overclaim verbs** (avoid): "demonstrate", "dictated by", "drive", "is characterized by", "was demonstrated"

### Brand term
- **"mirror regulation"** — Must appear in Discussion P1 sentence 2 and in the Conclusions final sentence. Describes conserved architecture deployed in opposite directions by disease context.

### Discussion structure (5 paragraphs, strictly segregated)
- P1: Mirror regulation definition + literature comparison
- P2: ATF4/MYC mechanism
- P3: Methodology (pathway-level vs gene-level resolution, hub gene independence)
- P4: Clinical boundary (TGS not prognostic)
- P5: Limitations → separate paragraph for broader significance

### Abstract vs Conclusions distinction
- Abstract = "What did we find?" (no mirror regulation)
- Conclusions = "Why does it matter?" (mirror regulation closure)

## Parameter Consistency (critical for reviewer scrutiny)

- **GSE57338 WGCNA**: 3,000 most variable genes, β=12, R²=0.92, manual stepwise method (not `blockwiseModules`)
- **GSE141198 WGCNA**: 5,003 genes (all expressed after filtering), β=4, R²=0.84, `blockwiseModules` one-step
- **Soft threshold figures**: Both S4A and S4B use `SFT.R.sq >= 0.85` for red coloring, with dashed line at 0.85. Gene counts in figure scripts must match WGCNA scripts.
- **ssGSEA**: 81 pathways total (48 Hallmark + 1 KEGG ribosome + 32 Reactome translation), of which 33 are translation/ribosome-related
- **TGS**: 6 genes (EEF1A1, FAU, RPL39, RPL3, RPL32, RPL41) — RPS28 excluded due to opposite cross-disease directionality

## Dataset Inclusion Criteria

All datasets must satisfy: (i) disease + control samples within same study, (ii) sample size ≥ 100, (iii) raw/normalized expression matrices publicly available.

| Dataset | Disease | n | Platform |
|---------|---------|---|----------|
| GSE57338 | HF | 313 LV samples | GPL11532 Affymetrix HuGene 1.1 ST |
| TCGA-LIHC | HCC | 424 (371 tumor + 50 normal) | Illumina RNA-seq |
| GSE141198 | HCC | 148 tumors (94 events) | RNA-seq |
| GSE14520 | HCC | 221 tumors | GPL3921 Affymetrix U133A |
| GSE76427 | HCC | 115 tumors | GPL10558 Illumina HumanHT-12 V4.0 |
