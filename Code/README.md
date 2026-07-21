# Code Repository: Translation Mirror Regulation in HF and HCC

This repository contains the complete analysis pipeline for:

**"Conserved Translational Co-expression Programs Exhibit Mirror Regulation Across Heart Failure and Hepatocellular Carcinoma"**

---

## Directory Structure

```
Code/
├── WGCNA/
│   ├── 00_download_GSE141198.R          # Download GSE141198 from GEO
│   ├── 00b_parse_GSE141198.R            # Parse GSE141198 expression matrix
│   ├── 00c_load_readcounts.R            # Load raw read counts
│   ├── 00d_ensembl_to_symbol.R          # Convert Ensembl IDs to gene symbols
│   ├── 00e_fix_vst.R                    # VST normalization
│   ├── 01_HCC_WGCNA.R                   # GSE141198 signed WGCNA (β=4, R²=0.84)
│   ├── 02_HF_WGCNA.R                    # GSE57338 signed WGCNA (β=12, R²=0.79)
│   └── 03_HF_WGCNA_parameter_selection.R # Soft-threshold selection for HF
├── ssGSEA/
│   └── 01_ssGSEA_cross_disease.R        # Parallel ssGSEA: 81 pathways, TCGA-LIHC + GSE57338
├── TF_analysis/
│   ├── 01_TF_TGS_correlation.R          # TF-TGS Spearman correlation (19 TFs × TCGA-LIHC)
│   ├── 02_direction_consistency.R       # Cross-disease hub gene direction analysis
│   ├── 03_TGS_RPL39_validation.R        # TGS/RPL39 validation in GSE14520 + GSE76427
│   └── 04_TGS_GSE141198_validation.R    # TGS survival validation in GSE141198
└── Figure_generation/
    ├── 01_main_figures.R                # Figures 2–5 (WGCNA, ssGSEA, TGS, TF)
    ├── 02_supplementary_figures.R       # Figures S1–S4
    ├── 03_fig1_study_design.R           # Figure 1: Study design overview
    ├── 04_fig6_mechanistic_model.R      # Figure 6: Mirror regulation framework
    ├── 05_graphical_abstract.R          # Graphical Abstract
    ├── 06_build_manuscript_docx.py      # Build submission-ready .docx from Markdown
    └── 07_table_S1_module_comparison.R  # Table S1: Module comparison + GO enrichment
```

## Execution Order

A complete reproduction run proceeds in this sequence:

```
01. Data Preprocessing
     → Code/WGCNA/00_download_GSE141198.R
     → Code/WGCNA/00b_parse_GSE141198.R
     → Code/WGCNA/00c_load_readcounts.R
     → Code/WGCNA/00d_ensembl_to_symbol.R
     → Code/WGCNA/00e_fix_vst.R

02. WGCNA (Co-expression Network Construction)
     → Code/WGCNA/03_HF_WGCNA_parameter_selection.R   # β selection
     → Code/WGCNA/02_HF_WGCNA.R                        # HF black module (99 genes)
     → Code/WGCNA/01_HCC_WGCNA.R                       # HCC blue module (1,315 genes)

03. ssGSEA (Pathway Activity Scoring)
     → Code/ssGSEA/01_ssGSEA_cross_disease.R           # 81 pathways × 2 diseases

04. TF Analysis (Upstream Regulator Inference)
     → Code/TF_analysis/01_TF_TGS_correlation.R         # ATF4/ISR → TGS
     → Code/TF_analysis/02_direction_consistency.R      # Hub gene cross-disease direction
     → Code/TF_analysis/03_TGS_RPL39_validation.R       # External cohort validation
     → Code/TF_analysis/04_TGS_GSE141198_validation.R   # GSE141198 survival

05. Figures & Tables
     → Code/Figure_generation/01_main_figures.R
     → Code/Figure_generation/02_supplementary_figures.R
     → Code/Figure_generation/07_table_S1_module_comparison.R
     → Code/Figure_generation/06_build_manuscript_docx.py   # Final .docx assembly
```

## Dependencies

**R ≥ 4.3.0** with packages:

| Package | Version used | Purpose |
|---------|-------------|---------|
| WGCNA | 1.74 | Signed co-expression network construction |
| DESeq2 | 1.52.0 | TCGA-LIHC differential expression |
| GSVA | 2.6.2 | ssGSEA pathway scoring |
| msigdbr | 26.1.0 | MSigDB gene set retrieval |
| clusterProfiler | 4.20.0 | GO enrichment analysis |
| survival | 3.8-6 | Kaplan–Meier survival analysis |
| survminer | 0.5.2 | Survival visualization |
| limma | 3.68.4 | Microarray differential expression |
| org.Hs.eg.db | 3.23.1 | Gene annotation database |
| metafor | 5.0.1 | Meta-analysis |
| ggplot2 | ≥ 3.5.0 | Figure generation |

**Python ≥ 3.9** with:
- `python-docx` — document assembly

## Data Sources

All datasets are publicly available:
- **GSE57338**: GEO (HF myocardium, n = 313)
- **GSE141198**: GEO (Taiwan HCC, n = 148)
- **GSE14520**: GEO (HCC validation, n = 244)
- **GSE76427**: GEO (HCC validation, n = 115)
- **TCGA-LIHC**: GDC Data Portal (HCC, clinical n = 377, expression n = 424)

## Key Parameters

| Parameter | GSE57338 (HF) | GSE141198 (HCC) |
|-----------|---------------|-----------------|
| Soft threshold (β) | 12 | 4 |
| Signed R² | 0.79 | 0.84 |
| Translation module | Black (99 genes) | Blue (1,315 genes) |
| Hub genes in module | 4/7 | 4/7 |

## Reproducibility

- Random seed fixed at `set.seed(42)` for all permutation tests (10,000 iterations)
- Bootstrap CIs use 2,000 resamples
- Session information recorded in `sessionInfo.txt`

## Contact

For questions regarding the analysis code, please open an issue on the GitHub repository:
https://github.com/zxy048/translation-mirror-hf-hcc
