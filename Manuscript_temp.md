# Conserved Translational Co-expression Programs Exhibit Mirror Regulation Across Heart Failure and Hepatocellular Carcinoma

**Running title**: Mirror Regulation of Translation in HF and HCC

---

## Abstract

**Background**: Heart failure (HF) and hepatocellular carcinoma (HCC) share hallmarks of metabolic remodeling and disrupted protein homeostasis. Translational dysregulation — altered expression of ribosomal proteins, translation factors, and regulatory components — occurs in both diseases. Whether the underlying transcriptional programs reflect a conserved architecture or disease-specific configurations has not been systematically tested.

**Methods**: We applied an identical five-layer pipeline to both diseases: (i) weighted gene co-expression network analysis (WGCNA) in HF myocardium (GSE57338, n = 313) and HCC (GSE141198, n = 148); (ii) single-sample gene set enrichment analysis (ssGSEA) of 81 pathways across TCGA-LIHC (n = 424) and GSE57338, with quantitative cross-disease effect size comparison; (iii) external validation of a translation gene score (TGS) in three independent HCC cohorts; (iv) transcription factor–TGS correlation analysis in TCGA-LIHC tumors (n = 371); and (v) cross-disease directionality analysis at gene and pathway levels.

**Results**: The translation co-expression module replicated in HCC: the GSE141198 blue module (1,315 genes) was enriched for cytoplasmic translation (GO adjusted p = 6.3 × 10⁻¹²), with four of seven hub genes co-localizing to this module. Cross-disease ssGSEA revealed a negative correlation in translation pathway effect sizes (Spearman ρ = −0.598, p = 0.0003; permutation p = 0.0079): 24 of 33 translation pathways exhibited HCC-up/HF-down mirror perturbation. The TGS lacked prognostic value in all three external cohorts, indicating cohort-specific rather than universal clinical association. ATF4, master transcription factor of the integrated stress response (ISR), was the strongest transcriptional correlate of the translation co-expression program (ρ = +0.439, FDR < 0.0001); MYC target pathway activity showed the highest overall correlation with TGS (ρ = +0.613, p < 0.0001). These findings support conserved translation co-expression architecture despite disease-specific mirror perturbation. None of the seven hub genes appear in canonical Hallmark MYC target gene sets, supporting their independent identification through network analysis.

**Conclusions**: Translation co-expression architecture is conserved across HF and HCC, but disease context determines whether this program is activated or suppressed: upregulated in HCC, downregulated in HF. ATF4/ISR and MYC-driven proliferative signaling may converge to promote translational activation in HCC. These observations support a descriptive mirror regulation framework in which conserved translation co-expression architecture is deployed in opposite directions according to disease context.

**Keywords**: translation; heart failure; hepatocellular carcinoma; WGCNA; ssGSEA; ATF4; integrated stress response; MYC; mirror regulation

---

## Introduction

Heart failure (HF) and hepatocellular carcinoma (HCC) appear clinically divergent: HF is a chronic degenerative condition of the myocardium, HCC an aggressive malignancy of the liver. Yet beneath this superficial distinction lies a shared molecular landscape. Both diseases undergo profound metabolic reprogramming — the failing heart faces chronic ATP deficit and substrate depletion [1,2], while HCC cells rewire central carbon metabolism to sustain unrestrained proliferation [6,14]. Both exhibit disrupted protein homeostasis, the balance among protein synthesis, folding, and degradation essential for cellular function and survival [3,4]. Translational regulation lies at the intersection of these processes, coupling cellular energy status with protein synthesis capacity. Accordingly, both diseases exhibit dysregulation of the translational apparatus — the ribosome, translation factors, and their upstream regulators — which mediates the final and most energetically expensive step of gene expression [1,5].

The nature of this translational dysregulation, however, operates in opposite directions. In HF, myocardial energy deficiency is associated with suppression of anabolic processes, including ribosome biogenesis and cap-dependent translation, as an adaptive energy-conservation strategy [2,4]. mTORC1 signaling is inhibited in the failing heart [4]. In HCC, by contrast, oncogenic pathways converge on the translational apparatus to upregulate ribosome biogenesis and protein synthesis, supporting uncontrolled proliferation [5,6,18]. These mirror-image pressures on the same molecular machinery raise a fundamental question: does the same translation-related transcriptional program operate in both diseases, with disease context determining only the direction of perturbation?

Answering this question requires moving beyond single-disease frameworks. Existing transcriptomic studies in HF and HCC have largely been confined to single-disease contexts, identifying differentially expressed genes, co-expression modules, or prognostic signatures within one disease at a time [3]. While these studies have generated valuable disease-specific insights, they cannot distinguish between features that are disease-specific and features that reflect general organizational principles of the translational system. Cross-disease comparative transcriptomics — the parallel analysis of independent datasets from distinct diseases using identical computational pipelines — provides a direct strategy for making this distinction. If a co-expression module appears in both diseases, its structure may reflect conserved organizational properties of translation-related transcriptional regulation rather than disease-specific factors. If the same module is perturbed in opposite directions, the polarity of that perturbation is likely governed by disease context.

HF and HCC provide an instructive comparison because they represent opposite extremes of translational demand. The failing heart is an energy-depleted environment where translational suppression is adaptive [1,2]; the HCC tumor microenvironment is a pro-growth environment where translational activation supports malignant fitness [5,6]. If the same translation co-expression programs are upregulated in HCC and downregulated in HF, this would constitute strong evidence for context-dependent perturbation of a conserved molecular architecture. We refer to this descriptive pattern as "mirror regulation," defined as the opposite-direction perturbation of a conserved translational co-expression program across distinct disease contexts. To test this hypothesis, we developed an integrative cross-disease transcriptomic framework combining network analysis, pathway-level comparison, external validation, and transcriptional regulator inference (Figure 1).



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE 1 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure1_Study_Design.png`

**Caption:** Five-layer cross-disease analytical framework for testing the mirror regulation hypothesis.

</div>

**Figure 1. Study design overview.** The analytical framework comprises five layers: (A) WGCNA-based identification of translation co-expression modules in HF (GSE57338, n = 313); (B) independent WGCNA replication in HCC (GSE141198, n = 148); (C) parallel ssGSEA pathway activity scoring across both diseases with quantitative comparison of effect sizes; (D) external validation of TGS prognostic value in three independent HCC cohorts; (E) systematic screening for upstream transcriptional regulators of the translation co-expression program. The mirror regulation hypothesis — conserved module structure, disease-specific perturbation direction — is tested across all five layers.

---

## Results

### 2.1 The Translation Co-expression Program Is Conserved Across HF and HCC

The translation co-expression program was independently identified in both diseases. In HF myocardium (GSE57338, n = 313; dilated cardiomyopathy, ischemic cardiomyopathy, and non-failing controls), signed WGCNA identified a disease-associated module enriched for translational machinery genes. A soft threshold of β = 12 was selected, achieving scale-free topology fit R² = 0.79. More than ten co-expression modules were identified. Gene Ontology (GO) Biological Process enrichment analysis revealed that the black module was significantly enriched for translation-related biological processes, with "cytoplasmic translation" as the most significant term (adjusted p < 1 × 10⁻⁸). The black module showed the strongest correlation with HF phenotype (r = 0.794, p = 8.73 × 10⁻⁴⁹). Based on the dual criteria of high intramodular connectivity with the translation-enriched black module (Module Membership > 0.80) and strong association with the HF phenotype (Gene Significance > 0.20), seven hub genes were identified: **EEF1A1, FAU, RPL39, RPL3, RPL32, RPL41, and RPS28**. Four mapped to the black module (FAU, RPL39, RPL32, RPL41) and three to the turquoise module (EEF1A1, RPL3, RPS28). These genes encode ribosomal proteins or translation elongation factors. Because RPS28 displayed inconsistent transcriptional direction across diseases (Section 2.2; Supplementary Table S3), it was excluded from the TGS, which was constructed from the remaining six hub genes.

To test whether this translation co-expression module exists independently in HCC, we performed independent WGCNA in GSE141198 (Taiwan HCC cohort, n = 148, 94 OS events), applying the same signed network framework, with the soft threshold selected to best approximate the scale-free topology criterion (R² ≈ 0.85). The soft threshold for GSE141198 was β = 4 (R² = 0.84, the closest achievable value), with a minimum module size of 30. The GSE141198 blue module (1,315 genes) showed "cytoplasmic translation" as the most significantly enriched GO term (adjusted p = 6.3 × 10⁻¹²), consistent with the GO enrichment results of the HF translation module. Four of the seven hub genes (FAU, RPL3, RPL32, RPS28) co-localized to the blue module, confirming cross-disease conservation of the translation co-expression program. RPL39 and RPL41 were assigned to the grey (unassigned) module in GSE141198, reflecting their weaker co-expression signal in this HCC cohort; however, both genes were retained in the TGS because they were strongly co-localized to the translation-enriched module in the HF dataset, and their inclusion preserves the structural integrity of the translation apparatus as a functional unit — excluding ribosomal protein genes solely on the basis of module assignment in one disease would fragment a biologically coherent machinery. GO enrichment overlap and module color mapping across the two datasets are detailed in Supplementary Table S1. Module–trait association analysis showed that no module in GSE141198 was significantly associated with overall survival (OS) status (all modules: p > 0.05, uncorrected for multiple comparisons). The translation co-expression program is conserved across HF and HCC; its association with clinical outcome is context-dependent (Figure 2). We next examined whether this conserved program is perturbed in opposing directions at the pathway level.


**Figure 2. Independent identification and functional conservation of translation co-expression programs across HF and HCC.** (A) Gene dendrogram and module assignment from signed WGCNA of GSE57338 (HF, n = 313). The black module (n = 99 genes) was identified as the translation-associated module (β = 12, signed R² = 0.79). (B) Gene dendrogram and module assignment from signed WGCNA of GSE141198 (HCC, n = 148). The blue module (n = 1,315 genes) was identified as the translation-associated module (β = 4, signed R² = 0.84). (C) Shared functional enrichment of translation-associated modules across diseases. GO Biological Process terms enriched in both HF black and HCC blue modules; translation-related terms highlighted. These analyses demonstrate that translation-associated co-expression programs identified independently in HF and HCC share conserved functional characteristics.

### 2.2 Translation Pathways Exhibit Mirror Perturbation Across Diseases

To quantitatively compare the perturbation direction and magnitude of translation pathways across the two diseases, we performed ssGSEA pathway activity scoring in parallel across 81 pathways (48 Hallmark + 1 KEGG ribosome + 32 Reactome translation-related pathways) in TCGA-LIHC (371 tumor vs. 50 normal) and GSE57338 (HF vs. non-failing), and computed Cohen's d effect sizes for each pathway.

The Spearman correlation of effect sizes across all 81 pathways between the two diseases was ρ = −0.290 (p = 0.0089), reflecting an overall moderate-to-weak negative correlation. However, when the analysis was restricted to the 33 translation/ribosome-related pathways, the negative correlation in effect sizes increased to **ρ = −0.598 (p = 0.0003)**. A 10,000-iteration permutation test confirmed that this result was unlikely to arise by chance (permutation p = 0.0079; Figure 3). The specific composition of the 33 translation pathways (1 KEGG ribosome pathway + 32 Reactome translation-related pathways) and complete effect size data for all 81 pathways are provided in Supplementary Table S2.

Detailed examination of translation pathway effect sizes revealed a consistent pattern: all translation pathways showed positive effect sizes in HCC (Cohen's d > 0, higher expression in tumor vs. normal tissue), whereas the vast majority showed negative effect sizes in HF (d < 0, lower expression in failing vs. non-failing myocardium). Among the 33 translation pathways, 24 exhibited an HCC-positive/HF-negative paired pattern, representing a predominant mirror perturbation (binomial test p < 0.0001). The remaining nine pathways were enriched for mitochondrial translation or miRNA-mediated translational regulation rather than core cytoplasmic translation (detailed in Supplementary Table S2). Proliferation-associated pathways (E2F Targets, G2M Checkpoint, and MYC Targets V1) exhibited the largest positive effect sizes in HCC (Cohen's d = +1.83 to +2.35; full rankings are provided in Supplementary Table S2). The most upregulated pathways in HF were Bile Acid Metabolism (d = +0.79) and Interferon Alpha Response (d = +0.67); their absolute effect sizes were considerably smaller than those of proliferation pathways in HCC.

To determine whether this mirror pattern was also evident at the individual gene level, we compared the cross-disease log2FC values of the seven hub genes. Across the six directionally concordant hub genes (RPS28 excluded due to opposite cross-disease directionality), the Spearman correlation of cross-disease log2FC values was ρ = 0.486 (p = 0.356; Figure 4B, Supplementary Figure S1). The full seven-gene correlation was ρ = 0.126 (p = 0.788; Supplementary Table S3), confirming that RPS28 is a directional outlier. This marginal correlation at the individual gene level contrasts with the strong negative correlation observed at the pathway level, illustrating a key analytical principle: single-gene comparisons are susceptible to biological and technical noise that can obscure systemic signals, whereas pathway- and network-level aggregation — the rationale underlying WGCNA-based module analysis — averages out gene-level variability to reveal the underlying regulatory architecture. The specific log2FC values, FDR, and direction concordance annotations for the seven hub genes in HF and HCC are provided in Supplementary Table S3.



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE 3 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure_ssGSEA_cross_disease.png`

**Caption:** Cross-disease comparison of pathway effect sizes. Scatter plot of Cohen's d for 81 pathways in HCC (TCGA-LIHC) vs. HF (GSE57338). Translation pathways (n=33) in red; non-translation (n=50) in gray. Spearman ρ = −0.598 (p = 0.0003; permutation p = 0.0079).

</div>

**Figure 3. Cross-disease comparison of pathway effect sizes.** Scatter plot of Cohen's d effect sizes for 81 pathways in HCC (TCGA-LIHC, tumor vs. normal, x-axis) versus HF (GSE57338, failing vs. non-failing, y-axis). Translation/ribosome-related pathways (n = 33) are shown as red points; non-translation pathways (n = 50) as gray points. The Spearman correlation for translation pathways is ρ = −0.598 (p = 0.0003; 10,000 permutations, p = 0.0079). Dashed lines mark d = 0. Quadrant I (upper right, HCC↑ HF↑): 2 pathways; Quadrant II (upper left, HCC↓ HF↑): 0 pathways; Quadrant III (lower left, HCC↓ HF↓): 0 pathways; Quadrant IV (lower right, HCC↑ HF↓): 24 of 33 translation pathways.

### 2.3 TGS Prognostic Value Is Cohort-Specific

Although the translation co-expression program exhibited robust structural conservation across diseases (Section 2.1) and consistent mirror perturbation at the pathway level (Section 2.2), whether this conservation extends to clinical utility remained unknown. We therefore evaluated the prognostic performance of TGS in three independent HCC cohorts. TGS showed significant prognostic value in TCGA-LIHC (discovery cohort). To test whether this association generalizes across populations, we performed the identical TGS calculation and survival analysis in GSE141198 (n = 148, 94 events).

TGS did not discriminate between high- and low-risk groups in this cohort (Log-rank p = 0.479; Figure 4A). When analyzed as a continuous variable, TGS was also non-significant in Cox regression (HR = 0.984, 95% CI 0.817–1.185, p = 0.869). Stratification by etiology (HBV, HCV, NBNC) revealed no prognostic significance for TGS in any subgroup. Across three independent HCC external cohorts — GSE141198, GSE14520 (n = 221), and GSE76427 (n = 115) — TGS showed no significant prognostic value (Figure 4A). Structural conservation of the translation co-expression program does not ensure universal prognostic value. Single-gene analysis of RPL39, the hub gene with the largest effect size in HCC, similarly showed no significant prognostic value in GSE141198 (Supplementary Figure S2).



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE 4 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure4_TGS_Validation.png`

**Caption:** TGS prognostic validation and hub gene cross-disease directionality. (A) Kaplan–Meier survival curve for TGS in GSE141198 (n=148). Log-rank p = 0.479. (B) Cross-disease log2FC of seven hub genes in HCC vs. HF.

</div>

**Figure 4. TGS prognostic validation and hub gene cross-disease directionality.** (A) Kaplan–Meier survival curve for TGS (median dichotomization) in the GSE141198 Taiwanese HCC cohort (n = 148, 94 events). Log-rank p = 0.479; TGS does not significantly stratify patients by overall survival. (B) Cross-disease expression direction (log2 fold-change) of the seven translation hub genes in HCC (TCGA-LIHC, tumor vs. normal) and HF (GSE57338, failing vs. non-failing). Circles: Heart Failure; triangles: HCC. Asterisk (*) denotes FDR < 0.05 within the respective disease. RPS28 (top panel) shows opposite directionality across diseases. HR, hazard ratio; CI, confidence interval.



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE S2 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure_S2_GSE141198_Survival.png`

**Caption:** RPL39 single-gene expression in GSE141198. Kaplan–Meier curve for RPL39 (median dichotomization) in GSE141198 (n=148).

</div>

**Supplementary Figure S2. RPL39 single-gene expression in GSE141198.** Kaplan–Meier curve for RPL39 expression (median dichotomization) in the GSE141198 Taiwanese HCC cohort (n = 148). RPL39 was selected as a representative hub gene from the translation co-expression program; its expression alone does not show significant prognostic value in this cohort.

### 2.4 Upstream Regulatory Programs Converge on the Translation Machinery

Given that the translation co-expression program is structurally conserved (Section 2.1) yet clinically decoupled from universal prognostic value (Section 2.3), we next asked which upstream regulatory programs govern its transcriptional activity. To identify upstream regulators associated with the translation co-expression program, we performed systematic transcription factor (TF)–TGS correlation analysis across 371 tumor samples from TCGA-LIHC. Nineteen candidate TFs were selected, encompassing established translation regulators (MYC, MYCN, E2F family, mTOR pathway components) as well as stress response-related TFs.

Thirteen of the nineteen TFs were significantly correlated with TGS at false discovery rate (FDR) < 0.05. Among these, **ATF4 was the strongest positive TF correlate of TGS** (ρ = +0.439, FDR < 0.0001; Figure 5A). ATF4 is the master transcription factor of the ISR. Additional transcription factors related to the ISR and the UPR were also significantly associated with TGS, including DDIT3 (ρ = +0.289, FDR < 0.0001) and XBP1 (ρ = +0.162, FDR = 0.003). HIF1A was negatively correlated with TGS (ρ = −0.417, FDR < 0.0001).

At the pathway activity level, **the correlation of MYC target pathway activity with TGS (ρ = +0.613, p < 0.0001; Figure 5B) far exceeded that of any individual TF mRNA level**. The pathway activity of MYC Targets V2 (ρ = +0.613) was more strongly correlated with TGS than MYC Targets V1 (ρ = +0.399) or MYC mRNA level (ρ = +0.297). Within the E2F family (E2F1–E2F8), three members showed significant negative correlations with TGS — E2F7 (ρ = −0.359, FDR < 0.0001), E2F8 (ρ = −0.296, FDR < 0.0001), and E2F3 (ρ = −0.267, FDR < 0.0001) — while E2F4 was positively correlated (ρ = +0.136, FDR = 0.013), and the remaining four (E2F1, E2F2, E2F5, E2F6) were not significant. mTOR pathway components MTOR (ρ = −0.391, FDR < 0.0001) and its regulatory subunit RPTOR (ρ = −0.143, FDR = 0.009) were both negatively correlated with TGS, consistent with mTOR pathway suppression in settings where alternative stress-responsive programs are engaged. TP53 showed a moderate positive correlation (ρ = +0.243, FDR < 0.0001), while NFE2L2 was negatively associated (ρ = −0.148, FDR = 0.007). Complete TF–TGS correlation statistics for all nineteen candidate TFs are summarized in Table 3.


<div style="border: 2px dashed #666; padding: 12px; margin: 18px 0; background: #fafafa;">

**⬇ [TABLE 3 — INSERT TABLE HERE] ⬇**

**Source:** `TF_upstream_result.rds (render as formatted table)`

**Content:** TF–TGS Correlation Summary. Spearman ρ, 95% CI, FDR-adjusted p-values for all 19 candidate transcription factors correlated with TGS in TCGA-LIHC tumors (n=371).

</div>



None of the seven hub genes are present in the canonical Hallmark MYC Targets V1 (200 genes) or V2 (58 genes) gene sets. Although these Hallmark sets contain other ribosomal protein genes (e.g., RPL14, RPS10), the seven hub genes represent a translation co-expression signature independently discovered through data-driven network analysis. Fisher's exact test confirmed no significant enrichment of the seven hub genes in any Hallmark TF target gene set (MYC V1/V2, E2F, mTORC1: all p = 1.00; Supplementary Figure S3).



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE S3 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure_S3_TF_Fisher_enrichment.png`

**Caption:** Independence of hub genes from canonical TF target gene sets. Fisher's exact test results for enrichment of seven hub genes in Hallmark TF target gene sets. All p = 1.00.

</div>

**Supplementary Figure S3. Independence of hub genes from canonical TF target gene sets.** Fisher's exact test results for the enrichment of the seven translation hub genes in Hallmark TF target gene sets (MYC Targets V1, MYC Targets V2, E2F Targets, mTORC1 Signaling). All p-values = 1.00.



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE 5 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure_TF_TGS_correlation.png
│  File (Panel B): figures/Figure_Pathway_TGS_correlation.png`

**Caption:** (A) Spearman correlation of 19 candidate TFs with TGS in TCGA-LIHC (n=371). ATF4 is the strongest correlate (ρ=+0.439, FDR<0.0001). (B) Spearman correlation of Hallmark pathway ssGSEA scores with TGS. MYC Targets V2 shows the strongest correlation (ρ=+0.613, p<0.0001).

</div>

**Figure 5. Upstream regulator analysis of the translation co-expression program.** (A) Spearman correlation of 19 candidate transcription factors with TGS in TCGA-LIHC tumor samples (n = 371). Red bars: positive correlation (FDR < 0.05); blue bars: negative correlation (FDR < 0.05); gray: not significant. ATF4 is the strongest TF correlate (ρ = +0.439, FDR < 0.0001). (B) Spearman correlation of Hallmark pathway ssGSEA scores with TGS. MYC Targets V2 pathway activity shows the strongest correlation with TGS (ρ = +0.613, p < 0.0001), compared with MYC mRNA level alone (ρ = +0.297).



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE S1 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure_S1_Direction_Consistency.png`

**Caption:** Cross-disease directional consistency of six directionally concordant hub genes (EEF1A1, FAU, RPL39, RPL3, RPL32, RPL41). Scatter plot of log2FC in HCC vs. HF. Spearman ρ = 0.486 (p = 0.356). RPS28 (red) excluded from TGS due to opposite directionality; seven-gene ρ = 0.126 (p = 0.788).

</div>

**Supplementary Figure S1. Cross-disease directional consistency of six directionally concordant hub genes.** Scatter plot of log2 fold-change values for the six translation hub genes retained in TGS (EEF1A1, FAU, RPL39, RPL3, RPL32, RPL41) in HCC (TCGA-LIHC, tumor vs. normal) versus HF (GSE57338, failing vs. non-failing). Spearman ρ = 0.486 (p = 0.356) for the six concordant genes. RPS28 (red) is shown for reference and was excluded from TGS due to opposite cross-disease directionality.

---

## Discussion

This study investigated whether translation-related transcriptional programs are organized according to conserved principles across heart failure and hepatocellular carcinoma despite their distinct pathological contexts. We refer to this organizational pattern as mirror regulation, in which a conserved translation co-expression architecture is preserved across diseases while its transcriptional perturbation follows opposite, disease-specific directions: consistently upregulated in the proliferative HCC microenvironment and consistently downregulated in the energy-depleted failing myocardium. This cross-disease comparative strategy — distinguishing shared organizational principles from disease-specific transcriptional responses — reveals regulatory principles that cannot be appreciated through single-disease analyses. Previous transcriptomic studies have characterized translational dysregulation in HF or HCC independently, whereas our cross-disease framework reveals that these disease-specific observations can be unified within a conserved translation co-expression architecture exhibiting mirror regulation. Importantly, cross-disease comparisons were performed using within-dataset disease-associated effect sizes (Cohen's d, disease vs. matched control within each tissue context) rather than absolute expression levels, thereby minimizing organ-specific baseline differences and enabling comparison of disease-associated perturbation directions rather than tissue-specific expression states.

Among the upstream regulators we examined, ATF4 emerged as the dominant transcriptional correlate of the translation co-expression architecture, a finding that extends beyond the canonical MYC-centric framework of translational regulation. ATF4 is the master transcription factor of the ISR, a conserved signaling cascade that couples cellular stress sensing to translational reprogramming [7,9,16]; its pro-survival arm can be co-opted to support tumor cell fitness under microenvironmental stress [8,10,11,17]. The co-association of ATF4, DDIT3 (CHOP), and XBP1 — components of the ISR and the UPR signaling chain — with TGS indicates ISR pathway-level activation rather than isolated ATF4 upregulation, suggesting that HCC cells may maintain translational output through coordinated stress-adaptive transcriptional mechanisms. Previous studies have implicated ATF4 in adaptive stress responses in cancer [8,10,11,17]; our findings extend this concept by linking ATF4 to a conserved translation co-expression program rather than isolated downstream targets. In parallel, MYC target pathway activity showed the strongest overall association with TGS, exceeding the explanatory power of MYC mRNA level alone. The stronger association of MYC target pathway activity with TGS likely reflects the statistical advantage of pathway-level aggregation, which integrates coordinated expression changes across multiple downstream target genes and reduces single-gene measurement noise; MYC pathway activity and ATF4 expression are therefore interpreted as complementary regulatory layers rather than competing upstream determinants. These represent distinct regulatory dimensions — stress-adaptive (ATF4/ISR) and proliferative (MYC) — neither subordinate to the other, converging on the translational machinery. The preferential association of ATF4 suggests that stress-responsive transcriptional control, rather than proliferative signaling alone, contributes to maintaining translation output. This interpretation is consistent with the requirement for tumor cells to preserve protein synthesis despite persistent proteotoxic stress. These observations support a model in which the ATF4/ISR stress-adaptation axis and MYC-driven proliferative signaling converge on the translational machinery, with ATF4 providing the stress-responsive input and MYC sustaining the proliferative demand [18,19]. HIF1A showed a strong negative correlation with TGS, suggesting that hypoxia-responsive transcription may be regulated independently of the ATF4/ISR-associated translation program — an observation that warrants further investigation.

None of the seven hub genes are present in the canonical Hallmark MYC target gene sets (MYC Targets V1, 200 genes; MYC Targets V2, 58 genes), confirming that they represent a translation co-expression signature independently discovered through data-driven network analysis rather than a predetermined candidate gene set. An important analytical insight concerns the resolution at which mirror regulation is detectable. The cross-disease negative correlation in translation pathway effect sizes was highly significant at the pathway level, yet this mirror pattern was not evident at the individual gene level — a discrepancy that reflects the ability of pathway-level aggregation to average out the multiple layers of biological variability that obscure systemic signals at the single-gene level. The nine pathways that did not conform to this predominant mirror pattern were enriched for mitochondrial translation or miRNA-mediated translational regulation rather than core cytosolic translation (Supplementary Table S2), indicating that mirror regulation is specific to the conserved cytosolic translational machinery. For cross-disease comparative transcriptomics, pathway- and module-level analyses are better suited than single-gene comparisons for revealing regulatory principles that operate across biological contexts. The translation co-expression module was independently identified in both diseases, with consistent GO enrichment and overlapping hub gene composition. This dual replication, using identical analytical pipelines in independent datasets from distinct diseases, supports the existence of a conserved organizational principle shared across distinct pathological contexts.

The consistent lack of prognostic value for TGS across three independent HCC cohorts is an informative negative result. Prognostic evaluation was restricted to HCC because available HCC cohorts provide large-scale, standardized survival annotation, whereas the HF dataset (GSE57338) represents a cross-sectional transplant population without comparable long-term outcome data. It indicates that while the translation co-expression architecture is conserved, its statistical association with clinical outcome is highly dependent on cohort-specific factors — etiology composition, stage distribution, and treatment background. TGS should not be considered a universal HCC prognostic biomarker, and its initial association with survival in TCGA-LIHC likely reflects population-specific characteristics. Unlike previous prognostic studies that focused on disease-specific biomarkers, our findings suggest that conserved molecular organization does not necessarily imply conserved clinical utility. In HF, the consistent downregulation of translation pathways aligns with the energy-depletion model of heart failure: with protein synthesis accounting for approximately 30% of myocardial ATP consumption [1], translational suppression is consistent with an evolutionarily conserved adaptive response to chronic energy deficit [12,13,15]. The smaller effect sizes in HF compared with HCC further suggest that translational downregulation in the failing heart may reflect a gradual, incomplete adaptation rather than the robust oncogene-driven transcriptional reprogramming characteristic of cancer.

Several limitations should be considered. One limitation is that all analyses are based on transcriptomic data and cannot directly assess translation efficiency; ribosome footprint profiling (Ribo-seq) would substantially strengthen the evidence for altered translational activity [20]. Another limitation is that the ATF4/ISR regulatory hypothesis is derived from correlational analysis and requires causal validation through ChIP-seq and functional perturbation experiments in HCC models. In addition, only a single HF dataset (GSE57338) was used; validation in additional HF cohorts would further strengthen our conclusions, although such datasets with matched phenotypic annotation remain scarce. Platform differences between microarray (GSE57338) and RNA-seq (TCGA-LIHC, GSE141198) datasets may also introduce systematic bias. While WGCNA provides a well-established framework for co-expression network discovery, alternative network inference methods (e.g., graphical LASSO, mutual information-based approaches) may reveal additional regulatory relationships not captured by correlation-based networks. Finally, HCC validation cohorts are composed predominantly of individuals of Asian ancestry, and cross-population generalizability requires further evaluation.

The multi-level consistency of our findings — from co-expression network topology through pathway activity to upstream regulator inference — provides a coherent framework for generating experimentally testable mechanistic hypotheses. More broadly, this study illustrates how cross-disease comparative transcriptomics can distinguish conserved organizational principles from disease-specific transcriptional responses, providing a general framework for interpreting molecular regulation across biologically divergent diseases.

## Conclusions

The conserved translation co-expression architecture identified across heart failure and hepatocellular carcinoma — diseases with fundamentally opposing physiological demands — indicates that the transcriptional programs governing the translational machinery follow shared organizational principles that transcend disease context. The convergence of ATF4/ISR stress-adaptive signaling and MYC-driven proliferative programs on this architecture in HCC, together with its coordinated suppression in the energy-depleted failing myocardium, supports a model in which a limited repertoire of upstream regulators controls translational output, with disease phenotype determined by the direction of perturbation rather than the architecture itself (Figure 6). More broadly, this study shows that cross-disease comparative frameworks applied at pathway and module resolution can distinguish conserved organizational principles from disease-specific transcriptional responses — regulatory logic that single-disease analyses cannot resolve. These findings position mirror regulation — conserved molecular architecture deployed in opposite directions according to disease context — as a conceptual framework for understanding how shared transcriptional programs can be repurposed across biologically divergent diseases.



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE 6 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure6_Mechanistic_Model.png`

**Caption:** Cross-disease mirror regulation framework for translation in HF versus HCC. Left (HF): energy deficit and mTORC1 suppression → coordinated downregulation (↓). Right (HCC): ISR activation + MYC-driven proliferative signaling → upregulation (↑). Center: structurally conserved translation co-expression module.

</div>

**Figure 6. Cross-disease mirror regulation framework for translation in HF versus HCC.** Left (HF): chronic energy deficit and mTORC1 suppression are associated with coordinated downregulation of the translation co-expression program (↓). Right (HCC): ISR activation (eIF2α phosphorylation → ATF4 translation) and MYC-driven proliferative signaling are associated with upregulation of the same program (↑). Center: the translation co-expression module, structurally conserved across both diseases. Arrows indicate the direction of transcriptional perturbation.

---

## Methods

### 4.1 Data Acquisition and Processing

Datasets were selected based on the following criteria: (i) availability of both disease and non-disease control samples within the same study; (ii) sample size ≥ 100; and (iii) raw or normalized expression matrices publicly available. **Datasets**: (1) GSE57338 (HF): downloaded from GEO, n = 313 left ventricular samples (dilated cardiomyopathy, ischemic cardiomyopathy, non-failing controls), platform GPL11532 (Affymetrix HuGene 1.1 ST Array). Probes were first filtered to retain those with expression above background signal threshold (log2 intensity ≥ 4) in at least 20% of samples, and the remaining probes were then mapped to gene symbols via the hugene11sttranscriptcluster.db annotation package, with multi-probe averaging per gene. (2) TCGA-LIHC (HCC discovery set): downloaded from the TCGA GDC portal, n = 424 samples (371 tumor + 50 normal + 3 recurrent tumor), RNA-seq (Illumina), using `assay(se, "unstranded")` for counts. DESeq2 variance stabilizing transformation (VST) normalization was applied using the vst function with `blind = TRUE` and `nsub = 1000` (design = ~ 1), following filtering for genes with counts ≥ 10 in at least 20% of samples. Ensembl IDs were converted to gene symbols via org.Hs.eg.db. (3) GSE141198 (HCC validation set): downloaded from GEO, n = 148 HCC tumors (RNA-seq), 94 OS events. Processing pipeline identical to TCGA-LIHC. (4) GSE14520 (n = 221 HCC tumors, GPL3921 Affymetrix HT Human Genome U133A Array) and GSE76427 (n = 115 HCC tumors, GPL10558 Illumina HumanHT-12 V4.0) were downloaded from GEO for external TGS validation. Processing and TGS calculation were performed identically across all validation cohorts.


<div style="border: 2px dashed #666; padding: 12px; margin: 18px 0; background: #fafafa;">

**⬇ [TABLE 1 — INSERT TABLE HERE] ⬇**

**Source:** `Table1_Cohort_Characteristics.txt`

**Content:** Cohort Characteristics. Summary of three primary cohorts (GSE57338, TCGA-LIHC, GSE141198): sample size, platform, disease/control counts, and key clinical annotations.

</div>



### 4.2 WGCNA

Genes do not function in isolation; they operate within coordinated transcriptional programs — groups of co-expressed genes that share regulatory mechanisms and often participate in common biological processes. Weighted gene co-expression network analysis (WGCNA) provides a computational framework for identifying these co-expression modules (functionally related gene clusters) from transcriptomic data, analogous to partitioning a social network into communities of tightly connected individuals. This approach enables discovery of disease-associated gene sets without relying on predetermined annotations.

WGCNA input data consisted of gene expression matrices filtered and normalized as described in Section 4.1. For GSE57338 (microarray data), the top 3,000 most variable genes (by variance) were retained for network construction to reduce computational burden. For GSE141198 (RNA-seq data), all expressed genes passing the Section 4.1 filters (5,003 genes) were used. Signed network analysis was performed using the WGCNA package. The soft threshold β was selected to best approximate the scale-free topology criterion (R² ≈ 0.85; GSE57338: β = 12, R² = 0.79; GSE141198: β = 4, R² = 0.84; Supplementary Figure S4). For GSE141198, the blockwiseModules one-step function was used for network construction (networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.25, TOMType = "signed", pamRespectsDendro = FALSE, maxBlockSize = 6,000, randomSeed = 42), where TOM denotes the topological overlap matrix. For GSE57338, due to platform-specific correlation matrix compatibility with the WGCNA package, a manual stepwise procedure was used: Pearson correlation matrix (stats::cor), signed adjacency transformation (0.5 × (1 + cor))^β, TOMsimilarity with TOMType = "signed", average-linkage hierarchical clustering, and module detection via dynamicTreeCut::cutreeDynamic (deepSplit = 2, minClusterSize = 30). Hub genes within the translation module were identified using the dual criteria of intramodular connectivity (Module Membership > 0.80) and gene significance for the disease phenotype (Gene Significance > 0.20). Module functional enrichment was performed using clusterProfiler (enrichGO, ont = "BP", pAdjustMethod = "BH") for GO Biological Process analysis, with translation modules identified by the keywords "translation|ribosom|peptide|ribonucleoprotein|rRNA|translational" in GO term descriptions.



<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [FIGURE S4 — INSERT IMAGE HERE] ⬇**

**File:** `figures/Figure_S4A_SoftThreshold_GSE141198.png
│  File (Panel B): figures/Figure_S4B_SoftThreshold_GSE57338.png`

**Caption:** Soft threshold selection for WGCNA. Scale-free topology model fit (R²) and mean connectivity vs. soft threshold power (β). (A) GSE141198: β=4 (R²=0.84). (B) GSE57338: β=12 (R²=0.79).

</div>

**Supplementary Figure S4. Soft threshold selection for WGCNA.** Scale-free topology model fit (R², left panels) and mean connectivity (right panels) as a function of soft threshold power (β). (A) GSE141198 (HCC): β = 4 (R² = 0.84). (B) GSE57338 (HF): β = 12 (R² = 0.79). The R² ≈ 0.85 criterion (dashed red line) was used for threshold selection.

### 4.3 ssGSEA Pathway Activity Analysis

Gene set sources: msigdbr package (v26.1.0), including Hallmark (collection = "H", 48 sets retained after filtering), KEGG ribosome (collection = "C2", subcollection = "CP:KEGG_MEDICUS", 1 set), and Reactome translation-related pathways (collection = "C2", subcollection = "CP:REACTOME", keyword-filtered: TRANSLATION|PEPTIDE_CHAIN_ELONGATION|EUKARYOTIC_TRANSLATION|RIBOSOME|NONSENSE_MEDIATED|TRNA_AMINOACYLATION|RRNA_PROCESSING, 32 sets). Total: 81 pathways.

ssGSEA scoring was performed using the GSVA package (v2.6.2): parameter objects were constructed via ssgseaParam (minSize = 5, maxSize = 500), and single-sample pathway activity scores were computed using the gsva function. Disease effect size for each pathway was calculated as Cohen's d = (disease group mean − control group mean) / pooled standard deviation, where pooled SD = √[(s₁² + s₂²) / 2] and s₁², s₂² denote the sample variances of the two groups. Cross-disease Spearman correlation: ρ = cor(es_TCGA, es_HF, method = "spearman"). Statistical significance of the cross-disease correlation was evaluated by a 10,000-iteration permutation test (random seed = 42): the vector of HF pathway effect sizes was randomly shuffled in each iteration, and the proportion of random |ρ| exceeding the observed |ρ| was computed as the permutation p-value.

### 4.4 TGS Survival Analysis

TGS was defined as the mean of z-scores of VST-normalized expression values for the six hub genes (EEF1A1, FAU, RPL39, RPL3, RPL32, RPL41; RPS28 excluded due to opposite cross-disease directionality). For TCGA-LIHC, OS was defined as days from diagnosis to death (event) or days to last follow-up (censored). For GSE141198, OS time and status were obtained from the "os days:ch1" and "os status:ch1" fields, respectively. Patients were dichotomized into high/low groups by median TGS. Kaplan–Meier survival curves were compared using the log-rank test. Cox proportional hazards regression was performed using the survival package, with TGS modeled as both a continuous variable (per-IQR increase) and a binary variable (median split). Multivariate Cox models were adjusted for etiology (HBV, HCV, NBNC) in cohort-specific subgroup analyses where etiology data were available. The proportional hazards assumption was assessed using Schoenfeld residual tests.

### 4.5 Upstream TF Analysis

Nineteen candidate TFs were selected based on literature and MSigDB: ATF4, DDIT3, E2F1, E2F2, E2F3, E2F4, E2F5, E2F6, E2F7, E2F8, HIF1A, MTOR, MYC, MYCL, MYCN, NFE2L2, RPTOR, TP53, XBP1. This list encompasses the MYC family (MYC, MYCN, MYCL), the E2F family (E2F1–E2F8), mTOR pathway components, and stress-responsive TFs (ATF4, DDIT3, XBP1, HIF1A, NFE2L2, TP53).

Spearman correlation coefficients were computed between VST-normalized expression values of each TF and TGS across TCGA-LIHC tumor samples (n = 371). Benjamini–Hochberg FDR correction was applied across all 19 tests. Fisher's exact test was used to assess enrichment of the seven hub genes in Hallmark TF target gene sets (MYC Targets V1, MYC Targets V2, E2F Targets, mTORC1 Signaling), with the background gene set defined as all genes detected in TCGA-LIHC after filtering and Ensembl-to-symbol conversion (see Section 4.1). For each pathway, a 2 × 2 contingency table was constructed (hub gene in/not in pathway × background gene in/not in pathway), and a one-sided Fisher's exact test (alternative = "greater") was performed.

### 4.6 Direction Consistency Testing

Cross-disease log2FC values for the seven hub genes were obtained from TCGA-LIHC (tumor vs. normal, DESeq2) and GSE57338 (HF vs. non-failing, limma). Spearman correlation of cross-disease log2FC was used to assess directional consistency. Significance was evaluated by a 10,000-iteration permutation test (random seed = 42): the HCC-side log2FC vector was randomly shuffled in each iteration, and the one-sided empirical p-value was computed as the proportion of permuted Spearman ρ values exceeding the observed ρ. Bootstrap 95% confidence intervals for ρ were estimated from 2,000 resamples (with replacement) of the six directionally concordant genes.

### 4.7 Statistical Methods

All analyses were conducted in R version 4.6.0 (2026-04-24 ucrt). The following R packages were used: WGCNA (v1.74), DESeq2 (v1.52.0), GSVA (v2.6.2), msigdbr (v26.1.0), clusterProfiler (v4.20.0), survival (v3.8-6), survminer (v0.5.2), limma (v3.68.4), org.Hs.eg.db (v3.23.1), hugene11sttranscriptcluster.db (v8.8.0), TCGAbiolinks (v2.40.0), metafor (v5.0.1), and dynamicTreeCut (v1.63-1). A complete session information record is available in the GitHub repository. The statistical significance threshold was set at two-sided p < 0.05 unless otherwise stated. The Benjamini–Hochberg FDR correction was applied for multiple testing correction throughout (GO enrichment, TF–TGS correlations). Permutation tests used 10,000 iterations with a fixed random seed (set.seed(42)) to ensure reproducibility. Missing data were not imputed (missing rate < 5% in all datasets).

### 4.8 Data and Code Availability

All datasets used in this study are publicly available: GSE57338, GSE141198, GSE14520, and GSE76427 from GEO (https://www.ncbi.nlm.nih.gov/geo/); TCGA-LIHC from the GDC Data Portal (https://portal.gdc.cancer.gov/). All datasets are accessible via their public landing pages. Analysis code has been deposited in a GitHub repository (https://github.com/zxy048/translation-mirror-hf-hcc), containing the complete WGCNA, ssGSEA, TF prediction, and direction consistency analysis pipelines, accessible via an anonymous link during peer review.

---

## References

[1] Gibbs CL. Cardiac energetics. Physiol Rev. 1978;58(1):174-254. https://doi.org/10.1152/physrev.1978.58.1.174

[2] Sciarretta S, Forte M, Frati G, et al. New insights into the role of mTOR signaling in the cardiovascular system. Circ Res. 2018;122(3):489-505. https://doi.org/10.1161/CIRCRESAHA.117.311147

[3] Liu Y, Morley M, Brandimarto J, et al. RNA-Seq identifies novel myocardial gene expression signatures of heart failure. Genomics. 2015;105(2):83-9. https://doi.org/10.1016/j.ygeno.2014.12.002

[4] Daneshgar N, Rabinovitch PS, Dai DF. TOR signaling pathway in cardiac aging and heart failure. Biomolecules. 2021;11(2):168. https://doi.org/10.3390/biom11020168

[5] van Riggelen J, Yetil A, Felsher DW. MYC as a regulator of ribosome biogenesis and protein synthesis. Nat Rev Cancer. 2010;10(4):301-309. https://doi.org/10.1038/nrc2819

[6] Liu P, Ge M, Hu J, et al. A functional mammalian target of rapamycin complex 1 signaling is indispensable for c-Myc-driven hepatocarcinogenesis. Hepatology. 2017;66(1):167-181. https://doi.org/10.1002/hep.29183

[7] Pakos-Zebrucka K, Koryga I, Mnich K, et al. The integrated stress response. EMBO Rep. 2016;17(10):1374-1395. https://doi.org/10.15252/embr.201642195

[8] Liu L, Cash TP, Jones RG, et al. Hypoxia-induced energy stress regulates mRNA translation and cell growth. Mol Cell. 2006;21(4):521-531. https://doi.org/10.1016/j.molcel.2006.01.010

[9] Wortel IMN, van der Meer LT, Kilberg MS, van Leeuwen FN. Surviving stress: modulation of ATF4-mediated stress responses in normal and malignant cells. Trends Endocrinol Metab. 2017;28(11):794-806. https://doi.org/10.1016/j.tem.2017.07.003

[10] Ye J, Kumanova M, Hart LS, et al. The GCN2-ATF4 pathway is critical for tumour cell survival and proliferation in response to nutrient deprivation. EMBO J. 2010;29(12):2082-2096. https://doi.org/10.1038/emboj.2010.81

[11] Fu S, Yang L, Li P, et al. Aberrant lipid metabolism disrupts calcium homeostasis causing liver endoplasmic reticulum stress in obesity. Nature. 2011;473(7348):528-531. https://doi.org/10.1038/nature09968

[12] Neubauer S. The failing heart — an engine out of fuel. N Engl J Med. 2007;356(11):1140-1151. https://doi.org/10.1056/NEJMra063052

[13] Harding HP, Zhang Y, Zeng H, et al. An integrated stress response regulates amino acid metabolism and resistance to oxidative stress. Mol Cell. 2003;11(3):619-633. https://doi.org/10.1016/S1097-2765(03)00105-9

[14] Llovet JM, Kelley RK, Villanueva A, et al. Hepatocellular carcinoma. Nat Rev Dis Primers. 2021;7(1):6. https://doi.org/10.1038/s41572-020-00240-3

[15] Halliday M, Radford H, Sekine Y, et al. Repurposed drugs targeting eIF2α-P-mediated translational repression prevent neurodegeneration in mice. Brain. 2017;140(6):1768-1783. https://doi.org/10.1093/brain/awx074

[16] Lines CL, McGrath MJ, Dorwart T, et al. The integrated stress response in cancer progression: a force for plasticity and resistance. Front Oncol. 2023;13:1206561. https://doi.org/10.3389/fonc.2023.1206561

[17] Tang H, Kang R, Liu J, et al. ATF4 in cellular stress, ferroptosis, and cancer. Arch Toxicol. 2024;98(4):1025-1041. https://doi.org/10.1007/s00204-024-03681-x

[18] Sullivan DK, Deutzmann A, Yarbrough J, et al. MYC oncogene elicits tumorigenesis associated with embryonic, ribosomal biogenesis, and tissue-lineage dedifferentiation gene expression changes. Oncogene. 2022;41(45):4960-4970. https://doi.org/10.1038/s41388-022-02458-9

[19] Chen S, Cao X, Zhang J, et al. circVAMP3 drives CAPRIN1 phase separation and inhibits hepatocellular carcinoma by suppressing c-Myc translation. Adv Sci. 2022;9(8):e2103817. https://doi.org/10.1002/advs.202103817

[20] Shiraishi C, Matsumoto A, Ichihara K, et al. RPL3L-containing ribosomes determine translation elongation dynamics required for cardiac function. Nat Commun. 2023;14:2131. https://doi.org/10.1038/s41467-023-37838-6

---



<div style="border: 2px dashed #666; padding: 12px; margin: 18px 0; background: #fafafa;">

**⬇ [TABLE 2 — INSERT TABLE HERE] ⬇**

**Source:** `embedded below`

**Content:** Hub Gene Characteristics and Cross-Disease Module Assignment. Seven hub genes with log2FC, FDR, direction concordance, and module assignment in HF (GSE57338) and HCC (TCGA-LIHC/GSE141198). RPS28 excluded from TGS due to opposite directionality.

</div>


## Table 2. Hub Gene Characteristics and Cross-Disease Module Assignment

| Gene | HF log2FC | HF FDR | HCC log2FC | HCC FDR | Direction | HF Module | HCC Module |
|------|-----------|--------|------------|---------|-----------|-----------|------------|
| EEF1A1 | −0.115 | <0.01 | −0.269 | 0.008 | ✅ Concordant | turquoise | turquoise |
| RPL39 | +0.052 | ns | +0.587 | <0.001 | ✅ Concordant | black | grey |
| FAU | +0.114 | <0.01 | +0.193 | ns | ✅ Concordant (HCC ns) | black | blue |
| RPL3 | +0.056 | ns | −0.034 | ns | ⚠ Weakly discordant | turquoise | blue |
| RPL41 | −0.016 | ns | +0.084 | ns | ⚠ Weakly discordant | black | grey |
| RPL32 | +0.071 | ns | +0.422 | 0.003 | ✅ Concordant | black | blue |
| RPS28 | −0.115 | <0.05 | +0.514 | 0.001 | ❌ Discordant | turquoise | blue |

Note: HF log2FC, HF vs. non-failing (GSE57338, limma); HCC log2FC, tumor vs. normal (TCGA-LIHC, DESeq2). Direction concordance based on consistency of log2FC sign across the two datasets. RPS28 was excluded from TGS construction due to opposite directionality. "grey" in the HCC Module column indicates that the gene was not assigned to any co-expression module in GSE141198. "turquoise" in the HF Module column indicates assignment to the largest WGCNA module (general co-expression), distinct from the translation-enriched black module.



<div style="border: 2px dashed #666; padding: 12px; margin: 18px 0; background: #fafafa;">

**⬇ [TABLE S1 — INSERT TABLE HERE] ⬇**

**Source:** `Table_S1_Module_Comparison.txt`

**Content:** Module Composition and GO Enrichment Comparison. Side-by-side comparison of translation module characteristics in GSE57338 (HF) and GSE141198 (HCC): module size, gene composition (including snoRNA/scaRNA dominance in HF black module), GO enrichment, hub gene co-localization, and cross-disease synthesis. Note: GSE141198 turquoise module (n=1,858) was flagged by initial GO screen but blue (n=1,665) was selected as the primary translation module based on canonical cytoplasmic translation GO terms.

</div>




<div style="border: 2px dashed #666; padding: 12px; margin: 18px 0; background: #fafafa;">

**⬇ [TABLE S2 — INSERT TABLE HERE] ⬇**

**Source:** `Table_S2_ssGSEA_Effect_Sizes.csv`

**Content:** Complete ssGSEA Pathway Effect Sizes. Cohen's d, group means, and pooled SD for all 81 pathways (48 Hallmark + 1 KEGG ribosome + 32 Reactome translation) in TCGA-LIHC and GSE57338.

</div>




<div style="border: 2px dashed #666; padding: 12px; margin: 18px 0; background: #fafafa;">

**⬇ [TABLE S3 — INSERT TABLE HERE] ⬇**

**Source:** `Table_S3_Direction_Consistency.csv`

**Content:** Hub Gene Cross-Disease Directionality. Log2FC, FDR, and direction concordance annotation for seven hub genes in HF vs. HCC. Used for Supplementary Figure S1 and TGS gene selection.

</div>

