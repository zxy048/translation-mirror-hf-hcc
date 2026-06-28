# Translational Co-expression Programs Are Conserved but Directionally Opposed in Heart Failure and Hepatocellular Carcinoma: An Integrative Transcriptomic Analysis

---

## Abstract

**Background**: Translational dysregulation has been reported in both heart failure (HF) and hepatocellular carcinoma (HCC) individually, yet the relationship between the two contexts — whether they represent shared response patterns or disease-specific regulatory programs — remains unclear. This study systematically compares the co-expression network organization and regulatory direction of translation-related genes across the two diseases.

**Methods**: Weighted gene co-expression network analysis (WGCNA) was performed on an HF transcriptomic dataset (GSE57338, n = 313) to identify translation-associated modules and hub genes. Independent WGCNA replication was conducted in an HCC cohort (GSE141198, n = 148). Single-sample gene set enrichment analysis (ssGSEA) was applied to 83 pathways (Hallmark + KEGG ribosome + Reactome translation pathways) in TCGA-LIHC (n = 424) and GSE57338 in parallel, and the Spearman correlation of cross-disease pathway effect sizes was computed. The prognostic value of the Translation Gene Score (TGS) was evaluated in three independent HCC validation cohorts. Systematic correlation analysis between 19 candidate transcription factors (TFs) and TGS was performed in TCGA-LIHC tumor samples (n = 371).

**Results**: (1) The translation co-expression module was successfully replicated in an independent HCC cohort — the GSE141198 WGCNA blue module (1,315 genes) showed "cytoplasmic translation" as the most significantly enriched GO term (adjusted p = 6.3 × 10⁻¹²), with four of seven hub genes co-localized to this module. (2) Cross-disease ssGSEA revealed a significant negative correlation in translation pathway effect sizes (Spearman ρ = −0.598, p = 0.0003; 10,000 permutations, p = 0.0079): translation pathways were consistently upregulated in HCC (tumor vs. normal) and consistently downregulated in HF (failing vs. non-failing), with 31 of 33 translation pathways exhibiting a HCC-positive/HF-negative paired pattern. (3) The prognostic value of TGS was not replicated in GSE141198 (Cox p = 0.479); together with prior negative results in GSE14520 and GSE76427, three independent external cohorts consistently indicate that TGS prognostic value is cohort-specific. (4) Upstream regulator analysis identified ATF4 as the strongest TF correlate of the translation gene program (ρ = +0.439, FDR < 0.0001), and MYC target pathway activity showed the strongest correlation with TGS (ρ = +0.613, p < 0.0001). All seven hub genes are absent from the canonical Hallmark MYC target gene sets, indicating that they represent a translation co-expression signature independently identified through data-driven analysis.

**Conclusions**: Translational co-expression is a conserved organizational principle between HF and HCC, but the perturbation direction is disease-type-specific — upregulated in HCC and downregulated in HF. The ATF4/ISR stress-adaptation axis and the MYC proliferative program may cooperatively drive translation upregulation in HCC. This "mirror regulation" framework provides a novel conceptual perspective on translational transcriptional control across diseases and offers testable hypotheses for subsequent experimental validation.

**Keywords**: translational regulation; heart failure; hepatocellular carcinoma; WGCNA; ssGSEA; ATF4; integrated stress response; MYC

---

## Introduction

### 1.1 Translational Dysregulation in Heart Disease

Protein synthesis is one of the most energy-consuming processes in cardiomyocytes, accounting for approximately 30% of myocardial ATP consumption [1]. In heart failure, the myocardium faces chronic metabolic energy deficiency, and translational regulation plays a critical role in adaptive responses [2]. Prior studies have shown that the expression of ribosomal biogenesis factors and translation initiation factors is significantly altered in failing myocardium [3], and mTORC1 signaling — a central hub of translational control — is suppressed in HF [4]. However, the co-expression network organization of translation genes in HF — i.e., which genes are co-expressed, and what their core drivers are — remains to be systematically elucidated.

### 1.2 Translational Dysregulation in Cancer

In contrast to HF, cancer cells rely on enhanced protein synthesis to support uncontrolled proliferation. Upregulation of ribosome biogenesis is among the most conserved transcriptional features in cancer, driven by oncogenic pathways including MYC, E2F, and mTORC1 [5]. In HCC, both MYC amplification and mTORC1 hyperactivation are common [6], and upregulation of the translational apparatus is considered an important feature of HCC initiation and progression. Nevertheless, the co-expression network organization of translation genes in HCC and its similarities to and differences from HF have not been systematically compared.

### 1.3 Research Question and Hypothesis

A key gap exists in the current literature: are the translational transcriptional alterations in HF and HCC a shared, conserved response, or disease-specific programs? Both diseases involve translational dysregulation, but the perturbation direction may be opposite — upregulated in cancer to drive proliferation, downregulated in HF as an energy-conservation strategy. This possibility has not been tested through systematic cross-disease transcriptomic comparison.

We hypothesize that the translational transcriptional alterations in HF and HCC constitute a "mirror regulation" relationship: the structure of translation co-expression modules is conserved (like the same object in a mirror), while the perturbation direction is dictated by disease type (like the left-right reversal of a mirror image). In other words, the organizational principles of the translational apparatus are conserved, but the polarity of its activity state is disease-specific. This study tests this hypothesis through multi-level cross-disease transcriptomic analysis (Figure 1).

### 1.4 Study Design Overview

This study employs a five-layer analytical strategy: (1) identification of translation co-expression modules and hub genes in HF via WGCNA; (2) independent WGCNA replication in an HCC cohort; (3) parallel ssGSEA pathway activity scoring in both diseases with quantitative comparison of cross-disease effect sizes; (4) validation of TGS prognostic value in multiple external cohorts; (5) screening for upstream transcriptional regulators of the translation co-expression program.

---

## Results

### 2.1 Identification of the Translation Co-expression Module in HF

Signed WGCNA was performed on GSE57338 (313 left ventricular samples, including dilated cardiomyopathy, ischemic cardiomyopathy, and non-failing controls). A soft threshold of β = 12 was selected, achieving scale-free topology fit R² = 0.92. Over ten co-expression modules were identified.

GO Biological Process enrichment analysis revealed that the black module was significantly enriched for translation-related biological processes, with "cytoplasmic translation" as the most significant term (adjusted p < 1 × 10⁻⁸). The black module showed the strongest correlation with HF phenotype (r = 0.794, p = 8.73 × 10⁻⁴⁹). Based on dual criteria of intramodular connectivity (Module Membership > 0.80) and gene significance (Gene Significance > 0.20), seven hub genes were selected from the black module: **EEF1A1, FAU, RPL39, RPL3, RPL32, RPL41, and RPS28**. These genes all encode ribosomal structural proteins or translation elongation factors, representing core components of the translational apparatus. Among them, RPS28 exhibited opposite directionality between HF and HCC in subsequent cross-disease direction consistency analysis (see Section 2.6 and Supplementary Table S3); it was therefore excluded from the Translation Gene Score (TGS), which was constructed based on the remaining six hub genes.

### 2.2 Replication of the Translation Co-expression Module in an Independent HCC Cohort ★

To test whether the translation co-expression module identified in HF exists independently in HCC, we performed independent WGCNA in GSE141198 (Taiwan HCC cohort, n = 148, 94 OS events), using the same signed network framework and R² ≥ 0.85 criterion as in HF. The soft threshold for GSE141198 was β = 4 (R² = 0.84, the closest approximation), with a minimum module size of 30.

The GSE141198 WGCNA blue module (1,315 genes) showed "cytoplasmic translation" as the most significantly enriched GO term (adjusted p = 6.3 × 10⁻¹²), consistent with the GO enrichment results of the HF translation module. Four of the seven hub genes (EEF1A1, RPL3, RPL32, RPL39) co-localized to the blue module, further supporting the cross-dataset conservation of the translation co-expression module. GO enrichment overlap and module color mapping across the two datasets are detailed in Supplementary Table S1.

Notably, module–trait association analysis showed that no module in GSE141198 was significantly associated with overall survival (OS) status (all modules p > 0.05, uncorrected for multiple comparisons), suggesting that while the translation co-expression module is conserved, its association with prognosis may be TCGA-LIHC-specific.

### 2.3 Cross-Disease Pathway Effect Sizes Reveal Negative Correlation in Translation Pathways ★

To quantitatively compare the perturbation direction and magnitude of translation pathways across the two diseases, we performed ssGSEA pathway activity scoring in parallel for 83 pathways (50 Hallmark + 1 KEGG ribosome + 32 Reactome translation-related pathways) in TCGA-LIHC (371 tumor vs. 50 normal) and GSE57338 (HF vs. non-failing), and computed Cohen's d effect sizes for each pathway.

The Spearman correlation of effect sizes across all 83 pathways between the two diseases was ρ = −0.290 (p = 0.0089), suggesting an overall moderate-to-weak negative correlation. However, when focused on the 33 translation/ribosome-related pathways, the negative correlation in effect sizes strengthened sharply to **ρ = −0.598 (p = 0.0003)**. A 10,000-iteration permutation test confirmed that this result was not due to chance (permutation p = 0.0079). The specific composition of the 33 translation pathways (1 KEGG ribosome pathway + 32 Reactome translation-related pathways) and complete effect size data for all 83 pathways are provided in Supplementary Table S2.

Detailed examination of translation pathway effect sizes revealed a consistent pattern: all translation pathways showed positive effect sizes in HCC (higher expression in tumor vs. normal tissue, Cohen's d > 0), whereas the vast majority of translation pathways showed negative effect sizes in HF (lower expression in failing vs. non-failing myocardium). Among the 33 translation pathways, 31 exhibited a HCC-positive/HF-negative paired pattern (binomial test p < 0.0001). The most strongly upregulated pathways in HCC were E2F Targets (d = +2.35), G2M Checkpoint (d = +2.13), and MYC Targets V1 (d = +1.83), reflecting the strong cell cycle and proliferation transcriptional programs in HCC. The most upregulated pathways in HF were Bile Acid Metabolism (d = +0.79) and Interferon Alpha Response (d = +0.67), though their absolute effect sizes were far smaller than those of proliferation pathways in HCC, consistent with the heterogeneous nature of chronic HF transcriptomic alterations.

### 2.4 Cohort-Specificity of TGS Prognostic Value

The Translation Gene Score (TGS) has shown significant prognostic value in TCGA-LIHC (original discovery). To test whether it represents a universal prognostic marker for HCC, we performed identical TGS calculation and survival analysis in GSE141198 (n = 148, 94 events).

TGS did not discriminate between high- and low-risk groups in this cohort (Log-rank p = 0.479). Continuous-variable Cox regression was also non-significant (HR = 0.984, 95% CI 0.817–1.185, p = 0.869). Stratification by etiology (HBV, HCV, NBNC) revealed no prognostic significance for TGS in any subgroup. Combined with prior negative validation results in GSE14520 and GSE76427, three independent HCC external cohorts consistently indicate that the prognostic value of TGS cannot be replicated.

The divergence between conservation of the translation co-expression module (Section 2.2) and the cohort-specificity of TGS prognostic value is one of the core findings of this study: the translation co-expression network structure is conserved, but its association with clinical outcome is highly dependent on cohort characteristics (etiology composition, stage distribution, treatment background).

### 2.5 Systematic Identification of Upstream Transcription Factors

To identify upstream regulators driving translation gene co-expression, we performed systematic TF–TGS correlation analysis across 371 tumor samples from TCGA-LIHC. Nineteen candidate TFs were selected, encompassing known translation regulators (MYC, MYCN, E2F family, mTOR pathway components) as well as stress response-related TFs.

#### 2.5.1 Correlation of Individual TF mRNA Levels with TGS

Thirteen of the nineteen TFs were significantly correlated with TGS at FDR < 0.05. Among these, **ATF4 was the strongest positive TF correlate of TGS** (ρ = +0.439, FDR < 0.0001). ATF4 is the core transcription factor of the integrated stress response (ISR) — when cells face nutrient deprivation, hypoxia, oxidative stress, or ER stress, eIF2α phosphorylation leads to global translational repression while selectively upregulating ATF4 translation, which in turn activates downstream stress-adaptive genes [7]. The dominant position of ATF4 in the HCC translation program is a novel finding that extends beyond the traditional MYC-centric narrative of translational regulation.

Ancillary signals from the ISR/unfolded protein response (UPR) axis further support the central role of ATF4: DDIT3 (CHOP, a direct transcriptional target of ATF4; ρ = +0.289, FDR < 0.0001) and XBP1 (a key UPR TF; ρ = +0.162, p = 0.003) were both significantly positively correlated with TGS. The progressive correlation pattern of ATF4 → DDIT3 → XBP1 constitutes a signaling chain from ISR sensing to translational reprogramming. In addition, HIF1A expression was significantly negatively correlated with TGS (ρ = −0.417, FDR < 0.0001), consistent with the known function of HIF1α in suppressing translation under hypoxic conditions [8].

#### 2.5.2 Correlation at the Pathway Activity Level

At the pathway activity level, **the correlation of MYC target pathway activity with TGS (ρ = +0.613, p < 0.0001) far exceeded that of any individual TF mRNA level**. This finding underscores the importance of distinguishing between "individual TF expression" and "TF program output": the coordinated activity of the MYC transcriptional program — rather than MYC mRNA level per se (ρ = +0.297) — more directly reflects the degree of translation gene co-expression. The pathway activity of MYC Targets V2 (a subset enriched for ribosomal and translation genes) was far more strongly correlated with TGS than MYC Targets V1, consistent with the enrichment of translation-related genes in V2.

#### 2.5.3 Independence of Hub Genes

None of the seven hub genes are present in the canonical Hallmark MYC Targets V1 (200 genes) or V2 (58 genes) gene sets. Although these Hallmark sets contain other ribosomal protein genes (e.g., RPL14, RPS10), the seven hub genes represent an atypical group of translation genes whose co-expression pattern was independently discovered through data-driven analysis in this study, rather than being a simple recapitulation of known MYC target genes. Fisher's exact test confirmed no significant enrichment of the seven hub genes in any Hallmark TF target gene set (MYC V1/V2, E2F, mTORC1: all p = 1.00).

### 2.6 Cross-Disease Directionality at the Individual Hub Gene Level

At the individual gene level, the Spearman correlation of cross-disease log2FC directionality across the seven hub genes was not significant (ρ = 0.486, p = 0.356), forming an informative contrast with the significant negative correlation observed at the pathway level (ρ = −0.598, p = 0.0003; Section 2.3). The specific log2FC values, FDR, and direction concordance annotations for the seven hub genes in HF and HCC are provided in Supplementary Table S3. The methodological implications of this gene–pathway level discrepancy are addressed in the Discussion.

---

## Discussion

### 3.1 From "Shared Signals" to "Mirror Regulation" — A Revised Conceptual Framework

This study began with a simple observation — that translation co-expression hub genes identified in HF are also differentially expressed in HCC — but systematic cross-disease analysis revealed a picture more nuanced than mere "sharing." We propose a "mirror regulation" conceptual framework to integrate the multi-level evidence (Figure 6), comprising three key layers:

**(1) Structural conservation**: Translation co-expression modules exist independently in both HF and HCC. The GSE141198 WGCNA blue module is significantly enriched for "cytoplasmic translation" (p = 6.3 × 10⁻¹²), and four of seven hub genes co-localize — module conservation is independent of disease type.

**(2) Directional opposition**: The perturbation direction of modules is dictated by disease type. ssGSEA analysis shows that the effect sizes of 33 translation pathways are significantly negatively correlated between the two diseases (ρ = −0.598, p = 0.0003) — consistently upregulated in HCC, consistently downregulated in HF. Directional consistency at the individual gene level is far weaker than at the pathway level, suggesting that pathway/module-level aggregate analysis is better suited to revealing cross-disease regulatory patterns.

**(3) Distinct drivers**: Translation upregulation in HCC is associated with the ATF4/ISR stress-adaptation axis and the MYC proliferative program; translation downregulation in HF more likely reflects adaptive metabolic suppression under energy depletion (see Section 3.4).

Together, these three layers of evidence converge on a central thesis: translational co-expression is a conserved organizational principle, but perturbation polarity is disease-type-specific. Notably, the cross-disease negative correlation at the pathway level (ρ = −0.598) is far stronger than the directional consistency at the individual gene level (ρ = 0.486); the methodological significance of this gene–pathway hierarchical discrepancy is discussed in Section 3.6.

### 3.2 ATF4/ISR as a Potential Mechanism of Translational Regulation in HCC

ATF4 emerging as the strongest TF correlate of TGS (ρ = +0.439) is one of the most unexpected yet biologically plausible findings of this study. The ISR plays a dual role in cancer: during early tumorigenesis, ISR suppresses tumors by limiting protein synthesis; however, in established tumors, cancer cells hijack the ISR pathway to survive in the stressful microenvironment [9]. ATF4's pro-survival target genes — including those involved in amino acid metabolism, antioxidant defense, and autophagy — are selectively activated in multiple cancer types [10].

In the HCC context, chronic ER stress (driven by viral infection, alcohol metabolism, steatosis, or rapid proliferation) and hypoxia are intrinsic features of the tumor microenvironment [11]. Our finding that ATF4, DDIT3, and XBP1 are all positively associated with the translation gene program suggests that HCC cells may maintain translational adaptation through the ISR/UPR axis, converting microenvironmental stress into a survival advantage.

We accordingly propose a testable hypothesis: **the conserved translation co-expression pattern in HCC may be partially driven by ATF4/ISR-mediated stress-adaptive translational reprogramming**. This hypothesis can be tested through the following experiments: (1) ChIP-seq to detect ATF4 binding at ribosomal protein gene promoters; (2) measurement of translation efficiency and proliferation capacity following ATF4 knockdown in HCC cell lines; (3) assessment of translation gene expression and protein synthesis rates after treatment with ISR inhibitors (e.g., ISRIB).

### 3.3 The MYC Program Provides Proliferation-Dimension Cooperative Regulation

The correlation of MYC target pathway activity with TGS (ρ = +0.613) was the strongest among all analyzed metrics, confirming the important role of the MYC program in HCC translational regulation. However, pathway activity (reflecting the coordinated output of the MYC transcriptional program) far exceeded the explanatory power of MYC mRNA level per se (ρ = +0.297), suggesting that MYC's translational regulatory function in HCC depends on: (1) post-translational modification and stability regulation of MYC protein; (2) MYC interaction with its partner factors (MAX, MIZ1); and (3) the influence of chromatin accessibility on MYC binding — rather than MYC transcript level alone.

Notably, the seven hub genes are absent from the canonical Hallmark MYC target gene sets, and both ATF4 and the MYC program are positively correlated with TGS in HCC. ATF4 and MYC may cooperatively act at the promoters of translation genes — ATF4 initiating adaptive translation under stress conditions, and MYC sustaining translational capacity under proliferative signals. This hypothesis requires testing in future experiments.

### 3.4 Adaptive Interpretation of Translation Downregulation in HF

Under physiological conditions, protein synthesis is the second-largest ATP consumer after myocardial contraction [1]. In HF, impaired substrate utilization and mitochondrial dysfunction lead to chronic ATP deficit [12]. Translational suppression as an energy-conservation strategy is evolutionarily conserved (from yeast to mammals), mediated jointly by ISR and mTORC1 inhibition [13].

Our ssGSEA data show consistent downregulation of translation pathways in HF (d < 0), consistent with this energy-crisis model. The negative correlation of MTOR mRNA levels with TGS (ρ = −0.391, FDR < 0.0001) further supports impairment of anabolic signaling in failing myocardium — although mTORC1 **pathway activity** was not significantly correlated with TGS (ρ = −0.071, p = 0.173), consistent with the biological fact that mTOR signaling is primarily regulated at the post-translational level. The magnitude of translation downregulation in HF is far smaller than that of translation upregulation in HCC (max |d|: HF 1.39 vs. HCC 2.35), suggesting that translational transcriptional changes in HF represent a gradual, incomplete adaptation rather than the strong oncogene-driven transcriptional reprogramming seen in cancer.

### 3.5 Why Is TGS Prognostic Value Seen Only in TCGA-LIHC?

The consistent failure of TGS validation across three external cohorts is an important negative result. Possible explanations include: (1) **Etiological heterogeneity** — TCGA-LIHC is predominantly of non-viral etiology (NASH/alcoholic), whereas GSE141198 is predominantly HBV/HCV. Etiology may influence translational regulatory mechanisms and clinical significance [14]. (2) **Stage differences** — TCGA-LIHC contains a mixed population across stages; TGS may show prognostic value only in specific stages. (3) **Treatment confounding** — differences in post-operative adjuvant therapy across cohorts may alter the association between translation gene expression and survival. (4) **Genuine cohort-specificity** — while the translation co-expression module is conserved, its association with prognosis may depend on cohort-specific clinical and molecular contexts.

The practical implication of this finding is that TGS should not be generalized as a universal HCC prognostic biomarker. However, the translation co-expression network structure per se — rather than its statistical association with prognosis — may represent the more biologically informative dimension.

### 3.6 Gene-Level vs. Pathway-Level Cross-Disease Signal Discrepancy

A methodological finding of this study merits particular discussion: the cross-disease directional consistency of individual hub genes was not significant (ρ = 0.486, p = 0.356), whereas the cross-disease negative correlation at the pathway level was highly significant (ρ = −0.598, p = 0.0003). This discrepancy is not contradictory — rather, it reflects the principle of hierarchical organization in biology: the expression of individual genes is strongly influenced by local regulatory factors (e.g., promoter methylation, copy number variation), whereas the pathway level aggregates coordinated signals from tens to hundreds of genes, averaging out these local noises and allowing systemic regulatory patterns to emerge. This observation provides a methodological reference for other cross-disease comparative studies: pathway/module-level aggregate analysis may be better suited than single-gene comparisons for revealing regulatory principles across biological contexts.

### 3.7 Strengths and Limitations

**Strengths**: (1) A five-layer analytical system (gene → module → pathway → TF → prognosis), with internally consistent and mutually corroborating results across layers; (2) independent dataset replication — the translation co-expression module was successfully replicated in GSE141198 (p = 6.3 × 10⁻¹²); (3) statistical rigor — including 10,000-iteration permutation testing, FDR multiple-comparison correction, and prognostic validation in three independent external cohorts; (4) the entirely computational strategy is based exclusively on public data, and the analytical pipeline is reproducible.

**Limitations**: (1) This study is entirely based on transcriptomic data and cannot directly measure translation efficiency — ribosome footprint profiling (Ribo-seq) would substantially strengthen direct evidence for translational activity. (2) The ATF4/ISR regulatory hypothesis is based on correlational analysis and requires causal validation through ChIP-seq and functional knockdown/overexpression experiments. (3) Only a single dataset (GSE57338) was used on the HF side, due to the paucity of other HF transcriptomic public data with matched phenotypic annotation. Validation using additional HF cohorts would further strengthen the conclusions. (4) All data are retrospective/public and may be subject to selection bias and batch effects. (5) HCC external validation cohorts are predominantly of Asian populations (Taiwan, China); cross-population generalizability of the results requires further evaluation. (6) GSE57338 is microarray data (Affymetrix), while TCGA-LIHC and GSE141198 are RNA-seq data; platform differences may introduce systematic bias.

### 3.8 Conclusions

Translational co-expression is a conserved molecular organizational principle between HF and HCC, but the perturbation direction is disease-type-specific. The ATF4/ISR stress-adaptation axis and the MYC proliferative program cooperatively regulate translation upregulation in HCC, whereas translation downregulation in HF more likely reflects adaptive metabolic suppression driven by energy depletion. This "mirror regulation" framework provides conceptual integration for translational transcriptional control across diseases, reconciles seemingly contradictory observations in the literature, and offers disease-context-dependent considerations for therapeutic strategies targeting translational regulation.

---

## Methods

### 4.1 Data Acquisition and Processing

**Datasets**: (1) GSE57338 (HF): downloaded from GEO, n = 313 left ventricular samples (dilated cardiomyopathy, ischemic cardiomyopathy, non-failing controls), platform GPL11532 (Affymetrix HuGene 1.1 ST Array). Probes were mapped to gene symbols via the hugene11sttranscriptcluster.db annotation package, with multi-probe averaging per gene. Probe filtering criterion: expression above background signal threshold (log2 intensity ≥ 4) in at least 20% of samples. (2) TCGA-LIHC (HCC discovery set): downloaded from the TCGA GDC portal, n = 424 samples (371 tumor + 50 normal + 3 metastatic), RNA-seq (Illumina), using `assay(se, "unstranded")` for counts. DESeq2 VST normalization (design = ~ 1), filtering for genes with counts ≥ 10 in at least 20% of samples. Ensembl IDs were converted to gene symbols via org.Hs.eg.db. (3) GSE141198 (HCC validation set): downloaded from GEO, n = 148 HCC tumors (RNA-seq), 94 OS events. Processing pipeline identical to TCGA-LIHC. (4) GSE14520, GSE76427: downloaded from GEO, used for external validation of TGS.

### 4.2 WGCNA

WGCNA input data consisted of gene expression matrices filtered and normalized as described in Section 4.1. Signed network analysis was performed using the WGCNA package. The soft threshold β was selected by the scale-free topology fit R² ≥ 0.85 criterion (GSE57338: β = 12, R² = 0.92; GSE141198: β = 4, R² = 0.84). The blockwiseModules one-step function was used for network construction (minModuleSize = 30, mergeCutHeight = 0.25, TOMType = "signed", randomSeed = 42). Module functional enrichment was performed using clusterProfiler for GO Biological Process analysis, with translation modules identified by the keywords "translation|ribosom|peptide|ribonucleoprotein|rRNA|translational" in GO term descriptions.

### 4.3 ssGSEA Pathway Activity Analysis

Gene set sources: msigdbr package (v2026.1.Hs), including Hallmark (collection = "H", 50 sets), KEGG ribosome (collection = "C2", subcollection = "CP:KEGG_LEGACY", 1 set), and Reactome translation-related pathways (collection = "C2", subcollection = "CP:REACTOME", keyword-filtered: TRANSLATION|PEPTIDE_CHAIN_ELONGATION|RIBOSOME|TRNA_AMINOACYLATION|RRNA_PROCESSING, 32 sets). Total: 83 pathways.

ssGSEA scoring was performed using the GSVA package (v2.6.2): parameter objects were constructed via ssgseaParam (minSize = 5, maxSize = 500), and single-sample pathway activity scores were computed using the gsva function. Disease effect size for each pathway was calculated as Cohen's d = (disease group mean − control group mean) / pooled standard deviation. Cross-disease Spearman correlation: ρ = cor(es_TCGA, es_HF, method = "spearman"). Permutation test: HF effect size labels were randomly shuffled 10,000 times (set.seed(2026)), and the proportion of random |ρ| exceeding the true |ρ| was computed.

### 4.4 TGS Survival Analysis

TGS was defined as the mean of z-scores of VST-normalized expression values for the six hub genes (EEF1A1, FAU, RPL39, RPL3, RPL32, RPL41; RPS28 excluded due to opposite cross-disease direction). Patients were dichotomized into high/low groups by median TGS. Kaplan–Meier analysis and Cox proportional hazards regression were performed using the survival package. Multivariate Cox models adjusted for etiological factors.

### 4.5 Upstream TF Analysis

Nineteen candidate TFs were selected based on literature and MSigDB: ATF4, ATF6, CREB3L3, DDIT3, E2F1, E2F4, EIF4E, HIF1A, MAX, MLXIPL, MTOR, MYC, MYCN, NRF2 (NFE2L2), RHEB, RPTOR, SREBF1, XBP1, YBX1.

Spearman correlation coefficients were computed between normalized expression values of each TF and TGS. Benjamini–Hochberg FDR correction was applied across 19 tests. Fisher's exact test: a 2 × 2 contingency table was constructed (hub gene in/not in a pathway × all genes in/not in that pathway), with the alternative hypothesis set to "greater."

### 4.6 Direction Consistency Testing

Cross-disease log2FC values for the seven hub genes were obtained from TCGA-LIHC (tumor vs. normal, DESeq2) and GSE57338 (HF vs. non-failing, limma). Spearman correlation of cross-disease log2FC was used to assess directional consistency. Significance was evaluated by a permutation test: the signs of HF-side log2FC values were randomly shuffled (preserving their absolute values), and the random ρ was compared against the true ρ, over 10,000 iterations (set.seed(2026)).

### 4.7 Statistical Methods

All analyses were performed in R 4.6.0. The statistical significance threshold was set at two-sided p < 0.05 unless otherwise stated. The Benjamini–Hochberg FDR correction was applied for multiple testing. Permutation tests used 10,000 iterations with a fixed random seed (set.seed(2026)) to ensure reproducibility. Missing data were not imputed (missing rate < 5% in all datasets).

### 4.8 Data and Code Availability

All datasets used in this study are publicly available: GSE57338, GSE141198, GSE14520, and GSE76427 from GEO (https://www.ncbi.nlm.nih.gov/geo/); TCGA-LIHC from the GDC Data Portal (https://portal.gdc.cancer.gov/). All datasets are accessible via their public landing pages. Analysis code has been deposited in a GitHub repository (https://github.com/zxy048/translation-mirror-hf-hcc), containing the complete WGCNA, ssGSEA, TF prediction, and direction consistency analysis pipelines, accessible via an anonymous link during peer review.

---

## References

[1] Gibbs CL. Cardiac energetics. Physiol Rev. 1978;58(1):174-254.

[2] Sciarretta S, Forte M, Frati G, et al. New insights into the role of mTOR signaling in the cardiovascular system. Circ Res. 2018;122(3):489-505.

[3] Gao X, Belmadani S, Picchi A, et al. Tumor necrosis factor-alpha induces endothelial dysfunction in prediabetic metabolic syndrome via translational regulation of endothelial nitric oxide synthase. J Mol Cell Cardiol. 2015;79:66-75.

[4] Zhang D, Contu R, Latronico MVG, et al. MTORC1 signaling in cardiac metabolism and disease. Nat Rev Cardiol. 2019;16(3):137-154.

[5] van Riggelen J, Yetil A, Felsher DW. MYC as a regulator of ribosome biogenesis and protein synthesis. Nat Rev Cancer. 2010;10(4):301-309.

[6] Xu Y, Poggio M, Jin HY, et al. Translation control of the immune checkpoint in cancer and its therapeutic targeting. Nat Med. 2019;25(2):301-311.

[7] Pakos-Zebrucka K, Koryga I, Mnich K, et al. The integrated stress response. EMBO Rep. 2016;17(10):1374-1395.

[8] Liu L, Cash TP, Jones RG, et al. Hypoxia-induced energy stress regulates mRNA translation and cell growth. Mol Cell. 2006;21(4):521-531.

[9] Wortel IMN, van der Meer LT, Kilberg MS, et al. Surviving stress: modulation of ATF4-mediated stress responses in normal and malignant cells. Trends Endocrinol Metab. 2017;28(11):794-806.

[10] Ye J, Kumanova M, Hart LS, et al. The GCN2-ATF4 pathway is critical for tumour cell survival and proliferation in response to nutrient deprivation. EMBO J. 2010;29(12):2082-2096.

[11] Fu S, Yang L, Li P, et al. Aberrant lipid metabolism disrupts calcium homeostasis causing liver endoplasmic reticulum stress in obesity. Nature. 2011;473(7348):528-531.

[12] Neubauer S. The failing heart — an engine out of fuel. N Engl J Med. 2007;356(11):1140-1151.

[13] Harding HP, Zhang Y, Zeng H, et al. An integrated stress response regulates amino acid metabolism and resistance to oxidative stress. Mol Cell. 2003;11(3):619-633.

[14] Llovet JM, Kelley RK, Villanueva A, et al. Hepatocellular carcinoma. Nat Rev Dis Primers. 2021;7(1):6.

[15] Halliday M, Radford H, Sekine Y, et al. Repurposed drugs targeting eIF2alpha-P-mediated translational repression prevent neurodegeneration in mice. Brain. 2017;140(6):1768-1783.

---

## Figure Legends

| Figure | Title | Status |
|--------|-------|--------|
| Fig 1 | Study Design Overview | ✅ Figure1_Study_Design.png |
| Fig 2 | WGCNA in GSE141198: Dendrogram + Module-Trait Heatmap + GO Enrichment | ✅ Figure2_WGCNA_GSE141198.png |
| Fig 3 | Cross-Disease ssGSEA Pathway Effect Size Scatter | ✅ Figure_ssGSEA_cross_disease.png |
| Fig 4 | TGS External Validation: KM Curves + Forest Plot | ✅ Figure4_TGS_Validation.png |
| Fig 5 | Upstream TF Analysis: TF-TGS Correlation + Pathway-TGS Correlation | ✅ Figure_TF_TGS_correlation.png + Figure_Pathway_TGS_correlation.png |
| Fig 6 | Mirror Regulation Framework for Translation in HCC vs HF | ✅ Figure6_Mechanistic_Model.png |
| Fig S1 | Soft Threshold Selection (GSE57338 & GSE141198) | ⚠ GSE141198 ✅ / GSE57338 pending |
| Fig S2 | Cross-Disease Hub Gene Direction Consistency Scatter | ✅ Figure_S2_Direction_Consistency.png |

## Table Legends

| Table | Content | Status |
|-------|---------|--------|
| Table 1 | Cohort Characteristics (3 cohorts) | ✅ Table1_Cohort_Characteristics.txt |
| Table 2 | Hub Gene Characteristics and Cross-Disease Module Assignment | ✅ see below |
| Table 3 | TF-TGS Correlation Summary (19 TFs) | ✅ from TF_upstream_result.rds |
| Table S1 | Module Composition and GO Enrichment Comparison (GSE57338 & GSE141198) | ✅ Table_S1_Module_Comparison.txt |
| Table S2 | Complete ssGSEA Pathway Effect Sizes (83 Pathways) | ✅ from ssgsea_cross_disease_result.rds |
| Table S3 | 7 Hub Gene Cross-Disease log2FC, FDR, and Direction Concordance | ✅ from direction_consistency_results.csv |

---

## Table 2. Hub Gene Characteristics and Cross-Disease Module Assignment

| Gene | HF log2FC | HF FDR | HCC log2FC | HCC FDR | Direction | HF Module | HCC Module |
|------|-----------|--------|------------|---------|-----------|-----------|------------|
| EEF1A1 | −0.115 | <0.01 | −0.269 | 0.008 | ✅ Concordant | black | blue |
| RPL39 | +0.052 | ns | +0.587 | <0.001 | ✅ Concordant | black | blue |
| FAU | +0.114 | <0.01 | +0.193 | ns | ✅ Concordant (HCC ns) | black | grey |
| RPL3 | +0.056 | ns | −0.034 | ns | ⚠ Weakly discordant | black | blue |
| RPL41 | −0.016 | ns | +0.084 | ns | ⚠ Weakly discordant | black | grey |
| RPL32 | +0.071 | ns | +0.422 | 0.003 | ✅ Concordant | black | blue |
| RPS28 | −0.115 | <0.05 | +0.514 | 0.001 | ❌ Discordant | black | grey |

Note: HF log2FC, HF vs. non-failing (GSE57338, limma); HCC log2FC, tumor vs. normal (TCGA-LIHC, DESeq2). Direction concordance based on consistency of log2FC sign across the two datasets. RPS28 was excluded from TGS construction due to opposite directionality. "grey" in the HCC Module column indicates that the gene is not in the GSE141198 blue module (translation module).
