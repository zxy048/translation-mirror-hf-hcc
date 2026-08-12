# -*- coding: utf-8 -*-
"""
Convert Manuscript_JTM_v1.md to BBA-MBD format:
  - Title page (title, authors, affiliations, 6 keywords)
  - Summary (unstructured, 100-200 words)
  - Structure: Introduction -> Materials and methods -> Results -> Discussion
  - Numbered references in BBA style (non-superscript, initials+surname)
  - Supplementary file references
"""
import re

SRC = r"D:\R_projects\revision_analysis\Manuscript_JTM_v1.md"
DST = r"D:\R_projects\revision_analysis\Manuscript_BBA_v1.md"

with open(SRC, encoding="utf-8") as f:
    text = f.read()

# ---------- 1. References: BBA style (complete author lists, verified) ----------
BBA_REFS = [
    "[1] M. Mericskay, C.J. Zuurbier, G.D. Lopaschuk, H. Taegtmeyer, Cardiac intermediary metabolism in heart failure: substrate use, signalling roles and therapeutic targets, Nat. Rev. Cardiol. 22 (2025) 704-727.",
    "[2] M. Abdellatif, P.P. Rainer, S. Sedej, G. Kroemer, Hallmarks of cardiovascular ageing, Nat. Rev. Cardiol. 20 (2023) 754-777.",
    "[3] X. Lu, P. Paliogiannis, D.F. Calvisi, X. Chen, Role of the mammalian target of rapamycin pathway in liver cancer: from molecular genetics to targeted therapies, Hepatology 73 (2021) 49-61.",
    "[4] J.M. Llovet, R. Pinyol, R.K. Kelley, A. El-Khoueiry, H.L. Reeves, X.W. Wang, G.J. Gores, A. Villanueva, Molecular pathogenesis and systemic therapies for hepatocellular carcinoma, Nat. Cancer 3 (2022) 386-401.",
    "[5] J.D. Lanzer, R.O. Ramirez Flores, J. Linares Blanco, M. Steier, A.Y. Rangrez, N. Frey, J. Saez-Rodriguez, A cross-study transcriptional patient map of heart failure defines conserved multicellular coordination in cardiac remodeling, Nat. Commun. 16 (2025) 9659.",
    "[6] S. Xie, M. Xu, W. Deng, Q. Tang, Metabolic landscape in cardiac aging: insights into molecular biology and therapeutic implications, Signal Transduct. Target. Ther. 8 (2023) 114.",
    "[7] J. Brito Querido, I. Diaz-Lopez, V. Ramakrishnan, The molecular basis of translation initiation and its regulation in eukaryotes, Nat. Rev. Mol. Cell Biol. 25 (2024) 168-186.",
    "[8] L. Jiao, Y. Liu, X.-Y. Yu, X. Pan, Y. Zhang, J. Tu, Y.-H. Song, Y. Li, Ribosome biogenesis in disease: new players and therapeutic targets, Signal Transduct. Target. Ther. 8 (2023) 15.",
    "[9] X. Jia, X. He, C. Huang, J. Li, Z. Dong, K. Liu, Protein translation: biological processes and therapeutic strategies for human diseases, Signal Transduct. Target. Ther. 9 (2024) 44.",
    "[10] C. Goul, R. Peruzzo, R. Zoncu, The molecular basis of nutrient sensing and signalling by mTORC1 in metabolism regulation and disease, Nat. Rev. Mol. Cell Biol. 24 (2023) 857-875.",
    "[11] M. Bartish, M.J. Abraham, C. Goncalves, O. Larsson, C. Rolny, S.V. del Rincon, The role of eIF4F-driven mRNA translation in regulating the tumour microenvironment, Nat. Rev. Cancer 23 (2023) 408-425.",
    "[12] S.T. Jakobsen, R. Siersbaek, Transcriptional regulation by MYC: an emerging new model, Oncogene 44 (2025) 1-7.",
    "[13] M.I. Love, W. Huber, S. Anders, Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2, Genome Biol. 15 (2014) 550.",
    "[14] B. Zhang, S. Horvath, A general framework for weighted gene co-expression network analysis, Stat. Appl. Genet. Mol. Biol. 4 (2005) Article 17.",
    "[15] P. Langfelder, S. Horvath, WGCNA: an R package for weighted correlation network analysis, BMC Bioinformatics 9 (2008) 559.",
    "[16] T. Wu, E. Hu, S. Xu, M. Chen, P. Guo, Z. Dai, T. Feng, L. Zhou, W. Tang, L. Zhan, X. Fu, S. Liu, X. Bo, G. Yu, clusterProfiler 4.0: a universal enrichment tool for interpreting omics data, Innovation 2 (2021) 100141.",
    "[17] S. Hanzelmann, R. Castelo, J. Guinney, GSVA: gene set variation analysis for microarray and RNA-seq data, BMC Bioinformatics 14 (2013) 7.",
    "[18] D.A. Barbie, P. Tamayo, J.S. Boehm, S.Y. Kim, S.E. Moody, I.F. Dunn, A.C. Schinzel, P. Sandy, E. Meylan, C. Scholl, et al., Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1, Nature 462 (2009) 108-112.",
    "[19] J. Ritterhoff, R. Tian, Metabolic mechanisms in physiological and pathological cardiac hypertrophy: new paradigms and challenges, Nat. Rev. Cardiol. 20 (2023) 739-753.",
    "[20] M. Costa-Mattioli, P. Walter, The integrated stress response: from mechanism to disease, Science 368 (2020) eaat5314.",
    "[21] R.C. Wek, T.G. Anthony, K.A. Staschke, Surviving and adapting to stress: translational control and the integrated stress response, Antioxid. Redox Signal. 39 (2023) 351-373.",
    "[22] S. Ramalho, A. Dopler, W.J. Faller, Ribosome specialization in cancer: a spotlight on ribosomal proteins, NAR Cancer 4 (2024) zcae029.",
    "[23] H.J. Lu, N. Koju, R. Sheng, Mammalian integrated stress responses in stressed organelles and their functions, Acta Pharmacol. Sin. 45 (2024) 1095-1114.",
    "[24] G. Neill, G.R. Masson, A stay of execution: ATF4 regulation and potential outcomes for the integrated stress response, Front. Mol. Neurosci. 16 (2023) 1112253.",
    "[25] C.L. Lines, M.J. McGrath, T. Dorwart, C.S. Conn, The integrated stress response in cancer progression: a force for plasticity and resistance, Front. Oncol. 13 (2023) 1206561.",
    "[26] C. Shiraishi, A. Matsumoto, K. Ichihara, T. Yamamoto, T. Yokoyama, T. Mizoo, A. Hatano, M. Matsumoto, Y. Tanaka, E. Matsuura-Suzuki, S. Iwasaki, S. Matsushima, H. Tsutsui, K.I. Nakayama, RPL3L-containing ribosomes determine translation elongation dynamics required for cardiac function, Nat. Commun. 14 (2023) 2131.",
]

bba_ref_text = '\n\n'.join(BBA_REFS)

# ---------- 2. Extract body sections ----------
def extract_section(title):
    m = re.search(rf'^## {title}\n(.*?)(?=\n## |\Z)', text, re.S | re.M)
    return m.group(1).strip() if m else ""

results = extract_section("Results")
discussion = extract_section("Discussion")
conclusion = extract_section("Conclusion")
methods = extract_section("Methods")
declarations = extract_section("Declarations")

# Merge Conclusion into Discussion (BBA has no separate conclusion)
discussion = discussion + "\n\n" + conclusion

# ---------- 3. Extract Introduction (after Keywords block, before '## Results') ----------
intro_match = re.search(r'\*\*Keywords:\*\*.*?\n---\s*\n\n(.*?)(?=\n## Results)', text, re.S)
intro = intro_match.group(1).strip()

# ---------- 4. Fix figure legend format "**Fig. 1 |" -> "**Fig. 1.**" ----------
results_bba = re.sub(r'\*\*Fig\. (\d+) \|', r'**Fig. \1.**', results)

# ---------- 5. Fix Additional file references -> Supplementary ----------
results_bba = results_bba.replace('Additional file 3: Figure S3b', 'Supplementary Fig. S3b')
results_bba = results_bba.replace('Additional file 9: Table S5', 'Supplementary Table S5')
results_bba = results_bba.replace('Additional file 1: Figure S1', 'Supplementary Fig. S1')
results_bba = results_bba.replace('Additional file 5: Table S1', 'Supplementary Table S1')
results_bba = results_bba.replace('Additional file 4: Figure S4', 'Supplementary Fig. S4')
results_bba = results_bba.replace('Additional file 6: Table S2', 'Supplementary Table S2')
results_bba = results_bba.replace('Additional file 8: Table S4', 'Supplementary Table S4')
results_bba = results_bba.replace('Additional file 11: Table S7', 'Supplementary Table S7')
results_bba = results_bba.replace('Additional file 10: Table S6', 'Supplementary Table S6')
results_bba = results_bba.replace('Additional file 3: Figure S3', 'Supplementary Fig. S3')

methods_bba = methods
for n in range(1, 12):
    methods_bba = methods_bba.replace(f'Additional file {n}: Figure S', f'Supplementary Fig. S')
    methods_bba = methods_bba.replace(f'Additional file {n}: Table S', f'Supplementary Table S')
    methods_bba = methods_bba.replace(f'Additional file {n}', f'Supplementary file {n}')

# ---------- 6. Build BBA document ----------
summary = """Heart failure (HF) and hepatocellular carcinoma (HCC) exhibit dysregulation of translation-related gene expression, yet whether this reflects shared transcriptional organization or disease-specific remodeling remains unresolved. Using weighted gene co-expression network analysis (WGCNA) of HF myocardium (GSE57338, n = 313) and HCC tumors (GSE141198, n = 148), single-sample gene set enrichment analysis (ssGSEA) across 81 pathways, and transcription factor (TF)-translation-associated transcriptional score (TATS) correlation analysis, we identified translation-related modules in both diseases (HF green module, 227 genes; HCC blue module, 1,315 genes) with limited gene-level overlap (62/227 genes; Fisher's OR = 1.7, P = 0.039). Despite distinct network architecture, 24 of 33 translation-related pathways exhibited mirror perturbation between diseases (Spearman rho = -0.598, P = 0.0003), a pattern absent in same-organ disease controls (hypertrophic cardiomyopathy, rho = -0.036; cirrhosis, rho = +0.402). ATF4 was the strongest individual TF correlate with TATS (rho = +0.500, FDR < 0.0001), while MYC target pathway activity showed the highest overall association (rho = +0.753, P < 0.0001). These findings identify disease-context-dependent remodeling, in which shared functional themes emerge through distinct network configurations with opposite perturbation directions, as an organizing principle of translation-related transcriptional programs across diseases."""

# Drop old Conclusion heading content into Discussion? BBA has no separate conclusion.
# Keep Discussion as-is; ensure Conclusion content merged if needed.

bba_doc = f"""# Distinct networks drive mirror perturbation of translation-related transcription in heart failure and hepatocellular carcinoma

Yuhe Hong<sup>a</sup>, Yili Liu<sup>a</sup>, Xianqi Li<sup>a</sup>, Xin Wu<sup>b#</sup>, Junhong Wang<sup>a,c*</sup>

<sup>a</sup> Department of Cardiology, The First Affiliated Hospital of Nanjing Medical University, Nanjing, China

<sup>b</sup> Department of Reproduction Medicine, The First Affiliated Hospital of Nanjing Medical University, Nanjing, China

<sup>c</sup> Department of Cardiology, Liyang People's Hospital, Liyang, China

<sup>*</sup> Lead corresponding author. E-mail: wangjunhong@jsph.org.cn

<sup>#</sup> Co-corresponding author. E-mail: wuxin0220@aliyun.com

**Keywords:** Translation; Heart failure; Hepatocellular carcinoma; WGCNA; ssGSEA; ATF4

## Summary

{summary}

## 1. Introduction

{intro}

## 2. Materials and methods

{methods_bba}

## 3. Results

{results_bba}

## 4. Discussion

{discussion}

## References

{bba_ref_text}

## Declarations

{declarations}

## Supplementary data

Supplementary data to this article can be found online at [URL to be inserted].
Supplementary Table S1. HF green module gene list (227 genes) with module membership and gene significance statistics.
Supplementary Table S2. Complete ssGSEA pathway effect sizes (Cohen's d) for all 81 pathways in HF and HCC.
Supplementary Table S3. Cross-disease log2FC and FDR for translation-related module genes.
Supplementary Table S4. Translation-related pathways showing mirror perturbation (HF down, HCC up) with binomial test statistics.
Supplementary Table S5. GO Biological Process enrichment results for the HF green module (clusterProfiler, BH-adjusted).
Supplementary Table S6. Complete TF-TATS Spearman correlation statistics for all 19 candidate transcription factors.
Supplementary Table S7. Gene overlap analysis between TATS gene set and Hallmark pathway gene sets, with de-overlap sensitivity analysis.
Supplementary Fig. S1. Green module hub gene identification (MM vs. GS scatter plot).
Supplementary Fig. S2. TATS exploratory survival analysis by etiology subgroup.
Supplementary Fig. S3. Soft threshold selection for WGCNA. (a) GSE141198 (HCC, beta = 4). (b) GSE57338 (HF, beta = 17, signed R2 = 0.865).
Supplementary Fig. S4. Green module internal robustness and cross-disease coherence assessment.
"""

with open(DST, 'w', encoding='utf-8') as f:
    f.write(bba_doc)

print(f"Written: {DST}")
print(f"Reference count: {len(BBA_REFS)}")
print(f"Intro length: {len(intro)} chars")
print(f"Methods length: {len(methods_bba)} chars")
print(f"Results length: {len(results_bba)} chars")
print(f"Discussion length: {len(discussion)} chars")
