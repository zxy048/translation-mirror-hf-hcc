# -*- coding: utf-8 -*-
"""
Convert Manuscript_BBA_v1.md to Hepatology Communications (LWW/AASLD) format:
  - HCC-primary narrative (title, abstract, introduction lead with HCC)
  - Structured abstract: Background / Methods / Results / Conclusion (<= 275 words)
  - Section order: Intro -> Methods -> Results -> Discussion -> Declarations ->
                   Acknowledgements (2 required paragraphs) -> References -> Figure legends
  - AMA-style references, numbered by first citation; all 26 refs renumbered
    because the new introduction changes first-appearance order
  - Figure legends extracted from Results and moved to a closing section
"""
import re

SRC = r"D:\R_projects\revision_analysis\Manuscript_BBA_v1.md"
DST = r"D:\R_projects\revision_analysis\Manuscript_HC_v1.md"

with open(SRC, encoding="utf-8") as f:
    text = f.read()

# ------------------------------------------------------------------
# 1. AMA-style references in NEW first-appearance order.
#    Mapping: new index -> old (BBA) index, so citations in the body
#    can be renumbered old -> new.
# ------------------------------------------------------------------
NEW_REFS = [
    "[1] Lu X, Paliogiannis P, Calvisi DF, Chen X. Role of the mammalian target of rapamycin pathway in liver cancer: from molecular genetics to targeted therapies. Hepatology. 2021;73:49-61.",                                        # old 3
    "[2] Llovet JM, Pinyol R, Kelley RK, et al. Molecular pathogenesis and systemic therapies for hepatocellular carcinoma. Nat Cancer. 2022;3:386-401.",                                                              # old 4
    "[3] Brito Querido J, Diaz-Lopez I, Ramakrishnan V. The molecular basis of translation initiation and its regulation in eukaryotes. Nat Rev Mol Cell Biol. 2024;25:168-186.",                                      # old 7
    "[4] Jiao L, Liu Y, Yu XY, et al. Ribosome biogenesis in disease: new players and therapeutic targets. Signal Transduct Target Ther. 2023;8:15.",                                                                    # old 8
    "[5] Jia X, He X, Huang C, Li J, Dong Z, Liu K. Protein translation: biological processes and therapeutic strategies for human diseases. Signal Transduct Target Ther. 2024;9:44.",                                 # old 9
    "[6] Bartish M, Abraham MJ, Goncalves C, Larsson O, Rolny C, del Rincon SV. The role of eIF4F-driven mRNA translation in regulating the tumour microenvironment. Nat Rev Cancer. 2023;23:408-425.",                  # old 11
    "[7] Ramalho S, Dopler A, Faller WJ. Ribosome specialization in cancer: a spotlight on ribosomal proteins. NAR Cancer. 2024;4:zcae029.",                                                                            # old 22
    "[8] Jakobsen ST, Siersbaek R. Transcriptional regulation by MYC: an emerging new model. Oncogene. 2025;44:1-7.",                                                                                                   # old 12
    "[9] Mericskay M, Zuurbier CJ, Lopaschuk GD, Taegtmeyer H. Cardiac intermediary metabolism in heart failure: substrate use, signalling roles and therapeutic targets. Nat Rev Cardiol. 2025;22:704-727.",              # old 1
    "[10] Abdellatif M, Rainer PP, Sedej S, Kroemer G. Hallmarks of cardiovascular ageing. Nat Rev Cardiol. 2023;20:754-777.",                                                                                        # old 2
    "[11] Xie S, Xu M, Deng W, Tang Q. Metabolic landscape in cardiac aging: insights into molecular biology and therapeutic implications. Signal Transduct Target Ther. 2023;8:114.",                                  # old 6
    "[12] Goul C, Peruzzo R, Zoncu R. The molecular basis of nutrient sensing and signalling by mTORC1 in metabolism regulation and disease. Nat Rev Mol Cell Biol. 2023;24:857-875.",                                   # old 10
    "[13] Lanzer JD, Ramirez Flores RO, Linares Blanco J, et al. A cross-study transcriptional patient map of heart failure defines conserved multicellular coordination in cardiac remodeling. Nat Commun. 2025;16:9659.", # old 5
    "[14] Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15:550.",                                                                      # old 13
    "[15] Zhang B, Horvath S. A general framework for weighted gene co-expression network analysis. Stat Appl Genet Mol Biol. 2005;4:Article 17.",                                                                       # old 14
    "[16] Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559.",                                                                                      # old 15
    "[17] Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. Innovation (Camb). 2021;2:100141.",                                                                   # old 16
    "[18] Hanzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7.",                                                                       # old 17
    "[19] Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462:108-112.",                                                           # old 18
    "[20] Ritterhoff J, Tian R. Metabolic mechanisms in physiological and pathological cardiac hypertrophy: new paradigms and challenges. Nat Rev Cardiol. 2023;20:739-753.",                                            # old 19
    "[21] Costa-Mattioli M, Walter P. The integrated stress response: from mechanism to disease. Science. 2020;368:eaat5314.",                                                                                          # old 20
    "[22] Wek RC, Anthony TG, Staschke KA. Surviving and adapting to stress: translational control and the integrated stress response. Antioxid Redox Signal. 2023;39:351-373.",                                         # old 21
    "[23] Lu HJ, Koju N, Sheng R. Mammalian integrated stress responses in stressed organelles and their functions. Acta Pharmacol Sin. 2024;45:1095-1114.",                                                            # old 23
    "[24] Neill G, Masson GR. A stay of execution: ATF4 regulation and potential outcomes for the integrated stress response. Front Mol Neurosci. 2023;16:1112253.",                                                    # old 24
    "[25] Lines CL, McGrath MJ, Dorwart T, Conn CS. The integrated stress response in cancer progression: a force for plasticity and resistance. Front Oncol. 2023;13:1206561.",                                          # old 25
    "[26] Shiraishi C, Matsumoto A, Ichihara K, et al. RPL3L-containing ribosomes determine translation elongation dynamics required for cardiac function. Nat Commun. 2023;14:2131.",                                   # old 26
]

# old -> new mapping (old index, 1-based, as in BBA refs list)
OLD_TO_NEW = {1: 9, 2: 10, 3: 1, 4: 2, 5: 13, 6: 11, 7: 3, 8: 4, 9: 5, 10: 12,
              11: 6, 12: 8, 13: 14, 14: 15, 15: 16, 16: 17, 17: 18, 18: 19,
              19: 20, 20: 21, 21: 22, 22: 7, 23: 23, 24: 24, 25: 25, 26: 26}

def renumber_citations(s):
    """Rewrite [n], [n,m], [n-m] in-text citations to new numbering, sorted ascending."""
    def repl(m):
        inner = m.group(1)
        nums = []
        for part in inner.split(','):
            part = part.strip()
            if '-' in part:
                a, b = part.split('-')
                nums.extend(range(int(a), int(b) + 1))
            elif part:
                nums.append(int(part))
        mapped = sorted(OLD_TO_NEW[n] for n in nums if n in OLD_TO_NEW)
        if not mapped:
            return m.group(0)
        # compress consecutive runs: [1,2,3] -> [1-3]; keep simple, use comma list
        return '[' + ','.join(str(n) for n in mapped) + ']'
    return re.sub(r'\[([0-9,\-\s]+)\]', repl, s)


# ------------------------------------------------------------------
# 2. Extract body sections from BBA manuscript
# ------------------------------------------------------------------
def extract_section(title):
    m = re.search(rf'^## {title}\n(.*?)(?=\n## |\Z)', text, re.S | re.M)
    return m.group(1).strip() if m else ""

results = extract_section("3. Results")
discussion = extract_section("4. Discussion")
methods = extract_section("2. Materials and methods")
declarations = extract_section("Declarations")

# ------------------------------------------------------------------
# 3. Split figure legends out of Results
# ------------------------------------------------------------------
def split_paragraphs(block):
    # keep "###" headings attached to their following paragraph
    paras = block.split('\n\n')
    return paras

def extract_legends(results_block):
    """Return (text_without_legends, legends_list)."""
    paras = split_paragraphs(results_block)
    legends = []
    kept = []
    for p in paras:
        if re.match(r'^\*\*Fig\. \d+\.\*\*', p.strip()):
            legends.append(p.strip())
        else:
            kept.append(p)
    return '\n\n'.join(kept), legends

results_body, figure_legends = extract_legends(results)

# ------------------------------------------------------------------
# 4. New HCC-first introduction (citations already in new numbering)
# ------------------------------------------------------------------
intro = """Hepatocellular carcinoma (HCC) is an aggressive malignancy whose growth is sustained by unrestrained proliferation, which in turn demands a dramatically expanded capacity for protein synthesis[1,2]. Translation — the final and most energetically expensive step of gene expression — is central to this demand[3]. In HCC, oncogenic pathways converge on the translational apparatus to upregulate ribosome biogenesis and protein synthesis[4,5], and translational control is increasingly recognized as a targetable vulnerability in cancer[6,7]. MYC, a master regulator of both proliferation and translation, is among the most frequently activated oncogenes in HCC[8].

Translation-related transcriptional programs, however, are not uniformly activated across disease states. Heart failure (HF), a chronic degenerative condition of the energy-depleted myocardium, is characterized by suppression of anabolic processes, including ribosome biogenesis and cap-dependent translation, as a response associated with energy conservation[9,10]. mTORC1 signaling, a central anabolic regulator of translation, is inhibited in the failing heart[11,12]. These mirror-image pressures on the same molecular machinery raise a fundamental question: does the same translation-related transcriptional program operate in both diseases, with disease context determining only the direction of perturbation[13]?

Answering this question requires moving beyond single-disease frameworks. Existing transcriptomic studies of translation-related programs have largely been confined to single-disease contexts, identifying differentially expressed genes, co-expression modules, or prognostic signatures within one disease at a time[13]. While these studies have generated valuable disease-specific insights, they cannot distinguish features that are disease-specific from those reflecting general organizational principles of the translational system. Cross-disease comparative transcriptomics — the parallel analysis of independent datasets from distinct diseases using identical computational pipelines — provides a direct strategy. Here, we applied an integrative cross-disease framework to investigate how disease context shapes translation-related transcriptional programs in HCC, with three objectives: (i) identify translation-related co-expression modules in HCC and HF through systematic WGCNA-based module discovery; (ii) quantitatively compare pathway-level perturbation across diseases; and (iii) investigate transcriptional correlates of translation-related activity. Rather than assuming module conservation, we explicitly tested whether network organization is shared or disease-specific, and whether pathway-level perturbation patterns are consistent across independent cohorts."""

# ------------------------------------------------------------------
# 5. Discussion: HCC-anchored opening (rest unchanged, citations renumbered)
# ------------------------------------------------------------------
disc_paras = split_paragraphs(discussion)
# renumber citations in the whole discussion first, then adjust opening para
discussion_rn = renumber_citations(discussion)
disc_paras = split_paragraphs(discussion_rn)
opening = "This study investigated how translation-related transcriptional programs in hepatocellular carcinoma are organized, and whether this organization reflects conserved network architecture or disease-context-dependent remodeling revealed through comparison with heart failure — a disease of fundamentally opposing translational demand."
disc_paras[0] = opening
discussion_final = '\n\n'.join(disc_paras)

# ------------------------------------------------------------------
# 6. Methods and Results: renumber citations
# ------------------------------------------------------------------
methods_rn = renumber_citations(methods)
results_rn = renumber_citations(results_body)

# rename "Materials and methods" heading to "Methods" (already done at section level)

# ------------------------------------------------------------------
# 7. Declarations (adapted from BBA Declarations section)
# ------------------------------------------------------------------
# Reuse BBA declarations but keep HC-friendly subheadings; renumber citations
decl_rn = renumber_citations(declarations)
# Remove the "## Declarations" heading if present in extracted body
decl_rn = decl_rn.replace('## Declarations', '').strip()
# Keep subheadings as-is; the section is added below.

# ------------------------------------------------------------------
# 8. Title, authors, abstract (HCC-first), keywords
# ------------------------------------------------------------------
title = "Translation-related transcriptional programs in hepatocellular carcinoma: disease-context-dependent remodeling revealed by cross-disease comparison with heart failure"

authors = "Yuhe Hong<sup>a</sup>, Yili Liu<sup>a</sup>, Xianqi Li<sup>a</sup>, Xin Wu<sup>b#</sup>, Junhong Wang<sup>a,c*</sup>"

affiliations = """<sup>a</sup> Department of Cardiology, The First Affiliated Hospital of Nanjing Medical University, Nanjing, China

<sup>b</sup> Department of Reproduction Medicine, The First Affiliated Hospital of Nanjing Medical University, Nanjing, China

<sup>c</sup> Department of Cardiology, Liyang People's Hospital, Liyang, China

<sup>*</sup> Lead corresponding author. E-mail: wangjunhong@jsph.org.cn

<sup>#</sup> Co-corresponding author. E-mail: wuxin0220@aliyun.com"""

abstract = """**Background:** Hepatocellular carcinoma (HCC) requires sustained protein synthesis to support malignant proliferation, and translation-related transcriptional programs are consistently activated in HCC. Whether this activation reflects a conserved program shared with diseases of opposing translational demand, or disease-context-dependent remodeling, remains unresolved.

**Methods:** We performed weighted gene co-expression network analysis (WGCNA) of HCC tumor transcriptomes (GSE141198, *n* = 148; TCGA-LIHC, *n* = 421) and failing myocardium (GSE57338, *n* = 313), single-sample gene set enrichment analysis (ssGSEA) across 81 pathways, and transcription factor (TF)-translation-associated transcriptional score (TATS) correlation analysis. Same-organ disease controls (hypertrophic cardiomyopathy, GSE141910; cirrhosis, GSE89377) and an independent HF cohort (GSE116250) were analyzed to assess disease-context dependence.

**Results:** WGCNA identified translation-related modules in HCC (blue module, 1,315 genes, enriched for cytoplasmic translation) and HF (green module, 227 genes, enriched for ribosome biogenesis) with limited gene-level overlap (62/227 genes; Fisher's OR = 1.7, *P* = 0.039). Despite distinct network architecture, 24 of 33 translation-related pathways exhibited mirror perturbation between diseases (Spearman ρ = −0.598, *P* = 0.0003), a pattern absent in same-organ disease controls (hypertrophic cardiomyopathy, ρ = −0.036; cirrhosis, ρ = +0.402). ATF4 was the strongest individual TF correlate of TATS (ρ = +0.500, FDR < 0.0001), while MYC target pathway activity showed the highest overall association (ρ = +0.753, *P* < 0.0001).

**Conclusion:** Translation-related programs in HCC are organized through disease-specific network architecture yet converge on shared functional themes with opposite perturbation directions relative to HF. These findings identify disease-context-dependent remodeling, with proliferative (MYC) and stress-responsive (ATF4/integrated stress response) signaling as dominant associated regulatory dimensions, providing a framework for translation-targeted investigation in HCC."""

keywords = "Translation; Hepatocellular carcinoma; Heart failure; WGCNA; ssGSEA; ATF4; Integrated stress response; MYC"

# ------------------------------------------------------------------
# 9. Acknowledgements (HC requires two paragraphs)
# ------------------------------------------------------------------
acknowledgements = """**Assistance with the study:** We thank all TCGA, GEO, and MSigDB database researchers and patients involved in the datasets used in this study, for their willingness to share relevant data and their contributions to medical progress.

**Presentation of preliminary data:** none."""

# ------------------------------------------------------------------
# 10. Assemble HC manuscript
# ------------------------------------------------------------------
hc_doc = f"""# {title}

{authors}

{affiliations}

**Keywords:** {keywords}

## Abstract

{abstract}

## 1. Introduction

{intro}

## 2. Methods

{methods_rn}

## 3. Results

{results_rn}

## 4. Discussion

{discussion_final}

## Declarations

{decl_rn}

## Acknowledgements

{acknowledgements}

## References

{chr(10).join(NEW_REFS)}

## Figure legends

{chr(10).join(figure_legends)}
"""

with open(DST, 'w', encoding='utf-8') as f:
    f.write(hc_doc)

# ------------------------------------------------------------------
# 11. Report
# ------------------------------------------------------------------
def word_count(s):
    return len(s.split())

body_words = word_count(intro) + word_count(methods_rn) + word_count(results_rn) + word_count(discussion_final)
abstract_words = word_count(abstract)
print(f"Written: {DST}")
print(f"Reference count: {len(NEW_REFS)}")
print(f"Abstract words: {abstract_words} (HC limit: 275)")
print(f"Body words (Intro+Methods+Results+Discussion): {body_words} (HC limit: 5000 excl. refs)")
print(f"Figure legends extracted: {len(figure_legends)}")
print(f"Intro length: {len(intro)} chars")
print(f"Methods length: {len(methods_rn)} chars")
print(f"Results length: {len(results_rn)} chars")
print(f"Discussion length: {len(discussion_final)} chars")
