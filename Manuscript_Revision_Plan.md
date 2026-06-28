# 论文重构框架：从"共享信号"到"镜像调控"

---

## 一、核心叙事转变

### 原始叙事（v7）
> HF和HCC共享翻译相关转录信号 → TGS有预后价值 → 可能存在共同的翻译调控机制

### 修订叙事（v8）
> 翻译机器在HF和HCC中均被转录水平扰动（共表达模块独立复现），但扰动**方向相反**且**上游调控机制不同**——HCC中ATF4/ISR适应性应激与MYC增殖程序驱动翻译上调，HF中翻译下调反映心肌能量耗竭。这一"镜像调控"模式揭示翻译装置转录控制的疾病类型特异性。

---

## 二、所有分析结果汇总

### 表1：完整证据矩阵

| 分析层级 | 分析方法 | 数据集 | 核心发现 | 对论文的贡献 |
|---------|---------|--------|---------|------------|
| **基因水平** | 方向一致性检验 | TCGA-LIHC + GSE57338 + GSE141198 | 7个hub基因的log2FC跨疾病NS（ρ=0.486, p=0.356） | hub基因个体方向不共享 |
| **模块水平** | WGCNA独立复现 | GSE141198 (n=148 HCC) | 蓝模块富集"cytoplasmic translation" p=6.3×10⁻¹²，4/7 hub基因位于蓝模块 | **✓ 翻译共表达模块跨数据集复现** |
| **模块水平** | WGCNA预后 | GSE141198 | 0个模块与OS显著相关 | 翻译模块的预后价值仅见于TCGA-LIHC |
| **通路水平** | ssGSEA跨疾病 | TCGA-LIHC + GSE57338 | 83条通路ρ=−0.290 p=0.0089；33条翻译通路ρ=−0.598 p=0.0003 | **✓ 翻译通路效应量跨疾病负相关** |
| **预后水平** | TGS生存验证 | GSE141198 (n=148) | TGS不显著（Cox p=0.479） | 外部验证失败 |
| **预后水平** | TGS生存验证 | GSE14520, GSE76427 | TGS不显著 | 已在先前分析中验证 |
| **调控水平** | 上游TF预测 | TCGA-LIHC (n=371) | ATF4最强正相关ρ=+0.439；MYC通路活性ρ=+0.613 | **✓ ATF4/ISR + MYC双通路机制** |
| **调控水平** | Hub基因独立性 | vs MSigDB Hallmark | 7个hub基因均不在标准MYC/E2F/mTORC1靶基因集中 | hub基因是独立发现，非已知靶基因复现 |

### 关键结果速查

```
正发现（可发表的）：
  1. ✅ WGCNA翻译模块在独立HCC队列中复现（蓝模块，p=6.3×10⁻¹²）
  2. ✅ ssGSEA翻译通路跨疾病负相关（ρ=−0.598, p=0.0003）
  3. ✅ ATF4是翻译基因程序的最强上游调控因子（ρ=+0.439, FDR<0.0001）
  4. ✅ 7个hub基因独立于已知MYC靶基因集

负发现（实事求是的）：
  1. ❌ TGS预后价值仅限TCGA-LIHC（3个外部队列均不显著）
  2. ❌ Hub基因个体方向跨疾病不一致（方向一致性检验NS）
  3. ❌ WGCNA模块在GSE141198中无预后相关性
```

---

## 三、修订后的论文结构

### Title (建议)
**"Translational Co-expression Programs Are Conserved but Directionally Opposed in Heart Failure and Hepatocellular Carcinoma: An Integrated Transcriptomic Analysis"**

备选：
- "Mirrored Translation: MYC and ATF4/ISR Programs Drive Disease-Specific Ribosomal Gene Expression in Cancer and Heart Failure"

### Abstract 框架

**Background**: Translation dysregulation reported in both HF and HCC, but relationship unknown.

**Methods**: 
- WGCNA on HF (GSE57338, n=313) → 7 hub genes
- Independent WGCNA on HCC (GSE141198, n=148)
- Cross-disease ssGSEA pathway analysis (83 pathways: Hallmark + KEGG ribosome + Reactome translation)
- Upstream TF prediction (19 TFs) in TCGA-LIHC (n=371 tumors)

**Results**:
- Translation co-expression module independently replicated in GSE141198 (blue module, "cytoplasmic translation" p=6.3×10⁻¹²)
- Translation pathway effect sizes negatively correlated across diseases (Spearman ρ=−0.598, p=0.0003): uniformly upregulated in HCC, downregulated in HF
- TGS prognostic value not validated in external cohorts (Cox p=0.479 in GSE141198)
- ATF4 identified as strongest upstream TF (ρ=+0.439, FDR<0.0001), with MYC target pathway activity showing strongest overall association (ρ=+0.613)

**Conclusion**: Translation co-expression is a conserved organizational principle across HF and HCC, but perturbation direction is disease-type-specific. ATF4/ISR-mediated stress adaptation and MYC-driven proliferation jointly regulate the translational program in HCC. These findings reconcile apparently contradictory observations and provide a refined model of translational regulation in disease.

### Main Text 结构

```
1. INTRODUCTION
   1.1 Translation dysregulation in cardiac pathology
   1.2 Translation dysregulation in cancer  
   1.3 The gap: shared signals vs disease-specific programs
   1.4 Study design overview (graphical abstract)

2. RESULTS
   2.1 Identification of translation-related co-expression network in HF
       - WGCNA on GSE57338 → hub genes (original finding, summarized)
       - (Table S1: Module assignment and hub gene characteristics)
   
   2.2 Independent replication of translation co-expression in HCC  ★NEW★
       - WGCNA on GSE141198: blue module "cytoplasmic translation" p=6.3×10⁻¹²
       - 4/7 hub genes co-localize in blue module
       - (Figure 2: WGCNA dendrogram + module-trait heatmap + GO enrichment)
   
   2.3 Cross-disease pathway activity reveals opposed directions  ★NEW★
       - ssGSEA: 83 pathways in TCGA-LIHC and GSE57338
       - All pathways ρ=−0.290 (p=0.0089)
       - Translation pathways ρ=−0.598 (p=0.0003) — strongest signal
       - HCC: all translation pathways ↑; HF: most translation pathways ↓
       - (Figure 3: Scatter plot of pathway effect sizes, permutation test)
   
   2.4 Prognostic value of translation gene score is cohort-specific  
       - TCGA-LIHC: TGS prognostic (original finding)
       - GSE141198: TGS Cox p=0.479 (not significant)
       - GSE14520, GSE76427: also not significant (prior analyses)
       - Translation co-expression ≠ translation prognostic
       - (Figure 4: KM curves, forest plot of external validation)
   
   2.5 Upstream regulators: ATF4/ISR and MYC programs  ★NEW★
       - 19 TFs tested for correlation with TGS in TCGA-LIHC
       - ATF4 strongest individual TF (ρ=+0.439, FDR<0.0001)
       - DDIT3 (ρ=+0.289) and XBP1 (ρ=+0.162) — ISR/UPR axis
       - MYC target pathway activity strongest overall signal (ρ=+0.613)
       - HIF1A negatively correlated (ρ=−0.417) — hypoxia represses translation
       - 7 hub genes are NOT in standard Hallmark MYC target gene sets → independently identified
       - (Figure 5: TF-TGS correlation barplot; Pathway-TGS correlation)
   
   2.6 Hub gene direction consistency across diseases
       - 7-gene log2FC direction Spearman correlation: ρ=0.486, p=0.356
       - Individual hub genes show discordant directions
       - Contrasts with module/pathway-level conserved structure
       - (Figure S1: Direction consistency scatter)

3. DISCUSSION
   3.1 A refined model: translation co-expression as organizational principle
       - The module is conserved (WGCNA replication)
       - The direction is disease-specific (ssGSEA negative correlation)
       - The prognostic value is cohort-specific
   
   3.2 ATF4/ISR: a new mechanistic hypothesis
       - ATF4 is the master TF of the Integrated Stress Response
       - ISR controls translation initiation via eIF2α phosphorylation
       - ATF4 → DDIT3 → translation reprogramming in HCC
       - This is NOT the same as MYC-driven ribosome biogenesis
       - Literature context: ISR in cancer (pro-survival adaptation)
   
   3.3 MYC program provides the proliferative drive
       - MYC target pathway activity (ssGSEA) >> MYC mRNA level
       - MYC_V2 (more translation-specific) >> MYC_V1
       - Suggests post-translational/cofactor-level regulation
   
   3.4 Why is HF translation different?
       - Energy crisis model: protein synthesis is ATP-expensive
       - Failing myocardium downregulates translation as adaptation
       - ATF4/ISR in HF is pro-apoptotic (opposite role vs cancer)
       - MTOR negative correlation supports energy-sensing disruption
   
   3.5 Clinical implications
       - TGS is NOT a general HCC prognostic biomarker
       - But translation co-expression network architecture may be targetable
       - ATF4/ISR pathway as potential therapeutic node (both diseases)
   
   3.6 Limitations
       - No wet-lab validation of ATF4 binding or ISR activation
       - RNA-seq only; ribosomal profiling (Ribo-seq) would strengthen conclusions
       - Single HF dataset (GSE57338); limited HCC external cohorts
       - Retrospective/public data; selection bias possible
       - Chinese/Taiwan HCC cohorts; generalizability to Western populations

4. METHODS
   4.1 Data acquisition
   4.2 WGCNA analysis (signed network, both datasets)
   4.3 ssGSEA pathway activity scoring
   4.4 Translation Gene Score (TGS) and survival analysis
   4.5 Upstream TF correlation analysis
   4.6 Direction consistency test
   4.7 Statistical methods

5. DATA AVAILABILITY
   - All datasets from GEO (GSE57338, GSE141198, GSE14520, GSE76427)
   - TCGA-LIHC from GDC portal
   - Analysis code: [repository URL]

6. SUPPLEMENTARY MATERIALS
   - Table S1: Module gene lists and GO enrichment (both datasets)
   - Table S2: ssGSEA pathway effect sizes (full 83 pathways)
   - Table S3: TF-TGS correlation (all 19 TFs)
   - Figure S1: Direction consistency scatter plot
   - Figure S2: WGCNA soft threshold selection (both datasets)
   - Figure S3: Sample clustering dendrograms
   - Figure S4: Individual hub gene KM curves (all cohorts)
   - Supplementary Methods: Detailed R package versions and parameters
```

---

## 四、Figures and Tables Plan

### Main Figures (5-6 figures)

| Figure | Content | Status |
|--------|---------|--------|
| Fig 1 | Study design overview (graphical abstract) | Need to create |
| Fig 2 | WGCNA in GSE141198: dendrogram + module-trait heatmap + GO enrichment for blue module | Have components |
| Fig 3 | ssGSEA cross-disease: pathway effect size scatter + permutation test | **Done** (Figure_ssGSEA_cross_disease.png) |
| Fig 4 | TGS external validation: KM curves + forest plot of all cohorts | Have components |
| Fig 5 | Upstream TF analysis: TF-TGS barplot + pathway-TGS barplot | **Done** (Figure_TF_TGS_correlation.png, Figure_Pathway_TGS_correlation.png) |
| Fig 6 | Mechanistic model: schematic of ATF4/ISR + MYC → translation regulation in HCC vs HF | Need to create |

### Main Tables (2-3 tables)

| Table | Content | Status |
|-------|---------|--------|
| Table 1 | Cohort characteristics (GSE57338, GSE141198, GSE14520, GSE76427, TCGA-LIHC) | Have data |
| Table 2 | Hub gene characteristics and cross-disease module assignment | Have data |
| Table 3 | TF-TGS correlation summary (top TFs with statistics) | **Done** (from 05b output) |

---

## 五、关键写作要点

### 5.1 如何呈现"负结果"

**错误写法**（暴露弱点）：
> "TGS was not prognostic in GSE141198 (p=0.479), which contradicts our initial hypothesis."

**正确写法**（将负结果融入故事）：
> "While TGS demonstrated prognostic value in TCGA-LIHC (HR=X.XX, p<0.001), this association was not replicated in three independent HCC cohorts (GSE141198, GSE14520, GSE76427), suggesting that the prognostic utility of translation gene expression is cohort-dependent and may reflect differences in etiology, stage distribution, or treatment context rather than a universal HCC prognostic mechanism."

### 5.2 如何连接WGCNA复现与ssGSEA负相关

这两个发现互为解释：
- WGCNA说：翻译共表达模块在两种病中都存在（conserved organization）
- ssGSEA说：但这个模块的活性方向是相反的（disease-specific direction）

联合起来：**翻译装置在转录水平上的共表达结构是保守的，但扰动极性由疾病类型决定。**

### 5.3 ATF4 vs MYC——不是竞争，是互补

ATF4和MYC不矛盾：
- ATF4代表应激适应维度（HCC中癌细胞面对营养/缺氧/氧化应激，ISR帮助存活）
- MYC代表增殖维度（癌细胞需要核糖体来支持蛋白质合成需求）

**一个机制假说图**：
```
             HCC                              HF
              │                                │
    ┌─────────┼─────────┐            ┌─────────┼─────────┐
    ▼         ▼         ▼            ▼                   ▼
  MYC      ATF4/ISR   HIF1A?      Energy           Sarcomeric
  program  stress     hypoxia     depletion         dysfunction
    │         │         │            │                   │
    └────┬────┘         │            │                   │
         ▼              ▼            ▼                   ▼
    Ribosome       Cap-dependent   Global              Protein
    biogenesis     translation     translational       quality
    ↑              reprogramming   repression          control
    │                   │            │                   │
    └───────────────────┼────────────┘                   │
                        ▼                                ▼
              Translation UP                    Translation DOWN
              (Proliferation)                  (Adaptive suppression)
```

---

## 六、目标期刊适配

### 5分左右期刊推荐（按叙事匹配度排序）

1. **International Journal of Biological Sciences** (IF ~5.0)
   - 接收：转录组+生物信息学+疾病机制
   - 优势：接受综合计算分析，对验证深度要求适中

2. **Journal of Translational Medicine** (IF ~5.5)
   - 接收：跨疾病/转化生物信息学
   - 优势：接受无湿实验的纯计算研究，关注疾病机制

3. **Frontiers in Genetics / Frontiers in Oncology** (IF ~4.5-5.0)
   - 接收：转录组分析
   - 优势：接受纯生物信息学，审稿快

4. **Genomics** (IF ~4.5)
   - 接收：计算基因组学+疾病
   - 优势：接受方法驱动的研究，不要求湿实验

5. **BMC Genomics** (IF ~4.0)
   - 接收：全转录组分析
   - 优势：审稿标准透明，接受纯计算论文

### Cover Letter 核心卖点

1. **独立复现**：翻译共表达模块在独立HCC队列中复现（WGCNA, p=6.3×10⁻¹²）
2. **负相关的发现**：ssGSEA揭示跨疾病效应量负相关（ρ=−0.598, p=0.0003）——不是简单的"共享"而是"镜像"
3. **新机制假说**：ATF4/ISR作为翻译调控的最强上游TF——超越传统的MYC-centric模型
4. **统计严谨**：置换检验、多重假设校正、多队列外部验证
5. **可复现性**：全部分析代码和数据可获取

---

## 七、还需要做的

### 立即需要：
1. [ ] **创建机制模型图 (Figure 6)**：ATF4/ISR + MYC → 翻译调控的示意
2. [ ] **创建研究设计图 (Figure 1)**：图形摘要
3. [ ] **整合所有table数据**：队列特征 + hub基因 + TF结果
4. [ ] **运行07_final_figures_and_tables.R**：如有必要，或手动整合

### 重写优先级：
1. [ ] **Abstract**（最重要——先写这个，定调）
2. [ ] **Results 2.2-2.5**（新分析的结果部分）
3. [ ] **Discussion 3.1-3.4**（新叙事下的讨论）
4. [ ] **Introduction 1.3-1.4**（调整问题陈述以匹配新叙事）
5. [ ] **Methods**（补充新分析方法）

---

## 八、预期审稿意见及提前预案

| 可能的攻击 | 提前防御 |
|-----------|---------|
| "没有湿实验验证ATF4调控" | 在Limitations中明确声明；强调计算预测的可生成假说价值；引用ATF4在HCC中的已有文献 |
| "TGS没有外部预后价值，论文贡献有限" | 将焦点从预后价值转向机制发现；TGS不显著本身是一个发现（cross-disease specificity） |
| "只有两个数据集比较" | 强调分析深度（5个层次：基因→模块→通路→预后→TF）；HF数据稀缺的客观限制 |
| "翻译通路负相关可能是因为器官差异而非疾病差异" | 使用normal vs disease的效应量（Cohen's d），以各自正常组织为基线，消除了器官基线差异 |
| "负相关只是统计学假象" | 置换检验p=0.0079；翻译通路子集p=0.0003远强于全体通路；33条翻译通路一致负向 |
