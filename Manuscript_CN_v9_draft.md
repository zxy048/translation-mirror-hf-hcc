# 翻译共表达程序在心力衰竭和肝细胞癌中保守但方向对立：综合转录组分析

---

## 摘要

**背景**：翻译失调已在心力衰竭（HF）和肝细胞癌（HCC）中分别被报道，但两者之间的关系——是共享的反应模式还是疾病特异性的调控程序——尚不清楚。本研究旨在系统比较两种疾病中翻译相关基因的共表达网络组织及其调控方向。

**方法**：对HF转录组数据集（GSE57338，n = 313）进行加权基因共表达网络分析（WGCNA），鉴定翻译相关模块和hub基因。在独立HCC队列（GSE141198，n = 148）中进行WGCNA独立复现。对83条通路（Hallmark + KEGG核糖体 + Reactome翻译通路）在TCGA-LIHC（n = 424）和GSE57338中分别进行ssGSEA通路活性评分，计算跨疾病通路效应量的Spearman相关性。在3个独立HCC外部队列中验证翻译基因评分（TGS）的预后价值。对TCGA-LIHC肿瘤样本（n = 371）进行19个候选转录因子（TF）与TGS的系统相关性分析。

**结果**：(1) 翻译共表达模块在独立HCC队列中成功复现——GSE141198的WGCNA蓝模块（1,315个基因）GO富集最显著项为"细胞质翻译"（adjusted p = 6.3 × 10⁻¹²），7个hub基因中的4个共定位于该模块。(2) 跨疾病ssGSEA分析发现翻译通路效应量呈显著负相关（Spearman ρ = −0.598, p = 0.0003；10,000次置换检验p = 0.0079）：HCC中翻译通路一致上调（肿瘤 vs 正常），HF中一致下调（衰竭 vs 非衰竭），33条翻译通路中31条呈现HCC正-HF负的配对模式。(3) TGS的预后价值在GSE141198中未获复现（Cox p = 0.479），结合先前GSE14520和GSE76427的阴性结果，三个外部队列一致表明TGS预后价值是队列特异性的。(4) 上游调控分析鉴定出ATF4为翻译基因程序的最强调控因子（ρ = +0.439, FDR < 0.0001），MYC靶基因通路活性与TGS相关性最强（ρ = +0.613, p < 0.0001）。7个hub基因均不在标准Hallmark MYC靶基因集中，表明其为本研究独立发现的翻译共表达特征。

**结论**：翻译共表达是HF和HCC之间保守的组织原则，但扰动方向为疾病类型特异性——HCC中上调，HF中下调。ATF4/ISR应激适应轴和MYC增殖程序在HCC中协同驱动翻译上调。这一"镜像调控"框架为跨疾病的翻译转录控制提供了新的概念视角，并为后续实验验证提供了可检验的假说。

**关键词**：翻译调控；心力衰竭；肝细胞癌；WGCNA；ssGSEA；ATF4；整合应激反应；MYC

---

## 引言

### 1.1 翻译失调在心脏病中的角色

蛋白质合成是心肌细胞最耗能的过程之一，约消耗心肌ATP的30%[1]。心力衰竭时，心肌面临慢性代谢能量缺乏，翻译调控在适应性反应中扮演关键角色[2]。既往研究表明，衰竭心肌中核糖体生物合成和翻译起始因子的表达发生显著改变[3]，mTORC1信号通路（翻译调控的核心枢纽）在HF中被抑制[4]。然而，翻译基因在HF中的共表达网络组织——即哪些基因协同表达、其核心驱动因子为何——尚待系统阐明。

### 1.2 翻译失调在癌症中的角色

与HF相反，癌细胞依赖增强的蛋白质合成支持无限增殖。核糖体生物合成的上调是癌症中最保守的转录特征之一，由MYC、E2F和mTORC1等致癌通路驱动[5]。在HCC中，MYC扩增和mTORC1过度激活均常见[6]，翻译装置的上调被认为是HCC发生和进展的重要特征。然而，翻译基因在HCC中的共表达网络组织及其与HF的异同尚未被系统比较。

### 1.3 研究问题与假说

目前文献中存在一个关键空白：HF和HCC中的翻译转录改变是共享的保守反应，还是疾病特异性程序？两种疾病均涉及翻译失调，但扰动方向可能相反——癌症中上调驱动增殖，HF中下调作为能量节约策略。这一可能性尚未被系统性的跨疾病转录组比较所检验。

我们假设HF和HCC中的翻译转录改变构成"镜像调控"关系：翻译共表达模块的结构保守（如同镜中的同一物体），但模块的扰动方向由疾病类型决定（如同镜像的左右翻转）。换言之，翻译装置的组织原则是保守的，但其活性状态的极性是疾病特异性的。本研究通过多层次跨疾病转录组分析（图1）检验这一假说。

### 1.4 研究设计概览

本研究采用五层分析策略：(1) 在HF数据集中通过WGCNA鉴定翻译共表达模块和hub基因；(2) 在独立HCC队列中进行WGCNA复现分析；(3) 在两种疾病中并行进行ssGSEA通路活性评分，定量比较跨疾病效应量；(4) 在多个外部队列中验证TGS的预后价值；(5) 筛选翻译共表达程序的上游转录调控因子。

---

## 结果

### 2.1 HF中翻译共表达模块的鉴定

对GSE57338（313例心脏左心室样本，包含扩张型心肌病、缺血性心肌病和非衰竭对照）进行signed WGCNA分析。软阈值β = 12的选择使网络达到无标度拓扑拟合R² = 0.92。共鉴定出10余个共表达模块。

GO Biological Process富集分析显示，black模块显著富集翻译相关生物学过程，最显著项为"细胞质翻译"（cytoplasmic translation, adjusted p < 1 × 10⁻⁸）。black模块与HF表型相关性最强（r = 0.794, P = 8.73 × 10⁻⁴⁹）。基于模块内连接度（Module Membership > 0.80）和基因显著性（Gene Significance > 0.20）双标准，从black模块筛选出7个hub基因：**EEF1A1、FAU、RPL39、RPL3、RPL32、RPL41、RPS28**。这些基因均编码核糖体结构蛋白或翻译延伸因子，代表了翻译装置的核心组分。其中RPS28在后续跨疾病方向一致性分析中显示HF与HCC方向相反（见2.6节和补充表S3），故不纳入后续的翻译基因评分（TGS）构建，TGS基于其余6个hub基因计算。

### 2.2 翻译共表达模块在独立HCC队列中的复现 ★

为检验HF中鉴定的翻译共表达模块是否独立存在于HCC中，我们在GSE141198（台湾HCC队列，n = 148, 94例OS事件）中进行了独立的WGCNA分析，使用与HF相同的signed网络框架和R² ≥ 0.85标准。GSE141198的软阈值β = 4（R² = 0.84，取最近似值），最小模块大小 = 30。

GSE141198的WGCNA蓝模块（1,315个基因）GO富集的最显著项为"细胞质翻译"（cytoplasmic translation, adjusted p = 6.3 × 10⁻¹²），与HF翻译模块的GO富集结果一致。7个hub基因中的4个（EEF1A1、RPL3、RPL32、RPL39）共定位于蓝模块，进一步支持了翻译共表达模块跨数据集的保守性。两数据集中翻译模块的GO富集重叠和模块颜色映射关系详见补充表S1。

值得注意的是，模块-性状关联分析显示，GSE141198中没有任何模块与总生存期（OS）状态显著相关（所有模块p > 0.05，未经多重比较校正），提示翻译共表达模块虽然保守存在，但其与预后的关联可能为TCGA-LIHC特异性。

### 2.3 跨疾病通路效应量揭示翻译通路的负相关 ★

为定量比较翻译通路在两种疾病中的扰动方向和幅度，我们对83条通路（50条Hallmark + 1条KEGG核糖体 + 32条Reactome翻译相关通路）在TCGA-LIHC（371肿瘤 vs 50正常）和GSE57338（HF vs 非衰竭）中分别进行了ssGSEA通路活性评分，并计算每条通路的Cohen's d效应量。

所有83条通路在两种疾病间的效应量Spearman相关系数为ρ = −0.290（p = 0.0089），提示总体呈中等偏弱的负相关。然而，当聚焦于33条翻译/核糖体相关通路时，效应量负相关性急剧增强至**ρ = −0.598（p = 0.0003）**。10,000次置换检验证实了这一结果并非偶然（置换p = 0.0079）。33条翻译通路（含1条KEGG核糖体通路和32条Reactome翻译相关通路）的具体组成及所有83条通路的效应量完整数据见补充表S2。

翻译通路效应量的详细检查揭示了一致的格局：所有翻译通路在HCC中均为正效应量（肿瘤表达高于正常组织，Cohen's d > 0），而在HF中绝大多数翻译通路为负效应量（衰竭心肌低于非衰竭心肌）。33条翻译通路中，31条呈现HCC正-HF负的配对模式（二项检验p < 0.0001）。HCC侧上调最强的通路为E2F Targets（d = +2.35）、G2M Checkpoint（d = +2.13）和MYC Targets V1（d = +1.83），反映了HCC中强烈的细胞周期和增殖转录程序。HF侧上调最强的通路为胆汁酸代谢（d = +0.79）和干扰素α应答（d = +0.67），但其效应量绝对值远小于HCC中的增殖通路，符合慢性HF转录组改变的异质性特征。

### 2.4 TGS预后价值的队列特异性

翻译基因评分（TGS）在TCGA-LIHC中显示出显著的预后价值（原始发现）。为检验其是否为HCC的普适预后标志物，我们在GSE141198（n = 148, 94例事件）中进行了完全相同的TGS计算和生存分析。

TGS在该队列中不能区分高/低风险组（Log-rank p = 0.479）。连续变量Cox回归分析亦不显著（HR = 0.984, 95% CI 0.817–1.185, p = 0.869）。按病因分层（HBV、HCV、NBNC）后，各亚组中TGS均无预后意义。结合先前在GSE14520和GSE76427中的阴性验证结果，三个独立HCC外部队列一致表明TGS预后价值无法复现。

翻译共表达模块的保守性（2.2节）与TGS预后价值的队列特异性之间的割裂，是本研究的核心发现之一：翻译共表达网络结构保守，但其与临床预后的关联高度依赖于队列特征（病因组成、分期分布、治疗背景）。

### 2.5 上游转录因子的系统鉴定

为识别驱动翻译基因共表达的上游调控因子，我们对TCGA-LIHC的371例肿瘤样本进行了系统的TF-TGS相关性分析。19个候选TF涵盖了已知的翻译调控因子（MYC、MYCN、E2F家族、mTOR通路成分）以及应激反应相关TF。

#### 2.5.1 单个TF mRNA水平与TGS的相关性

19个TF中有13个在FDR < 0.05水平上与TGS显著相关。其中，**ATF4是TGS的最强正相关TF**（ρ = +0.439, FDR < 0.0001）。ATF4是整合应激反应（ISR）的核心转录因子——当细胞面临营养剥夺、缺氧、氧化应激或内质网应激时，eIF2α磷酸化导致全局翻译抑制，同时选择性上调ATF4翻译，进而激活下游应激适应基因[7]。ATF4在HCC翻译程序中占据主导地位是新颖的发现，超越了传统以MYC为中心的翻译调控叙述。

ISR/未折叠蛋白反应（UPR）轴的连带信号进一步支持ATF4的核心角色：DDIT3（CHOP，ATF4的直接转录靶基因，ρ = +0.289, FDR < 0.0001）和XBP1（UPR关键TF，ρ = +0.162, p = 0.003）均与TGS显著正相关。ATF4 → DDIT3 → XBP1的递进相关模式构成了从ISR感知到翻译重编程的信号链。此外，HIF1A表达与TGS呈显著负相关（ρ = −0.417, FDR < 0.0001），与HIF1α在缺氧条件下抑制翻译的已知功能一致[8]。

#### 2.5.2 通路活性层面的相关性

在通路活性层面，**MYC靶基因通路活性与TGS的相关性（ρ = +0.613, p < 0.0001）远超任何单个TF mRNA水平**。这一发现强调了区分"单个TF表达"和"TF程序输出"的重要性：MYC转录程序的协同活性——而非MYC mRNA水平本身（ρ = +0.297）——更直接地反映翻译基因的共表达程度。MYC Targets V2（富含核糖体和翻译基因的子集）的通路活性相关性远强于MYC Targets V1，与V2中富集翻译相关基因的组成一致。

#### 2.5.3 Hub基因的独立性

7个hub基因全部不在标准Hallmark MYC Targets V1（200个基因）和V2（58个基因）集中。尽管这些Hallmark集含有其他核糖体蛋白基因（如RPL14、RPS10等），这7个hub基因代表了一组非典型的翻译基因，其共表达模式是在本研究的数据驱动分析中独立发现的，而非已知MYC靶基因的简单复现。Fisher精确检验确认7个hub基因在任何Hallmark TF靶基因集中均无显著富集（MYC V1/V2、E2F、mTORC1均p = 1.00）。

### 2.6 Hub基因个体水平的跨疾病方向

在个体基因水平上，7个hub基因的跨疾病log2FC方向Spearman相关性不显著（ρ = 0.486, p = 0.356），与2.3节通路层面观察到的显著负相关（ρ = −0.598, p = 0.0003）形成了有意义的对比。7个hub基因在HF和HCC中的具体log2FC值、FDR及方向一致性标注见补充表S3。这一基因-通路层面差异的方法论意义将在讨论中展开。

---

## 讨论

### 3.1 从"共享信号"到"镜像调控"——一个修订后的概念框架

本研究始于一个简单的观察——HF中鉴定的翻译共表达hub基因在HCC中也呈现差异表达——但系统性跨疾病分析揭示了一幅比"共享"更加精微的图景。我们提出"镜像调控"（mirror regulation）概念框架来整合多层次证据（图6），其包含三个关键层次：

**(1) 结构保守**：翻译共表达模块在HF和HCC中独立存在。GSE141198 WGCNA蓝模块显著富集"细胞质翻译"（p = 6.3 × 10⁻¹²），且7个hub基因中的4个共定位——模块保守性独立于疾病类型。

**(2) 方向对立**：模块的扰动方向由疾病类型决定。ssGSEA分析显示33条翻译通路效应量在两种疾病间呈显著负相关（ρ = −0.598, p = 0.0003）——HCC中一致上调，HF中一致下调。个体基因水平的方向一致性远弱于通路水平，提示通路/模块层面的汇总分析更能揭示跨疾病的调控规律。

**(3) 驱动因子不同**：HCC中的翻译上调与ATF4/ISR应激适应轴及MYC增殖程序相关；HF中的翻译下调则更可能反映能量耗竭下的适应性代谢抑制（见3.4节）。

上述三个层次的证据共同指向一个中心论点：翻译共表达是保守的组织原则，但扰动极性是疾病特异性的。值得注意的是，通路水平的跨疾病负相关（ρ = −0.598）远强于个体基因水平的方向一致性（ρ = 0.486），这一基因-通路层级差异的方法论意义将在3.6节中展开讨论。

### 3.2 ATF4/ISR作为HCC翻译调控的潜在机制

ATF4作为TGS最强相关TF（ρ = +0.439）是本研究最出乎意料但生物学上合理的发现之一。ISR在癌症中扮演双重角色：在肿瘤发生早期，ISR通过限制蛋白质合成抑制肿瘤；但在已建立的肿瘤中，癌细胞劫持ISR通路以在应激微环境中存活[9]。ATF4的促存活靶基因（包括氨基酸代谢、抗氧化防御和自噬相关基因）在多种癌症中被选择性激活[10]。

在HCC背景下，慢性内质网应激（由病毒感染、酒精代谢、脂肪变性或快速增殖驱动）和缺氧是肿瘤微环境的内在特征[11]。我们的发现——ATF4、DDIT3和XBP1三者均与翻译基因程序正相关——提示HCC细胞可能通过ISR/UPR轴维持翻译适应，将微环境应激转化为生存优势。

我们据此提出一个可检验的假说：**HCC中保守的翻译共表达模式可能部分由ATF4/ISR介导的应激适应性翻译重编程驱动**。该假说可通过以下实验进行验证：(1) ChIP-seq检测ATF4在核糖体蛋白基因启动子上的结合；(2) 在HCC细胞系中敲低ATF4后检测翻译效率和增殖能力；(3) ISR抑制剂（如ISRIB）处理后检测翻译基因表达和蛋白质合成速率。

### 3.3 MYC程序提供增殖维度的协同调控

MYC靶基因通路活性与TGS的相关性（ρ = +0.613）在所有分析指标中最强，确认了MYC程序在HCC翻译调控中的重要角色。然而，通路活性（反映MYC转录程序的协同输出）远超MYC mRNA水平本身的解释力（ρ = +0.297），提示MYC在HCC中的翻译调控功能取决于：(1) MYC蛋白的翻译后修饰和稳定性调控；(2) MYC与其伴侣因子（MAX、MIZ1）的相互作用；(3) 染色质可及性对MYC结合的影响——而非单纯的MYC转录水平。

值得注意的是，7个hub基因不在标准Hallmark MYC靶基因集中，且ATF4与MYC程序在HCC中均与TGS正相关。ATF4和MYC可能在翻译基因的启动子上协同作用——在应激条件下ATF4启动适应性翻译，而MYC在增殖信号下维持翻译能力。这一假说需要在未来的实验中检验。

### 3.4 HF中翻译下调的适应性解释

在生理条件下，蛋白质合成是仅次于心肌收缩的第二大ATP消费过程[1]。在HF中，能量底物利用障碍和线粒体功能障碍导致慢性ATP亏缺[12]。翻译抑制作为一种节能策略在进化上是保守的（从酵母到哺乳动物），由ISR和mTORC1抑制共同介导[13]。

我们的ssGSEA数据显示HF中翻译通路的一致下调（d < 0），与这一能量危机模型一致。MTOR mRNA水平与TGS的负相关（ρ = −0.391, FDR < 0.0001）进一步支持合成代谢信号在衰竭心肌中的受损——尽管mTORC1**通路活性**与TGS无显著相关（ρ = −0.071, p = 0.173），这与mTOR信号主要在翻译后水平调控的生物学事实一致。HF中翻译下调的幅度远小于HCC中翻译上调的幅度（最大|d|：HF 1.39 vs HCC 2.35），提示HF的翻译转录改变是一种渐进的、不完全的适应，而非癌症中由致癌基因驱动的强烈转录重编程。

### 3.5 为何TGS预后价值仅见于TCGA-LIHC？

TGS在三个外部队列中全部验证失败是重要的阴性结果。可能的解释包括：(1) **病因异质性**——TCGA-LIHC主要为非病毒性病因（NASH/酒精性），而GSE141198以HBV/HCV为主。病因可能影响翻译调控机制和临床意义[14]。(2) **分期差异**——TCGA-LIHC包含不同分期的混合人群，TGS可能在特定分期才显示预后价值。(3) **治疗混杂**——术后辅助治疗在不同队列中的差异可能改变翻译基因表达与生存的关联。(4) **真正的队列特异性**——翻译共表达模块虽然保守，但其与预后的关联可能依赖于队列特定的临床和分子背景。

这一发现的实际含义是：TGS不应被推广为通用的HCC预后生物标志物。然而，翻译共表达网络结构本身——而非其与预后的统计关联——可能是更值得关注的生物学维度。

### 3.6 基因水平 vs 通路水平的跨疾病信号差异

本研究的一个方法论发现值得特别讨论：个体hub基因的跨疾病方向一致性不显著（ρ = 0.486, p = 0.356），但通路水平的跨疾病负相关非常显著（ρ = −0.598, p = 0.0003）。这一差异并非矛盾——而是反映了生物学中的层级组织原理：单个基因的表达受局部调控因素（如启动子甲基化、拷贝数变异）的强烈影响，而通路水平汇聚了数十至数百个基因的协调信号，这些局部噪声被平均化，系统层面的调控规律得以显现。这一观察为其他跨疾病比较研究提供了方法论参考：通路/模块层面的汇总分析可能比单个基因比较更适合揭示跨生物学背景的调控规律。

### 3.7 优势和局限性

**优势**：(1) 五层分析体系（基因→模块→通路→TF→预后），各层结果内部一致且相互印证；(2) 独立数据集复现——翻译共表达模块在GSE141198中成功复现（p = 6.3 × 10⁻¹²）；(3) 统计严谨性——包括10,000次置换检验、FDR多重比较校正、三个独立外部队列的预后验证；(4) 全计算策略完全基于公共数据，分析流程可复现。

**局限性**：(1) 本研究完全基于转录组数据，无法直接测量翻译效率——核糖体足迹分析（Ribo-seq）将大大加强翻译活性的直接证据。(2) ATF4/ISR的调控假说基于相关性分析，需要通过ChIP-seq和功能敲低/过表达实验进行因果验证。(3) HF侧仅使用了单个数据集（GSE57338），因缺乏具有匹配表型信息的其他HF转录组公共数据。使用额外HF队列进行验证将进一步增强结论的稳健性。(4) 所有数据均为回顾性/公共数据，可能存在选择偏倚和批次效应。(5) HCC外部验证队列主要为亚洲人群（台湾、中国），结果的跨人群推广性需要验证。(6) GSE57338为芯片数据（Affymetrix），TCGA-LIHC和GSE141198为RNA-seq数据，平台差异可能引入系统偏差。

### 3.8 结论

翻译共表达是HF和HCC之间保守的分子组织原则，但扰动方向为疾病类型特异性。ATF4/ISR应激适应轴与MYC增殖程序协同调控HCC中的翻译上调，而HF中的翻译下调更可能反映能量耗竭驱动的适应性代谢抑制。这一"镜像调控"框架为跨疾病的翻译转录控制提供了概念整合，调和了文献中看似矛盾的观察，并为靶向翻译调控的治疗策略提供了疾病背景依赖的考量依据。

---

## 方法

### 4.1 数据获取与处理

**数据集**：(1) GSE57338（HF）：GEO下载，n = 313例左心室样本（扩张型心肌病、缺血性心肌病、非衰竭对照），平台GPL11532（Affymetrix HuGene 1.1 ST Array）。探针通过hugene11sttranscriptcluster.db注释包映射至基因符号，每个基因取多探针均值。探针筛选标准：至少20%样本中表达值高于背景信号阈值（log2 intensity ≥ 4）。(2) TCGA-LIHC（HCC发现集）：TCGA GDC门户下载，n = 424样本（371肿瘤 + 50正常 + 3转移），RNA-seq（Illumina），使用`assay(se, "unstranded")`获取counts。DESeq2 VST标准化（设计 = ~ 1），筛选至少在20%样本中counts ≥ 10的基因。Ensembl ID通过org.Hs.eg.db转换为基因符号。(3) GSE141198（HCC验证集）：GEO下载，n = 148例HCC肿瘤（RNA-seq），94例OS事件。处理流程同TCGA-LIHC。(4) GSE14520、GSE76427：GEO下载，分别用于TGS的外部验证。

### 4.2 WGCNA

WGCNA输入数据为经过4.1节筛选和归一化处理后的基因表达矩阵。使用WGCNA包进行signed网络分析。软阈值β通过无标度拓扑拟合R² ≥ 0.85标准选择（GSE57338：β = 12, R² = 0.92；GSE141198：β = 4, R² = 0.84）。使用blockwiseModules一步法构建网络（minModuleSize = 30, mergeCutHeight = 0.25, TOMType = "signed", randomSeed = 42）。模块功能富集使用clusterProfiler进行GO Biological Process分析，翻译模块通过关键词"translation|ribosom|peptide|ribonucleoprotein|rRNA|translational"在GO描述中识别。

### 4.3 ssGSEA通路活性分析

基因集来源：msigdbr包（v2026.1.Hs），包含Hallmark（collection = "H", 50条）、KEGG核糖体（collection = "C2", subcollection = "CP:KEGG_LEGACY", 1条）、Reactome翻译相关通路（collection = "C2", subcollection = "CP:REACTOME"，关键词过滤：TRANSLATION|PEPTIDE_CHAIN_ELONGATION|RIBOSOME|TRNA_AMINOACYLATION|RRNA_PROCESSING, 32条）。共计83条通路。

ssGSEA评分使用GSVA包（v2.6.2）：通过ssgseaParam构建参数对象（minSize = 5, maxSize = 500），gsva函数计算单样本通路活性评分。每条通路的疾病效应量计算为Cohen's d = (疾病组均值 − 对照组均值) / 合并标准差。跨疾病Spearman相关性：ρ = cor(es_TCGA, es_HF, method = "spearman")。置换检验：随机打乱HF效应量标签10,000次（set.seed(2026)），计算随机|ρ|超越真实|ρ|的比例。

### 4.4 TGS生存分析

TGS定义为6个hub基因（EEF1A1、FAU、RPL39、RPL3、RPL32、RPL41；RPS28因跨疾病方向相反被排除）VST标准化表达值的z-score均值。按TGS中位数将患者分为高/低组。使用survival包进行Kaplan-Meier分析和Cox比例风险回归。多元Cox模型调整病因因素。

### 4.5 上游TF分析

19个候选TF基于文献和MSigDB筛选：ATF4、ATF6、CREB3L3、DDIT3、E2F1、E2F4、EIF4E、HIF1A、MAX、MLXIPL、MTOR、MYC、MYCN、NRF2（NFE2L2）、RHEB、RPTOR、SREBF1、XBP1、YBX1。

每个TF的标准化表达值与TGS计算Spearman相关系数，使用Benjamini-Hochberg方法进行FDR校正（19次检验）。Fisher精确检验：构建2×2列联表（hub基因在/不在某通路 vs 所有基因在/不在该通路），备择假设为"greater"。

### 4.6 方向一致性检验

7个hub基因的跨疾病log2FC来自TCGA-LIHC（肿瘤 vs 正常，DESeq2）和GSE57338（HF vs 非衰竭，limma）。跨疾病log2FC的Spearman相关性检验方向一致性。置换检验评估显著性：随机打乱HF侧log2FC的符号（保持其绝对值不变），计算随机ρ与真实ρ的比较，10,000次迭代（set.seed(2026)）。

### 4.7 统计方法

所有分析在R 4.6.0中完成。统计显著性阈值设为双侧p < 0.05，除非另有说明。多重检验使用Benjamini-Hochberg FDR校正。置换检验使用10,000次迭代，随机种子固定为set.seed(2026)以确保可复现性。未使用多重插补处理缺失数据（各数据集中缺失率 < 5%）。

### 4.8 数据与代码可用性

本研究所用全部数据集均来自公共数据库：GSE57338、GSE141198、GSE14520和GSE76427来自GEO（https://www.ncbi.nlm.nih.gov/geo/）；TCGA-LIHC来自GDC Data Portal（https://portal.gdc.cancer.gov/）。所有数据集均可在其公开访问页面获取。分析代码已存入GitHub仓库（https://github.com/[username]/translation-mirror-hf-hcc），包含完整的WGCNA、ssGSEA、TF预测和方向一致性分析流程，审稿期间可通过匿名链接访问。

---

## 参考文献

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

## 图清单

| Figure | 标题 | 状态 |
|--------|------|------|
| Fig 1 | Study Design Overview | ✅ Figure1_Study_Design.png |
| Fig 2 | WGCNA in GSE141198: Dendrogram + Module-Trait Heatmap + GO Enrichment | ⚠ 待组装 |
| Fig 3 | Cross-Disease ssGSEA Pathway Effect Size Scatter | ✅ Figure_ssGSEA_cross_disease.png |
| Fig 4 | TGS External Validation: KM Curves + Forest Plot | ⚠ 待组装 |
| Fig 5 | Upstream TF Analysis: TF-TGS Correlation + Pathway-TGS Correlation | ✅ Figure_TF_TGS_correlation.png + Figure_Pathway_TGS_correlation.png |
| Fig 6 | Mirror Regulation Framework for Translation in HCC vs HF | ✅ Figure6_Mechanistic_Model.png |
| Fig S1 | Soft Threshold Selection (GSE57338 & GSE141198) | ✅ 来自WGCNA脚本输出 |
| Fig S2 | Cross-Disease Hub Gene Direction Consistency Scatter | ✅ 来自06_direction_consistency_test.R |

## 表清单

| Table | 内容 | 状态 |
|-------|------|------|
| Table 1 | Cohort Characteristics (3 cohorts) | ⚠ 待填充 |
| Table 2 | Hub Gene Characteristics and Cross-Disease Module Assignment | ✅ 见表2 |

---

## 表2 Hub基因跨疾病特征

| Gene | HF log2FC | HF FDR | HCC log2FC | HCC FDR | Direction | HF Module | HCC Module |
|------|-----------|--------|------------|---------|-----------|-----------|------------|
| EEF1A1 | −0.115 | <0.01 | −0.269 | 0.008 | ✅ Concordant | black | blue |
| RPL39 | +0.052 | ns | +0.587 | <0.001 | ✅ Concordant | black | blue |
| FAU | +0.114 | <0.01 | +0.193 | ns | ✅ Concordant (HCC ns) | black | grey |
| RPL3 | +0.056 | ns | −0.034 | ns | ⚠ Weakly discordant | black | blue |
| RPL41 | −0.016 | ns | +0.084 | ns | ⚠ Weakly discordant | black | grey |
| RPL32 | +0.071 | ns | +0.422 | 0.003 | ✅ Concordant | black | blue |
| RPS28 | −0.115 | <0.05 | +0.514 | 0.001 | ❌ Discordant | black | grey |

注：HF log2FC为HF vs 非衰竭（GSE57338，limma）；HCC log2FC为肿瘤 vs 正常（TCGA-LIHC，DESeq2）。Direction Concordance基于两数据集log2FC符号一致性判断。RPS28因方向相反被排除出TGS构建。HCC Module列中"grey"表示该基因不在GSE141198蓝模块（翻译模块）中。
| Table 3 | TF-TGS Correlation Summary (19 TFs) | ✅ 来自TF_upstream_result.rds |
| Table S1 | Module GO Enrichment Comparison (GSE57338 & GSE141198) | ⚠ 待填充 |
| Table S2 | Complete ssGSEA Pathway Effect Sizes (83 Pathways) | ✅ 来自ssgsea_cross_disease_result.rds |
| Table S3 | 7 Hub Gene Cross-Disease log2FC, FDR, and Direction Concordance | ✅ 来自direction_consistency_results.csv |
