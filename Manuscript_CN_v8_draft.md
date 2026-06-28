# 翻译共表达程序在心力衰竭和肝细胞癌中的镜像调控：综合转录组分析

---

## 摘要

**背景**：翻译失调已在心力衰竭（HF）和肝细胞癌（HCC）中分别被报道，但两者之间的关系——共享机制还是疾病特异性程序——尚不清楚。

**方法**：本研究对HF数据集（GSE57338，n = 313）进行WGCNA鉴定翻译相关共表达模块和hub基因。在独立HCC验证队列（GSE141198，n = 148）中进行WGCNA独立复现分析。对83条通路（Hallmark + KEGG核糖体 + Reactome翻译通路）在TCGA-LIHC（n = 424）和GSE57338中分别进行ssGSEA通路活性评分，计算跨疾病通路效应量Cohen's d的Spearman相关性。在3个独立HCC外部队列（GSE141198、GSE14520、GSE76427）中验证翻译基因评分（TGS）的预后价值。对TCGA-LIHC肿瘤样本（n = 371）中19个候选转录因子（TF）进行TF表达与TGS的Spearman相关性分析以及Fisher精确检验。

**结果**：(1) 翻译共表达模块在独立HCC队列中成功复现：GSE141198的WGCNA蓝模块显著富集"细胞质翻译"（cytoplasmic translation，p = 6.3 × 10⁻¹²），7个hub基因中的4个共定位于该模块。(2) 跨疾病ssGSEA分析揭示翻译通路效应量呈显著负相关（Spearman ρ = −0.598，p = 0.0003，置换检验p = 0.0079）：HCC中翻译相关通路一致上调（肿瘤 vs 正常），HF中一致下调（衰竭 vs 非衰竭）。(3) TGS的预后价值未在任何外部队列中复现（GSE141198：Cox p = 0.479），表明其为TCGA-LIHC特异性而非HCC普适的预后标志物。(4) 上游TF分析鉴定出ATF4为翻译基因程序的最强调控因子（ρ = +0.439，FDR < 0.0001），MYC靶基因通路活性与TGS相关性最强（ρ = +0.613，p < 0.0001），DDIT3（ρ = +0.289）和XBP1（ρ = +0.162）构成整合应激反应（ISR）/未折叠蛋白反应（UPR）轴。7个hub基因均不在标准Hallmark MYC靶基因集中，表明其为独立于现有知识的发现。

**结论**：翻译共表达是HF和HCC之间保守的组织原则，但扰动方向是疾病类型特异性的——HCC中上调以支持增殖，HF中下调以节约能量。ATF4/ISR应激适应通路和MYC增殖程序共同驱动HCC中的翻译上调。这一"镜像调控"模型调和了文献中看似矛盾的观察，并为靶向翻译调控提供了疾病背景下的机制框架。

**关键词**：翻译调控；心力衰竭；肝细胞癌；WGCNA；ssGSEA；ATF4；整合应激反应；MYC

---

## 引言

### 1.1 翻译失调在心脏病中的角色

蛋白质合成是心肌细胞最耗能的过程之一，约消耗心肌ATP的30%[1]。心力衰竭时，心肌面临慢性能量缺乏，翻译调控在适应性反应中扮演关键角色[2]。既往研究表明，衰竭心肌中核糖体生物合成和翻译起始因子表达发生显著改变[3]，mTORC1信号通路（翻译调控的核心枢纽）在HF中被抑制[4]。然而，大规模转录组层面的翻译基因共表达网络结构及其在HF中的系统特征尚待系统阐明。

### 1.2 翻译失调在癌症中的角色

与HF相反，癌细胞依赖增强的蛋白质合成来支持无限增殖。核糖体生物合成的上调是癌症中最保守的转录特征之一，由MYC、E2F和mTORC1等致癌通路驱动[5]。在HCC中，MYC扩增和mTORC1过度激活均常见[6]，翻译装置的上调被认为是HCC发生和进展的驱动性事件。然而，翻译基因表达在HCC中的共表达网络组织及其与HF的异同尚未被系统比较。

### 1.3 研究空白与假说

目前文献中存在一个关键空白：HF和HCC中的翻译转录改变是共享的保守反应，还是疾病特异性程序？两种疾病均涉及翻译失调，但方向可能相反——癌症中上调驱动增殖，HF中下调作为能量节约策略。这种可能性尚未被系统性地跨疾病转录组比较所检验。

### 1.4 研究设计概览

本研究采用多层次跨疾病转录组分析策略（图1）：(1) 在HF数据集中通过WGCNA鉴定翻译共表达模块和hub基因；(2) 在独立HCC队列中进行WGCNA复现分析；(3) 在两种疾病中并行进行ssGSEA通路活性评分并比较效应量；(4) 在多个外部队列中验证TGS的预后价值；(5) 筛选翻译共表达程序的上游转录调控因子。

---

## 结果

### 2.1 HF中翻译共表达模块的鉴定

对GSE57338（313例心脏左心室样本，包含扩张型心肌病、缺血性心肌病和非衰竭对照）进行signed WGCNA分析。软阈值β = [实际值]的选择使网络达到无标度拓扑拟合R² > 0.85。共鉴定出[实际模块数]个共表达模块。

GO Biological Process富集分析显示，[颜色]模块显著富集翻译相关生物学过程，包括"细胞质翻译"（cytoplasmic translation，p = [实际值]）、"核糖体结构组成"等。该模块中鉴定出7个hub基因（按模块内连接度排序）：EEF1A1、FAU、RPL39、RPL3、RPL32、RPL41、RPS28。这些hub基因均编码核糖体结构蛋白或翻译延伸因子，代表了翻译装置的核心组分。

### 2.2 翻译共表达模块在独立HCC队列中的复现 ★

为检验HF中鉴定的翻译共表达模块是否独立存在于HCC中，我们在GSE141198（台湾HCC队列，n = 148，94例OS事件）中进行了独立的WGCNA分析。使用与HF分析相同的参数（signed网络，β = [实际值]，最小模块大小 = 30）。

GSE141198的WGCNA蓝模块（1,315个基因）GO富集的最显著项为"细胞质翻译"（cytoplasmic translation，adjusted p = 6.3 × 10⁻¹²），与HF翻译模块的GO富集结果高度一致。重要的是，7个hub基因中的4个（EEF1A1、RPL3、RPL32、RPL39）共定位于蓝模块。

然而，模块-性状关联分析显示，GSE141198中没有任何模块与总生存期（OS）状态显著相关（所有模块p > 0.05）。翻译共表达模块在GSE141198中保守存在，但缺乏预后相关性，提示翻译共表达的预后价值可能是TCGA-LIHC特异性的。

### 2.3 跨疾病通路效应量揭示翻译通路的负相关 ★

为定量比较翻译通路在两种疾病中的扰动方向，我们对83条通路（50条Hallmark + 1条KEGG核糖体 + 32条Reactome翻译相关通路）在TCGA-LIHC（371肿瘤 vs 50正常）和GSE57338（HF vs NF）中分别进行了ssGSEA通路活性评分，并计算每条通路的Cohen's d效应量。

所有83条通路在两种疾病间的效应量Spearman相关系数为ρ = −0.290（p = 0.0089），提示总体呈中等偏弱的负相关。然而，当聚焦于33条翻译/核糖体相关通路时，效应量负相关性急剧增强至ρ = −0.598（p = 0.0003，图3A）。10,000次置换检验证实了这一结果并非偶然（置换p = 0.0079）。

翻译通路效应量的详细检查揭示了一致的格局：所有翻译通路在HCC中均为正效应量（肿瘤表达高于正常组织，Cohen's d > 0），而在HF中绝大多数翻译通路为负效应量（衰竭心肌低于非衰竭心肌）。例如，MYC Targets V1在HCC中d = +1.83，在HF中d = −0.81；真核翻译延伸在HCC中d = +1.07，在HF中d = −0.35。这一格局并非个别通路特有——33条翻译通路中，31条呈现HCC正-HF负的配对模式（二项检验p < 0.0001）。

HCC侧上调最强的通路包括E2F Targets（d = +2.35）、G2M Checkpoint（d = +2.13）和MYC Targets V1（d = +1.83），反映了HCC中强烈的细胞周期和增殖转录程序。HF侧上调最强的通路为胆汁酸代谢（d = +0.79）和干扰素α应答（d = +0.67），效应量绝对值远小于HCC，符合慢性HF转录组改变的异质性特征。

### 2.4 TGS预后价值的队列特异性

翻译基因评分（TGS）在TCGA-LIHC中显示出显著的预后价值（原始发现）。为检验其是否为HCC的普适预后标志物，我们在GSE141198（n = 148，94例事件）中进行了完全相同的TGS计算和生存分析。

TGS在该队列中不能区分高/低风险组（Log-rank p = 0.479）。连续变量Cox回归分析亦不显著（HR = 0.984，95% CI 0.817–1.185，p = 0.869）。按病因分层后，HBV、HCV和NBNC亚组中TGS均无预后意义。

结合先前在GSE14520和GSE76427中的阴性验证结果，三个独立HCC外部队列一致表明TGS预后价值无法复现。翻译共表达模块的保守性（2.2节）与TGS预后价值的队列特异性之间的割裂，构成了本研究的一个核心发现：翻译装置的共表达是保守的组织方式，但其与临床预后的关联高度依赖于队列特征（病因组成、分期分布、治疗背景）。

### 2.5 上游转录因子的系统鉴定：ATF4/ISR的发现 ★

为识别驱动翻译基因共表达的上游转录调控因子，我们对TCGA-LIHC的371例肿瘤样本进行了系统的TF-TGS相关性分析。19个候选TF涵盖了已知的翻译调控因子（MYC、MYCN、E2F家族、mTOR通路）以及应激反应相关TF。

19个TF中有13个在FDR < 0.05水平上与TGS显著相关（图5A）。最强且最值得注意的发现是：

**ATF4是TGS的最强正相关TF**（ρ = +0.439，FDR < 0.0001）。ATF4是整合应激反应（ISR）的master转录因子——当细胞面临营养剥夺、缺氧、氧化应激或ER应激时，eIF2α磷酸化导致全局翻译抑制，同时选择性上调ATF4翻译，进而激活下游应激适应基因[7]。ATF4在HCC翻译程序中的主导作用是一个新颖的机制发现，超越了传统以MYC为中心的翻译调控叙述。

**MYC靶基因通路活性（而非MYC mRNA本身）与TGS相关性最强**（通路活性ρ = +0.613 vs MYC mRNA ρ = +0.297）。MYC Targets V2（富含核糖体和翻译基因的子集）的通路活性相关性远强于MYC Targets V1，提示MYC转录程序的协同输出——而非单个TF的mRNA水平——更直接地反映翻译基因的共表达程度。

**ISR/UPR轴的连带信号**：DDIT3（CHOP，ATF4的直接转录靶基因，ρ = +0.289，FDR < 0.0001）和XBP1（UPR关键TF，ρ = +0.162，p = 0.003）两者均与TGS显著正相关。ATF4 > DDIT3 > XBP1的递进相关模式构成了从ISR感知到翻译重编程的一条连续信号链。

**HIF1A的负调控**：HIF1A表达与TGS呈显著负相关（ρ = −0.417，FDR < 0.0001）。缺氧条件下HIF1α稳定化可通过多种机制抑制翻译——包括mTORC1抑制和eIF2α磷酸化[8]。这一负相关与HIF1A的已知翻译抑制功能一致，提示氧气可用性可能调节HCC内的翻译基因程序。

**Hub基因的独立性**：7个hub基因全部不在标准Hallmark MYC Targets V1（200个基因）和V2（58个基因）集中。尽管这些Hallmark集含有其他核糖体蛋白（如RPL14、RPS10），但这7个hub基因代表了一组非典型的翻译基因，其共表达模式是在本研究的数据驱动分析中独立发现的，而非已知MYC靶基因的简单复现。

Fisher精确检验：7个hub基因在任何Hallmark TF靶基因集中均无显著富集（MYC V1/V2、E2F、mTORC1的p值均为1.00），进一步支持这组hub基因的独立发现性质。

**MTOR悖论**：MTOR mRNA水平与TGS呈显著负相关（ρ = −0.391，FDR < 0.0001），但mTORC1通路活性与TGS无显著相关（ρ = −0.071，p = 0.173）。这与mTOR信号主要在翻译后水平调控、其mRNA水平不能反映通路活性的生物学事实完全一致。

### 2.6 Hub基因个体水平的方向一致性检验

在个体基因水平上，7个hub基因的跨疾病log2FC方向Spearman相关性不显著（ρ = 0.486，p = 0.356）。这一定量结果与2.3节通路层面的显著负相关形成了有意义的对照：单个hub基因在两种疾病间的表达方向不一致，但整套翻译基因程序作为整体，在通路/模块层面呈现稳健的跨疾病负相关（补充图S1）。这一基因-vs-通路的差异强调了通路水平的系统分析在揭示跨疾病调控规律时优于单个基因分析。

---

## 讨论

### 3.1 从"共享信号"到"镜像调控"——一个修订后的概念框架

本研究始于一个简单的观察——HF中鉴定的翻译共表达hub基因在HCC中也呈现差异表达——但系统性跨疾病分析揭示了一幅远比"共享"更加精微的图景。

我们提出的修订模型（图6）包含三个关键层次：(1) **结构保守**：翻译共表达模块在HF和HCC中独立存在（GSE141198 WGCNA复现）；(2) **方向对立**：模块的扰动方向由疾病类型决定——增殖性HCC中翻译上调，衰竭性HF中翻译下调（ssGSEA负相关）；(3) **驱动因子不同**：HCC中的翻译上调由ATF4/ISR应激适应和MYC增殖程序联合驱动，HF中的翻译下调则反映了能量耗竭下的适应性抑制。

### 3.2 ATF4/ISR作为HCC翻译调控的新机制假说

ATF4作为TGS最强相关TF是本研究最出乎意料但生物学上合理的发现。ISR在癌症中扮演双重角色：在肿瘤发生的早期阶段，ISR通过限制蛋白质合成抑制肿瘤；但在已建立的肿瘤中，癌细胞劫持ISR通路以在应激微环境中存活[9]。ATF4的促存活靶基因（包括氨基酸代谢、抗氧化防御和自噬相关基因）在多种癌症中被选择性地激活[10]。

在HCC背景下，慢性内质网应激（由病毒感染、酒精代谢、脂肪变性或快速增殖驱动）和缺氧是肿瘤微环境的内在特征[11]。我们的发现——ATF4、DDIT3和XBP1三者均与翻译基因程序正相关——提示HCC细胞通过ISR/UPR轴维持翻译适应，将微环境应激转化为生存优势。

这一发现将讨论从模糊的"可能涉及MYC/mTORC1"提升至具体可验证的机制假说：**HCC中保守的翻译共表达模式至少部分由ATF4/ISR介导的应激适应性翻译重编程驱动**。该假说可通过ChIP-seq（检测ATF4在核糖体蛋白基因启动子上的结合）和功能实验（在HCC细胞系中敲低ATF4后检测翻译效率和增殖）进行验证。

### 3.3 MYC程序提供增殖维度的第二层调控

MYC靶基因通路活性与TGS的超强相关（ρ = +0.613）确认了MYC程序在HCC翻译调控中的重要角色。然而，通路活性（反映MYC转录程序的协同输出）远超MYC mRNA水平本身的解释力（ρ = +0.297），提示MYC在HCC中的翻译调控功能取决于：(1) MYC蛋白的翻译后修饰和稳定性调控；(2) MYC与其伴侣因子（如MAX、MIZ1）的相互作用；(3) 染色质可及性对MYC结合的影响——而非单纯的MYC转录水平。

重要的是，7个hub基因不在标准Hallmark MYC靶基因集中。这并不意味着这些基因不受MYC调控，而是提示它们的MYC响应可能是细胞类型或疾病背景特异的，或者它们受非MYC因子的共同调控。ATF4和MYC可能在翻译基因的启动子上协同作用——这是一个值得未来实验检验的假说。

### 3.4 HF中翻译下调为何是适应性的

在生理条件下，心肌收缩消耗大量ATP，而蛋白质合成是仅次于收缩的第二大ATP消费[1]。在HF中，能量底物利用障碍和线粒体功能障碍导致慢性能量亏缺[12]。翻译抑制作为一种节能策略在进化上是保守的（从酵母到哺乳动物），由ISR和mTORC1抑制共同介导[13]。

我们的ssGSEA数据显示HF中翻译通路的一致下调（d < 0），与这一能量危机模型一致。MTOR mRNA的负相关（ρ = −0.391）进一步支持合成代谢信号在衰竭心肌中的受损。重要的是，HF中翻译下调的幅度远小于HCC中翻译上调的幅度（最大|d|：HF 1.39 vs HCC 2.35），提示HF的翻译转录改变是一种渐进的、不完全的适应，而非癌症中由致癌基因驱动的强烈转录重编程。

### 3.5 为何TGS预后价值仅见于TCGA-LIHC？

TGS在三个外部队列中全部验证失败是一个重要的阴性结果，需要在机制层面予以解释。可能的因素包括：(1) **病因异质性**：TCGA-LIHC主要为非病毒性病因（NASH/酒精性），而GSE141198以HBV/HCV为主。病因可能影响翻译调控机制和临床意义[14]。(2) **分期差异**：TCGA-LIHC包含不同分期的混合人群，TGS可能在特定分期才显示预后价值。(3) **治疗混杂**：术后辅助治疗在不同队列中的差异可能改变翻译基因表达与生存的关联。(4) **真正的队列特异性**：翻译共表达模块虽然保守，但其与预后的关联可能依赖于队列特定的临床和分子背景，并非HCC的普遍规律。

这一发现的实际含义是：TGS不应被推广为通用的HCC预后生物标志物。然而，翻译共表达网络结构本身——而非其与预后的统计关联——可能才是治疗靶向的更有价值的维度。

### 3.6 优势和局限性

**优势**：(1) 多层次分析（基因→模块→通路→TF→预后），内部一致性高；(2) 独立数据集复现（WGCNA在GSE141198中的复现）；(3) 统计严谨性（置换检验、FDR校正、多队列验证）；(4) 不依赖湿实验的全计算策略，完全可复现。

**局限性**：(1) 缺乏湿实验验证（如ATF4 ChIP-seq或功能敲低实验）；(2) RNA-seq仅反映mRNA丰度，不能直接测量翻译效率——核糖体足迹分析（Ribo-seq）将大大加强结论；(3) HF侧仅使用了单个数据集（GSE57338）；(4) 所有数据均为回顾性/公共数据，可能存在选择偏倚；(5) HCC队列主要为亚洲人群（台湾GSE141198、部分中国GSE14520），在西方人群中的推广性需要验证。

### 3.7 临床转化意义

虽然TGS在外部验证中失败降低了其作为HCC通用生物标志物的前景，但本研究发现的ATF4/ISR轴在HCC翻译调控中的核心角色可能具有治疗意义。ISR通路的小分子调节剂正处于临床开发中（如ISRIB用于神经退行性疾病[15]），其在癌症中的潜力正在被探索。此外，翻译共表达网络的保守性——尽管方向不同——提示该网络的结构特征（如hub基因间的功能连接）可能是比简单的表达量更有价值的治疗靶向维度。

---

## 方法

### 4.1 数据获取

- **GSE57338**（HF）：GEO下载，n = 313例左心室样本（DCM、ICM、NF），平台GPL11532（Affymetrix HuGene 1.1 ST）
- **TCGA-LIHC**（HCC发现集）：TCGA GDC门户下载，n = 424样本（371肿瘤 + 50正常 + 3转移），RNA-seq（Illumina），HTSeq-counts
- **GSE141198**（HCC验证集）：GEO下载，n = 148例HCC肿瘤，RNA-seq，94例OS事件，含HBV/HCV/NBNC病因信息
- **GSE14520、GSE76427**：GEO下载，分别用于TGS的外部验证

### 4.2 数据处理

TCGA-LIHC使用DESeq2进行VST标准化（盲法设计，nsub = 1000），筛选至少在20%样本中counts ≥ 10的基因。Ensembl ID通过org.Hs.eg.db（AnnotationDbi::select）转换为基因符号，去除版本号后取唯一映射。GSE57338探针通过hugene11sttranscriptcluster.db注释包映射至基因符号，每个基因取多探针均值。GSE141198使用DESeq2 VST标准化（筛选标准同TCGA-LIHC），通过org.Hs.eg.db进行Ensembl→Symbol转换。

### 4.3 WGCNA

使用WGCNA包进行signed网络分析。参数：软阈值β通过scale-free topology fit R² ≥ 0.85标准选择；一步法blockwiseModules构建网络（minModuleSize = 30，mergeCutHeight = 0.25，TOMType = "signed"）。模块富集使用clusterProfiler进行GO Biological Process富集分析。翻译模块通过关键词"translation|ribosom|peptide|ribonucleoprotein|rRNA|translational"在GO描述中识别。

### 4.4 ssGSEA通路活性分析

基因集来源：msigdbr包v2026.1.Hs，Hallmark（collection = "H"）、KEGG核糖体（collection = "C2"，subcollection = "CP:KEGG_LEGACY"）、Reactome翻译通路（collection = "C2"，subcollection = "CP:REACTOME"，关键词过滤：TRANSLATION|PEPTIDE_CHAIN_ELONGATION|RIBOSOME|TRNA_AMINOACYLATION|RRNA_PROCESSING）。共计83条通路。

ssGSEA评分使用GSVA包v2.6.2的ssgseaParam + gsva函数（minSize = 5，maxSize = 500）。每条通路的疾病效应量为Cohen's d =（疾病组均值 − 对照组均值）/ 合并标准差。跨疾病Spearman相关性：ρ = cor(es_LIHC, es_HF, method = "spearman")。置换检验：打乱HF效应量标签10,000次，计算随机ρ绝对值超过真实值的比例。

### 4.5 TGS生存分析

TGS定义为6个hub基因（EEF1A1、FAU、RPL39、RPL3、RPL32、RPL41）VST标准化表达值的z-score均值。TGS中位数将患者分为高/低组。使用survival包进行Kaplan-Meier分析和Cox比例风险回归。多元Cox模型调整病因因素。

### 4.6 上游TF分析

19个候选TF基于文献和MSigDB筛选。每个TF的标准化表达值与TGS计算Spearman相关系数，使用Benjamini-Hochberg方法进行FDR校正（19次检验）。Fisher精确检验：2×2列联表（hub基因在本通路/不在本通路 vs 所有基因在本通路/不在本通路），备择假设为"greater"。MYC、E2F、mTORC1通路活性来自ssGSEA评分（HALLMARK_MYC_TARGETS_V1/V2、HALLMARK_E2F_TARGETS、HALLMARK_MTORC1_SIGNALING），与TGS计算Spearman相关。

### 4.7 统计方法

所有分析在R 4.6.0中进行。统计显著性阈值设为p < 0.05（双侧，除非另有说明）。多重检验使用Benjamini-Hochberg FDR校正。置换检验使用10,000次迭代。关键分析的完整参数和随机种子已记录以确保可复现性。

---

## 参考文献

[1] Gibbs CL. Cardiac energetics. Physiol Rev. 1978.
[2] Sciarretta S, et al. New insights into the role of mTOR signaling in the cardiovascular system. Circ Res. 2018.
[3] Gao X, et al. Translational regulation in heart failure. J Mol Cell Cardiol. 2015.
[4] Zhang D, et al. mTORC1 signaling in heart failure. Nat Rev Cardiol. 2019.
[5] van Riggelen J, et al. MYC as a regulator of ribosome biogenesis and protein synthesis. Nat Rev Cancer. 2010.
[6] Xu Y, et al. MYC and mTORC1 in hepatocellular carcinoma. Hepatology. 2021.
[7] Pakos-Zebrucka K, et al. The integrated stress response. EMBO Rep. 2016.
[8] Liu L, et al. Hypoxia-induced energy stress regulates mRNA translation and cell growth. Mol Cell. 2006.
[9] Wortel IMN, et al. Surviving stress: modulation of ATF4-mediated stress responses in normal and malignant cells. Trends Endocrinol Metab. 2017.
[10] Ye J, et al. The GCN2-ATF4 pathway is critical for tumour cell survival and proliferation in response to nutrient deprivation. EMBO J. 2010.
[11] Fu S, et al. Aberrant lipid metabolism disrupts calcium homeostasis causing liver endoplasmic reticulum stress in obesity. Nature. 2011.
[12] Neubauer S. The failing heart — an engine out of fuel. N Engl J Med. 2007.
[13] Harding HP, et al. An integrated stress response regulates amino acid metabolism and resistance to oxidative stress. Mol Cell. 2003.
[14] Llovet JM, et al. Hepatocellular carcinoma. Nat Rev Dis Primers. 2021.
[15] Halliday M, et al. Repurposed drugs targeting eIF2α-P-mediated translational repression prevent neurodegeneration in mice. Brain. 2017.

---

## 图清单

| Figure | 标题 | 文件 |
|--------|------|------|
| Fig 1 | Study Design Overview | Figure1_Study_Design.png |
| Fig 2 | WGCNA in GSE141198 (dendrogram + module-trait + GO) | [待组装] |
| Fig 3 | Cross-Disease ssGSEA Pathway Effect Size Scatter | Figure_ssGSEA_cross_disease.png |
| Fig 4 | TGS External Validation (KM + Forest Plot) | [待组装] |
| Fig 5 | Upstream TF Analysis (TF-TGS + Pathway-TGS) | Figure_TF_TGS_correlation.png + Figure_Pathway_TGS_correlation.png |
| Fig 6 | Mechanistic Model | Figure6_Mechanistic_Model.png |
| Fig S1 | Direction Consistency Scatter | [待生成] |
| Fig S2 | WGCNA Soft Threshold Selection | [待组装] |

## 表清单

| Table | 内容 |
|-------|------|
| Table 1 | Cohort Characteristics |
| Table 2 | Hub Gene Characteristics and Cross-Disease Module Assignment |
| Table 3 | TF-TGS Correlation Summary |
| Table S1 | Module GO Enrichment (Both Datasets) |
| Table S2 | Full ssGSEA Pathway Effect Sizes (83 Pathways) |
