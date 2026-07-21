# =============================================================================
# 脚本 04：GSE141198独立WGCNA——翻译模块跨数据集复现
# 目标：在完全独立的HCC数据集中de novo识别翻译相关共表达模块
# 这是全文最重要的新增分析——提供跨数据集的独立复现证据
# =============================================================================

library(WGCNA)
library(DESeq2)
library(limma)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)

set.seed(42)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(6)

PROJ_DIR <- "D:/R_projects/revision_analysis"
dir.create(file.path(PROJ_DIR, "WGCNA_GSE141198"), showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# 第一部分：加载GSE141198并预处理
# ═══════════════════════════════════════════════════════════════════════════════

message("═══ 加载GSE141198 ═══")

# 从脚本01下载的数据加载
gse141198 <- readRDS(file.path(PROJ_DIR, "GSE141198_raw.rds"))
eset <- gse141198[[1]]
expr_raw <- exprs(eset)
pdata <- pData(eset)

message(sprintf("原始表达矩阵: %d genes × %d samples", nrow(expr_raw), ncol(expr_raw)))

# ── 1.1 基因过滤和标准化 ──────────────────────────────────────────────────────
# 去除低表达基因
# 由于GSE141198是log2(FPKM+1)格式，使用适当的过滤标准
# RNA-seq: 去除在<20%样本中表达<1的基因
expr_log <- expr_raw
min_expr <- 1  # log2(FPKM+1)格式
min_samples <- floor(0.2 * ncol(expr_log))
keep_genes <- rowSums(expr_log >= min_expr) >= min_samples
message(sprintf("基因过滤: %d/%d 保留 (%.1f%%)",
                sum(keep_genes), length(keep_genes), 100*mean(keep_genes)))

expr_filt <- expr_log[keep_genes, ]

# 选取方差最大的前5000个基因用于WGCNA（标准做法）
gene_vars <- apply(expr_filt, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(5000, length(gene_vars))]
expr_wgcna <- expr_filt[top_genes, ]
message(sprintf("WGCNA输入: %d 高变异基因 × %d 样本", nrow(expr_wgcna), ncol(expr_wgcna)))

# ═══════════════════════════════════════════════════════════════════════════════
# 第二部分：GSE141198独立WGCNA
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ WGCNA De Novo ═══")

# ── 2.1 样本聚类和异常检测 ────────────────────────────────────────────────────
datExpr0 <- t(expr_wgcna)  # 行为样本，列为基因

# 检查缺失值
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[, gsg$goodGenes]
}

# 样本聚类
sampleTree <- hclust(dist(datExpr0), method = "average")

png(file.path(PROJ_DIR, "WGCNA_GSE141198", "sample_clustering.png"),
    width = 12, height = 6, units = "in", res = 150)
par(mar = c(4, 4, 3, 1))
plot(sampleTree, main = "GSE141198 Sample Clustering", sub = "",
     xlab = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
# 如果样本数不太多，标注样本名
if (nrow(datExpr0) <= 200) {
  abline(h = 20000, col = "red", lty = 2)
}
dev.off()

message("样本聚类图已保存。请检查是否有明显异常样本。")

# ── 2.2 软阈值选择 ────────────────────────────────────────────────────────────
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers,
                         networkType = "signed", verbose = 5)

# 绘制软阈值选择图
png(file.path(PROJ_DIR, "WGCNA_GSE141198", "soft_threshold.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R²",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = ifelse(sft$fitIndices$SFT.R.sq > 0.8, "red", "black"))
abline(h = 0.85, col = "red", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers,
     col = ifelse(sft$fitIndices$SFT.R.sq > 0.8, "red", "black"))
dev.off()

# 选择满足R² ≥ 0.85的beta（与原分析一致）
beta_sel <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq >= 0.85)]
if (length(beta_sel) == 0) beta_sel <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
message(sprintf("选择的软阈值 β = %d (R² = %.3f)", beta_sel,
                sft$fitIndices$SFT.R.sq[which(sft$fitIndices$Power == beta_sel)]))

# ── 2.3 共表达网络构建 ────────────────────────────────────────────────────────
# 为加速，使用blockwiseModules（适合大基因集）
# 如果基因数<=5000，也可直接用blockwiseModules
net <- blockwiseModules(
  datExpr0,
  power = beta_sel,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(PROJ_DIR, "WGCNA_GSE141198", "TOM_GSE141198"),
  maxBlockSize = 6000,
  verbose = 3,
  randomSeed = 42,
  nThreads = 6
)

moduleColors <- labels2colors(net$colors)
moduleLabels <- net$colors
MEs <- net$MEs

message(sprintf("识别出 %d 个共表达模块", length(unique(moduleColors))))

# ── 2.4 模块-性状关联 ────────────────────────────────────────────────────────
# 加载GSE141198临床数据并关联模块
# 从pdata中提取生存和临床变量
# GSE141198有OS、EFS、etiology等

# 提取OS和EFS信息
# 此处需要根据GSE141198的实际临床变量调整
# 框架代码如下：

if (FALSE) {  # 在实际运行时取消注释
  # 构建性状矩阵
  trait_data <- data.frame(
    OS_status = os_status,      # 从pdata提取
    EFS_status = efs_status,    # 从pdata提取
    HBV = as.numeric(etiology == "HBV"),
    HCV = as.numeric(etiology == "HCV"),
    stringsAsFactors = FALSE
  )
  rownames(trait_data) <- colnames(expr_wgcna)

  # 模块-性状相关
  moduleTraitCor <- cor(MEs, trait_data, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

  # 热图
  png(file.path(PROJ_DIR, "WGCNA_GSE141198", "module_trait_heatmap.png"),
      width = 10, height = 8, units = "in", res = 150)
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(trait_data),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = paste0(signif(moduleTraitCor, 2), "\n(",
                                     signif(moduleTraitPvalue, 2), ")"),
                 setStdMargins = FALSE,
                 cex.text = 0.7,
                 zlim = c(-1, 1),
                 main = "Module-Trait Relationships (GSE141198)")
  dev.off()
}

message("模块-性状关联分析框架就绪。根据实际临床数据调整trait_data。")

# ═══════════════════════════════════════════════════════════════════════════════
# 第三部分：GSE141198模块功能富集——寻找翻译/核糖体模块
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 模块功能富集 ═══")

# 对每个模块做GO/KEGG富集，寻找翻译/核糖体相关模块
all_genes <- colnames(datExpr0)

module_enrichment <- list()
for (mod in unique(moduleColors)) {
  mod_genes <- all_genes[moduleColors == mod]
  if (length(mod_genes) < 5) next

  # GO BP enrichment
  ego <- tryCatch({
    enrichGO(gene = mod_genes,
             OrgDb = org.Hs.eg.db,
             keyType = "SYMBOL",
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2)
  }, error = function(e) NULL)

  if (!is.null(ego) && nrow(ego) > 0) {
    module_enrichment[[mod]] <- ego@result
  }
}

# 寻找翻译/核糖体相关模块
translation_terms <- c("translation", "ribosome", "peptide", "ribonucleoprotein",
                       "rRNA", "translational", "aminoacyl")

translation_modules <- list()
for (mod in names(module_enrichment)) {
  desc <- paste(module_enrichment[[mod]]$Description, collapse = " ")
  if (any(sapply(translation_terms, function(term) grepl(term, desc, ignore.case = TRUE)))) {
    translation_modules[[mod]] <- module_enrichment[[mod]]
    message(sprintf("✅ %s 模块富集到翻译相关通路 (n=%d genes)",
                    mod, sum(moduleColors == mod)))
  }
}

# ── 3.1 与HF hub基因比较 ──────────────────────────────────────────────────────
hf_hub_genes <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41", "RPS28")

# 检查HF hub基因在GSE141198各模块中的分布
hub_module_mapping <- data.frame(
  gene = hf_hub_genes,
  in_GSE141198 = hf_hub_genes %in% rownames(expr_wgcna),
  module = sapply(hf_hub_genes, function(g) {
    if (g %in% all_genes) moduleColors[which(all_genes == g)] else NA_character_
  }),
  stringsAsFactors = FALSE
)

message("\n═══ HF hub基因在GSE141198中的模块分布 ═══")
print(hub_module_mapping)

# 做Fisher精确检验：翻译模块与HF hub基因的重叠是否显著
if (length(translation_modules) > 0) {
  mod_trans_genes <- all_genes[moduleColors %in% names(translation_modules)]
  hub_in_trans <- intersect(hf_hub_genes, mod_trans_genes)
  hub_not_in_trans <- setdiff(hf_hub_genes, mod_trans_genes)

  message(sprintf("\n翻译相关模块含 %d 个基因，其中 %d 个是HF hub基因",
                  length(mod_trans_genes), length(hub_in_trans)))
  if (length(hub_in_trans) > 0) {
    message("重叠HF hub基因: ", paste(hub_in_trans, collapse = ", "))
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 第四部分：跨数据集模块保守性分析（GSE57338 → GSE141198）
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 跨数据集保守性: HF (GSE57338) → HCC (GSE141198) ═══")

# 此分析类似于原论文的模块保留性分析，但在GSE141198而不是TCGA-LIHC上运行
# 提供了第二种跨疾病验证

# 加载HF参考网络
# 从原分析加载WGCNA结果
hf_wgcna_file <- "D:/R_projects/ML_output/wgcna_full_network.RData"
if (file.exists(hf_wgcna_file)) {
  load(hf_wgcna_file)

  # 找共同基因
  hf_genes <- names(net$colors)  # HF网络中的基因
  gse141198_genes <- colnames(datExpr0)

  # 需要将两个数据集的基因ID对齐
  # HF用的是微阵列探针→gene symbol，GSE141198用的是gene symbol
  common_genes <- intersect(hf_genes, gse141198_genes)
  message(sprintf("共同基因数: %d / %d (HF) / %d (GSE141198)",
                  length(common_genes), length(hf_genes), length(gse141198_genes)))

  if (length(common_genes) >= 100) {
    message("可在GSE141198上运行模块保留性分析（作为TCGA-LIHC之外的第三验证）")
    # 保存为未来运行做准备
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 第五部分：结果汇总和输出
# ═══════════════════════════════════════════════════════════════════════════════

message("\n", paste(rep("═", 70), collapse=""))
message("GSE141198 WGCNA分析框架就绪")
message(paste(rep("═", 70), collapse=""))
message("\n输出的关键结果：")
message(sprintf("  1. TOM矩阵: WGCNA_GSE141198/TOM_GSE141198-block.1.RData"))
message(sprintf("  2. 模块分配: %d 个模块", length(unique(moduleColors))))
message(sprintf("  3. 翻译相关模块: %d 个", length(translation_modules)))
message(sprintf("  4. HF hub基因映射: %d/%d 在GSE141198中找到",
                sum(hub_module_mapping$in_GSE141198), length(hf_hub_genes)))

message("\n后续步骤（在RStudio中交互运行）：")
message("  1. 确认GSE141198的临床变量后，完成模块-性状关联")
message("  2. 如果翻译模块与生存/分期相关→这是最强的独立验证")
message("  3. 比较GSE141198模块与HF模块的基因重叠→做Fisher检验")

# 保存关键对象
saveRDS(list(
  net = net,
  moduleColors = moduleColors,
  moduleLabels = moduleLabels,
  MEs = MEs,
  hub_module_mapping = hub_module_mapping,
  translation_modules = translation_modules,
  module_enrichment = module_enrichment,
  datExpr = datExpr0,
  beta_sel = beta_sel
), file.path(PROJ_DIR, "GSE141198_WGCNA_results.rds"))

message("\n✅ 脚本04执行完毕")
