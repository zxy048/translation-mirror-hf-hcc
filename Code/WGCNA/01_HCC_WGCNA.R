# =============================================================================
# 脚本 04b：GSE141198独立WGCNA——翻译模块跨数据集复现
# 输入：5003 genes × 148 samples VST标准化
# =============================================================================

library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

set.seed(42)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(6)

PROJ_DIR <- "D:/R_projects/revision_analysis"
wgcna_dir <- file.path(PROJ_DIR, "WGCNA_GSE141198")
dir.create(wgcna_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. 加载数据 ──────────────────────────────────────────────────────────────
expr_wgcna <- readRDS(file.path(PROJ_DIR, "GSE141198_wgcna_input.rds"))
clinical <- readRDS(file.path(PROJ_DIR, "GSE141198_clinical.rds"))

message(sprintf("WGCNA输入: %d genes × %d samples", nrow(expr_wgcna), ncol(expr_wgcna)))

datExpr0 <- t(expr_wgcna)
message(sprintf("datExpr: %d samples × %d genes", nrow(datExpr0), ncol(datExpr0)))

# ── 2. 样本聚类和异常检测 ──────────────────────────────────────────────────
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  message("已移除异常样本/基因")
}

# 样本聚类树
sampleTree <- hclust(dist(datExpr0), method = "average")
png(file.path(wgcna_dir, "sample_clustering.png"),
    width = 14, height = 7, units = "in", res = 150)
par(mar = c(4, 5, 3, 1))
plot(sampleTree, main = "GSE141198 Sample Clustering (HCC, n=148)",
     sub = "", xlab = "", cex = 0.5)
abline(h = 35000, col = "red", lty = 2, lwd = 1.5)  # 异常样本阈值
dev.off()
message("样本聚类图已保存")

# ── 3. 软阈值选择 ──────────────────────────────────────────────────────────
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers,
                         networkType = "signed", verbose = 5)

png(file.path(wgcna_dir, "soft_threshold.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R²",
     main = "Scale Independence (GSE141198)", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = ifelse(sft$fitIndices$SFT.R.sq > 0.8, "red", "black"))
abline(h = 0.85, col = "red", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean Connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers,
     col = ifelse(sft$fitIndices$SFT.R.sq > 0.8, "red", "black"))
dev.off()

# 选择beta
beta_candidates <- which(sft$fitIndices$SFT.R.sq >= 0.85)
if (length(beta_candidates) > 0) {
  beta_sel <- sft$fitIndices$Power[min(beta_candidates)]
} else {
  beta_sel <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
}
message(sprintf("选择 β = %d (R² = %.3f, mean connectivity = %.1f)",
                beta_sel,
                sft$fitIndices$SFT.R.sq[which(sft$fitIndices$Power == beta_sel)],
                sft$fitIndices$mean.k.[which(sft$fitIndices$Power == beta_sel)]))

# ── 4. 一步法网络构建 ──────────────────────────────────────────────────────
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
  saveTOMFileBase = file.path(wgcna_dir, "TOM_GSE141198"),
  maxBlockSize = 6000,
  verbose = 3,
  randomSeed = 42,
  nThreads = 6
)

moduleColors <- labels2colors(net$colors)
moduleLabels <- net$colors
MEs <- net$MEs

n_mods <- length(unique(moduleColors))
message(sprintf("识别出 %d 个共表达模块", n_mods))
cat("模块大小分布:\n")
tbl <- sort(table(moduleColors), decreasing = TRUE)
print(tbl[1:min(10, length(tbl))])

# ── 5. 模块-性状关联 ───────────────────────────────────────────────────────
# 匹配临床数据
sample_names <- rownames(datExpr0)
clinical$os_status_num <- as.numeric(clinical$os_status)
clinical$os_time_num <- clinical$os_time

# 构建性状矩阵
trait_data <- data.frame(
  row.names = sample_names,
  OS_status = clinical[sample_names, "os_status_num"],
  OS_time = clinical[sample_names, "os_time_num"] / 30.44,
  stringsAsFactors = FALSE
)

# 添加病因（哑变量）
if ("etiology" %in% colnames(clinical)) {
  for (etio in c("HBV", "HCV", "NBNC")) {
    trait_data[[paste0("Etiology_", etio)]] <- as.numeric(
      clinical[sample_names, "etiology"] == etio)
  }
}

# 添加CTNNB1
if ("ctnnb1" %in% colnames(clinical)) {
  trait_data$CTNNB1_mut <- as.numeric(grepl("mut|MT",
                                     clinical[sample_names, "ctnnb1"],
                                     ignore.case = TRUE))
}

# 模块-性状相关
nSamples <- nrow(datExpr0)
MEs0 <- orderMEs(MEs)
moduleTraitCor <- cor(MEs0, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# 热图
png(file.path(wgcna_dir, "module_trait_heatmap.png"),
    width = 12, height = 9, units = "in", res = 150)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(trait_data),
               yLabels = names(MEs0),
               ySymbols = names(MEs0),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = signif(moduleTraitCor, 2),
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships (GSE141198, n=148)")
dev.off()

# 找出与OS最显著相关的模块
os_cor <- moduleTraitCor[, "OS_status"]
os_p <- moduleTraitPvalue[, "OS_status"]
os_sig_modules <- names(os_p)[os_p < 0.05]
message(sprintf("\n与OS状态显著相关的模块: %d", length(os_sig_modules)))
for (m in os_sig_modules) {
  message(sprintf("  %s: r=%.3f, p=%.4f (n_genes=%d)",
                  m, os_cor[m], os_p[m], sum(moduleColors == m)))
}

# ── 6. 模块功能富集——寻找翻译/核糖体模块 ──────────────────────────────────
message("\n═══ 模块GO/KEGG富集 ═══")

all_genes <- colnames(datExpr0)
translation_terms <- "translation|ribosom|peptide|ribonucleoprotein|rRNA|translational"

translation_modules <- list()

for (mod in unique(moduleColors)) {
  mod_genes <- all_genes[moduleColors == mod]
  if (length(mod_genes) < 5) next

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
    desc_text <- paste(ego@result$Description, collapse = " ")
    if (grepl(translation_terms, desc_text, ignore.case = TRUE)) {
      translation_modules[[mod]] <- ego@result
      # 找最重要的翻译相关term
      trans_rows <- grep(translation_terms, ego@result$Description,
                         ignore.case = TRUE)
      top_term <- ego@result[trans_rows[1], ]
      message(sprintf("✅ %s (n=%d genes): %s (p=%.1e)",
                      mod, sum(moduleColors == mod),
                      top_term$Description, top_term$p.adjust))
    }
  }
}

message(sprintf("\n翻译/核糖体相关模块: %d/%d", length(translation_modules), n_mods))

# ── 7. Hub基因在GSE141198模块中的分布 ──────────────────────────────────────
hub_genes <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41", "RPS28")
hub_module_map <- data.frame(
  Gene = hub_genes,
  GSE141198_Module = sapply(hub_genes, function(g) {
    if (g %in% all_genes) as.character(moduleColors[which(all_genes == g)]) else "NOT_FOUND"
  }),
  stringsAsFactors = FALSE
)

cat("\n═══ Hub基因 GSE141198模块分布 ═══\n")
print(hub_module_map)

# 检查hub基因是否富集在同一个模块中
hub_mods <- hub_module_map$GSE141198_Module[hub_module_map$GSE141198_Module != "NOT_FOUND"]
hub_mod_table <- table(hub_mods)
cat("\nHub基因模块共定位:\n")
print(hub_mod_table)

# Fisher检验：hub基因在翻译模块中的富集
if (length(translation_modules) > 0) {
  trans_mod_genes <- all_genes[moduleColors %in% names(translation_modules)]
  hub_in_trans <- sum(hub_module_map$GSE141198_Module %in% names(translation_modules))
  hub_total <- sum(hub_module_map$GSE141198_Module != "NOT_FOUND")

  if (hub_in_trans >= 2) {
    a <- hub_in_trans; b <- hub_total - hub_in_trans
    c <- length(trans_mod_genes) - hub_in_trans
    d <- length(all_genes) - a - b - c
    fisher_p <- fisher.test(matrix(c(a, b, c, d), nrow = 2),
                            alternative = "greater")$p.value
    message(sprintf("\nFisher精确检验: %d/%d hub基因在翻译模块中, p=%.4f",
                    hub_in_trans, hub_total, fisher_p))
  }
}

# ── 8. 保存结果 ─────────────────────────────────────────────────────────────
wgcna_result <- list(
  net = net,
  moduleColors = moduleColors,
  moduleLabels = moduleLabels,
  MEs = MEs,
  beta_sel = beta_sel,
  n_modules = n_mods,
  translation_modules = translation_modules,
  hub_module_map = hub_module_map,
  os_sig_modules = os_sig_modules,
  moduleTraitCor = moduleTraitCor,
  moduleTraitPvalue = moduleTraitPvalue
)

saveRDS(wgcna_result, file.path(PROJ_DIR, "GSE141198_WGCNA_result.rds"))

# ── 9. 最终判定 ─────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════\n")
cat("  GSE141198 独立WGCNA 结果摘要\n")
cat("─────────────────────────────────\n")
cat(sprintf("  模块总数: %d\n", n_mods))
cat(sprintf("  翻译模块数: %d\n", length(translation_modules)))
cat(sprintf("  Hub基因共定位模块: %s\n",
            paste(names(which.max(hub_mod_table)),
                  collapse=", ")))
cat(sprintf("  与OS显著相关模块: %d\n", length(os_sig_modules)))

if (length(translation_modules) > 0 && length(os_sig_modules) > 0) {
  # 检查翻译模块是否与OS相关
  trans_os <- intersect(names(translation_modules), os_sig_modules)
  if (length(trans_os) > 0) {
    cat(sprintf("  ✅ 翻译模块 %s 同时与OS显著相关！\n",
                paste(trans_os, collapse=", ")))
    cat("  → 这构成独立复现：HF翻译基因集在GSE141198中形成\n")
    cat("    翻译相关模块，且该模块具有预后价值\n")
  }
}

if (length(translation_modules) == 0) {
  cat("  ⚠ 未独立识别出翻译相关模块\n")
  cat("  → 叙事改为：hub基因方向一致性在基因水平存在，\n")
  cat("    但独立HCC数据集不能重现核糖体共表达模块结构\n")
}

cat("═══════════════════════════════════════════════\n")
