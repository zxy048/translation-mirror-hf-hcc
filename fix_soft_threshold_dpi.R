# fix_soft_threshold_dpi.R
# 修复: Figure_S1_SoftThreshold_GSE141198.png 从 150 DPI → 300 DPI
# 同时修复: 阈值线从 0.8 → 0.85 (匹配稿件)
# 运行时间: ~30秒

library(WGCNA)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# 1. 加载 WGCNA 输入数据
expr_wgcna <- readRDS(file.path(PROJ_DIR, "GSE141198_wgcna_input.rds"))
message(sprintf("WGCNA输入: %d genes × %d samples", nrow(expr_wgcna), ncol(expr_wgcna)))

datExpr0 <- t(expr_wgcna)

# 2. 移除异常样本/基因 (与04b脚本完全一致)
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  message("已移除异常样本/基因")
}

# 3. 软阈值计算
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers,
                         networkType = "signed", verbose = 5)

# 4. 输出 300 DPI 图
png(file.path(PROJ_DIR, "figures", "Figure_S4A_SoftThreshold_GSE141198.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))

# 左: Scale independence
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R²",
     main = "Scale Independence (GSE141198)", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = ifelse(sft$fitIndices$SFT.R.sq >= 0.85, "red", "black"))
abline(h = 0.85, col = "red", lty = 2)

# 右: Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean Connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers,
     col = ifelse(sft$fitIndices$SFT.R.sq >= 0.85, "red", "black"))
dev.off()

# 5. 报告
beta_candidates <- which(sft$fitIndices$SFT.R.sq >= 0.85)
if (length(beta_candidates) > 0) {
  beta_sel <- sft$fitIndices$Power[min(beta_candidates)]
} else {
  beta_sel <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
}
message(sprintf("β = %d (R² = %.3f)", beta_sel,
                sft$fitIndices$SFT.R.sq[which(sft$fitIndices$Power == beta_sel)]))

# 6. 验证分辨率
message("→ figures/Figure_S4A_SoftThreshold_GSE141198.png 已保存 (300 DPI)")
message("完成。请在 R 中运行: source('fix_soft_threshold_dpi.R')")
