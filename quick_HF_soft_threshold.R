# quick_HF_soft_threshold.R
# 最小脚本：GSE57338 软阈值选择 → Figure S1

library(WGCNA)
library(GEOquery)
library(hugene11sttranscriptcluster.db)
library(AnnotationDbi)
enableWGCNAThreads(2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# 1. 加载本地文件
gz_file <- file.path(PROJ_DIR, "GSE57338_series_matrix.txt.gz")
gse <- getGEO(filename = gz_file, getGPL = FALSE)
eset <- gse
exprs_raw <- Biobase::exprs(eset)
cat(sprintf("%d probes × %d samples\n", nrow(exprs_raw), ncol(exprs_raw)))

# 2. 探针→基因
probe_ids <- as.character(rownames(exprs_raw))
sym_map <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
  keys = probe_ids, keytype = "PROBEID", columns = "SYMBOL")
sym_map <- sym_map[!is.na(sym_map$SYMBOL) & sym_map$SYMBOL != "", ]

keep <- rowSums(exprs_raw >= 4) >= ncol(exprs_raw) * 0.2
exprs_f <- exprs_raw[keep, ]
sym_f <- sym_map[sym_map$PROBEID %in% rownames(exprs_f), ]
exprs_f <- exprs_f[sym_f$PROBEID, , drop = FALSE]

uniq_syms <- unique(sym_f$SYMBOL)
exprs_gene <- t(sapply(uniq_syms, function(s) {
  probes <- sym_f$PROBEID[sym_f$SYMBOL == s]
  if (length(probes) == 1) exprs_f[probes, ] else colMeans(exprs_f[probes, , drop = FALSE])
}))
cat(sprintf("Gene matrix: %d × %d\n", nrow(exprs_gene), ncol(exprs_gene)))

# 3. 取高变异基因 + 软阈值
datExpr0 <- t(exprs_gene)
vars <- apply(datExpr0, 2, var, na.rm = TRUE)
datExpr0 <- datExpr0[, order(vars, decreasing = TRUE)[1:min(3000, ncol(datExpr0))]]

powers <- 1:20
sft <- pickSoftThreshold(datExpr0, powerVector = powers,
                         networkType = "signed", verbose = 2)

# 4. 输出图
png(file.path(PROJ_DIR, "figures", "Figure_S4B_SoftThreshold_GSE57338.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R²",
     main = "Scale Independence (GSE57338)", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = ifelse(sft$fitIndices$SFT.R.sq >= 0.85, "red", "black"))
abline(h = 0.85, col = "red", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean Connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers,
     col = ifelse(sft$fitIndices$SFT.R.sq >= 0.85, "red", "black"))
dev.off()

# 报告
idx <- which(sft$fitIndices$SFT.R.sq >= 0.85)
beta_sel <- if (length(idx) > 0) sft$fitIndices$Power[min(idx)] else sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
cat(sprintf("\nβ = %d (R²=%.3f)\n", beta_sel,
            sft$fitIndices$SFT.R.sq[which(sft$fitIndices$Power == beta_sel)]))
cat("→ figures/Figure_S4B_SoftThreshold_GSE57338.png saved\n")
