# =============================================================================
# 脚本 01e：Ensembl ID → Gene Symbol 转换 + Hub基因定位
# GSE141198使用ENSGxxxx.x格式，需要转换为Symbol匹配hub基因
# =============================================================================

library(dplyr)
library(org.Hs.eg.db)
library(DESeq2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ── 1. 加载数据 ──────────────────────────────────────────────────────────────
rc_filt <- readRDS(file.path(PROJ_DIR, "GSE141198_counts_filt.rds"))
clinical <- readRDS(file.path(PROJ_DIR, "GSE141198_clinical.rds"))

# ── 2. 去掉Ensembl版本号 ────────────────────────────────────────────────────
# ENSG00000156508.15 → ENSG00000156508
ensembl_ids_full <- rownames(rc_filt)
ensembl_ids <- gsub("\\.\\d+$", "", ensembl_ids_full)
message(sprintf("Ensembl ID总数: %d", length(ensembl_ids)))
message(sprintf("唯一Ensembl ID: %d", length(unique(ensembl_ids))))

# 处理重复ID（多个版本对应同一基因）
# 策略：保留表达量最高的那个版本
expr_matrix <- as.matrix(rc_filt)
rownames(expr_matrix) <- ensembl_ids

# 每个Ensembl ID保留平均表达量最高的版本
dup_ids <- ensembl_ids[duplicated(ensembl_ids)]
if (length(dup_ids) > 0) {
  message(sprintf("重复Ensembl ID: %d个 (%d个基因有多个版本)",
                  length(dup_ids), length(unique(dup_ids))))
  # 对每个重复基因，计算均值，保留最高的
  unique_ids <- unique(ensembl_ids)
  expr_dedup <- matrix(NA, nrow = length(unique_ids), ncol = ncol(expr_matrix))
  rownames(expr_dedup) <- unique_ids
  colnames(expr_dedup) <- colnames(expr_matrix)
  for (id in unique_ids) {
    rows <- which(ensembl_ids == id)
    if (length(rows) == 1) {
      expr_dedup[id, ] <- expr_matrix[rows, ]
    } else {
      # 保留均值最高的版本
      row_means <- rowMeans(expr_matrix[rows, , drop = FALSE])
      best <- rows[which.max(row_means)]
      expr_dedup[id, ] <- expr_matrix[best, ]
    }
  }
  expr_matrix <- expr_dedup
  message(sprintf("去重后: %d 基因", nrow(expr_matrix)))
}

# ── 3. Ensembl → Symbol 转换 ─────────────────────────────────────────────────
message("\n═══ Ensembl → Symbol ═══")

# 获取所有Ensembl ID的Symbol映射
ensembl_only <- rownames(expr_matrix)
# 确保是纯Ensembl ID（不含版本号）
stopifnot(all(grepl("^ENSG", ensembl_only)))

# 使用org.Hs.eg.db做映射
library(AnnotationDbi)
symbol_map <- select(org.Hs.eg.db,
                    keys = ensembl_only,
                    columns = c("SYMBOL", "GENENAME"),
                    keytype = "ENSEMBL")

message(sprintf("映射结果: %d 行", nrow(symbol_map)))
message(sprintf("有Symbol的Ensembl ID: %d/%d (%.0f%%)",
                sum(!is.na(symbol_map$SYMBOL)), length(unique(symbol_map$ENSEMBL)),
                100 * mean(!is.na(symbol_map$SYMBOL))))

# 处理一对多映射（一个Ensembl→多个Symbol）——取第一个
symbol_map <- symbol_map %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# ── 4. 替换行名为Symbol ──────────────────────────────────────────────────────
# 仅保留有Symbol映射的基因
common_ensembl <- intersect(ensembl_only, symbol_map$ENSEMBL)
message(sprintf("可映射基因: %d/%d (%.0f%%)",
                length(common_ensembl), length(ensembl_only),
                100 * length(common_ensembl) / length(ensembl_only)))

expr_symbol <- expr_matrix[common_ensembl, ]
# 替换行名
idx <- match(common_ensembl, symbol_map$ENSEMBL)
rownames(expr_symbol) <- symbol_map$SYMBOL[idx]

# 如果还有重复Symbol，取均值
dup_sym <- duplicated(rownames(expr_symbol))
if (any(dup_sym)) {
  message(sprintf("处理%d个重复Symbol...", sum(dup_sym)))
  unique_sym <- unique(rownames(expr_symbol))
  expr_dedup2 <- matrix(NA, nrow = length(unique_sym), ncol = ncol(expr_symbol))
  rownames(expr_dedup2) <- unique_sym
  colnames(expr_dedup2) <- colnames(expr_symbol)
  for (s in unique_sym) {
    rows <- which(rownames(expr_symbol) == s)
    if (length(rows) == 1) {
      expr_dedup2[s, ] <- expr_symbol[rows, ]
    } else {
      expr_dedup2[s, ] <- colMeans(expr_symbol[rows, , drop = FALSE])
    }
  }
  expr_symbol <- expr_dedup2
}
message(sprintf("最终Symbol矩阵: %d genes × %d samples", nrow(expr_symbol), ncol(expr_symbol)))

# ── 5. Hub基因检查 ──────────────────────────────────────────────────────────
hub_genes <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41", "RPS28")
hub_in_data <- hub_genes %in% rownames(expr_symbol)

cat("\n═══ Hub基因在GSE141198中（Symbol转换后）═══\n")
for (i in seq_along(hub_genes)) {
  cat(sprintf("  %s: %s\n", hub_genes[i],
              ifelse(hub_in_data[i], "✅ 找到", "❌ 缺失")))
}
message(sprintf("\n%d/%d hub基因在GSE141198中", sum(hub_in_data), length(hub_genes)))

# 如果仍有缺失，打印这些基因的Ensembl ID
missing_genes <- hub_genes[!hub_in_data]
if (length(missing_genes) > 0) {
  message("\n缺失基因的Ensembl ID查询:")
  for (g in missing_genes) {
    res <- select(org.Hs.eg.db, keys = g, columns = "ENSEMBL",
                  keytype = "SYMBOL")
    cat(sprintf("  %s → %s\n", g, paste(res$ENSEMBL, collapse=", ")))
  }
}

# ── 6. DESeq2标准化（为WGCNA和下游分析准备）─────────────────────────────────
message("\n═══ DESeq2 VST标准化 ═══")

# 创建DESeq2对象
dds <- DESeqDataSetFromMatrix(
  countData = expr_symbol,
  colData = data.frame(row.names = colnames(expr_symbol), condition = 1),
  design = ~ 1
)

# 快速VST（blind设计，因无分组变量）
vsd <- vst(dds, blind = TRUE, nsub = min(1000, nrow(dds)))
expr_vst <- assay(vsd)

message(sprintf("VST标准化后: %d genes × %d samples", nrow(expr_vst), ncol(expr_vst)))

# ── 7. 保存 ──────────────────────────────────────────────────────────────────
saveRDS(expr_symbol, file.path(PROJ_DIR, "GSE141198_counts_symbol.rds"))
saveRDS(expr_vst, file.path(PROJ_DIR, "GSE141198_vst.rds"))
saveRDS(clinical, file.path(PROJ_DIR, "GSE141198_clinical.rds"))

# 为WGCNA选取最高变异基因
gene_vars <- apply(expr_vst, 1, var)
top5000 <- names(sort(gene_vars, decreasing = TRUE))[1:min(5000, length(gene_vars))]
expr_wgcna <- expr_vst[top5000, ]
# 确保hub基因在其中（即使方差不够大）
extra_genes <- setdiff(hub_genes[hub_in_data], top5000)
if (length(extra_genes) > 0) {
  expr_wgcna <- rbind(expr_wgcna, expr_vst[extra_genes, , drop = FALSE])
}
saveRDS(expr_wgcna, file.path(PROJ_DIR, "GSE141198_wgcna_input.rds"))

cat("\n══════════════════════════════════════════════════════\n")
cat("GSE141198数据准备完成！\n")
cat("──────────────────────────────────────────────────\n")
cat(sprintf("  样本: %d (含94个OS事件, 64%%)\n", ncol(expr_vst)))
cat(sprintf("  病因: HBV %d / HCV %d / NBNC %d\n",
            sum(clinical$etiology == "HBV"),
            sum(clinical$etiology == "HCV"),
            sum(clinical$etiology == "NBNC")))
cat(sprintf("  Hub基因: %d/7找到\n", sum(hub_in_data)))
cat(sprintf("  WGCNA输入: %d genes × %d samples\n",
            nrow(expr_wgcna), ncol(expr_wgcna)))
cat("══════════════════════════════════════════════════════\n")
