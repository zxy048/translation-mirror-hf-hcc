# =============================================================================
# 脚本 01f：修复VST错误 + 完成GSE141198最终预处理
# 7/7 hub基因全在数据中 ✅
# =============================================================================

library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# 从原始counts重新开始，更稳健的处理
rc_raw <- readRDS(file.path(PROJ_DIR, "GSE141198_counts_raw.rds"))
clinical <- readRDS(file.path(PROJ_DIR, "GSE141198_clinical.rds"))

message(sprintf("原始counts: %d genes × %d samples", nrow(rc_raw), ncol(rc_raw)))

# ── 1. 干净地转换Ensembl → Symbol ───────────────────────────────────────────
ensembl_full <- rownames(rc_raw)
ensembl_clean <- gsub("\\.\\d+$", "", ensembl_full)

# 使用AnnotationDbi精确映射
map_df <- select(org.Hs.eg.db,
                 keys = ensembl_clean,
                 columns = "SYMBOL",
                 keytype = "ENSEMBL")

# 只保留1:1映射
map_df <- map_df[!is.na(map_df$SYMBOL), ]
# 去除一个Ensembl→多个Symbol的情况
dup_ens <- duplicated(map_df$ENSEMBL)
map_df <- map_df[!dup_ens, ]
# 去除一个Symbol→多个Ensembl的情况
dup_sym <- duplicated(map_df$SYMBOL)
map_df <- map_df[!dup_sym, ]

message(sprintf("1:1映射: %d Ensembl → Symbol", nrow(map_df)))

# 找到可映射的基因
keep_rows <- ensembl_clean %in% map_df$ENSEMBL
message(sprintf("可映射基因数: %d/%d (%.0f%%)",
                sum(keep_rows), length(keep_rows), 100*mean(keep_rows)))

rc_mapped <- rc_raw[keep_rows, ]
ensembl_mapped <- ensembl_clean[keep_rows]

# 直接替换行名为Symbol（使用match确保顺序正确）
idx <- match(ensembl_mapped, map_df$ENSEMBL)
rownames(rc_mapped) <- map_df$SYMBOL[idx]

# 验证
hub_genes <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41", "RPS28")
hub_found <- hub_genes[hub_genes %in% rownames(rc_mapped)]
message(sprintf("Hub基因在映射后矩阵中: %d/7", length(hub_found)))

# ── 2. DESeq2过滤和VST ─────────────────────────────────────────────────────
# 使用DESeq2的推荐过滤
dds <- DESeqDataSetFromMatrix(
  countData = rc_mapped,
  colData = data.frame(
    row.names = colnames(rc_mapped),
    condition = factor(rep("A", ncol(rc_mapped)))
  ),
  design = ~ 1
)

# DESeq2推荐：保留至少在最小样本数中有>=10 reads的基因
smallestGroupSize <- ncol(rc_mapped)  # 没有分组，使用全部
keep <- rowSums(counts(dds) >= 10) >= floor(0.2 * ncol(rc_mapped))
dds <- dds[keep, ]
message(sprintf("DESeq2过滤后: %d genes", nrow(dds)))

# VST
vsd <- vst(dds, blind = TRUE, nsub = min(1000, nrow(dds)))
expr_vst <- assay(vsd)
message(sprintf("VST矩阵: %d genes × %d samples", nrow(expr_vst), ncol(expr_vst)))

# 再次确认hub基因（过滤后可能丢失了一些）
hub_final <- hub_genes[hub_genes %in% rownames(expr_vst)]
message(sprintf("Hub基因在VST矩阵中: %d/7: %s",
                length(hub_final), paste(hub_final, collapse=", ")))

# ── 3. 准备WGCNA输入 ────────────────────────────────────────────────────────
# 选5000个最高变异基因
rv <- rowVars(expr_vst)
names(rv) <- rownames(expr_vst)
top_genes <- names(sort(rv, decreasing = TRUE))[1:min(5000, length(rv))]

# 强制加入所有hub基因
top_genes <- unique(c(top_genes, hub_final))
top_genes <- intersect(top_genes, rownames(expr_vst))

expr_wgcna <- expr_vst[top_genes, ]
message(sprintf("WGCNA输入: %d genes × %d samples", nrow(expr_wgcna), ncol(expr_wgcna)))

# ── 4. 保存 ──────────────────────────────────────────────────────────────────
saveRDS(expr_vst, file.path(PROJ_DIR, "GSE141198_vst.rds"))
saveRDS(expr_wgcna, file.path(PROJ_DIR, "GSE141198_wgcna_input.rds"))
saveRDS(counts(dds), file.path(PROJ_DIR, "GSE141198_counts_clean.rds"))
saveRDS(clinical, file.path(PROJ_DIR, "GSE141198_clinical.rds"))

# ── 5. 最终检查清单 ─────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════\n")
cat("  ✅ GSE141198 数据准备完成！\n")
cat("────────────────────────────────────\n")
cat(sprintf("  VST矩阵: %d genes × %d samples\n", nrow(expr_vst), ncol(expr_vst)))
cat(sprintf("  WGCNA输入: %d genes × %d samples\n", nrow(expr_wgcna), ncol(expr_wgcna)))
cat(sprintf("  Hub基因: %d/7 ✅\n", length(hub_final)))
cat(sprintf("  缺失Hub基因: %s\n",
            paste(setdiff(hub_genes, hub_final), collapse=", ")))
cat(sprintf("  OS events: %d/148\n",
            sum(clinical$os_status == 1, na.rm = TRUE)))
cat("═══════════════════════════════════════\n")
