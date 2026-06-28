# =============================================================================
# 脚本 01b：GSE141198数据解析（在01下载后立即运行）
# 根据实际pdata结构精确提取临床变量
# =============================================================================

library(dplyr)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ── 1. 加载已下载数据 ────────────────────────────────────────────────────────
gse <- readRDS(file.path(PROJ_DIR, "GSE141198_raw.rds"))
eset <- gse[[1]]

# ── 2. 检查表达矩阵 ──────────────────────────────────────────────────────────
expr <- exprs(eset)
message(sprintf("表达矩阵: %d genes × %d samples", nrow(expr), ncol(expr)))

# 检查基因ID格式
rn <- rownames(expr)
cat("\n前20个基因名:\n")
print(head(rn, 20))
cat("\n后20个基因名:\n")
print(tail(rn, 20))

# 检查是否有重复基因名
dup_genes <- sum(duplicated(rn))
message(sprintf("重复基因名: %d", dup_genes))

# 检查表达值范围（判断是否为log2转换后）
cat("\n表达值统计:\n")
print(summary(as.vector(expr[1:min(1000, nrow(expr)), ])))

# ── 3. 提取临床变量 ──────────────────────────────────────────────────────────
pdata <- pData(eset)

# 根据实际列名提取
clinical <- data.frame(
  geo_id = pdata$geo_accession,
  sample_name = pdata$title,
  row.names = pdata$geo_accession,
  stringsAsFactors = FALSE
)

# 解析characteristics_ch1系列
parse_char <- function(x) {
  if (is.null(x)) return(rep(NA, nrow(pdata)))
  # 有些列格式是 "label: value"
  val <- gsub("^[^:]+:\\s*", "", x)
  val <- trimws(val)
  return(val)
}

# characteristics_ch1 到 ch1.5
clinical$ctnnb1_status <- parse_char(pdata[["ctnnb1 status:ch1"]])

# 注意：列名可能略有不同，测试实际列名
cat("\ncharacteristics相关列名:\n")
chars_cols <- grep("characteristics|ch1", colnames(pdata), value = TRUE)
print(chars_cols)

# 打印各characteristics列的前3个值以确认内容
for (col in chars_cols) {
  cat(sprintf("\n  %s:", col))
  cat(sprintf(" [%s]", paste(head(pdata[[col]], 3), collapse=" | ")))
}

# 直接读取OS和EFS列
clinical$os_days <- as.numeric(pdata[["os days:ch1"]])
clinical$os_event <- as.numeric(pdata[["os event:ch1"]])
clinical$efs_days <- as.numeric(pdata[["efs days:ch1"]])
clinical$efs_event <- as.numeric(pdata[["efs event:ch1"]])

# 从characteristics列解析etiology
# 先确认etiology在哪个characteristics列
for (col in chars_cols) {
  vals <- pdata[[col]]
  if (any(grepl("HBV|HCV|alcohol|NBNC|NASH|etiology", vals, ignore.case = TRUE))) {
    cat(sprintf("\n✅ Etiology found in: %s\n", col))
    clinical$etiology <- parse_char(vals)
    break
  }
}

# ── 4. 生存数据整理 ──────────────────────────────────────────────────────────
# 清理：去除缺失OS数据的样本
clinical$os_time <- clinical$os_days
clinical$os_status <- clinical$os_event

# 报告
message(sprintf("\n═══ 临床数据摘要 ═══"))
message(sprintf("总样本数: %d", nrow(clinical)))
message(sprintf("OS可用: %d (events: %d)",
                sum(!is.na(clinical$os_time)),
                sum(clinical$os_status == 1, na.rm = TRUE)))
message(sprintf("EFS可用: %d (events: %d)",
                sum(!is.na(clinical$efs_days)),
                sum(clinical$efs_event == 1, na.rm = TRUE)))

if ("etiology" %in% colnames(clinical)) {
  message(sprintf("病因分布:"))
  print(table(clinical$etiology))
}
if ("ctnnb1_status" %in% colnames(clinical)) {
  message(sprintf("CTNNB1突变分布:"))
  print(table(clinical$ctnnb1_status))
}

# ── 5. 表达矩阵预处理 ────────────────────────────────────────────────────────
# 处理重复基因名（取均值）
if (dup_genes > 0) {
  message("\n⚠ 处理重复基因名...")
  expr_dedup <- matrix(NA, nrow = length(unique(rn)), ncol = ncol(expr))
  rownames(expr_dedup) <- unique(rn)
  colnames(expr_dedup) <- colnames(expr)
  for (g in unique(rn)) {
    idx <- which(rn == g)
    if (length(idx) == 1) {
      expr_dedup[g, ] <- expr[idx, ]
    } else {
      expr_dedup[g, ] <- colMeans(expr[idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  expr <- expr_dedup
  message(sprintf("去重后: %d genes × %d samples", nrow(expr), ncol(expr)))
}

# ── 6. 基因过滤（为WGCNA准备）──────────────────────────────────────────────
# 检查是否为log2转换后数据
is_log2 <- max(expr, na.rm = TRUE) < 30 && min(expr, na.rm = TRUE) >= 0
message(sprintf("数据是否为log2转换: %s", is_log2))

if (is_log2) {
  # log2(FPKM+1)格式：过滤低表达
  min_samples <- floor(0.2 * ncol(expr))
  keep <- rowSums(expr > 0.5) >= min_samples
} else {
  # Counts格式：过滤低counts
  min_samples <- floor(0.2 * ncol(expr))
  keep <- rowSums(expr >= 10) >= min_samples
}
message(sprintf("基因过滤: %d/%d 保留 (%.1f%%)", sum(keep), length(keep), 100*mean(keep)))
expr_filt <- expr[keep, ]

# ── 7. 保存处理后的数据 ──────────────────────────────────────────────────────
saveRDS(expr_filt, file.path(PROJ_DIR, "GSE141198_expr_clean.rds"))
saveRDS(clinical, file.path(PROJ_DIR, "GSE141198_clinical.rds"))

# 检查翻译hub基因是否在数据中
hub_genes <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41", "RPS28")
hub_in_data <- hub_genes %in% rownames(expr_filt)
cat("\n═══ Hub基因在GSE141198中 ═══\n")
for (i in seq_along(hub_genes)) {
  cat(sprintf("  %s: %s\n", hub_genes[i], ifelse(hub_in_data[i], "✅ 找到", "❌ 缺失")))
}

message(sprintf("\n%d/%d hub基因在GSE141198中", sum(hub_in_data), length(hub_genes)))

# ── 8. 输出供后续分析使用 ────────────────────────────────────────────────────
message("\n✅ GSE141198数据预处理完成")
message("输出文件:")
message("  GSE141198_expr_clean.rds — 清洗后表达矩阵")
message("  GSE141198_clinical.rds  — 结构化的临床数据")

# 重要的检查清单：
cat("\n─── 请在继续前确认 ───\n")
cat("1. 基因ID类型是Symbol还是Ensembl？\n")
cat("2. OS events数是否足够？（≥20 events用于Cox回归）\n")
cat("3. etiology + CTNNB1信息是否完整？\n")
cat("4. hub基因（尤其是RPL39和EEF1A1）是否在表达矩阵中？\n")
