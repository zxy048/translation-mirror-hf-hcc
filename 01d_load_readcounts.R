# =============================================================================
# 脚本 01d：从Supplementary文件加载GSE141198 read counts
# 结果：94 events / 148 samples — 统计检验力充足！
# =============================================================================

library(dplyr)
library(DESeq2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ── 1. 复制readcount文件到工作目录 ───────────────────────────────────────────
tmp_file <- file.path(tempdir(), "GSE141198", "GSE141198_TLCN.subset1.readcount.txt.gz")
if (!file.exists(tmp_file)) {
  # 从GEO重新下载
  library(GEOquery)
  tmp_dir <- tempdir()
  getGEOSuppFiles("GSE141198", baseDir = tmp_dir)
  tmp_file <- file.path(tmp_dir, "GSE141198", "GSE141198_TLCN.subset1.readcount.txt.gz")
}

message("读取文件: ", tmp_file)
message("文件大小: ", round(file.info(tmp_file)$size / 1024^2, 1), " MB")

# ── 2. 读取readcount矩阵 ────────────────────────────────────────────────────
rc <- read.table(gzfile(tmp_file), header = TRUE, row.names = 1,
                 check.names = FALSE, stringsAsFactors = FALSE)

message(sprintf("Read count矩阵: %d genes × %d samples", nrow(rc), ncol(rc)))

# 检查列名格式
cat("\n前5个列名:\n")
print(head(colnames(rc), 5))

# 检查基因ID类型（前20个基因名）
cat("\n前20个基因名:\n")
print(head(rownames(rc), 20))

# 检查表达值范围（raw counts还是标准化后）
cat("\n表达值统计:\n")
print(summary(as.vector(as.matrix(rc[1:min(2000, nrow(rc)), 1:min(5, ncol(rc))]))))

# ── 3. 匹配列名到GEO sample ID ────────────────────────────────────────────
# GEO的supplementary文件可能用sample title而非GSM ID
gse <- readRDS(file.path(PROJ_DIR, "GSE141198_raw.rds"))
eset <- gse[[1]]
pdata <- pData(eset)

# 映射：sample title → GSM ID
sample_map <- data.frame(
  title = pdata$title,
  geo_id = pdata$geo_accession,
  stringsAsFactors = FALSE
)

# 检查readcount列名匹配
cat("\nReadcount列名 vs GEO title:\n")
cat("Readcount列名样本: ", head(colnames(rc), 3), "\n")
cat("GEO title样本:     ", head(sample_map$title, 3), "\n")

# 尝试匹配
matches <- match(colnames(rc), sample_map$title)
message(sprintf("列名匹配率: %d/%d (%.0f%%)",
                sum(!is.na(matches)), length(matches),
                100 * mean(!is.na(matches))))

if (sum(!is.na(matches)) < ncol(rc) * 0.5) {
  # 可能用GSM ID做列名，检查
  matches_gsm <- match(colnames(rc), sample_map$geo_id)
  message(sprintf("GSM ID匹配率: %d/%d",
                  sum(!is.na(matches_gsm)), ncol(rc)))
  if (sum(!is.na(matches_gsm)) > ncol(rc) * 0.5) {
    matches <- matches_gsm
  }
}

# 重命名列为GSM ID
if (sum(!is.na(matches)) >= ncol(rc) * 0.5) {
  colnames(rc)[!is.na(matches)] <- sample_map$geo_id[matches[!is.na(matches)]]
  message("✅ 列名已转换为GSM ID")
}

# ── 4. 创建临床数据 ──────────────────────────────────────────────────────────
clinical <- data.frame(
  geo_id = pdata$geo_accession,
  os_time = as.numeric(pdata[["os days:ch1"]]),
  os_status = as.numeric(pdata[["os event:ch1"]]),
  efs_time = as.numeric(pdata[["efs days:ch1"]]),
  efs_status = as.numeric(pdata[["efs event:ch1"]]),
  ctnnb1 = pdata[["ctnnb1 status:ch1"]],
  row.names = pdata$geo_accession,
  stringsAsFactors = FALSE
)

# 解析characteristics列获取etiology
chars_cols <- grep("characteristics", colnames(pdata), value = TRUE)
message("\n寻找etiology...")
for (col in chars_cols) {
  vals <- as.character(pdata[[col]])
  if (any(grepl("HBV|HCV|Alcohol|NBNC|etiology", vals, ignore.case = TRUE))) {
    clinical$etiology <- gsub("^[^:]+:\\s*", "", trimws(vals))
    message(sprintf("✅ Etiology found in: %s", col))
    cat("etiology分布:\n")
    print(table(clinical$etiology))
    break
  }
}
# 如果char列没找到，打印前几个值
if (!"etiology" %in% colnames(clinical)) {
  for (col in chars_cols) {
    cat(sprintf("\n%s 前5值:", col))
    print(head(pdata[[col]], 5))
  }
}

# ── 5. 基因翻译hub基因 ────────────────────────────────────────────────────
hub_genes <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41", "RPS28")
hub_in_rc <- hub_genes %in% rownames(rc)
cat("\n═══ Hub基因在GSE141198 readcounts中 ═══\n")
for (i in seq_along(hub_genes)) {
  cat(sprintf("  %s: %s\n", hub_genes[i],
              ifelse(hub_in_rc[i], "✅ 找到", "❌ 缺失")))
}
message(sprintf("\n%d/%d hub基因在GSE141198中", sum(hub_in_rc), length(hub_genes)))

# 如果基因名是其他格式（Ensembl），可能需要转换
# 打印hub基因可能的Ensembl ID
if (sum(hub_in_rc) == 0) {
  message("\n⚠ 所有hub基因缺失！基因ID可能是Ensembl格式")
  message("前20个基因名（确认格式）:")
  print(head(rownames(rc), 20))
}

# ── 6. 基因过滤（DESeq2标准）────────────────────────────────────────────────
# Raw counts: 保留在≥20%样本中≥10 reads的基因
min_samples <- floor(0.2 * ncol(rc))
keep <- rowSums(rc >= 10) >= min_samples
message(sprintf("\n基因过滤: %d/%d (%.1f%%) 在≥%d样本中≥10 reads",
                sum(keep), length(keep), 100*mean(keep), min_samples))

rc_filt <- rc[keep, ]

# ── 7. 保存 ──────────────────────────────────────────────────────────────────
saveRDS(rc_filt, file.path(PROJ_DIR, "GSE141198_counts_filt.rds"))
saveRDS(rc, file.path(PROJ_DIR, "GSE141198_counts_raw.rds"))
saveRDS(clinical, file.path(PROJ_DIR, "GSE141198_clinical.rds"))

message("\n═══ GSE141198数据总结 ═══")
message(sprintf("样本数: %d", ncol(rc_filt)))
message(sprintf("基因数: %d (过滤后)", nrow(rc_filt)))
message(sprintf("OS events: %d / %d (%.0f%%)",
                sum(clinical$os_status == 1, na.rm = TRUE),
                sum(!is.na(clinical$os_status)),
                100 * mean(clinical$os_status == 1, na.rm = TRUE)))
message(sprintf("Hub基因: %d/7 在表达矩阵中", sum(hub_in_rc)))
message("\n✅ 数据预处理完成！")
