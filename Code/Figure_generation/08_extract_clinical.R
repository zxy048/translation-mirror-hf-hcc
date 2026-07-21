# 09_extract_clinical_characteristics.R
# 从三个数据集的pdata提取临床基线特征，填充Table 1
library(GEOquery)
library(Biobase)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ═══════════════════════════════════════════════════════
# 1. GSE57338 (HF) —— 本地GEO系列矩阵
# ═══════════════════════════════════════════════════════
cat("═══ GSE57338 (HF) ═══\n")
gz_file <- file.path(PROJ_DIR, "GSE57338_series_matrix.txt.gz")
gse_hf <- getGEO(filename = gz_file, getGPL = FALSE)
pd_hf <- pData(gse_hf)

cat(sprintf("样本数: %d\n", nrow(pd_hf)))
cat("\n所有列名:\n")
print(colnames(pd_hf))

# 打印 characteristics 列的内容
chars_cols <- grep("characteristics|:ch1", colnames(pd_hf), value = TRUE)
cat(sprintf("\ncharacteristics列 (%d个):\n", length(chars_cols)))
for (cn in chars_cols) {
  vals <- unique(pd_hf[[cn]])
  cat(sprintf("\n[%s] (%d unique values):\n", cn, length(vals)))
  for (v in head(vals, 10)) {
    cat(sprintf("  '%s'\n", as.character(v)))
  }
}

# 提取关键变量
# 年龄
age_col <- grep("age", colnames(pd_hf), ignore.case = TRUE, value = TRUE)
cat(sprintf("\n年龄相关列: %s\n", paste(age_col, collapse=", ")))

# 性别
sex_col <- grep("sex|gender", colnames(pd_hf), ignore.case = TRUE, value = TRUE)
cat(sprintf("性别相关列: %s\n", paste(sex_col, collapse=", ")))

# 分组/诊断
grp_col <- grep("heart failure|disease|diagnosis|group|condition|status|cardiomyopathy",
                colnames(pd_hf), ignore.case = TRUE, value = TRUE)
cat(sprintf("疾病相关列: %s\n", paste(grp_col, collapse=", ")))

# ── 解析characteristics提取年龄性别 ──
extract_val <- function(x) {
  # 处理 "label: value" 或直接 "value"
  v <- as.character(x)
  if (all(grepl(":", v))) v <- gsub("^[^:]+:\\s*", "", v)
  trimws(v)
}

# 遍历characteristics找age和sex
for (cn in chars_cols) {
  vals <- unique(pd_hf[[cn]])
  vals_str <- paste(vals, collapse=" | ")
  if (grepl("year|age|month", vals_str, ignore.case = TRUE) &&
      any(grepl("^\\d+$", extract_val(vals)))) {
    cat(sprintf("\n🔴 年龄列候选: %s\n", cn))
    ages <- as.numeric(extract_val(pd_hf[[cn]]))
    cat(sprintf("  年龄: mean=%.1f, sd=%.1f, range=[%.0f-%.0f], NA=%d\n",
                mean(ages, na.rm=TRUE), sd(ages, na.rm=TRUE),
                min(ages, na.rm=TRUE), max(ages, na.rm=TRUE), sum(is.na(ages))))
  }
  if (grepl("^male$|^female$|^M$|^F$", vals_str, ignore.case = TRUE) ||
      grepl("sex|gender", cn, ignore.case = TRUE)) {
    cat(sprintf("\n🔴 性别列候选: %s\n", cn))
    print(table(extract_val(pd_hf[[cn]])))
  }
}

# 疾病亚型
for (cn in chars_cols) {
  vals <- unique(pd_hf[[cn]])
  vals_str <- paste(vals, collapse=" | ")
  if (grepl("DCM|ICM|dilated|ischemic|cardiomyopathy|NF|normal|failure|yes|no",
            vals_str, ignore.case = TRUE)) {
    cat(sprintf("\n🔴 疾病分类列: %s\n", cn))
    print(table(extract_val(pd_hf[[cn]])))
  }
}

# ═══════════════════════════════════════════════════════
# 2. TCGA-LIHC
# ═══════════════════════════════════════════════════════
cat("\n\n═══ TCGA-LIHC ═══\n")

# 尝试本地文件
tcga_local <- "D:/R_projects/lihc_tcga/lihc_tcga/data_clinical_patient.txt"
if (file.exists(tcga_local)) {
  cat("使用本地TCGA临床文件\n")
  clin_tcga <- read.delim(tcga_local, skip = 4, stringsAsFactors = FALSE)
  cat(sprintf("行数: %d, 列数: %d\n", nrow(clin_tcga), ncol(clin_tcga)))
  cat("列名:\n")
  print(head(colnames(clin_tcga), 30))

  # 关键变量
  # OS_MONTHS, OS_STATUS, AGE, SEX, AJCC_PATHOLOGIC_TUMOR_STAGE
  os_cols <- grep("OS|overall.survival|vital|death|follow", colnames(clin_tcga),
                  ignore.case = TRUE, value = TRUE)
  cat(sprintf("\nOS相关列: %s\n", paste(os_cols, collapse=", ")))

  age_col <- grep("age|diagnosis", colnames(clin_tcga), ignore.case = TRUE, value = TRUE)
  cat(sprintf("年龄列: %s\n", paste(age_col, collapse=", ")))

} else {
  cat("本地TCGA文件不存在，尝试TCGAbiolinks...\n")
  # 尝试用TCGAbiolinks
  if (requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    library(TCGAbiolinks)
    clin_tcga <- GDCquery_clinic("TCGA-LIHC", type = "clinical")
    cat(sprintf("获取到 %d 行\n", nrow(clin_tcga)))
    cat("列名:\n")
    print(colnames(clin_tcga))

    # 年龄
    if ("age_at_index" %in% colnames(clin_tcga)) {
      ages <- clin_tcga$age_at_index
      cat(sprintf("年龄: mean=%.1f, sd=%.1f\n", mean(ages, na.rm=TRUE), sd(ages, na.rm=TRUE)))
    }
    # 性别
    if ("gender" %in% colnames(clin_tcga)) {
      print(table(clin_tcga$gender))
    }
    # 分期
    if ("ajcc_pathologic_stage" %in% colnames(clin_tcga)) {
      print(table(clin_tcga$ajcc_pathologic_stage))
    }
    # 生存
    cat(sprintf("vital_status: %s\n", paste(table(clin_tcga$vital_status), collapse=", ")))
  } else {
    cat("⚠ TCGAbiolinks不可用，请从v6.1分析中获取TCGA-LIHC数据\n")
    cat("  → 脚本 02_TGS_RPL39_validation.R 中有TCGA数据提取代码\n")
  }
}

# ═══════════════════════════════════════════════════════
# 3. GSE141198 (HCC) —— 已有RDS
# ═══════════════════════════════════════════════════════
cat("\n\n═══ GSE141198 (HCC) ═══\n")
clinical_141198 <- readRDS(file.path(PROJ_DIR, "GSE141198_clinical.rds"))

cat(sprintf("样本数: %d\n", nrow(clinical_141198)))
cat("列名:\n")
print(colnames(clinical_141198))

# OS
cat(sprintf("\nOS events: %d/%d (%.0f%%)\n",
            sum(clinical_141198$os_status == 1, na.rm = TRUE),
            sum(!is.na(clinical_141198$os_status)),
            100 * mean(clinical_141198$os_status == 1, na.rm = TRUE)))
cat(sprintf("OS time (days): median=%.0f, range=[%.0f-%.0f]\n",
            median(clinical_141198$os_time, na.rm = TRUE),
            min(clinical_141198$os_time, na.rm = TRUE),
            max(clinical_141198$os_time, na.rm = TRUE)))

# 病因
if ("etiology" %in% colnames(clinical_141198)) {
  cat("\n病因分布:\n")
  print(table(clinical_141198$etiology))
}

# CTNNB1
if ("ctnnb1" %in% colnames(clinical_141198)) {
  cat("\nCTNNB1:\n")
  print(table(clinical_141198$ctnnb1))
}

cat("\n═══════════════════════════════════\n")
cat("运行完毕。将上述输出填入 Table 1。\n")
