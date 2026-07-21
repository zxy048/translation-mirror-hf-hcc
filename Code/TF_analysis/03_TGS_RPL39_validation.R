# =============================================================================
# 脚本 02：翻译基因评分（TGS）+ RPL39 多队列预后验证
# 目标：用多基因评分替代单基因EEF1A1，在4个独立队列中验证
# 依赖：已有数据全部在 D:/R_projects/ 下
# =============================================================================

library(DESeq2)
library(limma)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(GSVA)
library(metafor)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
dir.create(file.path(PROJ_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJ_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

# ── 0. 定义hub基因 ────────────────────────────────────────────────────────────
# 7个hub基因（排除方向相反的RPS28，保留6个方向一致的）
hub_genes_6 <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41")
hub_genes_all <- c(hub_genes_6, "RPS28")  # 含RPS28作为对照

# 候选基因方向（基于原稿Table 1）：
# EEF1A1: HF上调, HCC下调 → 方向不一致，但两病均下调（需核实）
# 实际上原稿说：EEF1A1在HF(-0.115)和HCC(-0.269)均下调，方向一致✅
# RPL39: 两病均上调，效应量最大✅
# FAU: 两病方向一致✅
# RPS28: HF下调(-0.115), HCC上调(+0.514) 方向相反❌

# 方向一致的基因为hub_genes_6（排除RPS28）
message("翻译基因评分使用6个方向一致的hub基因: ", paste(hub_genes_6, collapse=", "))

# ═══════════════════════════════════════════════════════════════════════════════
# 第一部分：TCGA-LIHC 中构建 TGS 并做预后分析
# ═══════════════════════════════════════════════════════════════════════════════

message("\n", paste(rep("═", 70), collapse=""))
message("Part 1: TCGA-LIHC TGS construction and survival analysis")
message(paste(rep("═", 70), collapse=""))

# ── 1.1 加载TCGA-LIHC数据 ─────────────────────────────────────────────────────
se <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")
# 提取肿瘤样本表达矩阵
if ("shortLetterCode" %in% colnames(colData(se))) {
  tumor_idx <- which(colData(se)$shortLetterCode == "TP")
} else if ("sample_type" %in% colnames(colData(se))) {
  tumor_idx <- which(colData(se)$sample_type == "Primary Tumor" |
                     grepl("Tumor", colData(se)$sample_type))
} else {
  # 回退：用barcode判断 (TCGA barcode中01A表示原发肿瘤)
  tumor_idx <- grep("-01A|-01B", colnames(se))
}

message(sprintf("TCGA-LIHC 肿瘤样本数: %d", length(tumor_idx)))
expr_tumor <- assay(se, "HTSeq - Counts")[, tumor_idx]  # 或相应assay名
# 检查assay名称
message("可用assay: ", paste(assayNames(se), collapse=", "))

# 加载临床数据
clinical <- read.delim("D:/R_projects/lihc_tcga/lihc_tcga/data_clinical_patient.txt",
                       skip = 4, header = TRUE, stringsAsFactors = FALSE)
# 去除第一行（数据类型描述行）
# 实际需要根据文件结构调整

# ── 1.2 计算TGS ───────────────────────────────────────────────────────────────
# 方法1：ssGSEA法（推荐，类似通路评分）
# 方法2：z-score均值法（更简单）

# 确保基因在表达矩阵中
genes_available <- intersect(hub_genes_6, rownames(expr_tumor))
message(sprintf("6个hub基因中找到 %d 个: %s",
                length(genes_available),
                paste(genes_available, collapse=", ")))

# 使用z-score法计算TGS（每个样本的标准化均值）
expr_hub <- expr_tumor[genes_available, , drop = FALSE]
# log2转换（如果是counts）
if (max(expr_hub) > 1000) {
  expr_hub_log <- log2(expr_hub + 1)
} else {
  expr_hub_log <- expr_hub
}
# 每个基因z-score标准化
expr_z <- t(scale(t(expr_hub_log)))
# TGS = 每个样本的均值
tgs_score <- colMeans(expr_z, na.rm = TRUE)

# ── 1.3 TCGA-LIHC生存分析：TGS ──────────────────────────────────────────────
# 加载TCGA-LIHC生存数据（从ML_output加载之前分析的结果）
# 此处需要匹配临床数据
# 替代方案：使用TCGAbiolinks
library(TCGAbiolinks)

clinical_lihc <- GDCquery_clinic("TCGA-LIHC", type = "clinical")
surv_data <- data.frame(
  barcode = clinical_lihc$bcr_patient_barcode,
  os_time = clinical_lihc$days_to_death,
  os_status = ifelse(clinical_lihc$vital_status == "Dead", 1, 0),
  stage = clinical_lihc$ajcc_pathologic_stage,
  age = clinical_lihc$age_at_diagnosis,
  stringsAsFactors = FALSE
)

# OS时间为NA的用last_follow_up填补
na_idx <- is.na(surv_data$os_time)
surv_data$os_time[na_idx] <- clinical_lihc$days_to_last_follow_up[na_idx]
surv_data$os_status[na_idx] <- 0

# 简化barcode匹配（取前12位）
surv_data$barcode_short <- substr(surv_data$barcode, 1, 12)

# 匹配表达数据的barcode
tumor_barcodes <- colnames(expr_tumor)
tumor_short <- substr(tumor_barcodes, 1, 12)

# 创建分析用数据框（需要实际匹配逻辑，这里提供框架）
message("临床-表达匹配完成后，运行以下分析...")

# 保存中间结果
saveRDS(tgs_score, file.path(PROJ_DIR, "TCGA_LIHC_TGS_score.rds"))
message("✅ Part 1 完成：TGS在TCGA-LIHC中计算完毕")

# ═══════════════════════════════════════════════════════════════════════════════
# 第二部分：多队列TGS预后验证框架
# ═══════════════════════════════════════════════════════════════════════════════

message("\n", paste(rep("═", 70), collapse=""))
message("Part 2: Multi-cohort TGS validation framework")
message(paste(rep("═", 70), collapse=""))

# ── 2.1 定义可复用的分析函数 ─────────────────────────────────────────────────

#' 计算TGS评分
#' @param expr 表达矩阵（gene x sample），应为log2标准化后
#' @param genes 基因列表
calculate_tgs <- function(expr, genes) {
  genes_avail <- intersect(genes, rownames(expr))
  if (length(genes_avail) < 3) {
    warning(sprintf("仅找到 %d/%d 个基因，TGS可靠性降低", length(genes_avail), length(genes)))
  }
  expr_sub <- expr[genes_avail, , drop = FALSE]
  expr_z <- t(scale(t(expr_sub)))
  tgs <- colMeans(expr_z, na.rm = TRUE)
  return(list(tgs = tgs, genes_found = genes_avail, genes_missing = setdiff(genes, genes_avail)))
}

#' KM + Cox分析
#' @param tgs TGS评分向量
#' @param surv 生存数据 (time, status)
#' @param label 数据集标签
analyze_tgs_survival <- function(tgs, surv, label = "") {
  # 确保对齐
  common <- intersect(names(tgs), rownames(surv))
  tgs <- tgs[common]
  surv <- surv[common, ]

  # 中位数分组
  group <- ifelse(tgs > median(tgs), "High", "Low")
  surv$group <- factor(group, levels = c("Low", "High"))

  # KM
  fit <- survfit(Surv(time, status) ~ group, data = surv)

  # Log-rank
  lr_test <- survdiff(Surv(time, status) ~ group, data = surv)
  lr_p <- 1 - pchisq(lr_test$chisq, df = 1)

  # Cox (单变量)
  cox1 <- coxph(Surv(time, status) ~ tgs, data = surv)
  cox_summary <- summary(cox1)

  # Cox (连续变量，每IQR)
  iqr_val <- IQR(tgs)
  surv$tgs_per_iqr <- tgs / iqr_val
  cox_iqr <- coxph(Surv(time, status) ~ tgs_per_iqr, data = surv)

  list(
    fit = fit,
    logrank_p = lr_p,
    cox_hr = cox_summary$conf.int[1, 1],
    cox_ci_lower = cox_summary$conf.int[1, 3],
    cox_ci_upper = cox_summary$conf.int[1, 4],
    cox_p = cox_summary$coefficients[1, 5],
    cox_iqr_hr = exp(coef(cox_iqr)[1]),
    cox_iqr_p = summary(cox_iqr)$coefficients[1, 5],
    n_total = nrow(surv),
    n_events = sum(surv$status),
    median_os_low = ifelse(any(surv$group == "Low"),
                           summary(fit)$table["group=Low", "median"],
                           NA),
    median_os_high = summary(fit)$table["group=High", "median"],
    label = label
  )
}

message("✅ 分析函数定义完毕")

# ── 2.2 结果汇总表（在后续步骤中填充）────────────────────────────────────────
tgs_results <- data.frame(
  Cohort = character(),
  N = integer(),
  Events = integer(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  logrank_p = numeric(),
  Median_Low = numeric(),
  Median_High = numeric(),
  Platform = character(),
  stringsAsFactors = FALSE
)

message("\n✅ Part 2 完成：框架就绪")
message("提示：在确认TCGA-LIHC生存数据对齐后，依次对4个队列运行analyze_tgs_survival()")

# ═══════════════════════════════════════════════════════════════════════════════
# 第三部分：RPL39 单独验证
# ═══════════════════════════════════════════════════════════════════════════════

message("\n", paste(rep("═", 70), collapse=""))
message("Part 3: RPL39 independent validation")
message(paste(rep("═", 70), collapse=""))

#' RPL39单基因预后分析（跨队列）
#' RPL39是7个hub基因中效应量最大的（HCC log2FC=+0.587）
analyze_gene_survival <- function(expr_matrix, gene_name, surv_data, label = "") {
  # 基因表达
  if (!gene_name %in% rownames(expr_matrix)) {
    warning(sprintf("%s 不在表达矩阵中", gene_name))
    return(NULL)
  }
  expr_val <- as.numeric(expr_matrix[gene_name, ])
  names(expr_val) <- colnames(expr_matrix)

  # 匹配
  common <- intersect(names(expr_val), rownames(surv_data))
  expr_val <- expr_val[common]
  surv <- surv_data[common, ]

  # 中位数分组
  group <- ifelse(expr_val > median(expr_val), "High", "Low")
  surv$group <- factor(group, levels = c("Low", "High"))

  # KM + Cox
  fit <- survfit(Surv(time, status) ~ group, data = surv)
  cox1 <- coxph(Surv(time, status) ~ expr_val, data = surv)
  s <- summary(cox1)

  list(
    fit = fit,
    hr = s$conf.int[1, 1],
    ci_lower = s$conf.int[1, 3],
    ci_upper = s$conf.int[1, 4],
    p = s$coefficients[1, 5],
    n = nrow(surv),
    events = sum(surv$status),
    label = label,
    gene = gene_name
  )
}

message("✅ Part 3 完成：RPL39分析函数就绪")

# ═══════════════════════════════════════════════════════════════════════════════
# 第四部分：Meta分析（四队列汇总）
# ═══════════════════════════════════════════════════════════════════════════════

message("\n", paste(rep("═", 70), collapse=""))
message("Part 4: Meta-analysis across 4 cohorts")
message(paste(rep("═", 70), collapse=""))

#' Meta分析森林图
#' @param results_list 各队列analyze_tgs_survival()返回的列表
plot_meta_forest <- function(results_list, output_path = NULL) {
  # 提取HR和CI
  cohorts <- sapply(results_list, `[[`, "label")
  hrs <- sapply(results_list, `[[`, "cox_hr")
  ci_low <- sapply(results_list, `[[`, "cox_ci_lower")
  ci_upp <- sapply(results_list, `[[`, "cox_ci_upper")
  ns <- sapply(results_list, `[[`, "n_total")
  events <- sapply(results_list, `[[`, "n_events")

  # 随机效应meta分析
  log_hr <- log(hrs)
  se_log_hr <- (log(ci_upp) - log(ci_low)) / (2 * 1.96)

  meta_res <- rma(yi = log_hr, sei = se_log_hr, method = "REML")

  # 森林图数据
  forest_data <- data.frame(
    cohort = cohorts,
    n = ns,
    events = events,
    hr = hrs,
    ci_lower = ci_low,
    ci_upper = ci_upp,
    weight = weights(meta_res),
    stringsAsFactors = FALSE
  )

  # 输出
  cat("\n═══ Meta分析结果 ═══\n")
  cat(sprintf("合并HR: %.3f (95%%CI %.3f-%.3f)\n",
              exp(meta_res$b), exp(meta_res$ci.lb), exp(meta_res$ci.ub)))
  cat(sprintf("p值: %.4f\n", meta_res$pval))
  cat(sprintf("异质性 I²: %.1f%%\n", meta_res$I2))
  cat(sprintf("异质性p值: %.4f\n", meta_res$QEp))

  # 简单森林图（ggplot）
  p <- ggplot(forest_data, aes(x = hr, y = reorder(cohort, hr))) +
    geom_point(aes(size = events), color = "steelblue") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    scale_x_log10() +
    labs(x = "Hazard Ratio (log scale)", y = "",
         title = "TGS Prognostic Meta-Analysis Across Cohorts",
         subtitle = sprintf("Pooled HR=%.2f, 95%%CI %.2f-%.2f, p=%.4f, I²=%.1f%%",
                           exp(meta_res$b), exp(meta_res$ci.lb),
                           exp(meta_res$ci.ub), meta_res$pval, meta_res$I2)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 5, dpi = 300)
    message("森林图已保存至: ", output_path)
  }

  return(list(meta = meta_res, forest_data = forest_data, plot = p))
}

message("✅ 所有分析函数定义完毕")
message("\n═══════════════════════════════════════════════════════════════")
message("脚本02执行完毕。后续步骤：")
message("  1. 在RStudio中调整TCGA-LIHC临床数据匹配逻辑")
message("  2. 对4个队列分别运行TGS分析")
message("  3. 运行Meta分析")
message("═══════════════════════════════════════════════════════════════")

saveRDS(list(
  hub_genes = hub_genes_6,
  calculate_tgs = calculate_tgs,
  analyze_tgs_survival = analyze_tgs_survival,
  analyze_gene_survival = analyze_gene_survival,
  plot_meta_forest = plot_meta_forest
), file.path(PROJ_DIR, "analysis_functions.rds"))
