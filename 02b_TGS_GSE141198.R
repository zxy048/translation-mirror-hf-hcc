# =============================================================================
# 脚本 02b：TGS在GSE141198中的预后验证
# 使用实际临床数据：OS (94 events), HBV/HCV/NBNC
# =============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
dir.create(file.path(PROJ_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

# ── 1. 加载数据 ──────────────────────────────────────────────────────────────
expr_vst <- readRDS(file.path(PROJ_DIR, "GSE141198_vst.rds"))
clinical <- readRDS(file.path(PROJ_DIR, "GSE141198_clinical.rds"))

hub_genes_6 <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41")
hub_genes_all <- c(hub_genes_6, "RPS28")

message("═══ TGS预后验证：GSE141198 (n=148, events=94) ═══")

# ── 2. 计算TGS ──────────────────────────────────────────────────────────────
genes_found <- intersect(hub_genes_6, rownames(expr_vst))
message(sprintf("TGS使用基因: %d/6: %s", length(genes_found),
                paste(genes_found, collapse=", ")))

expr_hub <- expr_vst[genes_found, , drop = FALSE]
expr_z <- t(scale(t(expr_hub)))
tgs_score <- colMeans(expr_z, na.rm = TRUE)
message(sprintf("TGS范围: %.3f ~ %.3f, median=%.3f", min(tgs_score), max(tgs_score), median(tgs_score)))

# ── 3. 匹配生存数据 ─────────────────────────────────────────────────────────
common_samples <- intersect(names(tgs_score), clinical$geo_id)
message(sprintf("TGS-生存匹配样本: %d", length(common_samples)))

surv_df <- data.frame(
  sample = common_samples,
  tgs = tgs_score[common_samples],
  os_time = clinical[common_samples, "os_time"],
  os_status = clinical[common_samples, "os_status"],
  etiology = clinical[common_samples, "etiology"],
  stringsAsFactors = FALSE
)

# 去除缺失值
surv_df <- surv_df[!is.na(surv_df$os_time) & !is.na(surv_df$os_status), ]
surv_df$os_time <- surv_df$os_time / 30.44  # 转换为月
message(sprintf("有效样本: %d, events: %d", nrow(surv_df), sum(surv_df$os_status)))

# ── 4. KM分析（TGS中位数分组）──────────────────────────────────────────────
surv_df$tgs_group <- ifelse(surv_df$tgs > median(surv_df$tgs), "High", "Low")
surv_df$tgs_group <- factor(surv_df$tgs_group, levels = c("Low", "High"))

fit_tgs <- survfit(Surv(os_time, os_status) ~ tgs_group, data = surv_df)
lr_p <- survdiff(Surv(os_time, os_status) ~ tgs_group, data = surv_df)
lr_pval <- 1 - pchisq(lr_p$chisq, df = 1)

message(sprintf("\n═══ KM结果 ═══"))
message(sprintf("Log-rank p = %.4f", lr_pval))
message(sprintf("Median OS Low: %.1f months", summary(fit_tgs)$table["tgs_group=Low", "median"]))
message(sprintf("Median OS High: %.1f months", summary(fit_tgs)$table["tgs_group=High", "median"]))

# KM图
p_km <- ggsurvplot(fit_tgs, data = surv_df,
                   pval = TRUE, pval.size = 5,
                   risk.table = TRUE, risk.table.height = 0.25,
                   palette = c("#377EB8", "#E41A1C"),
                   legend.labs = c("TGS Low", "TGS High"),
                   title = "GSE141198: TGS Prognostic Validation",
                   xlab = "Time (months)", ylab = "Overall Survival",
                   ggtheme = theme_minimal(base_size = 13))
ggsave(file.path(PROJ_DIR, "figures", "GSE141198_TGS_KM.png"),
       p_km$plot, width = 8, height = 6, dpi = 300)

# ── 5. Cox回归 ──────────────────────────────────────────────────────────────
# 连续变量
cox_cont <- coxph(Surv(os_time, os_status) ~ tgs, data = surv_df)
s_cont <- summary(cox_cont)

# 每IQR
iqr_val <- IQR(surv_df$tgs)
surv_df$tgs_per_iqr <- surv_df$tgs / iqr_val
cox_iqr <- coxph(Surv(os_time, os_status) ~ tgs_per_iqr, data = surv_df)
s_iqr <- summary(cox_iqr)

message(sprintf("\n═══ Cox回归 ═══"))
message(sprintf("连续变量TGS: HR=%.3f, 95%%CI %.3f-%.3f, p=%.4f",
                s_cont$conf.int[1], s_cont$conf.int[3], s_cont$conf.int[4],
                s_cont$coefficients[5]))
message(sprintf("每IQR TGS:    HR=%.3f, 95%%CI %.3f-%.3f, p=%.4f",
                exp(coef(cox_iqr)), exp(confint(cox_iqr)[1]),
                exp(confint(cox_iqr)[2]), s_iqr$coefficients[5]))

# 多变量（TGS + etiology）
if ("etiology" %in% colnames(surv_df) && length(unique(surv_df$etiology)) >= 2) {
  cox_mv <- coxph(Surv(os_time, os_status) ~ tgs_per_iqr + etiology, data = surv_df)
  message(sprintf("\n多变量 Cox (TGS + etiology):"))
  print(summary(cox_mv)$coefficients)
}

# ── 6. RPL39单基因验证 ──────────────────────────────────────────────────────
if ("RPL39" %in% rownames(expr_vst)) {
  rpl39_expr <- as.numeric(expr_vst["RPL39", common_samples])
  names(rpl39_expr) <- common_samples

  surv_df$rpl39 <- rpl39_expr[surv_df$sample]
  surv_df$rpl39_group <- ifelse(surv_df$rpl39 > median(surv_df$rpl39), "High", "Low")
  surv_df$rpl39_group <- factor(surv_df$rpl39_group, levels = c("Low", "High"))

  fit_rpl39 <- survfit(Surv(os_time, os_status) ~ rpl39_group, data = surv_df)
  lr_rpl39 <- survdiff(Surv(os_time, os_status) ~ rpl39_group, data = surv_df)
  lr_rpl39_p <- 1 - pchisq(lr_rpl39$chisq, df = 1)

  cox_rpl39 <- coxph(Surv(os_time, os_status) ~ rpl39, data = surv_df)
  s_rpl39 <- summary(cox_rpl39)

  message(sprintf("\n═══ RPL39 单基因验证 ═══"))
  message(sprintf("Log-rank p = %.4f", lr_rpl39_p))
  message(sprintf("Cox HR = %.3f (%.3f-%.3f), p = %.4f",
                  s_rpl39$conf.int[1], s_rpl39$conf.int[3],
                  s_rpl39$conf.int[4], s_rpl39$coefficients[5]))

  p_rpl39_km <- ggsurvplot(fit_rpl39, data = surv_df,
                           pval = TRUE, pval.size = 5,
                           palette = c("#377EB8", "#E41A1C"),
                           legend.labs = c("RPL39 Low", "RPL39 High"),
                           title = "GSE141198: RPL39 Validation",
                           xlab = "Time (months)", ylab = "Overall Survival",
                           ggtheme = theme_minimal(base_size = 13))
  ggsave(file.path(PROJ_DIR, "figures", "GSE141198_RPL39_KM.png"),
         p_rpl39_km$plot, width = 8, height = 6, dpi = 300)
}

# ── 7. 按病因分层的TGS亚组分析 ──────────────────────────────────────────────
if ("etiology" %in% colnames(surv_df)) {
  message("\n═══ TGS预后——病因亚组 ═══")
  for (etio in unique(surv_df$etiology)) {
    sub <- surv_df[surv_df$etiology == etio, ]
    if (nrow(sub) >= 20 && sum(sub$os_status) >= 10) {
      cox_sub <- coxph(Surv(os_time, os_status) ~ tgs, data = sub)
      s_sub <- summary(cox_sub)
      message(sprintf("  %s (n=%d, events=%d): HR=%.3f, p=%.4f",
                      etio, nrow(sub), sum(sub$os_status),
                      s_sub$conf.int[1], s_sub$coefficients[5]))
    }
  }
}

# ── 8. 保存结果 ──────────────────────────────────────────────────────────────
tgs_gse141198 <- list(
  cohort = "GSE141198",
  n = nrow(surv_df),
  events = sum(surv_df$os_status),
  tgs_logrank_p = lr_pval,
  tgs_cox_hr = s_cont$conf.int[1],
  tgs_cox_ci_lower = s_cont$conf.int[3],
  tgs_cox_ci_upper = s_cont$conf.int[4],
  tgs_cox_p = s_cont$coefficients[5],
  tgs_iqr_hr = exp(coef(cox_iqr)),
  tgs_iqr_p = s_iqr$coefficients[5],
  median_os_low = summary(fit_tgs)$table["tgs_group=Low", "median"],
  median_os_high = summary(fit_tgs)$table["tgs_group=High", "median"]
)

# 如果有RPL39结果也加入
if (exists("s_rpl39")) {
  tgs_gse141198$rpl39_logrank_p <- lr_rpl39_p
  tgs_gse141198$rpl39_cox_hr <- s_rpl39$conf.int[1]
  tgs_gse141198$rpl39_cox_p <- s_rpl39$coefficients[5]
}

saveRDS(tgs_gse141198, file.path(PROJ_DIR, "tgs_GSE141198_result.rds"))

# ── 9. 关键结论 ─────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════\n")
cat("  GSE141198 TGS 预后验证结果\n")
cat("─────────────────────────────────────\n")
cat(sprintf("  TGS Log-rank p:  %.4f\n", lr_pval))
cat(sprintf("  TGS Cox HR:      %.3f (%.3f-%.3f)\n",
            s_cont$conf.int[1], s_cont$conf.int[3], s_cont$conf.int[4]))
cat(sprintf("  TGS Cox p:       %.4f\n", s_cont$coefficients[5]))
if (exists("lr_rpl39_p")) {
  cat(sprintf("  RPL39 Log-rank p: %.4f\n", lr_rpl39_p))
  cat(sprintf("  RPL39 Cox HR:     %.3f, p=%.4f\n",
              s_rpl39$conf.int[1], s_rpl39$coefficients[5]))
}
cat("═══════════════════════════════════════════════\n")
