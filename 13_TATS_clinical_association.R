# =============================================================================
# 脚本 13：TATS (Translation-Associated Transcriptional Score) 计算与临床关联
# 基于 HF WGCNA 重跑的 green module (227 genes)
# 替代旧的 TGS（基于 6 hub genes）
# =============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(gridExtra)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
FIG_DIR  <- file.path(PROJ_DIR, "figures")
TBL_DIR  <- file.path(PROJ_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load green module genes ──────────────────────────────────────────────
wgcna <- readRDS(file.path(PROJ_DIR, "GSE57338_WGCNA_rerun.rds"))
primary_mod <- wgcna$primary_translation_module  # "green"
green_genes <- wgcna$gene_symbols[wgcna$moduleColors == primary_mod]
cat(sprintf("Green module: %d genes\n", length(green_genes)))
cat(sprintf("Top 10 hub genes (by |GS|*|MM|): %s\n",
            paste(wgcna$hub_genes_auto[1:10], collapse = ", ")))

# ── 2. TATS in GSE141198 ────────────────────────────────────────────────────
cat("\n═══ TATS: GSE141198 ═══\n")

expr_vst <- readRDS(file.path(PROJ_DIR, "GSE141198_vst.rds"))
clinical <- readRDS(file.path(PROJ_DIR, "GSE141198_clinical.rds"))

# Find green module genes present in GSE141198
green_in_gse141198 <- intersect(green_genes, rownames(expr_vst))
cat(sprintf("Green genes in GSE141198: %d/%d\n",
            length(green_in_gse141198), length(green_genes)))

# Compute TATS
expr_green <- expr_vst[green_in_gse141198, , drop = FALSE]
expr_z <- t(scale(t(expr_green)))
tats_score <- colMeans(expr_z, na.rm = TRUE)
cat(sprintf("TATS range: %.3f to %.3f, median=%.3f\n",
            min(tats_score), max(tats_score), median(tats_score)))

# Match survival data
common_samples <- intersect(names(tats_score), clinical$geo_id)
cat(sprintf("TATS-survival matched samples: %d\n", length(common_samples)))

surv_df <- data.frame(
  sample = common_samples,
  tats = tats_score[common_samples],
  os_time = clinical[common_samples, "os_time"],
  os_status = clinical[common_samples, "os_status"],
  etiology = clinical[common_samples, "etiology"],
  stringsAsFactors = FALSE
)
surv_df <- surv_df[!is.na(surv_df$os_time) & !is.na(surv_df$os_status), ]
surv_df$os_time <- surv_df$os_time / 30.44  # days to months
cat(sprintf("Valid samples: %d, events: %d\n", nrow(surv_df), sum(surv_df$os_status)))

# Median split
surv_df$tats_group <- ifelse(surv_df$tats > median(surv_df$tats), "High", "Low")
surv_df$tats_group <- factor(surv_df$tats_group, levels = c("Low", "High"))

# Kaplan-Meier
fit <- survfit(Surv(os_time, os_status) ~ tats_group, data = surv_df)
logrank_p <- survdiff(Surv(os_time, os_status) ~ tats_group, data = surv_df)$pvalue
cat(sprintf("Log-rank p = %.4f\n", logrank_p))

# Cox regression (continuous)
cox_cont <- coxph(Surv(os_time, os_status) ~ tats, data = surv_df)
cox_summary <- summary(cox_cont)
cat(sprintf("Cox continuous: HR=%.3f, 95%%CI %.3f-%.3f, p=%.4f\n",
            cox_summary$conf.int[1, "exp(coef)"],
            cox_summary$conf.int[1, "lower .95"],
            cox_summary$conf.int[1, "upper .95"],
            cox_summary$coefficients[1, "Pr(>|z|)"]))

# ── 3. TATS KM plot (GSE141198) ─────────────────────────────────────────────
p_km <- ggsurvplot(fit, data = surv_df,
  pval = TRUE, pval.size = 3.5,
  palette = c("#2166AC", "#B2182B"),
  xlab = "Overall Survival (months)",
  ylab = "Overall Survival Probability",
  title = sprintf("TATS in GSE141198 (n=%d, events=%d)", nrow(surv_df), sum(surv_df$os_status)),
  subtitle = sprintf("227-gene green module | Log-rank p=%.4f | Cox HR=%.3f (%.3f-%.3f)",
                     logrank_p,
                     cox_summary$conf.int[1, "exp(coef)"],
                     cox_summary$conf.int[1, "lower .95"],
                     cox_summary$conf.int[1, "upper .95"]),
  legend.title = "TATS",
  legend.labs = c("Low", "High"),
  risk.table = TRUE,
  risk.table.height = 0.2,
  ggtheme = theme_minimal(base_size = 11)
)
png(file.path(FIG_DIR, "Figure4A_TATS_GSE141198_KM.png"),
    width = 8, height = 6, units = "in", res = 300)
print(p_km)
dev.off()
cat("Figure 4A (TATS KM) saved.\n")

# ── 4. Etiology-stratified analysis ──────────────────────────────────────────
cat("\n--- Etiology-stratified TATS ---\n")
for (et in unique(surv_df$etiology)) {
  sub <- surv_df[surv_df$etiology == et, ]
  if (nrow(sub) < 20) next
  sub$group <- ifelse(sub$tats > median(sub$tats), "High", "Low")
  if (length(unique(sub$group)) < 2) next
  fit_sub <- survdiff(Surv(os_time, os_status) ~ group, data = sub)
  cat(sprintf("  %s (n=%d, events=%d): log-rank p=%.4f\n",
              et, nrow(sub), sum(sub$os_status), fit_sub$pvalue))
}

# ── 5. Compare TATS vs old TGS ──────────────────────────────────────────────
cat("\n--- Comparison: TATS vs old TGS ---\n")

# Old TGS (6 hub genes)
old_hub <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41")
old_found <- intersect(old_hub, rownames(expr_vst))
cat(sprintf("Old TGS genes found in GSE141198: %d/6: %s\n",
            length(old_found), paste(old_found, collapse = ", ")))

if (length(old_found) > 0) {
  expr_old <- expr_vst[old_found, , drop = FALSE]
  expr_old_z <- t(scale(t(expr_old)))
  tgs_score <- colMeans(expr_old_z, na.rm = TRUE)

  common_both <- intersect(intersect(names(tats_score), names(tgs_score)), clinical$geo_id)
  comp_df <- data.frame(
    TATS = tats_score[common_both],
    TGS = tgs_score[common_both],
    stringsAsFactors = FALSE
  )
  cor_tats_tgs <- cor(comp_df$TATS, comp_df$TGS, method = "spearman")
  cat(sprintf("Spearman TATS vs TGS: ρ = %.3f\n", cor_tats_tgs))
}

# ── 6. Summary table ─────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════\n")
cat("  TATS SUMMARY\n")
cat("═══════════════════════════════════════\n")
cat(sprintf("Module: %s (%d genes)\n", primary_mod, length(green_genes)))
cat(sprintf("Genes present in GSE141198: %d\n", length(green_in_gse141198)))
cat(sprintf("GSE141198: Log-rank p=%.4f, Cox continuous p=%.4f\n",
            logrank_p, cox_summary$coefficients[1, "Pr(>|z|)"]))
cat(sprintf("TATS vs TGS Spearman ρ = %.3f\n", cor_tats_tgs))
cat("\nConclusion: TATS is an exploratory measure of translation-related\n")
cat("transcriptional activity. It does NOT constitute a validated prognostic\n")
cat("biomarker and should be interpreted as module activity, not a clinical score.\n")

# Save TATS scores for downstream use
saveRDS(list(
  module = primary_mod,
  n_genes = length(green_genes),
  genes_used = green_in_gse141198,
  tats_scores = tats_score,
  gse141198 = list(
    logrank_p = logrank_p,
    cox_hr = cox_summary$conf.int[1, "exp(coef)"],
    cox_hr_lower = cox_summary$conf.int[1, "lower .95"],
    cox_hr_upper = cox_summary$conf.int[1, "upper .95"],
    cox_p = cox_summary$coefficients[1, "Pr(>|z|)"]
  ),
  tats_vs_tgs_rho = cor_tats_tgs
), file.path(PROJ_DIR, "TATS_results.rds"))

cat("\nTATS results saved to TATS_results.rds\n")
