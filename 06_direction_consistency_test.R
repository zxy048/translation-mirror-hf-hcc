# =============================================================================
# 脚本 06：7基因跨疾病方向一致性正式统计检验
# 目标：将描述性的"方向一致性观察"转化为正式的统计检验
# 为什么重要：这是全文最基础的证据，必须从"描述"升级为"检验"
# =============================================================================

library(ggplot2)
library(dplyr)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"

# ═══════════════════════════════════════════════════════════════════════════════
# 第一部分：7基因跨疾病log2FC方向一致性
# ═══════════════════════════════════════════════════════════════════════════════

message("═══ 7基因跨疾病方向一致性统计检验 ═══")

# 来自原论文Table 1的数值
cross_disease_data <- data.frame(
  Gene = c("EEF1A1", "RPL39", "RPL32", "FAU", "RPL41", "RPL3", "RPS28"),
  HF_log2FC = c(-0.115, 0.052, 0.071, 0.114, -0.016, 0.056, -0.115),
  HF_FDR_significant = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
  HCC_log2FC = c(-0.269, 0.587, 0.422, 0.193, 0.084, -0.034, 0.514),
  HCC_FDR_significant = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
  Direction_Consistent = c("Yes", "Yes", "Yes", "Yes", "Weak_no", "Weak_no", "NO"),
  stringsAsFactors = FALSE
)

# ── 1.1 Spearman相关：HF log2FC vs HCC log2FC ─────────────────────────────────
# 使用全部7个基因
cor_7 <- cor.test(cross_disease_data$HF_log2FC, cross_disease_data$HCC_log2FC,
                  method = "spearman")
message(sprintf("\n7基因 Spearman ρ = %.3f, p = %.4f", cor_7$estimate, cor_7$p.value))

# 排除RPS28（方向相反）后
data_6 <- subset(cross_disease_data, Gene != "RPS28")
cor_6 <- cor.test(data_6$HF_log2FC, data_6$HCC_log2FC, method = "spearman")
message(sprintf("6基因（排除RPS28） Spearman ρ = %.3f, p = %.4f",
                cor_6$estimate, cor_6$p.value))

# ── 1.2 置换检验：6个方向一致基因的log2FC相关是否显著偏离0 ──────────────────
n_perm <- 10000
perm_rhos <- numeric(n_perm)
for (i in 1:n_perm) {
  perm_rhos[i] <- cor(data_6$HF_log2FC, sample(data_6$HCC_log2FC),
                      method = "spearman")
}

perm_p <- mean(perm_rhos >= cor_6$estimate)  # 单侧（预期正相关）
message(sprintf("置换检验 (n=%d): ρ=%.3f, 经验p(单侧)=%.4f",
                n_perm, cor_6$estimate, perm_p))

# ── 1.3 二项检验：6个基因中5个方向一致的概率 ─────────────────────────────
n_consistent <- sum(data_6$Direction_Consistent == "Yes")  # 4个明确一致 + 2个弱不一致
# 更严格：仅计算清晰方向一致的
n_clear <- sum(data_6$Direction_Consistent == "Yes")  # EEF1A1, RPL39, RPL32, FAU = 4
binom_test <- binom.test(n_clear, n = 6, p = 0.5, alternative = "greater")
message(sprintf("二项检验: %d/6基因方向一致, p=%.4f (H0: p=0.5)",
                n_clear, binom_test$p.value))

# ═══════════════════════════════════════════════════════════════════════════════
# 第二部分：Bootstrap效应量置信区间
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ Bootstrap置信区间 ═══")

# 对6个基因的log2FC做bootstrap重抽样，估计跨疾病相关性的置信区间
n_boot <- 2000
boot_rhos <- numeric(n_boot)
for (i in 1:n_boot) {
  idx <- sample(1:6, size = 6, replace = TRUE)
  boot_rhos[i] <- cor(data_6$HF_log2FC[idx], data_6$HCC_log2FC[idx],
                      method = "spearman")
}
ci_95 <- quantile(boot_rhos, c(0.025, 0.975), na.rm = TRUE)
message(sprintf("Bootstrap 95%% CI for Spearman ρ: [%.3f, %.3f]", ci_95[1], ci_95[2]))

# ═══════════════════════════════════════════════════════════════════════════════
# 第三部分：可视化
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 跨疾病方向一致性图 ═══")

# ── 3.1 log2FC散点图 ──────────────────────────────────────────────────────────
p1 <- ggplot(cross_disease_data, aes(x = HF_log2FC, y = HCC_log2FC)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", alpha = 0.5) +
  geom_point(aes(color = Direction_Consistent, size = abs(HCC_log2FC)),
             alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = Gene),
    size = 4,
    box.padding = 0.5,
    max.overlaps = 10
  ) +
  scale_color_manual(
    values = c("Yes" = "#2CA02C", "Weak_no" = "#FF7F0E", "NO" = "#D62728"),
    name = "Direction Consistent"
  ) +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  labs(
    x = "HF log2FC (GSE57338)",
    y = "HCC log2FC (TCGA-LIHC)",
    title = "Cross-Disease Expression Direction Consistency",
    subtitle = sprintf("6 genes (excl. RPS28): Spearman ρ=%.3f, Permutation p=%.4f, Bootstrap 95%%CI [%.3f, %.3f]",
                      cor_6$estimate, perm_p, ci_95[1], ci_95[2]),
    caption = paste0(
      "Quadrant I: concordant upregulation; Quadrant III: concordant downregulation\n",
      "RPS28 excluded from analysis (discordant direction: HF↓ HCC↑)"
    )
  ) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#2CA02C", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#2CA02C", alpha = 0.05) +
  annotate("text", x = 0.3, y = max(cross_disease_data$HCC_log2FC) * 0.6,
           label = "Concordant ↑", color = "#2CA02C", size = 3.5, fontface = "italic") +
  annotate("text", x = -0.06, y = min(cross_disease_data$HCC_log2FC) * 0.6,
           label = "Concordant ↓", color = "#2CA02C", size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.caption = element_text(size = 9, color = "grey40", hjust = 0)
  )

ggsave(file.path(PROJ_DIR, "figures", "Figure_CrossDisease_DirectionTest.png"),
       p1, width = 9, height = 8, dpi = 300)
message("方向一致性散点图已保存")

# ── 3.2 跨疾病log2FC森林图 ────────────────────────────────────────────────────
# 展示每个基因在HF和HCC中的效应量和方向
plot_data <- cross_disease_data %>%
  tidyr::pivot_longer(
    cols = c(HF_log2FC, HCC_log2FC),
    names_to = "Disease",
    values_to = "log2FC"
  ) %>%
  mutate(
    Disease = ifelse(grepl("HF", Disease), "Heart Failure", "HCC"),
    sig_label = case_when(
      Disease == "HF" & HF_FDR_significant ~ "*",
      Disease == "HCC" & HCC_FDR_significant ~ "*",
      TRUE ~ ""
    )
  )

p2 <- ggplot(plot_data, aes(x = log2FC, y = Gene, color = Disease, shape = Disease)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", alpha = 0.5) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(aes(label = sig_label), vjust = 0.5, hjust = -1.5, size = 5, color = "black") +
  scale_color_manual(values = c("Heart Failure" = "#377EB8", "HCC" = "#E41A1C")) +
  scale_shape_manual(values = c("Heart Failure" = 16, "HCC" = 17)) +
  facet_wrap(~ Direction_Consistent, scales = "free_y", ncol = 1,
             labeller = labeller(Direction_Consistent = c(
               "Yes" = "✅ Direction Consistent",
               "Weak_no" = "⚠️ Weakly Inconsistent",
               "NO" = "❌ Direction Opposite"
             ))) +
  labs(
    x = "log2 Fold Change",
    y = "",
    title = "Cross-Disease Expression Direction of Hub Genes",
    subtitle = "* FDR < 0.05",
    caption = "Direction consistency assessed across HF (GSE57338) and HCC (TCGA-LIHC)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave(file.path(PROJ_DIR, "figures", "Figure_CrossDisease_Forest.png"),
       p2, width = 8, height = 8, dpi = 300)

# ═══════════════════════════════════════════════════════════════════════════════
# 第四部分：生成跨疾病一致性核心表格
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 核心结果表 ═══")

result_table <- data.frame(
  Analysis = c(
    "7-gene Spearman correlation (HF vs HCC log2FC)",
    "6-gene Spearman correlation (excluding RPS28)",
    "Permutation test (one-sided, 10,000 iter)",
    "Bootstrap 95% CI for ρ",
    "Binomial test (concordant direction)",
    "Concordant genes (out of 7)",
    "Discordant genes",
    "Weakly inconsistent genes"
  ),
  Result = c(
    sprintf("ρ = %.3f, p = %.4f", cor_7$estimate, cor_7$p.value),
    sprintf("ρ = %.3f, p = %.4f", cor_6$estimate, cor_6$p.value),
    sprintf("p = %.4f (expected positive correlation)", perm_p),
    sprintf("[%.3f, %.3f]", ci_95[1], ci_95[2]),
    sprintf("p = %.4f (H0: p=0.5)", binom_test$p.value),
    sprintf("%d (EEF1A1, RPL39, RPL32, FAU)", n_clear),
    "1 (RPS28: HF↓, HCC↑)",
    "2 (RPL41, RPL3: non-significant both)"
  ),
  stringsAsFactors = FALSE
)

print(result_table, row.names = FALSE)
write.csv(result_table, file.path(PROJ_DIR, "tables", "direction_consistency_results.csv"),
          row.names = FALSE)

message("\n✅ 脚本06完成")
message("关键发现：")
if (cor_6$estimate > 0.7 && perm_p < 0.05) {
  message("✅ 6个基因跨疾病方向一致性具有统计显著性——这是全文的核心发现之一")
} else {
  message("⚠ 方向一致性未达统计显著性，需在论文中诚实报告并讨论其局限性")
}
