# =============================================================================
# 修复 Fig 5A/5B：方框挡字 + 底部标签裁切
# - 扩大 x 轴范围容纳 ρ 标签
# - 增加底部边距
# - 调整柱形图宽度
# =============================================================================

library(ggplot2)
library(dplyr)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
FIG_DIR  <- file.path(PROJ_DIR, "figures")

# ── 加载数据 ──────────────────────────────────────────────────
tf_cors   <- readRDS(file.path(PROJ_DIR, "TF_TATS_correlation_results.rds"))
pathway_cors <- readRDS(file.path(PROJ_DIR, "Pathway_TATS_correlation_results.rds"))

# ═══════════════════════════════════════════════════════════════
# Figure 5A: TF-TATS correlation（左上）
# ═══════════════════════════════════════════════════════════════

cat("\n═══ Regenerating Figure 5A: TF-TATS ═══\n")

n_samples_tf <- 148  # GSE141198

tf_cors$TF <- factor(tf_cors$TF, levels = rev(tf_cors$TF))
tf_cors$Significant <- ifelse(tf_cors$FDR < 0.05, "FDR < 0.05", "NS")
tf_cors$Direction <- ifelse(tf_cors$Spearman_R > 0, "Positive", "Negative")

# 扩大 xlim 防止标签被切
x_min_tf <- min(tf_cors$Spearman_R)
x_max_tf <- max(tf_cors$Spearman_R)
x_pad_tf <- max(abs(c(x_min_tf, x_max_tf))) * 0.35  # 35% padding

p_tf <- ggplot(tf_cors, aes(x = Spearman_R, y = TF)) +
  geom_bar(stat = "identity", aes(fill = Direction), width = 0.65) +
  scale_fill_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +
  geom_text(aes(label = sprintf("ρ=%.3f%s", Spearman_R,
                                 ifelse(FDR < 0.05, "*", "")),
                hjust = ifelse(Spearman_R > 0, -0.15, 1.15)),
            size = 3.2, color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(
    title = "TF Expression vs. TATS",
    subtitle = sprintf("GSE141198 (n=%d) | 19 candidate TFs | * FDR < 0.05", n_samples_tf),
    x = "Spearman \u03c1 with TATS",
    y = ""
  ) +
  xlim(x_min_tf - x_pad_tf, x_max_tf + x_pad_tf) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 5),
    legend.position = "none"
  )

ggsave(file.path(FIG_DIR, "Figure_TF_TATS_correlation.png"),
       p_tf, width = 9, height = 6.5, dpi = 300)
cat("Figure 5A (TF-TATS) saved: 9x6.5 in, 300 dpi\n")

# ═══════════════════════════════════════════════════════════════
# Figure 5B: Pathway-TATS correlation（右上）
# ═══════════════════════════════════════════════════════════════

cat("\n═══ Regenerating Figure 5B: Pathway-TATS ═══\n")

n_samples_pw <- 148

pw_top <- head(pathway_cors, 20)
pw_top$Pathway <- factor(pw_top$Pathway, levels = rev(pw_top$Pathway))
pw_top$Direction <- ifelse(pw_top$Spearman_R > 0, "Positive", "Negative")

# 扩大 xlim
x_min_pw <- min(pw_top$Spearman_R)
x_max_pw <- max(pw_top$Spearman_R)
x_pad_pw <- max(abs(c(x_min_pw, x_max_pw))) * 0.35

# 缩短过长 pathway 名称用于 y 轴标签
shorten_label <- function(x, max_chars = 42) {
  x_char <- as.character(x)
  ifelse(nchar(x_char) > max_chars,
         paste0(substr(x_char, 1, max_chars - 3), "..."),
         x_char)
}
pw_top$PathwayShort <- shorten_label(pw_top$Pathway)
pw_top$PathwayShort <- factor(pw_top$PathwayShort,
                               levels = rev(levels(factor(pw_top$PathwayShort))))

p_pw <- ggplot(pw_top, aes(x = Spearman_R, y = PathwayShort)) +
  geom_bar(stat = "identity", aes(fill = Direction), width = 0.65) +
  scale_fill_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +
  geom_text(aes(label = sprintf("ρ=%.3f%s", Spearman_R,
                                 ifelse(FDR < 0.05, "*", "")),
                hjust = ifelse(Spearman_R > 0, -0.15, 1.15)),
            size = 3.0, color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(
    title = "Hallmark Pathway Activity vs. TATS",
    subtitle = sprintf("GSE141198 (n=%d) | Top 20 pathways | * FDR < 0.05", n_samples_pw),
    x = "Spearman \u03c1 with TATS",
    y = ""
  ) +
  xlim(x_min_pw - x_pad_pw, x_max_pw + x_pad_pw) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(t = 10, r = 15, b = 15, l = 5),
    legend.position = "none"
  )

ggsave(file.path(FIG_DIR, "Figure_Pathway_TATS_correlation.png"),
       p_pw, width = 10, height = 7.5, dpi = 300)
cat("Figure 5B (Pathway-TATS) saved: 10x7.5 in, 300 dpi\n")

cat("\n═══ Done ═══\n")
