# =============================================================================
# 生成 Figure 6: 机制模型——ATF4/ISR + MYC 双通路驱动翻译调控的镜像模型
# =============================================================================

library(ggplot2)
library(dplyr)
library(grid)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ── HCC侧（左侧）元素 ─────────────────────────────────────────────────────
hcc <- list(
  # 疾病标签
  disease = data.frame(x = -2, y = 9.2, label = "Hepatocellular\nCarcinoma", color = "#C62828"),
  # 上游TF
  myc  = data.frame(x = -3, y = 7, label = "MYC Program\n(Pathway ρ=+0.613)", fill = "#FFCDD2"),
  atf4 = data.frame(x = -1, y = 7, label = "ATF4 / ISR\n(TF ρ=+0.439)", fill = "#FFCDD2"),
  dd   = data.frame(x = -2, y = 6.2, label = "DDIT3 (+0.289)\nXBP1 (+0.162)", fill = "#FFE0B2", size = 3),
  # 核糖体
  ribo = data.frame(x = -2, y = 4.5, label = "Ribosome Biogenesis\n& Translation Initiation\n(7 Hub Genes Co-expressed)", fill = "#EF9A9A"),
  # 输出
  out  = data.frame(x = -2, y = 2.8, label = "TRANSLATION ↑\nProliferation Support", fill = "#E53935"),
  # 注释
  note = data.frame(x = -2, y = 1.3, label = "Pro-survival ISR adaptation\n+ MYC-driven proliferation", color = "#B71C1C")
)

# ── HF侧（右侧）元素 ──────────────────────────────────────────────────────
hf <- list(
  disease = data.frame(x = 2, y = 9.2, label = "Heart Failure\n(Dilated/Ischemic)", color = "#1565C0"),
  stress = data.frame(x = 2, y = 7, label = "Energy Depletion\nMetabolic Stress", fill = "#BBDEFB"),
  mTOR  = data.frame(x = 2, y = 5.8, label = "mTORC1 ↓\n(MTOR ρ=−0.391)", fill = "#BBDEFB"),
  ribo  = data.frame(x = 2, y = 4.5, label = "Ribosome & Translation\nGenes Downregulated\n(7 Hub Genes Discordant)", fill = "#90CAF9"),
  out   = data.frame(x = 2, y = 2.8, label = "TRANSLATION ↓\nAdaptive Energy\nConservation", fill = "#1E88E5"),
  note  = data.frame(x = 2, y = 1.3, label = "ATP-costly translation suppressed\nas cardiomyocyte survival strategy", color = "#0D47A1")
)

# ── 中间共享元素 ─────────────────────────────────────────────────────────
mid <- list(
  module = data.frame(x = 0, y = 3.7, label = "Co-expression Module\nConserved Across Diseases\n(WGCNA Replication\np=6.3×10⁻¹²)", fill = "#E1BEE7"),
  cross  = data.frame(x = 0, y = 5.8, label = "Cross-Disease Pathway\nEffect Size Negative\nCorrelation\n(ssGSEA ρ=−0.598\np=0.0003)", fill = "#CE93D8")
)

# ── 绘图函数 ──────────────────────────────────────────────────────────────
add_rect <- function(df, w = 1.0, h = 0.7, alpha = 0.9) {
  geom_rect(data = df,
    aes(xmin = x - w, xmax = x + w, ymin = y - h, ymax = y + h),
    fill = df$fill, color = "grey40", linewidth = 0.8, alpha = alpha)
}

add_text <- function(df, size = 3.8, fontface = "plain", color = "grey20") {
  geom_text(data = df,
    aes(x = x, y = y, label = label),
    size = size, lineheight = 0.9, fontface = fontface, color = color)
}

# ── 箭头数据 ──────────────────────────────────────────────────────────────
arrows <- data.frame(
  x    = c(-3, -1, -2, -2, -2, -2, -2, 2, 2, 2, 2, -2,  2, 0, 0),
  xend = c(-3, -1, -2, -2, -2, -2, -2, 2, 2, 2, 2,  2, -2, 0, 0),
  y    = c(7.6, 7.6, 6.9, 7.6, 5.2, 5.2, 3.5, 7.6, 6.5, 5.2, 3.5, 5.2, 5.2, 4.4, 6.5),
  yend = c(5.2, 5.2, 5.2, 5.2, 3.5, 3.5, 3.5, 5.2, 5.2, 3.5, 3.5, 3.5, 3.5, 4.4, 4.5),
  color = c("#C62828", "#C62828", "#C62828", "#C62828", "#C62828", "#C62828", "#C62828",
            "#1565C0", "#1565C0", "#1565C0", "#1565C0", "#1565C0", "#1565C0", "#7B1FA2", "#7B1FA2")
)

p <- ggplot() +

  # ── 疾病标题 ────────────────────────────────────────────────────────────
  annotate("text", x = -2, y = 9.8, label = "HCC", fontface = "bold",
           size = 7, color = "#C62828") +
  annotate("text", x = 2, y = 9.8, label = "HF", fontface = "bold",
           size = 7, color = "#1565C0") +
  annotate("text", x = 0, y = 10.2, label = "vs", size = 4, color = "grey50") +

  # ── 分隔线 ──────────────────────────────────────────────────────────────
  annotate("segment", x = 0, xend = 0, y = 0.5, yend = 9.3,
           linetype = "dashed", color = "grey70", linewidth = 0.8) +
  annotate("text", x = 0, y = 9.5, label = "MIRROR", size = 3, color = "grey60",
           fontface = "italic") +

  # ── HCC侧的框 ───────────────────────────────────────────────────────────
  add_rect(hcc$myc,  w = 1.3, h = 0.55) +
  add_rect(hcc$atf4, w = 1.3, h = 0.55) +
  add_rect(hcc$dd,   w = 1.1, h = 0.45) +
  add_rect(hcc$ribo, w = 1.6, h = 0.65) +
  add_rect(hcc$out,  w = 1.4, h = 0.55) +

  add_text(hcc$myc,  size = 3.2) +
  add_text(hcc$atf4, size = 3.2) +
  add_text(hcc$dd,   size = 2.8, color = "#BF360C") +
  add_text(hcc$ribo, size = 3.0, fontface = "bold") +
  add_text(hcc$out,  size = 3.5, fontface = "bold", color = "white") +
  add_text(hcc$note, size = 2.7, color = "#B71C1C", fontface = "italic") +

  # ── HF侧的框 ────────────────────────────────────────────────────────────
  add_rect(hf$stress, w = 1.3, h = 0.55) +
  add_rect(hf$mTOR,   w = 1.3, h = 0.55) +
  add_rect(hf$ribo,   w = 1.6, h = 0.65) +
  add_rect(hf$out,    w = 1.4, h = 0.55) +

  add_text(hf$stress, size = 3.2) +
  add_text(hf$mTOR,   size = 3.0, color = "#0D47A1") +
  add_text(hf$ribo,   size = 3.0, fontface = "bold") +
  add_text(hf$out,    size = 3.5, fontface = "bold", color = "white") +
  add_text(hf$note,   size = 2.7, color = "#0D47A1", fontface = "italic") +

  # ── 中间元素 ────────────────────────────────────────────────────────────
  add_rect(mid$cross,  w = 1.6, h = 0.65) +
  add_rect(mid$module, w = 1.6, h = 0.65) +

  add_text(mid$cross,  size = 3.0, fontface = "bold", color = "#4A148C") +
  add_text(mid$module, size = 3.0, fontface = "bold", color = "#4A148C") +

  # ── 箭头 ────────────────────────────────────────────────────────────────
  # HCC侧
  annotate("segment", x = -3, xend = -3, y = 6.45, yend = 5.15,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           color = "#C62828", linewidth = 1.0) +
  annotate("segment", x = -1, xend = -1, y = 6.45, yend = 5.15,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           color = "#C62828", linewidth = 1.0) +
  annotate("segment", x = -2, xend = -2, y = 5.75, yend = 5.15,
           arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
           color = "#E57373", linewidth = 0.7) +
  # merge箭头
  annotate("segment", x = -2, xend = -2, y = 8.0, yend = 7.55,
           arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
           color = "#C62828", linewidth = 0.8) +

  # HCC ribo → output
  annotate("segment", x = -2, xend = -2, y = 3.85, yend = 3.35,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           color = "#C62828", linewidth = 1.2) +

  # HF侧
  annotate("segment", x = 2, xend = 2, y = 6.45, yend = 5.15,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           color = "#1565C0", linewidth = 1.0) +
  annotate("segment", x = 2, xend = 2, y = 3.85, yend = 3.35,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           color = "#1565C0", linewidth = 1.2) +

  # 交叉箭头 → 中间面板
  annotate("segment", x = -0.5, xend = 0, y = 5.15, yend = 5.15,
           arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
           color = "#7B1FA2", linewidth = 0.7, linetype = "dotted") +
  annotate("segment", x = 0.5, xend = 0, y = 5.15, yend = 5.15,
           arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
           color = "#7B1FA2", linewidth = 0.7, linetype = "dotted") +

  # 中间 → 模块
  annotate("segment", x = 0, xend = 0, y = 5.15, yend = 4.35,
           arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
           color = "#7B1FA2", linewidth = 0.8) +

  # ── HIF1A 标注 ──────────────────────────────────────────────────────────
  annotate("text", x = -3.8, y = 6.2, label = "HIF1A\n(ρ=−0.417)", size = 2.8,
           color = "#E57373", fontface = "italic", hjust = 0) +
  annotate("segment", x = -3.5, xend = -3.0, y = 6.2, yend = 6.6,
           arrow = arrow(length = unit(0.06, "cm"), type = "closed"),
           color = "#E57373", linewidth = 0.5, linetype = "dotted") +

  # ── 图例标注 ─────────────────────────────────────────────────────────────
  annotate("rect", xmin = -3.5, xmax = -2.5, ymin = 0.4, ymax = 0.9,
           fill = "#EF9A9A", color = "grey40", linewidth = 0.5) +
  annotate("text", x = -3, y = 0.65, label = "Source\nFunction", size = 2.5, lineheight = 0.9) +

  annotate("rect", xmin = -1.0, xmax = 1.0, ymin = 0.4, ymax = 0.9,
           fill = "#CE93D8", color = "grey40", linewidth = 0.5) +
  annotate("text", x = 0, y = 0.65, label = "Cross-Disease\nIntegration", size = 2.5, lineheight = 0.9) +

  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0.4, ymax = 0.9,
           fill = "#90CAF9", color = "grey40", linewidth = 0.5) +
  annotate("text", x = 2, y = 0.65, label = "Disease\nEffect", size = 2.5, lineheight = 0.9) +

  # 标题
  annotate("text", x = 0, y = 11.0,
           label = "Mechanistic Model: ATF4/ISR + MYC Drive Translation Regulation\nin HCC and Heart Failure",
           fontface = "bold", size = 5.5, color = "grey20") +

  coord_cartesian(xlim = c(-4.5, 4.5), ylim = c(0, 11.5)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(PROJ_DIR, "figures", "Figure6_Mechanistic_Model.png"),
       p, width = 13, height = 11, dpi = 300)
message("Figure 6 (Mechanistic Model) saved")
