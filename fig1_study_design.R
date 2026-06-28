# =============================================================================
# 生成 Figure 1: 研究设计流程图 (Graphical Abstract)
# =============================================================================

library(ggplot2)
library(dplyr)
library(grid)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ── 定义流程图的节点和边 ─────────────────────────────────────────────────────

# 三层结构：Data → Analysis → Findings
nodes <- data.frame(
  id = 1:18,
  x = c(
    # Row 1: Data sources
    1, 4, 7,
    # Row 2: Analysis - HF track
    1.5, 1.5,
    # Row 2: Analysis - Cross-disease
    4, 4, 4,
    # Row 2: Analysis - HCC track
    6.5, 6.5, 6.5,
    # Row 3: Key findings
    1, 3, 5, 7,
    2, 4, 6
  ),
  y = c(
    # Row 1: Data sources
    10, 10, 10,
    # Row 2: Analysis - HF
    8, 7,
    # Row 2: Analysis - Cross
    8, 7, 6,
    # Row 2: Analysis - HCC
    8, 7, 6,
    # Row 3: Findings
    4.5, 4.5, 4.5, 4.5,
    3, 3, 3
  ),
  label = c(
    "GSE57338\nHF (n=313)", "TCGA-LIHC\nHCC (n=371)", "GSE141198\nHCC (n=148)",
    "WGCNA\n(signed network)", "7 Hub Genes\n(EEF1A1, FAU,\nRPL39, RPL3,\nRPL32, RPL41, RPS28)",
    "ssGSEA\n83 Pathways", "Direction\nConsistency\nTest", "Upstream TF\nPrediction\n(19 TFs)",
    "WGCNA\n(Independent)", "TGS Survival\nValidation", "ssGSEA\nPathway Activity",
    "Module\nConservation",
    "Cross-Disease\nMirror Effect",
    "Prognostic\nSpecificity",
    "ATF4/ISR +\nMYC Regulation",
    "✓ Translation module\nreplicated (p=6.3e-12)",
    "Translation pathways\nρ = −0.598 (p=0.0003)",
    "TGS NS in 3 external\ncohorts; ATF4 top TF\n(ρ=+0.439, FDR<0.0001)"
  ),
  type = c("data", "data", "data",
           "analysis", "output",
           "analysis", "analysis", "analysis",
           "analysis", "analysis", "analysis",
           "finding", "finding", "finding", "finding",
           "conclusion", "conclusion", "conclusion"),
  stringsAsFactors = FALSE
)

# 颜色映射
fill_colors <- c(
  data = "#E8EAF6",
  analysis = "#FFF3E0",
  output = "#FFE0B2",
  finding = "#C8E6C9",
  conclusion = "#E0F7FA"
)

nodes$fill <- fill_colors[nodes$type]
nodes$border <- c(
  data = "#3F51B5", analysis = "#FF9800", output = "#E65100",
  finding = "#4CAF50", conclusion = "#0097A7"
)[nodes$type]

# ── 边（箭头） ───────────────────────────────────────────────────────────────
edges <- data.frame(
  from = c(1, 1, 2, 2, 2, 3, 3, 3,
           4, 4, 6, 7, 8, 9, 10, 11,
           5, 5, 5,
           12, 13, 14, 15),
  to   = c(4, 6, 6, 7, 8, 9, 10, 11,
           5, 5, 12, 13, 14, 16, 16, 17,
           7, 8, 14,
           17, 16, 18, 18),
  stringsAsFactors = FALSE
)

# ── 绘图 ─────────────────────────────────────────────────────────────────────

p <- ggplot() +
  # 背景面板标题
  annotate("rect", xmin = 0.2, xmax = 7.8, ymin = 9.2, ymax = 10.8,
           fill = "#E8EAF6", alpha = 0.3, color = NA) +
  annotate("text", x = 4, y = 10.9, label = "Data Sources",
           fontface = "bold", size = 4.5, color = "#283593") +

  annotate("rect", xmin = 0.2, xmax = 7.8, ymin = 5.5, ymax = 8.8,
           fill = "#FFF3E0", alpha = 0.3, color = NA) +
  annotate("text", x = 4, y = 8.9, label = "Analysis Pipeline",
           fontface = "bold", size = 4.5, color = "#E65100") +

  annotate("rect", xmin = 0.2, xmax = 7.8, ymin = 1.5, ymax = 5.2,
           fill = "#C8E6C9", alpha = 0.2, color = NA) +
  annotate("text", x = 4, y = 5.3, label = "Key Findings",
           fontface = "bold", size = 4.5, color = "#2E7D32") +

  # 边（箭头）
  geom_segment(data = edges,
    aes(x = nodes$x[from], y = nodes$y[from] - 0.3,
        xend = nodes$x[to], yend = nodes$y[to] + 0.3),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    color = "grey60", linewidth = 0.5) +

  # 节点框
  geom_rect(data = nodes,
    aes(xmin = x - 0.85, xmax = x + 0.85,
        ymin = y - 0.49, ymax = y + 0.49),
    fill = nodes$fill, color = nodes$border, linewidth = 1.2, alpha = 0.9) +

  # 节点文本
  geom_text(data = nodes,
    aes(x = x, y = y, label = label),
    size = 2.8, lineheight = 0.9, color = "grey20") +

  # 关键连接标注
  annotate("text", x = 2.0, y = 8.6, label = "→",
           size = 5, color = "#FF9800", fontface = "bold") +
  annotate("text", x = 6.0, y = 8.6, label = "→",
           size = 5, color = "#FF9800", fontface = "bold") +

  # 标题
  annotate("text", x = 4, y = 11.5,
           label = "Study Design: Cross-Disease Translation Program Analysis",
           fontface = "bold", size = 5.5, color = "grey20") +

  coord_cartesian(xlim = c(0, 8), ylim = c(1.2, 12)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(PROJ_DIR, "figures", "Figure1_Study_Design.png"),
       p, width = 14, height = 11, dpi = 300)
message("Figure 1 (Study Design) saved")
