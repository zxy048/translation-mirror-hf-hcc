# =============================================================================
# Figure 1 v2: 研究设计 —— 左→右流程，数据→分析→发现→结论
# =============================================================================

library(ggplot2)
library(grid)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# 颜色方案
col_data   <- "#5C6BC0"   # 靛蓝 - 数据
col_ana    <- "#FF9800"   # 橙 - 分析
col_find   <- "#43A047"   # 绿 - 发现
col_conc   <- "#E53935"   # 红 - 结论
col_bg     <- "#FAFAFA"

p <- ggplot() +

  # ═══════════════════ Row 1: DATA SOURCES ═══════════════════
  annotate("rect", xmin = 0.5,  xmax = 11.5, ymin = 9.5, ymax = 12.8,
           fill = "#E8EAF6", alpha = 0.5, color = NA) +
  annotate("text", x = 0.8, y = 12.5, label = "DATA", fontface = "bold",
           size = 5.5, color = col_data, hjust = 0) +

  # 三个数据集
  annotate("rect", xmin = 1.2,  xmax = 4.0, ymin = 9.8, ymax = 12.2,
           fill = "white", color = col_data, linewidth = 1.5) +
  annotate("text", x = 2.6, y = 11.5, label = "GSE57338", fontface = "bold",
           size = 4.5, color = "#283593") +
  annotate("text", x = 2.6, y = 10.8, label = "Heart Failure\nn = 313 (DCM/ICM + NF)\nAffymetrix HuGene 1.1 ST",
           size = 3.0, lineheight = 0.9, color = "grey40") +

  annotate("rect", xmin = 4.8,  xmax = 7.6, ymin = 9.8, ymax = 12.2,
           fill = "white", color = col_data, linewidth = 1.5) +
  annotate("text", x = 6.2, y = 11.5, label = "TCGA-LIHC", fontface = "bold",
           size = 4.5, color = "#283593") +
  annotate("text", x = 6.2, y = 10.8, label = "HCC Discovery\nn = 424 (371 T + 50 N)\nRNA-seq (Illumina)",
           size = 3.0, lineheight = 0.9, color = "grey40") +

  annotate("rect", xmin = 8.4,  xmax = 11.2, ymin = 9.8, ymax = 12.2,
           fill = "white", color = col_data, linewidth = 1.5) +
  annotate("text", x = 9.8, y = 11.5, label = "GSE141198", fontface = "bold",
           size = 4.5, color = "#283593") +
  annotate("text", x = 9.8, y = 10.8, label = "HCC Validation\nn = 148 (RNA-seq)\n94 OS events",
           size = 3.0, lineheight = 0.9, color = "grey40") +

  # Data → Analysis 箭头
  annotate("segment", x = 2.6, xend = 2.6, y = 9.75, yend = 8.9,
           arrow = arrow(length = unit(0.1, "cm")), color = col_ana, linewidth = 0.8) +
  annotate("segment", x = 6.2, xend = 5.5, y = 9.75, yend = 8.9,
           arrow = arrow(length = unit(0.1, "cm")), color = col_ana, linewidth = 0.8) +
  annotate("segment", x = 6.2, xend = 6.9, y = 9.75, yend = 8.9,
           arrow = arrow(length = unit(0.1, "cm")), color = col_ana, linewidth = 0.8) +
  annotate("segment", x = 9.8, xend = 9.8, y = 9.75, yend = 8.9,
           arrow = arrow(length = unit(0.1, "cm")), color = col_ana, linewidth = 0.8) +

  # ═══════════════════ Row 2: ANALYSIS ═══════════════════
  annotate("rect", xmin = 0.5,  xmax = 11.5, ymin = 5.8, ymax = 8.8,
           fill = "#FFF3E0", alpha = 0.5, color = NA) +
  annotate("text", x = 0.8, y = 8.55, label = "ANALYSIS", fontface = "bold",
           size = 5.5, color = col_ana, hjust = 0) +

  # 四个分析框
  analyses <- data.frame(
    xmn = c(1.2, 3.6,  6.0,  8.4),
    xmx = c(3.2, 5.6,  8.0, 10.4),
    lbl = c("WGCNA\n(Signed Network)", "ssGSEA\nPathway Activity", "Survival\nValidation", "Upstream\nTF Prediction"),
    sub = c("Module discovery\nHub gene identification\nGO enrichment",
            "81 pathways\nHallmark + KEGG + Reactome\nCross-disease effect sizes",
            "TGS (6-gene score)\nKM + Cox regression\nGSE141198 / GSE14520 / GSE76427",
            "19 TF candidates\nTF-TGS correlation\nMYC/E2F/mTORC1/ATF4"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(analyses)) {
    p <- p +
      annotate("rect",
               xmin = analyses$xmn[i], xmax = analyses$xmx[i],
               ymin = 6.1, ymax = 8.6,
               fill = "white", color = col_ana, linewidth = 1.3) +
      annotate("text",
               x = (analyses$xmn[i] + analyses$xmx[i]) / 2,
               y = 8.1,
               label = analyses$lbl[i],
               fontface = "bold", size = 3.8, color = "#E65100", lineheight = 0.9) +
      annotate("text",
               x = (analyses$xmn[i] + analyses$xmx[i]) / 2,
               y = 6.8,
               label = analyses$sub[i],
               size = 2.5, color = "grey40", lineheight = 0.9)
  }

  # Analysis → Findings 箭头
  for (cx in c(2.2, 4.6, 7.0, 9.4)) {
    annotate("segment", x = cx, xend = cx, y = 6.05, yend = 5.2,
             arrow = arrow(length = unit(0.1, "cm")), color = col_find, linewidth = 0.8)
  }

  # ═══════════════════ Row 3: FINDINGS ═══════════════════
  annotate("rect", xmin = 0.5,  xmax = 11.5, ymin = 2.0, ymax = 5.1,
           fill = "#E8F5E9", alpha = 0.5, color = NA) +
  annotate("text", x = 0.8, y = 4.85, label = "FINDINGS", fontface = "bold",
           size = 5.5, color = col_find, hjust = 0) +

  findings <- data.frame(
    xm = c(2.2, 4.6, 7.0, 9.4),
    icon = c("✓", "!", "✗", "★"),
    title = c("Module Conservation",
              "Mirror Effect",
              "Prognostic Specificity",
              "ATF4/ISR + MYC Axis"),
    desc = c("Translation module\nreplicated in GSE141198\nblue module p=6.3×10⁻¹²\n4/7 hub genes co-localized",
             "Cross-disease pathway\neffect sizes ρ = −0.598\np = 0.0003\nHCC↑ vs HF↓",
             "TGS NOT prognostic\nin 3 independent cohorts\nCox p = 0.479 (GSE141198)\nCohort-specific marker",
             "ATF4 strongest TF\nρ = +0.439, FDR<0.0001\nMYC pathway ρ = +0.613\nISR/UPR axis implicated"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(findings)) {
    p <- p +
      annotate("rect",
               xmin = findings$xm[i] - 1.2, xmax = findings$xm[i] + 1.2,
               ymin = 2.3, ymax = 4.9,
               fill = "white", color = col_find, linewidth = 1.3) +
      annotate("text", x = findings$xm[i], y = 4.4,
               label = findings$icon[i],
               size = 10, color = col_find, fontface = "bold") +
      annotate("text", x = findings$xm[i], y = 3.7,
               label = findings$title[i],
               fontface = "bold", size = 3.5, color = "#2E7D32") +
      annotate("text", x = findings$xm[i], y = 2.8,
               label = findings$desc[i],
               size = 2.5, color = "grey40", lineheight = 0.85)
  }

  # ═══════════════════ BOTTOM: CONCLUSION BANNER ═══════════════════
  annotate("rect", xmin = 0.5, xmax = 11.5, ymin = 0.3, ymax = 1.7,
           fill = "#FFEBEE", alpha = 0.7, color = col_conc, linewidth = 1.5) +
  annotate("text", x = 6.0, y = 1.2,
           label = "Translational co-expression is a conserved organizational principle\nacross HF and HCC, but perturbation direction is disease-type-specific",
           fontface = "bold", size = 4.2, color = "#B71C1C", lineheight = 0.9) +
  annotate("text", x = 6.0, y = 0.6,
           label = "ATF4/ISR stress adaptation + MYC proliferative program jointly regulate the translational program in HCC",
           size = 3.3, color = "#C62828", lineheight = 0.9) +

  # ═══════════════════ TITLE ═══════════════════
  annotate("text", x = 6.0, y = 13.3,
           label = "Cross-Disease Transcriptomic Analysis of Translation Regulation",
           fontface = "bold", size = 6.5, color = "grey20") +

  coord_cartesian(xlim = c(0, 12), ylim = c(0, 13.8)) +
  theme_void() +
  theme(plot.background = element_rect(fill = col_bg, color = NA),
        plot.margin = margin(15, 15, 15, 15))

ggsave(file.path(PROJ_DIR, "figures", "Figure1_Study_Design.png"),
       p, width = 16, height = 12, dpi = 300)
message("Figure 1 v2 saved")

# ═══════════════════════════════════════════════════════════════
# Figure 6 v2: 机制模型图 — 完全对称的镜像布局
# ═══════════════════════════════════════════════════════════════

# ── 所有坐标手动定义，确保完全对称 ──

# 左侧 HCC 元素
hcc <- list(
  disease = list(x = -4, y = 11.5, lab = "Hepatocellular\nCarcinoma", col = "#C62828", sz = 5),
  tf1 = list(x = -5.5, y = 9.2, lab = "ATF4 / ISR\nTF ρ = +0.439***", col = "#FFCDD2", sz = 3.2),
  tf2 = list(x = -2.5, y = 9.2, lab = "MYC Program\nPathway ρ = +0.613***", col = "#FFCDD2", sz = 3.2),
  sub = list(x = -4, y = 8.0, lab = "DDIT3 (+0.289)  XBP1 (+0.162)\nHIF1A (−0.417) negative feedback", col = "#FFE0B2", sz = 2.7),
  ribo = list(x = -4, y = 6.2, lab = "Ribosome Biogenesis\n& Translation Initiation\nEEF1A1  FAU  RPL39  RPL3\nRPL32  RPL41  RPS28", col = "#EF9A9A", sz = 3.0),
  out = list(x = -4, y = 4.2, lab = "TRANSLATION ↑", col = "#E53935", sz = 4.5),
  annot = list(x = -4, y = 2.8, lab = "Pro-survival ISR adaptation\n+ MYC-driven proliferation\n→ HCC growth", col = "#B71C1C", sz = 2.8)
)

# 右侧 HF 元素（对称于x=0）
hf <- list(
  disease = list(x = 4, y = 11.5, lab = "Heart Failure\n(DCM / ICM)", col = "#1565C0", sz = 5),
  tf1 = list(x = 5.5, y = 9.2, lab = "Energy Depletion\nMetabolic Stress", col = "#BBDEFB", sz = 3.2),
  tf2 = list(x = 2.5, y = 9.2, lab = "mTORC1 Signaling\nMTOR ρ = −0.391***", col = "#BBDEFB", sz = 3.2),
  sub = list(x = 4, y = 8.0, lab = "ATP-costly protein synthesis\nsuppressed as adaptation", col = "#E3F2FD", sz = 2.7),
  ribo = list(x = 4, y = 6.2, lab = "Ribosome & Translation\nGenes Downregulated\n7 Hub Genes Discordant", col = "#90CAF9", sz = 3.0),
  out = list(x = 4, y = 4.2, lab = "TRANSLATION ↓", col = "#1E88E5", sz = 4.5),
  annot = list(x = 4, y = 2.8, lab = "Energy conservation\nCardiomyocyte survival\n→ HF adaptation", col = "#0D47A1", sz = 2.8)
)

# 中间共享
mid <- list(
  ssgsea = list(x = 0, y = 7.0, lab = "Pathway Effect Size\nNegative Correlation\nssGSEA  ρ = −0.598\np = 0.0003", col = "#E1BEE7", sz = 3.2),
  wgcna = list(x = 0, y = 5.0, lab = "Co-expression Module\nConserved Across Diseases\nWGCNA Replication\np = 6.3 × 10⁻¹²", col = "#CE93D8", sz = 3.2),
)

p2 <- ggplot() +

  # ── 标题 ──
  annotate("text", x = 0, y = 13.5,
           label = "Mirror Regulation of Translation in HCC and Heart Failure",
           fontface = "bold", size = 7, color = "grey20") +
  annotate("text", x = 0, y = 12.8,
           label = "Shared co-expression architecture, opposed perturbation direction, distinct upstream regulators",
           size = 4, color = "grey50") +

  # ── 中间分割线 ──
  annotate("segment", x = 0, xend = 0, y = 1.5, yend = 11.0,
           linetype = "dashed", color = "grey70", linewidth = 0.8) +
  annotate("text", x = 0, y = 11.3, label = "M I R R O R", size = 3.5,
           color = "grey60", fontface = "italic") +

  # ── 疾病标签 ──
  annotate("text", x = -4, y = 12.3, label = "HCC", fontface = "bold",
           size = 9, color = "#C62828") +
  annotate("text", x = 4, y = 12.3, label = "HF", fontface = "bold",
           size = 9, color = "#1565C0") +

  # ── HCC 元素 ──
  annotate("rect", xmin = -6.5, xmax = -1.5, ymin = 8.6, ymax = 9.8,
           fill = "#FFCDD2", color = "#C62828", linewidth = 1) +
  annotate("text", x = -5.5, y = 9.2, label = "ATF4 / ISR\nTF ρ = +0.439***",
           size = 3.2, lineheight = 0.9, fontface = "bold", color = "#B71C1C") +
  annotate("text", x = -2.5, y = 9.2, label = "MYC Program\nPathway ρ = +0.613***",
           size = 3.2, lineheight = 0.9, fontface = "bold", color = "#B71C1C") +

  annotate("rect", xmin = -5.5, xmax = -2.5, ymin = 7.4, ymax = 8.5,
           fill = "#FFE0B2", color = "#E65100", linewidth = 0.8) +
  annotate("text", x = -4, y = 7.95,
           label = "DDIT3 (+0.289)  XBP1 (+0.162)\nHIF1A (−0.417)  ISR/UPR axis",
           size = 2.8, lineheight = 0.9, color = "#BF360C") +

  annotate("rect", xmin = -6.0, xmax = -2.0, ymin = 5.4, ymax = 7.0,
           fill = "#EF9A9A", color = "#C62828", linewidth = 1.2) +
  annotate("text", x = -4, y = 6.2,
           label = "Ribosome Biogenesis\n& Translation Initiation\nEEF1A1  FAU  RPL39  RPL3\nRPL32  RPL41  RPS28",
           size = 3.2, lineheight = 0.9, fontface = "bold", color = "#4A0000") +

  annotate("rect", xmin = -5.5, xmax = -2.5, ymin = 3.5, ymax = 4.8,
           fill = "#E53935", color = "#B71C1C", linewidth = 1.5) +
  annotate("text", x = -4, y = 4.2, label = "TRANSLATION ↑",
           size = 5, fontface = "bold", color = "white") +

  annotate("text", x = -4, y = 2.8,
           label = "Pro-survival ISR adaptation\n+ MYC-driven proliferation",
           size = 2.8, lineheight = 0.9, color = "#B71C1C", fontface = "italic") +

  # ── HF 元素 ──
  annotate("rect", xmin = 1.5, xmax = 6.5, ymin = 8.6, ymax = 9.8,
           fill = "#BBDEFB", color = "#1565C0", linewidth = 1) +
  annotate("text", x = 5.5, y = 9.2, label = "Energy Depletion\nMetabolic Stress",
           size = 3.2, lineheight = 0.9, fontface = "bold", color = "#0D47A1") +
  annotate("text", x = 2.5, y = 9.2, label = "mTORC1 Signaling\nMTOR ρ = −0.391***",
           size = 3.2, lineheight = 0.9, fontface = "bold", color = "#0D47A1") +

  annotate("rect", xmin = 2.5, xmax = 5.5, ymin = 7.4, ymax = 8.5,
           fill = "#E3F2FD", color = "#0D47A1", linewidth = 0.8) +
  annotate("text", x = 4, y = 7.95,
           label = "ATP-costly protein synthesis\nsuppressed as adaptation",
           size = 2.8, lineheight = 0.9, color = "#0D47A1") +

  annotate("rect", xmin = 2.0, xmax = 6.0, ymin = 5.4, ymax = 7.0,
           fill = "#90CAF9", color = "#1565C0", linewidth = 1.2) +
  annotate("text", x = 4, y = 6.2,
           label = "Ribosome & Translation\nGenes Downregulated\n7 Hub Genes Discordant",
           size = 3.2, lineheight = 0.9, fontface = "bold", color = "#002040") +

  annotate("rect", xmin = 2.5, xmax = 5.5, ymin = 3.5, ymax = 4.8,
           fill = "#1E88E5", color = "#0D47A1", linewidth = 1.5) +
  annotate("text", x = 4, y = 4.2, label = "TRANSLATION ↓",
           size = 5, fontface = "bold", color = "white") +

  annotate("text", x = 4, y = 2.8,
           label = "Energy conservation\nCardiomyocyte survival strategy",
           size = 2.8, lineheight = 0.9, color = "#0D47A1", fontface = "italic") +

  # ── 共享中间面板 ──
  annotate("rect", xmin = -1.8, xmax = 1.8, ymin = 6.2, ymax = 7.8,
           fill = "#E1BEE7", color = "#7B1FA2", linewidth = 1.2) +
  annotate("text", x = 0, y = 7.0,
           label = "Pathway Effect Size\nNegative Correlation\nssGSEA  ρ = −0.598\np = 0.0003",
           size = 3.2, lineheight = 0.9, fontface = "bold", color = "#4A148C") +

  annotate("rect", xmin = -1.8, xmax = 1.8, ymin = 4.2, ymax = 5.8,
           fill = "#CE93D8", color = "#7B1FA2", linewidth = 1.2) +
  annotate("text", x = 0, y = 5.0,
           label = "Co-expression Module\nConserved Across Diseases\nWGCNA Replication\np = 6.3 × 10⁻¹²",
           size = 3.2, lineheight = 0.9, fontface = "bold", color = "#4A148C") +

  # ── 箭头 ──

  # HCC 箭头
  annotate("segment", x = -5.5, xend = -5.5, y = 8.5, yend = 7.9,
           arrow = arrow(length = unit(0.1, "cm")), color = "#C62828", linewidth = 1.0) +
  annotate("segment", x = -2.5, xend = -2.5, y = 8.5, yend = 7.9,
           arrow = arrow(length = unit(0.1, "cm")), color = "#C62828", linewidth = 1.0) +
  annotate("segment", x = -4, xend = -4, y = 7.3, yend = 7.1,
           arrow = arrow(length = unit(0.08, "cm")), color = "#E65100", linewidth = 0.7) +

  annotate("segment", x = -4, xend = -4, y = 5.3, yend = 4.9,
           arrow = arrow(length = unit(0.12, "cm")), color = "#C62828", linewidth = 1.3) +

  # HF 箭头
  annotate("segment", x = 5.5, xend = 5.5, y = 8.5, yend = 7.9,
           arrow = arrow(length = unit(0.1, "cm")), color = "#1565C0", linewidth = 1.0) +
  annotate("segment", x = 2.5, xend = 2.5, y = 8.5, yend = 7.9,
           arrow = arrow(length = unit(0.1, "cm")), color = "#1565C0", linewidth = 1.0) +
  annotate("segment", x = 4, xend = 4, y = 5.3, yend = 4.9,
           arrow = arrow(length = unit(0.12, "cm")), color = "#1565C0", linewidth = 1.3) +

  # 横向连接 → 中间面板
  annotate("segment", x = -2.0, xend = -1.9, y = 6.9, yend = 6.9,
           arrow = arrow(length = unit(0.06, "cm")), color = "#7B1FA2", linewidth = 0.6, linetype = "dotted") +
  annotate("segment", x = 2.0, xend = 1.9, y = 6.9, yend = 6.9,
           arrow = arrow(length = unit(0.06, "cm")), color = "#7B1FA2", linewidth = 0.6, linetype = "dotted") +

  # 中间箭头向下
  annotate("segment", x = 0, xend = 0, y = 6.1, yend = 5.9,
           arrow = arrow(length = unit(0.1, "cm")), color = "#7B1FA2", linewidth = 1.0) +

  # ── 图例 ──
  annotate("rect", xmin = -6, xmax = -4, ymin = 1.2, ymax = 1.8,
           fill = "#FFCDD2", color = "#C62828", linewidth = 0.5) +
  annotate("text", x = -5, y = 1.5, label = "Upstream\nRegulator", size = 2.5, lineheight = 0.9) +

  annotate("rect", xmin = -2, xmax = 0, ymin = 1.2, ymax = 1.8,
           fill = "#EF9A9A", color = "#C62828", linewidth = 0.5) +
  annotate("text", x = -1, y = 1.5, label = "Translation\nProgram", size = 2.5, lineheight = 0.9) +

  annotate("rect", xmin = 1, xmax = 3, ymin = 1.2, ymax = 1.8,
           fill = "#E1BEE7", color = "#7B1FA2", linewidth = 0.5) +
  annotate("text", x = 2, y = 1.5, label = "Cross-Disease\nIntegration", size = 2.5, lineheight = 0.9) +

  annotate("rect", xmin = 4, xmax = 6, ymin = 1.2, ymax = 1.8,
           fill = "#E53935", color = "#B71C1C", linewidth = 0.5) +
  annotate("text", x = 5, y = 1.5, label = "Physiological\nOutput", size = 2.5, lineheight = 0.9, color = "white") +

  coord_cartesian(xlim = c(-7.5, 7.5), ylim = c(1, 14)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "#FAFAFA", color = NA),
        plot.margin = margin(15, 15, 15, 15))

ggsave(file.path(PROJ_DIR, "figures", "Figure6_Mechanistic_Model.png"),
       p2, width = 18, height = 14, dpi = 300)
message("Figure 6 v2 saved")

cat("\n=== Both figures regenerated ===\n")
cat("D:/R_projects/revision_analysis/figures/Figure1_Study_Design.png\n")
cat("D:/R_projects/revision_analysis/figures/Figure6_Mechanistic_Model.png\n")
