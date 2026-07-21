# graphical_abstract.R
# Generate Graphical Abstract for BBA-MBD submission
# Concept: mirror regulation framework — conserved architecture, opposed direction, distinct upstream programs

library(ggplot2)
library(grid)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ── Color palette ──
red    <- "#C62828"
blue   <- "#1565C0"
purp   <- "#7B1FA2"
dark   <- "grey20"
mid    <- "grey50"
light  <- "grey85"

# ── Build with ggplot2 + theme_void ──

p <- ggplot() +

  # ═══════════════ TOP ROW: Disease labels + arrows ═══════════════

  # Section header
  annotate("text", x = 0.50, y = 0.97,
           label = "Different Diseases — Opposite Translational Demands",
           fontface = "bold", size = 5.5, color = dark) +

  # HCC label
  annotate("text", x = 0.26, y = 0.91,
           label = "HCC", fontface = "bold", size = 9, color = red) +
  annotate("text", x = 0.26, y = 0.87,
           label = "Hepatocellular Carcinoma", size = 4, color = red) +

  # HF label
  annotate("text", x = 0.74, y = 0.91,
           label = "HF", fontface = "bold", size = 9, color = blue) +
  annotate("text", x = 0.74, y = 0.87,
           label = "Heart Failure", size = 4, color = blue) +

  # Direction arrows
  annotate("text", x = 0.18, y = 0.82, label = "↑  UP", fontface = "bold",
           size = 5.5, color = red) +
  annotate("text", x = 0.82, y = 0.82, label = "↓  DOWN", fontface = "bold",
           size = 5.5, color = blue) +

  # Center vertical divider
  annotate("segment", x = 0.50, xend = 0.50, y = 0.17, yend = 0.78,
           linetype = "dashed", color = light, linewidth = 1.2) +

  # ═══════════════ MIDDLE: Conserved program + Mirror ═══════════════

  # Conserved program box (full width, spanning both sides)
  annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.63, ymax = 0.78,
           fill = "#F3E5F5", color = purp, linewidth = 1.8) +
  annotate("text", x = 0.50, y = 0.74,
           label = "Conserved Translational Co-expression Program",
           fontface = "bold", size = 5.5, color = "#4A148C") +
  annotate("text", x = 0.50, y = 0.69,
           label = "WGCNA: replicated module structure  |  ssGSEA: 24/33 translation pathways mirror (HCC↑ HF↓)",
           size = 3.3, color = "#6A1B9A") +
  annotate("text", x = 0.50, y = 0.65,
           label = "ρ = -0.598, p = 0.0003",
           size = 3.3, color = "#6A1B9A") +

  # Down arrow from conserved program
  annotate("segment", x = 0.50, xend = 0.50, y = 0.625, yend = 0.58,
           arrow = arrow(length = unit(0.08, "npc")),
           color = purp, linewidth = 1.8) +

  # Mirror regulation
  annotate("text", x = 0.50, y = 0.555,
           label = "Mirror Regulation",
           fontface = "bold", size = 6, color = purp) +
  annotate("text", x = 0.50, y = 0.525,
           label = "Descriptive framework: conserved architecture deployed in opposite directions",
           size = 3.5, color = mid) +

  # ═══════════════ CENTER ARROW ═══════════════
  annotate("segment", x = 0.50, xend = 0.50, y = 0.50, yend = 0.44,
           arrow = arrow(length = unit(0.08, "npc")),
           color = dark, linewidth = 1.5) +

  # ═══════════════ BOTTOM: Distinct upstream associations ═══════════════

  annotate("text", x = 0.50, y = 0.41,
           label = "Distinct Upstream Associations",
           fontface = "bold", size = 5, color = dark) +

  # HCC upstream box
  annotate("rect", xmin = 0.10, xmax = 0.44, ymin = 0.22, ymax = 0.38,
           fill = "#FFEBEE", color = red, linewidth = 1.2) +
  annotate("text", x = 0.27, y = 0.35,
           label = "HCC: Pro-growth Drivers", fontface = "bold", size = 4.5, color = "#B71C1C") +
  annotate("text", x = 0.27, y = 0.31,
           label = "ATF4 / ISR", fontface = "bold", size = 4, color = red) +
  annotate("text", x = 0.27, y = 0.28,
           label = "ρ = +0.439, FDR < 0.0001",
           size = 3.2, color = "#BF360C") +
  annotate("text", x = 0.27, y = 0.25,
           label = "MYC Targets V2", fontface = "bold", size = 4, color = red) +
  annotate("text", x = 0.27, y = 0.22,
           label = "ρ = +0.613, p < 0.0001",
           size = 3.2, color = "#BF360C") +

  # HF upstream box
  annotate("rect", xmin = 0.56, xmax = 0.90, ymin = 0.22, ymax = 0.38,
           fill = "#E3F2FD", color = blue, linewidth = 1.2) +
  annotate("text", x = 0.73, y = 0.35,
           label = "HF: Energy-stress Adaptation", fontface = "bold", size = 4.5, color = "#0D47A1") +
  annotate("text", x = 0.73, y = 0.31,
           label = "Chronic ATP Deficit", fontface = "bold", size = 4, color = blue) +
  annotate("text", x = 0.73, y = 0.28,
           label = "Adaptive suppression of anabolic processes",
           size = 3.2, color = "#1565C0") +
  annotate("text", x = 0.73, y = 0.25,
           label = "mTORC1 Downregulation", fontface = "bold", size = 4, color = blue) +
  annotate("text", x = 0.73, y = 0.22,
           label = "Coordinate transcriptional repression",
           size = 3.2, color = "#1565C0") +

  # ═══════════════ BOTTOM DECORATION ═══════════════

  annotate("text", x = 0.50, y = 0.12,
           label = "Key insight: Conserved molecular organization does not imply universal prognostic utility",
           size = 3.5, color = mid, fontface = "italic") +

  annotate("text", x = 0.50, y = 0.06,
           label = "Graphical Abstract  |  BBA – Molecular Basis of Disease",
           size = 3, color = light) +

  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

# ── Save ──
outfile <- file.path(PROJ_DIR, "figures", "Graphical_Abstract.png")
ggsave(outfile, p, width = 16, height = 12, dpi = 300)
cat(sprintf("→ Graphical Abstract saved: %s\n", outfile))
