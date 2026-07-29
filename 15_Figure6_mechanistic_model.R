# =============================================================================
# Script 15: Figure 5C — Mechanistic Model
# Disease-context-dependent remodeling of translation-related programs
# =============================================================================

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(png)

PROJ_DIR <- "D:/R_projects/revision_analysis"
FIG_DIR  <- file.path(PROJ_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5C: Mechanistic Model — Disease-context-dependent remodeling
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══ Figure 5C: Mechanistic Model ═══\n")

# Use grid graphics to create a conceptual diagram
# Increased height for bottom margin
png(file.path(FIG_DIR, "Figure5_Mechanistic_Model.png"),
    width = 14, height = 13.0, units = "in", res = 300)

grid.newpage()

# ── Title ──
grid.text(
  "Disease-Context-Dependent Remodeling of Translation-Related Transcriptional Programs",
  x = 0.02, y = 0.998, just = c(0, 1),
  gp = gpar(fontface = "bold", fontsize = 14, col = "black")
)

# ══════════════════════════════════════════════════════════════════
# TOP PANEL: Two diseases with opposite perturbation
# Expanded viewports: width 0.30→0.32, height 0.28→0.30
# ══════════════════════════════════════════════════════════════════

# --- Left: Heart Failure ---
pushViewport(viewport(x = 0.18, y = 0.83, width = 0.32, height = 0.30))

# HF outer box
grid.rect(gp = gpar(fill = "#2166AC20", col = "#2166AC", lwd = 2))
grid.text("HEART FAILURE", x = 0.5, y = 0.93,
          gp = gpar(fontface = "bold", fontsize = 12, col = "#2166AC"))
grid.text("Energy-depleted myocardium", x = 0.5, y = 0.82,
          gp = gpar(fontsize = 9, col = "grey30"))
grid.text("GSE57338 (n = 313)", x = 0.5, y = 0.73,
          gp = gpar(fontsize = 8, col = "grey50"))

# Green module inner box: taller + wider to avoid text clipping
grid.rect(x = 0.5, y = 0.53, width = 0.78, height = 0.27,
          gp = gpar(fill = "#228B2220", col = "#228B22", lwd = 1.5))
grid.text("Green Module", x = 0.5, y = 0.65,
          gp = gpar(fontface = "bold", fontsize = 9, col = "#228B22"))
grid.text("227 genes | r = -0.521", x = 0.5, y = 0.54,
          gp = gpar(fontsize = 8, col = "grey30"))
grid.text("Ribosome biogenesis", x = 0.5, y = 0.43,
          gp = gpar(fontsize = 7.5, col = "grey50"))

# Down arrow — shifted slightly
grid.text(expression("\u2193"), x = 0.5, y = 0.25,
          gp = gpar(fontsize = 20, col = "#2166AC"))
grid.text("Translation-related", x = 0.5, y = 0.16,
          gp = gpar(fontsize = 7.5, col = "#2166AC", fontface = "italic"))
grid.text("transcriptional suppression", x = 0.5, y = 0.08,
          gp = gpar(fontsize = 7.5, col = "#2166AC", fontface = "italic"))
popViewport()

# --- Center: Distinct organization vs Functional mirroring ---
pushViewport(viewport(x = 0.50, y = 0.83, width = 0.34, height = 0.30))

# Middle: comparison zone
grid.text("DISTINCT ORGANIZATION", x = 0.5, y = 0.95,
          gp = gpar(fontface = "bold", fontsize = 10, col = "grey40"))
grid.text("Distinct module architecture", x = 0.5, y = 0.86,
          gp = gpar(fontsize = 8, col = "grey50"))
grid.text("Different hub genes", x = 0.5, y = 0.78,
          gp = gpar(fontsize = 8, col = "grey50"))
grid.text("Gene overlap: 62/227", x = 0.5, y = 0.70,
          gp = gpar(fontsize = 8, col = "grey50"))
grid.text("OR = 1.7, p = 0.039", x = 0.5, y = 0.62,
          gp = gpar(fontsize = 8, col = "grey50"))

# Divider
grid.lines(x = c(0.1, 0.9), y = c(0.54, 0.54),
           gp = gpar(col = "grey70", lwd = 1, lty = "dashed"))

grid.text("FUNCTIONAL MIRRORING", x = 0.5, y = 0.46,
          gp = gpar(fontface = "bold", fontsize = 10, col = "grey40"))
grid.text(expression("24/33 mirror pathways"), x = 0.5, y = 0.36,
          gp = gpar(fontsize = 8, col = "grey30"))
grid.text(expression("Spearman " * rho * " = -0.598"), x = 0.5, y = 0.28,
          gp = gpar(fontsize = 8, col = "grey30"))
grid.text("Same functional themes", x = 0.5, y = 0.20,
          gp = gpar(fontsize = 8, col = "grey50"))
grid.text("Opposite directions", x = 0.5, y = 0.12,
          gp = gpar(fontsize = 8, col = "grey50"))
popViewport()

# --- Right: HCC ---
pushViewport(viewport(x = 0.82, y = 0.83, width = 0.32, height = 0.30))

# HCC outer box
grid.rect(gp = gpar(fill = "#B2182B10", col = "#B2182B", lwd = 2))
grid.text("HEPATOCELLULAR", x = 0.5, y = 0.93,
          gp = gpar(fontface = "bold", fontsize = 12, col = "#B2182B"))
grid.text("CARCINOMA", x = 0.5, y = 0.85,
          gp = gpar(fontface = "bold", fontsize = 12, col = "#B2182B"))
grid.text("Proliferative tumor", x = 0.5, y = 0.75,
          gp = gpar(fontsize = 9, col = "grey30"))
grid.text("GSE141198 (n = 148)", x = 0.5, y = 0.67,
          gp = gpar(fontsize = 8, col = "grey50"))

# Blue module inner box: taller + wider to avoid text clipping
grid.rect(x = 0.5, y = 0.48, width = 0.78, height = 0.27,
          gp = gpar(fill = "#2166AC20", col = "#2166AC", lwd = 1.5))
grid.text("Blue Module", x = 0.5, y = 0.60,
          gp = gpar(fontface = "bold", fontsize = 9, col = "#2166AC"))
grid.text("1,315 genes", x = 0.5, y = 0.49,
          gp = gpar(fontsize = 8, col = "grey30"))
grid.text("Cytoplasmic translation", x = 0.5, y = 0.38,
          gp = gpar(fontsize = 7.5, col = "grey50"))

# Up arrow — shifted slightly
grid.text(expression("\u2191"), x = 0.5, y = 0.22,
          gp = gpar(fontsize = 20, col = "#B2182B"))
grid.text("Translation-related", x = 0.5, y = 0.13,
          gp = gpar(fontsize = 7.5, col = "#B2182B", fontface = "italic"))
grid.text("transcriptional activation", x = 0.5, y = 0.05,
          gp = gpar(fontsize = 7.5, col = "#B2182B", fontface = "italic"))
popViewport()

# ══════════════════════════════════════════════════════════════════
# MIDDLE PANEL: Same-organ controls
# ══════════════════════════════════════════════════════════════════

grid.text("Disease-Context Dependence: Same-Organ Controls",
          x = 0.5, y = 0.54, gp = gpar(fontface = "bold", fontsize = 11, col = "black"))

# HCM control
pushViewport(viewport(x = 0.20, y = 0.42, width = 0.22, height = 0.18))
grid.rect(gp = gpar(fill = "grey95", col = "grey60", lwd = 1.5))
grid.text("Cardiac Control", x = 0.5, y = 0.85,
          gp = gpar(fontface = "bold", fontsize = 9, col = "grey30"))
grid.text("HCM vs. NF", x = 0.5, y = 0.65,
          gp = gpar(fontsize = 8, col = "grey40"))
grid.text(expression(rho * " = -0.036"), x = 0.5, y = 0.48,
          gp = gpar(fontsize = 9, col = "grey40"))
grid.text("No mirror perturbation", x = 0.5, y = 0.32,
          gp = gpar(fontsize = 8, col = "grey50"))
grid.text(expression("\u2192"), x = 0.5, y = 0.15,
          gp = gpar(fontsize = 10, col = "grey50"))
grid.text("Not a broad cardiac effect", x = 0.5, y = 0.05,
          gp = gpar(fontsize = 7, col = "grey50", fontface = "italic"))
popViewport()

# Liver control
pushViewport(viewport(x = 0.42, y = 0.42, width = 0.22, height = 0.18))
grid.rect(gp = gpar(fill = "grey95", col = "grey60", lwd = 1.5))
grid.text("Liver Control", x = 0.5, y = 0.85,
          gp = gpar(fontface = "bold", fontsize = 9, col = "grey30"))
grid.text("Cirrhosis vs. Normal", x = 0.5, y = 0.65,
          gp = gpar(fontsize = 8, col = "grey40"))
grid.text(expression(rho * " = +0.402"), x = 0.5, y = 0.48,
          gp = gpar(fontsize = 9, col = "grey40"))
grid.text("Positive correlation", x = 0.5, y = 0.32,
          gp = gpar(fontsize = 8, col = "grey50"))
grid.text(expression("\u2192"), x = 0.5, y = 0.15,
          gp = gpar(fontsize = 10, col = "grey50"))
grid.text("Progressive remodeling", x = 0.5, y = 0.05,
          gp = gpar(fontsize = 7, col = "grey50", fontface = "italic"))
popViewport()

# Cohort variability
pushViewport(viewport(x = 0.64, y = 0.42, width = 0.22, height = 0.18))
grid.rect(gp = gpar(fill = "grey95", col = "grey60", lwd = 1.5))
grid.text("Cohort Dependency", x = 0.5, y = 0.85,
          gp = gpar(fontface = "bold", fontsize = 9, col = "grey30"))
grid.text("GSE116250 (n = 64)", x = 0.5, y = 0.65,
          gp = gpar(fontsize = 8, col = "grey40"))
grid.text("14/33 mirror pathways", x = 0.5, y = 0.48,
          gp = gpar(fontsize = 9, col = "grey40"))
grid.text("Partial replication", x = 0.5, y = 0.32,
          gp = gpar(fontsize = 8, col = "grey50"))
grid.text(expression("\u2192"), x = 0.5, y = 0.15,
          gp = gpar(fontsize = 10, col = "grey50"))
grid.text("Severity/etiology dependent", x = 0.5, y = 0.05,
          gp = gpar(fontsize = 7, col = "grey50", fontface = "italic"))
popViewport()

# ══════════════════════════════════════════════════════════════════
# BOTTOM PANEL: Regulatory dimensions
# ══════════════════════════════════════════════════════════════════

grid.text("Upstream regulatory dimensions (HCC)",
          x = 0.5, y = 0.24, gp = gpar(fontface = "bold", fontsize = 11, col = "black"))

# ATF4/ISR
pushViewport(viewport(x = 0.3, y = 0.10, width = 0.35, height = 0.20))
grid.rect(gp = gpar(fill = "#E8F5E9", col = "#228B22", lwd = 1.5))
grid.text("Stress-Responsive Axis", x = 0.5, y = 0.90,
          gp = gpar(fontface = "bold", fontsize = 10, col = "#228B22"))
grid.text(expression(bold("ATF4/ISR") * " (" * rho * " = +0.500)"), x = 0.5, y = 0.72,
          gp = gpar(fontsize = 9, col = "black"))
grid.text("DDIT3 (CHOP): \u03c1 = +0.338", x = 0.5, y = 0.54,
          gp = gpar(fontsize = 8, col = "grey30"))
grid.text("XBP1: \u03c1 = +0.180", x = 0.5, y = 0.38,
          gp = gpar(fontsize = 8, col = "grey30"))
grid.text("Cellular stress associated with", x = 0.5, y = 0.25,
          gp = gpar(fontsize = 8, col = "grey50", fontface = "italic"))
grid.text("translation-related remodeling", x = 0.5, y = 0.17,
          gp = gpar(fontsize = 8, col = "grey50", fontface = "italic"))
popViewport()

# MYC
pushViewport(viewport(x = 0.7, y = 0.10, width = 0.35, height = 0.20))
grid.rect(gp = gpar(fill = "#FCE4EC", col = "#B2182B", lwd = 1.5))
grid.text("Proliferative Axis", x = 0.5, y = 0.90,
          gp = gpar(fontface = "bold", fontsize = 10, col = "#B2182B"))
grid.text(expression(bold("MYC Pathway") * " (" * rho * " = +0.753)"), x = 0.5, y = 0.72,
          gp = gpar(fontsize = 9, col = "black"))
grid.text("MYC Targets V1: \u03c1 = +0.459", x = 0.5, y = 0.54,
          gp = gpar(fontsize = 8, col = "grey30"))
grid.text("E2F Targets: \u03c1 = +0.409", x = 0.5, y = 0.38,
          gp = gpar(fontsize = 8, col = "grey30"))
grid.text("Proliferative demand associated", x = 0.5, y = 0.25,
          gp = gpar(fontsize = 8, col = "grey50", fontface = "italic"))
grid.text("with ribosome biogenesis", x = 0.5, y = 0.17,
          gp = gpar(fontsize = 8, col = "grey50", fontface = "italic"))
popViewport()

# Complementarity label — moved down away from regulatory boxes
grid.text("Complementary — not competing — regulatory dimensions",
          x = 0.5, y = 0.022, just = c(0.5, 0.5),
          gp = gpar(fontsize = 9, col = "grey40", fontface = "italic"))

dev.off()

cat(sprintf("Figure 5C saved: %s\n",
            file.path(FIG_DIR, "Figure5_Mechanistic_Model.png")))
cat("\n=== Figure 5C generation complete ===\n")
