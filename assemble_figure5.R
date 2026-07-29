# =============================================================================
# 拼合 Figure 5: 5A (TF-TATS) + 5B (Pathway-TATS) + 5C (Model)
# Layout: A | B (top row), C spanning full width (bottom)
# =============================================================================

library(ggplot2)
library(grid)
library(png)

PROJ_DIR <- "D:/R_projects/revision_analysis"
FIG_DIR  <- file.path(PROJ_DIR, "figures")

read_png_grob <- function(path) {
  if (!file.exists(path)) return(nullGrob())
  img <- readPNG(path)
  rasterGrob(img, interpolate = TRUE)
}

cat("\n═══ Assembling Combined Figure 5 ═══\n")

# Load panels
gA <- read_png_grob(file.path(FIG_DIR, "Figure_TF_TATS_correlation.png"))
gB <- read_png_grob(file.path(FIG_DIR, "Figure_Pathway_TATS_correlation.png"))
gC <- read_png_grob(file.path(FIG_DIR, "Figure5_Mechanistic_Model.png"))

# Combined figure: 2 rows
# Row 1: A (left) + B (right) — wider to avoid squeezing
# Row 2: C (full width) — more bottom margin
png(file.path(FIG_DIR, "Figure5_Combined.png"),
    width = 20, height = 19, units = "in", res = 300)

grid.newpage()

# Main title
grid.text(
  "Figure 5. Transcription factor associations, pathway activity, and conceptual model\nof translation-related transcriptional remodeling in HCC.",
  x = 0.02, y = 0.998, just = c(0, 1),
  gp = gpar(fontface = "bold", fontsize = 14, lineheight = 0.9)
)

# Layout: 2 rows with more breathing room
top_row_y    <- 0.975   # top of first row
row1_height  <- 0.46    # height for top row
row2_height  <- 0.42    # height for bottom row
gap          <- 0.025   # gap between rows

# Row 1: Panel A (left half) + Panel B (right half)
pushViewport(viewport(x = 0.02, y = top_row_y,
                       width = 0.96, height = row1_height,
                       just = c(0, 1)))
grid.text("A", x = 0.00, y = 1.00, just = c(0, 1),
          gp = gpar(fontface = "bold", fontsize = 14))
grid.text("B", x = 0.49, y = 1.00, just = c(0, 1),
          gp = gpar(fontface = "bold", fontsize = 14))

# A: left half (0 to 0.49)
pushViewport(viewport(x = 0.00, y = 0.02,
                       width = 0.485, height = 0.92,
                       just = c(0, 0)))
grid.draw(gA)
popViewport()

# B: right half (0.49 to 0.98)
pushViewport(viewport(x = 0.49, y = 0.02,
                       width = 0.485, height = 0.92,
                       just = c(0, 0)))
grid.draw(gB)
popViewport()

popViewport()  # end row 1

# Row 2: Panel C (full width)
row2_y <- top_row_y - row1_height - gap
pushViewport(viewport(x = 0.02, y = row2_y,
                       width = 0.96, height = row2_height,
                       just = c(0, 1)))
grid.text("C", x = 0.00, y = 1.00, just = c(0, 1),
          gp = gpar(fontface = "bold", fontsize = 14))

pushViewport(viewport(x = 0.00, y = 0.02,
                       width = 1.00, height = 0.92,
                       just = c(0, 0)))
grid.draw(gC)
popViewport()

popViewport()  # end row 2

# Caption at bottom
grid.text(
  "(A) Spearman correlation of 19 candidate transcription factors with TATS in GSE141198 tumor samples (n = 148). Red bars: positive correlation (FDR < 0.05); blue bars: negative correlation (FDR < 0.05); gray: not significant. ATF4 is the strongest individual TF correlate (\u03c1 = +0.500, FDR < 0.0001). (B) Spearman correlation of Hallmark pathway ssGSEA scores with TATS. MYC Targets V2 pathway activity shows the strongest overall association (\u03c1 = +0.753, p < 0.0001), exceeding MYC mRNA level alone (\u03c1 = +0.255). ATF4/ISR (stress-responsive) and MYC (proliferative) represent complementary transcriptional dimensions associated with translation-related transcriptional remodeling. (C) Conceptual model of disease-context-dependent remodeling.",
  x = 0.02, y = 0.005, just = c(0, 0),
  gp = gpar(fontsize = 8, col = "grey40", lineheight = 0.9)
)

dev.off()

cat(sprintf("Combined Figure 5 saved: %s\n",
            file.path(FIG_DIR, "Figure5_Combined.png")))
cat("═══ Done ═══\n")
