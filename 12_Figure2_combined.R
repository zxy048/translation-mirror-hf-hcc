# =============================================================================
# 脚本 12：Figure 2 合并版（Phase 1 — HF WGCNA 重跑后）
# Panel A: HF dendrogram + module color bar (green highlighted)
# Panel B: Module-trait correlation bar plot
# Panel C: Green module GO enrichment bubble
# Panel D: Cross-disease module gene overlap (HF green vs HCC modules)
# =============================================================================

library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

PROJ_DIR <- "D:/R_projects/revision_analysis"
FIG_DIR  <- file.path(PROJ_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load rerun results ──
wgcna <- readRDS(file.path(PROJ_DIR, "GSE57338_WGCNA_rerun.rds"))
mod_colors <- wgcna$moduleColors
primary_mod <- wgcna$primary_translation_module  # "green"

# ── Color palette for modules ──
all_mods <- unique(mod_colors)
mod_palette <- setNames(all_mods, all_mods)
# Highlight translation module
highlight_colors <- rep("grey80", length(all_mods))
names(highlight_colors) <- all_mods
highlight_colors[primary_mod] <- "#228B22"  # forest green for emphasis

cat(sprintf("Primary translation module: %s\n", primary_mod))
cat(sprintf("Module genes: %d\n", sum(mod_colors == primary_mod)))

# =============================================================================
# Panel A: Dendrogram with module colors
# =============================================================================
cat("\n--- Panel A: Dendrogram ---\n")

geneTree <- wgcna$geneTree

png(file.path(FIG_DIR, "Figure2A_dendrogram.png"),
    width = 10, height = 3.5, units = "in", res = 300)
par(mar = c(2, 4, 1.5, 1))
WGCNA::plotDendroAndColors(
  geneTree,
  cbind(mod_colors),
  c("Module"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  guideHang = 0.05,
  main = sprintf("HF WGCNA Dendrogram — Translation module: %s (%d genes)",
                 primary_mod, sum(mod_colors == primary_mod)),
  cex.main = 1.0,
  cex.colorLabels = 0.8
)
dev.off()
cat("Panel A saved.\n")

# =============================================================================
# Panel B: Module-Trait Correlation
# =============================================================================
cat("\n--- Panel B: Module-Trait ---\n")

mt <- wgcna$module_trait_cor
mt <- mt[order(abs(mt$Correlation), decreasing = TRUE), ]
mt$Module <- factor(mt$Module, levels = rev(mt$Module))
mt$is_translation <- mt$Module == primary_mod

p_b <- ggplot(mt, aes(x = Correlation, y = Module)) +
  geom_bar(aes(fill = is_translation), stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#228B22", "FALSE" = "grey70"), guide = "none") +
  geom_text(aes(label = sprintf("r=%.3f", Correlation),
                hjust = ifelse(Correlation > 0, -0.15, 1.15)),
            size = 2.5, color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(
    title = "Module–Trait (HF vs NF) Correlation",
    subtitle = sprintf("Signed WGCNA, %s=%d, R%s=%.3f | %s (translation): r=%.3f, p=%.2e",
                       expression(beta), wgcna$beta_sel,
                       expression(""^2), wgcna$beta_r2,
                       primary_mod,
                       mt$Correlation[mt$Module == primary_mod],
                       mt$P_value[mt$Module == primary_mod]),
    x = "Pearson Correlation with HF Status",
    y = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8)
  )

ggsave(file.path(FIG_DIR, "Figure2B_module_trait.png"),
       p_b, width = 7, height = 4.5, dpi = 300)
cat("Panel B saved.\n")

# =============================================================================
# Panel C: GO Enrichment (green module)
# =============================================================================
cat("\n--- Panel C: GO ---\n")

ego <- wgcna$module_go_results[[primary_mod]]
if (!is.null(ego)) {
  go_df <- as.data.frame(ego)
  go_df <- go_df[order(go_df$p.adjust), ]
  # Top 15 terms
  go_top <- head(go_df, 15)
  go_top$Description <- factor(go_top$Description,
                               levels = rev(go_top$Description))
  go_top$neg_log10_padj <- -log10(go_top$p.adjust)
  go_top$Count <- as.numeric(sub("/.*", "", go_top$GeneRatio))

  # Classify as translation-related or not
  trans_kw <- c("ribosom", "translat", "rRNA", "peptide", "ribonucleoprotein", "translational")
  go_top$is_trans <- grepl(paste(trans_kw, collapse = "|"),
                           go_top$Description, ignore.case = TRUE)

  p_c <- ggplot(go_top, aes(x = neg_log10_padj, y = Description)) +
    geom_point(aes(size = Count, color = is_trans), alpha = 0.85) +
    scale_color_manual(values = c("TRUE" = "#228B22", "FALSE" = "grey50"), guide = "none") +
    scale_size_continuous(range = c(2, 7), name = "Gene count") +
    labs(
      title = sprintf("GO: %s Module (%d genes)", primary_mod,
                      sum(mod_colors == primary_mod)),
      subtitle = "Top 15 GO Biological Process terms",
      x = expression(-log[10](p.adjust)),
      y = ""
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8)
    )
  ggsave(file.path(FIG_DIR, "Figure2C_GO.png"),
         p_c, width = 8, height = 5, dpi = 300)
  cat("Panel C saved.\n")
} else {
  cat("Panel C skipped: no GO results.\n")
}

# =============================================================================
# Panel D: Cross-disease Module Overlap
# =============================================================================
cat("\n--- Panel D: Cross-Disease ---\n")

cd <- wgcna$cross_disease
if (!is.null(cd)) {
  hcc_dist <- cd$hcc_module_distribution
  if (length(hcc_dist) > 0) {
    overlap_df <- data.frame(
      HCC_Module = names(hcc_dist),
      Genes = as.numeric(hcc_dist),
      stringsAsFactors = FALSE
    )
    overlap_df <- overlap_df[order(overlap_df$Genes, decreasing = TRUE), ]
    overlap_df$HCC_Module <- factor(overlap_df$HCC_Module,
                                    levels = rev(overlap_df$HCC_Module))

    # Flag HCC blue module (original HCC translation module)
    overlap_df$is_hcc_trans <- overlap_df$HCC_Module == "blue"

    # Add Fisher test annotation
    fisher_or <- cd$fisher_test$odds_ratio
    fisher_p <- cd$fisher_test$p_value

    p_d <- ggplot(overlap_df, aes(x = Genes, y = HCC_Module)) +
      geom_bar(aes(fill = is_hcc_trans), stat = "identity", width = 0.65) +
      scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "grey70"), guide = "none") +
      geom_text(aes(label = Genes, hjust = -0.2), size = 3, color = "grey30") +
      labs(
        title = sprintf("HF %s Genes in HCC WGCNA Modules", primary_mod),
        subtitle = sprintf("Fisher OR=%.1f, p=%.3f | %d/%d HF %s genes found in HCC WGCNA",
                           fisher_or, fisher_p,
                           length(cd$hf_trans_in_hcc),
                           sum(mod_colors == primary_mod),
                           primary_mod),
        x = "Number of overlapping genes",
        y = "HCC Module"
      ) +
      xlim(0, max(overlap_df$Genes) * 1.25) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8)
      )
    ggsave(file.path(FIG_DIR, "Figure2D_overlap.png"),
           p_d, width = 7, height = 5, dpi = 300)
    cat("Panel D saved.\n")
  }
}

# =============================================================================
# Combined Figure 2 (2×2 layout)
# =============================================================================
cat("\n--- Combining Figure 2 ---\n")

# Read panels back as raster grobs
library(png)
library(grid)

read_png_grob <- function(path) {
  if (!file.exists(path)) return(nullGrob())
  img <- readPNG(path)
  rasterGrob(img, interpolate = TRUE)
}

gA <- read_png_grob(file.path(FIG_DIR, "Figure2A_dendrogram.png"))
gB <- read_png_grob(file.path(FIG_DIR, "Figure2B_module_trait.png"))
gC <- read_png_grob(file.path(FIG_DIR, "Figure2C_GO.png"))
gD <- read_png_grob(file.path(FIG_DIR, "Figure2D_overlap.png"))

# Combined layout
png(file.path(FIG_DIR, "Figure2_Translation_Module_Identification.png"),
    width = 16, height = 12, units = "in", res = 300)

grid.newpage()
# Title
grid.text(
  "Figure 2. Identification of disease-context-dependent translation-related modules in HF and HCC.",
  x = 0.02, y = 0.99, just = c(0, 1),
  gp = gpar(fontface = "bold", fontsize = 14)
)

# Layout: 2 rows x 2 columns
# Row 1: A (dendrogram, wide) | B (module-trait)
# Row 2: C (GO)               | D (cross-disease)
pushViewport(viewport(layout = grid.layout(2, 2,
  heights = c(0.45, 0.55),
  widths  = c(0.55, 0.45))))

# Labels
label_gp <- gpar(fontface = "bold", fontsize = 14, col = "black")

# A: dendrogram
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("A", x = 0.02, y = 0.97, just = c(0, 1), gp = label_gp)
grid.draw(gA)
popViewport()

# B: module-trait
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.text("B", x = 0.02, y = 0.97, just = c(0, 1), gp = label_gp)
grid.draw(gB)
popViewport()

# C: GO
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.text("C", x = 0.02, y = 0.97, just = c(0, 1), gp = label_gp)
grid.draw(gC)
popViewport()

# D: cross-disease
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.text("D", x = 0.02, y = 0.97, just = c(0, 1), gp = label_gp)
grid.draw(gD)
popViewport()

# Legend at bottom
grid.text(
  sprintf("(A) HF WGCNA dendrogram. Green = translation-related module (%s, %d genes). (B) Module–trait (HF vs NF) correlations. (C) GO BP enrichment of the %s module. (D) Overlap of HF %s module genes with HCC WGCNA modules (GSE141198).",
          primary_mod, sum(mod_colors == primary_mod), primary_mod, primary_mod),
  x = 0.02, y = 0.005, just = c(0, 0),
  gp = gpar(fontsize = 8, col = "grey40")
)

dev.off()
cat(sprintf("Combined Figure 2 saved: %s\n",
            file.path(FIG_DIR, "Figure2_Translation_Module_Identification.png")))

cat("\n=== Figure 2 generation complete ===\n")
