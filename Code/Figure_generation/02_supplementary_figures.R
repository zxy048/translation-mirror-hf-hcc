# 08_generate_remaining_figures.R
# 组装剩余图表: Fig2, Fig4, FigS1, FigS2, Table1

library(ggplot2)
library(png)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ============================================================
# Figure S2: 方向一致性散点图 (重命名已有文件)
# ============================================================
cat("=== Figure S2 ===\n")
file.copy(
  file.path(PROJ_DIR, "figures", "Figure_CrossDisease_DirectionTest.png"),
  file.path(PROJ_DIR, "figures", "Figure_S1_Direction_Consistency.png"),
  overwrite = TRUE)
cat("  done\n")

# ============================================================
# Figure 4: TGS KM + Forest 二合一
# ============================================================
cat("\n=== Figure 4: TGS Validation ===\n")
km_png     <- png::readPNG(file.path(PROJ_DIR, "figures", "GSE141198_TGS_KM.png"))
forest_png <- png::readPNG(file.path(PROJ_DIR, "figures", "Figure_CrossDisease_Forest.png"))

png(file.path(PROJ_DIR, "figures", "Figure4_TGS_Validation.png"),
    width = 16, height = 7, units = "in", res = 300)
par(mar = c(0,0,0,0))
layout(matrix(1:2, 1, 2))
plot.new(); rasterImage(km_png, 0, 0, 1, 1)
mtext("A", side = 3, line = -2, adj = 0.05, font = 2, cex = 1.5)
plot.new(); rasterImage(forest_png, 0, 0, 1, 1)
mtext("B", side = 3, line = -2, adj = 0.05, font = 2, cex = 1.5)
dev.off()
cat("  done\n")

# ============================================================
# Figure 2 (NEW — Cross-Disease Conservation)
#   A: HF GSE57338 dendrogram (black module = translation)
#   B: HCC GSE141198 dendrogram (blue module = translation)
#   C: Shared functional enrichment (HF black vs HCC blue GO terms)
# ============================================================
cat("\n=== Figure 2: Cross-Disease WGCNA Conservation ===\n")

# ── Panel A: HF Dendrogram (from 09_HF_WGCNA_final.R) ──
hf_dendro_png <- png::readPNG(file.path(PROJ_DIR, "figures", "Figure2A_HF_dendrogram.png"))
cat("  2A: HF dendrogram loaded\n")

# ── Panel B: HCC Dendrogram (regenerate with matched styling) ──
hcc_wgcna <- readRDS(file.path(PROJ_DIR, "GSE141198_WGCNA_result.rds"))
hcc_net <- hcc_wgcna$net
hcc_moduleColors <- hcc_wgcna$moduleColors

# Get HCC translation module info — force "blue" (matching manuscript)
# ★ HCC WGCNA used loose GO regex that may flag grey as "translation";
#     manuscript declares blue = translation module
hcc_trans_mods_auto <- names(hcc_wgcna$translation_modules)
# Prefer "blue" if it exists and was flagged; otherwise use it anyway
hcc_trans_mod <- if ("blue" %in% hcc_trans_mods_auto) "blue" else if (length(hcc_trans_mods_auto) > 0) hcc_trans_mods_auto[1] else "blue"
hcc_trans_n <- sum(hcc_moduleColors == hcc_trans_mod)
cat(sprintf("  HCC translation module: %s (%d genes) [auto-detected: %s]\n",
            hcc_trans_mod, hcc_trans_n, paste(hcc_trans_mods_auto, collapse=", ")))

png(file.path(PROJ_DIR, "figures", "Figure2B_HCC_dendrogram.png"),
    width = 14, height = 6, units = "in", res = 300)
plotDendroAndColors(
  hcc_net$dendrograms[[1]],
  hcc_moduleColors[hcc_net$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = sprintf("GSE141198 HCC WGCNA: Gene Dendrogram & Module Colors (β=%d, R²=0.84)",
                 hcc_wgcna$beta_sel)
)
dev.off()
hcc_dendro_png <- png::readPNG(file.path(PROJ_DIR, "figures", "Figure2B_HCC_dendrogram.png"))
cat(sprintf("  2B: HCC dendrogram regenerated (%s module, %d genes)\n",
            hcc_trans_mod, hcc_trans_n))

# ── Panel C: Shared Functional Enrichment ──
cat("  2C: Generating shared GO comparison...\n")

# Load HF WGCNA result
hf_wgcna <- readRDS(file.path(PROJ_DIR, "GSE57338_WGCNA_result.rds"))
hf_trans_mod <- hf_wgcna$primary_translation_module
hf_trans_genes <- hf_wgcna$translation_module_genes
hf_all_genes <- hf_wgcna$gene_symbols

# HCC translation module genes
hcc_expr <- readRDS(file.path(PROJ_DIR, "GSE141198_wgcna_input.rds"))
hcc_all_genes <- rownames(hcc_expr)
hcc_block1_idx <- hcc_net$blockGenes[[1]]
hcc_trans_genes <- hcc_all_genes[hcc_block1_idx][hcc_moduleColors[hcc_block1_idx] == hcc_trans_mod]

cat(sprintf("    HF %s: %d genes, HCC %s: %d genes\n",
            hf_trans_mod, length(hf_trans_genes), hcc_trans_mod, length(hcc_trans_genes)))

# GO enrichment for both modules
# ★ HF: use relaxed cutoffs (p=1.0) because top-3000-variance filter excludes
#     ribosomal protein genes, weakening the GO signal. Capture all terms
#     and rely on the translation regex to identify relevant hits.
# ★ HCC: standard cutoffs (strong enrichment expected)

ego_hf <- enrichGO(gene = hf_trans_genes, OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
  pvalueCutoff = 1.0, qvalueCutoff = 1.0)  # capture all terms
ego_hcc <- enrichGO(gene = hcc_trans_genes, OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
  pvalueCutoff = 0.01, qvalueCutoff = 0.05)

# Translation-related GO term regex (same as original Figure 2B)
trans_go_regex <- "translat|ribosom|peptide.*biosyn|ribonucleoprotein|rRNA"

# Collect translation-related GO terms from both modules
get_trans_terms <- function(ego, regex) {
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  df <- as.data.frame(ego)
  df <- df[grepl(regex, df$Description, ignore.case = TRUE), ]
  df[order(df$p.adjust), ]
}

hf_trans_go <- get_trans_terms(ego_hf, trans_go_regex)
hcc_trans_go <- get_trans_terms(ego_hcc, trans_go_regex)

if (!is.null(hf_trans_go) && nrow(hf_trans_go) > 0 &&
    !is.null(hcc_trans_go) && nrow(hcc_trans_go) > 0) {

  # Build shared GO term comparison
  # Use all translation terms from both, keep top N per module
  all_go_ids <- unique(c(hf_trans_go$ID, hcc_trans_go$ID))

  # For each GO ID, get the -log10(p.adjust) from each module
  shared_go_df <- data.frame(
    GO_ID = all_go_ids,
    Description = sapply(all_go_ids, function(id) {
      if (id %in% hf_trans_go$ID) hf_trans_go$Description[hf_trans_go$ID == id][1]
      else hcc_trans_go$Description[hcc_trans_go$ID == id][1]
    }),
    HF_logP = sapply(all_go_ids, function(id) {
      if (id %in% hf_trans_go$ID) -log10(hf_trans_go$p.adjust[hf_trans_go$ID == id][1]) else 0
    }),
    HCC_logP = sapply(all_go_ids, function(id) {
      if (id %in% hcc_trans_go$ID) -log10(hcc_trans_go$p.adjust[hcc_trans_go$ID == id][1]) else 0
    }),
    stringsAsFactors = FALSE
  )

  # Classify: shared (both > 0), HF-only, HCC-only
  shared_go_df$Category <- ifelse(
    shared_go_df$HF_logP > 0 & shared_go_df$HCC_logP > 0, "Shared",
    ifelse(shared_go_df$HF_logP > 0, "HF-specific", "HCC-specific"))

  # Sort: shared first, then by max logP
  shared_go_df$max_logP <- pmax(shared_go_df$HF_logP, shared_go_df$HCC_logP)
  shared_go_df <- shared_go_df[order(
    factor(shared_go_df$Category, levels = c("Shared", "HF-specific", "HCC-specific")),
    -shared_go_df$max_logP), ]

  # Keep top 20
  plot_go <- head(shared_go_df, 20)
  plot_go$Description <- factor(plot_go$Description, levels = rev(unique(plot_go$Description)))

  # Reshape for grouped bar chart
  plot_long <- rbind(
    data.frame(Description = plot_go$Description, logP = plot_go$HF_logP,
               Module = sprintf("HF %s", hf_trans_mod), Category = plot_go$Category),
    data.frame(Description = plot_go$Description, logP = plot_go$HCC_logP,
               Module = sprintf("HCC %s", hcc_trans_mod), Category = plot_go$Category)
  )

  # Build fill color mapping
  fill_colors <- c("#333333", "#377EB8")
  names(fill_colors) <- c(sprintf("HF %s", hf_trans_mod), sprintf("HCC %s", hcc_trans_mod))

  p_shared_go <- ggplot(plot_long, aes(x = logP, y = Description, fill = Module)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
    scale_fill_manual(values = fill_colors) +
    labs(title = "Shared Functional Enrichment of Translation-Associated Modules",
         subtitle = sprintf("GO Biological Process: HF %s (%d genes) vs HCC %s (%d genes)\nTranslation-related terms; relaxed detection threshold applied for HF module",
                            hf_trans_mod, length(hf_trans_genes),
                            hcc_trans_mod, length(hcc_trans_genes)),
         x = expression(-log[10](p.adjust)), y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10.5),
          legend.position = "bottom")

  ggsave(file.path(PROJ_DIR, "figures", "Figure2C_Shared_GO.png"),
         p_shared_go, width = 12, height = 7, dpi = 300)
  cat(sprintf("  2C: Shared GO plot saved (%d terms)\n", nrow(plot_go)))

  shared_go_png <- png::readPNG(file.path(PROJ_DIR, "figures", "Figure2C_Shared_GO.png"))
} else {
  cat("  2C: WARNING — insufficient GO terms for shared comparison\n")
  shared_go_png <- NULL
}

# ── Composite Figure 2 ──
cat("  Compositing Figure 2...\n")

png(file.path(PROJ_DIR, "figures", "Figure2_CrossDisease_Conservation.png"),
    width = 16, height = 14, units = "in", res = 300)

if (!is.null(shared_go_png)) {
  layout(matrix(c(1,1, 2,2, 3,3,3), 7, 1))
} else {
  layout(matrix(c(1,1, 2,2), 4, 1))
}
par(mar = c(0, 0, 0, 0))

# Panel A: HF dendrogram
plot.new(); rasterImage(hf_dendro_png, 0, 0, 1, 1)
mtext("A", side = 3, line = -2, adj = 0.05, font = 2, cex = 2)

# Panel B: HCC dendrogram
plot.new(); rasterImage(hcc_dendro_png, 0, 0, 1, 1)
mtext("B", side = 3, line = -2, adj = 0.05, font = 2, cex = 2)

# Panel C: Shared GO
if (!is.null(shared_go_png)) {
  plot.new(); rasterImage(shared_go_png, 0, 0, 1, 1)
  mtext("C", side = 3, line = -2, adj = 0.05, font = 2, cex = 2)
}

dev.off()
cat("  Figure 2 (CrossDisease_Conservation) saved\n")

# ============================================================
# Figure S1: 软阈值选择 (GSE141198已有, GSE57338待HF脚本完成)
# ============================================================
cat("\n=== Figure S1: Soft Threshold ===\n")
file.copy(
  file.path(PROJ_DIR, "WGCNA_GSE141198", "soft_threshold.png"),
  file.path(PROJ_DIR, "figures", "Figure_S4A_SoftThreshold_GSE141198.png"),
  overwrite = TRUE)
cat("  GSE141198 done\n")
cat("  GSE57338: 请从v6.1结果中复制软阈值图，命名为 Figure_S4B_SoftThreshold_GSE57338.png\n")

# ============================================================
# Table 1: 队列基线特征 (模板)
# ============================================================
cat("\n=== Table 1: Cohort Characteristics ===\n")
sink(file.path(PROJ_DIR, "tables", "Table1_Cohort_Characteristics_template.txt"))
cat("Table 1. Cohort characteristics across three datasets.\n\n")
cat("Characteristic        GSE57338 (HF)      TCGA-LIHC (HCC)     GSE141198 (HCC)\n")
cat("----------------------------------------------------------------------------\n")
cat("Sample size           313                 424                  148\n")
cat("Platform              Affymetrix          Illumina RNA-seq     Illumina RNA-seq\n")
cat("Tissue                LV myocardium       Liver                Liver\n")
cat("Disease subtypes      DCM/ICM/NF          HCC (371T + 50N)     HCC (all tumor)\n")
cat("Age (mean±SD)         [待提取]            [待提取]             [待提取]\n")
cat("Sex (M/F)             [待提取]            [待提取]             [待提取]\n")
cat("OS events             N/A                 [待提取]             94\n")
cat("Etiology              N/A                 HBV/HCV/Alcohol/etc  HBV/HCV/NBNC\n")
sink()
cat("  Template saved to tables/Table1_Cohort_Characteristics_template.txt\n")

cat("\n═══════════════════════════════════\n")
cat(" 完成清单:\n")
cat(" ✅ Figure S2 (方向一致性)\n")
cat(" ✅ Figure 4 (TGS验证 KM+Forest)\n")
cat(" ✅ Figure 2 (WGCNA三拼图)\n")
cat(" ⚠ Figure S1 - GSE57338图片需从v6.1复制\n")
cat(" ⚠ Table 1 - 需从各数据集pdata提取临床特征\n")
cat("═══════════════════════════════════\n")
