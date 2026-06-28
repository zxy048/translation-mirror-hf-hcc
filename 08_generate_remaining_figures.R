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
  file.path(PROJ_DIR, "figures", "Figure_S2_Direction_Consistency.png"),
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
# Figure 2: GSE141198 WGCNA 三拼图
#   A: 树状图+模块色带
#   B: GO富集气泡图 (blue模块top15)
#   C: 模块-性状热图
# ============================================================
cat("\n=== Figure 2: WGCNA ===\n")

wgcna_res <- readRDS(file.path(PROJ_DIR, "GSE141198_WGCNA_result.rds"))
net <- wgcna_res$net
moduleColors <- wgcna_res$moduleColors

# A: 树状图
png(file.path(PROJ_DIR, "figures", "Figure2A_dendrogram.png"),
    width = 12, height = 5, units = "in", res = 300)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "GSE141198 WGCNA: Gene Dendrogram & Module Colors"
)
dev.off()
cat("  2A done\n")

# B: GO富集气泡图
# ★ net$blockGenes[[1]] 返回的是整数索引，不是基因名
# 从WGCNA输入数据获取真正的基因SYMBOL
expr_wgcna <- readRDS(file.path(PROJ_DIR, "GSE141198_wgcna_input.rds"))
all_genes_input <- rownames(expr_wgcna)                  # 真正的SYMBOL (nrow个)
block1_idx <- net$blockGenes[[1]]                         # block1的基因索引
moduleColors_block1 <- moduleColors[block1_idx]           # 按block1索引取模块颜色

# 维度安全检查
if (length(block1_idx) != length(all_genes_input)) {
  cat(sprintf("  ⚠ 维度不匹配: block1=%d, input_genes=%d (goodSamplesGenes可能移除了基因)\n",
              length(block1_idx), length(all_genes_input)))
}
# 从block1索引对应的输入基因中取SYMBOL
all_genes <- all_genes_input[block1_idx]
blue_genes <- all_genes[moduleColors_block1 == "blue"]
cat(sprintf("  Blue module: %d genes\n", length(blue_genes)))

ego <- enrichGO(
  gene = blue_genes, OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL", ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05
)

if (!is.null(ego) && nrow(ego@result) > 0) {
  go_df <- as.data.frame(ego)
  go_df <- go_df[order(go_df$p.adjust), ]

  # 标记翻译相关条目
  go_df$is_translation <- grepl("translat|ribosom|peptide.*biosyn|ribonucleoprotein|rRNA",
                                go_df$Description, ignore.case = TRUE)

  plot_df <- head(go_df, 15)
  plot_df$Description <- factor(plot_df$Description, levels = rev(unique(plot_df$Description)))

  p_go <- ggplot(plot_df, aes(x = Count, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "steelblue", name = "p.adjust",
                          trans = "log10") +
    labs(title = "GO Enrichment: Blue Module (GSE141198)",
         subtitle = sprintf("%d genes; top 15 terms", length(blue_genes)),
         x = "Gene Count", y = "") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")

  ggsave(file.path(PROJ_DIR, "figures", "Figure2B_GO_Bubble.png"),
         p_go, width = 10, height = 5.5, dpi = 300)
  cat("  2B done\n")
} else {
  cat("  2B skipped (no GO results)\n")
}

# C: 组合三拼图
dendro_png  <- png::readPNG(file.path(PROJ_DIR, "figures", "Figure2A_dendrogram.png"))
go_png      <- if (file.exists(file.path(PROJ_DIR, "figures", "Figure2B_GO_Bubble.png")))
                   png::readPNG(file.path(PROJ_DIR, "figures", "Figure2B_GO_Bubble.png")) else NULL
heatmap_png <- png::readPNG(file.path(PROJ_DIR, "WGCNA_GSE141198", "module_trait_heatmap.png"))

png(file.path(PROJ_DIR, "figures", "Figure2_WGCNA_GSE141198.png"),
    width = 14, height = 15, units = "in", res = 300)
layout(matrix(c(1,1, 2,2, 3,3), 6, 1))
par(mar = c(0, 0, 0, 0))
plot.new(); rasterImage(dendro_png, 0, 0, 1, 1)
mtext("A", side = 3, line = -2, adj = 0.05, font = 2, cex = 2)
if (!is.null(go_png)) {
  plot.new(); rasterImage(go_png, 0, 0, 1, 1)
  mtext("B", side = 3, line = -2, adj = 0.05, font = 2, cex = 2)
}
plot.new(); rasterImage(heatmap_png, 0, 0, 1, 1)
mtext("C", side = 3, line = -2, adj = 0.05, font = 2, cex = 2)
dev.off()
cat("  Figure 2 combined done\n")

# ============================================================
# Figure S1: 软阈值选择 (GSE141198已有, GSE57338待HF脚本完成)
# ============================================================
cat("\n=== Figure S1: Soft Threshold ===\n")
file.copy(
  file.path(PROJ_DIR, "WGCNA_GSE141198", "soft_threshold.png"),
  file.path(PROJ_DIR, "figures", "Figure_S1_SoftThreshold_GSE141198.png"),
  overwrite = TRUE)
cat("  GSE141198 done\n")
cat("  GSE57338: 请从v6.1结果中复制软阈值图，命名为 Figure_S1_SoftThreshold_GSE57338.png\n")

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
