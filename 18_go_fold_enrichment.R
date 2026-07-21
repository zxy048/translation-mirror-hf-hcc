# =============================================================================
# Script 18: GO Fold Enrichment — Supplementary Table S5
# Computes GO BP enrichment with fold enrichment for green module genes
# =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# Load WGCNA results
wgcna <- readRDS(file.path(PROJ_DIR, "GSE57338_WGCNA_rerun.rds"))

# Green module genes
green_idx <- which(wgcna$moduleColors == wgcna$primary_translation_module)
green_genes <- wgcna$gene_symbols[green_idx]
cat(sprintf("Green module: %d genes\n", length(green_genes)))

# Background: all genes in WGCNA
bg_genes <- wgcna$gene_symbols
cat(sprintf("Background: %d genes\n", length(bg_genes)))

# Run GO BP enrichment
ego <- enrichGO(
  gene          = green_genes,
  universe      = bg_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

cat(sprintf("GO terms enriched: %d\n", nrow(ego)))

# Compute fold enrichment: GeneRatio / BgRatio
# GeneRatio = k/n, BgRatio = M/N
# Fold enrichment = (k/n) / (M/N)
parse_ratio <- function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
}

gene_ratio <- sapply(ego@result$GeneRatio, parse_ratio)
bg_ratio   <- sapply(ego@result$BgRatio, parse_ratio)
fold_enrichment <- gene_ratio / bg_ratio

# Build output table
tbl_go <- data.frame(
  GO_ID           = ego@result$ID,
  Description     = ego@result$Description,
  GeneRatio       = ego@result$GeneRatio,
  BgRatio         = ego@result$BgRatio,
  Fold_Enrichment = round(fold_enrichment, 2),
  pvalue          = formatC(ego@result$pvalue, format = "e", digits = 2),
  p_adjust        = formatC(ego@result$p.adjust, format = "e", digits = 2),
  Count           = ego@result$Count,
  row.names       = NULL
)

# Sort by fold enrichment descending (for translation-related: high FE is expected)
tbl_go <- tbl_go[order(-tbl_go$Fold_Enrichment), ]

# Mark translation-related terms
trans_keywords <- "translation|ribosom|peptide|ribonucleoprotein|rRNA|translational"
tbl_go$Translation_Related <- ifelse(
  grepl(trans_keywords, tbl_go$Description, ignore.case = TRUE), "Yes", ""
)

cat(sprintf("\nTranslation-related terms: %d\n",
            sum(tbl_go$Translation_Related == "Yes")))

# Save
out_path <- file.path(PROJ_DIR, "Table_S5_GO_fold_enrichment.csv")
write.csv(tbl_go, out_path, row.names = FALSE)

cat(sprintf("\n=== Table S5 Preview (top 15 by fold enrichment) ===\n"))
print(head(tbl_go[, c("Description", "Fold_Enrichment", "p_adjust", "Count", "Translation_Related")], 15))

cat(sprintf("\nTable S5 saved: %s\n", out_path))
