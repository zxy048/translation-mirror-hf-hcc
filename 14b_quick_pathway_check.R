# Quick script to get MYC V1/V2 pathway-TATS correlations
# Avoids S4 object loading issues by using only simple RDS

library(msigdbr)

# These should be safe to load (matrices/data.frames, no S4)
expr_vst <- readRDS("D:/R_projects/revision_analysis/GSE141198_vst.rds")
cat("expr_vst loaded:", nrow(expr_vst), "x", ncol(expr_vst), "\n")

# Get green genes from the WGCNA RDS (this has S4 objects but might load)
wgcna <- readRDS("D:/R_projects/revision_analysis/GSE57338_WGCNA_rerun.rds")
primary_mod <- as.character(wgcna$primary_translation_module)
green_genes <- wgcna$gene_symbols[wgcna$moduleColors == primary_mod]
cat("Green genes:", length(green_genes), "\n")

# Compute TATS
green_in <- intersect(green_genes, rownames(expr_vst))
expr_green <- expr_vst[green_in, , drop = FALSE]
expr_z <- t(scale(t(expr_green)))
tats <- colMeans(expr_z, na.rm = TRUE)
cat("TATS computed:", length(tats), "samples\n")

# Get MYC target genes from MSigDB
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")

# Compute simple pathway score (mean z-score of pathway genes)
compute_pw_score <- function(pw_name) {
  pw_genes_all <- msig_h$gene_symbol[msig_h$gs_name == pw_name]
  pw_genes <- intersect(pw_genes_all, rownames(expr_vst))
  pw_expr <- expr_vst[pw_genes, , drop = FALSE]
  pw_z <- t(scale(t(pw_expr)))
  pw_score <- colMeans(pw_z, na.rm = TRUE)
  return(list(score = pw_score, n_genes = length(pw_genes)))
}

for (pw_name in c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
                   "HALLMARK_E2F_TARGETS", "HALLMARK_MTORC1_SIGNALING")) {
  res <- compute_pw_score(pw_name)
  common <- intersect(names(tats), names(res$score))
  ct <- cor.test(res$score[common], tats[common], method = "spearman")
  cat(sprintf("%-40s rho=%+.4f  p=%.6f  genes=%d\n",
              pw_name, ct$estimate, ct$p.value, res$n_genes))
}
