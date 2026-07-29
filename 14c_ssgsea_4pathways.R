# Quick 4-pathway ssGSEA check
library(GSVA)
library(msigdbr)

expr_vst <- readRDS("D:/R_projects/revision_analysis/GSE141198_vst.rds")
wgcna <- readRDS("D:/R_projects/revision_analysis/GSE57338_WGCNA_rerun.rds")
primary_mod <- as.character(wgcna$primary_translation_module)
green_genes <- wgcna$gene_symbols[wgcna$moduleColors == primary_mod]
green_in <- intersect(green_genes, rownames(expr_vst))
expr_green <- expr_vst[green_in, , drop = FALSE]
expr_z <- t(scale(t(expr_green)))
tats <- colMeans(expr_z, na.rm = TRUE)
cat(sprintf("TATS: %d samples\n", length(tats)))

msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
gs <- list()
target_pws <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
                 "HALLMARK_E2F_TARGETS", "HALLMARK_MTORC1_SIGNALING")
for (pn in target_pws) {
  gs[[pn]] <- msig_h$gene_symbol[msig_h$gs_name == pn]
  cat(sprintf("%s: %d genes\n", pn, length(gs[[pn]])))
}

param <- ssgseaParam(as.matrix(expr_vst), gs, minSize = 5, maxSize = 500)
scores <- gsva(param, verbose = FALSE)
cat(sprintf("ssGSEA: %d pathways x %d samples\n", nrow(scores), ncol(scores)))

common <- intersect(colnames(scores), names(tats))
cat(sprintf("Common samples: %d\n", length(common)))

for (pn in target_pws) {
  ct <- cor.test(as.numeric(scores[pn, common]), tats[common], method = "spearman")
  cat(sprintf("%-40s rho=%+.4f  p=%.6f\n", pn, ct$estimate, ct$p.value))
}
