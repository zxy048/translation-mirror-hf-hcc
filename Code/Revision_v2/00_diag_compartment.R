# 00_diag_compartment.R -- is the module mismatch a COMPARTMENT effect?
# HF designated module looks like the cytoplasmic ribosome (RPS/RPL).
# HCC designated module looks like OXPHOS + mitoribosome (COX/ATP5/NDUF + MRPL/MRPS).
# Test whether Tier A genes partition by compartment, and whether the HF platform
# is simply blind to the mitochondrial side (R01 lost every MT- gene).
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

r2 <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))
r3 <- readRDS(file.path(OUT, "intermediate", "R03_translation_module.rds"))
A  <- r3$tierA_genes

setsA <- resolve_sets(SETS_TIER_A)
mito_sets <- grep("MITOCHONDRI", names(setsA), value = TRUE)
cyto_sets <- setdiff(names(setsA), mito_sets)
cat("mito Tier A sets  (", length(mito_sets), "):", paste(mito_sets, collapse = "\n                      "), "\n")
cat("cyto Tier A sets  (", length(cyto_sets), "):", paste(cyto_sets, collapse = ", "), "\n")

mito_g <- unique(unlist(setsA[mito_sets]))
cyto_g <- setdiff(unique(unlist(setsA[cyto_sets])), mito_g)
cat(sprintf("\nTier A genes: total %d | mito-only %d | cyto-ish %d\n",
            length(A), length(mito_g), length(cyto_g)))

labHF  <- r2$HF$net$colors
labHCC <- r2$HCC$net$colors
mHF <- r3$desHF$module; mHCC <- r3$desHCC$module

cat(sprintf("\n=== compartment composition of the two designated modules ===\n"))
for (nm in c("HF", "HCC")) {
  mm <- if (nm == "HF") r3$trans_genes$HF else r3$trans_genes$HCC
  mt <- intersect(mm, mito_g); cy <- intersect(mm, cyto_g)
  cat(sprintf("%s module '%s' (n=%d): mito %d (%.0f%%) | cyto %d (%.0f%%) | neither %d\n",
              nm, if (nm == "HF") mHF else mHCC, length(mm),
              length(mt), 100 * length(mt) / length(mm),
              length(cy), 100 * length(cy) / length(mm),
              length(mm) - length(mt) - length(cy)))
}

cat("\n=== where does each compartment's gene set sit in each network? ===\n")
for (nm in c("HF", "HCC")) {
  lab <- if (nm == "HF") labHF else labHCC
  n_mod <- length(unique(lab))
  cat(sprintf("\n-- %s (%d genes in %d modules) --\n", nm, length(lab), n_mod))
  for (cmp in c("mito", "cyto")) {
    g  <- intersect(names(lab), if (cmp == "mito") mito_g else cyto_g)
    tb <- sort(table(lab[g]), decreasing = TRUE)
    top <- head(tb, 3)
    # what share of this compartment lands in the designated module?
    des <- if (nm == "HF") mHF else mHCC
    cat(sprintf("   %s genes in network: %3d | in designated module '%s': %3d (%.0f%%) | top modules: %s\n",
                cmp, length(g), des, sum(lab[g] == des, na.rm = TRUE),
                100 * sum(lab[g] == des, na.rm = TRUE) / length(g),
                paste(sprintf("%s=%d", names(top), as.integer(top)), collapse = " ")))
  }
}

cat("\n=== is the HF network simply missing mito genes? ===\n")
cat(sprintf("mito-compartment Tier A genes in U            : %d\n", length(intersect(r2$U, mito_g))))
cat(sprintf("  ... present in HF network                   : %d\n", length(intersect(names(labHF), mito_g))))
cat(sprintf("  ... present in HCC network                  : %d\n", length(intersect(names(labHCC), mito_g))))
mt_enc <- grep("^MT-", A, value = TRUE)
cat(sprintf("MT-encoded Tier A genes in U                  : %d  [%s]\n",
            length(mt_enc), paste(sort(mt_enc), collapse = ", ")))

cat("\n=== HCC magenta: which Tier A genes, by compartment ===\n")
mmh <- r3$trans_genes$HCC
cat("mito:", paste(sort(intersect(mmh, mito_g)), collapse = " "), "\n")
cat("cyto:", paste(sort(intersect(mmh, cyto_g)), collapse = " "), "\n")
cat("\n=== HF purple: which Tier A genes, by compartment ===\n")
mmf <- r3$trans_genes$HF
cat("mito:", paste(sort(intersect(mmf, mito_g)), collapse = " "), "\n")
cat("cyto:", paste(sort(intersect(mmf, cyto_g)), collapse = " "), "\n")
