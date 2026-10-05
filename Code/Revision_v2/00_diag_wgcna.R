# 00_diag_wgcna.R -- why is the WGCNA network degenerate?
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

d <- readRDS(file.path(OUT, "intermediate", "R01_universe.rds"))
U <- d$U
hf_dis <- d$ann_hf$gsm[d$ann_hf$group == "HF"]

# WGCNA orientation: rows = SAMPLES, columns = GENES. R01 stores genes x samples.
datHF  <- t(d$E_hf_full[U,  hf_dis,  drop = FALSE])
datHCC <- t(d$E_hcc_full[U, colnames(d$E_hcc_full), drop = FALSE])

cat("=== matrix sanity ===\n")
for (nm in c("datHF", "datHCC")) {
  m <- get(nm)
  cat(sprintf("%s: %d x %d  class=%s\n", nm, nrow(m), ncol(m), class(m)[1]))
  cat("   NA:", sum(is.na(m)), " NaN:", sum(is.nan(m)), " Inf:", sum(is.infinite(m)), "\n")
  v <- apply(m, 1, var, na.rm = TRUE)
  cat("   zero-variance rows:", sum(v == 0 | is.na(v)), "\n")
  cat("   var range:", sprintf("%.3e .. %.3e", min(v, na.rm = TRUE), max(v, na.rm = TRUE)), "\n")
  cat("   range:", sprintf("%.3f .. %.3f", min(m), max(m)), "\n")
}

cat("\n=== raw pickSoftThreshold output (HF, NO blockSize) ===\n")
disableWGCNAThreads()
sft_hf <- pickSoftThreshold(datHF, powerVector = WGCNA$powerVector,
                            networkType = "signed", verbose = 2)
cat("\ncolnames(fitIndices):", paste(colnames(sft_hf$fitIndices), collapse = " | "), "\n\n")
print(sft_hf$fitIndices)

cat("\n=== same for HCC ===\n")
sft_hcc <- pickSoftThreshold(datHCC, powerVector = WGCNA$powerVector,
                             networkType = "signed", verbose = 2)
print(sft_hcc$fitIndices)

cat("\n=== interpret columns directly ===\n")
for (nm in c("HF", "HCC")) {
  fi <- if (nm == "HF") sft_hf$fitIndices else sft_hcc$fitIndices
  cat("\n--", nm, "--\n")
  cat("col2 as-is          :", sprintf("%.3f", fi[, 2]), "\n")
  cat("sign(col3)          :", sprintf("%+d", sign(fi[, 3])), "\n")
  cat("-sign(col3)*col2    :", sprintf("%.3f", -sign(fi[, 3]) * fi[, 2]), "\n")
}
