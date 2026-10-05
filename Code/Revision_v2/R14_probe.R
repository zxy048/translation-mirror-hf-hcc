source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(matrixStats))

r2 <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))
r3 <- readRDS(file.path(OUT, "intermediate", "R03_translation_module.rds"))
cat("R02 names:", paste(names(r2), collapse = ", "), "\n")
cat("R03 names:", paste(names(r3), collapse = ", "), "\n")
cat("r2$U length:", length(r2$U), "\n")
for (t in c("HF", "HCC")) {
  cat("--", t, "--\n")
  cat("  r2[[t]] names:", paste(names(r2[[t]]), collapse = ", "), "\n")
  cat("  power:", r2[[t]]$power, "\n")
  cat("  net names:", paste(names(r2[[t]]$net), collapse = ", "), "\n")
  cat("  n colors:", length(r2[[t]]$net$colors), "\n")
  dt <- if (t == "HF") r2$datHF else r2$datHCC
  cat("  dat dim:", paste(dim(dt), collapse = " x "), "\n")
  gx <- if (t == "HF") r2$gxHF else r2$gxHCC
  cat("  gx dim:", paste(dim(gx), collapse = " x "), "\n")
  dk <- if (t == "HF") "desHF" else "desHCC"
  cat("  r3[[", dk, "]] names:", paste(names(r3[[dk]]), collapse = ", "), "\n")
  cat("  module:", r3[[dk]]$module, "\n")
}

adj_from_cor <- function(r, power) {
  A <- (0.5 * (1 + r))^power
  A[A < 0] <- 0; A[A > 1] <- 1
  A
}
standardise <- function(dat) {
  Z <- t(dat); Z <- Z - rowMeans(Z)
  s <- rowSds(Z); s[s == 0] <- 1; Z / s
}
for (t in c("HF", "HCC")) {
  dt <- if (t == "HF") r2$datHF else r2$datHCC
  pw <- r2[[t]]$power
  U  <- r2$U
  n  <- nrow(dt)
  set.seed(STATS$SEED)
  sub <- dt[, U, drop = FALSE][, sample(length(U), 1200), drop = FALSE]
  Zs  <- standardise(sub)
  A   <- adj_from_cor((Zs %*% t(Zs)) / (n - 1), pw)
  k   <- rowSums(A)
  ref  <- WGCNA::TOMsimilarity(A, TOMType = "signed")
  mine <- (A %*% t(A) - A) / (outer(k - 1, k - 1, pmin) + 1 - A)
  off  <- upper.tri(ref)
  cat(sprintf("%s: power=%d  off-diag max|diff|=%.3e  diag(ref)=%.4f diag(mine)=%.4f\n",
              t, pw, max(abs(mine[off] - ref[off])), ref[1, 1], mine[1, 1]))
}
cat("PROBE DONE\n")
