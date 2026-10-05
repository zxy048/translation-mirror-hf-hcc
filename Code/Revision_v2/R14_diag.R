# Diagnostic: is my adjacency the one WGCNA used, and are these modules cohesive
# by a metric that does not depend on my TOM implementation?
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(matrixStats))

r2 <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))
r3 <- readRDS(file.path(OUT, "intermediate", "R03_translation_module.rds"))
U  <- r2$U

adj_from_cor <- function(r, power) {
  A <- (0.5 * (1 + r))^power
  A[A < 0] <- 0; A[A > 1] <- 1
  A
}

for (tag in c("HF", "HCC")) {
  dat   <- if (tag == "HF") r2$datHF else r2$datHCC
  power <- r2[[tag]]$power
  mod   <- r3[[if (tag == "HF") "desHF" else "desHCC"]]$module
  colors <- r2[[tag]]$net$colors
  if (is.null(names(colors))) names(colors) <- U
  genes <- names(colors)[colors == mod]
  cat("\n==================", tag, "module", mod, ":", length(genes), "genes ==================\n")

  # ---- 1. does my adjacency reproduce WGCNA::adjacency? ----
  set.seed(STATS$SEED)
  sub <- c(genes, sample(setdiff(U, genes), 2000 - length(genes)))
  ds  <- dat[, sub, drop = FALSE]                      # samples x genes
  Aw  <- WGCNA::adjacency(ds, power = power, type = "signed")
  Zs  <- t(ds); Zs <- Zs - rowMeans(Zs)
  s   <- rowSds(Zs); s[s == 0] <- 1; Zs <- Zs / s
  Am  <- adj_from_cor((Zs %*% t(Zs)) / (nrow(ds) - 1), power)
  cat(sprintf("  adjacency check vs WGCNA::adjacency: max|diff| = %.3e\n",
              max(abs(Am - Aw))))

  # ---- 2. simple, implementation-independent cohesion: mean pairwise r ----
  n <- nrow(dat)
  mean_r <- function(idx) {
    Z <- t(dat[, U[idx], drop = FALSE])
    Z <- Z - rowMeans(Z); s2 <- rowSds(Z); s2[s2 == 0] <- 1; Z <- Z / s2
    R <- (Z %*% t(Z)) / (n - 1)
    mean(R[upper.tri(R)])
  }
  gi <- match(genes, U)
  obs_r <- mean_r(gi)
  set.seed(STATS$SEED)
  null_r <- replicate(1000, mean_r(sample(length(U), length(gi))))
  cat(sprintf("  mean pairwise r: module %.4f | random %.4f (sd %.4f) | Z = %.2f\n",
              obs_r, mean(null_r), sd(null_r), (obs_r - mean(null_r)) / sd(null_r)))

  # ---- 3. kME as stored by R02 ----
  kme <- if (tag == "HF") r2$kmeHF else r2$kmeHCC
  cat("  kme object: ", paste(dim(kme), collapse = " x "),
      " cols: ", paste(head(colnames(kme), 6), collapse = ","), " ...\n")
  mecol <- intersect(c(mod, paste0("ME", mod)), colnames(kme))
  if (length(mecol)) {
    kv <- kme[genes, mecol]
    cat(sprintf("  |kME| for module genes: mean %.3f  median %.3f  min %.3f\n",
                mean(abs(kv)), median(abs(kv)), min(abs(kv))))
  }

  # ---- 4. full-universe connectivity, module vs random ----
  Zs2 <- t(dat[, U, drop = FALSE]); Zs2 <- Zs2 - rowMeans(Zs2)
  s3 <- rowSds(Zs2); s3[s3 == 0] <- 1; Zs2 <- Zs2 / s3
  k <- numeric(length(U))
  for (st in seq(1, length(U), by = 400)) {
    e <- min(st + 399, length(U))
    k[st:e] <- rowSums(adj_from_cor((Zs2[st:e, , drop = FALSE] %*% t(Zs2)) / (n - 1), power))
  }
  cat(sprintf("  full-universe k: all genes mean %.1f median %.1f | module genes mean %.1f median %.1f\n",
              mean(k), median(k), mean(k[gi]), median(k[gi])))
}
cat("\nDIAG DONE\n")
