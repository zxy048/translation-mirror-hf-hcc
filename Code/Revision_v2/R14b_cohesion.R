#!/usr/bin/env Rscript
# =============================================================================
# R14b_cohesion.R -- within-disease, connectivity-free module cohesion.
#
# WHY THIS IS A SEPARATE SCRIPT RATHER THAN PART OF R14
# -----------------------------------------------------
# R14's matched and unmatched nulls each draw 10,000 random gene sets and
# recompute a full signed TOM for every draw, which is a ~13 h run. The
# statistic computed here needs no TOM at all: it is the mean pairwise Pearson
# correlation among a module's genes measured in that module's OWN expression
# matrix, against a uniform random gene-set null. That is seconds of work.
#
# It nevertheless has to be a deposited artefact, because the manuscript quotes
# it. An earlier draft cited "mean pairwise r = 0.3916 against a null of 0.0239
# (Z = 26.0)" from an ad-hoc diagnostic log (R14_diag_log.txt) that no committed
# script reproduced and no table carried. That is precisely the defect class this
# revision exists to remove -- a number in the prose with no traceable producer.
#
# WHY THE STATISTIC MATTERS
# -------------------------
# It is the reason the matched null is the primary one. TOM is normalised by
# connectivity, so a uniform random gene set is not an exchangeable null for it,
# and both designated modules sit far below the universe's connectivity:
#
#   HF  purple   153 genes -- module k mean  6.6 vs universe 15.2
#   HCC magenta  146 genes -- module k mean 21.1 vs universe 74.1
#
# Under the uniform TOM null the HCC module scores Z = -0.17 -- apparently not
# cohesive at all -- while this connectivity-free measure of the same property
# says it is strongly cohesive. Same genes, same expression matrix, opposite
# conclusions; the difference is entirely the normalisation. Reporting the
# connectivity-free result alongside the two TOM nulls is what shows the
# unmatched null to be the problem rather than the module.
#
# Both are written to Table_S36; S32 carries the two TOM nulls and S33 the
# cross-disease coherence, so this is the within-disease counterpart to S33.
#
# Seed and permutation count follow STATS$SEED and R14's coherence null so the
# two scripts' uniform nulls are directly comparable.
#
# Run: Rscript --vanilla R14b_cohesion.R
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(matrixStats))

log_msg("=== R14b cohesion start ===")

r2 <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))
r3 <- readRDS(file.path(OUT, "intermediate", "R03_translation_module.rds"))
U  <- r2$U

# Same default as R14's cross-disease coherence null (1,000). Overridable so the
# script can be smoke-tested without touching a production filename.
N_PERM <- as.integer(Sys.getenv("R14B_N_PERM", "1000"))
SMOKE  <- N_PERM != 1000
tn <- function(x) if (SMOKE) paste0("SMOKE_", x) else x
if (SMOKE) log_msg("*** SMOKE MODE: ", N_PERM, " permutations ***")

# Column-centred and scaled genes x samples, so every correlation is one gemm.
# Identical to R14's standardise(); duplicated by necessity because R14 does not
# expose it, and kept byte-identical so the two scripts cannot drift.
standardise <- function(dat) {
  Z <- t(dat)
  Z <- Z - rowMeans(Z)
  s <- rowSds(Z); s[s == 0] <- 1
  Z / s
}

# Mean pairwise correlation among the rows indexed by ii, upper triangle only.
mean_pairwise_r <- function(Zs, ii) {
  r <- (Zs[ii, , drop = FALSE] %*% t(Zs[ii, , drop = FALSE])) / (ncol(Zs) - 1)
  mean(r[upper.tri(r)])
}

run_disease <- function(tag) {
  res   <- r2[[tag]]
  dat   <- r2[[if (tag == "HF") "datHF" else "datHCC"]]
  colors <- res$net$colors
  mod   <- r3[[if (tag == "HF") "desHF" else "desHCC"]]$module
  if (is.null(names(colors))) names(colors) <- U
  genes <- names(colors)[colors == mod]
  m     <- length(genes)

  Zs  <- standardise(dat[, U, drop = FALSE])
  idx <- match(genes, U)
  stopifnot(!any(is.na(idx)))

  obs <- mean_pairwise_r(Zs, idx)

  set.seed(STATS$SEED)
  null <- vapply(seq_len(N_PERM),
                 function(i) mean_pairwise_r(Zs, sample(nrow(Zs), m)),
                 numeric(1))

  z <- (obs - mean(null)) / sd(null)
  # One-sided empirical P with the observed value counted, matching R14's
  # convention: (1 + #{null >= obs}) / (1 + n).
  p <- (1 + sum(null >= obs)) / (1 + N_PERM)

  log_msg(sprintf("  %s %s (%d genes): mean pairwise r = %.4f | null %.4f (sd %.4f) | Z = %.2f | P = %.4f",
                  tag, mod, m, obs, mean(null), sd(null), z, p))

  data.frame(disease = tag, module = mod, n_genes = m,
             n_samples = nrow(dat), n_genes_universe = nrow(Zs),
             obs_mean_pairwise_r = signif(obs, 4),
             null_mean_r = signif(mean(null), 4),
             null_sd_r = signif(sd(null), 4),
             Z = round(z, 2), p = signif(p, 3), n_perm = N_PERM,
             stringsAsFactors = FALSE)
}

tab <- do.call(rbind, lapply(c("HF", "HCC"), run_disease))

# The connectivity figures that motivate the matched null, carried in the same
# table so a reader can see the confound and its correction together. Read from
# S32 if it exists -- S32 is the authority for k, and re-deriving it here would
# create a second source that could disagree.
s32p <- file.path(OUT, "tables", tn("Table_S32_module_internal_robustness.csv"))
if (file.exists(s32p)) {
  s32 <- utils::read.csv(s32p)
  key <- paste(s32$disease, s32$module)
  j   <- match(paste(tab$disease, tab$module), key)
  tab$k_module_mean   <- s32$k_module_mean[j]
  tab$k_universe_mean <- s32$k_universe_mean[j]
  tab$k_module_over_universe <- round(s32$k_module_mean[j] / s32$k_universe_mean[j], 3)
} else {
  log_msg("  NOTE: S32 not present yet; k columns left out of S36")
}

write_table(tab, tn("Table_S36_within_disease_cohesion.csv"))
cat("\n=== R14b within-disease connectivity-free cohesion ===\n")
print(tab, row.names = FALSE)
log_msg("=== R14b done ===")
