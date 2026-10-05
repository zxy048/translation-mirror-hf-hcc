# =============================================================================
# R14_tom_robustness.R
# -----------------------------------------------------------------------------
# Reproduce v1's module-internal robustness analysis on the harmonised modules.
#
# v1 claimed (Manuscript_PLOSONE_v1.md:51,105,249) that the green module's mean
# intra-modular TOM, compared against 10,000 random gene sets of equal size,
# gave Z = 62.8, P < 0.0001, plus a jackknife CV of 1.6% and a cross-disease
# coherence Z = 2.9. v2's pipeline (R00-R13) computes none of it, so without
# this script the claim would have to be dropped or carried over unverified.
#
# TWO NULLS ARE COMPUTED, AND THE UNMATCHED ONE IS NOT THE PRIMARY.
#
#   unmatched : 10,000 uniformly random gene sets of size m -- this is exactly
#               v1's null, kept so the response letter can show what changes and
#               why. It is reported for comparison only.
#   matched   : 10,000 gene sets drawn so that each member is matched to its
#               module counterpart on full-universe connectivity k. THIS is the
#               primary null.
#
# Why matching is required. TOM is connectivity-normalised: its denominator is
# min(k_i,k_j) + 1 - a_ij. A uniformly random gene set is therefore not an
# exchangeable null unless the module's genes happen to have the same k
# distribution as the universe, and they do not. Measured on the harmonised
# networks (R14_diag_log.txt):
#
#   HF  purple   153 genes -- module k mean  6.6 vs universe mean 15.2
#   HCC magenta  146 genes -- module k mean 21.1 vs universe mean 74.1
#
# Both modules are unambiguously cohesive by a connectivity-free measure (mean
# pairwise r, Z = 41.0 and 26.0 respectively), yet under the unmatched null HCC
# magenta scores Z = -0.17. That is the confound, not an absence of structure:
# low-k genes have both a small numerator and a small denominator, so their TOM
# ratio is dragged back toward the random ratio. The matched null removes it.
#
# TOM FORMULA. As verified against WGCNA::TOMsimilarity on real data (off-
# diagonal exact to ~4e-16 for both diseases), with A the signed adjacency
# (0.5*(1+cor))^power, a_ii = 1, and k = rowSums(A) INCLUDING the diagonal:
#
#     TOM_ij = ( (A A^T)_ij - a_ij ) / ( min(k_i, k_j) - a_ij )
#
# and WGCNA then sets the diagonal of its returned TOM to 1.
#
# Both halves of this were guessed wrong on the first attempt, twice over, and
# both errors were plausible-looking rather than obviously broken:
#
#   numerator   guessed (A A^T)_ij + a_ij; wrong by ~0.83. A A^T sums over ALL
#               u, so sum_{u != i,j} a_iu a_uj = (A A^T)_ij - 2 a_ij, giving
#               (A A^T)_ij - a_ij.
#   denominator guessed min(k_i,k_j) + 1 - a_ij, i.e. the documented form with
#               the self-term left in k; wrong by exactly 1. The code that
#               matches uses min(k_i,k_j) - a_ij, equivalently
#               min(k_i-1, k_j-1) + 1 - a_ij, i.e. the documented formula with a
#               connectivity that EXCLUDES self.
#
# A wrong numerator and a wrong null move the observed statistic and the null in
# the same direction, so the Z-score can look perfectly reasonable while being
# wrong. Hence SELFTEST, which pins both the adjacency and the TOM against
# WGCNA's own constructors and STOPS the script on any mismatch. Do not
# "simplify" either expression without re-running it.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(matrixStats))

log_msg("=== R14 start ===")
r2  <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))
r3  <- readRDS(file.path(OUT, "intermediate", "R03_translation_module.rds"))
U   <- r2$U
G   <- length(U)

# Overridable so the script can be smoke-tested end to end before the real run.
# Production values are the defaults; SMOKE forces a filename prefix so a smoke
# run can never overwrite a production artefact.
env_int <- function(nm, default) {
  v <- Sys.getenv(nm, "")
  if (nzchar(v)) as.integer(v) else default
}
N_PERM_TOM <- env_int("R14_N_PERM_TOM", 10000)   # v1 used 10,000
N_PERM_COH <- env_int("R14_N_PERM_COH", 1000)    # v1 used 1,000
JK_ITER    <- env_int("R14_JK_ITER", 100)
JK_KEEP    <- 0.90
BLOCK      <- 400L    # rows per connectivity block
N_BINS     <- 20L     # connectivity strata for the matched null
SMOKE      <- N_PERM_TOM != 10000 || N_PERM_COH != 1000 || JK_ITER != 100
if (SMOKE) log_msg("*** SMOKE MODE: ", N_PERM_TOM, " TOM perms, ", N_PERM_COH,
                   " coherence perms, ", JK_ITER, " jackknife iters ***")

# ---------------------------------------------------------------------------
# core pieces
# ---------------------------------------------------------------------------
# x^p by repeated squaring. R's `^` calls R_pow per element and is the single
# slowest step in the null loop (1.5e6 elements x 10,000 iterations x 2 nulls).
# p = 14 becomes 3 squarings and 2 multiplies, with no clamp: the base is
# (1 + r)/2, which is in [0,1] for |r| <= 1, and correlation can only exceed 1
# by floating-point noise that the power annihilates anyway.
fast_pow <- function(x, p) {
  if (p == 1L) return(x)
  res <- NULL; base <- x; e <- p
  while (e > 0) {
    if (e %% 2L == 1L) res <- if (is.null(res)) base else res * base
    e <- e %/% 2L
    if (e > 0) base <- base * base
  }
  res
}
adj_from_cor <- function(r, power) fast_pow(0.5 * (1 + r), power)

# genes x samples, column-centred and scaled, so every correlation is one gemm
standardise <- function(dat) {
  Z <- t(dat)
  Z <- Z - rowMeans(Z)
  s <- rowSds(Z); s[s == 0] <- 1
  Z / s
}

# k_i = sum_u a_iu over ALL genes, blockwise over rows
full_connectivity <- function(Zs, power, n) {
  k <- numeric(nrow(Zs))
  for (s in seq(1, nrow(Zs), by = BLOCK)) {
    e <- min(s + BLOCK - 1L, nrow(Zs))
    r <- (Zs[s:e, , drop = FALSE] %*% t(Zs)) / (n - 1)
    k[s:e] <- rowSums(adj_from_cor(r, power))
    rm(r); gc(verbose = FALSE)
  }
  k
}

# mean intra-set TOM, upper triangle; u-sums and k over ALL genes in Zs.
# Inner loop of a 10,000-permutation test: the m x G correlation matrix is
# computed once and the set's own adjacency is sliced out of it. No gc() here --
# at m ~ 150 the temporary is ~12 MB and calling gc() 10,000 times costs more
# than the garbage it collects.
mean_intra_tom <- function(Zs, power, n, idx, k) {
  A   <- adj_from_cor((Zs[idx, , drop = FALSE] %*% t(Zs)) / (n - 1), power)
  a   <- A[, idx, drop = FALSE]                  # the set's own m x m adjacency
  ks  <- k[idx]
  tom <- (A %*% t(A) - a) / (outer(ks, ks, pmin) - a)
  mean(tom[upper.tri(tom)])
}

# Connectivity-matched sampling. Homogeneous k-bins over the universe; a module
# gene in bin b is replaced by a uniformly drawn member of bin b. Bins are wide
# enough (m/20 ~ 7 genes per bin) that matching is approximate by construction,
# which under-corrects rather than over-corrects -- the conservative direction.
make_bins <- function(k, nb = N_BINS) {
  br  <- unique(quantile(k, probs = seq(0, 1, length.out = nb + 1)))
  b   <- findInterval(k, br, all.inside = TRUE)
  members <- split(seq_along(k), b)
  list(bin = b, members = members)
}
sample_matched <- function(module_bins, bins, exclude) {
  out <- integer(length(module_bins))
  for (i in seq_along(module_bins)) {
    pool <- bins$members[[as.character(module_bins[i])]]
    pool <- pool[!pool %in% exclude]
    out[i] <- if (length(pool) == 0L) NA_integer_ else
      pool[if (length(pool) == 1L) 1L else sample.int(length(pool), 1L)]
  }
  out
}

# ---------------------------------------------------------------------------
# SELFTEST -- must reproduce WGCNA before any number is reported
# ---------------------------------------------------------------------------
selftest <- function(dat, power, tag) {
  n <- nrow(dat)
  set.seed(STATS$SEED)
  sub <- dat[, sample(ncol(dat), 1200), drop = FALSE]
  Zs  <- standardise(sub)
  r   <- (Zs %*% t(Zs)) / (n - 1)
  A   <- pmin(pmax(adj_from_cor(r, power), 0), 1)   # clamp only for WGCNA's checks
  k   <- rowSums(A)

  ref   <- WGCNA::TOMsimilarity(A, TOMType = "signed")
  mine  <- (A %*% t(A) - A) / (outer(k, k, pmin) - A)
  off   <- upper.tri(ref)
  d_tom <- max(abs(mine[off] - ref[off]))

  # adjacency itself, against WGCNA's own constructor
  aw    <- WGCNA::adjacency(sub, power = power, type = "signed")
  d_adj <- max(abs(A - aw))

  log_msg(sprintf("  SELFTEST %s: TOM off-diag max|diff| = %.3e ; adjacency max|diff| = %.3e",
                  tag, d_tom, d_adj))
  if (!is.finite(d_tom) || !is.finite(d_adj) || d_tom > 1e-10 || d_adj > 1e-10)
    stop(sprintf(paste0("SELFTEST FAILED for %s (TOM %.3e, adjacency %.3e). The ",
                        "hand-rolled network does not reproduce WGCNA, so neither ",
                        "the observed statistic nor the null below is trustworthy."),
                 tag, d_tom, d_adj))
  invisible(c(tom = d_tom, adj = d_adj))
}

# ---------------------------------------------------------------------------
# per disease
# ---------------------------------------------------------------------------
run_disease <- function(tag) {
  res    <- r2[[tag]]
  dat    <- r2[[if (tag == "HF") "datHF" else "datHCC"]]
  power  <- res$power
  colors <- res$net$colors
  mod    <- r3[[if (tag == "HF") "desHF" else "desHCC"]]$module
  if (is.null(names(colors))) names(colors) <- U   # blockwiseModules normally
  genes  <- names(colors)[colors == mod]           # names these; belt and braces
  n      <- nrow(dat)
  m      <- length(genes)

  log_msg("--- ", tag, ": module ", mod, " (", m, " genes), power ", power,
          ", n ", n, " samples ---")

  selftest(dat[, U, drop = FALSE], power, tag)

  Zs  <- standardise(dat[, U, drop = FALSE])
  k   <- full_connectivity(Zs, power, n)
  idx <- match(genes, U)
  stopifnot(!any(is.na(idx)))
  log_msg(sprintf("  full-universe k: all genes mean %.1f median %.1f | module mean %.1f median %.1f",
                  mean(k), median(k), mean(k[idx]), median(k[idx])))
  obs <- mean_intra_tom(Zs, power, n, idx, k)

  # --- unmatched null: v1's method, reported for comparison only ---
  set.seed(STATS$SEED)
  null_un <- vapply(seq_len(N_PERM_TOM), function(i) {
    if (i %% 2500 == 0) log_msg("    unmatched null ", i, "/", N_PERM_TOM)
    mean_intra_tom(Zs, power, n, sample(G, m), k)
  }, numeric(1))

  # --- connectivity-matched null: primary ---
  bins     <- make_bins(k)
  mod_bins <- bins$bin[idx]
  set.seed(STATS$SEED)
  null_mat <- vapply(seq_len(N_PERM_TOM), function(i) {
    if (i %% 2500 == 0) log_msg("    matched null ", i, "/", N_PERM_TOM)
    ii <- sample_matched(mod_bins, bins, idx)
    if (anyNA(ii)) return(NA_real_)
    mean_intra_tom(Zs, power, n, ii, k)
  }, numeric(1))
  n_na <- sum(is.na(null_mat))
  if (n_na > 0) log_msg("  WARNING: ", n_na, " matched sets had an empty bin and were dropped")
  null_mat <- null_mat[!is.na(null_mat)]

  z_un  <- (obs - mean(null_un))  / sd(null_un)
  p_un  <- (1 + sum(null_un  >= obs)) / (1 + length(null_un))
  z_mat <- (obs - mean(null_mat)) / sd(null_mat)
  p_mat <- (1 + sum(null_mat >= obs)) / (1 + length(null_mat))
  log_msg(sprintf("  UNMATCHED null: obs %.4f | null %.4f (sd %.4f) | Z = %.2f | P = %.4f",
                  obs, mean(null_un), sd(null_un), z_un, p_un))
  log_msg(sprintf("  MATCHED   null: obs %.4f | null %.4f (sd %.4f) | Z = %.2f | P = %.4f",
                  obs, mean(null_mat), sd(null_mat), z_mat, p_mat))

  # --- jackknife: drop 10% of module genes, recompute intra-modular connectivity
  kim <- function(ii) {
    a <- adj_from_cor((Zs[ii, , drop = FALSE] %*% t(Zs[ii, , drop = FALSE])) / (n - 1), power)
    rowSums(a)
  }
  set.seed(STATS$SEED)
  keep_n <- floor(JK_KEEP * m)
  jk <- replicate(JK_ITER, kim(sample(idx, keep_n)))
  mean_kim  <- colMeans(jk)
  cv_across <- sd(mean_kim) / abs(mean(mean_kim))
  # Per-gene CV is heavy-tailed: genes whose kIM is near zero give enormous
  # ratios and a plain mean over genes is dominated by them. Report the median,
  # and count the degenerate genes so the tail is visible rather than hidden.
  per_gene_cv <- apply(jk, 2, function(x) sd(x) / abs(mean(x)))
  cv_gene_med <- median(per_gene_cv)
  n_degen     <- sum(colMeans(jk) < 0.05)

  log_msg(sprintf("  jackknife: CV of mean kIM = %.2f%% (v1's 1.6%% is this quantity) ; median per-gene CV = %.2f%% ; %d/%d genes with kIM < 0.05",
                  100 * cv_across, 100 * cv_gene_med, n_degen, m))

  # --- cross-disease coherence: this module's genes in the OTHER matrix.
  # Mean pairwise r is not connectivity-normalised, so the uniform null is
  # exchangeable here and no matching is needed.
  other_tag <- if (tag == "HF") "HCC" else "HF"
  other_gx  <- r2[[if (other_tag == "HF") "gxHF" else "gxHCC"]]     # genes x samples
  Zo <- standardise(t(other_gx))
  present <- intersect(genes, rownames(other_gx))
  oi <- match(present, rownames(other_gx))
  n_other <- ncol(Zo)
  coh <- function(ii) {
    r <- (Zo[ii, , drop = FALSE] %*% t(Zo[ii, , drop = FALSE])) / (n_other - 1)
    mean(r[upper.tri(r)])
  }
  obs_coh <- coh(oi)
  set.seed(STATS$SEED)
  null_coh <- vapply(seq_len(N_PERM_COH),
                     function(i) coh(sample(nrow(Zo), length(oi))), numeric(1))
  z_coh <- (obs_coh - mean(null_coh)) / sd(null_coh)
  p_coh <- (1 + sum(null_coh >= obs_coh)) / (1 + N_PERM_COH)
  log_msg(sprintf("  coherence in %s: %d of %d module genes present; obs mean r = %.4f ; Z = %.2f ; P = %.4f",
                  other_tag, length(present), m, obs_coh, z_coh, p_coh))

  list(tag = tag, module = mod, n_genes = m, power = power, n_samples = n,
       k_universe_mean = mean(k), k_module_mean = mean(k[idx]),
       obs_tom = obs,
       null_unmatched_mean = mean(null_un), null_unmatched_sd = sd(null_un),
       z_tom_unmatched = z_un, p_tom_unmatched = p_un,
       null_matched_mean = mean(null_mat), null_matched_sd = sd(null_mat),
       z_tom_matched = z_mat, p_tom_matched = p_mat,
       n_perm_tom = N_PERM_TOM, n_bins = N_BINS,
       jk_cv_across = cv_across, jk_cv_gene_median = cv_gene_med,
       jk_n_degenerate = n_degen, jk_iter = JK_ITER, jk_keep = JK_KEEP,
       coh_other = other_tag, coh_n_present = length(present),
       coh_obs = obs_coh, coh_null_mean = mean(null_coh), coh_null_sd = sd(null_coh),
       coh_z = z_coh, coh_p = p_coh, n_perm_coh = N_PERM_COH)
}

resHF  <- run_disease("HF")
resHCC <- run_disease("HCC")

# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------
# Smoke runs must never be able to land on a production filename.
tn <- function(x) if (SMOKE) paste0("SMOKE_", x) else x

tab <- do.call(rbind, lapply(list(resHF, resHCC), function(x)
  data.frame(disease = x$tag, module = x$module, n_genes = x$n_genes,
             power = x$power, n_samples = x$n_samples,
             k_universe_mean = round(x$k_universe_mean, 1),
             k_module_mean = round(x$k_module_mean, 1),
             obs_mean_intra_TOM = signif(x$obs_tom, 4),
             unmatched_null_mean = signif(x$null_unmatched_mean, 4),
             unmatched_null_sd = signif(x$null_unmatched_sd, 4),
             # Z is NOT comparable between the two nulls: matching collapses the
             # null's variance (sd 1.8e-5 vs 3.5e-4), which inflates Z
             # mechanically. The obs/null mean ratio is the interpretable
             # effect size and the one to quote.
             obs_over_unmatched_null = round(x$obs_tom / x$null_unmatched_mean, 2),
             Z_unmatched = round(x$z_tom_unmatched, 2),
             P_unmatched = signif(x$p_tom_unmatched, 3),
             matched_null_mean = signif(x$null_matched_mean, 4),
             matched_null_sd = signif(x$null_matched_sd, 4),
             obs_over_matched_null = round(x$obs_tom / x$null_matched_mean, 2),
             Z_matched = round(x$z_tom_matched, 2),
             P_matched = signif(x$p_tom_matched, 3),
             n_perm = x$n_perm_tom, n_bins = x$n_bins,
             jackknife_CV_of_mean_kIM_pct = round(100 * x$jk_cv_across, 2),
             jackknife_median_pergene_CV_pct = round(100 * x$jk_cv_gene_median, 2),
             jackknife_n_degenerate = x$jk_n_degenerate,
             stringsAsFactors = FALSE)))
write_table(tab, tn("Table_S32_module_internal_robustness.csv"))

coh <- do.call(rbind, lapply(list(resHF, resHCC), function(x)
  data.frame(module_disease = x$tag, module = x$module, tested_in = x$coh_other,
             n_genes_module = x$n_genes, n_genes_present = x$coh_n_present,
             obs_mean_pairwise_r = signif(x$coh_obs, 4),
             null_mean_r = signif(x$coh_null_mean, 4),
             null_sd_r = signif(x$coh_null_sd, 4),
             Z = round(x$coh_z, 2), p = signif(x$coh_p, 3),
             n_perm = x$n_perm_coh, stringsAsFactors = FALSE)))
write_table(coh, tn("Table_S33_cross_disease_module_coherence.csv"))

cat("\n=== R14 module internal robustness ===\n"); print(tab, row.names = FALSE)
cat("\n=== R14 cross-disease coherence ===\n");    print(coh, row.names = FALSE)

saveRDS(list(HF = resHF, HCC = resHCC, table = tab, coherence = coh,
             params = list(N_PERM_TOM = N_PERM_TOM, N_PERM_COH = N_PERM_COH,
                           JK_ITER = JK_ITER, JK_KEEP = JK_KEEP, N_BINS = N_BINS,
                           universe_size = G, block = BLOCK, smoke = SMOKE)),
        file.path(OUT, "intermediate", tn("R14_tom_robustness.rds")))
log_msg("=== R14 done ===")
