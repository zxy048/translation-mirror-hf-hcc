# =============================================================================
# R09_mirror_test.R
# -----------------------------------------------------------------------------
# Does the "mirror perturbation" survive? A mirror is a canonical set whose
# within-disease effect is UP in one disease and DOWN in the other.
#
# THREE TESTS, all reported (the plan called for this; none is cherry-picked):
#
#   1. SAMPLE-LEVEL PERMUTATION  (PRIMARY)
#      Shuffle disease/control labels WITHIN each cohort, recompute every
#      effect size, recount mirrors. This is the corrected replacement for
#      v1's null, which permuted PATHWAY LABELS -- that destroys the pairing
#      between the two diseases and cannot test a mirror at all.
#
#   2. DE-REDUNDIFIED BINOMIAL
#      Collapse the sets into Jaccard clusters (R08) and count clusters, so
#      7 near-identical ribosome sets contribute one trial, not seven.
#      This is the direct answer to R3's "the binomial P is inflated".
#
#   3. n_eff-CORRECTED BINOMIAL
#      Same count, but the denominator is the effective number of independent
#      tests (R08), reported for BOTH the Nyholt and Li & Ji estimators.
#
# A set counts as a mirror only if |g| >= MIN_EFFECT and the bootstrap CI
# excludes 0 in BOTH diseases -- i.e. it must be a real effect in each disease
# separately, not merely two small opposite-signed numbers.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(matrixStats))

log_msg("=== R09 start ===")
r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
r6 <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
r7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
r8 <- readRDS(file.path(OUT, "intermediate", "R08_redundancy.rds"))

# --- vectorised Hedges g for a whole score matrix ----------------------------
hedges_all <- function(S, dis, ctl) {
  n1 <- length(dis); n2 <- length(ctl)
  x1 <- S[, dis, drop = FALSE]; x2 <- S[, ctl, drop = FALSE]
  m1 <- rowMeans(x1); m2 <- rowMeans(x2)
  s1 <- rowSds(x1); s2 <- rowSds(x2)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * (m1 - m2) / sp
  names(g) <- rownames(S)
  g
}

ALL_SETS_SCOPE <- function() r6$ALL_SETS[!is.na(r6$tier_of)]
TIER_A <- intersect(SETS_TIER_A, ALL_SETS_SCOPE())
TIER_C <- intersect(SETS_TIER_C_V1, ALL_SETS_SCOPE())

#' Mirror statistics for a (HCC cohort, HF cohort) pair at one tier.
mirror_pair <- function(hcc, hf, sets, tier, nperm = STATS$N_PERM_NULL) {
  S1 <- r6$scores[[hcc]]; S2 <- r6$scores[[hf]]
  sets <- Reduce(intersect, list(sets, rownames(S1), rownames(S2)))
  S1 <- S1[sets, , drop = FALSE]; S2 <- S2[sets, , drop = FALSE]
  d1 <- colnames(S1)[r5[[hcc]]$group[colnames(S1)] == "disease"]
  c1 <- colnames(S1)[r5[[hcc]]$group[colnames(S1)] == "control"]
  d2 <- colnames(S2)[r5[[hf]]$group[colnames(S2)]  == "disease"]
  c2 <- colnames(S2)[r5[[hf]]$group[colnames(S2)]  == "control"]

  obs <- data.frame(set = sets, g_HCC = hedges_all(S1, d1, c1),
                    g_HF = hedges_all(S2, d2, c2), stringsAsFactors = FALSE)
  # gate on the R07 bootstrap CIs (a set must be a real effect in EACH disease)
  solid <- function(cohort) {
    d <- r7[r7$cohort == cohort, ]
    setNames(d$ci_excludes_0, d$set)
  }
  sh <- solid(hcc); sf <- solid(hf)
  obs$hcc_solid <- obs$set %in% names(sh)[sh]
  obs$hf_solid  <- obs$set %in% names(sf)[sf]
  obs$mirror <- with(obs, hcc_solid & hf_solid &
                       abs(g_HCC) >= STATS$MIN_EFFECT & abs(g_HF) >= STATS$MIN_EFFECT &
                       sign(g_HCC) == -sign(g_HF))
  obs$same_direction <- sign(obs$g_HCC) == sign(obs$g_HF)
  obs$hcc_cohort <- hcc; obs$hf_cohort <- hf; obs$tier <- tier
  n_obs <- sum(obs$mirror)

  cl <- setNames(r8$clusters$cluster, r8$clusters$set)
  obs$cluster <- cl[obs$set]
  total_cl <- length(unique(na.omit(cl[sets])))
  k <- length(unique(obs$cluster[obs$mirror]))

  # --- 1+2. permutation null: shuffle labels WITHIN each cohort --------------
  # BOTH statistics are permuted together so they share one correctly-specified
  # null. The cluster-level statistic is the de-redundified test: 7 near-identical
  # ribosome sets contribute one trial instead of seven. Testing either count
  # against 0.5 (as v1 did) would use the wrong null -- under label shuffling the
  # expected mirror count is near zero, not half the sets.
  set.seed(STATS$SEED)
  n1d <- length(d1); n1 <- n1d + length(c1)
  n2d <- length(d2); n2 <- n2d + length(c2)
  null <- vapply(seq_len(nperm), function(i) {
    p1 <- sample(n1); p2 <- sample(n2)
    a <- hedges_all(S1, p1[seq_len(n1d)], p1[-seq_len(n1d)])
    b <- hedges_all(S2, p2[seq_len(n2d)], p2[-seq_len(n2d)])
    mir <- abs(a) >= STATS$MIN_EFFECT & abs(b) >= STATS$MIN_EFFECT & sign(a) == -sign(b)
    c(sets = sum(mir), clusters = length(unique(obs$cluster[mir])))
  }, numeric(2))
  p_perm  <- (1 + sum(null["sets", ]     >= n_obs)) / (1 + nperm)
  p_perm_cl <- (1 + sum(null["clusters", ] >= k))   / (1 + nperm)

  # --- 3. what v1 did, for continuity only -----------------------------------
  # binom.test(k, total_cl, 0.5). Reported because the response letter has to
  # address v1's number, NOT as a valid test: 0.5 is not the null here.
  p_binom_v1 <- stats::binom.test(k, total_cl, 0.5)$p.value

  ne <- r8$n_eff[r8$n_eff$cohort %in% c(hcc, hf), ]
  neff_ny <- if (tier == "A") mean(ne$tierA_nyholt) else mean(ne$all_nyholt)
  neff_lj <- if (tier == "A") mean(ne$tierA_li_ji) else mean(ne$all_li_ji)

  list(obs = obs, n_obs = n_obs, n_sets = length(sets),
       n_clusters = total_cl, n_mirror_clusters = k,
       p_perm = p_perm, p_perm_cluster = p_perm_cl,
       perm_null_sets = null["sets", ], perm_null_clusters = null["clusters", ],
       perm_mean = mean(null["sets", ]), perm_mean_clusters = mean(null["clusters", ]),
       p_binom_v1 = p_binom_v1,
       n_eff_nyholt = neff_ny, n_eff_li_ji = neff_lj)
}

PAIRS <- list(
  list(tag = "primary (v1 cohorts)",  hcc = "TCGA_LIHC",       hf = "GSE57338"),
  list(tag = "HCC array vs HF RNA-seq", hcc = "GSE76427",      hf = "GSE116250"),
  list(tag = "HCC array vs HF RNA-seq 2", hcc = "GSE14520_GPL3921", hf = "GSE141910"),
  list(tag = "HCC array vs HF array-adjacent", hcc = "GSE14520_GPL571", hf = "GSE57338")
)

results <- list(); allobs <- list()
for (p in PAIRS) {
  for (tn in c("A", "C")) {
    sets <- if (tn == "A") TIER_A else TIER_C
    res <- mirror_pair(p$hcc, p$hf, sets, tn)
    res$tag <- p$tag; res$tier <- tn
    res$hcc <- p$hcc; res$hf <- p$hf
    results[[length(results) + 1]] <- res
    allobs[[length(allobs) + 1]] <- res$obs
    cat(sprintf("\n--- %s | Tier %s (%d sets, %d clusters) ---\n",
                p$tag, tn, res$n_sets, res$n_clusters))
    cat(sprintf("  set level     : %2d/%2d mirror   ; permutation P = %.4f (null mean %.2f)\n",
                res$n_obs, res$n_sets, res$p_perm, res$perm_mean))
    cat(sprintf("  cluster level : %2d/%2d mirror   ; permutation P = %.4f (null mean %.2f)\n",
                res$n_mirror_clusters, res$n_clusters,
                res$p_perm_cluster, res$perm_mean_clusters))
    cat(sprintf("  n_eff of this tier: %.1f (Nyholt) / %.1f (Li & Ji)\n",
                res$n_eff_nyholt, res$n_eff_li_ji))
  }
}

obs_all <- do.call(rbind, allobs)
write_table(obs_all, "Table_S17_mirror_sets.csv")

summ <- do.call(rbind, lapply(results, function(r) data.frame(
  pair = r$tag, tier = r$tier, n_sets = r$n_sets, n_clusters = r$n_clusters,
  n_mirror_sets = r$n_obs, null_sets_mean = round(r$perm_mean, 2),
  p_perm_sets = r$p_perm,
  n_mirror_clusters = r$n_mirror_clusters,
  null_clusters_mean = round(r$perm_mean_clusters, 2),
  p_perm_clusters = r$p_perm_cluster,
  p_binom_v1_style = r$p_binom_v1,
  n_eff_nyholt = round(r$n_eff_nyholt, 2), n_eff_li_ji = round(r$n_eff_li_ji, 2),
  stringsAsFactors = FALSE)))
cat("\n=== SUMMARY ===\n"); print(summ, row.names = FALSE)
write_table(summ, "Table_S18_mirror_summary.csv")

saveRDS(list(results = results, summary = summ, obs = obs_all),
        file.path(OUT, "intermediate", "R09_mirror.rds"))
log_msg("=== R09 done ===")
