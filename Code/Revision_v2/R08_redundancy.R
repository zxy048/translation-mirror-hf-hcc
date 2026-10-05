# =============================================================================
# R08_redundancy.R
# -----------------------------------------------------------------------------
# The canonical translation sets overlap heavily (REACTOME_TRANSLATION alone
# contains most of KEGG_RIBOSOME). R3's objection was that treating 33 such sets
# as 33 independent trials inflates any binomial/count test. Two corrections:
#
#   1. Jaccard-distance hierarchical clustering -> a set of redundancy clusters,
#      with one representative each. Used for the de-redundified binomial test.
#   2. The effective number of independent tests n_eff, from the eigenvalues of
#      the correlation matrix of the per-set ssGSEA scores. TWO estimators are
#      reported rather than one, because they answer the same question by
#      different approximations and a reviewer should be able to see both:
#        Nyholt (2004):   n_eff = 1 + (m - 1) * (1 - Var(lambda)/m)
#        Li & Ji (2005):  n_eff = sum_i [ I(lambda_i >= 1) + (lambda_i - floor(lambda_i)) ]
#      Both are computed on ALL m eigenvalues (Li & Ji use the full spectrum,
#      not just the ones above 1).
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R08 start ===")
r6 <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
r7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
sets_all <- r6$ALL_SETS
tier <- setNames(r6$tier_of, sets_all)
sets_all <- sets_all[!is.na(tier)]

# --- 1. Jaccard distance between gene sets -----------------------------------
resolved <- resolve_sets(sets_all)
J <- jaccard_dist(resolved)
hc <- hclust(as.dist(J), method = "average")

for (h in c(0.9, 0.7, 0.5)) {
  cl <- cutree(hc, h = h)
  cat(sprintf("\n=== Jaccard cut %.1f -> %d clusters (from %d sets) ===\n",
              h, length(unique(cl)), length(cl)))
  for (k in unique(cl)) {
    nm <- names(cl)[cl == k]
    if (length(nm) > 1) cat(sprintf("  [%d] %s\n", k, paste(nm, collapse = "  ~  ")))
  }
}

# Primary cut for de-redundification: 0.7 (sets sharing >30% of genes are one
# test). Recorded in the table so the choice is auditable.
CUT <- 0.7
cl <- cutree(hc, h = CUT)
# representative = the set with the most genes actually measured
gm <- vapply(resolved, length, integer(1))
rep_of <- vapply(unique(cl), function(k) names(cl)[cl == k][which.max(gm[cl == k])],
                 character(1))
assign_df <- data.frame(
  set = names(cl), tier = tier[names(cl)], cluster = as.integer(cl),
  is_representative = names(cl) %in% rep_of,
  n_genes = as.integer(gm[names(cl)]),
  stringsAsFactors = FALSE)
assign_df <- assign_df[order(assign_df$cluster, -assign_df$n_genes), ]
print(assign_df, row.names = FALSE)
write_table(assign_df, "Table_S15_set_redundancy_clusters.csv")

# --- 2. effective number of tests --------------------------------------------
n_eff_from_scores <- function(sc) {
  m <- nrow(sc)
  if (m < 2) return(c(nyholt = m, li_ji = m, m = m))
  R <- stats::cor(t(sc), use = "p")
  R <- (R + t(R)) / 2
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  ev[ev < 0] <- 0
  nyholt <- 1 + (m - 1) * (1 - stats::var(ev) / m)
  li_ji  <- sum(ifelse(ev >= 1, 1, 0) + (ev - floor(ev)))
  c(nyholt = nyholt, li_ji = li_ji, m = m)
}

neff <- do.call(rbind, lapply(names(r6$scores), function(nm) {
  sc <- r6$scores[[nm]]
  A <- intersect(rownames(sc), SETS_TIER_A)
  all_ <- rownames(sc)
  a <- n_eff_from_scores(sc[A, , drop = FALSE])
  b <- n_eff_from_scores(sc[all_, , drop = FALSE])
  data.frame(cohort = nm,
             tierA_m = a[["m"]], tierA_nyholt = round(a[["nyholt"]], 2),
             tierA_li_ji = round(a[["li_ji"]], 2),
             all_m = b[["m"]], all_nyholt = round(b[["nyholt"]], 2),
             all_li_ji = round(b[["li_ji"]], 2),
             stringsAsFactors = FALSE)
}))
cat("\n=== effective number of independent tests ===\n")
print(neff, row.names = FALSE)
write_table(neff, "Table_S16_effective_tests.csv")

cat(sprintf("\nTier A has 15 nominal sets; across cohorts n_eff is about %.1f-%.1f (Nyholt) / %.1f-%.1f (Li & Ji).\n",
            min(neff$tierA_nyholt), max(neff$tierA_nyholt),
            min(neff$tierA_li_ji), max(neff$tierA_li_ji)))

saveRDS(list(clusters = assign_df, hclust = hc, jaccard = J, cut = CUT,
             representatives = rep_of, n_eff = neff),
        file.path(OUT, "intermediate", "R08_redundancy.rds"))
log_msg("=== R08 done ===")
