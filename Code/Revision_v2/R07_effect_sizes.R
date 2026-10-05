# =============================================================================
# R07_effect_sizes.R
# -----------------------------------------------------------------------------
# Within-disease disease-vs-control effect size per canonical set, per cohort.
#
# WHAT IS CORRECTED vs v1 (03_ssGSEA_cross_disease_pathway.R, lines ~215-218):
#   v1:  d = (m1 - m2) / sqrt((s1^2 + s2^2) / 2)      EQUAL-WEIGHT pooled SD
#   v2:  g = Hedges-corrected d from an n-WEIGHTED pooled SD
#                 s_p = sqrt(((n1-1)s1^2 + (n2-1)s2^2) / (n1+n2-2))
#                 J   = 1 - 3/(4*(n1+n2) - 9);  g = J * (m1-m2)/s_p
# The equal-weight form is biased whenever n1 != n2, which is the norm here
# (TCGA 371 vs 50; GSE116250 50 vs 14). The estimand is unchanged -- this is the
# same quantity, correctly estimated.
#
# SIGN CONVENTION: g > 0 means HIGHER IN DISEASE (tumour, or failing heart).
# Bootstrap CI resamples SAMPLES within each arm (not pathway labels -- that was
# a second v1 defect).
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R07 start ===")
r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
r6 <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
scores <- r6$scores
tier <- setNames(r6$tier_of, r6$ALL_SETS)

hedges_g <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 3 || n2 < 3) return(c(d = NA, g = NA, sp = NA, J = NA))
  s1 <- stats::sd(x); s2 <- stats::sd(y)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(c(d = NA, g = NA, sp = sp, J = NA))
  d  <- (mean(x) - mean(y)) / sp
  J  <- 1 - 3 / (4 * (n1 + n2) - 9)
  c(d = d, g = J * d, sp = sp, J = J)
}

boot_g <- function(x, y, nboot = STATS$N_BOOT) {
  set.seed(STATS$SEED)
  v <- replicate(nboot, {
    xb <- sample(x, replace = TRUE); yb <- sample(y, replace = TRUE)
    hedges_g(xb, yb)[["g"]]
  })
  v <- v[is.finite(v)]
  if (length(v) < 100) return(c(NA, NA))
  stats::quantile(v, c(0.025, 0.975), names = FALSE)
}

rows <- list()
for (nm in names(scores)) {
  co <- r5[[nm]]
  if (!all(c("disease", "control") %in% co$group)) {
    log_msg("--- ", nm, ": no control arm, skipped (no effect size) ---")
    next
  }
  sc <- scores[[nm]]
  g <- co$group[colnames(sc)]          # group is named by sample ID (R05 mk())
  if (anyNA(g)) stop(nm, ": ", sum(is.na(g)), " scored samples have no group label")
  dis <- colnames(sc)[g == "disease"]
  ctl <- colnames(sc)[g == "control"]
  log_msg("--- ", nm, ": ", length(dis), " disease vs ", length(ctl), " control ---")

  for (s in rownames(sc)) {
    x <- sc[s, dis]; y <- sc[s, ctl]
    h <- hedges_g(x, y)
    ci <- boot_g(x, y)
    tt <- tryCatch(stats::t.test(x, y), error = function(e) NULL)
    rows[[length(rows) + 1]] <- data.frame(
      cohort = nm, disease = co$disease, tech = co$tech, set = s,
      tier = tier[[s]],
      n_disease = length(dis), n_control = length(ctl),
      mean_disease = mean(x), mean_control = mean(y),
      sd_disease = stats::sd(x), sd_control = stats::sd(y),
      cohens_d_equalweight = (mean(x) - mean(y)) / sqrt((stats::sd(x)^2 + stats::sd(y)^2) / 2),
      pooled_sd_nweighted = h[["sp"]], hedges_J = h[["J"]],
      d_nweighted = h[["d"]], hedges_g = h[["g"]],
      ci_lo = ci[1], ci_hi = ci[2],
      p_ttest = if (is.null(tt)) NA_real_ else tt$p.value,
      stringsAsFactors = FALSE)
  }
}

es <- do.call(rbind, rows)
es$ci_excludes_0 <- with(es, is.finite(ci_lo) & is.finite(ci_hi) &
                           (ci_lo > 0 | ci_hi < 0))
es$meets_min_effect <- abs(es$hedges_g) >= STATS$MIN_EFFECT
es$mirror_candidate <- es$ci_excludes_0 & es$meets_min_effect
es$p_fdr <- ave(es$p_ttest, es$cohort, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))

write_table(es, "Table_S14_effect_sizes_within_disease.csv")

cat("\n=== how much does the correction move the estimate? ===\n")
delta <- es$hedges_g - es$cohens_d_equalweight
cat(sprintf("equal-weight d vs n-weighted Hedges g: mean |delta| = %.3f, max |delta| = %.3f\n",
            mean(abs(delta), na.rm = TRUE), max(abs(delta), na.rm = TRUE)))
cat(sprintf("cases where the two disagree on SIGN: %d / %d\n",
            sum(sign(es$cohens_d_equalweight) != sign(es$hedges_g), na.rm = TRUE),
            sum(is.finite(es$hedges_g))))

cat("\n=== Tier A effect sizes, per cohort (sorted by |g|) ===\n")
ta <- es[es$tier == "A", ]
for (nm in unique(ta$cohort)) {
  d <- ta[ta$cohort == nm, ]; d <- d[order(-abs(d$hedges_g)), ]
  cat(sprintf("\n-- %s (%s, %d vs %d) --\n", nm, d$disease[1], d$n_disease[1], d$n_control[1]))
  print(head(d[, c("set", "hedges_g", "ci_lo", "ci_hi", "p_fdr")], 6), row.names = FALSE)
}

saveRDS(es, file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
log_msg("=== R07 done: ", nrow(es), " set x cohort effect sizes ===")
