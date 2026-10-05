# =============================================================================
# R11_composition.R
# -----------------------------------------------------------------------------
# Is the cross-disease difference in translation-program activity just cell
# composition? A tumour is not a hepatocyte and a failing heart is not a
# cardiomyocyte, and ribosomal-protein genes track proliferating and
# biosynthetically active cells. This tests that rather than asserting it.
#
# TWO composition scores per signature, because they answer different questions:
#   mean_z : mean z-scored expression of the marker genes across samples.
#            Simple, but moves with the sample's overall expression level.
#   rank   : mean WITHIN-SAMPLE rank of the marker genes (ESTIMATE-style).
#            This is the better "fraction" proxy: it is invariant to library
#            size and to any monotone per-sample transform, so a sample scores
#            high on "hepatocyte" only if hepatocyte genes are high RELATIVE to
#            everything else in that same sample. The rank score is the primary
#            covariate in the models; mean_z is reported alongside so the two
#            proxies can be compared.
#
# The adjustment is a linear model of the set score on group + covariates. The
# reported effect is beta_group / sigma(model) -- fully standardised, so it is
# directly comparable across models and to Cohen's d, and the attenuation from
# the unadjusted to the adjusted model is readable straight off the table.
#
# Sensitivity: one-at-a-time models for EVERY marker panel (not just the
# primary three), because a reviewer is entitled to see that the choice of
# covariate set does not carry the conclusion.
#
# No new packages, no new data: markers are the hand-curated lists frozen in
# R00_config.R (MARKERS), scored on the cohort matrices R05 already loaded.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R11 start ===")
r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
r6 <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
# r6$tier_of is a bare character vector; R06 names it only positionally via
# ALL_SETS. Key it once here so every lookup below is by set name.
tier_lookup <- setNames(r6$tier_of, r6$ALL_SETS)

# Primary covariate set per disease. Proliferation and immune content are the
# two confounders common to both; the third is the disease-specific purity
# proxy (hepatocyte fraction for tumour, cardiomyocyte fraction for myocardium).
PRIMARY_COVARS <- list(
  HCC = c("hepatocyte", "proliferation", "immune_pan"),
  HF  = c("cardiomyocyte", "fibroblast", "immune_pan")
)

# --- composition scores -------------------------------------------------------
#' Mean within-sample rank of `genes`, per sample. ESTIMATE-style: ranks are
#' computed WITHIN each sample across all measured genes, so the score is a
#' relative-abundance (fraction-like) measure rather than an absolute level.
rank_score <- function(expr, genes) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 3) return(NULL)
  # genes x samples; rank() on a column ranks across genes within that sample.
  # ties.method="average" keeps the score continuous when arrays have ties.
  r <- apply(expr, 2, rank, ties.method = "average")
  list(score = colMeans(r[g, , drop = FALSE]),
       n_used = length(g), n_total = length(genes), used = g)
}

#' Mean z-scored expression of `genes`, per sample.
z_score <- function(expr, genes) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 3) return(NULL)
  z <- t(scale(t(expr[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

score_cohort <- function(co) {
  expr <- co$expr
  out <- list(); cov <- list()
  for (mk in names(MARKERS)) {
    rs <- rank_score(expr, MARKERS[[mk]])
    zs <- z_score(expr, MARKERS[[mk]])
    if (is.null(rs)) {
      log_msg("  ", co$cohort, ": marker panel '", mk, "' has <3 genes present, skipped")
      next
    }
    out[[mk]] <- rs$score
    # The exact absent genes are recorded, not just the count: a reviewer is
    # entitled to know that e.g. the HF discovery array cannot see the lymphoid
    # panel at all, because that bounds what "adjusted for immune content" means.
    cov[[mk]] <- data.frame(
      cohort = co$cohort, panel = mk, n_used = rs$n_used,
      n_total = rs$n_total, frac_used = round(rs$n_used / rs$n_total, 3),
      mean_rank = mean(rs$score), mean_z = mean(zs),
      missing = paste(setdiff(MARKERS[[mk]], rs$used), collapse = ","),
      stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  list(scores = as.data.frame(out), coverage = do.call(rbind, cov))
}

comp <- list(); cov_all <- list()
for (nm in names(r5)) {
  res <- score_cohort(r5[[nm]])
  if (is.null(res)) next
  comp[[nm]] <- res$scores
  cov_all[[length(cov_all) + 1]] <- res$coverage
  log_msg(sprintf("  %s: %d samples x %d panels", nm, nrow(res$scores), ncol(res$scores)))
}
coverage <- do.call(rbind, cov_all)
write_table(coverage, "Table_S23_composition_marker_coverage.csv")

# --- do the composition scores actually shift with disease? -------------------
hedges <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 3 || n2 < 3) return(NA_real_)
  sp <- sqrt(((n1 - 1) * stats::sd(x)^2 + (n2 - 1) * stats::sd(y)^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (1 - 3 / (4 * (n1 + n2) - 9)) * (mean(x) - mean(y)) / sp
}

shift <- list()
for (nm in names(comp)) {
  co <- r5[[nm]]
  cs <- comp[[nm]]
  grp <- co$group[rownames(cs)]
  if (anyNA(grp) || !all(c("disease", "control") %in% grp)) {
    log_msg("  ", nm, ": no two-arm comparison, composition shift not tested")
    next
  }
  d <- rownames(cs)[grp == "disease"]; cc <- rownames(cs)[grp == "control"]
  for (p in colnames(cs)) {
    tt <- tryCatch(stats::t.test(cs[d, p], cs[cc, p]), error = function(e) NULL)
    shift[[length(shift) + 1]] <- data.frame(
      cohort = nm, disease = co$disease, panel = p,
      n_disease = length(d), n_control = length(cc),
      mean_disease = mean(cs[d, p]), mean_control = mean(cs[cc, p]),
      hedges_g = hedges(cs[d, p], cs[cc, p]),
      p_ttest = if (is.null(tt)) NA_real_ else tt$p.value,
      stringsAsFactors = FALSE)
  }
}
shift <- do.call(rbind, shift)
shift$p_fdr <- ave(shift$p_ttest, shift$cohort, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))

cat("\n=== do composition scores shift with disease? (rank score) ===\n")
print(shift[order(shift$cohort, -abs(shift$hedges_g)),
            c("cohort", "panel", "hedges_g", "p_fdr")], row.names = FALSE)
write_table(shift, "Table_S24_composition_shift.csv")

# --- correlation between composition and the translation sets -----------------
#' VIF for the covariates of a design matrix (excluding the intercept and the
#' group term). >10 means the adjustment is collinear and should be distrusted.
vif_of <- function(X) {
  cn <- setdiff(colnames(X), c("(Intercept)", "grpdisease"))
  if (length(cn) < 2) return(setNames(rep(1, length(cn)), cn))
  vapply(cn, function(v) {
    others <- setdiff(cn, v)
    r2 <- summary(stats::lm(X[, v] ~ X[, others, drop = FALSE]))$r.squared
    if (r2 >= 1) Inf else 1 / (1 - r2)
  }, numeric(1))
}

#' Fit score ~ group (+ covariates). Returns the standardised group effect
#' beta/sigma, its p-value, and the model's VIFs.
fit_model <- function(y, grp, covs) {
  grp <- factor(grp, levels = c("control", "disease"))
  df <- data.frame(y = y, grp = grp)
  if (!is.null(covs) && ncol(covs)) {
    covs <- as.data.frame(covs)
    # standardise covariates so the intercept and any comparison across cohorts
    # are on a common scale; z-scoring does not change the group coefficient.
    covs[] <- lapply(covs, function(v) {
      s <- stats::sd(v); if (!is.finite(s) || s == 0) 0 else (v - mean(v)) / s
    })
    df <- cbind(df, covs)
  }
  m <- stats::lm(y ~ ., data = df)
  co <- summary(m)$coefficients
  if (!"grpdisease" %in% rownames(co))
    return(list(beta = NA_real_, p = NA_real_, sigma = NA_real_, vif = NULL, n = nrow(df)))
  list(beta = co["grpdisease", "Estimate"],
       p    = co["grpdisease", "Pr(>|t|)"],
       sigma = summary(m)$sigma,
       vif = vif_of(model.matrix(m)),
       n = nrow(df))
}

adj <- list(); corr_rows <- list()
for (nm in names(comp)) {
  co <- r5[[nm]]
  if (!nm %in% names(r6$scores)) next
  sc <- r6$scores[[nm]]
  cs <- comp[[nm]][colnames(sc), , drop = FALSE]
  grp <- co$group[colnames(sc)]
  if (anyNA(grp) || !all(c("disease", "control") %in% grp))
    { log_msg("  ", nm, ": single-arm, no adjustment possible"); next }
  if (anyNA(cs)) { log_msg("  ", nm, ": NA composition score, skipped"); next }

  dz <- co$disease
  prim <- intersect(PRIMARY_COVARS[[dz]], colnames(cs))
  if (!length(prim)) { log_msg("  ", nm, ": no primary covariate present, skipped"); next }

  # how strongly is each composition score confounded with the set scores?
  ta <- intersect(rownames(sc), SETS_TIER_A)
  for (s in ta) for (p in colnames(cs)) {
    ct <- suppressWarnings(stats::cor.test(sc[s, ], cs[[p]]))
    corr_rows[[length(corr_rows) + 1]] <- data.frame(
      cohort = nm, disease = dz, set = s, panel = p,
      r = unname(ct$estimate), p = ct$p.value, stringsAsFactors = FALSE)
  }

  for (s in rownames(sc)) {
    y <- sc[s, ]
    m0 <- fit_model(y, grp, NULL)
    m1 <- fit_model(y, grp, cs[, prim, drop = FALSE])
    if (is.null(m1$vif) || !all(is.finite(m1$vif)) || any(m1$vif > 10))
      log_msg("  NOTE ", nm, " / ", s, ": VIF > 10, adjustment is collinear")
    # one-at-a-time, every panel
    onz <- lapply(colnames(cs), function(p) fit_model(y, grp, cs[, p, drop = FALSE]))
    names(onz) <- colnames(cs)
    b_onz <- vapply(onz, function(z) z$beta, numeric(1))
    p_onz <- vapply(onz, function(z) z$p, numeric(1))
    worst <- names(which.max(abs(b_onz - m0$beta)))
    adj[[length(adj) + 1]] <- data.frame(
      cohort = nm, disease = dz, set = s, tier = unname(tier_lookup[s]),
      n = m0$n,
      g_unadjusted   = m0$beta / m0$sigma, p_unadjusted = m0$p,
      g_adjusted     = m1$beta / m1$sigma, p_adjusted   = m1$p,
      n_covariates   = length(prim),
      max_vif        = if (length(m1$vif)) max(m1$vif) else NA_real_,
      g_one_at_a_time_worst = b_onz[[worst]] / m1$sigma,
      worst_covariate = worst,
      p_one_at_a_time_worst = p_onz[[worst]],
      covariates = paste(prim, collapse = "+"),
      stringsAsFactors = FALSE)
  }
}
adj <- do.call(rbind, adj)
adj$attenuation <- 1 - adj$g_adjusted / adj$g_unadjusted
adj$sign_kept   <- sign(adj$g_adjusted) == sign(adj$g_unadjusted)
adj$survives    <- adj$sign_kept & abs(adj$g_adjusted) >= STATS$MIN_EFFECT &
                   adj$p_adjusted < STATS$ALPHA
adj$p_fdr_unadj <- ave(adj$p_unadjusted, adj$cohort, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
adj$p_fdr_adj   <- ave(adj$p_adjusted,   adj$cohort, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
write_table(adj, "Table_S25_composition_adjustment.csv")

corr <- do.call(rbind, corr_rows)
corr$p_fdr <- ave(corr$p, corr$cohort, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
write_table(corr, "Table_S26_composition_set_correlation.csv")

# --- summary ------------------------------------------------------------------
cat("\n=== Tier A: does the disease effect survive composition adjustment? ===\n")
ta <- adj[adj$tier == "A", ]
for (nm in unique(ta$cohort)) {
  d <- ta[ta$cohort == nm, ]
  cat(sprintf("\n-- %s (%s; covariates: %s) --\n", nm, d$disease[1], d$covariates[1]))
  d <- d[order(-abs(d$g_unadjusted)), ]
  print(d[, c("set", "g_unadjusted", "g_adjusted", "attenuation",
              "p_unadjusted", "p_adjusted", "survives")], row.names = FALSE)
  cat(sprintf("   survives adjustment: %d / %d sets\n", sum(d$survives), nrow(d)))
}
cat(sprintf("\nTier A across all cohorts: %d / %d set-cohort pairs survive\n",
            sum(ta$survives), nrow(ta)))
cat(sprintf("median attenuation of |g|: %.1f%%\n",
            100 * median(adj$attenuation, na.rm = TRUE)))

cat("\n=== strongest composition confounders of the Tier A sets ===\n")
for (nm in unique(corr$cohort)) {
  d <- corr[corr$cohort == nm, ]
  d <- d[order(-abs(d$r)), ]
  cat(sprintf("\n-- %s --\n", nm))
  print(head(d[, c("set", "panel", "r", "p_fdr")], 5), row.names = FALSE)
}

saveRDS(list(scores = comp, coverage = coverage, shift = shift,
             adjustment = adj, correlation = corr,
             primary_covars = PRIMARY_COVARS),
        file.path(OUT, "intermediate", "R11_composition.rds"))
log_msg("=== R11 done ===")
