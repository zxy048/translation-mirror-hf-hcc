# =============================================================================
# R10c_clinical_covariates.R
# -----------------------------------------------------------------------------
# R10 tested one reviewer explanation for the failed HF replication (aetiology).
# R11 tested another (tissue composition). R2 Major 3 asks additionally for
# "stratified or covariate-adjusted analyses by ... relevant clinical
# characteristics", and neither of those covers patient age or sex, which are
# the two clinical variables every one of these public cohorts actually
# publishes. This closes that gap, and it is the last of the reviewer's three
# explanations that can be tested with the deposited data.
#
# WHAT THIS DOES
#   For every cohort whose public annotation carries age AND sex, refit the
#   disease-vs-control effect on each scored set with age and sex in the model,
#   and compare against the same fit without them. The effect is reported as
#   beta_group / sigma -- fully standardised, the same convention R11 uses -- so
#   the adjusted effect is directly comparable with the unadjusted one and with
#   Hedges g. Both fits use the SAME complete-case sample set, so the contrast
#   is paired and any change is attributable to the covariates rather than to
#   the samples.
#
#   The summary that matters is at the bottom: the cross-disease effect-size
#   correlation and the mirror count, recomputed on ADJUSTED effects. If the
#   mirror survived by being an artefact of who was in the cohorts, it would
#   move here.
#
# WHAT IT CANNOT DO, AND WHY THAT IS STATED RATHER THAN SKIPPED
#   NYHA class and left-ventricular ejection fraction are the clinical variables
#   a cardiologist would most want, and they exist in NONE of the three HF
#   cohort series matrices deposited with GEO. That is a property of the public
#   data, not a choice made here; the script enumerates the available annotation
#   columns per cohort so the claim is checkable.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(SummarizedExperiment))

log_msg("=== R10c start ===")
r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
r6 <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
r7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
tier_lookup <- setNames(r6$tier_of, r6$ALL_SETS)

# --- 1. which cohorts can be adjusted at all ----------------------------------
# Enumerate the annotation columns each cohort exposes, so that "not adjust" is
# reported as a fact about the cohort rather than passed over in silence.
log_msg("--- annotation available per cohort ---")
for (nm in names(r5)) {
  a <- r5[[nm]]$ann
  log_msg("  ", nm, ": ",
          if (is.null(a)) "no annotation in the frozen cohort object"
          else paste(colnames(a), collapse = ", "))
}

# --- 2. pull age and sex for the four cohorts that have them -------------------
#' Normalise a sex vector to M/F. The series matrices disagree on case and on
#' which level comes first, and a silent factor-level mismatch would flip the
#' sign of the sex coefficient without any error.
norm_sex <- function(v) {
  s <- tolower(trimws(as.character(v)))
  out <- rep(NA_character_, length(s))
  out[s %in% c("m", "male", "1")] <- "M"
  out[s %in% c("f", "female", "2")] <- "F"
  out[s == ""] <- NA_character_
  out
}

num_age <- function(v) suppressWarnings(as.numeric(as.character(v)))

clin <- list()

# GSE57338 -- HF discovery array. ann is in the frozen object.
a <- r5$GSE57338$ann
clin$GSE57338 <- data.frame(
  sample = a$gsm,
  age = num_age(a$age),
  sex = norm_sex(a$gender),
  source = "series matrix: age, gender",
  stringsAsFactors = FALSE)

# GSE141910 -- HF RNA-seq validation. ann is in the frozen object.
a <- r5$GSE141910$ann
clin$GSE141910 <- data.frame(
  sample = a$gsm,
  age = num_age(a$age),
  sex = norm_sex(a$Sex),
  source = "series matrix: age, Sex",
  stringsAsFactors = FALSE)

# GSE116250 -- HF RNA-seq. R05 passed ann = NULL for this cohort even though the
# parser returns the annotation, so it is re-read here rather than by editing the
# frozen loader: the cohort matrix, its gene filter and its group labels are
# untouched, and only the phenotype columns are picked up.
sm116 <- parse_series_matrix(file.path(VALID_DIR, "GSE116250_series_matrix.txt.gz"))
con <- gzfile(file.path(VALID_DIR, "GSE116250_series_matrix.txt.gz"), "rt")
l <- readLines(con, warn = FALSE); close(con)
tl <- grep("^!Sample_title\t", l, value = TRUE)[1]
titles <- gsub('^"|"$', "", strsplit(tl, "\t", fixed = TRUE)[[1]][-1])
stopifnot(length(titles) == nrow(sm116$ann))
# R05 keys GSE116250's group labels positionally through the same title vector.
# Rebuild that join from the titles and assert it reproduces the labels R05
# stored: if the title vector were misaligned, age and sex would be attached to
# the wrong samples and nothing downstream would raise an error.
st <- setNames(sm116$ann$disease, titles)[colnames(r5$GSE116250$expr)]
derived_grp <- ifelse(grepl("non-?failing|normal", st, ignore.case = TRUE),
                      "control", "disease")
stopifnot(identical(unname(derived_grp),
                    unname(r5$GSE116250$group[colnames(r5$GSE116250$expr)])))
clin$GSE116250 <- data.frame(
  sample = titles,
  age = num_age(sm116$ann$age),
  sex = norm_sex(sm116$ann$Sex),
  source = "series matrix: age, Sex",
  stringsAsFactors = FALSE)
# the matrix is keyed by the same strings
stopifnot(all(colnames(r5$GSE116250$expr) %in% clin$GSE116250$sample))

# TCGA-LIHC -- HCC primary inference cohort. Clinical fields live in colData;
# key on barcode, which is what the cohort matrix's columns are.
se <- readRDS(file.path(DATA_RAW, "TCGA_LIHC_se.rds"))
cd <- as.data.frame(colData(se))
clin$TCGA_LIHC <- data.frame(
  sample = rownames(cd),
  age = num_age(cd$age_at_index),
  sex = norm_sex(cd$gender),
  source = "GDC clinical: age_at_index, gender",
  stringsAsFactors = FALSE)

# --- 3. coverage report -------------------------------------------------------
cov <- do.call(rbind, lapply(names(clin), function(nm) {
  d <- clin[[nm]]
  data.frame(cohort = nm, n_samples = nrow(d),
             n_age = sum(!is.na(d$age)), n_sex = sum(!is.na(d$sex)),
             age_range = sprintf("%g-%g", min(d$age, na.rm = TRUE), max(d$age, na.rm = TRUE)),
             n_female = sum(d$sex == "F", na.rm = TRUE),
             n_male = sum(d$sex == "M", na.rm = TRUE),
             source = d$source[1], stringsAsFactors = FALSE)
}))
log_msg("--- age / sex coverage ---")
print(cov, row.names = FALSE)

# --- 3b. why each remaining cohort is absent, checked rather than asserted ----
# "Not adjusted" has to be a statement about the deposited annotation, not an
# omission, so each excluded cohort is tested for the reason it is excluded.
log_msg("--- cohorts not adjusted, and the reason ---")
for (nm in setdiff(names(r5), names(clin))) {
  a <- r5[[nm]]$ann
  if (is.null(a)) {
    log_msg("  ", nm, ": cohort object carries no annotation at all")
    next
  }
  cols <- tolower(colnames(a))
  has_age <- any(grepl("age", cols)); has_sex <- any(grepl("gender|sex", cols))
  if (!has_age && !has_sex) {
    log_msg("  ", nm, ": no age or sex column in the series matrix (",
            paste(colnames(a), collapse = ", "), ")")
    next
  }
  # age/sex present -- are they present for BOTH arms? A covariate measured on
  # one arm only cannot enter a disease-vs-control model.
  acol <- colnames(a)[grep("age", cols)][1]
  grp <- r5[[nm]]$group[as.character(a$gsm)]
  tab <- table(arm = grp, measured = !is.na(num_age(a[[acol]])) & a[[acol]] != "")
  log_msg("  ", nm, ": ", acol, " present for arm x measured -> ",
          paste(apply(tab, 1, function(r) sprintf("%s:%d/%d", r[1], r[2], sum(r))),
                collapse = "  "))
}

# --- 4. the adjusted fit ------------------------------------------------------
#' Fit score ~ group (+ age + sex) on a complete-case sample set.
#' Returns the standardised group effect beta/sigma, its p-value, and the VIF of
#' each covariate. Both the unadjusted and the adjusted fit are computed on the
#' SAME samples so that the difference between them is the covariates alone.
fit_pair <- function(y, grp, age, sex) {
  keep <- is.finite(y) & !is.na(grp) & is.finite(age) & !is.na(sex)
  y <- y[keep]; grp <- grp[keep]; age <- age[keep]; sex <- sex[keep]
  if (length(unique(grp)) < 2 || sum(grp == "disease") < 3 || sum(grp == "control") < 3)
    return(NULL)
  g <- factor(grp, levels = c("control", "disease"))
  a <- as.numeric(scale(age))
  s <- factor(sex, levels = c("M", "F"))
  if (nlevels(droplevels(s)) < 2) return(NULL)   # sex invariant in this subset

  m0 <- stats::lm(y ~ g)
  m1 <- stats::lm(y ~ g + a + s)
  c0 <- summary(m0)$coefficients; c1 <- summary(m1)$coefficients
  if (!"gdisease" %in% rownames(c1)) return(NULL)

  # VIF of the two clinical covariates, from the design matrix of m1
  X <- stats::model.matrix(m1)[, c("a", "sF"), drop = FALSE]
  vif <- vapply(colnames(X), function(v) {
    r2 <- summary(stats::lm(X[, v] ~ X[, setdiff(colnames(X), v), drop = FALSE]))$r.squared
    if (r2 >= 1) Inf else 1 / (1 - r2)
  }, numeric(1))

  list(n = length(y),
       g_unadj = c0["gdisease", "Estimate"] / summary(m0)$sigma,
       p_unadj = c0["gdisease", "Pr(>|t|)"],
       g_adj   = c1["gdisease", "Estimate"] / summary(m1)$sigma,
       p_adj   = c1["gdisease", "Pr(>|t|)"],
       b_age   = c1["a", "Estimate"] / summary(m1)$sigma,
       p_age   = c1["a", "Pr(>|t|)"],
       b_sex   = c1["sF", "Estimate"] / summary(m1)$sigma,
       p_sex   = c1["sF", "Pr(>|t|)"],
       max_vif = max(vif), n_female = sum(sex == "F"))
}

#' Percentile bootstrap CI for both standardised group effects, resampling
#' within each group exactly as R07's boot_g does, so the adjusted and
#' unadjusted intervals are built the same way and the manuscript's mirror gate
#' (opposite signs, |g| >= 0.2, CI excluding zero in both diseases) can be
#' applied to the adjusted effects without changing the gate.
boot_std_effect <- function(y, grp, age, sex, nboot = STATS$N_BOOT) {
  set.seed(STATS$SEED)
  i_d <- which(grp == "disease"); i_c <- which(grp == "control")
  X0 <- cbind(1, as.numeric(grp == "disease"))
  X1 <- cbind(1, as.numeric(grp == "disease"), age, as.numeric(sex == "F"))
  fit1 <- function(X, idx) {
    Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
    q <- qr(Xb)
    if (q$rank < ncol(Xb)) return(NA_real_)
    b <- qr.coef(q, yb); r <- qr.resid(q, yb)
    s <- sqrt(sum(r^2) / (length(yb) - ncol(Xb)))
    if (!is.finite(s) || s == 0) return(NA_real_)
    unname(b[2] / s)
  }
  bu <- ba <- numeric(nboot)
  for (i in seq_len(nboot)) {
    idx <- c(sample(i_d, replace = TRUE), sample(i_c, replace = TRUE))
    bu[i] <- fit1(X0, idx); ba[i] <- fit1(X1, idx)
  }
  q <- function(v) { v <- v[is.finite(v)]
                     if (length(v) < 100) c(NA_real_, NA_real_)
                     else stats::quantile(v, c(0.025, 0.975), names = FALSE) }
  c(ci_lo_unadj = q(bu)[1], ci_hi_unadj = q(bu)[2],
    ci_lo_adj   = q(ba)[1], ci_hi_adj   = q(ba)[2])
}

adj <- list()
for (nm in names(clin)) {
  co <- r5[[nm]]; sc <- r6$scores[[nm]]
  if (is.null(co) || is.null(sc)) { log_msg("  ", nm, ": not in R06, skipped"); next }
  cl <- clin[[nm]]
  idx <- match(colnames(sc), cl$sample)
  if (anyNA(idx)) { log_msg("  ", nm, ": ", sum(is.na(idx)),
                            " sample(s) unmatched to the annotation, skipped"); next }
  age <- cl$age[idx]; sex <- cl$sex[idx]
  grp <- unname(co$group[colnames(sc)])
  if (anyNA(grp) || !all(c("disease", "control") %in% grp)) {
    log_msg("  ", nm, ": no two-arm contrast, skipped"); next
  }
  n_ok <- sum(is.finite(age) & !is.na(sex) & !is.na(grp))
  log_msg(sprintf("  %s: %d/%d samples with age and sex", nm, n_ok, ncol(sc)))

  for (s in rownames(sc)) {
    y <- as.numeric(sc[s, ])
    r <- fit_pair(y, grp, age, sex)
    if (is.null(r)) next
    # the complete-case set fit_pair used, so the interval belongs to the same
    # estimate and the same samples
    keep <- is.finite(y) & !is.na(grp) & is.finite(age) & !is.na(sex)
    ci <- boot_std_effect(y[keep], grp[keep], age[keep], sex[keep])
    adj[[length(adj) + 1]] <- data.frame(
      cohort = nm, disease = co$disease, set = s, tier = unname(tier_lookup[s]),
      n = r$n, n_female = r$n_female,
      g_unadjusted = r$g_unadj, ci_lo_unadj = ci[["ci_lo_unadj"]],
      ci_hi_unadj  = ci[["ci_hi_unadj"]], p_unadjusted = r$p_unadj,
      g_adjusted   = r$g_adj,   ci_lo_adj   = ci[["ci_lo_adj"]],
      ci_hi_adj    = ci[["ci_hi_adj"]],   p_adjusted   = r$p_adj,
      b_age = r$b_age, p_age = r$p_age, b_sex_female = r$b_sex, p_sex = r$p_sex,
      max_vif = r$max_vif,
      covariates = "age+sex", stringsAsFactors = FALSE)
  }
}
adj <- do.call(rbind, adj)
adj$attenuation <- 1 - adj$g_adjusted / adj$g_unadjusted
adj$sign_kept   <- sign(adj$g_adjusted) == sign(adj$g_unadjusted)
adj$survives    <- adj$sign_kept & abs(adj$g_adjusted) >= STATS$MIN_EFFECT &
                   adj$p_adjusted < STATS$ALPHA
adj$p_fdr_unadj <- ave(adj$p_unadjusted, adj$cohort, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
adj$p_fdr_adj   <- ave(adj$p_adjusted,   adj$cohort, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
adj <- adj[order(adj$cohort, adj$tier, -abs(adj$g_unadjusted)), ]
write_table(adj, "Table_S37_clinical_covariate_adjustment.csv")

# --- 5. does the standardised group effect agree with R07's Hedges g? ---------
# A sanity check on the scale, not a result: the unadjusted beta/sigma and R07's
# Hedges g estimate the same thing on the same samples, so if they disagreed the
# adjusted column could not be compared with the mirror analysis at all.
log_msg("--- scale check: unadjusted beta/sigma vs R07 Hedges g ---")
for (nm in unique(adj$cohort)) {
  d <- adj[adj$cohort == nm, ]
  m <- merge(d[, c("set", "g_unadjusted")], r7[r7$cohort == nm, c("set", "hedges_g")], by = "set")
  if (nrow(m) < 3) next
  log_msg(sprintf("  %s: rho = %.3f over %d sets", nm, cor(m$g_unadjusted, m$hedges_g), nrow(m)))
}

# --- 6. the summary: mirror correlation and count, before and after -----------
# The gate is the manuscript's, applied identically to both the unadjusted and
# the adjusted effects: opposite signs, |g| >= 0.2, and a bootstrap CI excluding
# zero in BOTH diseases. Using one gate on both sides is the whole point -- the
# question is whether the count moves when age and sex enter the model, not
# whether a differently-defined count happens to be larger.
mirror_summary <- function(hcc, hf, tier) {
  s <- intersect(names(hcc$g), names(hf$g))
  if (tier != "ALL") {
    keep <- vapply(s, function(x) identical(unname(tier_lookup[x]), tier), logical(1))
    s <- s[keep]
  }
  g1 <- hcc$g[s]; g2 <- hf$g[s]
  c1 <- hcc$excl[s]; c2 <- hf$excl[s]
  ok <- is.finite(g1) & is.finite(g2)
  s <- s[ok]; g1 <- g1[ok]; g2 <- g2[ok]; c1 <- c1[ok]; c2 <- c2[ok]
  list(n = length(s),
       rho = suppressWarnings(cor(g1, g2, method = "spearman")),
       mirror = sum(sign(g1) != sign(g2) & abs(g1) >= STATS$MIN_EFFECT &
                    abs(g2) >= STATS$MIN_EFFECT & c1 & c2))
}

pull <- function(nm, col_g, col_lo, col_hi) {
  d <- adj[adj$cohort == nm, ]
  list(g = setNames(d[[col_g]], d$set),
       excl = setNames(is.finite(d[[col_lo]]) & is.finite(d[[col_hi]]) &
                       (d[[col_lo]] > 0 | d[[col_hi]] < 0), d$set))
}
tc_un  <- pull("TCGA_LIHC", "g_unadjusted", "ci_lo_unadj", "ci_hi_unadj")
tc_a   <- pull("TCGA_LIHC", "g_adjusted",   "ci_lo_adj",   "ci_hi_adj")

rows <- list()
for (nm in setdiff(unique(adj$cohort), "TCGA_LIHC")) {
  hf_un <- pull(nm, "g_unadjusted", "ci_lo_unadj", "ci_hi_unadj")
  hf_a  <- pull(nm, "g_adjusted",   "ci_lo_adj",   "ci_hi_adj")
  for (ti in c("A", "ALL")) {
    a <- mirror_summary(tc_a,  hf_a,  ti)
    u <- mirror_summary(tc_un, hf_un, ti)
    rows[[length(rows) + 1]] <- data.frame(
      hf_cohort = nm, tier = ti, n_sets = a$n,
      rho_unadjusted = u$rho, mirror_unadjusted = u$mirror,
      rho_adjusted   = a$rho, mirror_adjusted   = a$mirror,
      gate = "opposite signs, |g| >= 0.2, CI excludes 0 in both diseases",
      stringsAsFactors = FALSE)
  }
}
mir <- do.call(rbind, rows)
write_table(mir, "Table_S38_age_sex_mirror_recheck.csv")

cat("\n=== TCGA-LIHC effect, adjusted vs unadjusted (Tier A) ===\n")
print(adj[adj$cohort == "TCGA_LIHC" & adj$tier == "A",
          c("set", "n", "g_unadjusted", "g_adjusted", "attenuation",
            "p_unadjusted", "p_adjusted", "survives")], row.names = FALSE)

cat("\n=== does the mirror survive age and sex adjustment? ===\n")
print(mir, row.names = FALSE)

cat("\n=== age and sex effects on the Tier A sets (is either a confounder?) ===\n")
ta <- adj[adj$tier == "A", ]
cat(sprintf("  |b_age|  median %.3f, max %.3f, sets with p_age < 0.05: %d/%d\n",
            median(abs(ta$b_age)), max(abs(ta$b_age)),
            sum(ta$p_age < 0.05), nrow(ta)))
cat(sprintf("  |b_sex|  median %.3f, max %.3f, sets with p_sex < 0.05: %d/%d\n",
            median(abs(ta$b_sex_female)), max(abs(ta$b_sex_female)),
            sum(ta$p_sex < 0.05), nrow(ta)))
cat(sprintf("  max VIF across all fits: %.2f\n", max(adj$max_vif, na.rm = TRUE)))
cat(sprintf("  Tier A set-cohort pairs: %d; sign kept %d; survives gate %d\n",
            nrow(ta), sum(ta$sign_kept), sum(ta$survives)))
cat(sprintf("  median attenuation of |g|: %.1f%%\n", 100 * median(ta$attenuation, na.rm = TRUE)))

saveRDS(list(coverage = cov, adjustment = adj, mirror_recheck = mir, clinical = clin),
        file.path(OUT, "intermediate", "R10c_clinical.rds"))
log_msg("=== R10c done ===")
