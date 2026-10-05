# =============================================================================
# R10_heterogeneity.R
# -----------------------------------------------------------------------------
# The cohorts do not agree (R07/R09 showed GSE141910 running opposite to the
# other HF cohorts). This quantifies that instead of hiding it.
#
#   - Q and I^2 per canonical set across the cohorts of each disease, on the
#     Hedges g scale (DerSimonian-Laird random effects).
#   - Within HF, the discovery cohort is stratified by aetiology (DCM vs ICM):
#     if the two aetiologies disagree, "HF" is not one thing and the HF side of
#     the mirror is aetiology-dependent.
#
# A set is called homogeneous only if I^2 < 50%; everything else is reported as
# heterogeneous with its I^2, because selective reporting of the concordant
# cohorts is exactly the failure mode this revision exists to fix.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R10 start ===")
r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
r6 <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
r7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
suppressPackageStartupMessages(library(matrixStats))

# --- DerSimonian-Laird Q / I^2 on the g scale ---------------------------------
dl <- function(g, se) {
  ok <- is.finite(g) & is.finite(se) & se > 0
  g <- g[ok]; se <- se[ok]; k <- length(g)
  if (k < 2) return(c(k = k, Q = NA, df = k - 1, p_Q = NA, I2 = NA,
                      g_pooled = if (k) g[1] else NA, se_pooled = if (k) se[1] else NA))
  w <- 1 / se^2
  gbar <- sum(w * g) / sum(w)
  Q <- sum(w * (g - gbar)^2)
  df <- k - 1
  p_Q <- stats::pchisq(Q, df, lower.tail = FALSE)
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / C)
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  wr <- 1 / (se^2 + tau2)
  c(k = k, Q = Q, df = df, p_Q = p_Q, I2 = I2,
    g_pooled = sum(wr * g) / sum(wr), se_pooled = sqrt(1 / sum(wr)))
}

# SE of Hedges g from the R07 bootstrap CI (width / (2*1.96))
r7$se_g <- (r7$ci_hi - r7$ci_lo) / (2 * 1.959964)

het <- do.call(rbind, lapply(split(r7, r7$disease), function(d) {
  do.call(rbind, lapply(split(d, d$set), function(s) {
    r <- dl(s$hedges_g, s$se_g)
    data.frame(disease = d$disease[1], set = s$set[1], tier = s$tier[1],
               n_cohorts = r[["k"]], Q = r[["Q"]], p_Q = r[["p_Q"]],
               I2_pct = round(r[["I2"]], 1),
               g_pooled = r[["g_pooled"]], se_pooled = r[["se_pooled"]],
               cohorts = paste(s$cohort, collapse = "; "),
               stringsAsFactors = FALSE)
  }))
}))
het <- het[order(het$disease, -abs(het$g_pooled)), ]
write_table(het, "Table_S21_cohort_heterogeneity.csv")

cat("=== heterogeneity across cohorts, per disease (Tier A) ===\n")
for (dz in c("HCC", "HF")) {
  d <- het[het$disease == dz & het$tier == "A", ]
  cat(sprintf("\n-- %s (%d cohorts) --\n", dz, d$n_cohorts[1]))
  print(head(d[, c("set", "n_cohorts", "g_pooled", "I2_pct", "p_Q")], 15),
        row.names = FALSE)
  cat(sprintf("   I2 < 50%% (homogeneous): %d / %d sets\n",
              sum(d$I2_pct < 50, na.rm = TRUE), nrow(d)))
}

# --- HF stratified by aetiology (DCM vs ICM) in the discovery cohort ---------
log_msg("--- GSE57338 stratified by aetiology ---")
ann <- r5$GSE57338$ann
sc  <- r6$scores$GSE57338
aet <- setNames(ann$etiology, ann$gsm)[colnames(sc)]
grp <- r5$GSE57338$group[colnames(sc)]
hedges_all <- function(S, a, b) {
  n1 <- length(a); n2 <- length(b)
  x1 <- S[, a, drop = FALSE]; x2 <- S[, b, drop = FALSE]
  sp <- sqrt(((n1 - 1) * rowSds(x1)^2 + (n2 - 1) * rowSds(x2)^2) / (n1 + n2 - 2))
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  setNames(J * (rowMeans(x1) - rowMeans(x2)) / sp, rownames(S))
}
cf <- colnames(sc)
nf <- cf[grp == "control"]
dcm <- cf[grp == "disease" & aet == "DCM"]
icm <- cf[grp == "disease" & aet == "ICM"]
log_msg("  NF=", length(nf), " DCM=", length(dcm), " ICM=", length(icm))

strat <- data.frame(
  set = rownames(sc),
  g_DCM_vs_NF = hedges_all(sc, dcm, nf),
  g_ICM_vs_NF = hedges_all(sc, icm, nf),
  stringsAsFactors = FALSE)
strat$tier <- setNames(r6$tier_of, r6$ALL_SETS)[strat$set]
strat$same_direction <- sign(strat$g_DCM_vs_NF) == sign(strat$g_ICM_vs_NF)
write_table(strat, "Table_S22_HF_aetiology_stratified.csv")

cat("\n=== GSE57338: aetiology-stratified effect sizes (Tier A) ===\n")
sa <- strat[strat$tier == "A", ]
print(sa[order(-abs(sa$g_DCM_vs_NF)), c("set", "g_DCM_vs_NF", "g_ICM_vs_NF",
                                        "same_direction")], row.names = FALSE)
cat(sprintf("\nTier A sets where DCM and ICM point the SAME way: %d / %d\n",
            sum(sa$same_direction), nrow(sa)))

saveRDS(list(heterogeneity = het, stratified = strat,
             groups = list(nf = nf, dcm = dcm, icm = icm)),
        file.path(OUT, "intermediate", "R10_heterogeneity.rds"))
log_msg("=== R10 done ===")
