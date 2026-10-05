# =============================================================================
# R10b_hf_aetiology.R
# -----------------------------------------------------------------------------
# WHY THIS IS SEPARATE FROM R10
# R10 stratified the HF DISCOVERY cohort (GSE57338: DCM 82 / ICM 95 / NF 136).
# Probing the cohort annotations turned up a second, larger problem that R10 did
# not cover:
#
#   GSE141910 is NOT an "HCM cohort". It is 166 DCM + 28 HCM + 6 PPCM vs 166
#   non-failing donors. Manuscript v1 used only its 28 HCM samples, versus the
#   same 166 donors, and described that subset as the study's "same-organ
#   cardiac control" (rho = -0.036). v2's R05 instead pooled all 200 diseased
#   samples into one "HF validation cohort".
#
# Those are two different uses of one dataset, and neither was checked against
# the data. This script checks: what does each aetiology actually do, and does
# the DCM-only subset (166 vs 166 -- a precisely balanced design, and the
# aetiology that dominates both HF cohorts) still contradict the discovery
# cohort?
#
# This is a stratification of analyses the reviewers already asked for
# (R2 #3, R3 #7), not a new analysis: every quantity here is a Hedges g on the
# same Tier A sets, computed the same way as R07.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(matrixStats))

log_msg("=== R10b start ===")
r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
r6 <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
r7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
tier_lookup <- setNames(r6$tier_of, r6$ALL_SETS)

hedges_all <- function(S, a, b) {
  n1 <- length(a); n2 <- length(b)
  x1 <- S[, a, drop = FALSE]; x2 <- S[, b, drop = FALSE]
  sp <- sqrt(((n1 - 1) * rowSds(x1)^2 + (n2 - 1) * rowSds(x2)^2) / (n1 + n2 - 2))
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  setNames(J * (rowMeans(x1) - rowMeans(x2)) / sp, rownames(S))
}

# --- 1. per-aetiology effect sizes in each HF cohort with annotations --------
strat_one <- function(nm) {
  co <- r5[[nm]]; sc <- r6$scores[[nm]]
  et <- co$ann$etiology
  if (is.null(et)) { log_msg("  ", nm, ": no etiology annotation, skipped"); return(NULL) }
  aet <- setNames(et, co$ann$gsm)[colnames(sc)]
  grp <- co$group[colnames(sc)]
  cf  <- colnames(sc)
  ctl <- cf[grp == "control"]
  if (!length(ctl)) { log_msg("  ", nm, ": no control arm"); return(NULL) }
  cat(sprintf("\n-- %s: %d control (donor/NF) --\n", nm, length(ctl)))
  print(sort(table(aet[grp == "disease"]), decreasing = TRUE))

  A <- intersect(rownames(sc), SETS_TIER_A)
  rows <- list()
  for (a in sort(unique(na.omit(aet[grp == "disease"])))) {
    k <- cf[grp == "disease" & !is.na(aet) & aet == a]
    if (length(k) < 5) {
      log_msg("  ", nm, " / ", a, ": only ", length(k), " samples, not stratified")
      next
    }
    g <- hedges_all(sc[A, , drop = FALSE], k, ctl)
    rows[[length(rows) + 1]] <- data.frame(
      cohort = nm, aetiology = a, n_aetiology = length(k), n_control = length(ctl),
      set = A, tier = unname(tier_lookup[A]), hedges_g = g, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

strat <- do.call(rbind, lapply(c("GSE57338", "GSE141910"), strat_one))
strat$direction <- ifelse(strat$hedges_g > 0, "up in disease", "down in disease")
write_table(strat, "Table_S22b_HF_aetiology_effect_sizes.csv")

# --- 2. does each subgroup still oppose HCC (TCGA-LIHC)? ----------------------
hcc_g <- r7[r7$cohort == "TCGA_LIHC" & r7$tier == "A", c("set", "hedges_g", "ci_excludes_0")]
names(hcc_g)[2] <- "g_HCC"

cmp <- do.call(rbind, lapply(split(strat, interaction(strat$cohort, strat$aetiology)),
  function(d) {
    m <- merge(d, hcc_g, by = "set")
    if (nrow(m) < 5) return(NULL)
    rho <- suppressWarnings(stats::cor(m$hedges_g, m$g_HCC, method = "spearman"))
    mir <- sum(sign(m$hedges_g) == -sign(m$g_HCC) &
                 abs(m$hedges_g) >= STATS$MIN_EFFECT & abs(m$g_HCC) >= STATS$MIN_EFFECT)
    # the same gate R09 uses: a real effect in the HCC arm and the right sign here
    gated <- m$ci_excludes_0 & abs(m$g_HCC) >= STATS$MIN_EFFECT
    data.frame(cohort = d$cohort[1], aetiology = d$aetiology[1],
               n_aetiology = d$n_aetiology[1], n_control = d$n_control[1],
               n_sets = nrow(m), spearman_rho_vs_HCC = rho,
               mirror_sets_abs_g_ge_020 = mir,
               mirror_sets_gated = sum(mir_g <- (sign(m$hedges_g) == -sign(m$g_HCC) &
                 abs(m$hedges_g) >= STATS$MIN_EFFECT)[gated]),
               n_gated = sum(gated),
               mean_g = mean(m$hedges_g), stringsAsFactors = FALSE)
  }))
cmp <- cmp[order(cmp$cohort, -abs(cmp$spearman_rho_vs_HCC)), ]
write_table(cmp, "Table_S22c_HF_aetiology_vs_HCC.csv")

cat("\n=== each HF aetiology subgroup vs TCGA-LIHC (Tier A) ===\n")
cat("(negative rho / many mirrors = the subgroup reproduces the discovery pattern)\n")
print(cmp, row.names = FALSE)

# --- 3. DCM in the two cohorts: do they agree with each other? ----------------
# The two GEO series annotate the same aetiology differently -- GSE57338 says
# "DCM", GSE141910 says "Dilated cardiomyopathy (DCM)". Stripping the
# parenthetical does NOT reconcile them (the stems are "DCM" vs "Dilated
# cardiomyopathy"); the acronym inside the parentheses is the common key. A
# join on the raw label silently returns zero rows -- it did, twice.
acronym <- function(x) ifelse(grepl("\\(", x),
                              sub(".*\\(([^)]*)\\).*", "\\1", x), x)
strat$acronym <- acronym(strat$aetiology)
d1 <- strat[strat$cohort == "GSE57338" & strat$acronym == "DCM", c("set", "hedges_g")]
d2 <- strat[strat$cohort == "GSE141910" & strat$acronym == "DCM", c("set", "hedges_g")]
stopifnot(nrow(d1) > 0, nrow(d2) > 0)
names(d1)[2] <- "g_GSE57338_DCM"; names(d2)[2] <- "g_GSE141910_DCM"
dd <- merge(d1, d2, by = "set")
rho_dcm <- suppressWarnings(stats::cor(dd$g_GSE57338_DCM, dd$g_GSE141910_DCM,
                                       method = "spearman"))
cat(sprintf("\n=== DCM vs DCM: %s vs %s ===\n", "GSE57338", "GSE141910"))
cat(sprintf("  %d Tier A sets; Spearman rho = %+.3f\n", nrow(dd), rho_dcm))
cat(sprintf("  same direction in both cohorts: %d / %d\n",
            sum(sign(dd$g_GSE57338_DCM) == sign(dd$g_GSE141910_DCM)), nrow(dd)))
print(dd[order(-abs(dd$g_GSE57338_DCM)), ], row.names = FALSE)
write_table(dd, "Table_S22d_DCM_across_cohorts.csv")

saveRDS(list(stratified = strat, vs_hcc = cmp, dcm_cross = dd,
             rho_dcm_cross = rho_dcm),
        file.path(OUT, "intermediate", "R10b_hf_aetiology.rds"))
log_msg("=== R10b done ===")
