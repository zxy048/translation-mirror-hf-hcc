# =============================================================================
# R30_paired_sensitivity.R
# -----------------------------------------------------------------------------
# STAGE 1 ONLY: verify the pairing metadata and write the paired-sample
# manifest. No effect sizes, bootstrap, or permutation are computed here.
#
# The paired sensitivity analysis is restricted to cohorts in which a patient
# can be shown, from deposited sample metadata, to contribute one disease and
# one control sample. Complete pairs only; unmatched samples are never forced
# into a pair.
#
# Pairing keys:
#   GSE14520 (both platforms): `Individual` ends in A (disease) or B (control);
#       the patient key is `Individual` with the trailing A/B removed.
#   GSE76427: `patient.id`.
#   TCGA-LIHC: barcode characters 1-12 (patient); characters 14-15 are the
#       sample type (01 = primary tumour, 11 = solid tissue normal).
#   GSE57338 / GSE116250 / GSE141910: unpaired case-control designs.
#   GSE141198: tumour-only, no control arm.
#
# Outputs:
#   Table_S40a_paired_sample_manifest.csv  (one row per analysed sample)
#   Table_S40b_pairing_summary.csv         (one row per cohort)
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R30 start: pairing manifest only (no statistics) ===")
r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))

PAIRED_COHORTS <- c(
  "GSE14520_GPL3921",
  "GSE14520_GPL571",
  "GSE76427",
  "TCGA_LIHC"
)

PAIRING_KEY <- c(
  GSE14520_GPL3921 = "Individual stripped of trailing A/B (A=disease, B=control)",
  GSE14520_GPL571 = "Individual stripped of trailing A/B (A=disease, B=control)",
  GSE76427 = "patient.id",
  TCGA_LIHC = "barcode chars 1-12; type chars 14-15 (01=disease, 11=control)",
  GSE57338 = "none (unpaired case-control)",
  GSE116250 = "none (unpaired case-control)",
  GSE141910 = "none (unpaired case-control)",
  GSE141198 = "none (tumour-only, no control arm)"
)

UNPAIRED_REASON <- c(
  GSE57338 = "unpaired case-control design; not part of paired sensitivity",
  GSE116250 = "unpaired case-control design; not part of paired sensitivity",
  GSE141910 = "unpaired case-control design; not part of paired sensitivity",
  GSE141198 = "no control arm; not part of paired sensitivity"
)

make_manifest <- function(nm) {
  co <- r5[[nm]]
  ids <- names(co$group)
  grp <- unname(co$group)
  stopifnot(length(ids) == length(grp), !anyNA(ids), !anyDuplicated(ids))
  stopifnot(all(grp %in% c("disease", "control")))

  paired <- nm %in% PAIRED_COHORTS
  key <- rep(NA_character_, length(ids))

  if (nm %in% c("GSE14520_GPL3921", "GSE14520_GPL571")) {
    a <- co$ann
    stopifnot(!is.null(a), all(c("gsm", "Individual") %in% colnames(a)))
    a <- a[match(ids, a$gsm), , drop = FALSE]
    stopifnot(!anyNA(a$gsm), !anyNA(a$Individual))
    if (any(!grepl("[AB]$", a$Individual))) {
      stop(nm, ": Individual value without a trailing A/B suffix")
    }
    if (any(grp == "disease" & grepl("B$", a$Individual)) ||
        any(grp == "control" & grepl("A$", a$Individual))) {
      stop(nm, ": A/B suffix does not agree with the disease/control group")
    }
    key <- sub("[AB]$", "", a$Individual)
  } else if (nm == "GSE76427") {
    a <- co$ann
    stopifnot(!is.null(a), all(c("gsm", "patient.id") %in% colnames(a)))
    a <- a[match(ids, a$gsm), , drop = FALSE]
    stopifnot(!anyNA(a$gsm), !anyNA(a$patient.id), all(nzchar(a$patient.id)))
    key <- a$patient.id
  } else if (nm == "TCGA_LIHC") {
    key <- substr(ids, 1, 12)
    type <- substr(ids, 14, 15)
    if (any(!type %in% c("01", "11"))) {
      stop("TCGA_LIHC: unexpected sample type code(s): ",
           paste(unique(type[!type %in% c("01", "11")]), collapse = ","))
    }
    if (any(type == "01" & grp != "disease") || any(type == "11" & grp != "control")) {
      stop("TCGA_LIHC: barcode sample type does not agree with the group")
    }
  }

  pair_complete <- rep(NA, length(ids))
  included <- rep(FALSE, length(ids))
  pair_id <- rep(NA_character_, length(ids))
  exclusion_reason <- rep(NA_character_, length(ids))

  if (paired) {
    tab <- table(key, grp)
    both <- rownames(tab)[tab[, "disease"] > 0 & tab[, "control"] > 0]
    bad <- both[tab[both, "disease"] != 1 | tab[both, "control"] != 1]
    if (length(bad)) {
      stop(nm, ": patient(s) with more than one sample in an arm: ",
           paste(bad, collapse = ","))
    }
    pair_complete <- key %in% both
    included <- pair_complete
    pair_id[included] <- paste(nm, key[included], sep = ":")
    exclusion_reason[!included] <- "no paired sample in the other arm"
  } else {
    exclusion_reason[] <- unname(UNPAIRED_REASON[[nm]])
  }

  data.frame(
    cohort = nm,
    disease = co$disease,
    sample_id = ids,
    group = grp,
    paired_design = paired,
    patient_key = key,
    pair_complete = pair_complete,
    included_in_paired = included,
    pair_id = pair_id,
    exclusion_reason = exclusion_reason,
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, lapply(names(r5), make_manifest))

# ------------------------------------------------------------------ summary --
summarise_cohort <- function(nm) {
  m <- manifest[manifest$cohort == nm, , drop = FALSE]
  paired <- unique(m$paired_design)
  stopifnot(length(paired) == 1)

  if (paired) {
    n_pairs <- length(unique(m$patient_key[m$pair_complete]))
    n_unique <- length(unique(m$patient_key))
    n_unmatched_d <- sum(m$group == "disease" & !m$pair_complete)
    n_unmatched_c <- sum(m$group == "control" & !m$pair_complete)
  } else {
    n_pairs <- 0L
    n_unique <- NA_integer_
    n_unmatched_d <- NA_integer_
    n_unmatched_c <- NA_integer_
  }

  data.frame(
    cohort = nm,
    disease = unique(m$disease),
    paired_design = paired,
    n_samples = nrow(m),
    n_disease = sum(m$group == "disease"),
    n_control = sum(m$group == "control"),
    n_unique_patient_keys = n_unique,
    n_complete_pairs = n_pairs,
    n_unmatched_disease = n_unmatched_d,
    n_unmatched_control = n_unmatched_c,
    n_samples_in_complete_pairs = sum(m$included_in_paired),
    n_samples_not_in_complete_pairs = sum(!m$included_in_paired),
    pairing_key = unname(PAIRING_KEY[[nm]]),
    stringsAsFactors = FALSE
  )
}

pair_summary <- do.call(rbind, lapply(names(r5), summarise_cohort))

# ------------------------------- cross-platform GSE14520 overlap check --------
k3921 <- unique(manifest$patient_key[
  manifest$cohort == "GSE14520_GPL3921" & !is.na(manifest$patient_key)])
k571 <- unique(manifest$patient_key[
  manifest$cohort == "GSE14520_GPL571" & !is.na(manifest$patient_key)])
overlap <- intersect(k3921, k571)
log_msg("GSE14520 cross-platform patient-key overlap: ", length(overlap))
if (length(overlap)) {
  log_msg("WARNING: shared keys: ", paste(head(overlap, 20), collapse = ","))
}

# ------------------------------------------------------------------- write ---
write_table(manifest, "Table_S40a_paired_sample_manifest.csv")
write_table(pair_summary, "Table_S40b_pairing_summary.csv")

log_msg("Pairing summary:")
for (i in seq_len(nrow(pair_summary))) {
  s <- pair_summary[i, ]
  log_msg(sprintf(
    "  %-20s paired=%-5s samples=%d (%d disease / %d control) pairs=%s unmatched=%s/%s",
    s$cohort, s$paired_design, s$n_samples, s$n_disease, s$n_control,
    ifelse(s$paired_design, as.character(s$n_complete_pairs), "NA"),
    ifelse(s$paired_design, as.character(s$n_unmatched_disease), "NA"),
    ifelse(s$paired_design, as.character(s$n_unmatched_control), "NA")
  ))
}

log_msg("=== R30 stage 1 complete: manifest written ===")

# =============================================================================
# STAGE 2: paired-design sensitivity statistics
# -----------------------------------------------------------------------------
# Effect size (FIXED BEFORE COMPUTATION):
#   d_i      = score_disease,i - score_control,i   (one value per patient pair)
#   d_z      = mean(d_i) / sd(d_i)                 (paired standardised mean difference)
#   J_p      = 1 - 3 / (4 * n_pairs - 5)           (small-sample correction)
#   g_paired = J_p * d_z
# The denominator is the SD of WITHIN-PAIR DIFFERENCES. This is NOT the
# independent-sample Hedges g denominator (n-weighted pooled between-group SD),
# and the two quantities are reported in separate columns and never combined.
# The |g| >= 0.20 mirror gate is NOT applied to g_paired.
#
# Resampling unit = complete patient pair.
#   Bootstrap: pairs resampled with replacement; percentile 95% CI.
#   Permutation: within each pair the disease/control labels are swapped
#   (sign flip of d_i). ONE sign matrix is generated per cohort and reused for
#   every Tier A set, so the cross-set correlation structure is preserved.
#   Two-sided P = (1 + #{|stat_perm| >= |stat_obs|}) / (N_PERM + 1).
#
# The paired subset is a sensitivity analysis on complete pairs only. It does
# not replace the primary independent-sample analysis, and no cross-disease
# correlation is computed here.
# =============================================================================
RUN_STAGE2 <- TRUE

if (RUN_STAGE2) {
  r6 <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
  r7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
  scores <- r6$scores

  N_BOOT <- STATS$N_BOOT
  N_PERM <- STATS$N_PERM_NULL

  # paired mean/SD/d_z for every column of a pairs x sets matrix
  paired_col_stats <- function(M) {
    n <- nrow(M)
    mu <- colMeans(M)
    ss <- colMeans(M^2) - mu^2
    sdv <- sqrt(pmax(ss, 0) * n / (n - 1))
    dz <- mu / sdv
    dz[!is.finite(dz)] <- NA_real_
    list(mean = mu, sd = sdv, dz = dz)
  }

  # independent-sample Hedges g, identical formula to R07 (used only to
  # describe the complete-pair subset analysed without pairing)
  hedges_g_unpaired <- function(x, y) {
    n1 <- length(x); n2 <- length(y)
    if (n1 < 3 || n2 < 3) return(NA_real_)
    sp <- sqrt(((n1 - 1) * stats::sd(x)^2 + (n2 - 1) * stats::sd(y)^2) /
                 (n1 + n2 - 2))
    if (!is.finite(sp) || sp == 0) return(NA_real_)
    d <- (mean(x) - mean(y)) / sp
    J <- 1 - 3 / (4 * (n1 + n2) - 9)
    J * d
  }

  run_paired_cohort <- function(nm, cohort_index) {
    co <- r5[[nm]]
    sc <- scores[[nm]]
    mm <- manifest[manifest$cohort == nm & manifest$included_in_paired, ,
                   drop = FALSE]
    pids <- unique(mm$pair_id)
    if (!length(pids)) stop(nm, ": no complete pairs")

    d_ids <- ctl_ids <- rep(NA_character_, length(pids))
    for (i in seq_along(pids)) {
      p <- mm[mm$pair_id == pids[i], , drop = FALSE]
      if (nrow(p) != 2L ||
          !setequal(p$group, c("disease", "control"))) {
        stop(nm, ": pair ", pids[i], " is not exactly one disease + one control")
      }
      d_ids[i]   <- p$sample_id[p$group == "disease"]
      ctl_ids[i] <- p$sample_id[p$group == "control"]
    }
    stopifnot(all(d_ids %in% colnames(sc)), all(ctl_ids %in% colnames(sc)))

    sets <- intersect(SETS_TIER_A, rownames(sc))
    log_msg("--- ", nm, ": ", length(pids), " complete pairs, ",
            length(sets), " Tier A sets scored ---")

    # sets x pairs -> pairs x sets
    D <- t(sc[sets, d_ids, drop = FALSE] - sc[sets, ctl_ids, drop = FALSE])

    # pairs usable for the shared sign matrix must be complete across all sets
    ok_pair <- rowSums(!is.finite(D)) == 0L
    if (any(!ok_pair)) {
      log_msg("  excluding ", sum(!ok_pair),
              " pair(s) with non-finite Tier A scores")
    }
    D <- D[ok_pair, , drop = FALSE]
    n_pairs <- nrow(D)
    if (n_pairs < 5) stop(nm, ": fewer than 5 analysable pairs")

    obs <- paired_col_stats(D)
    Jp <- 1 - 3 / (4 * n_pairs - 5)
    g_obs <- Jp * obs$dz

    # ------------------------------------------------ bootstrap over pairs ---
    set.seed(STATS$SEED + 100L * cohort_index)
    boot <- matrix(NA_real_, nrow = N_BOOT, ncol = length(sets))
    for (b in seq_len(N_BOOT)) {
      idx <- sample.int(n_pairs, n_pairs, replace = TRUE)
      boot[b, ] <- Jp * paired_col_stats(D[idx, , drop = FALSE])$dz
    }
    ci <- apply(boot, 2, function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 100) return(c(NA_real_, NA_real_))
      stats::quantile(v, c(0.025, 0.975), names = FALSE)
    })

    # --------------------------------- within-pair sign-flip permutation -----
    # one sign matrix, shared by every set (preserves cross-set correlation)
    set.seed(STATS$SEED + 100000L + 100L * cohort_index)
    signs <- matrix(sample(c(-1, 1), n_pairs * N_PERM, replace = TRUE),
                    nrow = N_PERM, ncol = n_pairs)
    perm <- matrix(NA_real_, nrow = N_PERM, ncol = length(sets))
    for (k in seq_len(N_PERM)) {
      Dk <- D * signs[k, ]
      perm[k, ] <- paired_col_stats(Dk)$dz
    }
    p_perm <- vapply(seq_along(sets), function(j) {
      ps <- perm[, j]
      ps <- ps[is.finite(ps)]
      if (!is.finite(obs$dz[j]) || !length(ps)) return(NA_real_)
      (1 + sum(abs(ps) >= abs(obs$dz[j]))) / (length(ps) + 1)
    }, numeric(1))

    p_fdr <- stats::p.adjust(p_perm, method = STATS$FDR_METHOD)

    # paired t-test as an audit column only (permutation is primary)
    p_t <- vapply(seq_along(sets), function(j) {
      if (!is.finite(obs$dz[j])) return(NA_real_)
      tryCatch(stats::t.test(D[, j], mu = 0)$p.value,
               error = function(e) NA_real_)
    }, numeric(1))

    # ------------------------- independent-sample comparators (no mixing) ----
    full <- r7[r7$cohort == nm & r7$set %in% sets, , drop = FALSE]
    full <- full[match(sets, full$set), , drop = FALSE]

    g_subset <- vapply(sets, function(s) {
      hedges_g_unpaired(sc[s, d_ids], sc[s, ctl_ids])
    }, numeric(1))

    res <- data.frame(
      cohort = nm,
      disease = co$disease,
      set = sets,
      tier = "A",
      n_pairs_manifest = length(pids),
      n_pairs_analysed = n_pairs,
      n_disease_samples_in_pairs = n_pairs,
      n_control_samples_in_pairs = n_pairs,
      mean_within_pair_diff = obs$mean,
      sd_within_pair_diff = obs$sd,
      d_z = obs$dz,
      hedges_J_paired = Jp,
      g_paired_z = g_obs,
      paired_ci_lo = ci[1, ],
      paired_ci_hi = ci[2, ],
      paired_ci_excludes_0 = is.finite(ci[1, ]) & is.finite(ci[2, ]) &
        (ci[1, ] > 0 | ci[2, ] < 0),
      p_perm_paired = p_perm,
      p_fdr_within_cohort = p_fdr,
      p_paired_t_audit = p_t,
      n_disease_full_cohort = full$n_disease,
      n_control_full_cohort = full$n_control,
      hedges_g_full_cohort = full$hedges_g,
      p_fdr_full_cohort = full$p_fdr,
      g_pair_subset_unpaired = g_subset,
      direction_matches_full_cohort =
        ifelse(is.finite(g_obs) & is.finite(full$hedges_g),
               sign(g_obs) == sign(full$hedges_g), NA),
      direction_matches_pair_subset =
        ifelse(is.finite(g_obs) & is.finite(g_subset),
               sign(g_obs) == sign(g_subset), NA),
      stringsAsFactors = FALSE
    )
    res <- res[order(res$p_perm_paired), , drop = FALSE]

    list(res = res, n_pairs = n_pairs, sets = sets,
         Jp = Jp, N_BOOT = N_BOOT, N_PERM = N_PERM)
  }

  stage2 <- lapply(seq_along(PAIRED_COHORTS), function(i)
    run_paired_cohort(PAIRED_COHORTS[i], i))
  names(stage2) <- PAIRED_COHORTS

  paired_tab <- do.call(rbind, lapply(stage2, `[[`, "res"))
  rownames(paired_tab) <- NULL

  # ------------------------------------------------------------- summary ----
  summarise_paired <- function(nm) {
    d <- stage2[[nm]]$res
    data.frame(
      cohort = nm,
      disease = d$disease[1],
      n_pairs_manifest = d$n_pairs_manifest[1],
      n_pairs_analysed = d$n_pairs_analysed[1],
      n_tierA_sets_tested = nrow(d),
      n_direction_matches_full_cohort = sum(d$direction_matches_full_cohort, na.rm = TRUE),
      n_direction_differs_from_full_cohort = sum(!d$direction_matches_full_cohort, na.rm = TRUE),
      n_direction_matches_pair_subset = sum(d$direction_matches_pair_subset, na.rm = TRUE),
      n_paired_ci_excludes_0 = sum(d$paired_ci_excludes_0, na.rm = TRUE),
      n_p_perm_nominal_lt_0.05 = sum(d$p_perm_paired < 0.05, na.rm = TRUE),
      n_p_fdr_lt_0.05 = sum(d$p_fdr_within_cohort < 0.05, na.rm = TRUE),
      median_g_paired_z = stats::median(d$g_paired_z, na.rm = TRUE),
      median_hedges_g_full_cohort = stats::median(d$hedges_g_full_cohort, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  paired_summary <- do.call(rbind, lapply(PAIRED_COHORTS, summarise_paired))

  methods_tab <- data.frame(
    parameter = c(
      "paired_unit", "effect_size", "standardisation_denominator",
      "hedges_correction", "bootstrap_unit", "bootstrap_replicates",
      "bootstrap_CI", "permutation_scheme", "shared_sign_matrix_across_sets",
      "permutation_replicates", "permutation_p", "multiplicity_family",
      "fdr_method", "unmatched_samples", "scope", "cross_disease_rho"
    ),
    value = c(
      "complete patient pair (one disease + one control sample)",
      "g_paired_z = J_p * mean(d_i) / sd(d_i); d_i = disease - control within pair",
      "SD of within-pair differences (NOT pooled between-group SD)",
      "J_p = 1 - 3/(4*n_pairs - 5)",
      "patient pairs, with replacement",
      as.character(N_BOOT),
      "percentile 95%",
      "within-pair disease/control label swap (sign flip of d_i)",
      "yes: one sign matrix per cohort reused for every Tier A set",
      as.character(N_PERM),
      "two-sided (1 + #{|stat_perm| >= |stat_obs|}) / (N_PERM + 1)",
      "Tier A sets within each paired cohort",
      STATS$FDR_METHOD,
      "excluded from paired analysis; never forced into a pair",
      "complete-pair subset sensitivity analysis only",
      "not computed"
    ),
    stringsAsFactors = FALSE
  )

  write_table(paired_tab, "Table_S40c_paired_effect_sizes.csv")
  write_table(paired_summary, "Table_S40d_paired_summary.csv")
  write_table(methods_tab, "Table_S40e_paired_methods.csv")

  saveRDS(list(results = paired_tab, summary = paired_summary,
               methods = methods_tab,
               checks = lapply(stage2, function(x)
                 list(n_pairs = x$n_pairs, sets = x$sets, Jp = x$Jp))),
          file.path(OUT, "intermediate", "R30_paired_stats.rds"))

  log_msg("Paired sensitivity summary:")
  for (i in seq_len(nrow(paired_summary))) {
    s <- paired_summary[i, ]
    log_msg(sprintf(
      "  %-18s pairs=%d/%d sets=%d dir_full=%d flip=%d dir_subset=%d CI!=0=%d nominal=%d FDR=%d",
      s$cohort, s$n_pairs_analysed, s$n_pairs_manifest, s$n_tierA_sets_tested,
      s$n_direction_matches_full_cohort, s$n_direction_differs_from_full_cohort,
      s$n_direction_matches_pair_subset, s$n_paired_ci_excludes_0,
      s$n_p_perm_nominal_lt_0.05, s$n_p_fdr_lt_0.05))
  }

  flips <- paired_tab[!is.na(paired_tab$direction_matches_full_cohort) &
                        !paired_tab$direction_matches_full_cohort, ,
                      drop = FALSE]
  if (nrow(flips)) {
    log_msg("Sign differs from the full-cohort independent-sample g:")
    for (i in seq_len(nrow(flips))) {
      log_msg(sprintf("  %s %s: g_full=%.3f g_paired=%.3f p_perm=%.4g",
                      flips$cohort[i], flips$set[i],
                      flips$hedges_g_full_cohort[i], flips$g_paired_z[i],
                      flips$p_perm_paired[i]))
    }
  } else {
    log_msg("No Tier A set changes sign relative to the full-cohort analysis.")
  }

  log_msg("=== R30 stage 2 complete: paired sensitivity outputs written ===")
}
