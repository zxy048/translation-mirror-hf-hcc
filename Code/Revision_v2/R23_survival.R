# =============================================================================
# R23_survival.R
# -----------------------------------------------------------------------------
# The translation score against overall survival -- the analysis the v2
# manuscript reports but that no v2 script performed.
#
# WHY THIS FILE EXISTS
# --------------------
# Manuscript_PLOSONE_v2.md makes a survival claim in the Results:
#
#   "TCS was not associated with overall survival in any cohort with usable
#    annotation, whether modelled as a continuous per-IQR variable or as a
#    median split: TCGA-LIHC (n = 421) and GSE141198 (n = 148), and no
#    association emerged in the aetiology-stratified analyses (S29 Table)."
#
# and S29's caption says it includes "the null survival associations".
#
# Neither was true of the pipeline. There is no Cox model, no log-rank test and
# no survival-time field anywhere in Code/Revision_v2/, and Table_S29 -- a
# disease-versus-control effect-size table with one row per score per cohort --
# has no survival column and does not contain GSE141198 at all (that cohort has
# no control arm, so it never enters a case-control table).
#
# This is not a detail missed in passing; it is the failure the v2 rebuild was
# commissioned to remove. The first submission shipped supplementary tables that
# no script could produce, and a binomial P-value with no implementation in the
# repository. A survival paragraph whose numbers cannot be regenerated is the
# same defect one revision later.
#
# The intent was recorded and then not carried out. In
# v2_output/notes/manuscript_v2_change_list.md the v1 survival section is listed
# as "[N][C] Recompute under TCS" -- and the score matters, because v1's numbers
# (GSE141198 log-rank P = 0.303, Cox HR = 1.68 [0.67, 4.21]) belong to the
# withdrawn 227-gene green-module TATS, computed by the v1-era script
# 13_TATS_clinical_association.R. TCS shares no gene list with that score: it is
# the mean z-scored ssGSEA over the 15 Tier A canonical sets. Carrying the old
# numbers forward, or asserting the old null without recomputing it, would be
# reporting v1's result under v2's name.
#
# SCOPE, AND WHY IT IS THIS WIDE
# ------------------------------
# The sentence says "any cohort with usable annotation", so the cohort list is
# not a choice -- it is whatever the deposited data support. Every cohort in the
# pipeline was checked:
#
#   HCC  TCGA-LIHC          colData of TCGA_LIHC_se.rds   vital_status + days
#   HCC  GSE141198          GSE141198_clinical.rds        os_time + os_status
#   HCC  GSE76427           series matrix characteristics event_os, duryears_os
#   HCC  GSE14520 (x2 GPL)  GSE14520_Extra_Supplement     Survival status/months
#   HF   GSE57338           -- no survival field deposited
#   HF   GSE116250          -- no survival field deposited
#   HF   GSE141910          -- no survival field deposited
#
# The three heart-failure cohorts carry only tissue/disease/sex/age/etiology, so
# the claim is scoped to the hepatic cohorts and the coverage block below prints
# that per cohort rather than leaving it as an assertion. GSE14520's two platforms
# are the same patients assayed twice and are reported as two rows, not as
# replication -- the pipeline treats them as separate sub-cohorts everywhere else
# for the same reason.
#
# Time units are NOT harmonised. TCGA and GSE141198 deposit days, GSE76427
# years, GSE14520 months. A Cox hazard ratio is invariant to a linear rescaling
# of the time axis, and every model here is fitted within one cohort, so
# rewriting months into days with a 30.44 multiplier would add a fabricated
# constant without changing a single estimate. The native unit is carried in a
# `time_unit` column instead, where a reader can see it.
#
# ADJUSTMENT, AND WHY PURITY AND PROLIFERATION ARE NOT THE SAME COVARIATE
# ----------------------------------------------------------------------
# The unadjusted test is positive in GSE14520, so the first question is what the
# score is standing in for. The two candidates are tumour purity and
# proliferation, and they do not carry the same meaning:
#
#   purity         a MEASUREMENT confounder. A bulk tumour transcriptome is a
#                  mixture of tumour and stroma, and a score built from highly
#                  expressed growth genes moves with the tumour-cell fraction
#                  whether or not it says anything about the tumour. Adjusting
#                  for it is required. The hepatocyte rank score is the proxy,
#                  the same one R11 uses for HCC.
#   proliferation  a MEDIATOR, on this paper's own account. The manuscript's
#                  thesis is that the HCC translation program is proliferation-
#                  coupled, so part of any prognostic signal is expected to run
#                  THROUGH proliferation. Adjusting for it therefore does not
#                  correct a bias; it answers a different and stricter question
#                  -- whether TCS carries prognostic information beyond cell
#                  cycling. A null under this adjustment is informative, but it
#                  is not evidence that the unadjusted association was spurious.
#
# Both are reported, separately and together with immune content, and the
# `covariates` column labels them so the two cannot be read as one.
#
# `p_logrank` is the UNadjusted log-rank on the median split and is repeated down
# every adjusted row of a cohort by design: survdiff has no covariate term, and
# the log-rank is reported as the model-free counterpart to the Cox estimate, not
# as a test of the adjusted model. `p_cox` is the one that moves with `covariates`.
#
# The covariates are the ESTIMATE-style within-sample rank scores from
# R11_composition.rds: no new data, no new package, and the same primary set
# (hepatocyte + proliferation + immune_pan) R11 already froze for HCC. They are
# z-scored before entry, as R11 does, so the TCS hazard ratio stays per-IQR and
# the covariate hazard ratios are per-SD.
#
# Every adjusted model is fitted on the samples the unadjusted model used. The
# composition scores carry no missing values in these cohorts, so the complete-
# case set is the full set and no cohort silently loses samples to the
# adjustment -- which is checked in the code rather than assumed, because a
# shorter n under adjustment reads as attenuation that is really attrition.
#
# The events-per-term floor is applied per model, at the conventional ten events
# per covariate. A three-covariate fit therefore needs 40 events, which keeps
# GSE76427 (23) out of the fully adjusted model while leaving it in the one-
# covariate fits. That exclusion is reported, not silent.
#
# WHAT THIS SCRIPT DOES NOT DO
# ----------------------------
# It does not adjust for stage or grade, and that is a decision rather than an
# omission. Stage is a mediator as much as a confounder here: TCS is
# proliferation-coupled and stage is downstream of proliferation, so adjusting
# for it would remove part of the exposure's own effect and turn a null into a
# foregone conclusion. Composition is adjusted for, above; stage and grade are
# not, and the limitation is recorded in the manuscript.
#
# It does not run event-free survival for GSE141198, which deposits efs_time and
# efs_status. The manuscript's claim is about overall survival; adding a second
# endpoint would be a new analysis, not a recomputation of the one being fixed.
#
# No new dataset is introduced. Every file read here is already in the repository
# and is already read by another v2 script.
#
# Run:  Rscript Code/Revision_v2/R23_survival.R
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(survival)
})

log_msg("=== R23 start ===")

r5  <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
r12 <- readRDS(file.path(OUT, "intermediate", "R12_tats_tf.rds"))
r11 <- readRDS(file.path(OUT, "intermediate", "R11_composition.rds"))
SCORES <- c("TCS", "ME_HF_weighted", "ME_HCC_weighted")

# hepatocyte is the purity proxy, proliferation the mediator, immune_pan the
# shared confounder; R11 froze exactly this vector as the HCC primary covariate
# set, and it is reused rather than re-chosen here so the survival adjustment and
# the case-control adjustment cannot disagree about what "adjusted for
# composition" means. The unadjusted fit is the default `covar_label = "none"`.
COV_SETS <- list(
  purity                        = "hepatocyte",
  proliferation                 = "proliferation",
  "purity+proliferation+immune" = r11$primary_covars$HCC)

#' VIF of each column of a numeric design matrix. >10 means the adjustment is
#' collinear and the coefficient should be distrusted rather than reported. The
#' exposure is included alongside the covariates because collinearity between TCS
#' and the covariate is what inflates the TCS standard error, not collinearity
#' among the covariates alone.
vif_of <- function(X) {
  cn <- colnames(X)
  if (length(cn) < 2) return(stats::setNames(rep(1, length(cn)), cn))
  vapply(cn, function(v) {
    r2 <- summary(stats::lm(X[, v] ~ X[, setdiff(cn, v), drop = FALSE]))$r.squared
    if (!is.finite(r2) || r2 >= 1) Inf else 1 / (1 - r2)
  }, numeric(1))
}

# =============================================================================
# 1. Survival annotation, one loader per source, keyed to the cohort matrix
# =============================================================================
# Every loader keys on the sample ID the cohort matrix itself uses. A positional
# join between an expression matrix and a clinical table is the mistake that
# produces a plausible result from mismatched samples, and it is not detectable
# from the output -- R10c keys TCGA on the barcode for the same reason.
#
# Each loader returns sample / time / event / stratum, and the caller matches
# with `match()` and reports the unmatched count. An unmatched sample is dropped
# loudly; a silently recycled shorter vector would shift every patient's
# survival onto their neighbour's expression.
num_col <- function(x) suppressWarnings(as.numeric(as.character(x)))

load_surv <- list()

# --- TCGA-LIHC: primary tumours only ----------------------------------------
# colData carries a vital_status for the solid-tissue-normal aliquots too (they
# belong to the same patients), so an unrestricted join would enter 50 control
# samples into an overall-survival analysis and report n = 421 -- the cohort
# size, not the analysable n. That is where the manuscript's "n = 421" comes
# from, and it is the tumour-arm n that belongs there.
se <- readRDS(file.path(DATA_RAW, "TCGA_LIHC_se.rds"))
cd <- as.data.frame(colData(se))
cd <- cd[intersect(colnames(r5$TCGA_LIHC$expr), rownames(cd)), , drop = FALSE]
cd <- cd[cd$sample_type == "Primary Tumor", , drop = FALSE]
load_surv$TCGA_LIHC <- data.frame(
  sample = rownames(cd),
  time = num_col(ifelse(cd$vital_status == "Dead", cd$days_to_death,
                        cd$days_to_last_follow_up)),
  event = cd$vital_status == "Dead", stratum = NA_character_,
  time_unit = "days", source = "TCGA_LIHC_se.rds colData: vital_status, days_to_death / days_to_last_follow_up",
  stringsAsFactors = FALSE)

# --- GSE141198: GEO OS fields, plus aetiology for the subgroup models --------
cl <- readRDS(file.path(PROJ, "GSE141198_clinical.rds"))
load_surv$GSE141198 <- data.frame(
  sample = cl$geo_id, time = num_col(cl$os_time), event = cl$os_status == 1,
  stratum = cl$etiology, time_unit = "days",
  source = "GSE141198_clinical.rds: os_time, os_status, etiology",
  stringsAsFactors = FALSE)

# --- GSE76427: OS in the series-matrix characteristics -----------------------
sm <- parse_series_matrix(file.path(DATA_RAW, "GSE76427_series_matrix.txt.gz"))$ann
load_surv$GSE76427 <- data.frame(
  sample = sm$gsm, time = num_col(sm$duryears_os), event = num_col(sm$event_os) == 1,
  stratum = NA_character_, time_unit = "years",
  source = "GSE76427_series_matrix.txt.gz: event_os, duryears_os",
  stringsAsFactors = FALSE)

# --- GSE14520: survival lives in the series' own supplementary table ---------
# The series matrices for this accessions carry only Tissue / Disease state /
# Individual, so the survival columns are read from GSE14520_Extra_Supplement,
# which v1 also used ("the OS definitions followed the GEO record annotations").
# Rows are patients; each patient is assayed on both platforms, so the same
# patient's survival enters the GPL3921 row via Affy_GSM and the GPL571 row via
# Agilent_GSM. Only tumour rows are kept.
sup <- utils::read.delim(gzfile(file.path(DATA_RAW, "GSE14520_Extra_Supplement.txt.gz")),
                         check.names = FALSE, stringsAsFactors = FALSE)
sup <- sup[trimws(sup[["Tissue Type"]]) == "Tumor", , drop = FALSE]
sup_t <- num_col(sup[["Survival months"]])
sup_e <- trimws(sup[["Survival status"]])
for (gpl in c("GPL3921", "GPL571")) {
  key <- if (gpl == "GPL3921") "Affy_GSM" else "Agilent_GSM"
  nm <- paste0("GSE14520_", gpl)
  load_surv[[nm]] <- data.frame(
    sample = trimws(sup[[key]]), time = sup_t, event = sup_e == "1",
    stratum = NA_character_, time_unit = "months",
    source = paste0("GSE14520_Extra_Supplement.txt.gz: ", key,
                    " -> Survival months, Survival status"),
    stringsAsFactors = FALSE)
}

# =============================================================================
# 2. Attach to the pipeline's own samples and report what could not be attached
# =============================================================================
# The coverage table is the answer to "any cohort with usable annotation". It
# prints, per cohort, how many samples joined, how many were dropped for a
# missing or non-positive time, and how many events remain -- so a reader can
# check the scope of the claim rather than take it on trust.
build <- list(); report <- list()
for (nm in names(r5)) {
  ids <- colnames(r5[[nm]]$expr)
  sv <- load_surv[[nm]]
  if (is.null(sv)) {
    log_msg("  ", nm, ": NO survival annotation in the deposited data -- not analysed")
    report[[length(report) + 1]] <- data.frame(
      cohort = nm, disease = r5[[nm]]$disease, analysed = FALSE,
      n_samples = length(ids), n_joined = 0L, n_dropped_time = 0L,
      n_analysed = 0L, n_events = 0L, time_unit = NA_character_,
      note = "no survival field deposited in this series",
      source = "none deposited", stringsAsFactors = FALSE)
    next
  }
  m <- match(ids, sv$sample)
  d <- data.frame(sample = ids, time = sv$time[m], event = sv$event[m],
                  stratum = sv$stratum[m], stringsAsFactors = FALSE)
  ok <- !is.na(d$time) & d$time > 0 & !is.na(d$event)
  d <- d[ok, , drop = FALSE]
  # A cohort is "analysed" only if something survived the join. An annotation
  # file that exists but whose sample IDs do not intersect this cohort is not
  # usable annotation, and calling it analysed with n = 0 would put the claim
  # "every cohort with usable annotation was tested" one cohort wider than the
  # computation. GSE14520_GPL571 is exactly that case: the supplement's Agilent
  # block is GSM358026-GSM358100, the pipeline's GPL571 cohort is
  # GSM362947-GSM363449, and the two do not intersect.
  if (nrow(d)) build[[nm]] <- d
  note <- if (!nrow(d) && !is.na(sum(!is.na(m))))
    "annotation present but no sample ID in this cohort matches" else ""
  report[[length(report) + 1]] <- data.frame(
    cohort = nm, disease = r5[[nm]]$disease, analysed = nrow(d) > 0,
    n_samples = length(ids), n_joined = sum(!is.na(m)),
    n_dropped_time = sum(!is.na(m)) - sum(ok),
    n_analysed = nrow(d), n_events = sum(d$event),
    time_unit = sv$time_unit[1], note = note, source = sv$source[1],
    stringsAsFactors = FALSE)
}
coverage <- do.call(rbind, report)
log_msg("--- survival annotation coverage ---")
print(coverage[, c("cohort", "disease", "analysed", "n_samples", "n_analysed",
                   "n_events", "time_unit")], row.names = FALSE)

# =============================================================================
# 3. Models
# =============================================================================
#' One score in one cohort: per-IQR Cox, median-split Cox, log-rank on the split.
#' Returns one row per model, so the table carries a `model` column and a reader
#' can see which quantity each P-value belongs to. A wide table with four
#' P-values and no such column is how a log-rank P ends up quoted as a Cox P.
fit_surv <- function(df, score, cohort, model, stratum = NA_character_,
                     covars = character(0), covar_label = "none") {
  d <- data.frame(time = df$time, event = df$event, x = df[[score]])
  if (length(covars)) d <- cbind(d, as.data.frame(df[, covars, drop = FALSE]))
  d <- d[stats::complete.cases(d) & is.finite(d$x) & is.finite(d$time) &
           d$time > 0, , drop = FALSE]
  # The guard is on EVENTS per TERM: a Cox coefficient is informed by how many
  # events it has, at a conventional floor of ten per covariate. With no
  # covariates the floor is ten events, which is what keeps the NBNC subgroup
  # (n = 19, 10 events) in the aetiology block -- v1 reported all three strata,
  # and dropping one silently would narrow the manuscript's "aetiology-stratified
  # analyses" without saying so. With three covariates it is forty, which is what
  # keeps the fully adjusted fit out of GSE76427 (23 events). n and n_events are
  # both in the table so a reader can apply a different rule.
  k <- 1 + length(covars)
  if (nrow(d) < 10 || sum(d$event) < 10 * k) return(NULL)

  if (model == "continuous_per_IQR") {
    iqr <- stats::IQR(d$x)
    if (!is.finite(iqr) || iqr <= 0) return(NULL)
    d$v <- d$x / iqr                     # the HR is per one IQR of the score
  } else {
    d$v <- as.numeric(d$x > stats::median(d$x))
  }
  # Covariates z-scored, as R11 does. The rank scores run into the thousands, so
  # an unscaled one would give a hazard ratio per rank-unit: unreadable, and not
  # comparable to the per-IQR exposure. Scaling the covariate leaves the exposure
  # term untouched.
  if (length(covars)) d[covars] <- lapply(d[covars], function(v) {
    s <- stats::sd(v); if (!is.finite(s) || s == 0) v - mean(v) else (v - mean(v)) / s })

  form <- stats::as.formula(paste("Surv(time, event) ~",
                                  paste(c("v", covars), collapse = " + ")))
  fit <- tryCatch(survival::coxph(form, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  ph <- tryCatch(survival::cox.zph(fit)$table["GLOBAL", "p"], error = function(e) NA_real_)
  lr <- if (model == "median_split")
    tryCatch(stats::pchisq(survival::survdiff(survival::Surv(time, event) ~ v,
                                              data = d)$chisq, df = 1,
                           lower.tail = FALSE), error = function(e) NA_real_) else NA_real_
  # Covariate hazard ratios travel with the row as text. In an adjusted model the
  # most interpretable number is often the covariate's own HR -- does
  # proliferation carry the signal TCS was credited with? -- and a table that
  # reports only the adjusted exposure effect hides the answer.
  hr_cov <- if (length(covars))
    paste(sprintf("%s %.2f", covars, s$conf.int[covars, "exp(coef)"]),
          collapse = "; ") else ""
  data.frame(
    cohort = cohort, stratum = stratum, score = score, model = model,
    covariates = covar_label, n_covariates = length(covars),
    n = nrow(d), n_events = sum(d$event),
    hr = unname(s$conf.int["v", "exp(coef)"]),
    ci_lo = unname(s$conf.int["v", "lower .95"]),
    ci_hi = unname(s$conf.int["v", "upper .95"]),
    p_cox = unname(s$coefficients["v", "Pr(>|z|)"]),
    p_logrank = lr, p_ph_global = ph,
    # Flagged rather than left to the reader to threshold, because a proportional-
    # hazards violation changes what the P-value means and the affected models are
    # not the same in every cohort. In TCGA the TCS-adjusted models are the ones
    # that fail (P = 0.0009 to 0.020), so their P-values are reported but must not
    # be relied on; in GSE14520, where the association is, PH holds (min P = 0.20).
    ph_ok = is.na(ph) | ph >= 0.05,
    max_vif = if (length(covars))
      max(vif_of(as.matrix(d[, c("v", covars), drop = FALSE]))) else 1,
    hr_covariates = hr_cov, stringsAsFactors = FALSE)
}

res <- list()
for (nm in names(build)) {
  d <- build[[nm]]
  sc <- r12$scores[[nm]]
  if (is.null(sc)) { log_msg("  ", nm, ": no scores -- skipped"); next }
  # merge() rather than a positional cbind: `sc` is keyed by its own rownames and
  # the two objects were built in different scripts, so their column order is not
  # guaranteed to agree.
  d <- merge(d, sc[, SCORES, drop = FALSE], by.x = "sample", by.y = 0,
             all.x = FALSE, sort = FALSE)

  # Composition covariates, keyed the same way. A sample lost here would make the
  # adjusted model a smaller cohort than the unadjusted one, and the reader would
  # see n fall with no explanation -- so the loss is reported rather than left to
  # be discovered in the n column.
  cp <- r11$scores[[nm]]
  if (!is.null(cp)) {
    want <- intersect(unique(unlist(COV_SETS)), colnames(cp))
    before <- nrow(d)
    d <- merge(d, cp[, want, drop = FALSE], by.x = "sample", by.y = 0,
               all.x = TRUE, sort = FALSE)
    lost <- sum(!stats::complete.cases(d[, want, drop = FALSE]))
    if (lost) log_msg("  NOTE ", nm, ": ", lost, " of ", before,
                      " samples carry no composition score -- adjusted models use fewer samples")
  }

  for (s in SCORES) for (mo in c("continuous_per_IQR", "median_split"))
    res[[length(res) + 1]] <- fit_surv(d, s, nm, mo)

  # Adjusted models for TCS only. The two module-eigengene scores are the
  # sensitivity arm of the survival block, and running three scores across three
  # covariate sets would triple the table to answer a question that is being
  # asked of the primary score. A covariate set that fails the events-per-term
  # floor returns NULL and is logged here rather than dropped in silence.
  if (!is.null(cp)) for (cl in names(COV_SETS)) {
    cv <- intersect(COV_SETS[[cl]], colnames(d))
    if (length(cv) < length(COV_SETS[[cl]])) {
      log_msg("  NOTE ", nm, ": covariate set '", cl, "' not fitted -- ",
              paste(setdiff(COV_SETS[[cl]], cv), collapse = ", "),
              " not scored in this cohort")
      next
    }
    for (mo in c("continuous_per_IQR", "median_split")) {
      row <- fit_surv(d, "TCS", nm, mo, covars = cv, covar_label = cl)
      if (is.null(row))
        log_msg("  ", nm, " / TCS / ", mo, " / ", cl, ": NOT fitted (needs ",
                10 * (1 + length(cv)), " events, has ", sum(d$event), ")")
      else res[[length(res) + 1]] <- row
    }
  }

  # Aetiology subgroups, where annotated. v1 ran these as subgroups rather than
  # as an adjusted model, and the distinction is kept: an adjustment term
  # estimates one common within-stratum slope, whereas the claim under test is
  # that no subgroup shows an association at all.
  for (g in unique(stats::na.omit(d$stratum))) {
    sub <- d[!is.na(d$stratum) & d$stratum == g, , drop = FALSE]
    row <- fit_surv(sub, "TCS", nm, "continuous_per_IQR", g)
    if (is.null(row)) {
      log_msg("  ", nm, " stratum ", g, ": NOT reported -- n = ", nrow(sub),
              ", events = ", sum(sub$event, na.rm = TRUE),
              " (below the >= 10-event floor)")
      next
    }
    res[[length(res) + 1]] <- row
  }
}

surv_tab <- do.call(rbind, res)
surv_tab <- merge(surv_tab, coverage[, c("cohort", "time_unit")], by = "cohort",
                  all.x = TRUE, sort = FALSE)
# BH within cohort over the reported tests, matching every other table here
# (STATS$FDR_METHOD). It is reported for completeness, not for inference: the
# primary comparison is TCS in the two largest cohorts, and an FDR across a
# score x model x subgroup grid answers a question nobody asked.
surv_tab$p_fdr <- ave(surv_tab$p_cox, surv_tab$cohort,
                      FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
surv_tab <- surv_tab[order(surv_tab$cohort, surv_tab$stratum, surv_tab$score,
                           surv_tab$model, surv_tab$covariates), ]
write_table(surv_tab, "Table_S39_translation_score_survival.csv")

# The coverage block is a table in its own right rather than a diagnostic. It is
# the evidence for "every cohort with usable annotation", and it is what shows the
# three heart-failure cohorts are absent because their series deposit no survival
# field, not because they were left out. It carries its own S-number because
# R19_build_SI_pdf.py renders only the first file matching `Table_<label>_*.csv`,
# so two files under one caption would silently drop the second.
write_table(coverage, "Table_S39b_survival_annotation_coverage.csv")

log_msg("--- TCS against overall survival ---")
print(surv_tab[surv_tab$score == "TCS",
               c("cohort", "stratum", "model", "covariates", "n", "n_events",
                 "hr", "ci_lo", "ci_hi", "p_cox", "p_logrank", "p_ph_global")],
      row.names = FALSE)

# The adjustment block on its own, unadjusted first within each cohort, because
# the question it answers is whether the unadjusted estimate moves -- and that is
# read by comparing rows, which a table sorted any other way obscures.
log_msg("--- TCS, continuous per-IQR, by covariate set (whole-cohort models) ---")
pr <- surv_tab[surv_tab$score == "TCS" &
               surv_tab$model == "continuous_per_IQR" &
               is.na(surv_tab$stratum), ]
pr <- pr[order(pr$cohort, match(pr$covariates, c("none", names(COV_SETS)))), ]
print(pr[, c("cohort", "covariates", "n", "n_events", "hr", "ci_lo", "ci_hi",
             "p_cox", "ph_ok", "max_vif", "hr_covariates")], row.names = FALSE)

log_msg("--- any TCS test below P = 0.05? ---")
tcs <- surv_tab[surv_tab$score == "TCS", ]
for (lab in c("none", names(COV_SETS))) {
  z <- tcs[tcs$covariates == lab, ]
  if (!nrow(z)) next
  cat(sprintf("  %-28s %2d tests; smallest P = %.4f at %s / %s\n", lab, nrow(z),
              min(z$p_cox), z$cohort[which.min(z$p_cox)],
              z$model[which.min(z$p_cox)]))
}

saveRDS(list(table = surv_tab, coverage = coverage, surv = build),
        file.path(OUT, "intermediate", "R23_survival.rds"))
log_msg("=== R23 done ===")
