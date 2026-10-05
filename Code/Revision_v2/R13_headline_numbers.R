# =============================================================================
# R13_headline_numbers.R
# -----------------------------------------------------------------------------
# One place that prints every number the manuscript v2 and the Response to
# Reviewers quote, read from the saved R07-R12 objects rather than retyped.
#
# This exists because v1's central failure was numeric: the text asserted
# quantities that no script computed ("binomial test P < 0.0001", "all 33
# pathways positive in HCC"). Every figure that goes into v2 must be
# traceable to an object on disk, so they are assembled here and written to
# one file that the writing stage reads.
#
# DELIBERATELY DEFENSIVE. This script does not assume the shape of the R08/R09/
# R10/R12 objects -- it reports what each one actually contains and then prints
# the parts it recognises. A summary that dies on a guessed column name would
# be worse than no summary, because the failure would look like missing results.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R13 start ===")
rd <- function(f) readRDS(file.path(OUT, "intermediate", f))
obj <- list()
for (f in c("R03_translation_module.rds", "R06_ssgsea.rds", "R07_effect_sizes.rds",
            "R08_redundancy.rds", "R09_mirror.rds", "R10_heterogeneity.rds",
            "R11_composition.rds", "R12_tats_tf.rds", "R10b_hf_aetiology.rds")) {
  k <- sub("\\.rds$", "", f)
  obj[[k]] <- tryCatch(rd(f), error = function(e) {
    log_msg("  !! could not read ", f, ": ", conditionMessage(e)); NULL })
}
r6 <- obj$R06_ssgsea; r7 <- obj$R07_effect_sizes
r8 <- obj$R08_redundancy; r9 <- obj$R09_mirror
r10 <- obj$R10_heterogeneity; r11 <- obj$R11_composition
r12 <- obj$R12_tats_tf; r10b <- obj$R10b_hf_aetiology

con <- file(file.path(OUT, "tables", "v2_headline_numbers.txt"), open = "wt")
sink(con, split = TRUE)   # on screen AND to the file the writing stage reads
on.exit({ sink(); close(con) }, add = TRUE)

L <- function(...) cat(..., "\n", sep = "")
H <- function(x) { L(""); L("### ", x); L(strrep("-", 78)) }
# describe an object we have not assumed the shape of
describe <- function(x, nm) {
  if (is.null(x)) { L("  ", nm, ": <absent>"); return(invisible()) }
  if (is.data.frame(x)) L(sprintf("  %-22s data.frame %d x %d", nm, nrow(x), ncol(x)))
  else if (is.list(x))  L(sprintf("  %-22s list: %s", nm,
                                  paste(names(x) %||% seq_along(x), collapse = ", ")))
  else L(sprintf("  %-22s %s", nm, class(x)[1]))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
show_df <- function(x, nm, maxrows = 90) {
  if (is.null(x) || !is.data.frame(x)) return(invisible())
  L(""); L("-- ", nm, " --")
  print(utils::head(as.data.frame(x), maxrows), row.names = FALSE)
  if (nrow(x) > maxrows) L("  ... (", nrow(x) - maxrows, " more rows)")
}
# A wide table is still worth reporting if it has few COLUMNS: the heterogeneity
# table is 82 x 10 and is exactly what the response letter needs. Truncating on
# row count alone silently dropped it.
show_wide <- function(x, nm) {
  if (is.null(x) || !is.data.frame(x)) return(invisible())
  L(""); L("-- ", nm, " -- (", nrow(x), " x ", ncol(x), ")")
  if (ncol(x) <= 12) print(as.data.frame(x), row.names = FALSE)
  else print(utils::head(as.data.frame(x)[, seq_len(10)], 20), row.names = FALSE)
}

# ---------------------------------------------------------------- 01 cohorts --
H("01  cohort inventory")
if (!is.null(r6$scores)) {
  L("  cohort              samples  sets  tierA")
  for (nm in names(r6$scores)) {
    S <- r6$scores[[nm]]
    L(sprintf("  %-18s %6d %5d %6d", nm, ncol(S), nrow(S),
              sum(rownames(S) %in% SETS_TIER_A)))
  }
} else describe(r6, "R06 object")

# ------------------------------------------------------------ 02 effect sizes --
H("02  Tier A effect sizes per cohort (R07)")
if (!is.null(r7) && is.data.frame(r7) && "tier" %in% names(r7)) {
  a7 <- r7[r7$tier == "A", ]
  L("  cohort            n_sets  |g|>=0.20  CI excl 0     mean_g     min_g     max_g")
  for (nm in unique(a7$cohort)) {
    d <- a7[a7$cohort == nm, ]
    L(sprintf("  %-18s %5d %10d %9d  %+9.3f %+9.3f %+9.3f", nm, nrow(d),
              sum(abs(d$hedges_g) >= STATS$MIN_EFFECT),
              sum(d$ci_excludes_0, na.rm = TRUE),
              mean(d$hedges_g), min(d$hedges_g), max(d$hedges_g)))
  }
  L("")
  L("  cols: ", paste(colnames(r7), collapse = ", "))
} else { L("  R07 shape:"); describe(r7, "R07") }

# -------------------------------------------------------------- 03 redundancy --
H("03  redundancy / effective tests (R08)")
if (is.list(r8)) for (nm in names(r8)) {
  x <- r8[[nm]]
  if (is.data.frame(x) && nrow(x) <= 25) show_df(x, nm)
  else if (is.numeric(x) && length(x) <= 25) { L("  ", nm, ": ", paste(round(x, 3), collapse = ", ")) }
  else describe(x, nm)
} else describe(r8, "R08 object")

# ------------------------------------------------------------- 04 mirror test --
H("04  the mirror test (R09)  --  the manuscript's central claim")
if (is.list(r9)) for (nm in names(r9)) {
  x <- r9[[nm]]
  if (is.data.frame(x) && nrow(x) <= 40) show_df(x, nm)
  else describe(x, nm)
} else describe(r9, "R09 object")

# --------------------------------------------------------- 05 heterogeneity ---
H("05  between-cohort heterogeneity (R10)")
if (is.list(r10)) for (nm in names(r10)) {
  x <- r10[[nm]]
  if (is.data.frame(x)) show_wide(x, nm) else describe(x, nm)
} else describe(r10, "R10 object")

# ------------------------------------------------------- 06 composition adj ---
H("06  composition adjustment, Tier A (R11)")
if (!is.null(r11$adjustment)) {
  a <- r11$adjustment; a <- a[a$tier == "A", ]
  L(sprintf("  set-cohort pairs %d; survive %d (%.1f%%); median attenuation %.1f%%; max VIF %.2f",
            nrow(a), sum(a$survives), 100 * mean(a$survives),
            100 * median(a$attenuation, na.rm = TRUE), max(a$max_vif, na.rm = TRUE)))
  for (nm in unique(a$cohort))
    L(sprintf("    %-18s %2d / %2d survive", nm,
              sum(a$survives[a$cohort == nm]), sum(a$cohort == nm)))
} else describe(r11, "R11 object")

# -------------------------------------------------------- 07 aetiology (R10b) --
H("07  HF aetiology stratification (R10b)")
show_df(r10b$vs_hcc, "each HF subgroup vs TCGA-LIHC")
if (!is.null(r10b$dcm_cross)) {
  L("")
  L(sprintf("  DCM vs DCM across cohorts: rho = %+.3f; same direction %d / %d",
            r10b$rho_dcm_cross,
            sum(sign(r10b$dcm_cross$g_GSE57338_DCM) ==
                  sign(r10b$dcm_cross$g_GSE141910_DCM)), nrow(r10b$dcm_cross)))
}

# ------------------------------------------------------------ 08 TATS and TF --
H("08  TATS replacement and TF models (R12)")
if (is.list(r12)) for (nm in names(r12)) {
  x <- r12[[nm]]
  if (is.data.frame(x)) { L(""); L("-- ", nm, " -- ", nrow(x), " rows"); print(utils::head(x, 8), row.names = FALSE) }
  else describe(x, nm)
} else describe(r12, "R12 object")

# ------------------------------------------- 09 v1 numbers that must not return
H("09  v1 claims that v2 must NOT reproduce")
L("  'binomial test P < 0.0001' for 24/33 -> never implemented in any script.")
L("     The two-sided exact binomial for 24/33 at p=0.5 is P = ",
  formatC(stats::binom.test(24, 33, 0.5)$p.value, format = "f", digits = 4), ".")
L("  'all 33 translation-related pathways showed positive effect sizes in HCC'")
L("     -> tables/Table_S2_ssGSEA_Effect_Sizes.csv contains 6 negative values.")
r3 <- obj$R03_translation_module
if (!is.null(r3)) { L(""); L("  v2 designated modules (R03):"); describe(r3, "R03 object") }

log_msg("=== R13 done ===")
