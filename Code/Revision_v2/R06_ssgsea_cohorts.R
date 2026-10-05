# =============================================================================
# R06_ssgsea_cohorts.R
# -----------------------------------------------------------------------------
# ssGSEA score per sample for every canonical set, in every cohort.
#
# Runs ONE ssGSEA per cohort over the UNION of Tier A, Tier B and Tier C (Tier C
# is v1's contaminated 33 and is kept only for continuity reporting), then
# subsets into the three tier views. Running the union once instead of three
# times is not just faster: it guarantees the three tiers are scored from the
# identical ranked matrix, so a Tier A vs Tier C difference can never be an
# artifact of a re-run.
#
# Gene sets are intersected with each cohort's measured genes. A set left with
# < minSize genes is dropped for that cohort and the drop is recorded -- v1
# never reported this, which is how a "33-pathway" analysis silently became a
# different number of pathways per cohort.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(GSVA))

log_msg("=== R06 start ===")
cohorts <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))

ALL_SETS <- unique(c(SETS_TIER_A, SETS_TIER_B_EXTRA, SETS_TIER_C_V1))
log_msg("unique canonical sets to score: ", length(ALL_SETS))

sets_all <- resolve_sets(ALL_SETS)

scores   <- list()
coverage <- list()

for (nm in names(cohorts)) {
  co <- cohorts[[nm]]
  gx <- co$expr
  log_msg("--- ", nm, " (", nrow(gx), " genes x ", ncol(gx), " samples) ---")

  sets <- lapply(sets_all, function(g) intersect(g, rownames(gx)))
  n_before <- vapply(sets, length, integer(1))
  keep <- n_before >= 5
  dropped <- names(sets)[!keep]
  sets <- sets[keep]
  if (length(dropped))
    log_msg("  dropped (", length(dropped), " sets, <5 genes measured): ",
            paste(dropped, collapse = ", "))

  p <- GSVA::ssgseaParam(as.matrix(gx), sets, minSize = 5, maxSize = 500,
                         normalize = TRUE)
  sc <- GSVA::gsva(p, verbose = FALSE)
  log_msg("  scored ", nrow(sc), " sets x ", ncol(sc), " samples")

  scores[[nm]] <- sc
  coverage[[nm]] <- data.frame(
    cohort = nm, set = names(sets_all), n_genes_in_cohort = n_before,
    scored = keep, stringsAsFactors = FALSE)
}

cov <- do.call(rbind, coverage)
write_table(cov, "Table_S12_set_coverage_per_cohort.csv")

# --- per-tier views -----------------------------------------------------------
tier_of <- function(set) {
  ifelse(set %in% SETS_TIER_A, "A",
  ifelse(set %in% SETS_TIER_B_EXTRA, "B_extra",
  ifelse(set %in% SETS_TIER_C_V1, "C_v1", NA_character_)))
}
sizes <- do.call(rbind, lapply(names(scores), function(nm) {
  s <- scores[[nm]]
  data.frame(cohort = nm, n_sets_total = nrow(s),
             n_tierA = sum(rownames(s) %in% SETS_TIER_A),
             n_tierB = sum(rownames(s) %in% c(SETS_TIER_A, SETS_TIER_B_EXTRA)),
             n_tierC = sum(rownames(s) %in% SETS_TIER_C_V1),
             stringsAsFactors = FALSE)
}))
cat("\n=== sets actually scored per cohort ===\n"); print(sizes, row.names = FALSE)
write_table(sizes, "Table_S13_scored_set_counts.csv")

saveRDS(list(scores = scores, coverage = cov, sizes = sizes,
             ALL_SETS = ALL_SETS, tier_of = tier_of(ALL_SETS)),
        file.path(OUT, "intermediate", "R06_ssgsea.rds"))
log_msg("=== R06 done ===")
