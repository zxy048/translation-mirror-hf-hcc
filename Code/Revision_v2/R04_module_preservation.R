# =============================================================================
# R04_module_preservation.R
# -----------------------------------------------------------------------------
# Are the co-expression modules preserved ACROSS diseases, or is module
# architecture itself disease-specific? This is the direct test of the
# manuscript's title claim and is now a primary analysis.
#
# modulePreservation() is run BOTH WAYS -- HF modules tested in the HCC data and
# HCC modules tested in the HF data -- so a claim of "not preserved" cannot be
# an artifact of choosing one direction. v1 asserted that cross-platform
# differences precluded this analysis; R2 #4 correctly rejected that argument,
# so it is done here.
#
# Standard reading of Zsummary: >10 strong preservation, 2-10 weak-moderate,
# <2 no evidence. Zsummary is a permutation Z, so it is reported alongside the
# Bonferroni log p-values that modulePreservation also computes.
#
# RUN TIME: ~5 h for both directions at 200 permutations on this machine.
# The raw `pres` object is therefore written to disk IMMEDIATELY after the call
# and before any reshaping: a bug in the formatting code must not be able to
# destroy the computation. (It did once. That is why this comment exists.)
#
# ACCESSOR SHAPE (verified, not assumed -- see 00_probe_pres2.R):
#   pres$preservation$Z[["ref.HF"]][["inColumnsAlsoPresentIn.HCC"]]  -> data.frame
#     rownames = module colour, cols = moduleSize, Zsummary.pres,
#     Zdensity.pres, Zconnectivity.pres, medianRank.pres
#   pres$preservation$observed[[...]] and $log.pBonf[[...]] are likewise
#   data.frames. $q is NULL: q-values are not computed by this call.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R04 start ===")
r2 <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))
datHF <- r2$datHF; datHCC <- r2$datHCC
stopifnot(identical(colnames(datHF), colnames(datHCC)))

multiExpr <- list(HF  = list(data = datHF),
                  HCC = list(data = datHCC))
colorList <- list(HF  = r2$HF$net$colors,
                  HCC = r2$HCC$net$colors)

log_msg("HF: ", nrow(datHF), " samples x ", ncol(datHF), " genes; ",
        length(unique(colorList$HF)), " modules")
log_msg("HCC: ", nrow(datHCC), " samples x ", ncol(datHCC), " genes; ",
        length(unique(colorList$HCC)), " modules")

# The formatting below runs for ~5 h of upstream compute, so it is exercisable
# on its own: `Rscript R04_module_preservation.R` with R04_SMOKE=1 runs 2
# permutations and writes SMOKE-suffixed outputs, letting the reshaping be
# tested end to end in ~3 minutes. SMOKE results are never used in the paper;
# the flag is logged and stamped into every table.
SMOKE  <- nzchar(Sys.getenv("R04_SMOKE"))
n_perm <- if (SMOKE) 2L else STATS$N_PERM_PRESERV
# Every artefact this script writes goes through tn(), so a smoke run cannot
# overwrite a real table with 2-permutation numbers.
tn <- function(x) if (SMOKE) paste0("SMOKE_", x) else x
if (SMOKE) log_msg("*** SMOKE MODE: nPermutations = 2, results NOT for the paper ***")

disableWGCNAThreads()
t0 <- Sys.time()
pres <- WGCNA::modulePreservation(
  multiExpr, colorList,
  referenceNetworks = c(1, 2),
  testNetworks      = list(2, 1),   # HF ref -> HCC test ; HCC ref -> HF test
  nPermutations     = n_perm,
  randomSeed        = STATS$SEED,
  networkType       = WGCNA$networkType,
  verbose           = 3, indent = 0)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
log_msg(sprintf("modulePreservation finished in %.1f min", elapsed))

# --- SAVE THE COMPUTATION FIRST ----------------------------------------------
raw_file <- file.path(OUT, "intermediate",
                      if (SMOKE) "R04_preservation_raw_SMOKE.rds"
                      else "R04_preservation_raw.rds")
saveRDS(list(pres = pres, nPermutations = n_perm,
             elapsed_min = elapsed, config = CONFIG_VERSION, smoke = SMOKE),
        raw_file)
log_msg("raw modulePreservation object saved: ", basename(raw_file))

# --- flatten: named accessors, keyed by the direction actually run ------------
DIRS <- list(
  list(ref = "ref.HF",  in_test = "inColumnsAlsoPresentIn.HCC",
       reference = "HF",  test = "HCC"),
  list(ref = "ref.HCC", in_test = "inColumnsAlsoPresentIn.HF",
       reference = "HCC", test = "HF")
)

#' Pull one preservation data.frame and tag it with the direction.
grab <- function(slot, d, what = "") {
  x <- pres$preservation[[slot]][[d$ref]][[d$in_test]]
  if (is.null(x)) {
    log_msg("  !! ", slot, " missing for ", d$reference, " -> ", d$test)
    return(NULL)
  }
  x <- as.data.frame(x)
  cbind(reference = d$reference, test = d$test, module = rownames(x), x,
        stringsAsFactors = FALSE, row.names = NULL)
}

Z  <- do.call(rbind, lapply(DIRS, function(d) grab("Z", d)))
OB <- do.call(rbind, lapply(DIRS, function(d) grab("observed", d)))
PB <- do.call(rbind, lapply(DIRS, function(d) grab("log.pBonf", d)))

# one combined table: the Z statistics, the observed statistics they summarise,
# and the Bonferroni-corrected log p-value for the summary statistic.
key <- function(d) paste(d$reference, d$test, d$module)
m <- merge(Z, OB, by = c("reference", "test", "module"), suffixes = c("", ".obs"))
pb_keep <- PB[, c("reference", "test", "module", "log.p.Bonfsummary.pres",
                  "log.p.Bonfdensity.pres", "log.p.Bonfconnectivity.pres")]
pres_tab <- merge(m, pb_keep, by = c("reference", "test", "module"), all.x = TRUE)
pres_tab <- pres_tab[pres_tab$module != "gold", ]   # "gold" = all-network pseudo-module

# Benjamini-Hochberg within each direction, on the summary p-value. The
# Bonferroni log p is what modulePreservation reports; converting to BH makes it
# comparable with the FDR convention used everywhere else in this revision.
pres_tab$p_summary <- 10^(pres_tab$log.p.Bonfsummary.pres)
pres_tab$p_summary_bh <- ave(pres_tab$p_summary,
                             interaction(pres_tab$reference, pres_tab$test),
                             FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
pres_tab$preservation_class <- cut(pres_tab$Zsummary.pres,
                                   breaks = c(-Inf, 2, 10, Inf),
                                   labels = c("none (<2)", "weak-moderate (2-10)",
                                              "strong (>10)"))
pres_tab <- pres_tab[order(pres_tab$reference, -pres_tab$Zsummary.pres), ]
write_table(pres_tab, tn("Table_S19_module_preservation.csv"))

cat("\n=== modulePreservation (Zsummary: >10 strong, 2-10 weak, <2 none) ===\n")
for (rn in unique(pres_tab$reference)) {
  d <- pres_tab[pres_tab$reference == rn, ]
  cat(sprintf("\n-- %s modules tested in the other disease (%d modules) --\n",
              rn, nrow(d)))
  print(head(d[, c("module", "moduleSize", "Zsummary.pres", "Zdensity.pres",
                   "Zconnectivity.pres", "medianRank.pres", "p_summary_bh")], 12),
        row.names = FALSE)
}
cat(sprintf("\nmodules with Zsummary > 10 : %d / %d\n",
            sum(pres_tab$Zsummary.pres > 10), nrow(pres_tab)))
cat(sprintf("modules with Zsummary <  2 : %d / %d\n",
            sum(pres_tab$Zsummary.pres < 2), nrow(pres_tab)))
cat(sprintf("modules with BH p < 0.05  : %d / %d\n",
            sum(pres_tab$p_summary_bh < 0.05, na.rm = TRUE), nrow(pres_tab)))

cat("\n=== the two designated translation modules, in both directions ===\n")
r3 <- readRDS(file.path(OUT, "intermediate", "R03_translation_module.rds"))
des <- rbind(data.frame(reference = "HF",  module = r3$desHF$module),
             data.frame(reference = "HCC", module = r3$desHCC$module))
dd <- merge(des, pres_tab, by = c("reference", "module"))
print(dd[, c("reference", "module", "test", "moduleSize", "Zsummary.pres",
             "medianRank.pres", "p_summary_bh")], row.names = FALSE)
write_table(dd, tn("Table_S20_designated_module_preservation.csv"))

# --- what the designated-module result means for the title claim --------------
# A module that is NOT preserved across diseases is the evidence for
# disease-context-dependent organisation; a preserved one would undercut it.
# Recorded explicitly so the manuscript cannot overstate either way.
des_summary <- do.call(rbind, lapply(seq_len(nrow(dd)), function(i) {
  z <- dd$Zsummary.pres[i]
  data.frame(reference = dd$reference[i], module = dd$module[i],
             test = dd$test[i], Zsummary = z,
             verdict = if (!is.finite(z)) "undefined"
                       else if (z > 10) "strongly preserved"
                       else if (z >= 2) "weakly preserved"
                       else "not preserved",
             stringsAsFactors = FALSE)
}))
print(des_summary, row.names = FALSE)
write_table(des_summary, tn("Table_S20b_designated_module_verdict.csv"))

# tn() here too. This saveRDS was the one artefact that bypassed it, so a smoke
# run silently wrote 2-permutation numbers to the production path -- exactly the
# trap the comment above tn() claims to close. Do not reintroduce a bare name.
saveRDS(list(table = pres_tab, designated = dd, verdict = des_summary,
             nPermutations = n_perm, elapsed_min = elapsed,
             raw_file = raw_file, smoke = SMOKE),
        file.path(OUT, "intermediate", tn("R04_preservation.rds")))
log_msg("=== R04 done ===")
