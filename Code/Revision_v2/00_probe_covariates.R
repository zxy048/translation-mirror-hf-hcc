# =============================================================================
# 00_probe_covariates.R  (diagnostic, not part of the numbered pipeline)
# -----------------------------------------------------------------------------
# R2 Major 3 asks for adjustment for clinical covariates. Before writing any
# script that claims to adjust for age and sex, establish which cohorts actually
# carry those fields, and in what form. Print the annotation columns of every
# cohort in the frozen R05 object; do not assume from the GEO landing page.
#
# Output goes to a file, not the console: the Windows console is GBK and the
# series-matrix headers can be non-ASCII.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
sink(file.path(OUT, "logs", "00_probe_covariates.txt"))

cat("cohorts:", paste(names(r5), collapse = ", "), "\n\n")
for (nm in names(r5)) {
  co <- r5[[nm]]
  cat(strrep("=", 70), "\n", nm, "  n =", ncol(co$expr), " genes =", nrow(co$expr), "\n", sep = "")
  cat("tech:", co$tech, " subset:", co$subset_label, "\n")
  a <- co$ann
  if (is.null(a)) { cat("  ann: NULL\n\n"); next }
  cat("  ann dim:", nrow(a), "x", ncol(a), "\n")
  cat("  ann cols:", paste(colnames(a), collapse = " | "), "\n")
  for (cl in colnames(a)) {
    v <- a[[cl]]
    if (is.numeric(v)) {
      cat(sprintf("    %-18s numeric  n_nonNA=%d  range=[%g, %g]  median=%g\n",
                  cl, sum(!is.na(v)), min(v, na.rm = TRUE), max(v, na.rm = TRUE),
                  median(v, na.rm = TRUE)))
    } else {
      tt <- sort(table(v, useNA = "ifany"), decreasing = TRUE)
      show <- paste(sprintf("%s(%d)", names(tt), as.integer(tt))[seq_len(min(8, length(tt)))],
                    collapse = ", ")
      cat(sprintf("    %-18s factor   n_levels=%d  top: %s\n", cl, length(tt), show))
    }
  }
  cat("\n")
}
sink()
cat("written:", file.path(OUT, "logs", "00_probe_covariates.txt"), "\n")
