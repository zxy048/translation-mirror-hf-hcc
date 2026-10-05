# =============================================================================
# 00_probe_covariates2.R  (diagnostic)
# -----------------------------------------------------------------------------
# R05 passed `ann = NULL` for GSE116250 and TCGA_LIHC, so the frozen cohort
# object does not expose their clinical annotation even though the underlying
# files carry it. Before deciding which cohorts an age/sex sensitivity analysis
# can cover, look at the annotation as the parser actually returns it -- the
# question is not "does GEO list an age column" but "does the column survive
# parsing with values, and for how many of the samples that enter the model".
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
sink(file.path(OUT, "logs", "00_probe_covariates2.txt"))

cat("### GSE116250 series matrix annotation\n")
sm <- parse_series_matrix(file.path(VALID_DIR, "GSE116250_series_matrix.txt.gz"))
cat("cols:", paste(colnames(sm$ann), collapse = " | "), "\n")
for (cl in colnames(sm$ann)) {
  v <- sm$ann[[cl]]
  if (is.numeric(v)) cat(sprintf("  %-14s numeric n_nonNA=%d range=[%g,%g]\n", cl,
                                 sum(!is.na(v)), min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
  else cat(sprintf("  %-14s chr n_nonNA=%d vals: %s\n", cl, sum(!is.na(v) & v != ""),
                   paste(head(sort(unique(v)), 12), collapse = ", ")))
}

cat("\n### TCGA-LIHC colData columns mentioning age / gender / sex / stage\n")
se <- readRDS(file.path(DATA_RAW, "TCGA_LIHC_se.rds"))
cd <- as.data.frame(colData(se))
hit <- grep("age|gender|sex|stage|grade|vital|survival|days", colnames(cd), ignore.case = TRUE, value = TRUE)
for (cl in hit) {
  v <- cd[[cl]]
  if (is.numeric(v)) cat(sprintf("  %-26s numeric n_nonNA=%d range=[%g,%g]\n", cl,
                                 sum(!is.na(v)), min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
  else cat(sprintf("  %-26s chr n_nonNA=%d top: %s\n", cl, sum(!is.na(v) & v != ""),
                   paste(head(sort(table(v), decreasing = TRUE), 6), collapse = " ")))
}
cat("\n  all colData cols:", paste(colnames(cd), collapse = " | "), "\n")

cat("\n### GSE141198 / GSE57338 sanity\n")
r5 <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
for (nm in c("GSE57338", "GSE141910")) {
  a <- r5[[nm]]$ann
  cat(sprintf("\n-- %s (n=%d) --\n", nm, nrow(a)))
  for (cl in colnames(a)) {
    v <- a[[cl]]
    nv <- suppressWarnings(as.numeric(as.character(v)))
    cat(sprintf("  %-16s nonNA=%d  numeric-castable=%d  head: %s\n", cl,
                sum(!is.na(v) & v != ""), sum(!is.na(nv)),
                paste(head(v, 5), collapse = ",")))
  }
}
sink()
cat("written\n")
