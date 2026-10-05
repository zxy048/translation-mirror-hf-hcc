# 00_probe.R -- reconnaissance of data structures before writing R01.
# Read-only: loads and reports, writes nothing.

suppressPackageStartupMessages({
  library(msigdbr)
})

PROJ <- "D:/R_projects/revision_analysis"
safe_range <- function(o) {
  n <- suppressWarnings(tryCatch(as.numeric(unlist(o)), error = function(e) NA))
  if (all(is.na(n))) return("non-numeric")
  sprintf("%.3f .. %.3f", min(n, na.rm = TRUE), max(n, na.rm = TRUE))
}

cat("=== GSE141198 objects ===\n")
for (f in c("GSE141198_counts_raw.rds","GSE141198_counts_clean.rds","GSE141198_counts_filt.rds",
            "GSE141198_vst.rds","GSE141198_wgcna_input.rds","GSE141198_clinical.rds",
            "GSE141198_pdata.rds")) {
  p <- file.path(PROJ, f)
  if (!file.exists(p)) { cat(sprintf("  %-32s MISSING\n", f)); next }
  o <- readRDS(p)
  if (is.matrix(o) || is.data.frame(o)) {
    cat(sprintf("  %-32s %-11s dim=%d x %d  range=%s\n",
                f, class(o)[1], nrow(o), ncol(o), safe_range(o)))
    cat(sprintf("      rownames[1:3]: %s\n", paste(head(rownames(o),3), collapse=", ")))
    cat(sprintf("      colnames     : %s\n", paste(colnames(o), collapse=", ")))
  } else {
    cat(sprintf("  %-32s %s len=%d\n", f, class(o)[1], length(o)))
    if (is.list(o)) cat(sprintf("      names: %s\n", paste(names(o), collapse=", ")))
  }
}

cat("\n=== GSE57338 series matrix headers ===\n")
sm <- file.path(PROJ, "GSE57338_series_matrix.txt.gz")
con <- gzfile(sm, "rt"); lines <- readLines(con, n = 300); close(con)
for (h in grep("^!Sample_", lines, value = TRUE)) {
  cat("  ", substr(h, 1, 260), "\n")
}

cat("\n=== TCGA-LIHC availability ===\n")
for (p in c("D:/R_projects/TCGA_LIHC_se.rds", "D:/R_projects/GDCdata")) {
  cat(sprintf("  %-40s %s\n", p, file.exists(p)))
}
d <- "D:/R_projects/GDCdata"
if (dir.exists(d)) cat("   GDCdata contents:", paste(head(list.files(d), 10), collapse=", "), "\n")
