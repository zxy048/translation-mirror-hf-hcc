# 00_probe_covariates3.R (diagnostic) -- TCGA-LIHC colData inventory
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(SummarizedExperiment))
se <- readRDS(file.path(DATA_RAW, "TCGA_LIHC_se.rds"))
cd <- as.data.frame(colData(se))
out <- paste("ALL COLS:", paste(colnames(cd), collapse = " | "))
for (cl in colnames(cd)) {
  v <- cd[[cl]]
  if (is.numeric(v)) {
    out <- c(out, sprintf("  %-28s numeric n=%d range=[%g,%g]", cl,
                          sum(!is.na(v)), min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
  } else if (length(unique(v)) <= 12) {
    out <- c(out, sprintf("  %-28s chr n=%d: %s", cl, sum(!is.na(v) & v != ""),
                          paste(names(sort(table(v), decreasing = TRUE)), collapse = " ")))
  } else {
    out <- c(out, sprintf("  %-28s chr n=%d (%d levels)", cl,
                          sum(!is.na(v) & v != ""), length(unique(v))))
  }
}
writeLines(out, file.path(OUT, "logs", "00_probe_covariates3.txt"))
cat("done\n")
