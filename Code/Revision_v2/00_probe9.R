# 00_probe9.R -- structures of the validation cohorts (light: headers only).
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
V <- VALID_DIR

peek <- function(path, n = 4, width = 160) {
  if (!file.exists(path)) { cat("  MISSING:", path, "\n"); return(invisible()) }
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  l <- readLines(con, n = n, warn = FALSE); close(con)
  cat(paste(substr(l, 1, width), collapse = "\n"), "\n")
}

cat("========== GSE116250 (HF, RNA-seq) ==========\n")
cat("-- rpkm.txt.gz --\n"); peek(file.path(V, "GSE116250_rpkm.txt.gz"), 3, 200)
cat("\n-- series matrix (!Sample lines) --\n")
con <- gzfile(file.path(V, "GSE116250_series_matrix.txt.gz"), "rt")
l <- readLines(con, n = 200, warn = FALSE); close(con)
cat(paste(substr(grep("^!Sample_(title|source_name|characteristics)", l, value = TRUE), 1, 200),
          collapse = "\n"), "\n")
cat("table_begin at:", grep("^!series_matrix_table_begin", l), "\n")

cat("\n========== GSE141910 (HF, RNA-seq) ==========\n")
cat("-- series matrix !Sample lines --\n")
con <- gzfile(file.path(V, "GSE141910_series_matrix.txt.gz"), "rt")
l <- readLines(con, n = 200, warn = FALSE); close(con)
cat(paste(substr(grep("^!Sample_(title|source_name|characteristics)", l, value = TRUE), 1, 220),
          collapse = "\n"), "\n")
cat("-- per-sample raw file example --\n")
f <- list.files(file.path(V, "GSE141910"), full.names = TRUE)[1]
cat("file:", basename(f), "\n")
peek(f, 3, 160)
cat("n raw files:", length(list.files(file.path(V, "GSE141910"))), "\n")

cat("\n========== GSE14520 (HCC, array) ==========\n")
for (g in c("GPL3921", "GPL571")) {
  p <- file.path(DATA_RAW, sprintf("GSE14520-%s_series_matrix.txt.gz", g))
  con <- gzfile(p, "rt"); l <- readLines(con, n = 200, warn = FALSE); close(con)
  tb <- grep("^!series_matrix_table_begin", l)
  cat(sprintf("-- %s: table_begin=%s --\n", g, if (length(tb)) tb else "?"))
  cat(paste(substr(grep("^!Sample_(title|characteristics)", l, value = TRUE)[1:4], 1, 200),
            collapse = "\n"), "\n")
  if (length(tb)) cat("  first data row:", substr(l[tb + 2], 1, 140), "\n")
}

cat("\n========== GSE76427 (HCC, array GPL10558) ==========\n")
con <- gzfile(file.path(DATA_RAW, "GSE76427_series_matrix.txt.gz"), "rt")
l <- readLines(con, n = 200, warn = FALSE); close(con)
tb <- grep("^!series_matrix_table_begin", l)
cat("table_begin:", tb, "\n")
cat(paste(substr(grep("^!Sample_(title|characteristics)", l, value = TRUE)[1:4], 1, 200),
          collapse = "\n"), "\n")
if (length(tb)) cat("first data row:", substr(l[tb + 2], 1, 140), "\n")

cat("\n========== TCGA-LIHC SummarizedExperiment ==========\n")
suppressPackageStartupMessages(library(SummarizedExperiment))
se <- readRDS(file.path(DATA_RAW, "TCGA_LIHC_se.rds"))
cd <- as.data.frame(colData(se))
cat("dim:", paste(dim(se), collapse = " x "), "\n")
cat("sample_type:\n"); print(table(cd$sample_type))
cat("tumor_descriptor:\n"); print(table(cd$tumor_descriptor))
cat("assay 'unstranded' range:", sprintf("%.1f .. %.1f",
    min(assay(se, "unstranded")), max(assay(se, "unstranded"))), "\n")
cat("has rownames(ensembl):", head(rownames(se), 3), "\n")
cat("rowData cols:", paste(colnames(rowData(se)), collapse = ", "), "\n")
