# 诊断GSE57338系列矩阵文件内容
library(GEOquery)
library(Biobase)

# 下载系列矩阵（不下载GPL）
gse <- getGEO(filename = "D:/R_projects/GSE57338_series_matrix.txt.gz", getGPL = FALSE)

cat("=== Expression matrix dim ===\n")
cat(dim(exprs(gse)), "\n")
cat("First 5 probes:", paste(head(rownames(exprs(gse)), 5), collapse=", "), "\n")

cat("\n=== fData columns ===\n")
fdata <- fData(gse)
cat(ncol(fdata), "columns:\n")
for (cn in colnames(fdata)) {
  vals <- head(fdata[[cn]], 3)
  cat(sprintf("  %s: %s\n", cn, paste(as.character(vals), collapse=" | ")))
}

cat("\n=== Phenotype data columns ===\n")
pd <- pData(gse)
for (cn in colnames(pd)) {
  vals <- unique(head(pd[[cn]], 5))
  cat(sprintf("  %s: %s\n", cn, paste(as.character(vals), collapse=" | ")))
}

cat("\n=== Phenotype characteristics ===\n")
for (cn in grep("characteristics", colnames(pd), value=TRUE)) {
  vals <- unique(pd[[cn]])
  cat(sprintf("  %s:\n", cn))
  for (v in head(vals, 10)) cat(sprintf("    %s\n", as.character(v)))
}

cat("\n=== Platform ID ===\n")
cat("GPL:", annotation(gse), "\n")
cat("Platform:", gse@annotation, "\n")
