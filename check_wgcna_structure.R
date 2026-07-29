# Quick check of WGCNA RDS structure
w <- readRDS("D:/R_projects/revision_analysis/GSE57338_WGCNA_rerun.rds")
cat("=== cross_disease ===\n")
cat("class:", class(w$cross_disease), "\n")
if (is.list(w$cross_disease)) {
  cat("names:", paste(names(w$cross_disease), collapse=", "), "\n")
  for (n in names(w$cross_disease)) {
    x <- w$cross_disease[[n]]
    if (is.vector(x)) {
      cat(sprintf("  %s: length=%d\n", n, length(x)))
    } else if (is.matrix(x) || is.data.frame(x)) {
      cat(sprintf("  %s: dim=(%d,%d)\n", n, nrow(x), ncol(x)))
    } else {
      cat(sprintf("  %s: class=%s\n", n, class(x)))
    }
  }
}

cat("\n=== TOM ===\n")
cat("class:", class(w$TOM), "\n")
if (is.matrix(w$TOM)) {
  cat("dim:", dim(w$TOM), "\n")
}

cat("\n=== gene_symbols ===\n")
cat("length:", length(w$gene_symbols), "\n")
cat("head:", head(w$gene_symbols, 10), "\n")

cat("\n=== moduleColors ===\n")
cat("length:", length(w$moduleColors), "\n")
cat("table:\n")
print(table(w$moduleColors))

cat("\n=== primary_translation_module ===\n")
cat("value:", w$primary_translation_module, "\n")
cat("class:", class(w$primary_translation_module), "\n")
