pk <- c("org.Hs.eg.db","AnnotationDbi","GSVA","limma","SummarizedExperiment","hugene11sttranscriptcluster.db","data.table","matrixStats","cluster","fields")
for (p in pk) cat(sprintf("%-38s %s\n", p, if (requireNamespace(p, quietly=TRUE)) "OK" else "MISSING"))
cat("\n-- generic soft parser in R01? --\n")
