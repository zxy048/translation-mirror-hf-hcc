library(SummarizedExperiment)
library(DESeq2)

se <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")
counts_lihc <- assay(se, "unstranded")

cat("class:", class(counts_lihc), "\n")
cat("type:", typeof(counts_lihc), "\n")
cat("dim:", dim(counts_lihc), "\n")

# Check rownames
rn <- rownames(counts_lihc)
cat("rownames[1:5]:", paste(head(rn, 5), collapse=", "), "\n")
cat("duplicated rownames:", sum(duplicated(rn)), "\n")
cat("NA rownames:", sum(is.na(rn)), "\n")

# Check for sparse matrix / DelayedMatrix
if (is(counts_lihc, "DelayedMatrix")) {
  cat("Is DelayedMatrix - converting to dense\n")
  counts_dense <- as.matrix(counts_lihc)
  cat("dense class:", class(counts_dense), "\n")
} else {
  counts_dense <- as.matrix(counts_lihc)
}

# Try DESeq2 with dense matrix
dds <- DESeqDataSetFromMatrix(
  countData = counts_dense,
  colData = data.frame(row.names = colnames(counts_dense), cond = factor(1)),
  design = ~ 1)
cat("DESeqDataSet created successfully\n")
cat("dds dim:", dim(dds), "\n")
