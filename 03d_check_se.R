library(SummarizedExperiment)
se <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")
cat("class:", class(se), "\n")
cat("assayNames:", paste(assayNames(se), collapse=", "), "\n")
cat("dim:", paste(dim(se), collapse=" x "), "\n")
# First few rownames and colnames
cat("First 3 genes:", paste(head(rownames(se), 3), collapse=", "), "\n")
cat("First 3 samples:", paste(head(colnames(se), 3), collapse=", "), "\n")
# Check expression range
a1 <- assay(se, 1)
cat("Expr range:", range(a1, na.rm=TRUE), "\n")
cat("Expr[1:3,1:3]:\n")
print(a1[1:3, 1:3])
