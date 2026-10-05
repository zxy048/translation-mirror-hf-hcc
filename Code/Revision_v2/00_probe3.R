# 00_probe3.R -- discover the msigdbr 26.1.0 API, then check canonical names.
suppressPackageStartupMessages({ library(msigdbr) })

cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n\n")

cat("=== msigdbr_collections() ===\n")
cc <- tryCatch(msigdbr_collections(), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
if (!is.null(cc)) { print(as.data.frame(cc)) }

cat("\n=== formals(msigdbr) ===\n")
print(names(formals(msigdbr)))
cat("\n=== formals(msigdbr_collections) ===\n")
print(names(formals(msigdbr_collections)))

cat("\n=== try a call ===\n")
x <- tryCatch(msigdbr(species = "Homo sapiens"),
              error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
if (!is.null(x)) {
  cat("dim:", paste(dim(x), collapse=" x "), "\n")
  cat("cols:", paste(colnames(x), collapse=", "), "\n")
  cat("gs_collection values:\n"); print(table(x$gs_collection))
  cat("gs_subcollection head:\n"); print(head(unique(x$gs_subcollection), 30))
}
