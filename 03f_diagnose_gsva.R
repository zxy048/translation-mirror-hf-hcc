# 诊断GSVA 2.6.2 API
library(GSVA)
cat("GSVA version:", as.character(packageVersion("GSVA")), "\n\n")

# Check available functions
cat("=== gsva function args ===\n")
print(args(gsva))
cat("\n")

# Check if parametric objects exist
cat("gsvaParam exists:", exists("gsvaParam"), "\n")
cat("ssgseaParam exists:", exists("ssgseaParam"), "\n")
cat("gsva is S4 generic:", isS4(gsva), "\n")

# Test new API with a tiny example
if (exists("gsvaParam")) {
  cat("\n=== gsvaParam args ===\n")
  print(args(gsvaParam))

  # Test with tiny data
  set.seed(42)
  tiny_expr <- matrix(rnorm(1000), nrow=100, ncol=10,
                       dimnames=list(paste0("gene", 1:100), paste0("samp", 1:10)))
  tiny_gs <- list(set1=paste0("gene", 1:20), set2=paste0("gene", 21:40))
  tiny_gs2 <- list(set1=paste0("gene", 1:20))

  cat("\n=== Test gsvaParam ===\n")
  param_test <- tryCatch({
    gsvaParam(exprData = tiny_expr, geneSets = tiny_gs, minSize = 5, maxSize = 500)
    cat("Success!\n")
  }, error = function(e) cat("Error:", e$message, "\n"))

  # Try only with kcdf
  cat("\n=== Test gsvaParam with kcdf ===\n")
  param_test2 <- tryCatch({
    gsvaParam(exprData = tiny_expr, geneSets = tiny_gs, minSize = 5, maxSize = 500, kcdf = "Gaussian")
    cat("Success!\n")
  }, error = function(e) cat("Error:", e$message, "\n"))

  # What class does gsva expect?
  cat("\n=== showMethods(gsva) ===\n")
  tryCatch({
    print(showMethods("gsva"))
  }, error = function(e) cat("Error:", e$message, "\n"))
}
