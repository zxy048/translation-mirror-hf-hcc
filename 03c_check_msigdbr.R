library(msigdbr)

# Try different approaches to find collection names
cat("=== Try 1: msigdbr_collections() ===\n")
tryCatch({
  cols <- msigdbr_collections()
  print(head(cols, 30))
}, error = function(e) cat("Error:", e$message, "\n"))

cat("\n=== Try 2: Search for Hallmark genes directly ===\n")
# Just try getting hallmark genes without specifying collection
tryCatch({
  h_df <- msigdbr(species = "Homo sapiens")
  cats <- unique(h_df$gs_cat)
  cat("Available gs_cat values:\n")
  print(cats)
  cat("\nAvailable gs_subcat values:\n")
  print(unique(h_df$gs_subcat))
}, error = function(e) cat("Error:", e$message, "\n"))
