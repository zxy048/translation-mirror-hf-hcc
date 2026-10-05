source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
r4 <- readRDS(file.path(OUT, "intermediate", "R04_preservation.rds"))
tb <- r4$table
tb$class <- ifelse(tb$Zsummary.pres > 10, "strong",
            ifelse(tb$Zsummary.pres > 2, "weak-mod", "none"))

cat("=== rows per direction, and is grey present? ===\n")
for (ref in unique(tb$reference)) {
  s <- tb[tb$reference == ref, ]
  cat(sprintf("\nreference=%s  n=%d  modules: %s\n", ref, nrow(s),
              paste(sort(s$module), collapse = ", ")))
  cat("  grey present:", any(tolower(s$module) == "grey"), "\n")
}

cat("\n=== class breakdown INCLUDING grey ===\n")
print(table(tb$reference, tb$class))

cat("\n=== class breakdown EXCLUDING grey (named modules only) ===\n")
nb <- tb[tolower(tb$module) != "grey", ]
print(table(nb$reference, nb$class))
cat("\ncounts per direction:\n"); print(table(nb$reference))

cat("\n=== the designated modules ===\n")
print(r4$verdict, row.names = FALSE)
print(r4$designated[, c("reference", "module", "test", "moduleSize",
                        "Zsummary.pres", "medianRank.pres")], row.names = FALSE)

cat("\n=== ranked position of designated module among NAMED modules ===\n")
for (ref in unique(nb$reference)) {
  s <- nb[nb$reference == ref, ]
  s <- s[order(-s$Zsummary.pres), ]
  des <- if (ref == "HF") DES_HF else DES_HCC
  i <- which(tolower(s$module) == des)
  cat(sprintf("  reference=%s: %s is rank %d of %d named modules (Z=%.2f)\n",
              ref, des, i, nrow(s), s$Zsummary.pres[i]))
}
cat("\nPROBE DONE\n")
