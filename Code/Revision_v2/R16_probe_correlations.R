# Diagnostic only: cross-disease Spearman rho for the v2 manuscript's comparison
# pairs, so the Results text can quote exact numbers that match R07/R09.
# Produces no shipping artefact.
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
es <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
w <- reshape(es[, c("cohort", "set", "tier", "hedges_g")],
             idvar = c("set", "tier"), timevar = "cohort", direction = "wide")
names(w) <- sub("^hedges_g[.]", "", names(w))

cat("cols:", paste(names(w), collapse = ", "), "\n\n")
pairs <- list(c("TCGA_LIHC", "GSE57338"), c("TCGA_LIHC", "GSE116250"),
              c("TCGA_LIHC", "GSE141910"), c("GSE14520_GPL3921", "GSE57338"))
for (p in pairs) {
  for (tr in list("A", c("A", "B_extra", "C_v1"))) {
    if (!all(p %in% names(w))) next
    d <- w[w[["tier"]] %in% tr, ]
    d <- d[is.finite(d[[p[1]]]) & is.finite(d[[p[2]]]), ]
    cat(sprintf("%-18s vs %-10s %-5s n=%2d rho=%+.3f\n", p[1], p[2],
                if (length(tr) == 1) "TierA" else "All", nrow(d),
                cor(d[[p[1]]], d[[p[2]]], method = "spearman")))
  }
}

cat("\n--- every HCC cohort vs GSE57338, Tier A ---\n")
for (cc in c("TCGA_LIHC", "GSE141198", "GSE14520_GPL3921", "GSE14520_GPL571", "GSE76427")) {
  if (!cc %in% names(w)) next
  d <- w[w[["tier"]] == "A", ]
  d <- d[is.finite(d[[cc]]) & is.finite(d[["GSE57338"]]), ]
  cat(sprintf("  %-18s n=%2d rho=%+.3f\n", cc, nrow(d),
              cor(d[[cc]], d[["GSE57338"]], method = "spearman")))
}

r15 <- readRDS(file.path(OUT, "intermediate", "R15_cirrhosis_control.rds"))
g15 <- r15[["wide"]]
cat("\n--- GSE89377 (R15) ---\n")
for (tr in list("A", c("A", "B_extra", "C_v1"))) {
  a <- g15[g15[["tier"]] %in% tr, c("set", "Cirrhosis", "HCC")]
  b <- w[w[["tier"]] %in% tr, c("set", "TCGA_LIHC", "GSE57338")]
  m <- merge(a, b, by = "set")
  cat(sprintf("  %-6s n=%2d | cirrhosis vs TCGA rho=%+.3f | GSE89377-HCC vs TCGA rho=%+.3f | GSE89377-HCC vs GSE57338 rho=%+.3f\n",
              if (length(tr) == 1) "TierA" else "All", nrow(m),
              cor(m[["Cirrhosis"]], m[["TCGA_LIHC"]], method = "spearman"),
              cor(m[["HCC"]], m[["TCGA_LIHC"]], method = "spearman"),
              cor(m[["HCC"]], m[["GSE57338"]], method = "spearman")))
}
cat("PROBE DONE\n")
