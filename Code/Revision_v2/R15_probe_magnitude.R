# Is GSE89377's HCC arm compressed relative to the other HCC cohorts? If its
# effect sizes are systematically smaller, "cirrhosis > HCC" inside GSE89377
# would be a cohort-scaling artefact rather than a biological statement.
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
es7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
r15 <- readRDS(file.path(OUT, "intermediate", "R15_cirrhosis_control.rds"))

ta <- intersect(SETS_TIER_A, unique(es7$set))
hcc7 <- es7[es7$disease == "HCC" & es7$set %in% ta, c("cohort", "set", "hedges_g", "n_disease", "n_control")]
w7 <- reshape(hcc7, idvar = "set", timevar = "cohort", direction = "wide")
names(w7) <- sub("^hedges_g\\.", "", names(w7))

g15 <- r15$wide[r15$wide$tier == "A", ]
m <- merge(g15[, c("set", "Cirrhosis", "HCC")], w7, by = "set")
names(m)[names(m) == "HCC"] <- "GSE89377_HCC"

cohorts <- setdiff(names(m), c("set", "Cirrhosis", "GSE89377_HCC"))
cat("\n=== |g| for the HCC arm, GSE89377 vs each HCC cohort (Tier A) ===\n")
for (cc in cohorts) {
  v <- m[[cc]]
  ok <- is.finite(v) & is.finite(m$GSE89377_HCC)
  cat(sprintf("  %-18s n=%2d | mean|g| %.2f vs GSE89377 %.2f | median|g| %.2f vs %.2f\n",
              cc, sum(ok),
              mean(abs(v[ok])), mean(abs(m$GSE89377_HCC[ok])),
              median(abs(v[ok])), median(abs(m$GSE89377_HCC[ok]))))
}
cat(sprintf("\n  mean|g| across the four external HCC cohorts vs GSE89377 HCC arm: %.2f vs %.2f\n",
            mean(abs(as.matrix(m[, cohorts])), na.rm = TRUE),
            mean(abs(m$GSE89377_HCC), na.rm = TRUE)))

cat("\n=== side by side: cirrhosis, HCC within GSE89377, and HCC elsewhere ===\n")
show <- m[order(-abs(m$GSE89377_HCC)), ]
show$set <- substr(show$set, 1, 46)
print(show, row.names = FALSE, digits = 2)

# If GSE89377 compresses, its HCC arm should be smaller than its OWN cirrhosis
# arm by more than the external cohorts differ from each other.
cat("\n=== within-GSE89377 comparison, and the spread among external cohorts ===\n")
cat(sprintf("  GSE89377: mean|g| cirrhosis %.2f vs HCC %.2f\n",
            mean(abs(m$Cirrhosis)), mean(abs(m$GSE89377_HCC))))
ext <- as.matrix(m[, cohorts])
cat(sprintf("  external HCC cohorts: mean|g| range across cohorts %.2f - %.2f\n",
            min(colMeans(abs(ext), na.rm = TRUE)), max(colMeans(abs(ext), na.rm = TRUE))))
cat(sprintf("  paired: cirrhosis > HCC in %d / %d Tier A sets\n",
            sum(abs(m$Cirrhosis) > abs(m$GSE89377_HCC)),
            nrow(m)))
cat("PROBE DONE\n")
