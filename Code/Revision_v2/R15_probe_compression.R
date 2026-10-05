# Is GSE89377's smaller |g| a cohort-wide (platform / pipeline) scaling effect,
# or specific to the translation sets? If it is cohort-wide, the internal
# cirrhosis-vs-HCC comparison is still interpretable; if it is specific to the
# translation sets, that comparison is confounded and must not be asserted.
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
es7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
r15 <- readRDS(file.path(OUT, "intermediate", "R15_cirrhosis_control.rds"))

g15 <- r15$wide
hcc7 <- es7[es7$disease == "HCC", c("cohort", "set", "hedges_g")]
w7 <- reshape(hcc7, idvar = "set", timevar = "cohort", direction = "wide")
names(w7) <- sub("^hedges_g\\.", "", names(w7))
names(w7)[names(w7) == "HCC" & FALSE] <- NA  # no-op; cohorts are named explicitly

m <- merge(g15[, c("set", "tier", "Cirrhosis", "HCC")], w7, by = "set")
names(m)[names(m) == "HCC"] <- "GSE89377_HCC"
ext <- setdiff(names(w7), "set")

cat("=== ALL shared sets (any tier), n =", nrow(m), "===\n")
cat(sprintf("  GSE89377 HCC arm mean|g| = %.3f\n", mean(abs(m$GSE89377_HCC), na.rm = TRUE)))
for (cc in ext)
  cat(sprintf("  %-18s mean|g| = %.3f | ratio to GSE89377 = %.2fx\n",
              cc, mean(abs(m[[cc]]), na.rm = TRUE),
              mean(abs(m[[cc]]), na.rm = TRUE) / mean(abs(m$GSE89377_HCC), na.rm = TRUE)))

cat("\n=== by tier: is the deficit uniform across tiers? ===\n")
for (tr in unique(m$tier)) {
  s <- m[m$tier == tr, ]
  if (!nrow(s)) next
  cat(sprintf("  tier %-8s n=%2d | GSE89377 %.2f | TCGA %.2f | GSE76427 %.2f | GPL3921 %.2f\n",
              tr, nrow(s),
              mean(abs(s$GSE89377_HCC), na.rm = TRUE),
              mean(abs(s$TCGA_LIHC), na.rm = TRUE),
              mean(abs(s$GSE76427), na.rm = TRUE),
              mean(abs(s$GSE14520_GPL3921), na.rm = TRUE)))
}

# The decisive comparison: the ratio (GSE89377 HCC / TCGA HCC) computed over
# Tier A sets versus over non-Tier-A sets. If the two ratios agree, the deficit
# is cohort-wide scaling; the internal cirrhosis comparison then stands.
rA <- mean(abs(m$GSE89377_HCC[m$tier == "A"]), na.rm = TRUE) /
      mean(abs(m$TCGA_LIHC[m$tier == "A"]), na.rm = TRUE)
rO <- mean(abs(m$GSE89377_HCC[m$tier != "A"]), na.rm = TRUE) /
      mean(abs(m$TCGA_LIHC[m$tier != "A"]), na.rm = TRUE)
cat(sprintf("\n  GSE89377/TCGA |g| ratio: Tier A %.2f vs non-Tier-A %.2f  -> %s\n",
            rA, rO,
            if (abs(rA - rO) < 0.25) "UNIFORM (cohort-wide scaling)"
            else "NON-UNIFORM (translation sets specifically attenuated)"))

cat("\n=== TG subgroups: is the pooled HCC arm attenuated by heterogeneity? ===\n")
cat("  (need per-TG scores; see note)\n")
sc <- r15$scores; stg <- r15$stage
for (s in c("KEGG_RIBOSOME", "REACTOME_TRANSLATION", "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")) {
  if (!s %in% rownames(sc)) next
  ctl <- sc[s, names(stg)[stg == "Normal"]]
  parts <- vapply(c("HCC", "eHCC", "Cirrhosis", "DN_high"), function(a) {
    x <- sc[s, names(stg)[stg == a]]
    n1 <- length(x); n2 <- length(ctl)
    sp <- sqrt(((n1-1)*sd(x)^2 + (n2-1)*sd(ctl)^2) / (n1+n2-2))
    round((mean(x) - mean(ctl)) / sp, 3)
  }, numeric(1))
  cat(sprintf("  %-42s HCC %.2f | eHCC %.2f | Cirrh %.2f | DN_high %.2f\n",
              substr(s, 1, 42), parts[1], parts[2], parts[3], parts[4]))
}
cat("PROBE DONE\n")
