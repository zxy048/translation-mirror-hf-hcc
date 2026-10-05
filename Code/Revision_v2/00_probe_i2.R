t <- read.csv("D:/R_projects/revision_analysis/v2_output/tables/Table_S21_cohort_heterogeneity.csv",
              stringsAsFactors = FALSE)
a <- t[t$tier == "A", ]
for (d in c("HF", "HCC")) {
  s <- a[a$disease == d, ]
  cat(sprintf("%s tier A: n=%d  median I2=%.2f  median pQ=%.3g  frac pQ<0.05 = %d/%d = %.0f%%\n",
              d, nrow(s), median(s$I2_pct), median(s$p_Q),
              sum(s$p_Q < 0.05), nrow(s), 100 * mean(s$p_Q < 0.05)))
}
cat(sprintf("pooled over both diseases, tier A: frac pQ<0.05 = %d/%d = %.0f%%\n",
            sum(a$p_Q < 0.05), nrow(a), 100 * mean(a$p_Q < 0.05)))
cat(sprintf("all tiers, both diseases: median I2 HF=%.2f HCC=%.2f\n",
            median(t$I2_pct[t$disease == "HF"]), median(t$I2_pct[t$disease == "HCC"])))
