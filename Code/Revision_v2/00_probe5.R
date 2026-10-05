# 00_probe5.R -- read the v1 Table S2 translation rows for an honest cross-walk.
t <- read.csv("D:/R_projects/revision_analysis/tables/Table_S2_ssGSEA_Effect_Sizes.csv",
              stringsAsFactors = FALSE)
cat("dim:", dim(t), "\n")
cat("cols:", paste(colnames(t), collapse = " | "), "\n\n")
cat("Category counts:\n"); print(table(t$Category))

tr <- t[t$Category == "Translation/Ribosome", ]
cat("\n=== Translation/Ribosome rows (", nrow(tr), ") ===\n", sep = "")
neg <- tr$Pathway[tr$Cohens_d_HCC < 0]
cat("NEGATIVE d_HCC:", length(neg), "\n")
for (i in seq_len(nrow(tr))) {
  cat(sprintf("%2d  dHCC=%+.3f dHF=%+.3f  %-58s | %s\n",
              i, tr$Cohens_d_HCC[i], tr$Cohens_d_HF[i],
              substr(tr$Pathway[i], 1, 58), tr$Short_Label[i]))
}
cat("\nquadrant counts (sign dHCC, sign dHF):\n")
print(table(sign(tr$Cohens_d_HCC), sign(tr$Cohens_d_HF)))
