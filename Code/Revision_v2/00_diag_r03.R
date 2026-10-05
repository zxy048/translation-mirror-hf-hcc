# 00_diag_r03.R -- what ARE the two designated translation modules?
# The cross-disease gene overlap came out at 5/294 (Jaccard 0.017). Before
# reporting that, establish whether it is biology (different programs) or an
# artifact of the designation rule (each module being a different KIND of set).
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

r2 <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))
r3 <- readRDS(file.path(OUT, "intermediate", "R03_translation_module.rds"))
A  <- r3$tierA_genes

labHF  <- r2$HF$net$colors
labHCC <- r2$HCC$net$colors

mm <- list(HF = r3$trans_genes$HF, HCC = r3$trans_genes$HCC)

cat("=== 1. overlap composition ===\n")
ov <- intersect(mm$HF, mm$HCC)
cat("shared genes (", length(ov), "):", paste(sort(ov), collapse = ", "), "\n")
cat("shared that are Tier A:", paste(sort(intersect(ov, A)), collapse = ", "), "\n")

cat("\n=== 2. Tier A genes inside each designated module ===\n")
aHF <- intersect(mm$HF, A); aHCC <- intersect(mm$HCC, A)
cat(sprintf("HF  module: %d/%d Tier A (%.0f%%)\n", length(aHF), length(mm$HF),
            100 * length(aHF) / length(mm$HF)))
cat(sprintf("HCC module: %d/%d Tier A (%.0f%%)\n", length(aHCC), length(mm$HCC),
            100 * length(aHCC) / length(mm$HCC)))
cat(sprintf("Tier A shared between the two modules: %d  [%s]\n",
            length(intersect(aHF, aHCC)), paste(sort(intersect(aHF, aHCC)), collapse = ", ")))
cat(sprintf("Tier A total in U: %d ; in HF module: %d ; in HCC module: %d\n",
            length(A), length(aHF), length(aHCC)))

cat("\n=== 3. where do the Tier A genes live in each network? ===\n")
for (nm in c("HF", "HCC")) {
  lab <- if (nm == "HF") labHF else labHCC
  ta  <- intersect(names(lab), A)
  t <- sort(table(lab[ta]), decreasing = TRUE)
  cat(sprintf("\n-- %s: %d/%d Tier A genes present in the network --\n", nm, length(ta), length(A)))
  print(head(data.frame(module = names(t), n_tierA = as.integer(t),
                        pct_of_module = round(100 * as.integer(t) / as.integer(table(lab)[names(t)])), 1),
             10), row.names = FALSE)
}

cat("\n=== 4. top genes by |kME| in each designated module ===\n")
for (nm in c("HF", "HCC")) {
  k <- if (nm == "HF") r2$kmeHF else r2$kmeHCC
  m <- if (nm == "HF") r3$desHF$module else r3$desHCC$module
  d <- k[k$module == m, ]
  d <- d[order(-abs(d$kME)), ]
  cat(sprintf("\n-- %s module '%s' (n=%d); top 25 by |kME| --\n", nm, m, nrow(d)))
  cat(paste(head(d$gene, 25), collapse = " "), "\n")
  cat("   marked * = Tier A\n")
  cat(paste(ifelse(head(d$gene, 40) %in% A, paste0(head(d$gene, 40), "*"), head(d$gene, 40)),
            collapse = " "), "\n")
}

cat("\n=== 5. do the two modules connect? cross-tab of Tier A membership ===\n")
tab <- table(TierA_HF = ifelse(ov <- NULL, "", ""))
allg <- union(mm$HF, mm$HCC)
df <- data.frame(
  gene = allg,
  in_HF  = allg %in% mm$HF,
  in_HCC = allg %in% mm$HCC,
  is_TierA = allg %in% A, stringsAsFactors = FALSE)
print(table(df$in_HF, df$in_HCC, dnn = c("in_HF_module", "in_HCC_module")))
cat("\namong Tier A genes:\n")
print(table(df$in_HF[df$is_TierA], df$in_HCC[df$is_TierA],
            dnn = c("in_HF_module", "in_HCC_module")))
