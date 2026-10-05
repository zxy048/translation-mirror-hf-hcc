source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
r2 <- readRDS(file.path(OUT,"intermediate","R02_wgcna.rds"))
disableWGCNAThreads()
multiExpr <- list(HF = list(data=r2$datHF), HCC = list(data=r2$datHCC))
colorList <- list(HF = r2$HF$net$colors, HCC = r2$HCC$net$colors)
pres <- WGCNA::modulePreservation(
  multiExpr, colorList, referenceNetworks=c(1,2), testNetworks=list(2,1),
  nPermutations=2, randomSeed=STATS$SEED, networkType=WGCNA$networkType,
  verbose=0, indent=0)
cat("\n=== top-level names ===\n"); print(names(pres))
cat("\n=== preservation names ===\n"); print(names(pres$preservation))
cat("\n=== class of Z ===\n"); print(class(pres$preservation$Z))
cat("length Z:", length(pres$preservation$Z), "\n")
cat("names Z:", paste(names(pres$preservation$Z), collapse=" | "), "\n")
z1 <- pres$preservation$Z[[1]]
cat("\nclass Z[[1]]:", class(z1), " length:", length(z1), "\n")
cat("names Z[[1]]:", paste(names(z1), collapse=" | "), "\n")
cat("\nclass Z[[1]][[1]]:", class(pres$preservation$Z[[1]][[1]]), "\n")
cat("dim:", paste(dim(pres$preservation$Z[[1]][[1]]), collapse="x"), "\n")
cat("colnames:", paste(colnames(pres$preservation$Z[[1]][[1]]), collapse=" | "), "\n")
cat("\n--- head ---\n"); print(head(pres$preservation$Z[[1]][[1]], 4))
cat("\n=== observed names ===\n"); print(names(pres$preservation$Observed))
cat("class Observed[[1]][[1]]:", class(pres$preservation$Observed[[1]][[1]]), "\n")
cat("colnames Observed:", paste(colnames(pres$preservation$Observed[[1]][[1]]), collapse=" | "), "\n")
