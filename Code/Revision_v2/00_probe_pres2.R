source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
r2 <- readRDS(file.path(OUT,"intermediate","R02_wgcna.rds"))
disableWGCNAThreads()
multiExpr <- list(HF = list(data=r2$datHF), HCC = list(data=r2$datHCC))
colorList <- list(HF = r2$HF$net$colors, HCC = r2$HCC$net$colors)
pres <- WGCNA::modulePreservation(
  multiExpr, colorList, referenceNetworks=c(1,2), testNetworks=list(2,1),
  nPermutations=2, randomSeed=STATS$SEED, networkType=WGCNA$networkType,
  verbose=0, indent=0)
show <- function(lbl, x) {
  cat("\n### ", lbl, " class=", paste(class(x), collapse="/"), "\n", sep="")
  if (is.data.frame(x) || is.matrix(x)) {
    cat("dim:", paste(dim(x), collapse="x"), "\n")
    cat("cols:", paste(colnames(x), collapse=" | "), "\n")
    cat("rows:", paste(rownames(x), collapse=","), "\n")
    print(utils::head(x, 3))
  } else { print(utils::str(x, max.level=2)) }
}
show("Z$ref.HF$inColumnsAlsoPresentIn.HCC", pres$preservation$Z[["ref.HF"]][["inColumnsAlsoPresentIn.HCC"]])
show("Z$ref.HCC$inColumnsAlsoPresentIn.HF", pres$preservation$Z[["ref.HCC"]][["inColumnsAlsoPresentIn.HF"]])
show("observed$ref.HF$...HCC", pres$preservation$observed[["ref.HF"]][["inColumnsAlsoPresentIn.HCC"]])
show("log.pBonf$ref.HF$...HCC", pres$preservation$log.pBonf[["ref.HF"]][["inColumnsAlsoPresentIn.HCC"]])
show("q$ref.HF$...HCC", pres$preservation$q[["ref.HF"]][["inColumnsAlsoPresentIn.HCC"]])
