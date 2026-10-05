source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
r2 <- readRDS(file.path(OUT,"intermediate","R02_wgcna.rds")); r3 <- readRDS(file.path(OUT,"intermediate","R03_translation_module.rds"))
A <- r3$tierA_genes; U <- r2$U
mt <- grep("^MT-", A, value = TRUE)
cat("MT- genes in Tier A canonical list:", length(mt), "\n")
cat("  present in U :", length(intersect(mt,U)), " [", paste(sort(intersect(mt,U)),collapse=", "), "]\n")
cat("  LOST from U  :", length(setdiff(mt,U)), " [", paste(sort(setdiff(mt,U)),collapse=", "), "]\n")
# nuclear-encoded mitoribosome (MRPL/MRPS) coverage
mr <- grep("^(MRPL|MRPS)", A, value = TRUE)
cat("\nMRPL/MRPS in Tier A:", length(mr), "| in U:", length(intersect(mr,U)),
    "| lost:", length(setdiff(mr,U)), "\n")
cat("  lost:", paste(sort(setdiff(mr,U)),collapse=", "), "\n")
