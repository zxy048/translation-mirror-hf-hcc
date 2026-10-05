source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
r5 <- readRDS(file.path(OUT,"intermediate","R05_cohorts.rds"))
for (nm in c("GSE141910","GSE57338","GSE116250")) {
  co <- r5[[nm]]
  cat("\n==== ", nm, " (", ncol(co$expr), " samples) ====\n", sep="")
  cat("ann cols:", paste(colnames(co$ann), collapse=", "), "\n")
  et <- co$ann$etiology
  if (!is.null(et)) {
    tb <- sort(table(et, useNA="ifany"), decreasing=TRUE)
    print(tb)
    g <- co$group[co$ann$gsm]
    cat("\n-- group x etiology --\n")
    print(table(group=g, etiology=et, useNA="ifany"))
  } else cat("(no etiology column)\n")
}
