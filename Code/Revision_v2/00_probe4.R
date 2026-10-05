# 00_probe4.R -- list the REAL translation-related gs_name values in each
# relevant subcollection of msigdbr 26.1.0.
suppressPackageStartupMessages({ library(msigdbr) })

pull <- function(coll, sub = NULL) {
  a <- list(db_species = "HS", species = "Homo sapiens", collection = coll)
  if (!is.null(sub)) a$subcollection <- sub
  d <- do.call(msigdbr, a)
  data.frame(sub = if (is.null(sub)) coll else sub,
             gs_name = unique(d$gs_name),
             n_gene  = as.integer(table(d$gs_name)[unique(d$gs_name)]),
             stringsAsFactors = FALSE)
}

sets <- rbind(
  pull("H"),
  pull("C2", "CP:REACTOME"),
  pull("C2", "CP:KEGG_LEGACY"),
  pull("C2", "CP:KEGG_MEDICUS")
)
cat("total sets pulled:", nrow(sets), "\n\n")

pat <- "TRANSLAT|RIBOSOM|RRNA|AMINOACYL|TRNA|INITIATION_FACTOR|ELONGATION_FACTOR"
tr <- sets[grepl(pat, sets$gs_name, ignore.case = TRUE), ]
tr <- tr[order(tr$sub, tr$gs_name), ]
cat("=== translation-related sets (", nrow(tr), ") ===\n", sep = "")
for (i in seq_len(nrow(tr))) {
  cat(sprintf("  %-16s n=%-5d %s\n", tr$sub[i], tr$n_gene[i], tr$gs_name[i]))
}

cat("\n=== Hallmark: MYC / mTORC1 / UPR / E2F / G2M (the exclusion list) ===\n")
h <- sets[sets$sub == "H", ]
for (nm in c("HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2","HALLMARK_MTORC1_SIGNALING",
             "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")) {
  r <- h[h$gs_name == nm, ]
  cat(sprintf("  %-40s %s\n", nm, if (nrow(r)) paste0("OK n=", r$n_gene) else "MISSING"))
}

cat("\n=== does KEGG_RIBOSOME exist, and under which subcollection? ===\n")
for (nm in c("KEGG_RIBOSOME","KEGG_AMINOACYL_TRNA_BIOSYNTHESIS","KEGG_RNA_POLYMERASE")) {
  r <- sets[sets$gs_name == nm, ]
  cat(sprintf("  %-36s %s\n", nm,
      if (nrow(r)) paste(sprintf("%s(n=%d)", r$sub, r$n_gene), collapse="; ") else "MISSING"))
}
