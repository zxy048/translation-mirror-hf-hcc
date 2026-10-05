# 00_probe2.R -- verify TCGA availability and that every canonical gene-set name
# we intend to freeze actually exists in this msigdbr version (26.1.0).

suppressPackageStartupMessages({ library(msigdbr) })

cat("=== TCGA-LIHC on disk? ===\n")
for (p in c("D:/R_projects/TCGA_LIHC_se.rds",
            "D:/R_projects/revision_analysis/TCGA_LIHC_se.rds",
            "D:/R_projects/GDCdata")) {
  cat(sprintf("  %-52s %s\n", p, file.exists(p)))
}
if (file.exists("D:/R_projects/TCGA_LIHC_se.rds")) {
  o <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")
  cat("  class:", class(o)[1], "\n")
  if (is.list(o)) cat("  names:", paste(names(o), collapse=", "), "\n")
  if (requireNamespace("SummarizedExperiment", quietly=TRUE) &&
      methods::is(o, "SummarizedExperiment")) {
    cat("  dim:", paste(dim(o), collapse=" x "), "\n")
    cat("  assays:", paste(SummarizedExperiment::assayNames(o), collapse=", "), "\n")
    cat("  colData cols:", paste(head(colnames(SummarizedExperiment::colData(o)), 20), collapse=", "), "\n")
  }
}

cat("\n=== msigdbr 26.1.0 -- which canonical names exist? ===\n")
CANON <- list(
  KEGG = c("KEGG_RIBOSOME", "KEGG_AMINOACYL_TRNA_BIOSYNTHESIS", "KEGG_RNA_POLYMERASE"),
  REACTOME = c(
    "REACTOME_TRANSLATION",
    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
    "REACTOME_EUKARYOTIC_TRANSLATION_TERMINATION",
    "REACTOME_CAP_DEPENDENT_TRANSLATION_INITIATION",
    "REACTOME_FORMATION_OF_THE_TERNARY_COMPLEX_AND_SUBSEQUENTLY_THE_43S_COMPLEX",
    "REACTOME_TRNA_AMINOACYLATION",
    "REACTOME_CYTOSOLIC_TRNA_AMINOACYLATION",
    "REACTOME_RRNA_PROCESSING",
    "REACTOME_RRNA_PROCESSING_IN_THE_NUCLEUS_AND_CYTOSOL",
    "REACTOME_RIBOSOME_BIOGENESIS",
    "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
    "REACTOME_3_UTR_MEDIATED_TRANSLATIONAL_REGULATION",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
    "REACTOME_MITOCHONDRIAL_TRANSLATION",
    "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S",
    "REACTOME_GTP_HYDROLYSIS_AND_JOINING_OF_THE_60S_RIBOSOMAL_SUBUNIT",
    "REACTOME_L13A_MEDIATED_TRANSLATIONAL_TERMINATION",
    "REACTOME_REGULATION_OF_EXPRESSION_OF_SLIT_AND_ROBO_GENES"
  ),
  HALLMARK_EXCLUDE = c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
                       "HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                       "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
)

all_sets <- tryCatch({
  a <- msigdbr(species = "Homo sapiens", collection = "CP")
  b <- msigdbr(species = "Homo sapiens", collection = "H")
  rbind(
    data.frame(gs_name = unique(a$gs_name), coll = "CP"),
    data.frame(gs_name = unique(b$gs_name), coll = "H")
  )
}, error = function(e) { cat("  msigdbr ERROR:", conditionMessage(e), "\n"); NULL })

if (!is.null(all_sets)) {
  cat("  total distinct set names:", nrow(all_sets), "\n")
  for (grp in names(CANON)) {
    cat(sprintf("\n  --- %s ---\n", grp))
    for (nm in CANON[[grp]]) {
      hit <- all_sets$gs_name[all_sets$gs_name == nm]
      if (length(hit)) {
        cat(sprintf("    OK      %s\n", nm))
      } else {
        # fuzzy: find close candidates
        cand <- grep(sub("^(KEGG|REACTOME|HALLMARK)_", "", nm), all_sets$gs_name,
                     value = TRUE, ignore.case = TRUE)
        cand <- setdiff(cand, nm)
        cat(sprintf("    MISSING %s   -> candidates: %s\n", nm,
                    paste(head(cand, 4), collapse=" | ")))
      }
    }
  }
}

cat("\n=== all Reactome translation-ish names available ===\n")
if (!is.null(all_sets)) {
  tr <- grep("TRANSLAT|RIBOSOM|RRNA|AMINOACYL|TRNA", all_sets$gs_name, value = TRUE)
  cat(paste(sort(tr), collapse = "\n"), "\n")
}
