# Structure probe for the v2 figure generator (R17). No shipping artefact.
# Purpose: know exactly what each intermediate contains so R17 never
# subscript-guesses (the R15 bug class).
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

d <- function(nm) {
  f <- file.path(OUT, "intermediate", paste0(nm, ".rds"))
  if (!file.exists(f)) { cat(sprintf("\n### %s : MISSING\n", nm)); return(invisible()) }
  x <- readRDS(f)
  cat(sprintf("\n### %s : class=%s\n", nm, paste(class(x), collapse = "/")))
  if (is.list(x)) {
    for (n in names(x)) {
      v <- x[[n]]
      cat(sprintf("  $%-22s %-18s ", n, paste(class(v), collapse = "/")))
      if (is.data.frame(v)) {
        cat(sprintf("dim %dx%d | cols: %s\n", nrow(v), ncol(v),
                    paste(head(names(v), 14), collapse = ", ")))
      } else if (is.matrix(v)) {
        cat(sprintf("matrix %dx%d | rn: %s | cn: %s\n", nrow(v), ncol(v),
                    paste(head(rownames(v), 3), collapse = ","),
                    paste(head(colnames(v), 6), collapse = ",")))
      } else if (is.atomic(v)) {
        cat(sprintf("len %d | %s\n", length(v),
                    paste(head(as.character(v), 6), collapse = ", ")))
      } else cat("(other)\n")
    }
  }
}

for (nm in c("R02_wgcna", "R03_translation_module", "R04_preservation",
             "R06_ssgsea", "R07_effect_sizes", "R08_redundancy", "R09_mirror",
             "R10_heterogeneity", "R10b_hf_aetiology", "R12_tats_tf",
             "R15_cirrhosis_control")) d(nm)

cat("\n### R01_universe\n")
u <- readRDS(file.path(OUT, "intermediate", "R01_universe.rds"))
cat("  class:", class(u), "| names:", paste(names(u), collapse = ", "), "\n")
if (is.list(u)) for (n in names(u)) {
  v <- u[[n]]
  cat(sprintf("  $%-18s %-14s %s\n", n, paste(class(v), collapse = "/"),
              if (is.atomic(v)) paste0("len ", length(v), " | ", paste(head(as.character(v), 4), collapse = ", ")) else ""))
}

cat("\n### SETS_TIER_A\n")
cat("  n =", length(SETS_TIER_A), "\n")
cat("  ", paste(SETS_TIER_A, collapse = "\n   "), "\n")
cat("\nPROBE DONE\n")
