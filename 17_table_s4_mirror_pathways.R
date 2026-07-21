# =============================================================================
# Script 17: Generate Supplementary Table S4 — 24/33 Mirror Pathways
# =============================================================================

PROJ_DIR <- "D:/R_projects/revision_analysis"

ssgsea <- readRDS(file.path(PROJ_DIR, "ssgsea_cross_disease_result.rds"))
plot_df <- ssgsea$plot_df

# Filter translation pathways
trans_df <- plot_df[plot_df$is_translation == TRUE, ]

# Determine mirror direction: HCC↑ (es_hcc > 0) and HF↓ (es_hf < 0)
trans_df$mirror <- trans_df$es_hcc > 0 & trans_df$es_hf < 0
trans_df$direction_hcc <- ifelse(trans_df$es_hcc > 0, "Up", "Down")
trans_df$direction_hf  <- ifelse(trans_df$es_hf > 0, "Up", "Down")

# Determine source database
trans_df$source <- ifelse(grepl("^HALLMARK_", trans_df$pathway), "Hallmark",
                   ifelse(grepl("^KEGG_", trans_df$pathway), "KEGG",
                   ifelse(grepl("^REACTOME_", trans_df$pathway), "Reactome", "Other")))

# Clean pathway names
trans_df$pathway_clean <- gsub("^HALLMARK_|^KEGG_|^REACTOME_", "", trans_df$pathway)
trans_df$pathway_clean <- gsub("_", " ", trans_df$pathway_clean)

# Order by mirror status then effect size
trans_df <- trans_df[order(-trans_df$mirror, -abs(trans_df$es_hcc - trans_df$es_hf)), ]

# Build Table S4
tbl_s4 <- data.frame(
  Pathway = trans_df$pathway_clean,
  Source = trans_df$source,
  HCC_Effect_d = round(trans_df$es_hcc, 3),
  HCC_Direction = trans_df$direction_hcc,
  HF_Effect_d = round(trans_df$es_hf, 3),
  HF_Direction = trans_df$direction_hf,
  Mirror = ifelse(trans_df$mirror, "Yes", "No"),
  row.names = NULL
)

cat(sprintf("Total translation pathways: %d\n", nrow(tbl_s4)))
cat(sprintf("Mirror pathways (HCC Up, HF Down): %d\n", sum(tbl_s4$Mirror == "Yes")))
cat(sprintf("Non-mirror pathways: %d\n", sum(tbl_s4$Mirror == "No")))

# Save as CSV
write.csv(tbl_s4, file.path(PROJ_DIR, "Table_S4_mirror_pathways.csv"), row.names = FALSE)

# Print preview
cat("\n=== Table S4 Preview (first 10 rows) ===\n")
print(head(tbl_s4, 10))

cat("\n=== Non-mirror pathways ===\n")
print(tbl_s4[tbl_s4$Mirror == "No", ])

cat(sprintf("\nTable S4 saved: %s\n", file.path(PROJ_DIR, "Table_S4_mirror_pathways.csv")))
