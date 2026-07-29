# Minimal script: Generate module-trait heatmap for GSE57338 HF WGCNA
# Uses cached MEs from original WGCNA run (v3: directly use stored MEs)
# Avoids recomputing MEs → guarantees r=0.794 matches original result

library(ggplot2)

PROJ_DIR <- "D:/R_projects/revision_analysis"
FIG_DIR  <- file.path(PROJ_DIR, "figures")

# ── Load cached WGCNA result ──
wgcna <- readRDS(file.path(PROJ_DIR, "GSE57338_WGCNA_result.rds"))
MEs <- wgcna$MEs
mod_names <- gsub("^ME", "", colnames(MEs))
cat(sprintf("Loaded MEs: %d samples x %d modules\n", nrow(MEs), ncol(MEs)))
cat(sprintf("Modules: %s\n", paste(mod_names, collapse = ", ")))

# ── Read phenotype from GEO series matrix (base R, no getGEO) ──
cat("Reading GSE57338 phenotype...\n")
gz_path <- file.path(PROJ_DIR, "GSE57338_series_matrix.txt.gz")

con <- gzfile(gz_path, "r")
all_lines <- readLines(con)
close(con)

# Find the "heart failure: yes/no" characteristic
# GSE57338 has 5 characteristics lines, all tagged !Sample_characteristics_ch1.
# Each line = one characteristic dimension across all samples.
hf_status <- NULL
char_lines <- all_lines[grepl("Sample_characteristics_ch1", all_lines, fixed = TRUE)]
for (cl in char_lines) {
  parts <- strsplit(cl, "\t")[[1]][-1]  # drop tag
  parts <- gsub('"', '', parts)
  if (any(grepl("heart.failure|heart failure", parts, ignore.case = TRUE))) {
    hf_status <- parts
    cat(sprintf("Found HF status in characteristics (line with values: %s)\n",
                paste(unique(parts)[1:min(2, length(unique(parts)))], collapse=", ")))
    break
  }
}

if (is.null(hf_status)) stop("Could not find heart failure status in series matrix")

cat(sprintf("Samples: %d, Unique: %s\n", length(hf_status),
            paste(unique(hf_status), collapse = ", ")))

# Binary HF indicator
is_hf <- grepl(": *yes|heart failure.*yes", hf_status, ignore.case = TRUE)
hf_binary <- ifelse(is_hf, 1, 0)
cat(sprintf("HF=%d, NF=%d\n", sum(is_hf), sum(!is_hf)))

# ── Compute module-trait correlations ──
cat("Computing module-trait correlations...\n")
mod_trait_cor <- cor(MEs, hf_binary, use = "pairwise.complete.obs")
mod_trait_p <- apply(MEs, 2, function(me) {
  cor.test(me, hf_binary, method = "pearson")$p.value
})

cor_df <- data.frame(
  Module = mod_names,
  Correlation = as.numeric(mod_trait_cor),
  P_value = mod_trait_p,
  stringsAsFactors = FALSE
)
cor_df <- cor_df[order(abs(cor_df$Correlation), decreasing = TRUE), ]
cor_df$Module <- factor(cor_df$Module, levels = cor_df$Module)

# Flag black module
cor_df$is_black <- cor_df$Module == "black"

black_r <- cor_df$Correlation[cor_df$Module == "black"]
black_p <- cor_df$P_value[cor_df$Module == "black"]
cat(sprintf("Black module: r=%.3f, p=%.2e\n", black_r, black_p))

# ── Bar plot ──
p <- ggplot(cor_df, aes(x = Correlation, y = Module)) +
  geom_bar(aes(fill = is_black), stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#C62828", "FALSE" = "grey60"), guide = "none") +
  geom_text(aes(label = sprintf("r=%.2f", Correlation),
                hjust = ifelse(Correlation > 0, -0.1, 1.1)),
            size = 2.8, color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "GSE57338 HF: Module–Trait (HF vs NF) Correlation",
    subtitle = paste0("Signed WGCNA, beta=12 | Black module (ribosome/translation) r=",
                      round(black_r, 3), ", p=",
                      format(black_p, digits = 2, scientific = TRUE)),
    x = "Pearson Correlation with HF Status",
    y = "Module"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "Figure2_ModuleTrait_Heatmap.png"),
       p, width = 8, height = max(5, nrow(cor_df) * 0.38), dpi = 300,
       limitsize = FALSE)
cat("Done: Figure2_ModuleTrait_Heatmap.png\n")
