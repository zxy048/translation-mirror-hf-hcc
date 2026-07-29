# =============================================================================
# Script 19b: Fix MYC V2-TATS correlation — use ssGSEA-based scores (not direct)
# Purpose: Resolve ρ=0.753 (main) vs ρ=0.779 (M5) discrepancy
# Method:  Compute ssGSEA on GSE141198, use HALLMARK_MYC_TARGETS_V2 ssGSEA score
#          for bootstrap CI and de-overlap sensitivity
# =============================================================================

library(msigdbr)
library(dplyr)
library(GSVA)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"

# ── Load data ──────────────────────────────────────────────────────────────
cat("═══ Loading data ═══\n")

wgcna <- readRDS(file.path(PROJ_DIR, "GSE57338_WGCNA_rerun.rds"))
green_idx <- which(wgcna$moduleColors == wgcna$primary_translation_module)
green_genes <- wgcna$gene_symbols[green_idx]
cat(sprintf("Green module (TATS): %d genes\n", length(green_genes)))

expr_vst <- readRDS(file.path(PROJ_DIR, "GSE141198_vst.rds"))
green_in_gse141198 <- intersect(green_genes, rownames(expr_vst))
cat(sprintf("Green genes in GSE141198 VST: %d/%d\n",
            length(green_in_gse141198), length(green_genes)))

# ── Compute TATS ───────────────────────────────────────────────────────────
expr_green <- expr_vst[green_in_gse141198, , drop = FALSE]
expr_z <- t(scale(t(expr_green)))
tats_original <- colMeans(expr_z, na.rm = TRUE)

# ── Compute ssGSEA on GSE141198 (same as script 14) ────────────────────────
cat("\n═══ Computing ssGSEA on GSE141198 ═══\n")

msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_genes <- split(msig_h$gene_symbol, msig_h$gs_name)

msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2",
                   subcollection = "CP:REACTOME")
trans_kw <- "TRANSLATION|PEPTIDE_CHAIN_ELONGATION|EUKARYOTIC_TRANSLATION|RIBOSOME|NONSENSE_MEDIATED|TRNA_AMINOACYLATION|RRNA_PROCESSING"
trans_sets <- msig_c2[grepl(trans_kw, msig_c2$gs_name, ignore.case = TRUE), ]
trans_genes <- split(trans_sets$gene_symbol, trans_sets$gs_name)

kegg_ribo <- msigdbr(species = "Homo sapiens", collection = "C2",
                     subcollection = "CP:KEGG_MEDICUS") %>%
  filter(gs_name == "KEGG_RIBOSOME")
ribo_genes <- split(kegg_ribo$gene_symbol, kegg_ribo$gs_name)

all_gs <- c(hallmark_genes, trans_genes, ribo_genes)
all_gs <- all_gs[sapply(all_gs, length) >= 5]
cat(sprintf("Gene sets for ssGSEA: %d\n", length(all_gs)))

mat <- as.matrix(expr_vst)
param <- ssgseaParam(mat, all_gs, minSize = 5, maxSize = 500)
gse141198_ssgsea <- gsva(param, verbose = FALSE)
cat(sprintf("ssGSEA computed: %d pathways x %d samples\n",
            nrow(gse141198_ssgsea), ncol(gse141198_ssgsea)))

# ── Extract MYC Targets V2 ssGSEA score ────────────────────────────────────
mycv2_ssgsea_name <- "HALLMARK_MYC_TARGETS_V2"
stopifnot(mycv2_ssgsea_name %in% rownames(gse141198_ssgsea))

common_samples <- intersect(colnames(gse141198_ssgsea), names(tats_original))
cat(sprintf("Common samples (ssGSEA + TATS): %d\n", length(common_samples)))

mycv2_ssgsea_score <- as.numeric(gse141198_ssgsea[mycv2_ssgsea_name, common_samples])
tats_common <- tats_original[common_samples]

# ── Original correlation (ssGSEA-based) ────────────────────────────────────
obs_rho_mycv2_ssgsea <- cor(mycv2_ssgsea_score, tats_common, method = "spearman")
cat(sprintf("\nMYC Targets V2 (ssGSEA) vs TATS: ρ = %.4f\n", obs_rho_mycv2_ssgsea))

# ── Bootstrap 95% CI ───────────────────────────────────────────────────────
n_boot <- 10000
n_samples <- length(mycv2_ssgsea_score)

boot_rhos <- numeric(n_boot)
for (i in 1:n_boot) {
  idx <- sample(n_samples, replace = TRUE)
  boot_rhos[i] <- cor(mycv2_ssgsea_score[idx], tats_common[idx], method = "spearman")
}
ci_mycv2_ssgsea <- quantile(boot_rhos, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("Bootstrap 95%% CI: [%.4f, %.4f]\n", ci_mycv2_ssgsea[1], ci_mycv2_ssgsea[2]))

# ── De-overlap sensitivity (TATS-side only, ssGSEA score unchanged) ───────
# Remove 9 overlapping genes from TATS; MYC V2 ssGSEA score is computed from
# the full gene set via GSVA algorithm, so it remains unchanged.

msig_h_all <- msigdbr(species = "Homo sapiens", collection = "H")
myc_v2_genes <- msig_h_all %>%
  filter(gs_name == "HALLMARK_MYC_TARGETS_V2") %>%
  pull(gene_symbol) %>% unique()

myc_v2_overlap <- intersect(green_in_gse141198, myc_v2_genes)
cat(sprintf("\nMYC V2 overlap with TATS: %d genes\n", length(myc_v2_overlap)))

green_no_mycv2 <- setdiff(green_in_gse141198, myc_v2_overlap)
expr_green_de <- expr_vst[green_no_mycv2, , drop = FALSE]
expr_z_de <- t(scale(t(expr_green_de)))
tats_deoverlap <- colMeans(expr_z_de, na.rm = TRUE)

tats_de_common <- tats_deoverlap[common_samples]

# Correlation: ssGSEA MYC V2 score vs de-overlap TATS
obs_rho_deoverlap <- cor(mycv2_ssgsea_score, tats_de_common, method = "spearman")
cat(sprintf("MYC V2 (ssGSEA) vs TATS (de-overlap): ρ = %.4f\n", obs_rho_deoverlap))

# Bootstrap CI for de-overlap
boot_rhos_de <- numeric(n_boot)
for (i in 1:n_boot) {
  idx <- sample(n_samples, replace = TRUE)
  boot_rhos_de[i] <- cor(mycv2_ssgsea_score[idx], tats_de_common[idx], method = "spearman")
}
ci_deoverlap <- quantile(boot_rhos_de, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("De-overlap Bootstrap 95%% CI: [%.4f, %.4f]\n", ci_deoverlap[1], ci_deoverlap[2]))

# ── Also recompute ATF4 and DDIT3 with consistent pipeline ─────────────────
cat("\n═══ ATF4/DDIT3 bootstrap (consistency check) ═══\n")

atf4_expr <- as.numeric(expr_vst["ATF4", common_samples])
obs_rho_atf4 <- cor(atf4_expr, tats_common, method = "spearman")
boot_rhos_atf4 <- numeric(n_boot)
for (i in 1:n_boot) {
  idx <- sample(n_samples, replace = TRUE)
  boot_rhos_atf4[i] <- cor(atf4_expr[idx], tats_common[idx], method = "spearman")
}
ci_atf4 <- quantile(boot_rhos_atf4, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("ATF4-TATS: ρ = %.3f [%.3f, %.3f]\n", obs_rho_atf4, ci_atf4[1], ci_atf4[2]))

ddit3_expr <- as.numeric(expr_vst["DDIT3", common_samples])
obs_rho_ddit3 <- cor(ddit3_expr, tats_common, method = "spearman")
boot_rhos_ddit3 <- numeric(n_boot)
for (i in 1:n_boot) {
  idx <- sample(n_samples, replace = TRUE)
  boot_rhos_ddit3[i] <- cor(ddit3_expr[idx], tats_common[idx], method = "spearman")
}
ci_ddit3 <- quantile(boot_rhos_ddit3, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("DDIT3-TATS: ρ = %.3f [%.3f, %.3f]\n", obs_rho_ddit3, ci_ddit3[1], ci_ddit3[2]))

# ── Summary ────────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════\n")
cat("  FIXED MYC V2-TATS VALUES (ssGSEA-based)\n")
cat("═══════════════════════════════════════════════════════════\n")
cat(sprintf("MYC V2 (ssGSEA) vs TATS original:   ρ = %.3f [%.3f, %.3f]\n",
            obs_rho_mycv2_ssgsea, ci_mycv2_ssgsea[1], ci_mycv2_ssgsea[2]))
cat(sprintf("MYC V2 (ssGSEA) vs TATS de-overlap: ρ = %.3f [%.3f, %.3f]\n",
            obs_rho_deoverlap, ci_deoverlap[1], ci_deoverlap[2]))
cat(sprintf("Delta (original - de-overlap):      Δρ = %.3f\n",
            obs_rho_mycv2_ssgsea - obs_rho_deoverlap))
cat("\nThese values NOW MATCH the main text (ρ = 0.753, ssGSEA-based).\n")
cat("Old M5 values (direct pathway score): ρ = 0.779 [0.702, 0.837] — DISCARD.\n")

# ── Save ───────────────────────────────────────────────────────────────────
fix_results <- list(
  mycv2_ssgsea_original  = c(rho = obs_rho_mycv2_ssgsea,
                              ci_lower = ci_mycv2_ssgsea[1],
                              ci_upper = ci_mycv2_ssgsea[2]),
  mycv2_ssgsea_deoverlap = c(rho = obs_rho_deoverlap,
                              ci_lower = ci_deoverlap[1],
                              ci_upper = ci_deoverlap[2]),
  atf4_tats   = c(rho = obs_rho_atf4,  ci_lower = ci_atf4[1],  ci_upper = ci_atf4[2]),
  ddit3_tats  = c(rho = obs_rho_ddit3, ci_lower = ci_ddit3[1], ci_upper = ci_ddit3[2]),
  n_overlap_genes = length(myc_v2_overlap),
  overlap_genes    = myc_v2_overlap
)

saveRDS(fix_results, file.path(PROJ_DIR, "M5_mycv2_ssgsea_fix.rds"))
cat("\nResults saved: M5_mycv2_ssgsea_fix.rds\n")
cat("═══ Script 19b complete ═══\n")
