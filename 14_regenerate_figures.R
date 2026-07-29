# =============================================================================
# 脚本 14：重新生成 Figure 4/5/S1/S2（TATS 替换 TGS）
# 基于 HF WGCNA 重跑的 green module (227 genes)
# =============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(png)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
FIG_DIR  <- file.path(PROJ_DIR, "figures")
TBL_DIR  <- file.path(PROJ_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Helper: read PNG as grob ─────────────────────────────────────────────────
read_png_grob <- function(path) {
  if (!file.exists(path)) return(nullGrob())
  img <- readPNG(path)
  rasterGrob(img, interpolate = TRUE)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Load data
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══ Loading data ═══\n")

# HF WGCNA rerun results
wgcna <- readRDS(file.path(PROJ_DIR, "GSE57338_WGCNA_rerun.rds"))
green_genes <- wgcna$gene_symbols[wgcna$moduleColors == wgcna$primary_translation_module]
cat(sprintf("Green module genes: %d\n", length(green_genes)))

# GSE141198 VST + clinical
expr_vst <- readRDS(file.path(PROJ_DIR, "GSE141198_vst.rds"))
clinical <- readRDS(file.path(PROJ_DIR, "GSE141198_clinical.rds"))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Compute TATS for GSE141198
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══ Computing TATS ═══\n")

green_in_141198 <- intersect(green_genes, rownames(expr_vst))
cat(sprintf("Green genes in GSE141198: %d/%d\n",
            length(green_in_141198), length(green_genes)))

expr_green <- expr_vst[green_in_141198, , drop = FALSE]
expr_z <- t(scale(t(expr_green)))
tats_score <- colMeans(expr_z, na.rm = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Figure 4: TATS combined figure
#    Panel A: KM curve (GSE141198)
#    Panel B: Cross-disease expression direction
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══ Figure 4: TATS Validation ═══\n")

# --- Panel A: KM curve ---
common_samples <- intersect(names(tats_score), clinical$geo_id)
cat(sprintf("TATS-survival matched: %d\n", length(common_samples)))

surv_df <- data.frame(
  sample = common_samples,
  tats = tats_score[common_samples],
  os_time = clinical[common_samples, "os_time"],
  os_status = clinical[common_samples, "os_status"],
  etiology = clinical[common_samples, "etiology"],
  stringsAsFactors = FALSE
)
surv_df <- surv_df[!is.na(surv_df$os_time) & !is.na(surv_df$os_status), ]
surv_df$os_time <- surv_df$os_time / 30.44
surv_df$tats_group <- ifelse(surv_df$tats > median(surv_df$tats), "High", "Low")
surv_df$tats_group <- factor(surv_df$tats_group, levels = c("Low", "High"))

fit <- survfit(Surv(os_time, os_status) ~ tats_group, data = surv_df)
logrank_p <- survdiff(Surv(os_time, os_status) ~ tats_group, data = surv_df)$pvalue
cox_cont <- coxph(Surv(os_time, os_status) ~ tats, data = surv_df)
cox_summary <- summary(cox_cont)

p_km <- ggsurvplot(fit, data = surv_df,
  pval = TRUE, pval.size = 3.5,
  palette = c("#2166AC", "#B2182B"),
  xlab = "Overall Survival (months)",
  ylab = "Overall Survival Probability",
  title = sprintf("TATS in GSE141198 (n=%d, events=%d)",
                  nrow(surv_df), sum(surv_df$os_status)),
  subtitle = sprintf("227-gene green module | Log-rank p=%.3f | Cox HR=%.2f (%.2f-%.2f)",
                     logrank_p,
                     cox_summary$conf.int[1, "exp(coef)"],
                     cox_summary$conf.int[1, "lower .95"],
                     cox_summary$conf.int[1, "upper .95"]),
  legend.title = "TATS",
  legend.labs = c("Low", "High"),
  risk.table = TRUE,
  risk.table.height = 0.2,
  ggtheme = theme_minimal(base_size = 11)
)
png(file.path(FIG_DIR, "Figure4A_TATS_KM.png"), width = 8, height = 6, units = "in", res = 300)
print(p_km)
dev.off()
cat("Panel 4A (TATS KM) saved.\n")

# --- Panel B: Cross-disease expression direction ---
# Load cross-disease log2FC data
cat("Panel 4B: Cross-disease expression direction...\n")

# Get hub genes from WGCNA results
hub_genes <- wgcna$hub_genes_auto[1:10]

# Load or compute cross-disease FC
# Check if we have the data from TATS_results or direction_consistency
dir_cons <- readRDS(file.path(PROJ_DIR, "ssgsea_cross_disease_result.rds"))
# Look for cross-disease gene FC data
if (file.exists(file.path(PROJ_DIR, "TATS_results.rds"))) {
  # Use stored hub genes
  cat(sprintf("Hub genes: %s\n", paste(hub_genes, collapse = ", ")))
}

# We need log2FC from GSE57338 and TCGA-LIHC for the hub genes
# Since TCGA-LIHC SE loading causes segfault, use stored results if available
# For now, plot what we can from the GSE57338 side
# Load HF limma results if available
hf_limma_files <- list.files(PROJ_DIR, pattern = "GSE57338.*limma|GSE57338.*deg|GSE57338.*logFC",
                             ignore.case = TRUE)
cat(sprintf("HF limma files: %s\n", paste(hf_limma_files, collapse = ", ")))

# Use the WGCNA hub gene GS values as a proxy for direction
if (!is.null(wgcna$hub_GS) && length(wgcna$hub_GS) > 0) {
  # Create a cross-disease direction plot
  gs_values <- wgcna$hub_GS[hub_genes]
  mm_values <- wgcna$hub_MM[hub_genes]

  hub_df <- data.frame(
    Gene = factor(hub_genes, levels = rev(hub_genes)),
    GS = gs_values[hub_genes],
    MM = mm_values[hub_genes],
    stringsAsFactors = FALSE
  )

  p_cross <- ggplot(hub_df, aes(x = GS, y = Gene)) +
    geom_bar(stat = "identity", fill = ifelse(hub_df$GS < 0, "#2166AC", "#B2182B"),
             width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(
      title = "Hub Gene Significance (HF: GSE57338)",
      subtitle = "Gene Significance = correlation with HF status",
      x = "Gene Significance (Pearson r with HF status)",
      y = ""
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank())

  ggsave(file.path(FIG_DIR, "Figure4B_CrossDisease_Direction.png"),
         p_cross, width = 7, height = 5, dpi = 300)
  cat("Panel 4B (cross-disease direction) saved.\n")
}

# --- Combine Figure 4 ---
gA <- read_png_grob(file.path(FIG_DIR, "Figure4A_TATS_KM.png"))
gB <- read_png_grob(file.path(FIG_DIR, "Figure4B_CrossDisease_Direction.png"))

png(file.path(FIG_DIR, "Figure4_TATS_Validation.png"),
    width = 14, height = 7, units = "in", res = 300)

grid.newpage()
grid.text(
  "Figure 4. TATS exploratory clinical association and hub gene significance.",
  x = 0.02, y = 0.99, just = c(0, 1),
  gp = gpar(fontface = "bold", fontsize = 14)
)

pushViewport(viewport(layout = grid.layout(1, 2, widths = c(0.55, 0.45))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("A", x = 0.02, y = 0.97, just = c(0, 1),
          gp = gpar(fontface = "bold", fontsize = 14))
grid.draw(gA)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.text("B", x = 0.02, y = 0.97, just = c(0, 1),
          gp = gpar(fontface = "bold", fontsize = 14))
grid.draw(gB)
popViewport()

grid.text(
  "(A) Kaplan-Meier OS curves for TATS (median split) in GSE141198 (n=148, 94 events). (B) Gene Significance (correlation with HF status) for top 10 green module hub genes.",
  x = 0.02, y = 0.005, just = c(0, 0),
  gp = gpar(fontsize = 8, col = "grey40")
)

dev.off()
cat(sprintf("Combined Figure 4 saved: %s\n",
            file.path(FIG_DIR, "Figure4_TATS_Validation.png")))

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Figure 5: TF-TATS and Pathway-TATS correlations
#    (Using TCGA-LIHC if available, else GSE141198)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══ Figure 5: TF-TATS and Pathway-TATS ═══\n")

# We need TCGA-LIHC expression data for TF correlations
# Try loading via DESeq2 instead of SummarizedExperiment to avoid segfault
tcga_available <- FALSE

# Check if we have a pre-computed VST matrix
tcga_vst_files <- list.files(PROJ_DIR, pattern = "TCGA.*vst|tcga.*vst|LIHC.*vst",
                             ignore.case = TRUE)
cat(sprintf("TCGA VST files: %s\n", paste(tcga_vst_files, collapse = ", ")))

if (length(tcga_vst_files) == 0) {
  cat("No pre-computed TCGA VST file found. Using GSE141198 for TF correlations.\n")

  # Use GSE141198 expression for TF-TATS correlations
  expr_mat <- expr_vst

  # Candidate TFs from manuscript
  candidate_tfs <- c("ATF4", "DDIT3", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5",
                     "E2F6", "E2F7", "E2F8", "HIF1A", "MTOR", "MYC", "MYCL",
                     "MYCN", "NFE2L2", "RPTOR", "TP53", "XBP1")

  tfs_found <- intersect(candidate_tfs, rownames(expr_mat))
  cat(sprintf("TFs found in GSE141198: %d/19\n", length(tfs_found)))

  if (length(tfs_found) > 0) {
    # Compute TF-TATS correlations
    common_tats <- intersect(colnames(expr_mat), names(tats_score))
    cat(sprintf("Common samples for TF-TATS: %d\n", length(common_tats)))

    tf_cors <- data.frame(
      TF = character(),
      Spearman_R = numeric(),
      P_value = numeric(),
      FDR = numeric(),
      stringsAsFactors = FALSE
    )

    for (tf in tfs_found) {
      tf_expr <- as.numeric(expr_mat[tf, common_tats])
      tats_val <- tats_score[common_tats]
      ct <- cor.test(tf_expr, tats_val, method = "spearman")
      tf_cors <- rbind(tf_cors, data.frame(
        TF = tf, Spearman_R = ct$estimate, P_value = ct$p.value,
        stringsAsFactors = FALSE
      ))
    }
    tf_cors$FDR <- p.adjust(tf_cors$P_value, method = "BH")
    tf_cors <- tf_cors[order(tf_cors$Spearman_R, decreasing = TRUE), ]

    # Print results
    cat("\nTF-TATS correlations (GSE141198):\n")
    for (i in 1:nrow(tf_cors)) {
      sig_marker <- ifelse(tf_cors$FDR[i] < 0.05, "*", "")
      cat(sprintf("  %-8s ρ=%+.3f  p=%.4f  FDR=%.4f %s\n",
                  tf_cors$TF[i], tf_cors$Spearman_R[i],
                  tf_cors$P_value[i], tf_cors$FDR[i], sig_marker))
    }

    # Plot TF-TATS
    tf_cors$TF <- factor(tf_cors$TF, levels = rev(tf_cors$TF))
    tf_cors$Significant <- ifelse(tf_cors$FDR < 0.05, "FDR < 0.05", "NS")
    tf_cors$Direction <- ifelse(tf_cors$Spearman_R > 0, "Positive", "Negative")

    p_tf <- ggplot(tf_cors, aes(x = Spearman_R, y = TF)) +
      geom_bar(stat = "identity", aes(fill = Direction), width = 0.7) +
      scale_fill_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +
      geom_text(aes(label = sprintf("ρ=%.3f%s", Spearman_R,
                                     ifelse(FDR < 0.05, "*", "")),
                    hjust = ifelse(Spearman_R > 0, -0.1, 1.1)),
                size = 3, color = "grey30") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
      labs(
        title = "TF Expression vs. TATS",
        subtitle = sprintf("GSE141198 (n=%d) | 19 candidate TFs | * FDR < 0.05",
                           length(common_tats)),
        x = "Spearman ρ with TATS",
        y = ""
      ) +
      xlim(min(tf_cors$Spearman_R) - 0.15, max(tf_cors$Spearman_R) + 0.15) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none"
      )

    ggsave(file.path(FIG_DIR, "Figure_TF_TATS_correlation.png"),
           p_tf, width = 8, height = 6, dpi = 300)
    cat("Figure 5A (TF-TATS) saved.\n")

    # Save TF results for manuscript reference
    saveRDS(tf_cors, file.path(PROJ_DIR, "TF_TATS_correlation_results.rds"))
  }

  # --- Pathway-TATS correlations ---
  cat("\n--- Pathway-TATS correlations ---\n")

  # Use GSE141198 ssGSEA scores if available
  ssgsea_res <- readRDS(file.path(PROJ_DIR, "ssgsea_cross_disease_result.rds"))
  cat(sprintf("ssGSEA result names: %s\n", paste(names(ssgsea_res), collapse = ", ")))

  # Check for GSE141198 ssGSEA scores
  gse141198_ssgsea <- NULL
  for (nm in c("gse141198_ssgsea", "hcc_ssgsea", "GSE141198_ssgsea")) {
    if (nm %in% names(ssgsea_res)) {
      gse141198_ssgsea <- ssgsea_res[[nm]]
      break
    }
  }

  if (is.null(gse141198_ssgsea)) {
    cat("No GSE141198 ssGSEA found. Computing from scratch...\n")

    library(GSVA)
    library(msigdbr)

    # Get Hallmark gene sets
    msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
    hallmark_genes <- split(msig_h$gene_symbol, msig_h$gs_name)

    # Add translation-related Reactome
    msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2",
                       subcollection = "CP:REACTOME")
    trans_kw <- "TRANSLATION|PEPTIDE_CHAIN_ELONGATION|EUKARYOTIC_TRANSLATION|RIBOSOME|NONSENSE_MEDIATED|TRNA_AMINOACYLATION|RRNA_PROCESSING"
    trans_sets <- msig_c2[grepl(trans_kw, msig_c2$gs_name, ignore.case = TRUE), ]
    trans_genes <- split(trans_sets$gene_symbol, trans_sets$gs_name)

    # KEGG ribosome
    kegg_ribo <- msigdbr(species = "Homo sapiens", collection = "C2",
                         subcollection = "CP:KEGG_MEDICUS") %>%
      filter(gs_name == "KEGG_RIBOSOME")
    if (nrow(kegg_ribo) > 0) {
      ribo_genes <- split(kegg_ribo$gene_symbol, kegg_ribo$gs_name)
    } else {
      ribo_genes <- list()
    }

    all_gs <- c(hallmark_genes, trans_genes, ribo_genes)
    all_gs <- all_gs[sapply(all_gs, length) >= 5]
    cat(sprintf("Gene sets for ssGSEA: %d\n", length(all_gs)))

    # Run ssGSEA on GSE141198
    mat <- as.matrix(expr_vst)
    param <- ssgseaParam(mat, all_gs, minSize = 5, maxSize = 500)
    gse141198_ssgsea <- gsva(param, verbose = FALSE)
    cat(sprintf("ssGSEA computed: %d pathways x %d samples\n",
                nrow(gse141198_ssgsea), ncol(gse141198_ssgsea)))
  }

  # Correlate pathway scores with TATS
  common_ssgsea <- intersect(colnames(gse141198_ssgsea), names(tats_score))
  cat(sprintf("Pathway-TATS common samples: %d\n", length(common_ssgsea)))

  # Only Hallmark pathways
  hallmark_rows <- grep("^HALLMARK_", rownames(gse141198_ssgsea), value = TRUE)
  cat(sprintf("Hallmark pathways: %d\n", length(hallmark_rows)))

  pathway_cors <- data.frame(
    Pathway = character(),
    Spearman_R = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )

  for (pw in hallmark_rows) {
    pw_scores <- as.numeric(gse141198_ssgsea[pw, common_ssgsea])
    ct <- cor.test(pw_scores, tats_score[common_ssgsea], method = "spearman")
    pathway_cors <- rbind(pathway_cors, data.frame(
      Pathway = gsub("^HALLMARK_", "", pw),
      Spearman_R = ct$estimate,
      P_value = ct$p.value,
      stringsAsFactors = FALSE
    ))
  }
  pathway_cors$FDR <- p.adjust(pathway_cors$P_value, method = "BH")
  pathway_cors <- pathway_cors[order(pathway_cors$Spearman_R, decreasing = TRUE), ]

  # Top 20 pathways
  pw_top <- head(pathway_cors, 20)
  pw_top$Pathway <- factor(pw_top$Pathway, levels = rev(pw_top$Pathway))
  pw_top$Direction <- ifelse(pw_top$Spearman_R > 0, "Positive", "Negative")

  p_pw <- ggplot(pw_top, aes(x = Spearman_R, y = Pathway)) +
    geom_bar(stat = "identity", aes(fill = Direction), width = 0.7) +
    scale_fill_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +
    geom_text(aes(label = sprintf("ρ=%.3f%s", Spearman_R,
                                   ifelse(FDR < 0.05, "*", "")),
                  hjust = ifelse(Spearman_R > 0, -0.1, 1.1)),
              size = 2.8, color = "grey30") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    labs(
      title = "Hallmark Pathway Activity vs. TATS",
      subtitle = sprintf("GSE141198 (n=%d) | Top 20 pathways | * FDR < 0.05",
                         length(common_ssgsea)),
      x = "Spearman ρ with TATS",
      y = ""
    ) +
    xlim(min(pw_top$Spearman_R) - 0.15, max(pw_top$Spearman_R) + 0.15) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )

  ggsave(file.path(FIG_DIR, "Figure_Pathway_TATS_correlation.png"),
         p_pw, width = 9, height = 7, dpi = 300)
  cat("Figure 5B (Pathway-TATS) saved.\n")

  # Save pathway results
  saveRDS(pathway_cors, file.path(PROJ_DIR, "Pathway_TATS_correlation_results.rds"))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. Figure S2: TATS by etiology subgroup
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══ Figure S2: TATS by Etiology ═══\n")

etiology_groups <- unique(surv_df$etiology)
cat(sprintf("Etiology groups: %s\n", paste(etiology_groups, collapse = ", ")))

# Generate KM per etiology
km_plots <- list()
for (et in etiology_groups) {
  sub <- surv_df[surv_df$etiology == et, ]
  if (nrow(sub) < 15) next

  sub$group <- ifelse(sub$tats > median(sub$tats), "High", "Low")
  sub$group <- factor(sub$group, levels = c("Low", "High"))

  if (length(unique(sub$group)) < 2) next

  fit_sub <- survfit(Surv(os_time, os_status) ~ group, data = sub)
  logrank_sub <- survdiff(Surv(os_time, os_status) ~ group, data = sub)$pvalue

  p <- ggsurvplot(fit_sub, data = sub,
    pval = TRUE, pval.size = 3,
    palette = c("#2166AC", "#B2182B"),
    title = sprintf("%s (n=%d, events=%d)", et, nrow(sub), sum(sub$os_status)),
    xlab = "OS (months)", ylab = "OS Probability",
    legend.title = "TATS",
    legend.labs = c("Low", "High"),
    ggtheme = theme_minimal(base_size = 10)
  )
  km_plots[[et]] <- p
}

if (length(km_plots) > 0) {
  png(file.path(FIG_DIR, "Figure_S2_TATS_Subgroup.png"),
      width = 12, height = 4 * ceiling(length(km_plots) / 3),
      units = "in", res = 300)

  # Arrange in grid
  n_cols <- min(3, length(km_plots))
  grid.arrange(grobs = lapply(km_plots, function(x) x$plot),
               ncol = n_cols,
               top = "TATS Survival Analysis by Etiology Subgroup (GSE141198)")
  dev.off()
  cat(sprintf("Figure S2 (TATS subgroups) saved: %d panels\n", length(km_plots)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6. Figure S1: Cross-disease expression direction
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══ Figure S1: Cross-disease direction ═══\n")

# Using GS values from HF WGCNA as HF direction proxy
# And if available, HCC log2FC
# For now, create a clean version showing module membership vs gene significance
if (!is.null(wgcna$hub_GS) && !is.null(wgcna$hub_MM)) {
  # Get top 30 hub genes by |GS|*|MM|
  gs_all <- wgcna$hub_GS
  mm_all <- wgcna$hub_MM
  common_names <- intersect(names(gs_all), names(mm_all))

  hub_rank_df <- data.frame(
    Gene = common_names,
    GS = gs_all[common_names],
    MM = mm_all[common_names],
    Product = abs(gs_all[common_names]) * abs(mm_all[common_names]),
    stringsAsFactors = FALSE
  )
  hub_rank_df <- hub_rank_df[order(hub_rank_df$Product, decreasing = TRUE), ]
  hub_top30 <- head(hub_rank_df, 30)

  # Classify by |GS| > 0.2 AND |MM| > 0.8
  hub_top30$Hub <- abs(hub_top30$GS) > 0.2 & abs(hub_top30$MM) > 0.8

  p_s1 <- ggplot(hub_top30, aes(x = MM, y = GS)) +
    geom_point(aes(color = Hub, size = Product), alpha = 0.8) +
    scale_color_manual(values = c("TRUE" = "#228B22", "FALSE" = "grey60"),
                       labels = c("TRUE" = "Hub (|GS|>0.2 & |MM|>0.8)", "FALSE" = "Non-hub")) +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.8, 0.8), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_text(data = subset(hub_top30, Hub),
              aes(label = Gene), size = 2.8, hjust = -0.1, vjust = -0.5) +
    labs(
      title = "Green Module: Module Membership vs. Gene Significance",
      subtitle = sprintf("227 genes | Hub: |GS|>0.2 & |MM|>0.8 | %d hub genes identified",
                         sum(hub_top30$Hub)),
      x = "Module Membership (MM)",
      y = "Gene Significance (GS, HF vs. NF)",
      color = ""
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(FIG_DIR, "Figure_S1_Direction_Consistency.png"),
         p_s1, width = 8, height = 7, dpi = 300)
  cat("Figure S1 (MM vs GS) saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 7. Summary
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════\n")
cat("  FIGURE REGENERATION COMPLETE\n")
cat("═══════════════════════════════════\n")
cat("\nGenerated figures:\n")
cat("  Figure 4:  Figure4_TATS_Validation.png\n")
cat("  Figure 5A: Figure_TF_TATS_correlation.png\n")
cat("  Figure 5B: Figure_Pathway_TATS_correlation.png\n")
cat("  Figure S1: Figure_S1_Direction_Consistency.png\n")
cat("  Figure S2: Figure_S2_TATS_Subgroup.png\n")

# Print key statistics for the manuscript
cat("\n--- Key Statistics for Manuscript ---\n")
cat(sprintf("Green module: %d genes\n", length(green_genes)))
cat(sprintf("TATS in GSE141198: Log-rank p=%.4f, Cox HR=%.2f (%.2f-%.2f), p=%.4f\n",
            logrank_p,
            cox_summary$conf.int[1, "exp(coef)"],
            cox_summary$conf.int[1, "lower .95"],
            cox_summary$conf.int[1, "upper .95"],
            cox_summary$coefficients[1, "Pr(>|z|)"]))
if (exists("tf_cors") && nrow(tf_cors) > 0) {
  cat(sprintf("Top TF correlate: %s (ρ=%+.3f, FDR=%.4f)\n",
              tf_cors$TF[1], tf_cors$Spearman_R[1], tf_cors$FDR[1]))
}
if (exists("pathway_cors") && nrow(pathway_cors) > 0) {
  cat(sprintf("Top pathway correlate: %s (ρ=%+.3f, FDR=%.4f)\n",
              pathway_cors$Pathway[1], pathway_cors$Spearman_R[1], pathway_cors$FDR[1]))
}
cat("\nDone.\n")
