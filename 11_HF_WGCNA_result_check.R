# =============================================================================
# 脚本 11：HF WGCNA 重跑结果检查
# 快速回答：有没有翻译模块？该走路线 1 还是路线 2？
# =============================================================================

PROJ_DIR <- "D:/R_projects/revision_analysis"

wgcna <- readRDS(file.path(PROJ_DIR, "GSE57338_WGCNA_rerun.rds"))

cat("============================================\n")
cat("  HF WGCNA RERUN — RESULT CHECK\n")
cat("============================================\n\n")

# ── 1. 基本参数 ──
cat(sprintf("Input genes: %d\n", wgcna$n_genes_input))
cat(sprintf("Beta: %d (R^2=%.3f)\n", wgcna$beta_sel, wgcna$beta_r2))
cat(sprintf("Modules: %d\n", wgcna$n_modules))
cat(sprintf("Ribosomal genes in input: %d/%d\n",
            wgcna$ribo_genes_in_input, wgcna$ribo_genes_total))

# ── 2. 翻译模块是否存在？──
primary_mod <- wgcna$primary_translation_module
cat(sprintf("\nPrimary translation module: %s\n",
            if (is.null(primary_mod)) "NONE" else primary_mod))

if (!is.null(primary_mod)) {
  cat(sprintf("Translation module size: %d genes\n",
              sum(wgcna$moduleColors == primary_mod)))
  cat(sprintf("Translation score: %.1f\n",
              wgcna$module_translation_score[primary_mod]))
}

# ── 3. Module-trait correlations ──
cat("\n--- Module-Trait Correlations ---\n")
print(wgcna$module_trait_cor[, c("Module", "Correlation", "P_value", "Gene_Count")])

# ── 4. Hub gene audit ──
cat("\n--- Manuscript 7 Hub Genes ---\n")
print(wgcna$hub_manuscript_audit)

# ── 5. Auto-detected hub genes ──
cat("\n--- Auto-detected Hub Genes ---\n")
if (length(wgcna$hub_genes_auto) > 0) {
  for (g in wgcna$hub_genes_auto) {
    cat(sprintf("  %s: GS=%.3f, MM=%.3f\n", g,
                wgcna$hub_GS[g], wgcna$hub_MM[g]))
  }
} else {
  cat("  None (no translation module detected)\n")
}

# ── 6. Ribosomal gene distribution across modules ──
cat("\n--- Ribosomal Gene Distribution ---\n")
ribo_dist <- wgcna$ribo_genes_per_module
ribo_dist <- ribo_dist[order(ribo_dist, decreasing = TRUE)]
print(ribo_dist[ribo_dist > 0])

# ── 7. Cross-disease comparison ──
if (!is.null(wgcna$cross_disease)) {
  cat("\n--- Cross-Disease ---\n")
  cd <- wgcna$cross_disease
  cat(sprintf("HF translation genes in HCC WGCNA: %d\n", length(cd$hf_trans_in_hcc)))
  cat(sprintf("Fisher: OR=%.1f, p=%.2e\n",
              cd$fisher_test$odds_ratio, cd$fisher_test$p_value))
  cat("HCC module distribution:\n")
  print(cd$hcc_module_distribution)
}

# ── 8. GO top terms ──
if (!is.null(primary_mod) && primary_mod %in% names(wgcna$module_go_results)) {
  cat("\n--- Top 5 GO Terms for", primary_mod, "---\n")
  ego <- wgcna$module_go_results[[primary_mod]]
  top5 <- head(as.data.frame(ego)[order(as.data.frame(ego)$p.adjust), ], 5)
  for (i in 1:nrow(top5)) {
    cat(sprintf("  %s | p.adj=%.2e | %s\n",
                top5$ID[i], top5$p.adjust[i], top5$Description[i]))
  }
}

# ── 9. 结论 ──
cat("\n============================================\n")
cat("  CONCLUSION\n")
cat("============================================\n")

if (!is.null(primary_mod)) {
  trans_size <- sum(wgcna$moduleColors == primary_mod)
  trans_r <- wgcna$module_trait_cor$Correlation[
    wgcna$module_trait_cor$Module == primary_mod]

  cat(sprintf("Translation module FOUND: %s (%d genes, r=%.3f)\n",
              primary_mod, trans_size, trans_r))
  cat("\nROUTE 1: Keep module story, revise manuscript.\n")
  cat("  - Update all figures with new module assignments\n")
  cat("  - Replace 7 hub genes with auto-detected hubs\n")
  cat("  - Re-run TGS/TMS, downstream analyses\n")
} else {
  cat("NO translation module found.\n")
  cat("\nROUTE 2: Switch to pathway-level mirror perturbation story.\n")
  cat("  - Remove 'conserved co-expression program' claims\n")
  cat("  - ssGSEA results become the core finding\n")
  cat("  - WGCNA provides disease-specific network context only\n")
  cat("  - TGS → TMS or remove altogether\n")
}

cat("\n============================================\n")
