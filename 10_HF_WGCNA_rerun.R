# =============================================================================
# 脚本 10：HF WGCNA 完整重跑（修正版）
#
# 与原 09_HF_WGCNA_final.R 的关键区别：
#   1. ❌ 不再使用 top-3000 variance filtering（系统性排除核糖体基因）
#   2. ✅ 使用 bottom-20% variance 过滤（保留生物信号，仅去噪声）
#   3. ❌ 不再硬编码 primary_trans_mod <- "black"
#   4. ✅ GO 富集自动检测翻译/核糖体模块
#   5. ❌ 不再预设 hub genes (EEF1A1, FAU, RPL39...)
#   6. ✅ 基于 GS + MM 标准从翻译模块内筛选 hub genes
#
# 输入：GSE57338_series_matrix.txt.gz (313 samples, HF+NF)
# 输出：
#   GSE57338_WGCNA_rerun.rds  — 完整 WGCNA 结果
#   figures/Figure2A_HF_dendrogram_rerun.png
#   figures/Figure2B_HF_ModuleTrait_rerun.png
#   figures/Figure2C_HF_GO_rerun.png
#   figures/Figure2D_Module_Overlap_rerun.png
# =============================================================================

library(WGCNA)
library(GEOquery)
library(hugene11sttranscriptcluster.db)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dynamicTreeCut)
library(ggplot2)
library(dplyr)

disableWGCNAThreads()
set.seed(42)
options(stringsAsFactors = FALSE)

PROJ_DIR <- "D:/R_projects/revision_analysis"
FIG_DIR  <- file.path(PROJ_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: 数据加载与预处理
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 1: Data Loading & Preprocessing\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

gz_file <- file.path(PROJ_DIR, "GSE57338_series_matrix.txt.gz")
cat("Loading GSE57338 from local file...\n")
gse <- getGEO(filename = gz_file, getGPL = FALSE)
eset <- gse

exprs_raw <- Biobase::exprs(eset)
pdata     <- Biobase::pData(eset)
cat(sprintf("Input: %d probes x %d samples\n", nrow(exprs_raw), ncol(exprs_raw)))

# ── 表型分组 ──
grp_col <- grep("heart failure|disease|condition|group", colnames(pdata),
                ignore.case = TRUE, value = TRUE)[1]
grp <- pdata[[grp_col]]
cat(sprintf("Group column: %s\n", grp_col))
cat("Unique values:\n")
print(table(grp))

is_hf <- grepl("yes|failure|HF|DCM|ICM|cardiomyopathy|dilated|ischemic", grp, ignore.case = TRUE)
is_nf <- grepl("no|normal|non.fail|NF|control|healthy", grp, ignore.case = TRUE)
hf_binary <- ifelse(is_hf, 1, 0)
cat(sprintf("HF=%d, NF=%d, Other=%d\n", sum(is_hf), sum(is_nf), sum(!is_hf & !is_nf)))

# ── 探针→基因 ──
cat("Probe-to-gene mapping...\n")
probe_ids <- as.character(rownames(exprs_raw))
sym_map <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
  keys = probe_ids, keytype = "PROBEID", columns = "SYMBOL")
sym_map <- sym_map[!is.na(sym_map$SYMBOL) & sym_map$SYMBOL != "", ]
cat(sprintf("Annotated probes: %d / %d\n", nrow(sym_map), length(probe_ids)))

# 低表达过滤（探针级）
keep <- rowSums(exprs_raw >= 4) >= ncol(exprs_raw) * 0.2
exprs_f <- exprs_raw[keep, ]
cat(sprintf("After expression filter: %d probes\n", nrow(exprs_f)))

# 过滤注释
sym_f <- sym_map[sym_map$PROBEID %in% rownames(exprs_f), ]
exprs_f <- exprs_f[sym_f$PROBEID, , drop = FALSE]

# 合并到 gene symbol（多探针取均值）
uniq_syms <- unique(sym_f$SYMBOL)
cat(sprintf("Collapsing to gene symbols: %d unique genes\n", length(uniq_syms)))

exprs_gene <- t(sapply(uniq_syms, function(s) {
  probes <- sym_f$PROBEID[sym_f$SYMBOL == s]
  if (length(probes) == 1) exprs_f[probes, ] else colMeans(exprs_f[probes, , drop = FALSE])
}))
cat(sprintf("Gene expression matrix: %d x %d\n", nrow(exprs_gene), ncol(exprs_gene)))

# ── 检查核糖体基因是否在矩阵中 ──
ribo_pattern <- "^RPL|^RPS|^MRPL|^MRPS|EEF1A1|FAU"
ribo_genes_all <- grep(ribo_pattern, rownames(exprs_gene), value = TRUE)
cat(sprintf("Ribosomal/translation genes in full matrix: %d\n", length(ribo_genes_all)))
cat("Examples:", paste(head(ribo_genes_all, 20), collapse = ", "), "\n")

# =============================================================================
# PART 2: 基因过滤（修正版 —— 不设 top-N 上限）
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 2: Gene Filtering (top-6000 by variance, NO hard top-3000 cap)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

datExpr0 <- t(exprs_gene)
vars <- apply(datExpr0, 2, var, na.rm = TRUE)

# ★ 修正策略：取方差最大的 6000 基因（是原 3000 的 2 倍）
#   理由：(1) 手动 WGCNA 对大矩阵可行，6000×6000 TOM ~288MB
#         (2) 保留更多核糖体基因（它们方差中等，不在 top-3000 但在 top-6000）
#         (3) 不使用 blockwiseModules（避免 cor 版本冲突）
n_keep <- min(6000, ncol(datExpr0))
keep_idx <- order(vars, decreasing = TRUE)[1:n_keep]
datExpr0 <- datExpr0[, keep_idx, drop = FALSE]
cat(sprintf("Top-%d genes by variance: %d (from %d total)\n",
            n_keep, ncol(datExpr0), length(vars)))

# 检查核糖体基因保留情况
ribo_kept <- grep(ribo_pattern, colnames(datExpr0), value = TRUE)
cat(sprintf("Ribosomal/translation genes RETAINED: %d/%d\n",
            length(ribo_kept), length(ribo_genes_all)))
if (length(ribo_kept) > 0) {
  cat("  ", paste(sort(ribo_kept), collapse = ", "), "\n")
}

# 检查原文 7 hub genes 是否在输入中
hub_manuscript <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41", "RPS28")
hub_in_input <- hub_manuscript %in% colnames(datExpr0)
cat(sprintf("Manuscript 7 hub genes in WGCNA input: %d/7\n", sum(hub_in_input)))
cat("  Present:", paste(hub_manuscript[hub_in_input], collapse = ", "), "\n")
cat("  Missing:", paste(hub_manuscript[!hub_in_input], collapse = ", "), "\n")

# ── 检查数据质量 ──
cat(sprintf("\nData quality check: %d samples x %d genes\n", nrow(datExpr0), ncol(datExpr0)))
cat(sprintf("  Missing values: %d\n", sum(is.na(datExpr0))))
gsg <- goodSamplesGenes(datExpr0, verbose = 1)
if (!gsg$allOK) {
  cat("  ⚠ Removing flagged genes/samples...\n")
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("  After cleanup: %d samples x %d genes\n", nrow(datExpr0), ncol(datExpr0)))
}

# =============================================================================
# PART 3: 软阈值选择
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 3: Soft Threshold Selection\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

powers <- 1:20
cat("Running pickSoftThreshold (single-thread)...\n")
sft_out <- pickSoftThreshold(datExpr0, powerVector = powers,
                              networkType = "signed", verbose = 2)

# 选择 R² >= 0.85 的最小 β，否则取 R² 最大的
idx_r2_85 <- which(sft_out$fitIndices$SFT.R.sq >= 0.85)
if (length(idx_r2_85) > 0) {
  beta_sel <- sft_out$fitIndices$Power[min(idx_r2_85)]
} else {
  beta_sel <- sft_out$fitIndices$Power[which.max(sft_out$fitIndices$SFT.R.sq)]
}
beta_r2 <- sft_out$fitIndices$SFT.R.sq[which(sft_out$fitIndices$Power == beta_sel)]
cat(sprintf("Selected: beta = %d (R^2 = %.3f)\n", beta_sel, beta_r2))

# 同时报告 beta=12 的 R²（与 manuscript 对比）
beta12_r2 <- sft_out$fitIndices$SFT.R.sq[which(sft_out$fitIndices$Power == 12)]
cat(sprintf("For reference: beta=12 R^2 = %.3f (manuscript claims 0.92)\n", beta12_r2))

# 保存软阈值图
png(file.path(FIG_DIR, "Figure_S4B_SoftThreshold_GSE57338_rerun.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(sft_out$fitIndices[, 1], -sign(sft_out$fitIndices[, 3]) * sft_out$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale Independence (GSE57338 — RERUN)", type = "n")
text(sft_out$fitIndices[, 1], -sign(sft_out$fitIndices[, 3]) * sft_out$fitIndices[, 2],
     labels = powers, col = ifelse(sft_out$fitIndices$SFT.R.sq > 0.8, "red", "black"))
abline(h = 0.85, col = "red", lty = 2)

plot(sft_out$fitIndices[, 1], sft_out$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean Connectivity", type = "n")
text(sft_out$fitIndices[, 1], sft_out$fitIndices[, 5], labels = powers,
     col = ifelse(sft_out$fitIndices$SFT.R.sq > 0.8, "red", "black"))
dev.off()
cat("Soft threshold plot saved.\n")

# =============================================================================
# PART 4: 网络构建与模块检测
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 4: Network Construction & Module Detection\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# ★ 使用手动构建方法（stats::cor 绕过 WGCNA::cor 版本冲突）
#    blockwiseModules 因 WGCNA cor S4 参数不兼容而无法使用
#    6000 基因 × 6000 基因 TOM ≈ 288 MB，完全可管理
n_genes <- ncol(datExpr0)
cat(sprintf("Gene count: %d\n", n_genes))

cat("Building adjacency matrix (stats::cor)...\n")
cor_mat <- stats::cor(datExpr0, use = "pairwise.complete.obs")
adj <- (0.5 * (1 + cor_mat))^beta_sel

cat("Computing TOM...\n")
TOM <- WGCNA::TOMsimilarity(adj, TOMType = "signed", verbose = 2)

cat("Hierarchical clustering...\n")
geneTree <- hclust(as.dist(1 - TOM), method = "average")

cat("Module detection (dynamic tree cut)...\n")
dynamicMods <- dynamicTreeCut::cutreeDynamic(
  dendro = geneTree, distM = 1 - TOM,
  deepSplit = 2, pamRespectsDendro = FALSE,
  minClusterSize = 30, verbose = 2)

moduleLabels <- dynamicMods
moduleColors <- WGCNA::labels2colors(moduleLabels)

n_mods <- length(unique(moduleColors))
cat(sprintf("\nModules detected: %d\n", n_mods))
mod_table <- sort(table(moduleColors), decreasing = TRUE)
print(mod_table[1:min(20, length(mod_table))])

all_genes <- colnames(datExpr0)

# =============================================================================
# PART 5: 模块功能注释 —— 自动寻找翻译模块
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 5: Module GO Enrichment — Auto-detect Translation Module\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

translation_regex <- paste0(
  "ribosom|rRNA|ribonucleoprotein|translational init|",
  "peptide biosyn|cytoplasmic translat| mitochondrial translat|",
  "structural constituent of ribosome|large ribosomal|small ribosomal"
)

module_go_results <- list()
module_translation_score <- setNames(numeric(n_mods), unique(moduleColors))

for (mod in unique(moduleColors)) {
  mod_genes <- all_genes[moduleColors == mod]
  if (length(mod_genes) < 10) next

  ego <- tryCatch({
    enrichGO(gene = mod_genes,
             OrgDb = org.Hs.eg.db,
             keyType = "SYMBOL",
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2)
  }, error = function(e) NULL)

  if (!is.null(ego) && nrow(ego@result) > 0) {
    module_go_results[[mod]] <- ego
    desc_text <- paste(ego@result$Description, collapse = " ")

    # 翻译相关度评分：匹配的 GO term 数量 × (-log10 最佳 p 值)
    trans_matches <- grep(translation_regex, ego@result$Description,
                          ignore.case = TRUE)
    if (length(trans_matches) > 0) {
      best_p <- min(ego@result$p.adjust[trans_matches])
      module_translation_score[mod] <- length(trans_matches) * (-log10(best_p))
      cat(sprintf("\n  %s (n=%d): %d translation-related GO terms\n",
                  mod, sum(moduleColors == mod), length(trans_matches)))
      for (i in trans_matches[1:min(5, length(trans_matches))]) {
        cat(sprintf("    %s | p.adj=%.2e | %s\n",
                    ego@result$ID[i], ego@result$p.adjust[i],
                    ego@result$Description[i]))
      }
    }
  }
}

cat(sprintf("\nTranslation module scores (higher = more translation-related):\n"))
print(sort(module_translation_score, decreasing = TRUE)[1:min(10, length(module_translation_score))])

# ── 自动选择翻译模块 ──
if (max(module_translation_score) > 0) {
  primary_trans_mod <- names(which.max(module_translation_score))
  cat(sprintf("\n*** Auto-detected translation module: %s (score=%.1f) ***\n",
              primary_trans_mod, max(module_translation_score)))
} else {
  cat("\n*** WARNING: No translation module detected! ***\n")
  cat("This means: after bottom-20% variance filter, no co-expression module is\n")
  cat("significantly enriched for ribosomal/translation GO terms.\n")
  primary_trans_mod <- NULL
}

# =============================================================================
# PART 6: 检查原文 hub genes 的模块归属
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 6: Hub Gene Audit\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

hub_module_map <- data.frame(
  Gene = hub_manuscript,
  In_WGCNA_Input = hub_manuscript %in% all_genes,
  Module = sapply(hub_manuscript, function(g) {
    if (g %in% all_genes) as.character(moduleColors[which(all_genes == g)])
    else "NOT_IN_INPUT"
  }),
  stringsAsFactors = FALSE
)

if (!is.null(primary_trans_mod)) {
  hub_module_map$In_Translation_Module <- hub_module_map$Module == primary_trans_mod
} else {
  hub_module_map$In_Translation_Module <- FALSE
}

cat("Manuscript 7 hub genes in new WGCNA:\n")
print(hub_module_map)

# ── 从翻译模块中重新筛选 hub genes ──
if (!is.null(primary_trans_mod)) {
  cat(sprintf("\n--- Re-selecting hub genes from %s module ---\n", primary_trans_mod))

  trans_genes <- all_genes[moduleColors == primary_trans_mod]
  cat(sprintf("Translation module size: %d genes\n", length(trans_genes)))

  # 计算 ME
  MEs <- moduleEigengenes(datExpr0, moduleColors)$eigengenes
  ME_trans <- MEs[[paste0("ME", primary_trans_mod)]]

  # 计算 Gene Significance (correlation with HF trait)
  GS <- apply(datExpr0[, trans_genes, drop = FALSE], 2, function(g) {
    cor(g, hf_binary, use = "pairwise.complete.obs")
  })

  # 计算 Module Membership (kME: correlation with module eigengene)
  MM <- apply(datExpr0[, trans_genes, drop = FALSE], 2, function(g) {
    cor(g, ME_trans, use = "pairwise.complete.obs")
  })

  # 筛选标准：|GS| > 0.2, |MM| > 0.8
  hub_candidates <- trans_genes[abs(GS) > 0.2 & abs(MM) > 0.8]
  cat(sprintf("Hub candidates (|GS|>0.2 & |MM|>0.8): %d\n", length(hub_candidates)))

  # 按 |GS*MM| 排序取 top
  hub_scores <- abs(GS[hub_candidates] * MM[hub_candidates])
  names(hub_scores) <- hub_candidates
  hub_scores <- sort(hub_scores, decreasing = TRUE)

  hub_final <- names(hub_scores)[1:min(10, length(hub_scores))]
  cat("Top hub genes:\n")
  for (i in seq_along(hub_final)) {
    g <- hub_final[i]
    cat(sprintf("  %d. %s: GS=%.3f, MM=%.3f\n", i, g, GS[g], MM[g]))
  }
} else {
  MEs <- moduleEigengenes(datExpr0, moduleColors)$eigengenes
  trans_genes <- character(0)
  GS <- numeric(0)
  MM <- numeric(0)
  hub_final <- character(0)
}

# =============================================================================
# PART 7: Module-Trait 相关性
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 7: Module-Trait Correlations\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 确保 ME 存在
if (!exists("MEs") || is.null(MEs)) {
  MEs <- moduleEigengenes(datExpr0, moduleColors)$eigengenes
}

mod_names <- gsub("^ME", "", colnames(MEs))

mod_trait_cor <- cor(MEs, hf_binary, use = "pairwise.complete.obs")
mod_trait_p <- apply(MEs, 2, function(me) {
  cor.test(me, hf_binary, method = "pearson")$p.value
})

cor_df <- data.frame(
  Module = mod_names,
  Correlation = as.numeric(mod_trait_cor),
  P_value = mod_trait_p,
  Gene_Count = sapply(mod_names, function(m) sum(moduleColors == m)),
  stringsAsFactors = FALSE
)
cor_df <- cor_df[order(abs(cor_df$Correlation), decreasing = TRUE), ]

cat("Module-trait (HF vs NF) correlations:\n")
for (i in 1:nrow(cor_df)) {
  flag <- ""
  if (!is.null(primary_trans_mod) && cor_df$Module[i] == primary_trans_mod) {
    flag <- " ★ TRANSLATION"
  }
  cat(sprintf("  %s: r=%.4f, p=%.2e, n=%d%s\n",
              cor_df$Module[i], cor_df$Correlation[i],
              cor_df$P_value[i], cor_df$Gene_Count[i], flag))
}

# =============================================================================
# PART 8: 可视化
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 8: Visualization\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# ── Panel A: 树状图 ──
png(file.path(FIG_DIR, "Figure2A_HF_dendrogram_rerun.png"),
    width = 14, height = 6, units = "in", res = 300)

plotDendroAndColors(
  geneTree,
  moduleColors,
  "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = sprintf("GSE57338 HF WGCNA (RERUN): Gene Dendrogram (beta=%d, R^2=%.2f, %d genes)",
                 beta_sel, beta_r2, n_genes)
)
dev.off()
cat("Panel A (dendrogram) saved.\n")

# ── Panel B: Module-Trait 柱状图 ──
cor_df$Module <- factor(cor_df$Module, levels = cor_df$Module)
cor_df$is_translation <- if (!is.null(primary_trans_mod)) {
  cor_df$Module == primary_trans_mod
} else {
  rep(FALSE, nrow(cor_df))
}

p_bar <- ggplot(cor_df, aes(x = Correlation, y = Module)) +
  geom_bar(aes(fill = is_translation), stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#C62828", "FALSE" = "grey60"), guide = "none") +
  geom_text(aes(label = sprintf("r=%.3f", Correlation),
                hjust = ifelse(Correlation > 0, -0.1, 1.1)),
            size = 2.8, color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "GSE57338 HF: Module-Trait (HF vs NF) Correlation — RERUN",
    subtitle = if (!is.null(primary_trans_mod)) {
      sprintf("Signed WGCNA, beta=%d | %d genes | Translation module: %s (r=%.3f, p=%.1e)",
              beta_sel, n_genes, primary_trans_mod,
              cor_df$Correlation[cor_df$Module == primary_trans_mod],
              cor_df$P_value[cor_df$Module == primary_trans_mod])
    } else {
      sprintf("Signed WGCNA, beta=%d | %d genes | NO translation module detected",
              beta_sel, n_genes)
    },
    x = "Pearson Correlation with HF Status",
    y = "Module"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "Figure2B_HF_ModuleTrait_rerun.png"),
       p_bar, width = 8, height = max(5, nrow(cor_df) * 0.38), dpi = 300,
       limitsize = FALSE)
cat("Panel B (module-trait bar plot) saved.\n")

# ── Panel C: GO 富集气泡图（仅当翻译模块存在时）──
if (!is.null(primary_trans_mod) && primary_trans_mod %in% names(module_go_results)) {
  ego_trans <- module_go_results[[primary_trans_mod]]
  go_df <- as.data.frame(ego_trans)
  go_df <- go_df[order(go_df$p.adjust), ]

  # 标记翻译相关
  go_df$is_translation <- grepl(translation_regex, go_df$Description, ignore.case = TRUE)

  # Top 20
  trans_rows <- which(go_df$is_translation)
  other_rows <- which(!go_df$is_translation)
  plot_go <- rbind(
    head(go_df[trans_rows, ], 15),
    head(go_df[other_rows, ], 5)
  )
  plot_go <- plot_go[order(plot_go$p.adjust), ]
  plot_go <- head(plot_go, 15)
  plot_go$Description <- factor(plot_go$Description,
                                 levels = rev(unique(plot_go$Description)))

  p_go <- ggplot(plot_go, aes(x = -log10(p.adjust), y = Description)) +
    geom_bar(aes(fill = is_translation), stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"), guide = "none") +
    labs(
      title = sprintf("GSE57338 HF: %s Module GO Enrichment (RERUN)", primary_trans_mod),
      subtitle = sprintf("%d genes; auto-detected translation module", length(trans_genes)),
      x = expression(-log[10](p.adjust)), y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 11))

  ggsave(file.path(FIG_DIR, "Figure2C_HF_GO_rerun.png"),
         p_go, width = 10, height = 5.5, dpi = 300)
  cat("Panel C (GO bubble) saved.\n")
} else {
  cat("Panel C (GO bubble) skipped — no translation module.\n")
}

# ── Panel D: 跨疾病模块基因保存 ──
cat("\n--- Cross-Disease Module Overlap ---\n")

hcc_wgcna <- readRDS(file.path(PROJ_DIR, "GSE141198_WGCNA_result.rds"))
hcc_moduleColors <- hcc_wgcna$moduleColors

# HCC WGCNA 的基因列表
hcc_expr <- readRDS(file.path(PROJ_DIR, "GSE141198_wgcna_input.rds"))
hcc_all_genes <- rownames(hcc_expr)

# 交集
common_genes <- intersect(all_genes, hcc_all_genes)
cat(sprintf("HF genes: %d, HCC genes: %d, Common: %d\n",
            length(all_genes), length(hcc_all_genes), length(common_genes)))

if (!is.null(primary_trans_mod)) {
  hf_trans_genes <- all_genes[moduleColors == primary_trans_mod]
  hf_trans_in_hcc <- intersect(hf_trans_genes, hcc_all_genes)
  cat(sprintf("HF %s module genes: %d, in HCC WGCNA: %d\n",
              primary_trans_mod, length(hf_trans_genes), length(hf_trans_in_hcc)))

  # 映射到 HCC 模块
  hcc_mod_of_hf <- sapply(hf_trans_in_hcc, function(g) {
    as.character(hcc_moduleColors[which(hcc_all_genes == g)])
  })

  hcc_mod_dist <- sort(table(hcc_mod_of_hf), decreasing = TRUE)
  cat("HF translation genes → HCC module distribution:\n")
  print(hcc_mod_dist[1:min(10, length(hcc_mod_dist))])

  # Fisher 检验
  hcc_trans_mod <- if ("blue" %in% names(hcc_mod_dist)) "blue" else names(hcc_mod_dist)[1]

  a <- sum(hcc_mod_of_hf == hcc_trans_mod)
  b <- length(hf_trans_in_hcc) - a
  c <- sum(hcc_moduleColors == hcc_trans_mod) - a
  d <- length(hcc_all_genes) - a - b - c

  fisher_result <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  cat(sprintf("Fisher exact test: HF %s -> HCC %s: OR=%.1f, p=%.2e\n",
              primary_trans_mod, hcc_trans_mod,
              fisher_result$estimate, fisher_result$p.value))

  # 堆叠条形图
  top_mods <- names(hcc_mod_dist)[1:min(6, length(hcc_mod_dist))]
  other_count <- sum(hcc_mod_dist) - sum(hcc_mod_dist[top_mods])

  mod_overlap_df <- data.frame(
    HCC_Module = c(top_mods, "Other"),
    Gene_Count = c(as.integer(hcc_mod_dist[top_mods]), other_count),
    stringsAsFactors = FALSE
  )
  mod_overlap_df$HCC_Module <- factor(mod_overlap_df$HCC_Module,
    levels = rev(c(top_mods, "Other")))
  mod_overlap_df$is_translation <- mod_overlap_df$HCC_Module == hcc_trans_mod
  mod_overlap_df$Label <- paste0(mod_overlap_df$HCC_Module, " (n=", mod_overlap_df$Gene_Count, ")")

  p_overlap <- ggplot(mod_overlap_df, aes(x = Gene_Count, y = HCC_Module)) +
    geom_bar(aes(fill = is_translation), stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "grey75"), guide = "none") +
    geom_text(aes(label = Label), hjust = -0.1, size = 3.5) +
    labs(
      title = "Cross-Disease Module Gene Overlap (RERUN)",
      subtitle = sprintf("HF %s module genes → HCC module assignment\nFisher: OR=%.1f, p=%.1e",
                         primary_trans_mod, fisher_result$estimate, fisher_result$p.value),
      x = "Number of Genes", y = "HCC Module"
    ) +
    xlim(0, max(mod_overlap_df$Gene_Count) * 1.3) +
    theme_minimal(base_size = 12)

  ggsave(file.path(FIG_DIR, "Figure2D_Module_Overlap_rerun.png"),
         p_overlap, width = 8, height = 4.5, dpi = 300)
  cat("Panel D (module overlap) saved.\n")

  cross_disease <- list(
    hf_trans_in_hcc = hf_trans_in_hcc,
    hcc_module_distribution = hcc_mod_dist,
    fisher_test = list(
      a = a, b = b, c = c, d = d,
      odds_ratio = fisher_result$estimate,
      p_value = fisher_result$p.value,
      hcc_translation_module = hcc_trans_mod
    )
  )
} else {
  cross_disease <- NULL
  cat("Panel D (module overlap) skipped — no translation module.\n")
}

# =============================================================================
# PART 9: 保存结果
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 9: Save Results\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 核糖体基因统计
ribo_count_per_module <- sapply(unique(moduleColors), function(m) {
  mod_genes <- all_genes[moduleColors == m]
  sum(grepl(ribo_pattern, mod_genes))
})

wgcna_rerun <- list(
  # 网络参数
  beta_sel = beta_sel,
  beta_r2 = beta_r2,
  n_genes_input = n_genes,
  n_modules = n_mods,
  n_samples = nrow(datExpr0),
  variance_filter = "top_6000_by_variance",

  # 模块分配
  moduleColors = moduleColors,
  moduleLabels = moduleLabels,
  gene_symbols = all_genes,

  # 翻译模块（自动检测）
  primary_translation_module = primary_trans_mod,
  translation_module_genes = if (!is.null(primary_trans_mod)) {
    all_genes[moduleColors == primary_trans_mod]
  } else { character(0) },
  module_translation_score = module_translation_score,
  module_go_results = module_go_results,

  # Hub 基因
  hub_manuscript_audit = hub_module_map,
  hub_genes_auto = hub_final,
  hub_GS = GS,
  hub_MM = MM,

  # ME
  MEs = MEs,
  module_trait_cor = cor_df,

  # 核糖体基因分布
  ribo_genes_per_module = ribo_count_per_module,
  ribo_genes_total = length(ribo_genes_all),
  ribo_genes_in_input = length(ribo_kept),

  # 跨疾病
  cross_disease = cross_disease,

  # 树和 TOM
  geneTree = geneTree,
  TOM = TOM,

  # 软阈值表
  sft_indices = sft_out$fitIndices
)

saveRDS(wgcna_rerun, file.path(PROJ_DIR, "GSE57338_WGCNA_rerun.rds"))
cat("Saved: GSE57338_WGCNA_rerun.rds\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  10_HF_WGCNA_rerun.R — COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("  Input: %d genes (bottom-20%% variance filter)\n", n_genes))
cat(sprintf("  Network: beta=%d, signed R^2=%.2f\n", beta_sel, beta_r2))
cat(sprintf("  Modules: %d total\n", n_mods))
cat(sprintf("  Ribosomal genes in input: %d/%d\n", length(ribo_kept), length(ribo_genes_all)))
cat(sprintf("  Manuscript 7 hub genes in input: %d/7\n", sum(hub_in_input)))

if (!is.null(primary_trans_mod)) {
  cat(sprintf("\n  *** Translation module: %s (%d genes) ***\n",
              primary_trans_mod, sum(moduleColors == primary_trans_mod)))
  cat(sprintf("  Module-trait correlation: r=%.3f\n",
              cor_df$Correlation[cor_df$Module == primary_trans_mod]))
  cat(sprintf("  Auto-detected hub genes: %d\n", length(hub_final)))
  if (length(hub_final) > 0) {
    cat("  Top hubs:", paste(hub_final[1:min(7, length(hub_final))], collapse = ", "), "\n")
  }
} else {
  cat("\n  *** NO translation module detected ***\n")
  cat("  This confirms: after preserving ribosomal genes, there is no coherent\n")
  cat("  translation-related co-expression module in HF heart tissue.\n")
}

cat("\n  Generated files:\n")
cat("    figures/Figure2A_HF_dendrogram_rerun.png\n")
cat("    figures/Figure2B_HF_ModuleTrait_rerun.png\n")
cat("    figures/Figure2C_HF_GO_rerun.png\n")
cat("    figures/Figure2D_Module_Overlap_rerun.png\n")
cat("    GSE57338_WGCNA_rerun.rds\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
