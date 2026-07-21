# =============================================================================
# 脚本 05b：上游转录因子预测（重写版）
# 目标：
#   1. hub基因在已知TF靶基因集中的富集（Fisher精确检验）
#   2. 核糖体全体基因在MYC/E2F/mTORC1靶基因中的富集
#   3. TCGA-LIHC中TF表达与TGS的相关性
#   4. MYC/E2F/mTORC1通路活性与TGS的相关性比较
# =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(GSVA)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
dir.create(file.path(PROJ_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# Part 1: 准备基因集
# ═══════════════════════════════════════════════════════════════════════════════

message("═══ 准备基因集 ═══")

hub_genes <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41", "RPS28")

# KEGG核糖体全部基因（新版msigdbr API）
kegg_ribo <- msigdbr(species = "Homo sapiens", collection = "C2",
                     subcollection = "CP:KEGG_LEGACY") %>%
  filter(gs_name == "KEGG_RIBOSOME") %>%
  pull(gene_symbol) %>% unique()
message(sprintf("KEGG核糖体基因: %d", length(kegg_ribo)))

# Hallmark通路
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")

myc_v1 <- msig_h %>% filter(gs_name == "HALLMARK_MYC_TARGETS_V1") %>% pull(gene_symbol) %>% unique()
myc_v2 <- msig_h %>% filter(gs_name == "HALLMARK_MYC_TARGETS_V2") %>% pull(gene_symbol) %>% unique()
e2f    <- msig_h %>% filter(gs_name == "HALLMARK_E2F_TARGETS")    %>% pull(gene_symbol) %>% unique()
mtorc1 <- msig_h %>% filter(gs_name == "HALLMARK_MTORC1_SIGNALING") %>% pull(gene_symbol) %>% unique()

message(sprintf("MYC Targets V1: %d genes", length(myc_v1)))
message(sprintf("MYC Targets V2: %d genes", length(myc_v2)))
message(sprintf("E2F Targets:    %d genes", length(e2f)))
message(sprintf("mTORC1:         %d genes", length(mtorc1)))

# ═══════════════════════════════════════════════════════════════════════════════
# Part 2: Hub基因在TF靶基因集中的重叠
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ Hub基因与TF靶基因集重叠 ═══")

# 获取数据集中的所有基因作为背景
se <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")
# Ensembl IDs → Symbol（需要映射后才能比对）
library(AnnotationDbi)
ensembl_ids <- sub("\\.\\d+$", "", rownames(se))
sym_map <- AnnotationDbi::select(org.Hs.eg.db,
  keys = ensembl_ids, keytype = "ENSEMBL", columns = "SYMBOL")
sym_map <- sym_map[!is.na(sym_map$SYMBOL), ]
bg_genes <- unique(sym_map$SYMBOL)
bg_size <- length(bg_genes)
message(sprintf("背景基因集: %d genes", bg_size))

# hub基因存在性
hub_found <- hub_genes[hub_genes %in% bg_genes]
message(sprintf("Hub基因存在于数据集: %d/7: %s", length(hub_found), paste(hub_found, collapse=", ")))

# 重叠检查
overlap_check <- function(hub_set, tf_set, tf_name) {
  overlap <- intersect(hub_set, tf_set)
  message(sprintf("  %s: %d/7 — %s", tf_name, length(overlap), paste(overlap, collapse=", ")))
  return(overlap)
}

hub_myc_v1 <- overlap_check(hub_found, myc_v1, "MYC V1")
hub_myc_v2 <- overlap_check(hub_found, myc_v2, "MYC V2")
hub_e2f    <- overlap_check(hub_found, e2f,    "E2F")
hub_mtorc1 <- overlap_check(hub_found, mtorc1, "mTORC1")

# Fisher精确检验
run_fisher <- function(hub_in_pathway, pathway_genes, bg, pathway_name) {
  a <- length(hub_in_pathway)
  b <- length(hub_found) - a
  c <- length(intersect(bg, pathway_genes)) - a
  d <- bg_size - a - b - c

  if (c < 0 || d < 0) {
    message(sprintf("  %s: Fisher参数无效", pathway_name))
    return(NULL)
  }

  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  message(sprintf("  %s: OR=%.2f, p=%.4f (%d/%d hub genes)",
                  pathway_name, ft$estimate, ft$p.value, a, length(hub_found)))
  return(ft)
}

message("\n═══ Fisher精确检验: Hub基因富集 ═══")
fisher_myc_v1 <- run_fisher(hub_myc_v1, myc_v1, bg_genes, "MYC V1")
fisher_myc_v2 <- run_fisher(hub_myc_v2, myc_v2, bg_genes, "MYC V2")
fisher_e2f    <- run_fisher(hub_e2f,    e2f,    bg_genes, "E2F")
fisher_mtorc1 <- run_fisher(hub_mtorc1, mtorc1, bg_genes, "mTORC1")

# ═══════════════════════════════════════════════════════════════════════════════
# Part 3: 全体核糖体基因在TF靶基因中的富集
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 全体核糖体基因富集 ═══")

ribo_in_bg <- intersect(kegg_ribo, bg_genes)
message(sprintf("核糖体基因在数据集中: %d/%d", length(ribo_in_bg), length(kegg_ribo)))

if (length(ribo_in_bg) >= 30) {
  ribo_in_myc1 <- intersect(ribo_in_bg, myc_v1)
  ribo_in_myc2 <- intersect(ribo_in_bg, myc_v2)
  ribo_in_e2f  <- intersect(ribo_in_bg, e2f)
  ribo_in_mtor <- intersect(ribo_in_bg, mtorc1)

  run_fisher_kegg <- function(kegg_set, tf_set, tf_name) {
    a <- length(kegg_set)
    b <- length(ribo_in_bg) - a
    c <- length(intersect(bg_genes, tf_set)) - a
    d <- bg_size - a - b - c

    if (c < 0 || d < 0) {
      message(sprintf("  %s: Fisher参数无效", tf_name))
      return(NULL)
    }

    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    message(sprintf("  %s: OR=%.2f, p=%.2e (%d/%d ribo genes in pathway)",
                    tf_name, ft$estimate, ft$p.value, a, length(ribo_in_bg)))
    return(ft)
  }

  f_ribo_myc1 <- run_fisher_kegg(ribo_in_myc1, myc_v1, "MYC V1")
  f_ribo_myc2 <- run_fisher_kegg(ribo_in_myc2, myc_v2, "MYC V2")
  f_ribo_e2f  <- run_fisher_kegg(ribo_in_e2f,  e2f,    "E2F")
  f_ribo_mtor <- run_fisher_kegg(ribo_in_mtor, mtorc1, "mTORC1")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Part 4: TCGA-LIHC中TF表达与TGS的相关性
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ TCGA-LIHC: TF表达 vs TGS ═══")

# 加载TCGA-LIHC VST表达（复用脚本03b的数据准备逻辑）
se <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")
counts_lihc <- assay(se, "unstranded")
n_samples <- ncol(counts_lihc)

# 仅肿瘤样本
tumor_cols <- grep("-01A|-01B", colnames(counts_lihc))
counts_tumor <- counts_lihc[, tumor_cols]

dds <- DESeqDataSetFromMatrix(
  countData = counts_tumor,
  colData = S4Vectors::DataFrame(row.names = colnames(counts_tumor)),
  design = ~ 1)
keep <- rowSums(counts(dds) >= 10) >= floor(0.2 * ncol(dds))
dds <- dds[keep, ]
vsd <- vst(dds, blind = TRUE, nsub = 1000)
expr_tumor_vst <- assay(vsd)

# Ensembl → Symbol
ens_ids <- sub("\\.\\d+$", "", rownames(expr_tumor_vst))
s_map <- AnnotationDbi::select(org.Hs.eg.db,
  keys = ens_ids, keytype = "ENSEMBL", columns = "SYMBOL")
s_map <- s_map[!is.na(s_map$SYMBOL) & !duplicated(s_map$ENSEMBL), ]
keep_rows <- match(s_map$ENSEMBL, ens_ids)
expr_tumor_vst <- expr_tumor_vst[keep_rows, ]
rownames(expr_tumor_vst) <- s_map$SYMBOL
expr_tumor_vst <- expr_tumor_vst[!duplicated(rownames(expr_tumor_vst)), ]
message(sprintf("TCGA-LIHC tumor: %d genes × %d samples", nrow(expr_tumor_vst), ncol(expr_tumor_vst)))

# 计算TGS
hub_6 <- intersect(c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41"), rownames(expr_tumor_vst))
expr_hub <- expr_tumor_vst[hub_6, , drop = FALSE]
expr_z <- t(scale(t(expr_hub)))
tgs_score <- colMeans(expr_z, na.rm = TRUE)

# 候选TF列表
candidate_tfs <- c("MYC", "MYCN", "MYCL", "E2F1", "E2F2", "E2F3",
                   "E2F4", "E2F5", "E2F6", "E2F7", "E2F8",
                   "MTOR", "RPTOR", "TP53", "HIF1A", "NFE2L2",
                   "ATF4", "XBP1", "DDIT3")

# TF-TGS相关性
tf_cor_results <- data.frame()
for (tf in candidate_tfs) {
  if (!tf %in% rownames(expr_tumor_vst)) next

  tf_expr <- as.numeric(expr_tumor_vst[tf, ])
  ct <- cor.test(tf_expr, tgs_score, method = "spearman")

  tf_cor_results <- rbind(tf_cor_results, data.frame(
    TF = tf,
    rho = round(ct$estimate, 4),
    p_value = ct$p.value,
    p_adj = NA_real_,
    stringsAsFactors = FALSE
  ))
}

tf_cor_results$p_adj <- p.adjust(tf_cor_results$p_value, method = "BH")
tf_cor_results <- tf_cor_results[order(abs(tf_cor_results$rho), decreasing = TRUE), ]

cat("\n═══ TF表达与TGS相关性 ═══\n")
print(tf_cor_results, row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# Part 5: MYC/E2F/mTORC1通路活性与TGS的相关性比较（ssGSEA通路活性）
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 关键通路活性 vs TGS ═══")

# 加载ssGSEA结果（从脚本03b）
ssgsea_res <- readRDS(file.path(PROJ_DIR, "ssgsea_cross_disease_result.rds"))
ssgsea_lihc <- ssgsea_res$ssgsea_lihc

# 仅肿瘤样本的ssGSEA
tumor_cols_ssgsea <- grep("-01A|-01B", colnames(ssgsea_lihc))
ssgsea_tumor <- ssgsea_lihc[, tumor_cols_ssgsea, drop = FALSE]

# 对比上游通路
upstream_pw <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
                 "HALLMARK_E2F_TARGETS", "HALLMARK_MTORC1_SIGNALING")
existing_pw <- intersect(upstream_pw, rownames(ssgsea_tumor))
message(sprintf("可用上游通路: %d/%d", length(existing_pw), length(upstream_pw)))

# 匹配样本
common_samps <- intersect(colnames(ssgsea_tumor), names(tgs_score))
message(sprintf("共同样本: %d", length(common_samps)))

pw_cor_results <- data.frame()
for (pw in existing_pw) {
  pw_score <- as.numeric(ssgsea_tumor[pw, common_samps])
  ct <- cor.test(tgs_score[common_samps], pw_score, method = "spearman")

  pw_cor_results <- rbind(pw_cor_results, data.frame(
    pathway = gsub("HALLMARK_", "", pw),
    rho = round(ct$estimate, 4),
    p_value = ct$p.value,
    stringsAsFactors = FALSE
  ))
}

pw_cor_results <- pw_cor_results[order(abs(pw_cor_results$rho), decreasing = TRUE), ]

cat("\n═══ 上游通路活性与TGS相关性 ═══\n")
print(pw_cor_results, row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# Part 6: 可视化
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 可视化 ═══")

# 6.1 TF-TGS相关性条形图
if (nrow(tf_cor_results) > 0) {
  plot_tf <- tf_cor_results[1:min(15, nrow(tf_cor_results)), ]
  plot_tf$TF <- factor(plot_tf$TF, levels = rev(plot_tf$TF))
  plot_tf$sig <- ifelse(plot_tf$p_adj < 0.05, "FDR<0.05",
                  ifelse(plot_tf$p_value < 0.05, "p<0.05", "NS"))

  p_tf <- ggplot(plot_tf, aes(x = TF, y = rho, fill = sig)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("FDR<0.05" = "#E41A1C", "p<0.05" = "#FF7F00", "NS" = "#999999")) +
    coord_flip() +
    labs(title = "TF Expression vs Translation Gene Score (TGS) in TCGA-LIHC",
         subtitle = paste0("Spearman correlation, n=", length(tgs_score), " tumors"),
         x = "", y = expression(rho)) +
    theme_minimal(base_size = 12)

  ggsave(file.path(PROJ_DIR, "figures", "Figure_TF_TGS_correlation.png"),
         p_tf, width = 9, height = 5, dpi = 300)
  message("TF-TGS相关图已保存")
}

# 6.2 上游通路 vs TGS
if (nrow(pw_cor_results) > 0) {
  plot_pw <- pw_cor_results
  plot_pw$pathway <- factor(plot_pw$pathway, levels = rev(plot_pw$pathway))
  plot_pw$sig <- ifelse(plot_pw$p_value < 0.05, "p<0.05", "NS")

  p_pw <- ggplot(plot_pw, aes(x = pathway, y = rho, fill = sig)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("p<0.05" = "#E41A1C", "NS" = "#999999")) +
    coord_flip() +
    labs(title = "Pathway Activity vs TGS in TCGA-LIHC",
         subtitle = paste0("ssGSEA pathway scores, n=", length(common_samps), " tumors"),
         x = "", y = expression(rho)) +
    theme_minimal(base_size = 12)

  ggsave(file.path(PROJ_DIR, "figures", "Figure_Pathway_TGS_correlation.png"),
         p_pw, width = 8, height = 4, dpi = 300)
  message("通路-TGS相关图已保存")
}

# 6.3 Fisher检验结果汇总图
fisher_results <- data.frame(
  level = rep(c("7 Hub Genes", "All Ribosome (KEGG)"), each = 4),
  pathway = rep(c("MYC V1", "MYC V2", "E2F", "mTORC1"), 2),
  OR = NA_real_, p = NA_real_, stringsAsFactors = FALSE
)

# 填充（如果Fisher检验成功）
for (i in 1:nrow(fisher_results)) {
  if (fisher_results$level[i] == "7 Hub Genes") {
    res <- switch(fisher_results$pathway[i],
                  "MYC V1" = fisher_myc_v1, "MYC V2" = fisher_myc_v2,
                  "E2F" = fisher_e2f, "mTORC1" = fisher_mtorc1)
  } else {
    res <- switch(fisher_results$pathway[i],
                  "MYC V1" = if (exists("f_ribo_myc1")) f_ribo_myc1 else NULL,
                  "MYC V2" = if (exists("f_ribo_myc2")) f_ribo_myc2 else NULL,
                  "E2F"    = if (exists("f_ribo_e2f"))  f_ribo_e2f  else NULL,
                  "mTORC1" = if (exists("f_ribo_mtor")) f_ribo_mtor else NULL)
  }
  if (!is.null(res)) {
    fisher_results$OR[i] <- res$estimate
    fisher_results$p[i] <- res$p.value
  }
}

fisher_results <- fisher_results[!is.na(fisher_results$OR), ]
if (nrow(fisher_results) > 0) {
  fisher_results$log10p <- -log10(fisher_results$p)
  fisher_results$pathway <- factor(fisher_results$pathway,
    levels = c("MYC V1", "MYC V2", "E2F", "mTORC1"))

  # 将 Hub 基因 p=1.00 的无富集结果加一个最小可视高度（避免柱子为零）
  hub_mask <- fisher_results$level == "7 Hub Genes"
  fisher_results$bar_height <- fisher_results$log10p
  fisher_results$bar_height[hub_mask] <- pmax(fisher_results$log10p[hub_mask], 0.15)

  # 为 Hub 基因准备独立标签（p=1.00 时 -log10=0，标签放柱子顶部）
  fisher_results$label_y <- fisher_results$bar_height + 0.5
  fisher_results$label <- ifelse(hub_mask & fisher_results$p > 0.99,
    sprintf("p=1.00\n(no enrichment)"),
    sprintf("OR=%.1f\np=%.1e", fisher_results$OR, fisher_results$p))

  p_fisher <- ggplot(fisher_results,
    aes(x = pathway, y = bar_height, fill = level)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    annotate("text", x = 4.5, y = -log10(0.05) + 0.3,
             label = expression(p==0.05), size = 3, color = "grey50", hjust = 1) +
    geom_text(aes(label = label, y = label_y,
                  vjust = 0), position = position_dodge(width = 0.8),
              size = 3, lineheight = 0.9) +
    scale_fill_manual(values = c("7 Hub Genes" = "#E41A1C",
                                  "All Ribosome (KEGG)" = "#377EB8")) +
    labs(title = "Fisher Enrichment: Translation Genes in TF Target Sets",
         subtitle = "Hub genes show no enrichment in canonical TF target sets (p = 1.00)",
         x = "TF Target Gene Set", y = expression(-log[10](p)),
         fill = "") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(PROJ_DIR, "figures", "Figure_TF_Fisher_enrichment.png"),
         p_fisher, width = 8, height = 6, dpi = 300)
  message("Fisher富集图已保存")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Part 7: 保存 + 摘要
# ═══════════════════════════════════════════════════════════════════════════════

tf_result <- list(
  hub_genes = hub_genes,
  fisher_hub = list(MYC_V1 = fisher_myc_v1, MYC_V2 = fisher_myc_v2,
                    E2F = fisher_e2f, mTORC1 = fisher_mtorc1),
  fisher_ribosome = list(MYC_V1 = f_ribo_myc1, MYC_V2 = f_ribo_myc2,
                         E2F = f_ribo_e2f, mTORC1 = f_ribo_mtor),
  tf_tgs_cor = tf_cor_results,
  pathway_tgs_cor = pw_cor_results,
  ribo_genes_in_bg = ribo_in_bg
)
saveRDS(tf_result, file.path(PROJ_DIR, "TF_upstream_result.rds"))

cat("\n═══════════════════════════════════════════════\n")
cat("  上游TF预测分析完成\n")
cat("─────────────────────────────────\n")

# 总结关键发现
cat("\n【关键发现】\n")

# 最强的TF
if (nrow(tf_cor_results) > 0) {
  sig_tfs <- tf_cor_results[tf_cor_results$p_adj < 0.05, ]
  nominal_tfs <- tf_cor_results[tf_cor_results$p_value < 0.05 & tf_cor_results$p_adj >= 0.05, ]

  if (nrow(sig_tfs) > 0) {
    cat(sprintf("  FDR显著TF: %d (%s)\n", nrow(sig_tfs),
                paste(sig_tfs$TF, collapse = ", ")))
  }
  if (nrow(nominal_tfs) > 0) {
    cat(sprintf("  名义显著TF: %d (%s)\n", nrow(nominal_tfs),
                paste(nominal_tfs$TF, collapse = ", ")))
  }
}

# Fisher检验结论
for (nm in names(tf_result$fisher_hub)) {
  f <- tf_result$fisher_hub[[nm]]
  if (!is.null(f) && f$p.value < 0.05) {
    cat(sprintf("  ✅ Hub genes significantly enriched in %s (OR=%.1f, p=%.4f)\n",
                nm, f$estimate, f$p.value))
  } else if (!is.null(f)) {
    cat(sprintf("  ❌ Hub genes NOT enriched in %s (OR=%.1f, p=%.4f)\n",
                nm, f$estimate, f$p.value))
  }
}

# 上游通路匹配
cat("\n【上游通路活性 vs TGS】\n")
for (i in 1:nrow(pw_cor_results)) {
  cat(sprintf("  %s: ρ=%.3f, p=%.4f\n",
              pw_cor_results$pathway[i], pw_cor_results$rho[i], pw_cor_results$p_value[i]))
}

cat("═══════════════════════════════════════════════\n")
