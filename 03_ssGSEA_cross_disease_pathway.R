# =============================================================================
# 脚本 03：ssGSEA跨疾病通路活性效应量比较
# 目标：用单样本通路评分+跨疾病效应量相关分析，替代原有不显著的GSEA重叠检验
# 核心创新点：比较HF和HCC中通路的效应量方向相关性
# =============================================================================

library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(msigdbr)  # MSigDB基因集

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"

# ═══════════════════════════════════════════════════════════════════════════════
# 第一部分：加载MSigDB基因集
# ═══════════════════════════════════════════════════════════════════════════════

message("═══ 加载基因集 ═══")

# Hallmark + KEGG翻译相关 + Reactome翻译通路
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_sets <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)

# KEGG 核糖体
msig_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
kegg_ribosome <- msig_kegg %>%
  filter(gs_name == "KEGG_RIBOSOME") %>%
  pull(gene_symbol) %>%
  unique()

# Reactome 翻译相关通路
msig_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
reactome_translation <- msig_reactome %>%
  filter(grepl("TRANSLATION|PEPTIDE_CHAIN_ELONGATION|EUKARYOTIC_TRANSLATION|RIBOSOME|NONSENSE_MEDIATED_DECAY|SELENOAMINO_ACID|TRNA_AMINOACYLATION|RRNA_PROCESSING",
               gs_name, ignore.case = TRUE))

reactome_gene_sets <- split(reactome_translation$gene_symbol,
                            reactome_translation$gs_name)

# 合并所有基因集
all_gene_sets <- c(hallmark_sets, reactome_gene_sets)
# 添加单独的核糖体通路
all_gene_sets[["KEGG_RIBOSOME"]] <- kegg_ribosome

# 额外添加翻译相关通路
msig_kegg_all <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
kegg_sets <- split(msig_kegg_all$gene_symbol, msig_kegg_all$gs_name)
# 选择翻译相关KEGG通路
translation_kegg <- grep("RIBOSOME|AMINOACYL|PROTEIN_EXPORT|RNA_POLYMERASE|SPLICEOSOME",
                          names(kegg_sets), value = TRUE)
for (nm in translation_kegg) {
  all_gene_sets[[nm]] <- kegg_sets[[nm]]
}

message(sprintf("总基因集数: %d", length(all_gene_sets)))
message(sprintf("Hallmark: %d", length(hallmark_sets)))
message(sprintf("翻译相关通路: KEGG_RIBOSOME (%d genes) + Reactome翻译 (%d sets)",
                length(kegg_ribosome), length(reactome_gene_sets)))

# ═══════════════════════════════════════════════════════════════════════════════
# 第二部分：TCGA-LIHC ssGSEA
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ TCGA-LIHC ssGSEA ═══")

# 加载TCGA-LIHC肿瘤表达数据
se <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")

# 获取表达矩阵（使用log2转换后的counts或TPM）
# 如果可用，使用vst或rlog
if ("vst" %in% assayNames(se)) {
  expr_lihc <- assay(se, "vst")
} else if ("HTSeq - Counts" %in% assayNames(se)) {
  # 用vst快速转换
  dds <- DESeqDataSet(se, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  # 仅对>10样本中>5 reads的基因
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds <- dds[keep, ]
  vsd <- vst(dds, blind = TRUE, nsub = 1000)
  expr_lihc <- assay(vsd)
} else {
  expr_lihc <- assay(se, 1)
}

# 仅保留肿瘤样本
if ("shortLetterCode" %in% colnames(colData(se))) {
  tumor_idx <- which(colData(se)$shortLetterCode == "TP")
  normal_idx <- which(colData(se)$shortLetterCode == "NT")
} else {
  tumor_idx <- grep("-01A|-01B", colnames(se))
  normal_idx <- grep("-11A|-11B", colnames(se))
}

message(sprintf("TCGA-LIHC: %d tumor, %d normal", length(tumor_idx), length(normal_idx)))

# ssGSEA
ssgsea_lihc <- gsva(expr_lihc, all_gene_sets,
                    method = "ssgsea",
                    kcdf = "Gaussian",
                    min.sz = 5,
                    max.sz = 500,
                    verbose = TRUE)

message(sprintf("ssGSEA完成: %d pathways × %d samples", nrow(ssgsea_lihc), ncol(ssgsea_lihc)))

# 计算肿瘤 vs 正常的效应量
# 使用所有样本（肿瘤+正常）
if (length(normal_idx) >= 5) {
  # 有足够正常样本：计算Cohen's d
  calc_cohens_d <- function(ssgsea_mat, tumor_cols, normal_cols) {
    apply(ssgsea_mat, 1, function(x) {
      t_mean <- mean(x[tumor_cols], na.rm = TRUE)
      n_mean <- mean(x[normal_cols], na.rm = TRUE)
      t_sd <- sd(x[tumor_cols], na.rm = TRUE)
      n_sd <- sd(x[normal_cols], na.rm = TRUE)
      pooled_sd <- sqrt((t_sd^2 + n_sd^2) / 2)
      d <- (t_mean - n_mean) / pooled_sd
      return(d)
    })
  }
  li_hc_es <- calc_cohens_d(ssgsea_lihc, tumor_idx, normal_idx)
} else {
  # 无足够正常样本：仅描述肿瘤内通路分数
  li_hc_es <- rowMeans(ssgsea_lihc[, tumor_idx], na.rm = TRUE)
  message("⚠ 正常样本不足(n<5)，使用肿瘤内均值作为描述性指标")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 第三部分：GSE57338 ssGSEA
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ GSE57338 (HF) ssGSEA ═══")

# 加载HF表达数据
# GSE57338是微阵列数据，需要通过GEOquery解析
# 使用已处理好的数据
hf_expr_file <- "D:/R_projects/GSE57338_series_matrix.txt.gz"

if (file.exists(hf_expr_file)) {
  # 读取series matrix
  library(GEOquery)
  gse_hf <- getGEO(filename = hf_expr_file, getGPL = TRUE)
  eset_hf <- gse_hf

  # 寻找HF和正常对照样本
  pd_hf <- pData(eset_hf)

  # 识别表型列（DCM/ICM/NF等）
  # GSE57338包含：NF (non-failing), DCM, ICM
  # 寻找phenotype列
  pheno_cols <- grep("diagnosis|disease|phenotype|group|condition|tissue|status",
                     colnames(pd_hf), ignore.case = TRUE, value = TRUE)

  if (length(pheno_cols) > 0) {
    cat("可能的表型列:", pheno_cols, "\n")
    # 打印唯一值帮助判断
    for (col in pheno_cols[1:min(3, length(pheno_cols))]) {
      cat(sprintf("  %s: %s\n", col, paste(unique(pd_hf[, col])[1:min(4, length(unique(pd_hf[, col])))], collapse=", ")))
    }
  }

  expr_hf <- exprs(eset_hf)
  message(sprintf("HF表达矩阵: %d genes × %d samples", nrow(expr_hf), ncol(expr_hf)))

  # 确保基因名是gene symbol（可能需要从探针ID转换）
  # GPL11532是Affymetrix Human Gene 1.1 ST Array
  # 特征数据中通常包含gene symbol
  fdata <- fData(eset_hf)
  symbol_col <- grep("symbol|gene_symbol|gene_assignment", colnames(fdata),
                     ignore.case = TRUE, value = TRUE)
  if (length(symbol_col) > 0) {
    message("使用fData列作为gene symbol: ", symbol_col[1])
  }

  # ssGSEA
  ssgsea_hf <- gsva(expr_hf, all_gene_sets,
                    method = "ssgsea",
                    kcdf = "Gaussian",
                    min.sz = 5,
                    max.sz = 500,
                    verbose = TRUE)

  # 计算HF vs NF效应量
  # 需要识别HF和NF样本
  # 假设能从pheno中找到
  hf_es <- rowMeans(ssgsea_hf, na.rm = TRUE)  # 占位，需替换

  message(sprintf("HF ssGSEA完成: %d pathways × %d samples", nrow(ssgsea_hf), ncol(ssgsea_hf)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 第四部分：跨疾病通路效应量相关分析（核心分析）
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 跨疾病效应量相关分析 ═══")

# 注：此部分在ssGSEA完成并正确分组后运行
# 以下提供分析框架

analyze_cross_disease_pathway <- function(ssgsea_lihc, ssgsea_hf,
                                          li_hc_groups, hf_groups,
                                          output_dir = NULL) {
  # li_hc_groups: TCGA样本分组 (tumor/normal)
  # hf_groups: HF样本分组 (HF/NF)

  # 1. 计算各通路的效应量
  calc_effect_size <- function(ssgsea, groups_1, groups_2) {
    apply(ssgsea, 1, function(x) {
      m1 <- mean(x[groups_1], na.rm = TRUE)
      m2 <- mean(x[groups_2], na.rm = TRUE)
      s1 <- sd(x[groups_1], na.rm = TRUE)
      s2 <- sd(x[groups_2], na.rm = TRUE)
      pooled_sd <- sqrt((s1^2 + s2^2) / 2)
      d <- (m1 - m2) / pooled_sd
      return(d)
    })
  }

  # 仅用Hallmark通路做主要分析（避免冗余）
  hallmark_names <- names(hallmark_sets)
  common_pathways <- intersect(hallmark_names,
                               intersect(rownames(ssgsea_lihc), rownames(ssgsea_hf)))

  es_lihc <- calc_effect_size(ssgsea_lihc[common_pathways, ], li_hc_groups[[1]], li_hc_groups[[2]])
  es_hf <- calc_effect_size(ssgsea_hf[common_pathways, ], hf_groups[[1]], hf_groups[[2]])

  # 2. Spearman相关
  cor_test <- cor.test(es_lihc, es_hf, method = "spearman")
  message(sprintf("\n跨疾病通路效应量 Spearman ρ = %.3f, p = %.4f",
                  cor_test$estimate, cor_test$p.value))

  # 3. 置换检验（打乱通路标签，1000次）
  n_perm <- 1000
  perm_rhos <- numeric(n_perm)
  for (i in 1:n_perm) {
    perm_rhos[i] <- cor(es_lihc, sample(es_hf), method = "spearman")
  }
  perm_p <- mean(abs(perm_rhos) >= abs(cor_test$estimate))

  message(sprintf("置换检验 (n=%d): 经验p = %.4f", n_perm, perm_p))

  # 4. 识别翻译相关通路
  translation_pathways <- grep("RIBOSOM|TRANSLAT|PEPTIDE|AMINOACYL|TRNA|RRNA|EIF|EEF|MYC_TARGET",
                               common_pathways, value = TRUE, ignore.case = TRUE)

  # 5. 构建结果表
  pathway_results <- data.frame(
    pathway = common_pathways,
    es_hcc = es_lihc,
    es_hf = es_hf,
    is_translation = common_pathways %in% translation_pathways,
    stringsAsFactors = FALSE
  )

  # 6. 绘制效应量散点图
  pathway_results$category <- ifelse(pathway_results$is_translation,
                                     "Translation-related", "Other Hallmark")

  p <- ggplot(pathway_results, aes(x = es_hcc, y = es_hf)) +
    geom_point(aes(color = category), size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "darkgrey", alpha = 0.3) +
    ggrepel::geom_text_repel(
      data = subset(pathway_results, is_translation | abs(es_hcc) > 2 | abs(es_hf) > 2),
      aes(label = gsub("HALLMARK_|KEGG_|REACTOME_", "", pathway)),
      size = 3, max.overlaps = 20
    ) +
    scale_color_manual(values = c("Translation-related" = "#E41A1C",
                                  "Other Hallmark" = "#377EB8")) +
    geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) +
    labs(
      x = "Effect Size in HCC (Tumor vs Normal, Cohen's d)",
      y = "Effect Size in HF (HF vs NF, Cohen's d)",
      title = "Cross-Disease Pathway Activity Comparison",
      subtitle = sprintf("Spearman ρ = %.3f, Permutation p = %.4f (n=%d)",
                        cor_test$estimate, perm_p, n_perm),
      caption = sprintf("N=%d common Hallmark pathways", length(common_pathways))
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())

  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, "Figure_CrossDisease_PathwayES.png"),
           p, width = 10, height = 8, dpi = 300)
  }

  # 7. 翻译相关通路详细表
  trans_pathway_detail <- pathway_results %>%
    filter(is_translation | grepl("MYC|MTOR|E2F", pathway)) %>%
    arrange(desc(abs(es_hcc) + abs(es_hf)))

  print(trans_pathway_detail)

  return(list(
    pathway_results = pathway_results,
    cor_test = cor_test,
    perm_p = perm_p,
    perm_rhos = perm_rhos,
    trans_detail = trans_pathway_detail,
    plot = p
  ))
}

message("\n✅ 脚本03完成：ssGSEA分析框架就绪")
message("后续步骤：")
message("  1. 确认GSE57338的HF/NF分组后，重新运行ssGSEA")
message("  2. 运行analyze_cross_disease_pathway()获取核心结果")
message("  3. 如果跨疾病Spearman ρ显著，这是全文最强的通路层面证据")

saveRDS(list(
  all_gene_sets = all_gene_sets,
  hallmark_sets = hallmark_sets,
  analyze_cross_disease_pathway = analyze_cross_disease_pathway
), file.path(PROJ_DIR, "ssgsea_framework.rds"))
