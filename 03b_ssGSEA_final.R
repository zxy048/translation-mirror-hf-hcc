# =============================================================================
# 脚本 03b：ssGSEA跨疾病通路活性比较（使用实际数据 + 新版msigdbr API）
# =============================================================================

library(GSVA)
library(msigdbr)
library(limma)
library(ggplot2)
library(dplyr)
library(GEOquery)
library(DESeq2)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
dir.create(file.path(PROJ_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════
# Part 1: 准备基因集（新版msigdbr API：collection="H"大写）
# ══════════════════════════════════════════════════════════════════════
message("═══ 准备基因集 ═══")

msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(msig_h$gene_symbol, msig_h$gs_name)
message(sprintf("Hallmark基因集: %d", length(hallmark_sets)))

# KEGG核糖体
kegg_df <- msigdbr(species = "Homo sapiens", collection = "C2",
                   subcollection = "CP:KEGG_MEDICUS")
kegg_ribo <- kegg_df %>%
  filter(grepl("RIBOSOME", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>% unique()
hallmark_sets[["KEGG_RIBOSOME"]] <- kegg_ribo

# Reactome翻译通路
reac_df <- msigdbr(species = "Homo sapiens", collection = "C2",
                   subcollection = "CP:REACTOME")
reac_trans <- reac_df %>%
  filter(grepl("TRANSLATION|PEPTIDE_CHAIN_ELONGATION|EUKARYOTIC_TRANSLATION|RIBOSOME|NONSENSE_MEDIATED|TRNA_AMINOACYLATION|RRNA_PROCESSING",
               gs_name, ignore.case = TRUE))
for (nm in unique(reac_trans$gs_name)) {
  hallmark_sets[[nm]] <- reac_trans$gene_symbol[reac_trans$gs_name == nm]
}
message(sprintf("总基因集: %d (Hallmark + KEGG Ribosome + Reactome Translation)",
                length(hallmark_sets)))

# ══════════════════════════════════════════════════════════════════════
# Part 2: TCGA-LIHC ssGSEA
# ══════════════════════════════════════════════════════════════════════
message("\n═══ TCGA-LIHC ssGSEA ═══")

library(org.Hs.eg.db)
library(AnnotationDbi)

se <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")

# 使用unstranded assay（原始counts）
counts_lihc <- assay(se, "unstranded")
message(sprintf("Raw counts: %d genes × %d samples", nrow(counts_lihc), ncol(counts_lihc)))

# DESeq2 VST标准化
n_samples <- ncol(counts_lihc)
col_data <- S4Vectors::DataFrame(row.names = colnames(counts_lihc))
dds <- DESeqDataSetFromMatrix(
  countData = counts_lihc,
  colData = col_data,
  design = ~ 1)
keep <- rowSums(counts(dds) >= 10) >= floor(0.2 * ncol(dds))
dds <- dds[keep, ]
message(sprintf("Filtered: %d genes (min 10 counts in >=20%% samples)", nrow(dds)))
vsd <- vst(dds, blind = TRUE, nsub = 1000)
expr_vst <- assay(vsd)

# Ensembl ID → Gene Symbol 转换（msigdbr基因集用Symbol）
ensembl_ids <- sub("\\.\\d+$", "", rownames(expr_vst))  # 去除版本号
symbol_map <- AnnotationDbi::select(org.Hs.eg.db,
  keys = ensembl_ids, keytype = "ENSEMBL", columns = "SYMBOL")
symbol_map <- symbol_map[!is.na(symbol_map$SYMBOL) & !duplicated(symbol_map$ENSEMBL), ]
# 去重Symbol（保留表达量最高的）
keep_rows <- match(symbol_map$ENSEMBL, ensembl_ids)
expr_vst <- expr_vst[keep_rows, ]
rownames(expr_vst) <- symbol_map$SYMBOL
expr_vst <- expr_vst[!duplicated(rownames(expr_vst)), ]
message(sprintf("Ensembl→Symbol: %d genes retained", nrow(expr_vst)))

expr_lihc <- expr_vst

tumor_idx <- grep("-01A|-01B", colnames(expr_lihc))
normal_idx <- grep("-11A|-11B", colnames(expr_lihc))
message(sprintf("TCGA-LIHC: %d tumor, %d normal", length(tumor_idx), length(normal_idx)))

param_lihc <- ssgseaParam(exprData = expr_lihc, geneSets = hallmark_sets,
                          minSize = 5, maxSize = 500)
ssgsea_lihc <- gsva(param_lihc, verbose = FALSE)
message(sprintf("ssGSEA完成: %d pathways × %d samples", nrow(ssgsea_lihc), ncol(ssgsea_lihc)))

lihc_es <- apply(ssgsea_lihc, 1, function(x) {
  t_m <- mean(x[tumor_idx], na.rm = TRUE)
  n_m <- mean(x[normal_idx], na.rm = TRUE)
  ps <- sqrt((sd(x[tumor_idx])^2 + sd(x[normal_idx])^2) / 2)
  if (ps == 0) 0 else (t_m - n_m) / ps
})

# ══════════════════════════════════════════════════════════════════════
# Part 3: GSE57338 (HF) ssGSEA
# ══════════════════════════════════════════════════════════════════════
message("\n═══ GSE57338 (HF) ssGSEA ═══")

gse_hf <- getGEO(filename = "D:/R_projects/GSE57338_series_matrix.txt.gz", getGPL = FALSE)
expr_hf <- exprs(gse_hf)
pd_hf <- pData(gse_hf)

# 探针→Symbol映射：使用Bioconductor注释包
# GPL11532 = Affymetrix Human Gene 1.1 ST Array
# 探针ID是transcript cluster ID（纯数字，如7892501）
suppressPackageStartupMessages({
  library(hugene11sttranscriptcluster.db)
})

probe_ids <- rownames(expr_hf)
message(sprintf("总探针: %d", length(probe_ids)))

# 获取探针→Symbol映射
sym_list <- tryCatch({
  AnnotationDbi::select(hugene11sttranscriptcluster.db,
    keys = probe_ids, keytype = "PROBEID", columns = "SYMBOL")
}, error = function(e) {
  message("select with PROBEID failed, trying keys as is")
  # 探针可能是纯数字，需要特殊处理
  AnnotationDbi::select(hugene11sttranscriptcluster.db,
    keys = as.character(probe_ids), keytype = "PROBEID", columns = "SYMBOL")
})

sym_list <- sym_list[!is.na(sym_list$SYMBOL) & sym_list$SYMBOL != "", ]
sym_list <- sym_list[!duplicated(sym_list$PROBEID), ]
message(sprintf("注释探针: %d/%d", nrow(sym_list), length(probe_ids)))

# 映射
common_probes <- intersect(rownames(expr_hf), sym_list$PROBEID)
expr_hf_mapped <- expr_hf[common_probes, , drop = FALSE]
rownames(expr_hf_mapped) <- sym_list$SYMBOL[match(common_probes, sym_list$PROBEID)]

# 按Symbol去重（取均值）
uniq_syms <- unique(rownames(expr_hf_mapped))
expr_hf_sym <- t(sapply(uniq_syms, function(g) {
  rows <- which(rownames(expr_hf_mapped) == g)
  if (length(rows) == 1) expr_hf_mapped[rows, ]
  else colMeans(expr_hf_mapped[rows, , drop = FALSE], na.rm = TRUE)
}))
if (!is.matrix(expr_hf_sym)) expr_hf_sym <- as.matrix(expr_hf_sym)
rownames(expr_hf_sym) <- uniq_syms
message(sprintf("探针→Symbol: %d genes", nrow(expr_hf_sym)))

# 分组：characteristics_ch1.1 = "heart failure: yes" vs "heart failure: no"
hf_col <- grep("heart.failure|heart failure", colnames(pd_hf),
               ignore.case = TRUE, value = TRUE)[1]
if (is.na(hf_col)) {
  # fallback: use characteristics_ch1.2 disease status
  hf_col <- grep("disease.status|disease status", colnames(pd_hf),
                 ignore.case = TRUE, value = TRUE)[1]
}
message(sprintf("分组列: %s", hf_col))
message(sprintf("分组值: %s", paste(unique(pd_hf[[hf_col]]), collapse=", ")))

hf_idx <- which(grepl("yes|DCM|ICM|dilated", pd_hf[[hf_col]], ignore.case = TRUE))
nf_idx <- which(grepl("no|non-failing", pd_hf[[hf_col]], ignore.case = TRUE))
message(sprintf("HF: %d, NF: %d", length(hf_idx), length(nf_idx)))

param_hf <- ssgseaParam(exprData = expr_hf_sym, geneSets = hallmark_sets,
                      minSize = 5, maxSize = 500)
ssgsea_hf <- gsva(param_hf, verbose = FALSE)
message(sprintf("HF ssGSEA完成: %d pathways × %d samples", nrow(ssgsea_hf), ncol(ssgsea_hf)))

hf_es <- apply(ssgsea_hf, 1, function(x) {
  h_m <- mean(x[hf_idx], na.rm = TRUE)
  n_m <- mean(x[nf_idx], na.rm = TRUE)
  ps <- sqrt((sd(x[hf_idx])^2 + sd(x[nf_idx])^2) / 2)
  if (ps == 0) 0 else (h_m - n_m) / ps
})

# ══════════════════════════════════════════════════════════════════════
# Part 4: 跨疾病通路效应量相关分析
# ══════════════════════════════════════════════════════════════════════
message("\n═══ 跨疾病通路效应量相关 ═══")

common_pw <- intersect(names(lihc_es), names(hf_es))
message(sprintf("共同通路: %d", length(common_pw)))

es_lihc_v <- lihc_es[common_pw]
es_hf_v <- hf_es[common_pw]

cor_all <- cor.test(es_lihc_v, es_hf_v, method = "spearman")
message(sprintf("所有通路 Spearman ρ = %.3f, p = %.4f", cor_all$estimate, cor_all$p.value))

# 翻译相关通路
trans_pw <- grep("RIBOSOM|TRANSLAT|PEPTIDE|AMINOACYL|RRNA|EIF|EEF|MYC_TARGET|MTORC1",
                 common_pw, value = TRUE, ignore.case = TRUE)
message(sprintf("翻译/核糖体相关通路: %d/%d", length(trans_pw), length(common_pw)))
for (pw in trans_pw) {
  message(sprintf("  %s: HCC d=%.2f, HF d=%.2f",
                  gsub("HALLMARK_|KEGG_|REACTOME_", "", pw),
                  es_lihc_v[pw], es_hf_v[pw]))
}

if (length(trans_pw) >= 3) {
  cor_trans <- cor.test(es_lihc_v[trans_pw], es_hf_v[trans_pw], method = "spearman")
  message(sprintf("翻译通路 Spearman ρ = %.3f, p = %.4f",
                  cor_trans$estimate, cor_trans$p.value))
}

# 置换检验
n_perm <- 10000
perm_rhos <- numeric(n_perm)
for (i in 1:n_perm) {
  perm_rhos[i] <- cor(es_lihc_v, sample(es_hf_v), method = "spearman")
}
perm_p <- mean(abs(perm_rhos) >= abs(cor_all$estimate))
message(sprintf("置换检验 p = %.4f", perm_p))

# ══════════════════════════════════════════════════════════════════════
# Part 5: 可视化
# ══════════════════════════════════════════════════════════════════════

is_trans <- common_pw %in% trans_pw

plot_df <- data.frame(
  pathway = common_pw,
  es_hcc = es_lihc_v,
  es_hf = es_hf_v,
  is_translation = is_trans,
  label = gsub("HALLMARK_|KEGG_|REACTOME_", "", common_pw),
  stringsAsFactors = FALSE
)

p1 <- ggplot(plot_df, aes(x = es_hcc, y = es_hf)) +
  geom_point(aes(color = is_translation, size = is_translation), alpha = 0.7) +
  geom_abline(slope = -1, intercept = 0, linetype = "dotted", color = "grey50", alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50", alpha = 0.25, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  ggrepel::geom_text_repel(
    data = subset(plot_df, is_translation),
    aes(label = label), size = 2.8, max.overlaps = 20, color = "#E41A1C"
  ) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#4DAF4A"),
                     labels = c("TRUE" = "Translation-related", "FALSE" = "Other")) +
  scale_size_manual(values = c("TRUE" = 3.5, "FALSE" = 1.8), guide = "none") +
  annotate("text", x = max(es_lihc_v) * 0.75, y = min(es_hf_v) * 0.9,
           label = paste0("All pathways: ρ=", round(cor_all$estimate, 3),
                          ", p=", round(cor_all$p.value, 4),
                          "\nPermutation p=", round(perm_p, 4)),
           size = 3.5, hjust = 0, color = "grey30") +
  labs(x = "Effect Size: HCC Tumor vs Normal (Cohen's d)",
       y = "Effect Size: HF vs NF (Cohen's d)",
       title = "Cross-Disease Pathway Activity Comparison (ssGSEA)",
       subtitle = paste0(length(trans_pw), " translation-related pathways highlighted")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(PROJ_DIR, "figures", "Figure_ssGSEA_cross_disease.png"),
       p1, width = 10, height = 8, dpi = 300)
message("\n✅ 图片已保存")

# ══════════════════════════════════════════════════════════════════════
# Part 6: 保存
# ══════════════════════════════════════════════════════════════════════
result <- list(
  cor_all = cor_all, cor_trans = if (exists("cor_trans")) cor_trans else NULL,
  perm_p = perm_p, trans_pathways = trans_pw,
  es_lihc = lihc_es, es_hf = hf_es,
  ssgsea_lihc = ssgsea_lihc, ssgsea_hf = ssgsea_hf, plot_df = plot_df
)
saveRDS(result, file.path(PROJ_DIR, "ssgsea_cross_disease_result.rds"))

cat("\n═══════════════════════════════════════════════\n")
cat("  ssGSEA 跨疾病通路分析完成\n")
cat(sprintf("  所有通路: ρ=%.3f, p=%.4f, 置换p=%.4f\n",
            cor_all$estimate, cor_all$p.value, perm_p))
if (exists("cor_trans")) {
  cat(sprintf("  翻译通路: ρ=%.3f, p=%.4f\n",
              cor_trans$estimate, cor_trans$p.value))
}
cat(sprintf("  翻译通路数: %d\n", length(trans_pw)))
cat("═══════════════════════════════════════════════\n")
