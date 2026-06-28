# 提取 HF WGCNA 参数 —— v4（手动分步法，绕过 WGCNA cor 兼容性问题）

library(WGCNA)
library(GEOquery)
library(hugene11sttranscriptcluster.db)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dynamicTreeCut)
enableWGCNAThreads(2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ── 1. 加载本地文件 ──
gz_file <- file.path(PROJ_DIR, "GSE57338_series_matrix.txt.gz")
cat("加载本地文件...\n")
gse <- getGEO(filename = gz_file, getGPL = FALSE)
eset <- gse

exprs_raw <- Biobase::exprs(eset)
pdata <- Biobase::pData(eset)
cat(sprintf("✅ %d 探针 × %d 样本\n", nrow(exprs_raw), ncol(exprs_raw)))

# ── 2. 表型分组 ──
grp_col <- grep("heart failure|disease|condition|group", colnames(pdata),
                ignore.case = TRUE, value = TRUE)[1]
grp <- pdata[[grp_col]]
is_hf <- grepl("yes|failure|HF|DCM|ICM|cardiomyopathy", grp, ignore.case = TRUE)
is_nf <- grepl("no|normal|non.fail|NF|control|healthy", grp, ignore.case = TRUE)
cat(sprintf("HF=%d, NF=%d, Other=%d\n", sum(is_hf), sum(is_nf), sum(!is_hf & !is_nf)))

# ── 3. 探针→基因 ──
cat("探针→基因...\n")
probe_ids <- as.character(rownames(exprs_raw))
sym_map <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
  keys = probe_ids, keytype = "PROBEID", columns = "SYMBOL")
sym_map <- sym_map[!is.na(sym_map$SYMBOL) & sym_map$SYMBOL != "", ]

keep <- rowSums(exprs_raw >= 4) >= ncol(exprs_raw) * 0.2
exprs_f <- exprs_raw[keep, ]
sym_f <- sym_map[sym_map$PROBEID %in% rownames(exprs_f), ]
exprs_f <- exprs_f[sym_f$PROBEID, , drop = FALSE]
uniq_syms <- unique(sym_f$SYMBOL)

exprs_gene <- t(sapply(uniq_syms, function(s) {
  probes <- sym_f$PROBEID[sym_f$SYMBOL == s]
  if (length(probes) == 1) exprs_f[probes, ] else colMeans(exprs_f[probes, , drop = FALSE])
}))
cat(sprintf("基因矩阵: %d × %d\n", nrow(exprs_gene), ncol(exprs_gene)))

# ── 4. 手动 WGCNA（不用 blockwiseModules）──
cat("\n═══ 手动 WGCNA ═══\n")
datExpr0 <- t(exprs_gene)

# 选方差最大的 5000 基因
vars <- apply(datExpr0, 2, var, na.rm = TRUE)
datExpr0 <- datExpr0[, order(vars, decreasing = TRUE)[1:min(3000, ncol(datExpr0))]]

# 软阈值选择
powers <- 1:20
sft_out <- pickSoftThreshold(datExpr0, powerVector = powers,
                              networkType = "signed", verbose = 2)

idx <- which(sft_out$fitIndices$SFT.R.sq >= 0.85)
beta_sel <- if (length(idx) > 0) sft_out$fitIndices$Power[min(idx)] else sft_out$fitIndices$Power[which.max(sft_out$fitIndices$SFT.R.sq)]
cat(sprintf("✅ β = %d (R²=%.3f)\n", beta_sel,
            sft_out$fitIndices$SFT.R.sq[which(sft_out$fitIndices$Power == beta_sel)]))

# ★ 手动构建 signed 邻接矩阵（用 stats::cor 绕过 WGCNA::cor 兼容性问题）
cat("构建邻接矩阵...\n")
cor_mat <- stats::cor(datExpr0, use = "pairwise.complete.obs")
adj <- (0.5 * (1 + cor_mat))^beta_sel  # signed adjacency

# TOM
cat("计算 TOM...\n")
TOM <- WGCNA::TOMsimilarity(adj, TOMType = "signed", verbose = 2)

# 聚类
cat("层次聚类...\n")
geneTree <- hclust(as.dist(1 - TOM), method = "average")

# 动态剪切
cat("模块检测...\n")
dynamicMods <- dynamicTreeCut::cutreeDynamic(
  dendro = geneTree, distM = 1 - TOM,
  deepSplit = 2, pamRespectsDendro = FALSE,
  minClusterSize = 30, verbose = 2)

moduleLabels <- dynamicMods
moduleColors <- WGCNA::labels2colors(moduleLabels)
n_mods <- length(unique(moduleColors))
cat(sprintf("✅ 模块数 = %d\n", n_mods))
print(sort(table(moduleColors), decreasing = TRUE))

# ── 5. GO 富集找翻译模块 ──
cat("\n═══ GO富集找翻译模块 ═══\n")
all_genes <- colnames(datExpr0)

for (mod in unique(moduleColors)) {
  mod_genes <- all_genes[moduleColors == mod]
  if (length(mod_genes) < 10) next
  ego <- enrichGO(gene = mod_genes, OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  if (!is.null(ego) && nrow(ego@result) > 0) {
    trans <- grepl("translat|ribosom|peptide.*biosyn|ribonucleoprotein|rRNA",
                   ego@result$Description, ignore.case = TRUE)
    if (any(trans)) {
      cat(sprintf("\n🔴 %s (%d genes) ★ 翻译模块!\n", mod, length(mod_genes)))
      top <- head(ego@result[trans, ], 6)
      for (i in seq_len(nrow(top)))
        cat(sprintf("  %s | p=%.2e | %s\n", top$ID[i], top$p.adjust[i], top$Description[i]))
    }
  }
}

# ── 6. Hub基因 ──
cat("\n═══ Hub基因分布 ═══\n")
for (g in c("EEF1A1","FAU","RPL39","RPL3","RPL32","RPL41","RPS28")) {
  if (g %in% all_genes)
    cat(sprintf("  %s → %s\n", g, moduleColors[which(all_genes == g)]))
  else
    cat(sprintf("  %s → 不在\n", g))
}

cat(sprintf("\n═══════════════════\n 汇总: β=%d | 模块数=%d\n═══════════════════\n",
            beta_sel, n_mods))
