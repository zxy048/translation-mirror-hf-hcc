# 10_generate_table_S1_GO_comparison.R
# Table S1: GSE57338 black module vs GSE141198 blue module GO enrichment comparison
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# ═══════════════════════════════════════════════════════
# 1. GSE57338 HF black module
# ═══════════════════════════════════════════════════════
cat("═══ 1. GSE57338 Black Module ═══\n")

# 加载 v6.1 WGCNA 结果
hf_file <- "D:/R_projects/ML_output/wgcna_full_network.RData"
load(hf_file)

# 查看加载了什么对象
cat("Loaded objects:\n")
print(ls())

# 找到模块颜色和基因名
# 通常 blockwiseModules 输出包含: net, moduleColors, MEs 等
if (exists("moduleColors")) {
  cat(sprintf("\nmoduleColors length: %d\n", length(moduleColors)))
  cat("Module distribution:\n")
  print(sort(table(moduleColors), decreasing = TRUE))

  # 找 black 模块
  black_genes <- names(moduleColors)[moduleColors == "black"]
  cat(sprintf("\nBlack module: %d genes\n", length(black_genes)))
  cat("First 20 black genes:\n")
  print(head(black_genes, 20))
} else if (exists("net")) {
  cat("\nUsing net$colors...\n")
  cols <- labels2colors(net$colors)
  cat(sprintf("net$colors length: %d\n", length(cols)))
  print(sort(table(cols), decreasing = TRUE))

  black_idx <- which(cols == "black")
  cat(sprintf("\nBlack module: %d genes\n", length(black_idx)))

  # 基因名可能来自 net 的其他组件
  if (!is.null(names(net$colors))) {
    black_genes <- names(net$colors)[black_idx]
  } else if (!is.null(names(cols))) {
    black_genes <- names(cols)[black_idx]
  } else {
    cat("⚠ 基因名未命名，检查 net 结构...\n")
    cat("net 组件:\n")
    print(names(net))
  }
}

# 检查基因名格式
if (exists("black_genes") && length(black_genes) > 0) {
  cat(sprintf("\nBlack genes sample:\n"))
  print(head(black_genes, 30))
  cat(sprintf("Gene ID format check (first 10): %s\n",
              paste(head(black_genes, 10), collapse=", ")))
}

# ═══════════════════════════════════════════════════════
# 2. 检查是否需要ID转换
# ═══════════════════════════════════════════════════════
cat("\n═══ 2. Gene ID Check ═══\n")

# HF 数据是 Affymetrix 探针 → Symbol，所以应该是 SYMBOL
# 但也可能是探针ID
if (exists("black_genes")) {
  # 测试前几个基因是否能被 enrichGO 识别
  test_genes <- head(black_genes, 50)

  # 检查是否为 SYMBOL
  sym_check <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = test_genes,
                          columns = "SYMBOL", keytype = "SYMBOL")
  )
  cat(sprintf("SYMBOL→SYMBOL 映射: %d/%d\n",
              sum(!is.na(sym_check$SYMBOL)), length(test_genes)))

  # 如果不是 SYMBOL，检查是否为 ENSEMBL 或 ENTREZ
  if (sum(!is.na(sym_check$SYMBOL)) < 5) {
    cat("可能是其他ID类型，尝试 ENTREZ...\n")
    ent_check <- suppressMessages(
      AnnotationDbi::select(org.Hs.eg.db, keys = test_genes,
                            columns = "SYMBOL", keytype = "ENTREZID")
    )
    cat(sprintf("ENTREZ→SYMBOL 映射: %d/%d\n",
                sum(!is.na(ent_check$SYMBOL)), length(test_genes)))
  }
}

# ═══════════════════════════════════════════════════════
# 3. GSE141198 Blue Module GO (复用 Figure 2B 逻辑)
# ═══════════════════════════════════════════════════════
cat("\n═══ 3. GSE141198 Blue Module GO ═══\n")

expr_wgcna_141198 <- readRDS(file.path(PROJ_DIR, "GSE141198_wgcna_input.rds"))
wgcna_res_141198 <- readRDS(file.path(PROJ_DIR, "GSE141198_WGCNA_result.rds"))

all_genes_141198 <- rownames(expr_wgcna_141198)
block1_idx_141198 <- wgcna_res_141198$net$blockGenes[[1]]
moduleColors_141198 <- wgcna_res_141198$moduleColors
mc_block1 <- moduleColors_141198[block1_idx_141198]
blue_genes <- all_genes_141198[block1_idx_141198][mc_block1 == "blue"]
cat(sprintf("GSE141198 Blue module: %d genes\n", length(blue_genes)))

ego_blue <- enrichGO(
  gene = blue_genes, OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL", ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2
)

if (!is.null(ego_blue) && nrow(ego_blue@result) > 0) {
  cat(sprintf("GSE141198 Blue GO terms: %d\n", nrow(ego_blue@result)))
  # 标记翻译相关
  blue_df <- as.data.frame(ego_blue)
  blue_df$is_translation <- grepl("translat|ribosom|peptide.*biosyn|ribonucleoprotein|rRNA|translational",
                                   blue_df$Description, ignore.case = TRUE)
  cat(sprintf("  Translation-related: %d terms\n", sum(blue_df$is_translation)))
}

# ═══════════════════════════════════════════════════════
# 4. 生成 Table S1
# ═══════════════════════════════════════════════════════
cat("\n═══ 4. Table S1 Generation ═══\n")

# 等步骤1确认 black_genes 的ID类型后再继续
cat("\n请确认上述输出：\n")
cat("  1. HF black module 基因数是否 > 0？\n")
cat("  2. 基因ID格式是什么（SYMBOL / ENTREZ / 探针ID）？\n")
cat("  3. 如果是探针ID，会报后面补充转换代码\n")
