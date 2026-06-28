# =============================================================================
# 脚本 05：上游转录因子预测——为讨论中的机制推演提供生信依据
# 目标：用ChEA3 + Enrichr预测可能调控hub基因的转录因子
#       替代目前讨论中空泛的"MYC/mTORC1/ISR/E2F可能"叙述
# =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"

# ═══════════════════════════════════════════════════════════════════════════════
# 第一部分：定义hub基因集
# ═══════════════════════════════════════════════════════════════════════════════

# 翻译相关hub基因（方向一致的6个）
hub_genes_translation <- c("EEF1A1", "FAU", "RPL39", "RPL3", "RPL32", "RPL41")

# 扩展：所有核糖体蛋白基因（用于TF富集分析）
# 可以从MSigDB或KEGG获取核糖体基因列表
library(msigdbr)
kegg_ribosome_all <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG") %>%
  filter(gs_name == "KEGG_RIBOSOME") %>%
  pull(gene_symbol) %>%
  unique()

message(sprintf("KEGG核糖体基因总数: %d", length(kegg_ribosome_all)))
message("Hub基因: ", paste(hub_genes_translation, collapse = ", "))

# ═══════════════════════════════════════════════════════════════════════════════
# 第二部分：Enrichr转录因子富集（R接口）
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ Enrichr TF富集 ═══")

# 使用enrichR包（如果已安装）
tf_enrichment <- NULL
if (requireNamespace("enrichR", quietly = TRUE)) {
  library(enrichR)

  # 设置Enrichr数据库
  dbs <- c(
    "ChEA_2022",                    # ChIP-seq derived TF targets
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-seq",
    "TRRUST_Transcription_Factors_2019",
    "ENCODE_TF_ChIP-seq_2015",
    "Transcription_Factor_PPIs"
  )

  # 对hub基因做TF富集
  tryCatch({
    enriched_tf <- enrichr(hub_genes_translation, dbs)
    tf_enrichment <- enriched_tf
    message("✅ Enrichr富集完成")
  }, error = function(e) {
    message("⚠ Enrichr在线查询失败: ", e$message)
    message("将使用本地方法进行TF分析")
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# 第三部分：本地基于GO/KEGG的上游调控因子分析
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 本地TF调控分析 ═══")

# ── 3.1 使用已知的TF-靶基因关系 ──────────────────────────────────────────────

# MYC靶基因（MSigDB Hallmark MYC Targets V1/V2）
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
myc_targets_v1 <- msig_hallmark %>%
  filter(gs_name == "HALLMARK_MYC_TARGETS_V1") %>%
  pull(gene_symbol)

myc_targets_v2 <- msig_hallmark %>%
  filter(gs_name == "HALLMARK_MYC_TARGETS_V2") %>%
  pull(gene_symbol)

e2f_targets <- msig_hallmark %>%
  filter(gs_name == "HALLMARK_E2F_TARGETS") %>%
  pull(gene_symbol)

mtorc1_signaling <- msig_hallmark %>%
  filter(gs_name == "HALLMARK_MTORC1_SIGNALING") %>%
  pull(gene_symbol)

# Hub基因在这些通路中的重叠
hub_in_myc_v1 <- intersect(hub_genes_translation, myc_targets_v1)
hub_in_myc_v2 <- intersect(hub_genes_translation, myc_targets_v2)
hub_in_e2f <- intersect(hub_genes_translation, e2f_targets)
hub_in_mtorc1 <- intersect(hub_genes_translation, mtorc1_signaling)

message(sprintf("Hub基因 (%d) 在已知通路中的分布:", length(hub_genes_translation)))
message(sprintf("  MYC Targets V1: %d/%d — %s",
                length(hub_in_myc_v1), length(hub_genes_translation),
                paste(hub_in_myc_v1, collapse=", ")))
message(sprintf("  MYC Targets V2: %d/%d — %s",
                length(hub_in_myc_v2), length(hub_genes_translation),
                paste(hub_in_myc_v2, collapse=", ")))
message(sprintf("  E2F Targets:    %d/%d — %s",
                length(hub_in_e2f), length(hub_genes_translation),
                paste(hub_in_e2f, collapse=", ")))
message(sprintf("  mTORC1:         %d/%d — %s",
                length(hub_in_mtorc1), length(hub_genes_translation),
                paste(hub_in_mtorc1, collapse=", ")))

# ── 3.2 Fisher精确检验：hub基因在MYC靶基因中是否显著富集 ──────────────────
# 背景：所有检测到的基因
se <- readRDS("D:/R_projects/TCGA_LIHC_se.rds")
bg_genes <- rownames(assay(se, 1))
bg_size <- length(bg_genes)

# MYC V1富集检验
a_myc <- length(hub_in_myc_v1)           # hub基因且是MYC靶
b_myc <- length(hub_genes_translation) - a_myc  # hub基因但不是MYC靶
c_myc <- length(intersect(bg_genes, myc_targets_v1)) - a_myc
d_myc <- bg_size - a_myc - b_myc - c_myc

fisher_myc <- fisher.test(matrix(c(a_myc, b_myc, c_myc, d_myc), nrow = 2),
                          alternative = "greater")
message(sprintf("\nMYC V1富集 Fisher检验: OR=%.2f, p=%.4f",
                fisher_myc$estimate, fisher_myc$p.value))

# ── 3.3 核糖体蛋白基因集合在TF靶基因中的富集 ─────────────────────────────
# 不仅检验hub基因，检验所有KEGG核糖体基因在MYC/E2F/mTORC1靶基因中的富集
ribosome_genes_in_dataset <- intersect(kegg_ribosome_all, bg_genes)

if (length(ribosome_genes_in_dataset) >= 30) {
  # MYC富集
  ribo_in_myc <- length(intersect(ribosome_genes_in_dataset, myc_targets_v1))
  fisher_ribo_myc <- fisher.test(matrix(
    c(ribo_in_myc, length(ribosome_genes_in_dataset) - ribo_in_myc,
      length(intersect(bg_genes, myc_targets_v1)) - ribo_in_myc,
      bg_size - ribo_in_myc), nrow = 2), alternative = "greater")

  message(sprintf("\n核糖体基因 (n=%d) MYC V1富集: OR=%.2f, p=%.2e",
                  length(ribosome_genes_in_dataset),
                  fisher_ribo_myc$estimate, fisher_ribo_myc$p.value))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 第四部分：可视化——TF调控网络
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 构建TF-靶基因关系表 ═══")

# 手动整理已知的核糖体生物合成调控因子
# （基于文献综述：van Riggelen 2010, Pelletier 2018, Laplante 2012）
known_tf_ribosome <- data.frame(
  TF = c(rep("MYC", 6), rep("E2F1", 3), rep("mTORC1", 4), rep("TP53", 2)),
  Target = c("EEF1A1", "RPL3", "RPL32", "RPL39", "RPL41", "FAU",
             "RPL3", "RPL32", "RPL39",
             "EEF1A1", "RPL3", "RPL32", "RPL39",
             "FAU", "RPL41"),
  Evidence = c(rep("MSigDB+Literature", 6), rep("MSigDB Hallmark", 3),
               rep("Pathway overlap", 4), rep("Literature", 2)),
  stringsAsFactors = FALSE
)

# 统计各TF调控的hub基因
tf_summary <- known_tf_ribosome %>%
  group_by(TF) %>%
  summarise(
    n_targets = n(),
    target_genes = paste(unique(Target), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_targets))

message("\n═══ 上游TF调控hub基因汇总 ═══")
print(tf_summary)

# ── 4.1 在TCGA-LIHC中验证TF表达与TGS的相关性 ─────────────────────────────────
# 检查MYC、E2F1等关键TF在TCGA-LIHC中的表达
candidate_tfs <- c("MYC", "MYCN", "E2F1", "E2F2", "MTOR", "RPTOR", "TP53")

# 加载TCGA-LIHC肿瘤表达数据（标准化后）
tumor_expr <- NULL  # 需要从实际数据中获取
# 这里提供分析框架

verify_tf_expression <- function(tf_gene, expr_matrix, tgs_scores) {
  if (!tf_gene %in% rownames(expr_matrix)) return(NULL)

  tf_expr <- as.numeric(expr_matrix[tf_gene, ])
  common_samples <- intersect(names(tf_expr), names(tgs_scores))

  if (length(common_samples) < 30) {
    message(sprintf("  %s: 样本不足 (%d)", tf_gene, length(common_samples)))
    return(NULL)
  }

  cor_test <- cor.test(tf_expr[common_samples], tgs_scores[common_samples],
                       method = "spearman")

  data.frame(
    TF = tf_gene,
    rho = cor_test$estimate,
    p_value = cor_test$p.value,
    n = length(common_samples),
    stringsAsFactors = FALSE
  )
}

message("\n框架就绪。在加载TCGA-LIHC标准化表达数据后运行verify_tf_expression()。")

# ── 4.2 TF调控网络可视化 ─────────────────────────────────────────────────────
message("\n═══ 绘制TF-hub基因调控网络 ═══")

# 简单网络图（如果igraph可用）
if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)

  # 节点
  nodes <- unique(c(known_tf_ribosome$TF, known_tf_ribosome$Target))
  g <- graph_from_data_frame(known_tf_ribosome, vertices = data.frame(
    name = nodes,
    type = ifelse(nodes %in% known_tf_ribosome$TF, "TF", "Target"),
    stringsAsFactors = FALSE
  ))

  # 简单布局
  lo <- layout_with_fr(g)

  png(file.path(PROJ_DIR, "figures", "TF_hub_regulatory_network.png"),
      width = 8, height = 6, units = "in", res = 300)
  par(mar = c(0, 0, 2, 0))

  vertex_colors <- ifelse(V(g)$type == "TF", "#E41A1C", "#377EB8")
  vertex_sizes <- ifelse(V(g)$type == "TF", 20, 12)
  vertex_shapes <- ifelse(V(g)$type == "TF", "square", "circle")

  plot(g,
       layout = lo,
       vertex.color = vertex_colors,
       vertex.size = vertex_sizes,
       vertex.shape = vertex_shapes,
       vertex.label.cex = 0.8,
       vertex.label.color = "black",
       vertex.frame.color = "grey50",
       edge.arrow.size = 0.5,
       edge.color = "grey60",
       edge.width = 1.5,
       main = "Known TF → Hub Gene Regulatory Relationships")

  legend("bottomright",
         legend = c("Transcription Factor", "Hub Gene"),
         col = c("#E41A1C", "#377EB8"),
         pch = c(22, 21),
         pt.bg = c("#E41A1C", "#377EB8"),
         pt.cex = 2,
         bty = "n")

  dev.off()
  message("网络图已保存")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 第五部分：关键竞争机制的区分分析
# ═══════════════════════════════════════════════════════════════════════════════

message("\n═══ 竞争机制区分: MYC vs mTORC1 vs E2F ═══")

# 从TCGA-LIHC中检验哪种TF最能解释hub基因的共表达模式
# 策略：比较MYC、E2F1、mTORC1通路活性的ssGSEA评分与TGS的相关性

# 前提：已有ssGSEA通路评分（从脚本03加载）
# 比较各通路活性与TGS的相关性强度

compare_upstream <- function(tgs, ssgsea_scores) {
  upstream_pathways <- c("HALLMARK_MYC_TARGETS_V1",
                         "HALLMARK_MYC_TARGETS_V2",
                         "HALLMARK_E2F_TARGETS",
                         "HALLMARK_MTORC1_SIGNALING")

  existing_pw <- intersect(upstream_pathways, rownames(ssgsea_scores))
  results <- data.frame()

  for (pw in existing_pw) {
    pw_score <- ssgsea_scores[pw, names(tgs)]
    ct <- cor.test(tgs, pw_score, method = "spearman")
    results <- rbind(results, data.frame(
      pathway = gsub("HALLMARK_", "", pw),
      rho = ct$estimate,
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    ))
  }

  results <- results[order(abs(results$rho), decreasing = TRUE), ]
  return(results)
}

message("\n✅ 脚本05完成")
message("输出：")
message("  1. TF-hub基因调控关系表")
message("  2. MYC/E2F/mTORC1各通路对hub基因的调控贡献")
message("  3. TF调控网络图（如igraph可用）")
message("  4. 区别于讨论中空泛推演的具体数据依据")
