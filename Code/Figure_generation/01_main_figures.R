# =============================================================================
# 脚本 07：最终图表整合——为修订稿生成4张核心新图 + 2张新表
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # 拼图

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
fig_dir <- file.path(PROJ_DIR, "figures")
tbl_dir <- file.path(PROJ_DIR, "tables")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure X1（新增）：跨疾病翻译基因方向一致性 + 三数据集WGCNA复现
# ═══════════════════════════════════════════════════════════════════════════════
message("═══ Figure X1: Cross-disease translational gene convergence ═══")

# 此图整合：
# Panel A: 7基因跨疾病log2FC散点图（来自脚本06）
# Panel B: GSE141198独立WGCNA模块富集（来自脚本04）
# Panel C: 三数据集翻译模块基因重叠Venn图

# ═══════════════════════════════════════════════════════════════════════════════
# Figure X2（新增）：TGS跨四队列预后Meta分析
# ═══════════════════════════════════════════════════════════════════════════════
message("═══ Figure X2: TGS prognostic meta-analysis ═══")

# 此图整合：
# Panel A-D: 四队列KM曲线（TCGA-LIHC, GSE14520, GSE76427, GSE141198）
# Panel E: Meta分析森林图（随机效应模型）
# Panel F: 按病因分层的亚组分析

# ═══════════════════════════════════════════════════════════════════════════════
# Figure X3（新增）：ssGSEA跨疾病通路活性比较
# ═══════════════════════════════════════════════════════════════════════════════
message("═══ Figure X3: Cross-disease pathway activity comparison ═══")

# 此图整合：
# Panel A: Hallmark通路效应量散点图（HF vs HCC，标注翻译相关通路）
# Panel B: 翻译相关通路ssGSEA分数箱线图（两病并排比较）
# Panel C: 置换检验null分布 + 观察值

# ═══════════════════════════════════════════════════════════════════════════════
# Figure X4（新增）：上游TF预测 + 调控网络
# ═══════════════════════════════════════════════════════════════════════════════
message("═══ Figure X4: Upstream TF prediction ═══")

# 此图整合：
# Panel A: MYC/mTORC1/E2F通路与TGS的相关性barplot
# Panel B: TF-hub基因调控网络
# Panel C: MYC表达与TGS的散点图（TCGA-LIHC验证）

# ═══════════════════════════════════════════════════════════════════════════════
# Table X1（新增）：多队列TGS预后验证结果汇总
# ═══════════════════════════════════════════════════════════════════════════════
message("═══ Table X1: Multi-cohort TGS validation ═══")

create_table_x1 <- function(results_list) {
  # results_list: 每个队列的TGS分析结果列表
  tab <- data.frame(
    Cohort = c("TCGA-LIHC", "GSE14520", "GSE76427", "GSE141198"),
    N = c(371, 221, 115, 148),
    Platform = c("RNA-seq", "Microarray", "RNA-seq", "RNA-seq"),
    Events = NA_integer_,
    HR_per_IQR = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    Cox_p = NA_real_,
    Logrank_p = NA_real_,
    Concordance = NA_real_,
    stringsAsFactors = FALSE
  )
  return(tab)
}

# Table X2（新增）：跨疾病通路SSGSEA效应量
# Table 汇总每个翻译相关通路在HF和HCC中的效应量、方向和一致性

message("\n✅ 最终图表框架就绪")
message("在实际数据填充后，运行本脚本生成最终稿图表")

# ═══════════════════════════════════════════════════════════════════════════════
# 附录：新旧图表对照
# ═══════════════════════════════════════════════════════════════════════════════
cat("
───────────────────────────────────────────────────────────────
            修订稿图表对照表
───────────────────────────────────────────────────────────────
Figure 1:  研究设计流程图（保留，更新加入GSE141198）
Figure 2:  WGCNA树状图+模块性状热图（保留）
Figure 3:  【新增】跨疾病基因方向一致性 + 三数据集WGCNA复现
Figure 4:  hub基因跨疾病表达森林图（保留，更新数据）
Figure 5:  【新增】TGS预后Meta分析（四队列）
Figure 6:  【新增】ssGSEA跨疾病通路活性比较
Figure 7:  【新增】上游TF预测网络
Figure 8:  KM生存曲线（保留，增加TGS分组）
───────────────────────────────────────────────────────────────
Sup Fig 1-7: 保留原有补充图
Sup Fig 8:  模块保留性（移入正文Figure 3）
Sup Fig 9:  GSE141198独立WGCNA详细结果
Sup Fig 10: GSE116174 TGS验证（额外验证队列）
Sup Fig 11: 置换检验详细结果
───────────────────────────────────────────────────────────────
")
