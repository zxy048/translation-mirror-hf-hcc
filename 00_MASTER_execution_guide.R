# =============================================================================
# MASTER执行指南：论文增强分析完整流程
# 期刊目标：~5分生物信息学/医学基因组学期刊
# 运行顺序：按脚本编号依次执行
# =============================================================================

# 在RStudio中逐一source以下脚本：

# ── 第0步：环境检查 ──────────────────────────────────────────────────────────
cat("
╔══════════════════════════════════════════════════════════════════╗
║   HF-HCC 翻译转录信号论文——增强分析执行指南                      ║
║   目标期刊：5分SCI (BMC Medical Genomics/Frontiers系列)           ║
╚══════════════════════════════════════════════════════════════════╝

执行顺序：
  01 → 下载GSE141198并交互确认数据结构
  06 → 【可独立运行】方向一致性统计检验（最快获得新结果）
  02 → TGS + RPL39多队列验证
  03 → ssGSEA跨疾病通路比较
  04 → GSE141198独立WGCNA（核心新增分析）
  05 → 上游TF预测
  07 → 最终汇总和图表整合

建议执行时间线：
  Day 1: 01 + 06 + 02（获得初步结果）
  Day 2: 03 + 04（核心分析）
  Day 3: 05 + 07（收尾）
\n")

# ── 环境检查 ──────────────────────────────────────────────────────────────────
PROJ_DIR <- "D:/R_projects/revision_analysis"
setwd(PROJ_DIR)

source("D:/R_projects/revision_analysis/01_download_GSE141198.R")
# ⚠️ 运行完01后，必须在RStudio中执行以下交互检查：
#    View(pdata)   — 查看临床变量
#    head(expr_matrix) — 确认基因ID类型
#    table(pdata$characteristics_ch1.??) — 确认分组变量

# ── 06：方向一致性检验（最快获得新结果，可立即运行）──────────────────────
message("\n▶ 即将运行：06_direction_consistency_test.R")
message("  预计耗时：< 1分钟")
source("D:/R_projects/revision_analysis/06_direction_consistency_test.R")

# ── 02：TGS + RPL39验证（需要先匹配TCGA-LIHC临床数据）─────────────────────
message("\n▶ 即将运行：02_TGS_RPL39_validation.R")
message("  预计耗时：10-20分钟（含4队列分析）")
# ⚠️ 需要先在TCGA-LIHC中完成表达-预后数据匹配
source("D:/R_projects/revision_analysis/02_TGS_RPL39_validation.R")

# ── 03：ssGSEA通路比较（需要先加载HF和HCC表达数据）───────────────────────
message("\n▶ 即将运行：03_ssGSEA_cross_disease_pathway.R")
message("  预计耗时：30-60分钟（ssGSEA计算密集）")
source("D:/R_projects/revision_analysis/03_ssGSEA_cross_disease_pathway.R")

# ── 04：GSE141198独立WGCNA（最重要，仅当01完成数据确认后运行）────────────
message("\n▶ 即将运行：04_GSE141198_WGCNA.R")
message("  预计耗时：1-3小时（WGCNA TOM矩阵计算密集）")
source("D:/R_projects/revision_analysis/04_GSE141198_WGCNA.R")

# ── 05：上游TF预测 ──────────────────────────────────────────────────────────
message("\n▶ 即将运行：05_TF_upstream_prediction.R")
message("  预计耗时：5-15分钟")
source("D:/R_projects/revision_analysis/05_TF_upstream_prediction.R")

# ── 完成 ─────────────────────────────────────────────────────────────────────
cat("
╔══════════════════════════════════════════════════════════════════╗
║   ✅ 所有分析脚本已执行完毕                                       ║
║                                                                  ║
║   输出目录：                                                      ║
║     figures/ — 所有图片                                          ║
║     tables/  — 所有结果表                                        ║
║     WGCNA_GSE141198/ — GSE141198 WGCNA结果                       ║
║                                                                  ║
║   下一步：                                                        ║
║     1. 检查各分析输出结果                                         ║
║     2. 根据实际结果调整论文叙事                                  ║
║     3. 运行 07_final_figures.R 生成最终版图表                     ║
╚══════════════════════════════════════════════════════════════════╝
\n")
