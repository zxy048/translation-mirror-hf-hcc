# =============================================================================
# 脚本 01：下载 GSE141198（HCC RNA-seq, n=148, 台湾队列）
# 用途：作为第三个独立HCC数据集，用于de novo WGCNA和TGS预后验证
# 运行环境：R 4.6.0, GEOquery
# =============================================================================

library(GEOquery)
library(dplyr)
library(tidyr)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
dir.create(PROJ_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. 下载 GSE141198 ─────────────────────────────────────────────────────────
message("═══ 下载 GSE141198 ═══")

# 下载series matrix（表达矩阵+临床注释）
gse <- getGEO("GSE141198", GSEMatrix = TRUE, getGPL = TRUE)
eset <- gse[[1]]  # ExpressionSet

# 提取表达矩阵
# GSE141198 数据以log2(FPKM+1)或类似格式提供
expr_matrix <- exprs(eset)
message(sprintf("表达矩阵维度: %d genes × %d samples", nrow(expr_matrix), ncol(expr_matrix)))

# 提取表型数据
pdata <- pData(eset)

# ── 2. 提取关键临床变量 ──────────────────────────────────────────────────────
message("\n═══ 提取临床变量 ═══")

# 检查列名
cat("可用临床变量:\n")
print(colnames(pdata))

# 提取关键变量（GSE141198的实际列名可能略有不同）
# 常见模式：characteristics_ch1 系列列
# 关键信息通常包括：OS time, OS status, etiology, stage 等

# 尝试提取生存数据
surv_cols <- grep("overall survival|os\\b|survival|follow.?up|event|status|dead|death",
                  colnames(pdata), ignore.case = TRUE, value = TRUE)
cat("\n可能的生存相关列:\n")
print(surv_cols)

# 打印前几行查看数据结构
cat("\n表型数据前3列的前5行:\n")
print(head(pdata[, 1:min(8, ncol(pdata))], 3))

# ── 3. 标准化临床变量 ────────────────────────────────────────────────────────
# 注：GSE141198的临床数据需要通过characteristics_ch1提取
# 以下代码在RStudio中交互运行时会显示实际列名

# ── 4. 保存原始数据 ──────────────────────────────────────────────────────────
saveRDS(gse, file.path(PROJ_DIR, "GSE141198_raw.rds"))
saveRDS(expr_matrix, file.path(PROJ_DIR, "GSE141198_expr.rds"))
saveRDS(pdata, file.path(PROJ_DIR, "GSE141198_pdata.rds"))

message("\n✅ GSE141198下载完成。请在RStudio中查看pdata列名后继续。")
message("文件保存至: ", PROJ_DIR)

# ── 5. 基因ID映射准备（用于后续WGCNA）────────────────────────────────────────
# GSE141198可能使用gene symbol或Ensembl ID
# 检查前5个基因名确认ID类型
cat("\n前20个基因名（确认ID类型）:\n")
print(head(rownames(expr_matrix), 20))

# 保存session信息
writeLines(capture.output(sessionInfo()),
           file.path(PROJ_DIR, "sessionInfo_01_download.txt"))
