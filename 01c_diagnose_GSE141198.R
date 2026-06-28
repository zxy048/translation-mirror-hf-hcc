# =============================================================================
# 脚本 01c：诊断GSE141198数据结构
# =============================================================================

PROJ_DIR <- "D:/R_projects/revision_analysis"

library(GEOquery)

# 方法1：重新下载，但这次不假设GSE Matrix格式
gse_list <- getGEO("GSE141198", GSEMatrix = TRUE, getGPL = FALSE)
message("GSE对象数量: ", length(gse_list))
message("GSE对象类: ", class(gse_list))
message("第一个元素类: ", class(gse_list[[1]]))

eset <- gse_list[[1]]
message("\nExpressionSet结构:")
message("  assayData元素: ", paste(names(assayData(eset)), collapse=", "))
message("  exprs维度: ", nrow(exprs(eset)), " × ", ncol(exprs(eset)))

# 检查是否有其他assay
print(assayData(eset))

# 检查feature data
fdata <- fData(eset)
message("\nfData列:")
print(colnames(fdata))
message("fData行数: ", nrow(fdata))

# 检查是否有表达数据在fData中
if (nrow(fdata) > 0 && ncol(fdata) > 0) {
  cat("\nfData前3列的前5行:\n")
  print(head(fdata[, 1:min(5, ncol(fdata))], 5))
}

# 方法2：直接从GEO series文件获取（如果Matrix不可用）
# 尝试用getGEO(filename=...)方式
message("\n═══ 尝试替代方案 ═══")

# 下载原始series matrix文件
tmp_dir <- tempdir()
gse_raw <- getGEOSuppFiles("GSE141198", baseDir = tmp_dir, makeDirectory = TRUE)
message("Supplementary文件:")
print(list.files(file.path(tmp_dir, "GSE141198"), full.names = TRUE))

# 方法3：直接检查是否在assayData的exprs之外有数据
message("\n═══ 深度检查assayData ═══")
ad <- assayData(eset)
message("assayData存储模式: ", storageMode(ad))
ls_ad <- ls(ad)
message("assayData环境中的对象: ", paste(ls_ad, collapse=", "))

# 逐个检查
for (nm in ls_ad) {
  obj <- get(nm, envir = ad)
  message(sprintf("  %s: class=%s, dim=%s", nm, class(obj)[1],
                  paste(dim(obj), collapse="×")))
}

# 检查phenoData
pd <- pData(eset)
message("\nphenoData行数: ", nrow(pd))
message("phenoData列数: ", ncol(pd))
cat("\nphenoData OS相关列的唯一值:\n")
cat("  os days:ch1 范围:", range(as.numeric(pd[["os days:ch1"]]), na.rm=TRUE), "\n")
cat("  os event:ch1 分布:\n")
print(table(pd[["os event:ch1"]]))
