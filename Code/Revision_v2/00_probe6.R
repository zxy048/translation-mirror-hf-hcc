# 00_probe6.R -- GSE57338 series-matrix structure: are rows probe IDs or symbols?
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
sm <- COHORTS$GSE57338$file
cat("file:", sm, " size:", file.size(sm), "\n\n")

con <- gzfile(sm, "rt")
hd  <- readLines(con, n = 60)
close(con)
cat("=== first 12 lines ===\n")
cat(paste(substr(hd[1:12], 1, 150), collapse = "\n"), "\n")

cat("\n=== lines that start with ! ===\n")
cat(paste(substr(grep("^!", hd, value = TRUE), 1, 120), collapse = "\n"), "\n")

# find the platform table
con <- gzfile(sm, "rt")
chunk <- readLines(con, n = 200)
close(con)
tbl_start <- grep("^!platform_table_begin", chunk)
cat("\nplatform_table_begin at line:", if (length(tbl_start)) tbl_start else "NOT in first 200", "\n")
if (length(tbl_start)) {
  cat(paste(substr(chunk[(tbl_start+1):(tbl_start+6)], 1, 160), collapse = "\n"), "\n")
}

# same for the sample table
s_start <- grep("^!series_matrix_table_begin", chunk)
cat("\nseries_matrix_table_begin at line:", if (length(s_start)) s_start else "NOT in first 200", "\n")
if (length(s_start)) {
  cat(paste(substr(chunk[(s_start+1):(s_start+5)], 1, 200), collapse = "\n"), "\n")
}

# how many samples / how many lines total
con <- gzfile(sm, "rt"); all <- readLines(con); close(con)
cat("\ntotal lines:", length(all), "\n")
cat("table_begin line:", grep("^!series_matrix_table_begin", all), "\n")
cat("table_end   line:", grep("^!series_matrix_table_end", all), "\n")

cat("\n=== GPL11532.soft.gz: does it map transcript-cluster -> symbol? ===\n")
gpl <- "D:/R_projects/GPL11532.soft.gz"
cat("exists:", file.exists(gpl), "\n")
con <- gzfile(gpl, "rt"); h <- readLines(con, n = 40); close(con)
cat(paste(substr(h[1:30], 1, 130), collapse = "\n"), "\n")
