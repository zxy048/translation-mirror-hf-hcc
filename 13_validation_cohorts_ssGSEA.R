# =============================================================================
# 脚本 13：新增验证队列 ssGSEA 分析
# 模块1: GSE116250 HF独立验证
# 模块2: GSE141910 HCM 同器官心肌对照 + GSE89377 Cirrhosis 同器官肝脏对照
# =============================================================================

library(GSVA)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(GEOquery)

set.seed(42)
PROJ_DIR <- "D:/R_projects/revision_analysis"
DATA_DIR <- file.path(PROJ_DIR, "validation_data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJ_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

# More robust download settings for Windows
options(timeout = 600)
options(download.file.method = "wininet")

# Cache file for incremental saves
CACHE_RDS <- file.path(PROJ_DIR, "validation_cache.rds")

# ══════════════════════════════════════════════════════════════════════
# Part 0: 加载基因集 + 已有结果
# ══════════════════════════════════════════════════════════════════════
message("═══ 加载基因集 ═══")

msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(msig_h$gene_symbol, msig_h$gs_name)
message(sprintf("Hallmark: %d", length(hallmark_sets)))

kegg_df <- msigdbr(species = "Homo sapiens", collection = "C2",
                   subcollection = "CP:KEGG_MEDICUS")
kegg_ribo <- kegg_df %>%
  filter(grepl("RIBOSOME", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>% unique()
hallmark_sets[["KEGG_RIBOSOME"]] <- kegg_ribo

reac_df <- msigdbr(species = "Homo sapiens", collection = "C2",
                   subcollection = "CP:REACTOME")
reac_trans <- reac_df %>%
  filter(grepl("TRANSLATION|PEPTIDE_CHAIN_ELONGATION|EUKARYOTIC_TRANSLATION|RIBOSOME|NONSENSE_MEDIATED|TRNA_AMINOACYLATION|RRNA_PROCESSING",
               gs_name, ignore.case = TRUE))
for (nm in unique(reac_trans$gs_name)) {
  hallmark_sets[[nm]] <- reac_trans$gene_symbol[reac_trans$gs_name == nm]
}
message(sprintf("Total gene sets: %d", length(hallmark_sets)))

# 加载已有结果
prev <- readRDS(file.path(PROJ_DIR, "ssgsea_cross_disease_result.rds"))
es_hcc  <- prev$es_lihc   # TCGA-LIHC Cohen's d
es_hf   <- prev$es_hf     # GSE57338 Cohen's d
trans_pw <- prev$trans_pathways

# 33 translation pathway names
trans_pw_names <- intersect(names(es_hcc), trans_pw)
message(sprintf("Existing: HCC ES=%d, HF ES=%d, Translation pathways=%d",
                length(es_hcc), length(es_hf), length(trans_pw_names)))

# ══════════════════════════════════════════════════════════════════════
# Helper: compute Cohen's d from ssGSEA matrix
# ══════════════════════════════════════════════════════════════════════
compute_cohens_d <- function(ssgsea_mat, case_idx, ctrl_idx) {
  apply(ssgsea_mat, 1, function(x) {
    cm <- mean(x[case_idx], na.rm = TRUE)
    nm <- mean(x[ctrl_idx], na.rm = TRUE)
    ps <- sqrt((var(x[case_idx]) + var(x[ctrl_idx])) / 2)
    if (ps == 0) 0 else (cm - nm) / ps
  })
}

# Helper: cross-disease Spearman correlation + permutation
run_cross_disease_test <- function(es_new, es_hcc, trans_pw_vec, label = "") {
  common <- intersect(names(es_new), names(es_hcc))
  es_a <- es_hcc[common]
  es_b <- es_new[common]

  trans_in_common <- intersect(common, trans_pw_vec)

  # All pathways
  cor_all <- cor.test(es_a, es_b, method = "spearman")

  # Translation pathways only
  cor_trans <- if (length(trans_in_common) >= 3) {
    cor.test(es_a[trans_in_common], es_b[trans_in_common], method = "spearman")
  } else NULL

  # Permutation test (translation pathways)
  n_perm <- 10000
  perm_rhos <- numeric(n_perm)
  for (i in 1:n_perm) {
    perm_rhos[i] <- cor(es_a[trans_in_common], sample(es_b[trans_in_common]),
                        method = "spearman")
  }
  perm_p <- if (!is.null(cor_trans)) {
    mean(abs(perm_rhos) >= abs(cor_trans$estimate))
  } else NA

  # Sign concordance: how many translation pathways show same direction?
  sign_concordance <- if (length(trans_in_common) > 0) {
    dir_a <- sign(es_a[trans_in_common])
    dir_b <- sign(es_b[trans_in_common])
    n_mirror <- sum(dir_a > 0 & dir_b < 0)  # HCC up, other down = mirror
    n_concordant <- sum(dir_a * dir_b > 0)   # same direction
    n_total <- length(trans_in_common)
    data.frame(n_total = n_total, n_mirror = n_mirror, n_concordant = n_concordant,
               pct_mirror = round(100 * n_mirror / n_total, 1))
  } else NULL

  list(cor_all = cor_all, cor_trans = cor_trans, perm_p = perm_p,
       sign_conc = sign_concordance, common_pathways = common,
       trans_in_common = trans_in_common, es_a = es_a, es_b = es_b,
       label = label)
}

# ══════════════════════════════════════════════════════════════════════
# Part 1: GSE116250 — HF Independent Validation (RNA-seq RPKM)
# ══════════════════════════════════════════════════════════════════════
message("\n═══ Part 1: GSE116250 HF Validation ═══")

gse116250_file <- file.path(DATA_DIR, "GSE116250_rpkm.txt.gz")
if (!file.exists(gse116250_file)) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116250/suppl/GSE116250_rpkm.txt.gz",
    gse116250_file, mode = "wb")
  message("Downloaded GSE116250 RPKM data")
}

rpkm_116250 <- read.table(gzfile(gse116250_file), header = TRUE,
                           check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
# Format: Gene (Ensembl) | Common_name (Symbol) | NF10 | NF11 | ...
# Use Common_name as gene symbol, keep highest expression per symbol
expr_cols_116250 <- grep("^NF|^DCM|^ICM|^HF", colnames(rpkm_116250), value = TRUE, ignore.case = TRUE)
if (length(expr_cols_116250) == 0) {
  # If no clear pattern, all columns except first 2 are expression
  expr_cols_116250 <- colnames(rpkm_116250)[-(1:2)]
}
message(sprintf("Expression columns: %d", length(expr_cols_116250)))

# Compute mean expression for dedup
expr_subset <- rpkm_116250[, expr_cols_116250, drop = FALSE]
expr_subset[] <- lapply(expr_subset, as.numeric)
rpkm_116250$mean_expr <- rowMeans(expr_subset, na.rm = TRUE)
rpkm_116250 <- rpkm_116250[order(rpkm_116250$mean_expr, decreasing = TRUE), ]
rpkm_116250 <- rpkm_116250[!duplicated(rpkm_116250$Common_name) & rpkm_116250$Common_name != "", ]
rownames(rpkm_116250) <- rpkm_116250$Common_name
# Keep only expression columns
expr_116250 <- as.matrix(log2(rpkm_116250[, expr_cols_116250, drop = FALSE] + 1))
rpkm_116250 <- NULL  # free memory
message(sprintf("GSE116250: %d genes × %d samples (log2 RPKM)", nrow(expr_116250), ncol(expr_116250)))

# Sample annotation from GEO - use locally cached series matrix
gse116250_sm <- file.path(DATA_DIR, "GSE116250_series_matrix.txt.gz")
suppressWarnings({
  gse116250_geo <- getGEO(filename = gse116250_sm, getGPL = FALSE)
})
pd_116250 <- pData(gse116250_geo)

# Inspect phenotype columns
message("Phenotype columns: ", paste(colnames(pd_116250), collapse=", "))
for (cc in grep("characteristics|diagnosis|title", colnames(pd_116250), ignore.case = TRUE, value = TRUE)) {
  message(sprintf("  %s => e.g. %s", cc, paste(head(unique(pd_116250[[cc]]), 5), collapse=" | ")))
}

# GSE116250 RPKM columns use descriptive names: NF1..NF14, CM1..CM50, ICM1..ICM13
# From the sample titles, we can directly classify by prefix
nf_cols_116250  <- grep("^NF\\d+", colnames(expr_116250), ignore.case = TRUE)
dcm_cols_116250 <- grep("^CM\\d+|^DCM\\d+", colnames(expr_116250), ignore.case = TRUE)
icm_cols_116250 <- grep("^ICM\\d+", colnames(expr_116250), ignore.case = TRUE)
hf_cols_116250 <- c(dcm_cols_116250, icm_cols_116250)

message(sprintf("GSE116250 by column prefix: NF=%d, DCM=%d, ICM=%d -> HF=%d",
                length(nf_cols_116250), length(dcm_cols_116250), length(icm_cols_116250),
                length(hf_cols_116250)))

if (length(hf_cols_116250) == 0 || length(nf_cols_116250) == 0) {
  # Fallback: use GEO phenotype to map
  message("Falling back to GEO phenotype mapping...")
  dx_col <- grep("diagnosis|disease|characteristic", colnames(pd_116250), ignore.case = TRUE, value = TRUE)[1]
  dx_values <- pd_116250[[dx_col]]
  message(sprintf("Diagnosis: %s", paste(unique(dx_values), collapse=" | ")))
  is_hf <- grepl("dilated|ischemic|cardiomyopathy|DCM|ICM", dx_values, ignore.case = TRUE)
  is_nf <- grepl("non.failing|nonfailing|control|normal|healthy", dx_values, ignore.case = TRUE)
  hf_acc <- pd_116250$geo_accession[is_hf]
  nf_acc <- pd_116250$geo_accession[is_nf]
  hf_cols_116250 <- which(colnames(expr_116250) %in% hf_acc)
  nf_cols_116250 <- which(colnames(expr_116250) %in% nf_acc)
  message(sprintf("GEO map: HF=%d cols, NF=%d cols", length(hf_cols_116250), length(nf_cols_116250)))
}

# Run ssGSEA
param_116250 <- ssgseaParam(exprData = expr_116250, geneSets = hallmark_sets,
                            minSize = 5, maxSize = 500)
ssgsea_116250 <- gsva(param_116250, verbose = FALSE)
es_116250 <- compute_cohens_d(ssgsea_116250, hf_cols_116250, nf_cols_116250)

message(sprintf("HF Cohen's d: median=%.3f, range=[%.3f, %.3f]",
                median(es_116250, na.rm = TRUE),
                min(es_116250, na.rm = TRUE), max(es_116250, na.rm = TRUE)))

result_116250 <- run_cross_disease_test(es_116250, es_hcc, trans_pw_names,
                                        "GSE116250 HF Validation")
message(sprintf("GSE116250 vs TCGA-LIHC: All ρ=%.3f, p=%.4f",
                result_116250$cor_all$estimate, result_116250$cor_all$p.value))
if (!is.null(result_116250$cor_trans)) {
  message(sprintf("  Translation: ρ=%.3f, p=%.4f, perm_p=%.4f",
                  result_116250$cor_trans$estimate,
                  result_116250$cor_trans$p.value,
                  result_116250$perm_p))
}
if (!is.null(result_116250$sign_conc)) {
  message(sprintf("  Sign conc: %d/%d mirror, %d concordant",
                  result_116250$sign_conc$n_mirror,
                  result_116250$sign_conc$n_total,
                  result_116250$sign_conc$n_concordant))
}

# ══════════════════════════════════════════════════════════════════════
# Part 2: GSE141910 — Heart Same-Organ Control (HCM vs NF, RNA-seq)
# Series Matrix for this RNA-seq dataset is empty; use supplementary CSV files
# ══════════════════════════════════════════════════════════════════════
message("\n═══ Part 2: GSE141910 HCM Heart Control ═══")

options(timeout = 300)  # 5 min timeout for large download

gse141910_tar <- file.path(DATA_DIR, "GSE141910_RAW.tar")
gse141910_dir <- file.path(DATA_DIR, "GSE141910")

if (!file.exists(gse141910_tar)) {
  tryCatch({
    download.file(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141910/suppl/GSE141910_RAW.tar",
      gse141910_tar, mode = "wb")
    message("Downloaded GSE141910 RAW.tar")
  }, error = function(e) {
    message("Download failed: ", e$message)
    message("Attempting curl download...")
    system(sprintf('curl -L -o "%s" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141910/suppl/GSE141910_RAW.tar" --max-time 300',
                   gse141910_tar))
  })
}

if (file.exists(gse141910_tar) && file.info(gse141910_tar)$size > 50000000) {
  dir.create(gse141910_dir, showWarnings = FALSE)
  untar(gse141910_tar, exdir = gse141910_dir)
  csv_files <- list.files(gse141910_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
  message(sprintf("GSE141910: %d CSV.gz files extracted", length(csv_files)))

  # Read each CSV: expected format is gene_id, expression_value
  # Build merged expression matrix
  sample1 <- read.csv(gzfile(csv_files[1]), header = TRUE, stringsAsFactors = FALSE, nrows = 3)
  message(sprintf("CSV columns: %s", paste(colnames(sample1), collapse=" | ")))

  expr_list <- list()
  gene_universe <- c()
  for (f in csv_files) {
    dat <- read.csv(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
    gene_col <- grep("gene|ensembl|symbol|id|Gene", colnames(dat), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(gene_col)) gene_col <- colnames(dat)[1]
    expr_col <- setdiff(colnames(dat), gene_col)[1]
    if (is.na(expr_col) || length(expr_col) == 0) expr_col <- colnames(dat)[2]

    sample_id <- gsub("\\.csv\\.gz$", "", basename(f))
    vals <- setNames(dat[[expr_col]], dat[[gene_col]])
    # Remove NAs and empty gene names
    vals <- vals[!is.na(names(vals)) & names(vals) != "" & !is.na(vals)]
    expr_list[[sample_id]] <- vals
    gene_universe <- union(gene_universe, names(vals))
  }

  message(sprintf("Gene universe: %d, Samples: %d", length(gene_universe), length(expr_list)))
  # Show sample gene IDs
  message(sprintf("Sample gene IDs: %s", paste(head(gene_universe, 10), collapse=", ")))

  # Build expression matrix
  expr_141910 <- matrix(NA, nrow = length(gene_universe), ncol = length(expr_list))
  rownames(expr_141910) <- gene_universe
  colnames(expr_141910) <- names(expr_list)
  for (nm in names(expr_list)) {
    genes <- names(expr_list[[nm]])
    expr_141910[genes, nm] <- expr_list[[nm]][genes]
  }
  # Remove rows with >50% NAs
  na_frac <- rowMeans(is.na(expr_141910))
  expr_141910 <- expr_141910[na_frac < 0.5, ]
  message(sprintf("Filtered: %d genes × %d samples", nrow(expr_141910), ncol(expr_141910)))

  # Map gene IDs to symbols. The X column likely contains Ensembl IDs or Entrez IDs
  # Try to detect and map
  library(org.Hs.eg.db)
  gene_ids <- rownames(expr_141910)

  # Check if Ensembl IDs (ENSG...) or Entrez (numeric)
  if (any(grepl("^ENSG", gene_ids))) {
    message("Detected Ensembl IDs, mapping to symbols...")
    gene_ids_clean <- sub("\\.\\d+$", "", gene_ids)
    sym_map <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_ids_clean,
                                     keytype = "ENSEMBL", columns = "SYMBOL")
    sym_map <- sym_map[!is.na(sym_map$SYMBOL) & !duplicated(sym_map$ENSEMBL), ]
    # Align
    keep <- match(sym_map$ENSEMBL, gene_ids_clean)
    expr_141910 <- expr_141910[keep, , drop = FALSE]
    rownames(expr_141910) <- sym_map$SYMBOL
  } else if (all(grepl("^\\d+$", gene_ids[1:10]))) {
    message("Detected Entrez IDs, mapping to symbols...")
    sym_map <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_ids,
                                     keytype = "ENTREZID", columns = "SYMBOL")
    sym_map <- sym_map[!is.na(sym_map$SYMBOL) & !duplicated(sym_map$ENTREZID), ]
    keep <- match(sym_map$ENTREZID, gene_ids)
    expr_141910 <- expr_141910[keep, , drop = FALSE]
    rownames(expr_141910) <- sym_map$SYMBOL
  } else {
    message("Unknown gene ID format, keeping as-is")
  }

  # Remove duplicates
  expr_141910 <- expr_141910[!duplicated(rownames(expr_141910)), ]
  message(sprintf("After gene symbol mapping: %d genes", nrow(expr_141910)))

  # Log2 if needed (raw counts >100 median)
  med_expr <- median(expr_141910, na.rm = TRUE)
  if (med_expr > 100) {
    expr_141910 <- log2(expr_141910 + 1)
    message(sprintf("Log2 transformed (median raw = %.1f)", med_expr))
  }
} else {
  message("GSE141910 RAW.tar not available, skipping this dataset")
  expr_141910 <- NULL
}

if (!is.null(expr_141910)) {
  # Get phenotype from locally cached series matrix
  gse141910_sm <- file.path(DATA_DIR, "GSE141910_series_matrix.txt.gz")
  suppressWarnings({
    gse141910_pd <- getGEO(filename = gse141910_sm, getGPL = FALSE)
  })
  pd_141910 <- pData(gse141910_pd)

  message("GSE141910 phenotype columns:")
  for (cc in grep("characteristics|diagnosis|title", colnames(pd_141910), ignore.case = TRUE, value = TRUE)) {
    message(sprintf("  %s => %s", cc, paste(head(unique(pd_141910[[cc]]), 4), collapse=" | ")))
  }

  # Classify: use characteristics_ch1 which has "etiology: ..."
  et_col <- grep("etiology|characteristics_ch1", colnames(pd_141910), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(et_col)) et_col <- grep("characteristics", colnames(pd_141910), ignore.case = TRUE, value = TRUE)[1]

  etiology <- pd_141910[[et_col]]
  message(sprintf("Using column '%s' for classification", et_col))

  is_hcm <- grepl("HCM|hypertrophic", etiology, ignore.case = TRUE)
  is_nf  <- grepl("Non.Failing|nonfailing|NF|control|healthy|donor", etiology, ignore.case = TRUE)
  is_other <- grepl("DCM|dilated|PPCM|peripartum", etiology, ignore.case = TRUE)

  hcm_acc <- pd_141910$geo_accession[is_hcm & !is_other]
  nf_acc  <- pd_141910$geo_accession[is_nf & !is_other]

  message(sprintf("GSE141910: HCM=%d, NF=%d, Other(DCM/PPCM)=%d",
                  length(hcm_acc), length(nf_acc), sum(is_other)))

  # Match: expression colnames = sample IDs from CSV (e.g., "C00039")
  # GEO accessions = "GSM4215858", title = "C00039"
  # We need to match CSV colnames to GEO accessions via titles
  title_to_acc <- setNames(pd_141910$geo_accession, pd_141910$title)

  hcm_titles <- pd_141910$title[is_hcm & !is_other]
  nf_titles  <- pd_141910$title[is_nf & !is_other]

  # CSV files: "GSM4215858_C00039.csv.gz" → sample_id="GSM4215858_C00039"
  # PD title: "C00039"
  # Need partial matching: check if PD title is a substring of colname
  hcm_cols <- which(sapply(colnames(expr_141910), function(x) any(sapply(hcm_titles, function(t) grepl(t, x, fixed = TRUE)))))
  nf_cols  <- which(sapply(colnames(expr_141910), function(x) any(sapply(nf_titles,  function(t) grepl(t, x, fixed = TRUE)))))

  message(sprintf("Expression colnames (first 5): %s", paste(head(colnames(expr_141910), 5), collapse=", ")))
  message(sprintf("HCM titles (first 5): %s", paste(head(hcm_titles, 5), collapse=", ")))

  message(sprintf("Aligned: HCM=%d, NF=%d cols", length(hcm_cols), length(nf_cols)))

  if (length(hcm_cols) >= 3 && length(nf_cols) >= 3) {
    param_141910 <- ssgseaParam(exprData = expr_141910, geneSets = hallmark_sets,
                                minSize = 5, maxSize = 500)
    ssgsea_141910 <- gsva(param_141910, verbose = FALSE)
    es_141910 <- compute_cohens_d(ssgsea_141910, hcm_cols, nf_cols)

    result_141910 <- run_cross_disease_test(es_141910, es_hcc, trans_pw_names,
                                            "GSE141910 HCM vs NF")
    message(sprintf("GSE141910 HCM vs TCGA-LIHC: All ρ=%.3f, p=%.4f",
                    result_141910$cor_all$estimate, result_141910$cor_all$p.value))
    if (!is.null(result_141910$cor_trans)) {
      message(sprintf("  Translation: ρ=%.3f, p=%.4f",
                      result_141910$cor_trans$estimate, result_141910$cor_trans$p.value))
    }
  } else {
    message("ERROR: Not enough HCM/NF samples in GSE141910")
    result_141910 <- NULL
    es_141910 <- NULL
  }
} else {
  result_141910 <- NULL
  es_141910 <- NULL
}

# ══════════════════════════════════════════════════════════════════════
# Part 3: GSE89377 — Liver Same-Organ Control (Cirrhosis vs Normal, Microarray)
# ══════════════════════════════════════════════════════════════════════
message("\n═══ Part 3: GSE89377 Cirrhosis Liver Control ═══")

# Use locally cached Series Matrix (processed data)
gse89377_sm <- file.path(DATA_DIR, "GSE89377_series_matrix.txt.gz")
suppressWarnings({
  eset_89377 <- getGEO(filename = gse89377_sm, getGPL = FALSE)
})
expr_89377_raw <- as.matrix(exprs(eset_89377))
pd_89377 <- pData(eset_89377)

message(sprintf("GSE89377 Series Matrix: %d probes × %d samples", nrow(expr_89377_raw), ncol(expr_89377_raw)))

# Inspect phenotype
message("GSE89377 phenotype:")
for (cc in grep("characteristics|title", colnames(pd_89377), ignore.case = TRUE, value = TRUE)) {
  message(sprintf("  %s => %s", cc, paste(head(unique(pd_89377[[cc]]), 4), collapse=" | ")))
}

# Classify: Normal (N-), Cirrhosis (CS-)
titles_89377 <- pd_89377$title
is_normal <- grepl("^N-|normal", titles_89377, ignore.case = TRUE)
is_cs    <- grepl("^CS-|cirrhosis|cirrhotic", titles_89377, ignore.case = TRUE)

message(sprintf("Normal: %d, Cirrhosis: %d", sum(is_normal), sum(is_cs)))

normal_cols <- which(is_normal)
cs_cols    <- which(is_cs)

message(sprintf("Aligned: Normal=%d, CS=%d columns", length(normal_cols), length(cs_cols)))

# Map Illumina probe ID → Gene Symbol using locally downloaded GPL annotation
gpl_file <- file.path(DATA_DIR, "GPL6947.annot.gz")
message("Parsing GPL6947 annotation...")

# Read the GEO annot file (GEO platform annotation format, tab-separated)
# Skip header lines, find platform_table_begin
annot_lines <- readLines(gzfile(gpl_file))
table_start <- grep("^!platform_table_begin", annot_lines)
annot_raw <- read.delim(gzfile(gpl_file), header = TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE, skip = table_start, sep = "\t")
message(sprintf("Annotation columns: %s", paste(colnames(annot_raw), collapse=" | ")))

# Find ID and Symbol columns
id_col <- grep("^ID$", colnames(annot_raw), ignore.case = TRUE, value = TRUE)[1]
gs_col <- grep("^Gene symbol$|^Symbol$|Gene.?symbol", colnames(annot_raw), ignore.case = TRUE, value = TRUE)[1]

if (!is.na(id_col) && !is.na(gs_col)) {
  message(sprintf("Using columns: ID='%s', Symbol='%s'", id_col, gs_col))
  probe_to_sym <- setNames(annot_raw[[gs_col]], annot_raw[[id_col]])
  probe_to_sym <- probe_to_sym[!is.na(probe_to_sym) & probe_to_sym != ""]
  message(sprintf("Probes with gene symbol: %d", length(probe_to_sym)))

  common <- intersect(rownames(expr_89377_raw), names(probe_to_sym))
  expr_89377 <- expr_89377_raw[common, , drop = FALSE]
  rownames(expr_89377) <- probe_to_sym[common]

  # Dedup
  uniq_g <- unique(rownames(expr_89377))
  expr_89377 <- t(sapply(uniq_g, function(g) {
    rows <- which(rownames(expr_89377) == g)
    if (length(rows) == 1) expr_89377[rows, , drop = TRUE]
    else colMeans(expr_89377[rows, , drop = FALSE], na.rm = TRUE)
  }))
  message(sprintf("After probe→gene: %d genes", nrow(expr_89377)))
} else {
  message("ERROR: Cannot find ID/Symbol columns in annotation")
  expr_89377 <- NULL
}

if (!is.null(expr_89377) && length(normal_cols) >= 3 && length(cs_cols) >= 3) {
  param_89377 <- ssgseaParam(exprData = expr_89377, geneSets = hallmark_sets,
                              minSize = 5, maxSize = 500)
  ssgsea_89377 <- gsva(param_89377, verbose = FALSE)
  es_89377 <- compute_cohens_d(ssgsea_89377, cs_cols, normal_cols)

  result_89377 <- run_cross_disease_test(es_89377, es_hcc, trans_pw_names,
                                          "GSE89377 Cirrhosis vs Normal")
  message(sprintf("GSE89377 Cirrhosis vs TCGA-LIHC: All ρ=%.3f, p=%.4f",
                  result_89377$cor_all$estimate, result_89377$cor_all$p.value))
  if (!is.null(result_89377$cor_trans)) {
    message(sprintf("  Translation: ρ=%.3f, p=%.4f, perm_p=%.4f",
                    result_89377$cor_trans$estimate,
                    result_89377$cor_trans$p.value,
                    result_89377$perm_p))
  }
  if (!is.null(result_89377$sign_conc)) {
    message(sprintf("  Sign conc: %d/%d mirror, %d concordant",
                    result_89377$sign_conc$n_mirror,
                    result_89377$sign_conc$n_total,
                    result_89377$sign_conc$n_concordant))
  }
} else {
  message("ERROR: Not enough Cirrhosis/Normal samples aligned in GSE89377")
  result_89377 <- NULL
}

# ══════════════════════════════════════════════════════════════════════
# Part 4: Also run GSE57338 (existing HF) cross-disease for reference
# ══════════════════════════════════════════════════════════════════════
message("\n═══ Part 4: GSE57338 Reference Correlation ═══")

result_57338 <- run_cross_disease_test(es_hf, es_hcc, trans_pw_names,
                                        "GSE57338 HF Discovery")
message(sprintf("GSE57338 vs TCGA-LIHC: All ρ=%.3f, Translation ρ=%.3f",
                result_57338$cor_all$estimate,
                result_57338$cor_trans$estimate))

# ══════════════════════════════════════════════════════════════════════
# Part 5: Combined Supplementary Figure
# ══════════════════════════════════════════════════════════════════════
message("\n═══ Part 5: Generate Combined Figure ═══")

build_plot_df <- function(result_obj, xlab_text) {
  if (is.null(result_obj)) return(NULL)
  is_trans <- result_obj$common_pathways %in% result_obj$trans_in_common
  data.frame(
    pathway = result_obj$common_pathways,
    es_hcc = result_obj$es_a,
    es_comp = result_obj$es_b,
    is_translation = is_trans,
    label = gsub("HALLMARK_|KEGG_|REACTOME_", "", result_obj$common_pathways),
    comparison = result_obj$label,
    stringsAsFactors = FALSE
  )
}

df_57338  <- build_plot_df(result_57338, "HF (GSE57338)")
df_116250 <- build_plot_df(result_116250, "HF Validation (GSE116250)")
df_141910 <- build_plot_df(result_141910, "HCM (GSE141910)")
df_89377  <- build_plot_df(result_89377, "Cirrhosis (GSE89377)")

all_plot_df <- do.call(rbind, Filter(Negate(is.null),
  list(df_57338, df_141910, df_89377, df_116250)))

# Add annotation
annot_df <- data.frame(
  comparison = c("GSE57338 HF Discovery", "GSE141910 HCM vs NF",
                 "GSE89377 Cirrhosis vs Normal", "GSE116250 HF Cohort"),
  rho_label = c(
    sprintf("ρ=%.3f, p=%.4f", result_57338$cor_trans$estimate, result_57338$cor_trans$p.value),
    if (!is.null(result_141910)) sprintf("ρ=%.3f, p=%.4f", result_141910$cor_trans$estimate,
                                          result_141910$cor_trans$p.value) else "N/A",
    if (!is.null(result_89377)) sprintf("ρ=%.3f, p=%.4f", result_89377$cor_trans$estimate,
                                         result_89377$cor_trans$p.value) else "N/A",
    sprintf("ρ=%.3f, p=%.4f", result_116250$cor_trans$estimate, result_116250$cor_trans$p.value)
  ),
  sign_label = c(
    sprintf("%d/%d mirror", result_57338$sign_conc$n_mirror, result_57338$sign_conc$n_total),
    if (!is.null(result_141910)) sprintf("%d/%d mirror", result_141910$sign_conc$n_mirror,
                                          result_141910$sign_conc$n_total) else "",
    if (!is.null(result_89377)) sprintf("%d/%d mirror", result_89377$sign_conc$n_mirror,
                                         result_89377$sign_conc$n_total) else "",
    sprintf("%d/%d mirror", result_116250$sign_conc$n_mirror, result_116250$sign_conc$n_total)
  ),
  category = c("HF-HCC", "Same-Organ Heart", "Same-Organ Liver", "HF-HCC")
)

all_plot_df$comparison <- factor(all_plot_df$comparison,
  levels = c("GSE57338 HF Discovery", "GSE141910 HCM vs NF",
             "GSE89377 Cirrhosis vs Normal", "GSE116250 HF Cohort"))

p_combined <- ggplot(all_plot_df, aes(x = es_hcc, y = es_comp)) +
  geom_point(aes(color = is_translation, size = is_translation), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50", alpha = 0.2, linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#999999"),
                     labels = c("TRUE" = "Translation", "FALSE" = "Other")) +
  scale_size_manual(values = c("TRUE" = 2.5, "FALSE" = 1.2), guide = "none") +
  geom_text(data = annot_df, aes(x = -Inf, y = Inf, label = rho_label),
            hjust = -0.1, vjust = 2, size = 3.2, inherit.aes = FALSE) +
  geom_text(data = annot_df, aes(x = -Inf, y = Inf, label = sign_label),
            hjust = -0.1, vjust = 4, size = 3.0, color = "grey40", inherit.aes = FALSE) +
  facet_wrap(~ comparison, ncol = 2, scales = "free") +
  labs(x = "HCC Tumor vs Normal (Cohen's d)",
       y = "Comparison Group Effect Size (Cohen's d)",
       title = "Disease Specificity of Cross-Disease Mirror Perturbation",
       subtitle = paste0("Translation pathways (n=", length(trans_pw_names),
                         "): red points. Mirror = HCC↑ & Comparison↓.")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 11))

ggsave(file.path(PROJ_DIR, "figures", "Figure_S5_Validation_SameOrgan_Controls.png"),
       p_combined, width = 14, height = 12, dpi = 300)
message("✅ Combined figure saved: Figure_S5_Validation_SameOrgan_Controls.png")

# ══════════════════════════════════════════════════════════════════════
# Part 6: Sign Concordance Summary Table
# ══════════════════════════════════════════════════════════════════════
message("\n═══ Part 6: Sign Concordance Summary ═══")

concordance_table <- data.frame(
  Comparison = c("HF Discovery (GSE57338)", "Heart Control (GSE141910 HCM)",
                 "Liver Control (GSE89377 Cirrhosis)", "HF Cohort (GSE116250)"),
  n_Translation_Pathways = c(
    result_57338$sign_conc$n_total,
    if (!is.null(result_141910)) result_141910$sign_conc$n_total else NA,
    if (!is.null(result_89377)) result_89377$sign_conc$n_total else NA,
    result_116250$sign_conc$n_total
  ),
  n_Mirror = c(
    result_57338$sign_conc$n_mirror,
    if (!is.null(result_141910)) result_141910$sign_conc$n_mirror else NA,
    if (!is.null(result_89377)) result_89377$sign_conc$n_mirror else NA,
    result_116250$sign_conc$n_mirror
  ),
  Spearman_rho = c(
    round(result_57338$cor_trans$estimate, 3),
    if (!is.null(result_141910)) round(result_141910$cor_trans$estimate, 3) else NA,
    if (!is.null(result_89377)) round(result_89377$cor_trans$estimate, 3) else NA,
    round(result_116250$cor_trans$estimate, 3)
  ),
  p_value = c(
    signif(result_57338$cor_trans$p.value, 3),
    if (!is.null(result_141910)) signif(result_141910$cor_trans$p.value, 3) else NA,
    if (!is.null(result_89377)) signif(result_89377$cor_trans$p.value, 3) else NA,
    signif(result_116250$cor_trans$p.value, 3)
  ),
  Permutation_p = c(
    signif(result_57338$perm_p, 3),
    if (!is.null(result_141910)) signif(result_141910$perm_p, 3) else NA,
    if (!is.null(result_89377)) signif(result_89377$perm_p, 3) else NA,
    signif(result_116250$perm_p, 3)
  ),
  Category = c("HF-HCC", "Same-Organ Heart", "Same-Organ Liver", "HF-HCC"),
  stringsAsFactors = FALSE
)

print(concordance_table)
write.csv(concordance_table, file.path(PROJ_DIR, "validation_sign_concordance.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════
# Part 7: Save All Results
# ══════════════════════════════════════════════════════════════════════
message("\n═══ Saving Results ═══")

validation_results <- list(
  es_hcc = es_hcc,
  es_hf_57338 = es_hf,
  es_hf_116250 = es_116250,
  es_hcm_141910 = if (exists("es_141910")) es_141910 else NULL,
  es_cirrhosis_89377 = if (exists("es_89377")) es_89377 else NULL,
  result_57338 = result_57338,
  result_116250 = result_116250,
  result_141910 = result_141910,
  result_89377 = result_89377,
  trans_pathways = trans_pw_names,
  concordance_table = concordance_table
)

saveRDS(validation_results, file.path(PROJ_DIR, "validation_cohorts_result.rds"))
message("✅ Results saved to validation_cohorts_result.rds")

cat("\n═══════════════════════════════════════════════\n")
cat("  Validation Analysis Complete\n")
cat("═══════════════════════════════════════════════\n")
cat(sprintf("  HF Discovery (GSE57338):    ρ=%.3f, %d/%d mirror\n",
            result_57338$cor_trans$estimate,
            result_57338$sign_conc$n_mirror, result_57338$sign_conc$n_total))
cat(sprintf("  HF Validation (GSE116250):  ρ=%.3f, %d/%d mirror\n",
            result_116250$cor_trans$estimate,
            result_116250$sign_conc$n_mirror, result_116250$sign_conc$n_total))
if (!is.null(result_141910)) {
  cat(sprintf("  Heart Control (GSE141910): ρ=%.3f, %d/%d mirror\n",
              result_141910$cor_trans$estimate,
              result_141910$sign_conc$n_mirror, result_141910$sign_conc$n_total))
}
if (!is.null(result_89377)) {
  cat(sprintf("  Liver Control (GSE89377):  ρ=%.3f, %d/%d mirror\n",
              result_89377$cor_trans$estimate,
              result_89377$sign_conc$n_mirror, result_89377$sign_conc$n_total))
}
cat("═══════════════════════════════════════════════\n")
