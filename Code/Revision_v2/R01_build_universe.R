# =============================================================================
# R01_build_universe.R
# -----------------------------------------------------------------------------
# Load both discovery platforms, apply a DECLARED expression filter (no variance
# pre-filter -- that is what systematically removed ribosomal protein genes in
# v1), map to symbols, and build the common gene universe U.
#
# Emits a full loss chain: features -> detected -> symbol-mapped -> in U,
# plus per-tier coverage of the canonical translation sets.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R01 start ===")

# I/O helpers (parse_series_matrix, parse_soft, collapse_by_mean, ...) come from
# R00_io.R, sourced via R00_config.R. They were first written here; extracting
# them is what keeps R01 and R05 from drifting apart.

# ===========================================================================
# 1. HF / GSE57338  (array, GPL11532)
# ===========================================================================
log_msg("--- HF: GSE57338 ---")
sm <- parse_series_matrix(COHORTS$GSE57338$file)
E_hf_raw <- sm$expr
ann_hf   <- sm$ann
log_msg("series matrix: ", nrow(E_hf_raw), " features x ", ncol(E_hf_raw), " samples")

stopifnot(ncol(E_hf_raw) == 313)
cat("\nHF sample annotation columns:", paste(setdiff(names(ann_hf), "gsm"), collapse = ", "), "\n")
if ("heart.failure" %in% names(ann_hf)) print(table(ann_hf$heart.failure))
if ("disease.status" %in% names(ann_hf)) print(table(ann_hf$disease.status))

# group assignment
ann_hf$group <- ifelse(grepl("^yes", ann_hf$heart.failure, ignore.case = TRUE), "HF",
                ifelse(grepl("^no",  ann_hf$heart.failure, ignore.case = TRUE), "NF", NA))
ann_hf$etiology <- ifelse(grepl("dilated", ann_hf$disease.status, ignore.case = TRUE), "DCM",
                   ifelse(grepl("ischemic", ann_hf$disease.status, ignore.case = TRUE), "ICM",
                   ifelse(grepl("non-?failing|normal", ann_hf$disease.status, ignore.case = TRUE), "NF", "other")))
cat("\nHF group x etiology:\n"); print(table(ann_hf$group, ann_hf$etiology, useNA = "ifany"))

n_hf_dis <- sum(ann_hf$group == "HF", na.rm = TRUE)
log_msg("HF disease samples: ", n_hf_dis, "  (config says ", COHORTS$GSE57338$n_disease, ")")

# probe -> symbol
suppressPackageStartupMessages(library(hugene11sttranscriptcluster.db))
sym_hf_all <- suppressWarnings(AnnotationDbi::mapIds(
  hugene11sttranscriptcluster.db, keys = rownames(E_hf_raw),
  column = "SYMBOL", keytype = "PROBEID", multiVals = "first"))

# --- expression filter (array): above the global median in >= 50% of samples
gmed <- stats::median(E_hf_raw, na.rm = TRUE)
det_hf <- rowSums(E_hf_raw > gmed, na.rm = TRUE) >= 0.5 * ncol(E_hf_raw)
log_msg(sprintf("array filter: global median=%.3f -> %d/%d probes detected (%.1f%%)",
                gmed, sum(det_hf), nrow(E_hf_raw), 100 * mean(det_hf)))

E_hf_det <- E_hf_raw[det_hf, , drop = FALSE]
sym_hf_det <- sym_hf_all[det_hf]
hf_coll <- collapse_by_mean(E_hf_det, sym_hf_det)
E_hf <- hf_coll$expr
log_msg("HF after collapse: ", nrow(E_hf), " genes")

# ===========================================================================
# 2. HCC / GSE141198  (RNA-seq, already symbol-mapped)
# ===========================================================================
log_msg("--- HCC: GSE141198 ---")
cnt <- readRDS(COHORTS$GSE141198$file)          # NB: config file field points at vst
vst <- readRDS(file.path(PROJ, "GSE141198_vst.rds"))
cnt <- readRDS(file.path(PROJ, "GSE141198_counts_clean.rds"))   # raw counts
stopifnot(identical(rownames(cnt), rownames(vst)))
log_msg("counts: ", nrow(cnt), " genes x ", ncol(cnt), " samples (vst agrees)")

# --- expression filter (RNA-seq): count >= 10 in >= 50% of samples
det_hcc <- rowSums(cnt >= 10) >= 0.5 * ncol(cnt)
log_msg(sprintf("RNA-seq filter: count>=10 in >=50%% -> %d/%d genes (%.1f%%)",
                sum(det_hcc), nrow(cnt), 100 * mean(det_hcc)))

E_hcc <- vst[det_hcc, , drop = FALSE]
hcc_coll <- collapse_by_mean(E_hcc, rownames(E_hcc))
E_hcc <- hcc_coll$expr
log_msg("HCC after collapse: ", nrow(E_hcc), " genes")

# ===========================================================================
# 3. common universe U
# ===========================================================================
U <- intersect(rownames(E_hf), rownames(E_hcc))
log_msg("=== universe U: ", length(U), " genes ===")

E_hf_U  <- E_hf[U, , drop = FALSE]
E_hcc_U <- E_hcc[U, , drop = FALSE]

# ===========================================================================
# 4. loss chain table
# ===========================================================================
loss <- rbind(
  data.frame(platform = "GPL11532 (HF array)", stage = c(
    "features on platform (matrix rows)",
    "detected by expression filter",
    "mapped to a gene symbol",
    "unique symbols after collapsing",
    "present in common universe U"),
    n = c(nrow(E_hf_raw), sum(det_hf), sum(!is.na(sym_hf_det)), nrow(E_hf), length(U))),
  data.frame(platform = "GSE141198 (HCC RNA-seq)", stage = c(
    "features on platform (symbol-mapped counts)",
    "detected by expression filter",
    "mapped to a gene symbol",
    "unique symbols after collapsing",
    "present in common universe U"),
    n = c(nrow(cnt), sum(det_hcc), nrow(E_hcc), nrow(E_hcc), length(U)))
)
print(loss)
write_table(loss, "Table_S1_gene_universe_loss_chain.csv")

# ===========================================================================
# 5. canonical-set coverage inside U  (the number the reviewers asked for)
# ===========================================================================
tiers <- list(TierA = SETS_TIER_A,
              TierB = c(SETS_TIER_A, SETS_TIER_B_EXTRA),
              TierC_v1 = SETS_TIER_C_V1)
cov_rows <- list()
for (tn in names(tiers)) {
  sets <- resolve_sets(tiers[[tn]])
  for (sn in names(sets)) {
    g <- sets[[sn]]
    cov_rows[[paste(tn, sn)]] <- data.frame(
      tier = tn, set = sn, n_genes_in_msigdb = length(g),
      n_in_HF_platform = sum(g %in% rownames(E_hf)),
      n_in_HCC_platform = sum(g %in% rownames(E_hcc)),
      n_in_common_U = sum(g %in% U),
      stringsAsFactors = FALSE)
  }
}
coverage <- do.call(rbind, cov_rows)
coverage$frac_in_U <- round(coverage$n_in_common_U / coverage$n_genes_in_msigdb, 3)
print(coverage[coverage$tier == "TierA", ], row.names = FALSE)
write_table(coverage, "Table_S2_canonical_set_coverage.csv")

# genes lost from Tier A specifically
A <- resolve_sets(SETS_TIER_A)
A_genes <- unique(unlist(A))
A_lost <- setdiff(A_genes, U)
log_msg("Tier A genes absent from U: ", length(A_lost), " / ", length(A_genes))

# Per-gene loss reason: which platform dropped it, and at which stage.
A_tab <- data.frame(
  gene            = A_genes,
  in_HF_platform  = A_genes %in% rownames(E_hf),
  in_HCC_platform = A_genes %in% rownames(E_hcc),
  in_common_U     = A_genes %in% U,
  stringsAsFactors = FALSE)
A_tab$loss_reason <- ifelse(
  A_tab$in_common_U, "retained",
  ifelse(!A_tab$in_HF_platform & !A_tab$in_HCC_platform, "absent from both platforms",
  ifelse(!A_tab$in_HF_platform, "lost on HF array (GPL11532)",
  ifelse(!A_tab$in_HCC_platform, "lost in HCC RNA-seq", "unexplained"))))
print(table(A_tab$loss_reason))
write_table(A_tab, "Table_S3_tierA_gene_loss.csv")
cat("\nTier A genes lost (", length(A_lost), "):\n", sep = "")
cat(paste(sort(A_lost), collapse = ", "), "\n\n")

# ===========================================================================
# 6. save
# ===========================================================================
saveRDS(list(expr = E_hf_U,  ann = ann_hf,  universe = U),
        file.path(OUT, "intermediate", "R01_HF.rds"))
saveRDS(list(expr = E_hcc_U, ann = data.frame(gsm = colnames(E_hcc_U)),
             universe = U),
        file.path(OUT, "intermediate", "R01_HCC.rds"))
saveRDS(list(E_hf_full = E_hf, E_hcc_full = E_hcc, U = U,
             ann_hf = ann_hf, loss = loss, coverage = coverage),
        file.path(OUT, "intermediate", "R01_universe.rds"))

log_msg("=== R01 done: |U| = ", length(U), " genes; ",
        "HF ", ncol(E_hf_U), " samples; HCC ", ncol(E_hcc_U), " samples ===")
