# =============================================================================
# R05_load_cohorts.R
# -----------------------------------------------------------------------------
# Load all seven cohorts into a common shape: genes x samples in SYMBOL space,
# with a two-level group label (disease / control) so every cohort can yield a
# within-disease disease-vs-control effect size.
#
# WHY within-disease effect sizes: v1's "mirror" compared Cohens_d_HCC (TCGA-LIHC
# tumour vs normal) against Cohens_d_HF (GSE57338 HF vs non-failing), computed in
# 03_ssGSEA_cross_disease_pathway.R with an EQUAL-WEIGHT pooled SD
# sqrt((s1^2+s2^2)/2). R07 recomputes the same quantity with an n-weighted pooled
# SD + Hedges correction -- the same estimand, correctly estimated.
#
# FILTERING POLICY: discovery cohorts keep R01's declared filter and the common
# universe U. Validation cohorts get only a MINIMAL detection filter, because
# their job is to test replication of a discovery effect; imposing a different
# gene universe on them would confound replication with gene-set change. This
# asymmetry is deliberate and is reported per cohort (n_genes, Tier A coverage).
#
# DISCOVERY (gene universe U):
#   GSE57338  HF  array  GPL11532     177 HF / 136 NF      [from R01]
#   GSE141198 HCC RNA-seq             148 tumour / 0 ctrl  [from R01]
#   -> GSE141198 has NO control arm; it is a network-discovery cohort only and
#      yields no within-disease effect size. HCC effects come from TCGA/GSE14520/
#      GSE76427 (the same cohorts v1 used for the HCC side).
# VALIDATION:
#   TCGA_LIHC HCC RNA-seq             371 tumour / 50 normal
#   GSE14520  HCC array  GPL3921+GPL571 (2 sub-cohorts, analysed separately)
#   GSE76427  HCC array  GPL10558
#   GSE116250 HF  RNA-seq             DCM+ICM / non-failing
#   GSE141910 HF  RNA-seq             HF / non-failing donor
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages({library(SummarizedExperiment)})

log_msg("=== R05 start ===")
cohorts <- list()

# --- helper: assemble the standard record ------------------------------------
mk <- function(expr, group, cohort, disease, tech, note = "",
               ann = NULL, subset_label = NA_character_) {
  stopifnot(ncol(expr) == length(group))
  if (is.null(rownames(expr)) || anyDuplicated(rownames(expr)))
    stop("rownames missing or duplicated for ", cohort)
  # Samples whose group could not be determined are DROPPED, explicitly and
  # with a log line -- never silently carried as an NA level into a group test.
  bad <- is.na(group) | group == ""
  if (any(bad)) {
    log_msg("  ", cohort, ": dropping ", sum(bad),
            " sample(s) with undetermined group")
    expr <- expr[, !bad, drop = FALSE]; group <- group[!bad]
  }
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  # group is ALWAYS named by sample ID. Callers that build it with ifelse()
  # produce an unnamed vector, and every downstream join keys on these names.
  names(group) <- colnames(expr)
  structure(list(expr = expr, group = group, cohort = cohort, disease = disease,
                 tech = tech, note = note, ann = ann,
                 subset_label = subset_label),
            class = "cohort")
}

# =============================================================================
# 1. HF discovery -- GSE57338 (reuse R01 exactly)
# =============================================================================
log_msg("--- GSE57338 (HF array, discovery) ---")
r1 <- readRDS(file.path(OUT, "intermediate", "R01_universe.rds"))
ann_hf <- r1$ann_hf
grp_hf <- ifelse(ann_hf$group == "HF", "disease",
          ifelse(ann_hf$group == "NF", "control", NA_character_))
names(grp_hf) <- ann_hf$gsm
E_hf <- r1$E_hf_full[, ann_hf$gsm, drop = FALSE]
cohorts$GSE57338 <- mk(E_hf, grp_hf[colnames(E_hf)], "GSE57338", "HF", "array",
                       note = "discovery; gene universe U",
                       ann = ann_hf, subset_label = "all")
log_msg("  ", nrow(E_hf), " genes x ", ncol(E_hf), " samples; ",
        paste(names(table(grp_hf, useNA = "ifany")), table(grp_hf, useNA = "ifany"),
              collapse = " "))

# =============================================================================
# 2. HCC discovery -- GSE141198 (tumour only)
# =============================================================================
log_msg("--- GSE141198 (HCC RNA-seq, discovery) ---")
E_hcc <- r1$E_hcc_full
grp_hcc <- setNames(rep("disease", ncol(E_hcc)), colnames(E_hcc))
cohorts$GSE141198 <- mk(E_hcc, grp_hcc, "GSE141198", "HCC", "RNA-seq",
                        note = "discovery; tumour only, no control arm -> no effect size",
                        subset_label = "all")
log_msg("  ", nrow(E_hcc), " genes x ", ncol(E_hcc), " samples (all tumour)")

# =============================================================================
# 3. TCGA-LIHC (HCC RNA-seq validation)
# =============================================================================
log_msg("--- TCGA-LIHC (HCC RNA-seq) ---")
se <- readRDS(file.path(DATA_RAW, "TCGA_LIHC_se.rds"))
cd <- as.data.frame(colData(se))
cnt <- assay(se, "unstranded")
sym <- rowData(se)$gene_name
keep_type <- cd$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
cnt <- cnt[, keep_type, drop = FALSE]
grp <- ifelse(cd$sample_type[keep_type] == "Primary Tumor", "disease", "control")
log_msg("  raw: ", nrow(cnt), " features x ", ncol(cnt), " samples; ",
        paste(names(table(grp)), table(grp), collapse = " "))
cnt_det <- cnt[filter_counts(cnt, 10), , drop = FALSE]
col <- collapse_by_mean(cnt_det, sym[filter_counts(cnt, 10)])
log_msg("  after count filter + collapse: ", nrow(col$expr), " genes")
cohorts$TCGA_LIHC <- mk(log2(col$expr + 1), grp, "TCGA_LIHC", "HCC", "RNA-seq",
                        note = "primary tumour vs solid tissue normal; log2(count+1)",
                        subset_label = "Primary Tumor / Solid Tissue Normal")

# =============================================================================
# 4. GSE116250 (HF RNA-seq validation) -- rpkm with an explicit symbol column
# =============================================================================
log_msg("--- GSE116250 (HF RNA-seq) ---")
rp <- read.delim(gzfile(file.path(VALID_DIR, "GSE116250_rpkm.txt.gz")),
                 check.names = FALSE)
sym_col <- rp$Common_name
expr <- as.matrix(rp[, setdiff(names(rp), c("Gene", "Common_name")), drop = FALSE])
sm116 <- parse_series_matrix(file.path(VALID_DIR, "GSE116250_series_matrix.txt.gz"))
# rpkm column names are sample titles; map title -> disease state via the matrix
ti <- setNames(sm116$ann$disease, sm116$ann$gsm)
# !Sample_title holds the same strings as the rpkm header ("NF1", "DCM8", ...)
con <- gzfile(file.path(VALID_DIR, "GSE116250_series_matrix.txt.gz"), "rt")
l <- readLines(con, warn = FALSE); close(con)
tl <- grep("^!Sample_title\t", l, value = TRUE)[1]
titles <- gsub('^"|"$', "", strsplit(tl, "\t", fixed = TRUE)[[1]][-1])
disease_state <- setNames(sm116$ann$disease, titles)
st <- disease_state[colnames(expr)]
st[is.na(st)] <- "non-failing"
grp116 <- ifelse(grepl("non-?failing|normal", st, ignore.case = TRUE), "control", "disease")
expr116 <- expr[filter_nonzero(expr, 0.5), , drop = FALSE]
col116 <- collapse_by_mean(expr116, sym_col[filter_nonzero(expr, 0.5)])
log_msg("  states: ", paste(names(table(st)), table(st), collapse = " | "))
log_msg("  after filter + collapse: ", nrow(col116$expr), " genes")
cohorts$GSE116250 <- mk(col116$expr, grp116, "GSE116250", "HF", "RNA-seq",
                        note = "DCM+ICM vs non-failing; rpkm",
                        subset_label = "all")

# =============================================================================
# 5. GSE14520 (HCC array validation) -- TWO platforms, kept separate
# =============================================================================
log_msg("--- GSE14520 (HCC array, GPL3921 + GPL571) ---")
for (gpl in c("GPL3921", "GPL571")) {
  key <- paste0("GSE14520_", gpl)
  sm <- parse_series_matrix(file.path(DATA_RAW,
          sprintf("GSE14520-%s_series_matrix.txt.gz", gpl)))
  soft <- parse_soft(file.path(DATA_RAW, sprintf("%s.soft.gz", gpl)))
  p2s <- soft_probe2symbol(soft)
  tissue <- sm$ann$Tissue
  grp <- ifelse(grepl("Non-Tumor", tissue, ignore.case = TRUE), "control",
         ifelse(grepl("Tumor", tissue, ignore.case = TRUE), "disease", NA_character_))
  det <- filter_array(sm$expr)
  col <- collapse_by_mean(sm$expr[det, , drop = FALSE], p2s[rownames(sm$expr)[det]])
  log_msg("  ", gpl, ": ", nrow(col$expr), " genes x ", ncol(col$expr),
          " samples; ", paste(names(table(grp, useNA = "ifany")),
                              table(grp, useNA = "ifany"), collapse = " "))
  cohorts[[key]] <- mk(col$expr, grp, key, "HCC", "array",
                       note = paste("tumour vs adjacent non-tumour;", gpl),
                       ann = sm$ann, subset_label = gpl)
}
rm(soft, p2s)

# =============================================================================
# 6. GSE76427 (HCC array validation, GPL10558)
# =============================================================================
log_msg("--- GSE76427 (HCC array, GPL10558) ---")
sm764 <- parse_series_matrix(file.path(DATA_RAW, "GSE76427_series_matrix.txt.gz"))
soft764 <- parse_soft(file.path(DATA_RAW, "GPL10558.soft.gz"))
p2s764 <- soft_probe2symbol(soft764)
tis <- sm764$ann$tissue
grp764 <- ifelse(grepl("adjacent|non-?tumou?r|normal", tis, ignore.case = TRUE), "control",
          ifelse(grepl("tumor|tumour|carcinoma", tis, ignore.case = TRUE), "disease",
                 NA_character_))
det764 <- filter_array(sm764$expr)
col764 <- collapse_by_mean(sm764$expr[det764, , drop = FALSE],
                           p2s764[rownames(sm764$expr)[det764]])
log_msg("  ", nrow(col764$expr), " genes x ", ncol(col764$expr), " samples; ",
        paste(names(table(grp764, useNA = "ifany")),
              table(grp764, useNA = "ifany"), collapse = " "))
cohorts$GSE76427 <- mk(col764$expr, grp764, "GSE76427", "HCC", "array",
                       note = "tumour vs adjacent normal; GPL10558",
                       ann = sm764$ann, subset_label = "all")

# =============================================================================
# 7. GSE141910 (HF RNA-seq validation) -- 366 per-sample CSVs
# =============================================================================
log_msg("--- GSE141910 (HF RNA-seq, per-sample CSVs) ---")
sm141 <- parse_series_matrix(file.path(VALID_DIR, "GSE141910_series_matrix.txt.gz"))
ety <- sm141$ann$etiology
grp141 <- ifelse(grepl("non-?failing|donor|normal", ety, ignore.case = TRUE), "control",
          ifelse(is.na(ety) | ety == "", NA_character_, "disease"))
names(grp141) <- sm141$ann$gsm
files <- list.files(file.path(VALID_DIR, "GSE141910"), full.names = TRUE)
log_msg("  per-sample files: ", length(files))
lst <- lapply(files, function(f) {
  d <- read.csv(gzfile(f), row.names = 1, check.names = FALSE)
  setNames(d[[1]], rownames(d))
})
# filename is GSM<id>_<title>.csv.gz -- key on the GSM (which is what the series
# matrix annotates), not on the title.
gsm141 <- sub("_.*$", "", basename(files))
names(lst) <- gsm141
allg <- Reduce(union, lapply(lst, names))
m141 <- vapply(lst, function(v) { out <- rep(NA_real_, length(allg))
                                  out[match(names(v), allg)] <- v; out },
               numeric(length(allg)))
rownames(m141) <- allg
ok <- gsm141 %in% sm141$ann$gsm
log_msg("  GSM match to series matrix: ", sum(ok), "/", length(gsm141))
m141 <- m141[, ok, drop = FALSE]; gsm141 <- gsm141[ok]
sym141 <- ensg2symbol(rownames(m141))
det141 <- filter_nonzero(m141, 0.5)
col141 <- collapse_by_mean(m141[det141, , drop = FALSE], sym141[det141])
grp141f <- grp141[colnames(col141$expr)]
log_msg("  ", nrow(col141$expr), " genes x ", ncol(col141$expr), " samples; ",
        paste(names(table(grp141f, useNA = "ifany")),
              table(grp141f, useNA = "ifany"), collapse = " "))
cohorts$GSE141910 <- mk(col141$expr, grp141f, "GSE141910", "HF", "RNA-seq",
                        note = "HF vs non-failing donor; per-sample CSV, ENSG->symbol",
                        ann = sm141$ann, subset_label = "all")

# =============================================================================
# 8. Coverage report -- how much of each canonical tier survives per cohort
# =============================================================================
log_msg("--- coverage per cohort ---")
A <- unique(unlist(resolve_sets(SETS_TIER_A)))
cov <- do.call(rbind, lapply(names(cohorts), function(nm) {
  co <- cohorts[[nm]]
  g <- rownames(co$expr)
  data.frame(cohort = nm, disease = co$disease, tech = co$tech,
             n_genes = length(g),
             n_samples = ncol(co$expr),
             n_disease = sum(co$group == "disease", na.rm = TRUE),
             n_control = sum(co$group == "control", na.rm = TRUE),
             n_unlabeled = sum(is.na(co$group)),
             tierA_covered = length(intersect(g, A)),
             tierA_total = length(A),
             stringsAsFactors = FALSE)
}))
cov$tierA_pct <- round(100 * cov$tierA_covered / cov$tierA_total, 1)
print(cov, row.names = FALSE)
write_table(cov, "Table_S11_cohort_coverage.csv")

if (any(cov$n_control == 0 & cov$n_unlabeled == 0))
  log_msg("NOTE: cohorts with no control arm cannot yield an effect size: ",
          paste(cov$cohort[cov$n_control == 0 & cov$n_unlabeled == 0], collapse = ", "))

saveRDS(cohorts, file.path(OUT, "intermediate", "R05_cohorts.rds"))
log_msg("=== R05 done: ", length(cohorts), " cohorts ===")
