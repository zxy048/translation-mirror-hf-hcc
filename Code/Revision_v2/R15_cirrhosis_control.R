# =============================================================================
# R15_cirrhosis_control.R
# -----------------------------------------------------------------------------
# GSE89377 as a SAME-ORGAN, NON-MALIGNANT disease control for the HCC side.
#
# WHY THIS EXISTS. R2 Major 3 and R4 #1 both press the same point: the HCC arm
# of v1's mirror compares tumour against ADJACENT NORMAL LIVER, so "translation
# programs go up in HCC" may be a liver-disease effect -- inflammation, fibrosis,
# regenerative proliferation -- rather than anything cancer-specific. v1 had no
# way to separate the two, because every HCC cohort it used contrasts tumour with
# non-tumour tissue from the same organ.
#
# GSE89377 (GPL6947, Illumina HumanHT-12 V3) is a multistep hepatocarcinogenesis
# series with an explicit non-malignant arm, so it can separate them directly:
#
#   N  Normal (13) | FL chronic hepatitis, low grade (8) | FH chronic hepatitis,
#   high grade (12) | CS cirrhosis (12) | DL dysplastic nodule, low grade (11) |
#   DH dysplastic nodule, high grade (11) | eHCC early HCC (5) |
#   TG1/TG2/TG3 HCC (9/12/14)                                    total 107
#
# Cirrhosis is the decisive arm: cirrhotic liver is diseased but not malignant.
# If the Tier A translation sets move the same way in cirrhosis as in HCC, the
# HCC arm of the mirror is a liver-disease effect and the manuscript must say so.
# If they do not, the cancer-specific reading survives with actual evidence
# behind it instead of an assumption. Either answer is publishable; the
# unfalsifiable version is what the reviewers objected to.
#
# The ordered spectrum also lets the transition point be located rather than
# asserted: a monotone rise from Normal through cirrhosis to dysplasia to HCC is
# a different claim from a step change at the malignant transition.
#
# NORMALISATION. The series ships raw intensities with per-probe detection
# p-values. Probes are kept when detected (p < 0.05) in >= 50% of samples, then
# log2 + quantile normalised. This is a minimal, standard Illumina pipeline; no
# cohort-specific gene universe is imposed, because this cohort exists to test
# replication of a discovery effect (same policy as R05's validation cohorts).
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(GSVA))

log_msg("=== R15 start ===")

# ---------------------------------------------------------------------------
# 1. sample annotation
# ---------------------------------------------------------------------------
sm <- parse_series_matrix(file.path(VALID_DIR, "GSE89377_series_matrix.txt.gz"))
ph <- sm$ann$phenotype
gsm <- sm$ann$gsm
log_msg("  series matrix: ", length(gsm), " samples; phenotype levels:")
print(table(ph, useNA = "ifany"))

# The source data misspells "early" as "earyl"; match both rather than fix the
# data. Stage order is the biological progression, used for the trend column.
STAGES <- c("Normal", "CH_low", "CH_high", "Cirrhosis",
            "DN_low", "DN_high", "eHCC", "HCC")
stage_of <- function(p) {
  p <- tolower(p)
  ifelse(grepl("normal", p),                          "Normal",
  ifelse(grepl("hepatitis with low", p),              "CH_low",
  ifelse(grepl("hepatitis with high", p),             "CH_high",
  ifelse(grepl("cirrhosis", p),                       "Cirrhosis",
  ifelse(grepl("dysplastic nodules with low", p),     "DN_low",
  ifelse(grepl("dysplastic nodules with high", p),    "DN_high",
  ifelse(grepl("earyl|early hepatocell", p),          "eHCC",
  ifelse(grepl("hepatocellular carcinoma \\(tg", p),  "HCC",
         NA_character_))))))))
}
stage <- stage_of(ph)
if (anyNA(stage)) {
  log_msg("  WARNING: ", sum(is.na(stage)), " samples with unrecognised phenotype, dropped: ",
          paste(unique(ph[is.na(stage)]), collapse = " | "))
}
log_msg("  stage counts: ", paste(names(table(stage)), table(stage), sep = "=", collapse = " "))

# ---------------------------------------------------------------------------
# 2. raw intensities + detection p
# ---------------------------------------------------------------------------
con <- gzfile(file.path(VALID_DIR, "GSE89377_non-normalized.txt.gz"), "rt")
hdr <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
raw <- read.delim(con, header = FALSE, col.names = hdr, check.names = FALSE,
                  stringsAsFactors = FALSE)
close(con)
log_msg("  raw: ", nrow(raw), " probes x ", ncol(raw), " columns")

# columns alternate expression / detection p; the header names them SAMPLE k
idcol  <- hdr[1]
expr_c <- grep("^SAMPLE [0-9]+$", hdr)
det_c  <- grep("^Detection Pval", hdr)
stopifnot(length(expr_c) == length(gsm), length(det_c) == length(gsm),
          length(expr_c) == ncol(raw) - 1 - length(det_c))

# SAMPLE k maps to the k-th column of the series matrix. That mapping is the
# single assumption this whole script rests on -- get it wrong and every sample
# is labelled with someone else's disease stage, silently. So it is checked
# rather than assumed: the title prefixes (N-, FL-, CS-, TG1- ...) encode the
# same stage as the phenotype characteristic, and the two must agree sample for
# sample. A mismatch means the column order differs and the script must stop.
samp_no <- as.integer(sub("^SAMPLE ", "", hdr[expr_c]))
ord <- order(samp_no)
gm  <- gsm[ord]

tg <- grep("^!Sample_title\t", readLines(gzfile(file.path(VALID_DIR,
        "GSE89377_series_matrix.txt.gz")), warn = FALSE), value = TRUE)[1]
titles <- gsub('^"|"$', "", strsplit(tg, "\t", fixed = TRUE)[[1]][-1])[ord]
stopifnot(length(titles) == length(gm), all(nzchar(titles)))

stage_from_title <- function(t) {
  p <- toupper(sub("-.*$", "", t))
  ifelse(p == "N",    "Normal",
  ifelse(p == "FL",   "CH_low",
  ifelse(p == "FH",   "CH_high",
  ifelse(p == "CS",   "Cirrhosis",
  ifelse(p == "DL",   "DN_low",
  ifelse(p == "DH",   "DN_high",
  ifelse(p == "EHCC", "eHCC",
  ifelse(grepl("^TG", p), "HCC", NA_character_))))))))
}
st_title <- stage_from_title(titles)
st_pheno <- stage_of(ph)[ord]
n_ok <- sum(st_title == st_pheno)
log_msg("  title-vs-phenotype agreement: ", n_ok, " / ", length(gm))
if (n_ok != length(gm)) {
  bad <- which(st_title != st_pheno)
  log_msg("  MISMATCH at ", length(bad), " sample(s), first few: ",
          paste(sprintf("%s(title=%s,pheno=%s)", gm[bad], st_title[bad], st_pheno[bad])[seq_len(min(5, length(bad)))],
                collapse = " ; "))
  stop("GSE89377: SAMPLE k does not map to series-matrix column k")
}

E <- as.matrix(raw[, expr_c[ord], drop = FALSE])
P <- as.matrix(raw[, det_c[ord], drop = FALSE])
rownames(E) <- rownames(P) <- raw[[idcol]]
colnames(E) <- colnames(P) <- gm
storage.mode(E) <- "double"; storage.mode(P) <- "double"

# ---------------------------------------------------------------------------
# 3. detection filter -> log2 -> quantile normalise -> probe -> symbol
# ---------------------------------------------------------------------------
keep <- rowMeans(P < 0.05, na.rm = TRUE) >= 0.5
log_msg("  detection filter: ", sum(keep), " / ", length(keep), " probes kept")
E <- log2(E[keep, , drop = FALSE] + 1)

qn <- function(m) {
  r <- apply(m, 2, rank, ties.method = "average")
  ref <- sort(rowMeans(apply(m, 2, sort)))
  out <- apply(r, 2, function(x) ref[round(x)])
  dimnames(out) <- dimnames(m); out
}
E <- qn(E)
log_msg("  after log2 + quantile: range ", sprintf("%.2f-%.2f", min(E), max(E)))

ann <- read.delim(gzfile(file.path(VALID_DIR, "GPL6947.annot.gz")), skip = 28,
                  check.names = FALSE, quote = "", stringsAsFactors = FALSE)
log_msg("  GPL6947 annot: ", nrow(ann), " probes; cols: ",
        paste(head(names(ann), 4), collapse = ", "))
sym <- setNames(ann[["Gene symbol"]], ann[["ID"]])
p2s <- sym[rownames(E)]
col <- collapse_by_mean(E, unname(p2s))
log_msg("  ", nrow(col$expr), " genes x ", ncol(col$expr), " samples after collapse")

gx  <- col$expr
grp_stage <- setNames(stage[ord], colnames(gx))
if (anyNA(grp_stage)) {
  ok <- !is.na(grp_stage); gx <- gx[, ok, drop = FALSE]; grp_stage <- grp_stage[ok]
}
log_msg("  final: ", nrow(gx), " genes x ", ncol(gx), " samples; ",
        paste(names(table(grp_stage)), table(grp_stage), sep = "=", collapse = " "))

# ---------------------------------------------------------------------------
# 4. ssGSEA -- same parameter object as R06, so scores are comparable
# ---------------------------------------------------------------------------
sets_all <- resolve_sets(unique(c(SETS_TIER_A, SETS_TIER_B_EXTRA, SETS_TIER_C_V1)))
sets <- lapply(sets_all, function(g) intersect(g, rownames(gx)))
n_before <- vapply(sets, length, integer(1))
sc <- GSVA::gsva(GSVA::ssgseaParam(as.matrix(gx), sets[keep <- n_before >= 5],
                                   minSize = 5, maxSize = 500, normalize = TRUE),
                 verbose = FALSE)
log_msg("  ssGSEA: ", nrow(sc), " sets x ", ncol(sc), " samples; dropped ",
        sum(!keep), " sets (<5 genes measured)")
tier <- setNames(ifelse(rownames(sc) %in% SETS_TIER_A, "A",
                 ifelse(rownames(sc) %in% SETS_TIER_B_EXTRA, "B_extra",
                 ifelse(rownames(sc) %in% SETS_TIER_C_V1, "C_v1", NA_character_))),
                 rownames(sc))
log_msg("  tiers scored: ", paste(names(table(tier)), table(tier), sep = "=", collapse = " "))

# ---------------------------------------------------------------------------
# 5. effect size of every stage vs Normal, with R07's estimator
# ---------------------------------------------------------------------------
hedges_g <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 3 || n2 < 3) return(c(d = NA, g = NA, sp = NA, J = NA))
  s1 <- stats::sd(x); s2 <- stats::sd(y)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(c(d = NA, g = NA, sp = sp, J = NA))
  d <- (mean(x) - mean(y)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  c(d = d, g = J * d, sp = sp, J = J)
}
boot_g <- function(x, y, nboot = STATS$N_BOOT) {
  set.seed(STATS$SEED)
  v <- replicate(nboot, {
    xb <- sample(x, replace = TRUE); yb <- sample(y, replace = TRUE)
    hedges_g(xb, yb)[["g"]]
  })
  v <- v[is.finite(v)]
  if (length(v) < 100) return(c(NA, NA))
  stats::quantile(v, c(0.025, 0.975), names = FALSE)
}

ctl <- colnames(sc)[grp_stage[colnames(sc)] == "Normal"]
stopifnot(length(ctl) >= 5)
log_msg("  reference arm: Normal, n = ", length(ctl))

# One arm per non-reference stage; "HCC" pools TG1+TG2+TG3 (stage_of maps all
# three TG prefixes onto it, so this is the same arm the loop above builds).
arms <- lapply(setdiff(STAGES, "Normal"), function(s) {
  colnames(sc)[grp_stage[colnames(sc)] == s]
})
names(arms) <- setdiff(STAGES, "Normal")

rows <- list()
for (a in names(arms)) {
  d <- arms[[a]]
  if (length(d) < 3) { log_msg("  ", a, ": n = ", length(d), " < 3, skipped"); next }
  for (s in rownames(sc)) {
    x <- sc[s, d]; y <- sc[s, ctl]
    h  <- hedges_g(x, y)
    ci <- boot_g(x, y)
    tt <- tryCatch(stats::t.test(x, y), error = function(e) NULL)
    rows[[length(rows) + 1]] <- data.frame(
      cohort = "GSE89377", set = s, tier = tier[[s]],
      stage = a, stage_index = match(a, STAGES),
      n_stage = length(d), n_control = length(ctl),
      mean_stage = mean(x), mean_control = mean(y),
      hedges_g = h[["g"]], ci_lo = ci[1], ci_hi = ci[2],
      p_ttest = if (is.null(tt)) NA_real_ else tt$p.value,
      stringsAsFactors = FALSE)
  }
  log_msg("  ", a, " (", length(d), " vs ", length(ctl), ") done")
}
es <- do.call(rbind, rows)
es$ci_excludes_0   <- with(es, is.finite(ci_lo) & is.finite(ci_hi) & (ci_lo > 0 | ci_hi < 0))
es$meets_min_effect <- abs(es$hedges_g) >= STATS$MIN_EFFECT
es$p_fdr <- ave(es$p_ttest, es$stage, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
write_table(es, "Table_S34_GSE89377_stage_spectrum_effect_sizes.csv")

# ---------------------------------------------------------------------------
# 6. the decisive question: is the HCC direction cancer-specific?
# ---------------------------------------------------------------------------
wide <- reshape(es[, c("set", "tier", "stage", "hedges_g")],
                idvar = c("set", "tier"), timevar = "stage", direction = "wide")
names(wide) <- sub("^hedges_g\\.", "", names(wide))
stopifnot("Cirrhosis" %in% names(wide), "HCC" %in% names(wide))

# Direction agreement, per tier.
agree <- function(a, b, tier_sel = NULL) {
  d <- wide
  if (!is.null(tier_sel)) d <- d[d$tier %in% tier_sel, ]
  d <- d[is.finite(d[[a]]) & is.finite(d[[b]]), ]
  both_up   <- sum(d[[a]] > 0 & d[[b]] > 0)
  both_down <- sum(d[[a]] < 0 & d[[b]] < 0)
  disc      <- sum(sign(d[[a]]) != sign(d[[b]]))
  rho <- suppressWarnings(stats::cor(d[[a]], d[[b]], method = "spearman"))
  list(n = nrow(d), both_up = both_up, both_down = both_down, discordant = disc,
       concordant_pct = round(100 * (both_up + both_down) / nrow(d), 1),
       rho = round(rho, 3))
}
log_msg("--- cirrhosis vs HCC direction agreement ---")
for (tr in list("A", c("A", "B_extra", "C_v1"))) {
  lbl <- if (length(tr) == 1) paste0("Tier ", tr) else "ALL tiers"
  a <- agree("Cirrhosis", "HCC", tr)
  log_msg(sprintf("  %-10s n=%2d | both-up %2d | both-down %2d | discordant %2d | rho=%+.3f",
                  lbl, a$n, a$both_up, a$both_down, a$discordant, a$rho))
}

wide$cirrhosis_same_direction_as_HCC <- sign(wide$Cirrhosis) == sign(wide$HCC)
wide$cirrhosis_magnitude_ratio <- round(abs(wide$Cirrhosis) / pmax(abs(wide$HCC), 1e-9), 2)
write_table(wide, "Table_S35_GSE89377_HCC_vs_cirrhosis_specificity.csv")

# The same question asked of the malignant transition rather than of cancer vs
# normal liver: does the change happen AT the transition (dysplasia -> eHCC ->
# HCC) or is it already present in non-malignant disease?
log_msg("--- eHCC vs HCC direction agreement (is the change at the transition?) ---")
for (tr in list("A", c("A", "B_extra", "C_v1"))) {
  lbl <- if (length(tr) == 1) paste0("Tier ", tr) else "ALL tiers"
  a <- agree("eHCC", "HCC", tr)
  log_msg(sprintf("  %-10s n=%2d | both-up %2d | both-down %2d | discordant %2d | rho=%+.3f",
                  lbl, a$n, a$both_up, a$both_down, a$discordant, a$rho))
}

cat("\n=== Tier A: does cirrhosis move the same way as HCC? ===\n")
ta <- wide[wide$tier == "A", c("set", "CH_low", "CH_high", "Cirrhosis",
                               "DN_low", "DN_high", "eHCC", "HCC",
                               "cirrhosis_same_direction_as_HCC")]
ta <- ta[order(-abs(ta$HCC)), ]
print(ta, row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# 7. cross-check: does GSE89377's HCC arm reproduce the other HCC cohorts?
# ---------------------------------------------------------------------------
r7 <- readRDS(file.path(OUT, "intermediate", "R07_effect_sizes.rds"))
hcc7 <- r7[r7$disease == "HCC" & r7$cohort != "GSE141198", c("cohort", "set", "hedges_g")]
wide7 <- reshape(hcc7, idvar = "set", timevar = "cohort", direction = "wide")
names(wide7) <- sub("^hedges_g\\.", "", names(wide7))
me <- merge(wide[wide$tier == "A", c("set", "HCC")], wide7, by = "set")
log_msg("--- GSE89377 HCC arm vs the HCC cohorts v1 already used (Tier A, n=", nrow(me), ") ---")
for (cc in setdiff(names(me), c("set", "HCC"))) {
  v <- me[[cc]]
  ok <- is.finite(v) & is.finite(me$HCC)
  if (sum(ok) < 5) next
  log_msg(sprintf("  vs %-14s n=%2d | sign agreement %.0f%% | Spearman rho=%.3f",
                  cc, sum(ok), 100 * mean(sign(v[ok]) == sign(me$HCC[ok])),
                  suppressWarnings(stats::cor(me$HCC[ok], v[ok], method = "spearman"))))
}

saveRDS(list(es = es, wide = wide, scores = sc, stage = grp_stage,
             crosscheck = me, stages = STAGES),
        file.path(OUT, "intermediate", "R15_cirrhosis_control.rds"))
log_msg("=== R15 done ===")
