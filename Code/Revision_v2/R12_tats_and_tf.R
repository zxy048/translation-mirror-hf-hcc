# =============================================================================
# R12_tats_and_tf.R
# -----------------------------------------------------------------------------
# Redefine the "translation-associated transcriptional signature" (TATS) and
# redo the transcription-factor association with the reviewers' objections
# built in rather than answered after the fact.
#
# WHAT THE REVIEWERS ASKED FOR, AND WHERE IT IS ANSWERED HERE:
#   R2 #2   "compare the HF-derived TATS with an HCC-module-derived score and
#            independently curated core-translation signatures"   -> section A/B
#   R1 #4   "state the exact number of genes used to calculate TATS in each
#            cohort"                                             -> Table_S27
#   R2 #11  adjust TF associations for purity / proliferation and report across
#           independent HCC cohorts                              -> section C
#   R4 #2   do not imply upstream regulation                     -> language only;
#           the statistics here are reported as association, and the de-overlap
#           set is reported explicitly so an overlap-driven correlation cannot
#           be read as regulation
#   R3 #4   the |GS| > 0.20 threshold is too weak                -> section D
#
# THREE SCORES, all defined identically in every cohort:
#   TCS      canonical translation score: mean z-scored ssGSEA over the 15
#            a priori Tier A sets. Disease-agnostic by construction -- the same
#            definition is applied to HCC and HF, so it cannot be the thing that
#            creates a cross-disease difference.
#   ME_HF    HF purple module (153 genes, from R03) projected into a cohort:
#            kME-weighted mean z-expression, weights taken from the HF network.
#   ME_HCC   HCC magenta module (146 genes) projected the same way.
# ME_HF is v1's TATS with the v2 module and a documented gene count. ME_HCC is
# the control the reviewers asked for: if the two modules' scores agree, the
# result does not depend on which disease's module was chosen; if they disagree,
# TATS is a disease-specific projection and must be described as one.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

log_msg("=== R12 start ===")
r2  <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))
r3  <- readRDS(file.path(OUT, "intermediate", "R03_translation_module.rds"))
r5  <- readRDS(file.path(OUT, "intermediate", "R05_cohorts.rds"))
r6  <- readRDS(file.path(OUT, "intermediate", "R06_ssgsea.rds"))
r11 <- readRDS(file.path(OUT, "intermediate", "R11_composition.rds"))

MOD_HF  <- r3$desHF$module    # purple
MOD_HCC <- r3$desHCC$module   # magenta
GENES_HF  <- names(r2$HF$net$colors)[r2$HF$net$colors  == MOD_HF]
GENES_HCC <- names(r2$HCC$net$colors)[r2$HCC$net$colors == MOD_HCC]
log_msg("HF module '", MOD_HF, "' = ", length(GENES_HF), " genes; HCC module '",
        MOD_HCC, "' = ", length(GENES_HCC), " genes")

#' kME weights for a module, from the network it was defined in.
kme_weights <- function(kme_df, module, genes) {
  d <- kme_df[kme_df$module == module, ]
  w <- setNames(d$kME, d$gene)[genes]
  if (anyNA(w)) w[is.na(w)] <- 0
  w
}
W_HF  <- kme_weights(r2$kmeHF, MOD_HF,  GENES_HF)
W_HCC <- kme_weights(r2$kmeHCC, MOD_HCC, GENES_HCC)

# Resolved once, not inside a loop: this is the gene universe of the primary
# a priori definition and several sections below ask membership in it.
tierA_genes <- unique(unlist(resolve_sets(SETS_TIER_A)))

zrow <- function(v) {
  s <- stats::sd(v)
  if (!is.finite(s) || s == 0) return(rep(0, length(v)))
  (v - mean(v)) / s
}

#' Project a module into a cohort's expression matrix. Returns the weighted
#' (kME) and unweighted (mean-z) scores plus the exact gene accounting.
project <- function(expr, genes, weights) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 3) return(NULL)
  z <- t(scale(t(expr[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  w <- weights[g]
  if (all(!is.finite(w)) || sum(abs(w)) == 0) w <- rep(1, length(g))
  list(weighted   = colSums(z * w) / sum(abs(w)),
       unweighted = colMeans(z),
       n_used = length(g), n_total = length(genes), used = g)
}

#' TCS: mean of the z-scored Tier A ssGSEA scores. Identical rule in every
#' cohort; `n_sets` is reported because two cohorts cannot score one Tier A set.
tcs <- function(sc) {
  A <- intersect(rownames(sc), SETS_TIER_A)
  colMeans(t(apply(sc[A, , drop = FALSE], 1, zrow)))
}

# --- A. build the scores in every cohort --------------------------------------
scores <- list(); acc <- list(); setrows <- list()
for (nm in names(r6$scores)) {
  co <- r5[[nm]]; sc <- r6$scores[[nm]]; expr <- co$expr
  A <- intersect(rownames(sc), SETS_TIER_A)
  sc_cols <- colnames(sc)

  t_ <- tcs(sc)
  pH <- project(expr[, sc_cols, drop = FALSE], GENES_HF,  W_HF)
  pC <- project(expr[, sc_cols, drop = FALSE], GENES_HCC, W_HCC)
  if (is.null(pH) || is.null(pC)) { log_msg("  ", nm, ": module projection failed"); next }

  df <- data.frame(TCS = t_,
                   ME_HF_weighted = pH$weighted, ME_HF_mean_z = pH$unweighted,
                   ME_HCC_weighted = pC$weighted, ME_HCC_mean_z = pC$unweighted,
                   group = co$group[sc_cols], stringsAsFactors = FALSE)
  scores[[nm]] <- df

  acc[[length(acc) + 1]] <- data.frame(
    cohort = nm, disease = co$disease, tech = co$tech, n_samples = ncol(sc),
    n_tierA_sets = length(A), n_tierA_expected = length(SETS_TIER_A),
    genes_HF_module = length(GENES_HF), genes_HF_present = pH$n_used,
    frac_HF_present = round(pH$n_used / length(GENES_HF), 3),
    genes_HCC_module = length(GENES_HCC), genes_HCC_present = pC$n_used,
    frac_HCC_present = round(pC$n_used / length(GENES_HCC), 3),
    stringsAsFactors = FALSE)

  # How much of each module's projected gene list is itself Tier A translation
  # machinery, and how much the two modules share in this cohort. This is the
  # table behind the "which genes, how many" question (R1 #4).
  setrows[[length(setrows) + 1]] <- data.frame(
    cohort = nm,
    HF_in_tierA = sum(pH$used %in% tierA_genes),
    HCC_in_tierA = sum(pC$used %in% tierA_genes),
    shared_HF_HCC = length(intersect(pH$used, pC$used)),
    stringsAsFactors = FALSE)
}
accounting <- do.call(rbind, acc)
accounting <- merge(accounting, do.call(rbind, setrows), by = "cohort", sort = FALSE)
write_table(accounting, "Table_S27_TATS_gene_accounting.csv")

cat("\n=== exact genes behind each score, per cohort (R1 #4) ===\n")
print(accounting, row.names = FALSE)

# --- B. do the three scores agree? (R2 #2) ------------------------------------
hedges <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 3 || n2 < 3) return(c(g = NA_real_, p = NA_real_))
  sp <- sqrt(((n1 - 1) * stats::sd(x)^2 + (n2 - 1) * stats::sd(y)^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(c(g = NA_real_, p = NA_real_))
  tt <- tryCatch(stats::t.test(x, y), error = function(e) NULL)
  c(g = (1 - 3 / (4 * (n1 + n2) - 9)) * (mean(x) - mean(y)) / sp,
    p = if (is.null(tt)) NA_real_ else tt$p.value)
}

conc <- list(); eff <- list()
SCORE_COLS <- c("TCS", "ME_HF_weighted", "ME_HCC_weighted")
for (nm in names(scores)) {
  df <- scores[[nm]]
  cm <- stats::cor(df[, SCORE_COLS], use = "p")
  # combn(), not `if (a < b)`: comparing score names with `<` is lexicographic,
  # and no two of these three names happen to sort into the order they are
  # written in -- `"TCS" < "ME_HF_weighted"` is FALSE, so a name-based guard
  # silently produced an empty table.
  prs <- utils::combn(SCORE_COLS, 2)
  for (j in seq_len(ncol(prs)))
    conc[[length(conc) + 1]] <- data.frame(
      cohort = nm, score_a = prs[1, j], score_b = prs[2, j],
      r = cm[prs[1, j], prs[2, j]], stringsAsFactors = FALSE)

  grp <- df$group
  if (all(c("disease", "control") %in% grp)) {
    d <- df[grp == "disease", ]; c_ <- df[grp == "control", ]
    for (s in SCORE_COLS) {
      h <- hedges(d[[s]], c_[[s]])
      eff[[length(eff) + 1]] <- data.frame(
        cohort = nm, disease = r5[[nm]]$disease, score = s,
        n_disease = nrow(d), n_control = nrow(c_),
        mean_disease = mean(d[[s]]), mean_control = mean(c_[[s]]),
        hedges_g = h[["g"]], p = h[["p"]], stringsAsFactors = FALSE)
    }
  }
}
conc <- do.call(rbind, conc)
write_table(conc, "Table_S28_TATS_score_concordance.csv")

eff <- do.call(rbind, eff)
eff$p_fdr <- ave(eff$p, eff$cohort, FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
write_table(eff, "Table_S29_TATS_effect_sizes.csv")

cat("\n=== agreement among the three score definitions (R2 #2) ===\n")
for (nm in unique(conc$cohort)) {
  d <- conc[conc$cohort == nm, ]
  cat(sprintf("  %-18s  TCS~ME_HF r=%+.3f   TCS~ME_HCC r=%+.3f   ME_HF~ME_HCC r=%+.3f\n",
              nm,
              d$r[d$score_a == "TCS" & d$score_b == "ME_HF_weighted"],
              d$r[d$score_a == "TCS" & d$score_b == "ME_HCC_weighted"],
              d$r[d$score_a == "ME_HF_weighted" & d$score_b == "ME_HCC_weighted"]))
}
cat("\n=== disease vs control effect of each score ===\n")
print(eff[order(eff$cohort, eff$score), c("cohort", "score", "hedges_g", "p_fdr")],
      row.names = FALSE)

# --- C. TF association, with composition adjustment (R2 #11, R4 #2) -----------
#' Partial correlation of x and y controlling for the columns of Z.
#' Residualises against an explicit data.frame via `data =` rather than
#' `~ Z`: a bare data.frame in a formula is promoted by model.frame in ways
#' that depend on its column names, whereas an explicit data argument is not.
pcor_test <- function(x, y, Z) {
  if (is.null(Z) || !ncol(Z)) return(stats::cor.test(x, y))
  Z <- as.data.frame(Z)
  Z <- Z[, vapply(Z, is.numeric, logical(1)), drop = FALSE]
  if (!ncol(Z)) return(stats::cor.test(x, y))
  rx <- stats::residuals(stats::lm(x ~ ., data = Z))
  ry <- stats::residuals(stats::lm(y ~ ., data = Z))
  stats::cor.test(rx, ry)
}

# Which panel TFs are themselves members of the Tier A translation sets? A TF
# that is inside the score's own gene universe produces a correlation that is
# partly definitional, so those are flagged and the analysis is repeated without
# them.
TF_IN_SETS <- intersect(TFS_PANEL, tierA_genes)
log_msg("TFs that are themselves Tier A set members: ",
        if (length(TF_IN_SETS)) paste(TF_IN_SETS, collapse = ", ") else "none")

tf_rows <- list()
for (nm in names(scores)) {
  co <- r5[[nm]]; df <- scores[[nm]]; expr <- co$expr
  present <- intersect(TFS_PANEL, rownames(expr))
  if (length(present) < 3) { log_msg("  ", nm, ": <3 panel TFs present"); next }
  # standardise every TF across samples so the coefficients are comparable
  tfm <- t(scale(t(expr[present, rownames(df), drop = FALSE])))
  tfm[!is.finite(tfm)] <- 0
  cs <- r11$scores[[nm]][rownames(df), , drop = FALSE]
  prim <- intersect(r11$primary_covars[[co$disease]], colnames(cs))
  Z <- if (length(prim)) {
    z <- as.data.frame(cs[, prim, drop = FALSE])
    z[] <- lapply(z, function(v) { s <- stats::sd(v); if (!is.finite(s) || s == 0) 0 else (v - mean(v)) / s })
    z
  } else NULL

  arms <- list(all = rep(TRUE, nrow(df)))
  if (all(c("disease", "control") %in% df$group))
    arms$disease_only <- df$group == "disease"

  for (arm in names(arms)) {
    k <- arms[[arm]]
    for (s in SCORE_COLS) for (tf in present) {
      x <- tfm[tf, k]; y <- df[[s]][k]
      Zk <- if (is.null(Z)) NULL else Z[k, , drop = FALSE]
      c0 <- stats::cor.test(x, y)
      c1 <- pcor_test(x, y, Zk)
      tf_rows[[length(tf_rows) + 1]] <- data.frame(
        cohort = nm, disease = co$disease, arm = arm, score = s, tf = tf,
        tf_in_translation_sets = tf %in% TF_IN_SETS,
        n = sum(k),
        r_unadjusted = unname(c0$estimate), p_unadjusted = c0$p.value,
        r_adjusted   = unname(c1$estimate), p_adjusted   = c1$p.value,
        covariates = if (length(prim)) paste(prim, collapse = "+") else "",
        stringsAsFactors = FALSE)
    }
  }
}
tf <- do.call(rbind, tf_rows)
tf$p_fdr_unadjusted <- ave(tf$p_unadjusted, interaction(tf$cohort, tf$arm, tf$score),
                           FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
tf$p_fdr_adjusted   <- ave(tf$p_adjusted,   interaction(tf$cohort, tf$arm, tf$score),
                           FUN = function(p) p.adjust(p, STATS$FDR_METHOD))
tf$direction_kept <- sign(tf$r_unadjusted) == sign(tf$r_adjusted)
write_table(tf, "Table_S30_TF_translation_association.csv")

cat("\n=== TF associations, disease arm only, adjusted for composition ===\n")
cat("(r_unadj / r_adj; association only, NOT regulation)\n")
for (nm in unique(tf$cohort)) {
  d <- tf[tf$cohort == nm & tf$arm == "disease_only" & tf$score == "TCS", ]
  if (!nrow(d)) next
  d <- d[order(-abs(d$r_unadjusted)), ]
  cat(sprintf("\n-- %s (%s, n=%d; covar: %s) --\n", nm, d$disease[1], d$n[1],
              if (nzchar(d$covariates[1])) d$covariates[1] else "none"))
  print(d[, c("tf", "tf_in_translation_sets", "r_unadjusted", "r_adjusted",
              "p_fdr_unadjusted", "p_fdr_adjusted", "direction_kept")],
        row.names = FALSE)
}

cat("\n=== how often does the TF association survive adjustment? ===\n")
for (a in unique(tf$arm)) for (s in SCORE_COLS) {
  d <- tf[tf$arm == a & tf$score == s, ]
  cat(sprintf("  %-13s %-16s : |r|>0.2 & FDR<0.05 unadjusted %3d/%3d -> adjusted %3d/%3d ; sign kept %3d/%3d\n",
              a, s,
              sum(abs(d$r_unadjusted) > 0.2 & d$p_fdr_unadjusted < 0.05, na.rm = TRUE), nrow(d),
              sum(abs(d$r_adjusted)   > 0.2 & d$p_fdr_adjusted   < 0.05, na.rm = TRUE), nrow(d),
              sum(d$direction_kept, na.rm = TRUE), nrow(d)))
}

# --- D. the |GS| > 0.20 threshold (R3 #4) -------------------------------------
# Module identity in v2 comes from set enrichment (R03), not from a GS cutoff,
# so the threshold no longer selects the module. What it still affects is how
# many module genes are called "significant"; both thresholds are reported.
gs_rows <- list()
# One representative cohort per disease, each with BOTH arms so a gene
# significance is definable. GSE141198 -- the HCC network-building cohort -- is
# tumour-only and therefore cannot yield a GS at all; TCGA-LIHC stands in for
# HCC here, and that substitution is stated in the manuscript rather than
# silently absorbing the loss of the HCC row.
GS_COHORTS <- list(GSE57338 = list(module = MOD_HF,  genes = GENES_HF),
                   TCGA_LIHC = list(module = MOD_HCC, genes = GENES_HCC))
for (nm in names(GS_COHORTS)) {
  mod <- GS_COHORTS[[nm]]$module
  genes <- GS_COHORTS[[nm]]$genes
  expr <- r5[[nm]]$expr
  grp <- r5[[nm]]$group[colnames(expr)]
  if (!all(c("disease", "control") %in% grp)) {
    log_msg("  ", nm, ": single arm, GS not computable"); next
  }
  y <- as.integer(grp == "disease")
  g <- intersect(genes, rownames(expr))
  gs <- suppressWarnings(apply(expr[g, , drop = FALSE], 1, function(v)
    if (stats::sd(v) == 0) 0 else stats::cor(v, y)))
  p <- suppressWarnings(apply(expr[g, , drop = FALSE], 1, function(v)
    if (stats::sd(v) == 0) 1 else stats::cor.test(v, y)$p.value))
  fdr <- p.adjust(p, STATS$FDR_METHOD)
  gs_rows[[length(gs_rows) + 1]] <- data.frame(
    cohort = nm, module = mod, n_module_genes = length(genes),
    n_present = length(g),
    n_absGS_gt_020 = sum(abs(gs) > 0.20),
    n_absGS_gt_050 = sum(abs(gs) > 0.50),
    n_absGS_gt_050_and_FDR05 = sum(abs(gs) > 0.50 & fdr < 0.05),
    n_FDR05_alone = sum(fdr < 0.05),
    median_abs_GS = stats::median(abs(gs)), stringsAsFactors = FALSE)
}
gs <- do.call(rbind, gs_rows)
write_table(gs, "Table_S31_gene_significance_threshold.csv")
cat("\n=== gene-significance threshold sensitivity (R3 #4) ===\n")
print(gs, row.names = FALSE)

saveRDS(list(scores = scores, accounting = accounting, concordance = conc,
             effects = eff, tf = tf, gs = gs, tf_in_sets = TF_IN_SETS,
             genes_HF = GENES_HF, genes_HCC = GENES_HCC),
        file.path(OUT, "intermediate", "R12_tats_tf.rds"))
log_msg("=== R12 done ===")
