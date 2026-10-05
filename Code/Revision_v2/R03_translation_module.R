# =============================================================================
# R03_translation_module.R
# -----------------------------------------------------------------------------
# Identify the translation module in each disease QUANTITATIVELY instead of by
# keyword matching on gene names (what manuscript v1 did, and what R3 objected
# to).
#
#   1. ssGSEA score per sample for every canonical translation set
#   2. correlate each set score with each module eigengene
#   3. translation module = argmax over modules of mean |r| across Tier A sets,
#      with a sample-level permutation p-value (max-statistic null, so the
#      p-value is corrected across the whole sets x modules matrix)
#   4. independent secondary check: hypergeometric over-representation of
#      Tier A genes inside each module
#
# ORIENTATION: R02 stores samples x genes (WGCNA); GSVA needs genes x samples,
# so the analysis uses r2$gxHF / r2$gxHCC.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
suppressPackageStartupMessages(library(GSVA))

log_msg("=== R03 start ===")
r2 <- readRDS(file.path(OUT, "intermediate", "R02_wgcna.rds"))

tiers <- list(TierA = SETS_TIER_A,
              TierB = c(SETS_TIER_A, SETS_TIER_B_EXTRA),
              TierC_v1 = SETS_TIER_C_V1)

score_sets <- function(gx, set_list) {
  sets <- resolve_sets(set_list)
  sets <- lapply(sets, function(g) intersect(g, rownames(gx)))
  sets <- sets[vapply(sets, length, integer(1)) >= 5]
  p <- GSVA::ssgseaParam(as.matrix(gx), sets, minSize = 5, maxSize = 500,
                         normalize = TRUE)
  GSVA::gsva(p, verbose = FALSE)      # returns sets x samples
}

analyse <- function(res, gx, tag) {
  me_all <- res$net$MEs
  me_all <- me_all[, colnames(me_all) != "MEgrey", drop = FALSE]
  # module labels are stored as colour names without the "ME" prefix
  colnames(me_all) <- sub("^ME", "", colnames(me_all))
  me <- me_all
  out <- list()
  for (tn in names(tiers)) {
    sc <- score_sets(gx, tiers[[tn]])
    rs <- rownames(sc)
    r  <- suppressWarnings(cor(t(sc), me, use = "p"))     # sets x modules

    # Sample-level permutation null. Two statistics per permutation:
    #  (a) global max |r|          -> per-set p, FWER across the whole matrix
    #  (b) max over modules of the mean |r| across sets
    #                              -> module-level p for the ranking statistic
    set.seed(STATS$SEED)
    nperm <- 1000
    perm <- replicate(nperm, {
      scp <- sc[, sample(ncol(sc)), drop = FALSE]
      rp  <- abs(suppressWarnings(cor(t(scp), me, use = "p")))
      c(max(rp), max(rowMeans(rp)))
    })
    perm_r   <- perm[1, ]
    perm_agg <- perm[2, ]

    agg <- colMeans(abs(r))                     # named by module
    df <- expand.grid(set = rs, module = colnames(r), stringsAsFactors = FALSE)
    df$r         <- as.vector(r)
    df$tier      <- tn
    df$disease   <- tag
    df$p_fwer    <- vapply(abs(df$r),
                           function(x) (1 + sum(perm_r >= x)) / (1 + nperm), numeric(1))
    df$module_agg <- as.numeric(agg[df$module])
    df$p_module   <- vapply(df$module,
                            function(m) (1 + sum(perm_agg >= agg[[m]])) / (1 + nperm),
                            numeric(1))
    out[[tn]] <- df
  }
  do.call(rbind, out)
}

corr_all <- rbind(analyse(r2$HF,  r2$gxHF,  "HF"),
                  analyse(r2$HCC, r2$gxHCC, "HCC"))
write_table(corr_all, "Table_S8_set_module_correlation.csv")

# --- criterion 1: eigengene-track ranking -------------------------------------
# Which module's eigengene best tracks canonical translation-set activity?
tracking <- function(corr, tag, tier = "TierA") {
  d <- corr[corr$disease == tag & corr$tier == tier, ]
  agg <- aggregate(abs(r) ~ module, data = d, FUN = mean)
  names(agg)[2] <- "mean_abs_r"
  p <- unique(d[, c("module", "p_module")])
  agg <- merge(agg, p, by = "module")
  agg <- agg[order(-agg$mean_abs_r), ]
  list(ranking = agg, module = agg$module[1],
       mean_abs_r = agg$mean_abs_r[1], p_module = agg$p_module[1],
       detail = d[d$module == agg$module[1], ][order(-abs(d[d$module == agg$module[1], "r"])), ])
}
trkHF  <- tracking(corr_all, "HF")
trkHCC <- tracking(corr_all, "HCC")

cat("\n=== Criterion 1: eigengene tracking (modules ranked by mean |r| to Tier A sets) ===\n")
cat("\n-- HF --\n");  print(head(trkHF$ranking,  8), row.names = FALSE)
cat("\n-- HCC --\n"); print(head(trkHCC$ranking, 8), row.names = FALSE)
cat(sprintf("\nHF  best-tracking module: %s (mean|r|=%.3f, perm p=%.4f)\n",
            trkHF$module,  trkHF$mean_abs_r,  trkHF$p_module))
cat(sprintf("HCC best-tracking module: %s (mean|r|=%.3f, perm p=%.4f)\n",
            trkHCC$module, trkHCC$mean_abs_r, trkHCC$p_module))

# --- criterion 2: hypergeometric over-representation (module membership) ------
A_genes <- unique(unlist(resolve_sets(SETS_TIER_A)))
hyper <- function(lab, hit, universe) {
  rows <- lapply(names(sort(table(lab), decreasing = TRUE)), function(m) {
    g <- names(lab)[lab == m]
    k <- length(intersect(g, hit)); n <- length(g)
    K <- length(intersect(universe, hit)); N <- length(universe)
    data.frame(module = m, n_module = n, n_overlap = k,
               expected = round(n * K / N, 1),
               fold = round((k / n) / (K / N), 2),
               p = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows); d$fdr <- p.adjust(d$p, "BH"); d[order(d$p), ]
}
hypHF  <- hyper(r2$HF$net$colors,  A_genes, r2$U);  hypHF$disease  <- "HF"
hypHCC <- hyper(r2$HCC$net$colors, A_genes, r2$U); hypHCC$disease <- "HCC"
write_table(rbind(hypHF, hypHCC), "Table_S9_translation_overrepresentation.csv")
cat("\n=== Criterion 2: hypergeometric over-representation of Tier A genes ===\n")
cat("\n-- HF --\n");  print(head(hypHF,  6), row.names = FALSE)
cat("\n-- HCC --\n"); print(head(hypHCC, 6), row.names = FALSE)

# --- designation ---------------------------------------------------------------
# Three candidate rules are computed SIDE BY SIDE so the choice is explicit and
# auditable rather than implicit in the code:
#
#   (a) raw significance  argmin p                  -- size-biased: a module with
#       2,000 genes wins simply by being large (this is how HCC "blue" with
#       fold 1.99 outranked "magenta" with fold 3.99)
#   (b) fold enrichment   argmax fold among FDR<alpha -- density, size-corrected
#   (c) eigengene track   argmax mean|r|              -- activity, not membership
#
# PRIMARY = (b), applied identically to both diseases. Rationale: the question
# "which module is the translation module" is a membership question, and fold
# enrichment is the only one of the three that is not confounded by module size.
# Because (a) is size-biased it would compare a 153-gene module in one disease
# against a 1,948-gene module in the other -- the exact non-equivalence that
# invalidated v1's cross-disease comparison. (c) is reported as a confirmatory
# criterion and disagreement is stated, not hidden.
designate <- function(hyp, trk, tag) {
  sig <- hyp[hyp$fdr < STATS$ALPHA, ]
  if (!nrow(sig))
    return(list(module = NA_character_, note = "no module enriched at FDR<alpha"))
  prim <- sig[which.max(sig$fold), ]          # (b) densest among significant
  by_p <- sig[which.min(sig$p), ]             # (a)
  i_prim <- match(prim$module, trk$ranking$module)
  list(
    module = prim$module, n_module = prim$n_module, n_overlap = prim$n_overlap,
    fold = prim$fold, p = prim$p, fdr = prim$fdr, expected = prim$expected,
    tracking_mean_abs_r = trk$ranking$mean_abs_r[i_prim],
    tracking_rank = i_prim,
    most_significant = by_p$module, most_significant_fold = by_p$fold,
    best_tracking = trk$module,
    agree_with_tracking = identical(prim$module, trk$module))
}
desHF  <- designate(hypHF,  trkHF,  "HF")
desHCC <- designate(hypHCC, trkHCC, "HCC")

cat("\n=== DESIGNATION (rule b: max fold among FDR<alpha, applied to both) ===\n")
for (nm in c("HF", "HCC")) {
  d <- if (nm == "HF") desHF else desHCC
  cat(sprintf("\n%s: translation module = %s\n", nm, d$module))
  cat(sprintf("    size          : %d genes\n", d$n_module))
  cat(sprintf("    Tier A members: %d (expected %.1f) -> fold %.2f, p=%.2e, FDR=%.2e\n",
              d$n_overlap, d$expected, d$fold, d$p, d$fdr))
  cat(sprintf("    tracking      : mean|r|=%.3f, rank %d of %d modules\n",
              d$tracking_mean_abs_r, d$tracking_rank, nrow(trkHF$ranking)))
  cat(sprintf("    rule agreement: by significance=%s (fold %.2f) | by tracking=%s%s\n",
              d$most_significant, d$most_significant_fold, d$best_tracking,
              if (d$agree_with_tracking) " [matches primary]" else " [differs from primary]"))
}

# --- translation module gene lists --------------------------------------------
trans_genes <- list(
  HF  = if (is.na(desHF$module))  character(0) else names(r2$HF$net$colors)[r2$HF$net$colors   == desHF$module],
  HCC = if (is.na(desHCC$module)) character(0) else names(r2$HCC$net$colors)[r2$HCC$net$colors == desHCC$module])
ov <- intersect(trans_genes$HF, trans_genes$HCC)
write_table(data.frame(
  disease = c(rep("HF", length(trans_genes$HF)), rep("HCC", length(trans_genes$HCC))),
  module  = c(rep(desHF$module,  length(trans_genes$HF)),
              rep(desHCC$module, length(trans_genes$HCC))),
  gene    = c(trans_genes$HF, trans_genes$HCC), stringsAsFactors = FALSE),
  "Table_S10_translation_module_genes.csv")

cat(sprintf("\ntranslation modules: HF=%d genes, HCC=%d genes, overlap=%d (Jaccard=%.3f)\n",
            length(trans_genes$HF), length(trans_genes$HCC), length(ov),
            if (length(union(trans_genes$HF, trans_genes$HCC)))
              length(ov) / length(union(trans_genes$HF, trans_genes$HCC)) else NaN))
cat("Tier A genes in HF module: ", sum(trans_genes$HF  %in% A_genes), "\n")
cat("Tier A genes in HCC module:", sum(trans_genes$HCC %in% A_genes), "\n")

saveRDS(list(trkHF = trkHF, trkHCC = trkHCC, desHF = desHF, desHCC = desHCC,
             corr = corr_all, hypHF = hypHF, hypHCC = hypHCC,
             trans_genes = trans_genes, overlap = ov, tierA_genes = A_genes),
        file.path(OUT, "intermediate", "R03_translation_module.rds"))
log_msg("=== R03 done ===")
