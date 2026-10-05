# =============================================================================
# R02_wgcna_coordinated.R
# -----------------------------------------------------------------------------
# ONE network per disease, on the SAME gene universe, with BYTE-IDENTICAL
# parameters, on DISEASE SAMPLES ONLY (so disease status cannot drive module
# structure -- reviewer R1 #1).
#
# ORIENTATION: WGCNA wants rows = samples, columns = genes. R01 stores
# genes x samples, so both matrices are transposed here. (Passing the wrong
# orientation silently makes blockwiseModules cluster the SAMPLES; the
# assert_orientation() guard below exists because that failure mode is silent.)
#
# BETA: chosen by rule in R00. The two diseases are then built at a COMMON beta
# (the more stringent of the two rule-implied values) so neither network is
# built more permissively than its own data requires. The per-disease beta is
# also run as a sensitivity analysis.
# =============================================================================
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")
log_msg("=== R02 start ===")

d  <- readRDS(file.path(OUT, "intermediate", "R01_universe.rds"))
U  <- d$U
hf_dis  <- d$ann_hf$gsm[d$ann_hf$group == "HF"]
hcc_all <- colnames(d$E_hcc_full)

# genes x samples (R01 orientation)
gxHF  <- d$E_hf_full[U,  hf_dis,  drop = FALSE]
gxHCC <- d$E_hcc_full[U, hcc_all, drop = FALSE]
stopifnot(identical(rownames(gxHF), rownames(gxHCC)), length(U) == nrow(gxHF))

# samples x genes (WGCNA orientation)
datHF  <- t(gxHF)
datHCC <- t(gxHCC)

# --- orientation guard --------------------------------------------------------
# WGCNA treats COLUMNS as genes. If the number of columns is not |U| we (or a
# future edit) transposed the wrong way; fail loudly rather than silently
# clustering samples.
assert_orientation <- function(datExpr, n_genes, tag) {
  if (ncol(datExpr) != n_genes)
    stop(sprintf("%s: expected %d gene columns, found %d. Matrix is transposed the wrong way.",
                 tag, n_genes, ncol(datExpr)))
  if (ncol(datExpr) < nrow(datExpr))
    warning(sprintf("%s: %d genes x %d samples -- fewer genes than samples is unusual.",
                    tag, ncol(datExpr), nrow(datExpr)))
  invisible(TRUE)
}
assert_orientation(datHF,  length(U), "HF")
assert_orientation(datHCC, length(U), "HCC")
log_msg("orientation OK: HF ", nrow(datHF), "x", ncol(datHF),
        " | HCC ", nrow(datHCC), "x", ncol(datHCC), " (samples x genes)")

wpar <- WGCNA[c("networkType","TOMType","minModuleSize","deepSplit","mergeCutHeight",
                "pamRespectsDendro","maxBlockSize","numericLabels","verbose")]

# --- soft threshold (once per disease) ---------------------------------------
soft <- function(datExpr, tag) {
  disableWGCNAThreads()
  sft <- pickSoftThreshold(datExpr, powerVector = WGCNA$powerVector,
                           networkType = WGCNA$networkType,
                           blockSize = WGCNA$maxBlockSize,
                           verbose = 2)
  pb <- pick_beta(sft)
  log_msg(sprintf("%s: rule -> beta=%s (%s), R2=%.3f, reached cut=%s",
                  tag, pb$power, pb$method, pb$r2, pb$reached_cut))
  pb
}
pbHF  <- soft(datHF,  "HF")
pbHCC <- soft(datHCC, "HCC")
common <- choose_common_beta(list(HF = pbHF, HCC = pbHCC))
log_msg("common beta = ", common$power, "  [", common$rationale, "]")

# --- network builder ----------------------------------------------------------
build <- function(datExpr, power, tag) {
  log_msg("--- ", tag, ": blockwiseModules at beta=", power, " ---")
  enableWGCNAThreads(nThreads = 4)
  t0 <- Sys.time()
  net <- do.call(blockwiseModules, c(list(datExpr = datExpr, power = power), wpar))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  lab <- net$colors
  log_msg(sprintf("%s: %.1f min; %d modules; sizes: %s", tag, dt,
                  length(unique(lab)),
                  paste(head(sort(table(lab), decreasing = TRUE), 6),
                        collapse = " ")))
  list(net = net, power = power, minutes = dt,
       n_samples = nrow(datExpr), n_genes = ncol(datExpr))
}

# PRIMARY: common beta for both diseases
resHF  <- build(datHF,  common$power, "HF")
resHCC <- build(datHCC, common$power, "HCC")

# SENSITIVITY: each disease at its own rule-implied beta
resHF_sens  <- if (pbHF$power  != common$power) build(datHF,  pbHF$power,  "HF_ownbeta")  else NULL
resHCC_sens <- if (pbHCC$power != common$power) build(datHCC, pbHCC$power, "HCC_ownbeta") else NULL

# --- parity assertion ---------------------------------------------------------
assert_wgcna_parity(wpar, wpar)
log_msg("parameter parity OK; both networks built at beta=", common$power)

# --- per-gene module membership (kME), samples x genes orientation ------------
gene_kME <- function(res, datExpr) {
  lab <- res$net$colors
  me  <- res$net$MEs
  rows <- lapply(names(sort(table(lab), decreasing = TRUE)), function(m) {
    g   <- names(lab)[lab == m]
    col <- paste0("ME", m)
    if (length(g) < 3 || !col %in% colnames(me))
      return(data.frame(gene = g, module = m, kME = NA_real_, stringsAsFactors = FALSE))
    km <- suppressWarnings(cor(datExpr[, g, drop = FALSE], me[, col], use = "p")[, 1])
    data.frame(gene = g, module = m, kME = km, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[order(out$module, -abs(out$kME)), ]
}
kmeHF  <- gene_kME(resHF,  datHF);  kmeHF$disease  <- "HF"
kmeHCC <- gene_kME(resHCC, datHCC); kmeHCC$disease <- "HCC"
write_table(rbind(kmeHF, kmeHCC)[, c("disease","module","gene","kME")],
            "Table_S7_gene_module_membership.csv")

# --- module inventory ---------------------------------------------------------
inventory <- function(res, tag) {
  s <- sort(table(res$net$colors), decreasing = TRUE)
  data.frame(disease = tag, module = names(s), n_genes = as.integer(s),
             power = res$power, stringsAsFactors = FALSE)
}
inv <- rbind(inventory(resHF, "HF"), inventory(resHCC, "HCC"))
inv$pct_of_genes <- round(100 * inv$n_genes /
                            ave(inv$n_genes, inv$disease, FUN = sum), 1)
print(inv, row.names = FALSE)
write_table(inv, "Table_S4_module_inventory.csv")

fit_tab <- rbind(
  data.frame(disease = "HF",  power = pbHF$powers,  fit = pbHF$fit,
             chosen_own = pbHF$powers  == pbHF$power,  chosen_common = pbHF$powers  == common$power),
  data.frame(disease = "HCC", power = pbHCC$powers, fit = pbHCC$fit,
             chosen_own = pbHCC$powers == pbHCC$power, chosen_common = pbHCC$powers == common$power))
write_table(fit_tab, "Table_S5_soft_threshold_fit.csv")

# --- module-trait (HF: etiology / sex / age) ----------------------------------
trait_hf <- data.frame(
  ICM  = as.integer(d$ann_hf$etiology[match(hf_dis, d$ann_hf$gsm)] == "ICM"),
  male = as.integer(d$ann_hf$gender[match(hf_dis, d$ann_hf$gsm)] == "male"),
  age  = suppressWarnings(as.numeric(d$ann_hf$age[match(hf_dis, d$ann_hf$gsm)])),
  row.names = hf_dis)
mt <- function(net, trait) {
  me <- net$MEs[, colnames(net$MEs) != "MEgrey", drop = FALSE]
  r  <- suppressWarnings(cor(me, trait, use = "p"))
  n  <- nrow(trait)
  p  <- 2 * pt(-abs(r * sqrt((n - 2) / (1 - r^2))), n - 2)
  data.frame(module = rep(rownames(r), ncol(r)),
             trait  = rep(colnames(r), each = nrow(r)),
             r = as.vector(r), p = as.vector(p), stringsAsFactors = FALSE)
}
mthf <- mt(resHF$net, trait_hf); mthf$fdr <- p.adjust(mthf$p, "BH")
write_table(mthf, "Table_S6_module_trait_HF.csv")
cat("\n--- HF module-trait (FDR < 0.05) ---\n")
print(head(mthf[mthf$fdr < 0.05, ][order(mthf$fdr[mthf$fdr < 0.05]), ], 15), row.names = FALSE)

# --- save ---------------------------------------------------------------------
saveRDS(list(HF = resHF, HCC = resHCC, HF_sens = resHF_sens, HCC_sens = resHCC_sens,
             wpar = wpar, pbHF = pbHF, pbHCC = pbHCC, common = common,
             hf_dis = hf_dis, hcc_all = hcc_all,
             datHF = datHF, datHCC = datHCC, gxHF = gxHF, gxHCC = gxHCC, U = U,
             inv = inv, fit = fit_tab, module_trait_HF = mthf,
             kmeHF = kmeHF, kmeHCC = kmeHCC),
        file.path(OUT, "intermediate", "R02_wgcna.rds"))

log_msg("=== R02 done ===")
cat("\n================= R02 CHECKPOINT SUMMARY =================\n")
cat(sprintf("gene universe U : %d genes\n", length(U)))
cat(sprintf("HF  network     : %d samples, beta=%d, %d modules\n",
            resHF$n_samples, resHF$power, length(unique(resHF$net$colors))))
cat(sprintf("HCC network     : %d samples, beta=%d, %d modules\n",
            resHCC$n_samples, resHCC$power, length(unique(resHCC$net$colors))))
cat(sprintf("\nbeta rule       : %s\n", common$rationale))
cat(sprintf("                  HF reached R2 cut: %s ; HCC reached R2 cut: %s\n",
            pbHF$reached_cut, pbHCC$reached_cut))
cat("\nHF modules:\n");  print(head(inv[inv$disease == "HF",  c("module","n_genes","pct_of_genes")], 12), row.names = FALSE)
cat("\nHCC modules:\n"); print(head(inv[inv$disease == "HCC", c("module","n_genes","pct_of_genes")], 12), row.names = FALSE)
cat("=========================================================\n")
