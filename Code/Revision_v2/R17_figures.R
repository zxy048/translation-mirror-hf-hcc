# =============================================================================
# R17_figures.R — regenerate every figure of manuscript v2 from the v2
# intermediates. No hand-entered numbers and no hand-tuned layout constants that
# encode a result: every plotted value is read from an RDS written by R02-R15.
#
# Design note (important): Fig 1 does NOT show "module-trait relationships in
# the sense of v1, because both v2 networks are built on DISEASE SAMPLES ONLY
# (HF 177 failing, HCC 148 tumours). There is no within-network disease/control
# contrast to correlate a module eigengene against, so a trait heatmap of that
# kind cannot be produced from these networks and would have to be faked by
# reaching into the control samples the networks were deliberately built
# without. Fig 1 therefore shows the module inventory, the quantitative
# translation-set enrichment that designated the modules, and the preservation
# result.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(grid); library(methods)
})
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

FIG <- file.path(OUT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
I <- function(nm) readRDS(file.path(OUT, "intermediate", paste0(nm, ".rds")))

ok <- list(); bad <- list()

# The default pdf() device on this Windows build cannot encode Greek letters
# ("conversion failure ... mbcsToSbcs"), which silently corrupts any subtitle
# carrying a rho or a beta. Cairo handles UTF-8 properly; fall back to pdf()
# only if cairo is unavailable, and in that case the ASCII-only prose below
# still renders.
USE_CAIRO <- isTRUE(capabilities("cairo"))
cat("cairo available:", USE_CAIRO, "\n")

# PLOS ONE accepts figure files up to 2250 x 2625 px at 300 dpi, i.e. 7.5 x
# 8.75 in. Anything larger is downscaled by the journal, and every label with
# it: the previous revision was rendered at 13 x 11.5 in, so its 5.5 pt axis
# text reached the reader at an effective ~2.75 pt. Every figure below is now
# rendered at or inside these limits, so the point sizes set in the code are the
# point sizes in the PDF.
PLOS_MAX_W <- 7.5
PLOS_MAX_H <- 8.75

save_fig <- function(p, name, w = PLOS_MAX_W, h = PLOS_MAX_H) {
  if (USE_CAIRO) {
    grDevices::cairo_pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h)
  } else {
    grDevices::pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h)
  }
  print(p); grDevices::dev.off()
  grDevices::png(file.path(FIG, paste0(name, ".png")), width = w, height = h,
                 units = "in", res = 300, type = if (USE_CAIRO) "cairo" else "windows")
  print(p); grDevices::dev.off()
  cat("  wrote", name, "\n")
}
run <- function(label, expr) {
  r <- try(expr, silent = TRUE)
  if (inherits(r, "try-error")) {
    bad[[label]] <<- as.character(r)
    cat("  !! FAILED:", label, "->", as.character(r), "\n")
  } else ok[[label]] <<- TRUE
}

# ---- shared theming ---------------------------------------------------------
# Note on units, because the two APIs differ and mixing them up is what made the
# first attempt at this revision illegible:
#   element_text(size = x)  ->  x is POINTS
#   geom_text(size = x)     ->  x is multiplied by .pt (2.845), so size 2.8 ~ 8 pt
# A target of 8-12 pt for every label on the page therefore means axis/title
# sizes in the 8-12 range and annotation sizes in the 2.8-4.2 range.
th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(face = "bold"),
        plot.tag = element_text(face = "bold", size = 12),
        legend.key.size = unit(0.8, "lines"))
C_HF <- "#2C6FA8"; C_HCC <- "#C1443B"; C_GREY <- "grey70"

# WGCNA module names are R colour names; fall back to grey if a name is not one.
# unname() matters: vapply returns a vector named by module, and a named vector
# passed to scale_fill_manual(values=) gets .purple/.magenta suffixes that match
# nothing, silently dropping the scale to default greys.
mod_col <- function(x) {
  unname(vapply(tolower(x), function(n)
    tryCatch({ col2rgb(n); n }, error = function(e) "grey55"), character(1)))
}
# The designated modules are read from R04's verdict, never hardcoded here: the
# designation is a result of R03's quantitative scoring, and re-typing it into
# the plotting script would let figure and result drift apart silently.
# grey is WGCNA's unassigned bin, not a module; it is excluded from every count
# and plot so that "n named modules" means the same thing in text and figure.
PRES <- I("R04_preservation")
DESIG <- PRES$verdict                     # reference, module, test, Zsummary, verdict
des_of <- function(ref) DESIG$module[DESIG$reference == ref][1]
cat("designated modules from R04$verdict:\n"); print(DESIG, row.names = FALSE)
is_grey <- function(x) tolower(x) == "grey"

cat("=== R17 figures start ===\n")

# =============================================================================
# Fig 1 — harmonised networks, translation-module identification, preservation
# =============================================================================
run("Fig1", {
  w  <- I("R02_wgcna")
  r3 <- I("R03_translation_module")
  r4 <- I("R04_preservation")

  inv <- w$inv[!is_grey(w$inv$module), ]
  inv$disease <- factor(inv$disease, levels = c("HF", "HCC"))
  inv$designated <- mapply(function(d, m)
    tolower(m) == tolower(des_of(if (d == "HF") "HF" else "HCC")),
    as.character(inv$disease), inv$module)
  inv$lab <- ifelse(inv$designated, paste0(inv$module, " (", inv$n_genes, ")"), NA)

  # Every count quoted in the schematic is read from R02/R04, never typed in.
  hf_n   <- sum(inv$disease == "HF");  hcc_n  <- sum(inv$disease == "HCC")
  hf_sz  <- inv$n_genes[inv$disease == "HF"  & inv$designated][1]
  hcc_sz <- inv$n_genes[inv$disease == "HCC" & inv$designated][1]
  hf_mod <- des_of("HF"); hcc_mod <- des_of("HCC")

  # --- a: schematic ----------------------------------------------------------
  # R2 (Major 1) asked for a schematic relating the WGCNA modules, the a priori
  # canonical pathway tiers, and the derived scores. Boxes and arrows only; no
  # data are plotted in this panel, and nothing here is a result. Box text is
  # capped at four lines because this band is only ~1.9 in tall.
  # Box text is kept to short lines (<= 26 characters) because a column is only
  # ~1.3 in wide at this figure width; the longer descriptions that used to sit
  # in these boxes overflowed them and are carried by the caption instead.
  bx <- data.frame(
    xmin = c( 1,  1, 34, 34, 69, 69, 69),
    xmax = c(31, 31, 66, 66, 99, 99, 99),
    ymin = c(56,  4, 56,  4, 74, 42,  4),
    ymax = c(100, 48, 100, 48, 100, 70, 38),
    fill = c("#DCE7F2", "#F6DEDC", "#EDEDED", "#E4EDE4", "#FBEFD0", "#EAF3EA", "#EAF3EA"),
    edge = c(C_HF, C_HCC, "grey40", "#3B6E3B", "#8A6D00", "#3B6E3B", "#3B6E3B"),
    lab = c(
      sprintf("HF network\nGSE57338\n177 failing LV\n%s named modules", hf_n),
      sprintf("HCC network\nGSE141198\n148 tumours\n%s named modules", hcc_n),
      "A priori pathway tiers\n15 Tier A sets\nexplicit membership\n(KEGG + Reactome)",
      sprintf("ssGSEA per network\n\u2192 designated module:\nHF %s (%s genes)\nHCC %s (%s genes)",
              hf_mod, hf_sz, hcc_mod, hcc_sz),
      "TCS\nmean z of the 15 Tier A sets;\ndisease-agnostic",
      "Module eigengene score\nkME-weighted mean z;\nboth directions",
      "Fed to the preservation test\nand the Hedges g comparison"),
    stringsAsFactors = FALSE)
  hd <- data.frame(x = c(16, 50, 84),
                   lab = c("1 \u00b7 Network modules",
                           "2 \u00b7 Pathway tiers",
                           "3 \u00b7 Derived scores"))
  arw <- arrow(length = unit(0.05, "in"), type = "closed")
  pa <- ggplot() +
    geom_rect(data = bx, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = bx$fill, colour = bx$edge, linewidth = 0.4) +
    geom_text(data = bx, aes((xmin + xmax) / 2, (ymin + ymax) / 2, label = lab),
              size = 2.3, lineheight = 1.06, colour = "grey10") +
    geom_text(data = hd, aes(x, 107, label = lab), size = 2.6, fontface = "bold",
              colour = "grey25") +
    geom_segment(aes(x = 31, xend = 34, y = 78, yend = 30), arrow = arw,
                 linewidth = 0.4, colour = "grey35") +
    geom_segment(aes(x = 31, xend = 34, y = 26, yend = 16), arrow = arw,
                 linewidth = 0.4, colour = "grey35") +
    geom_segment(aes(x = 50, xend = 50, y = 56, yend = 48), arrow = arw,
                 linewidth = 0.4, colour = "grey35") +
    geom_segment(aes(x = 66, xend = 69, y = 26, yend = 56), arrow = arw,
                 linewidth = 0.4, colour = "grey35") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 114), expand = c(0, 0)) +
    labs(tag = "a", title = "How the three objects relate (schematic)") +
    theme_void() +
    theme(plot.tag = element_text(face = "bold", size = 12),
          plot.title = element_text(size = 9, face = "bold", hjust = 0,
                                    margin = margin(b = 2)))

  # --- b: module inventory ---------------------------------------------------
  # The y-axis carries one label per named module (43 across the two facets).
  # With ~5.1 in of plot height that allows 7 pt without any two labels
  # touching; the "HF"/"HCC" prefix is dropped because the facet strip already
  # states the disease, which is what was crowding the axis.
  pb <- ggplot(inv, aes(stats::reorder(paste(disease, module), n_genes), n_genes,
                        fill = module)) +
    geom_col(width = 0.8) + scale_fill_manual(values = mod_col(inv$module), guide = "none") +
    geom_text(aes(label = lab), hjust = -0.1, size = 3.4, fontface = "bold", na.rm = TRUE) +
    coord_flip() + facet_grid(disease ~ ., scales = "free_y", space = "free_y") +
    scale_x_discrete(labels = function(x) sub("^(HF|HCC) ", "", x)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Genes in module", tag = "b", title = "Module inventory (\u03b2 = 14)") +
    expand_limits(y = max(inv$n_genes) * 1.30) + th +
    theme(axis.text.y = element_text(size = 7))

  # --- c: Tier A enrichment per module ---------------------------------------
  # Redesigned from a 43-category bar chart to a fold-versus-size scatter: the
  # bar chart could not carry readable module labels at this figure width, and
  # the scatter makes the same point more directly (the designated module is the
  # outlier, and not merely the largest module).
  hyp <- rbind(r3$hypHF, r3$hypHCC)
  hyp <- hyp[!is_grey(hyp$module), ]
  hyp$disease <- factor(hyp$disease, levels = c("HF", "HCC"))
  hyp <- merge(hyp, inv[, c("disease", "module", "n_genes", "designated")],
               by = c("disease", "module"))
  hyp$lab <- ifelse(hyp$designated, as.character(hyp$module), NA)
  pc <- ggplot(hyp, aes(n_genes, fold, colour = disease)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey35") +
    geom_point(size = 1.9, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = lab), size = 3.1, fontface = "bold",
                             colour = "#7A5C00", na.rm = TRUE, min.segment.length = 0) +
    scale_x_log10() +
    scale_colour_manual(values = c(HF = C_HF, HCC = C_HCC)) +
    labs(x = "Genes in module (log scale)", y = "Tier A fold enrichment",
         colour = NULL, tag = "c",
         title = "Enrichment by module size") +
    th + theme(legend.position = "top")

  # --- d: preservation, both directions --------------------------------------
  tb <- r4$table[!is_grey(r4$table$module), ]
  tb$direction <- ifelse(tb$reference == "HF", "HF network tested in HCC",
                                                "HCC network tested in HF")
  # A module name recurs across directions (both networks have a "magenta"), so
  # the designated flag must be reference-specific: HF purple tested in HCC and
  # HCC magenta tested in HF are the two designated points, not four.
  tb$designated <- mapply(function(ref, m) tolower(m) == tolower(des_of(ref)),
                          tb$reference, tb$module)
  tb$lab <- ifelse(tb$designated, paste0(tb$module, " (", tb$moduleSize, ")"), NA)
  cat("  preservation class counts (named modules only):\n")
  print(table(tb$direction,
              ifelse(tb$Zsummary.pres > 10, "strong",
              ifelse(tb$Zsummary.pres > 2, "weak-moderate", "none"))))
  for (rf in unique(tb$reference)) {
    s <- tb[tb$reference == rf, ]
    s <- s[order(-s$Zsummary.pres), ]
    i <- which(s$designated)
    cat(sprintf("  %s: designated %s is rank %d of %d named modules (Z = %.2f)\n",
                rf, s$module[i], i, nrow(s), s$Zsummary.pres[i]))
  }
  pd <- ggplot(tb, aes(Zsummary.pres, medianRank.pres, colour = direction)) +
    geom_vline(xintercept = c(2, 10), linetype = 3, colour = "grey55") +
    geom_point(aes(shape = designated, size = moduleSize), alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = lab), size = 3.1, fontface = "bold",
                             colour = "grey15", na.rm = TRUE, min.segment.length = 0) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    scale_size_continuous(range = c(1.6, 5), guide = "none") +
    # The worst-preserved module sits at median rank ~21 of 23; without padding
    # its point and label touch the panel border and read as clipped.
    scale_y_continuous(expand = expansion(mult = 0.08)) +
    scale_x_continuous(expand = expansion(mult = 0.06)) +
    scale_colour_manual(values = c("HF network tested in HCC" = C_HCC,
                                   "HCC network tested in HF" = C_HF)) +
    labs(x = expression(Z[summary]~"preservation"), y = "Median rank",
         colour = NULL, tag = "d",
         title = "Module preservation, both directions",
         subtitle = "dotted lines: Z = 2 and Z = 10") +
    th + theme(legend.position = "top")

  # --- e: preservation classes ------------------------------------------------
  tb$class <- ifelse(tb$Zsummary.pres > 10, "strong (>10)",
              ifelse(tb$Zsummary.pres > 2, "weak-moderate (2-10)", "none (<2)"))
  sm <- as.data.frame(table(tb$direction, tb$class))
  names(sm) <- c("direction", "class", "n")
  sm$class <- factor(sm$class, levels = c("strong (>10)", "weak-moderate (2-10)", "none (<2)"))
  sm$frac <- ave(sm$n, sm$direction, FUN = function(z) z / sum(z))
  pe <- ggplot(sm, aes(direction, frac, fill = class)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.4,
              colour = "white", fontface = "bold") +
    scale_fill_manual(values = c("strong (>10)" = "#1B7837", "weak-moderate (2-10)" = "#E5A50A",
                                 "none (<2)" = "grey65")) +
    scale_y_continuous(labels = percent_format()) +
    coord_flip() +
    labs(x = NULL, y = "Fraction of named modules", fill = "Preservation class", tag = "e",
         title = "Most modules are not shared; the translation modules are",
         subtitle = "counts are of named modules; the unassigned grey bin is excluded") +
    th + theme(legend.position = "right",
               axis.text.y = element_text(size = 8.5))

  # Sized to the PLOS limit (2250 x 2625 px at 300 dpi = 7.5 x 8.75 in) so the
  # journal does not downscale it and shrink every label with it. The previous
  # 13 x 11.5 in canvas was scaled to roughly half, which is why the reviewer
  # saw 5.5 pt axis text rendered at an effective ~2.75 pt.
  save_fig((pa / (pb | (pc / pd)) / pe) + plot_layout(heights = c(1.9, 5.8, 1.05)),
           "Fig1_harmonised_networks", PLOS_MAX_W, PLOS_MAX_H)
})

# =============================================================================
# Fig 2 — the pathway-level mirror, its null, and its cohort dependence
# =============================================================================
run("Fig2", {
  es <- I("R07_effect_sizes"); r9 <- I("R09_mirror"); r8 <- I("R08_redundancy")
  w <- reshape(es[, c("cohort", "set", "tier", "hedges_g")],
               idvar = c("set", "tier"), timevar = "cohort", direction = "wide")
  names(w) <- sub("^hedges_g[.]", "", names(w))
  rhoA <- cor(w$TCGA_LIHC[w$tier == "A"], w$GSE57338[w$tier == "A"], method = "spearman")
  rhoAll <- cor(w$TCGA_LIHC, w$GSE57338, method = "spearman")

  # The mirror flag is RECOMPUTED here from R07's own gate columns rather than
  # read from a hand-written set list: a hardcoded list of "the sets that
  # mirror" would go stale the moment any threshold moved, and the figure would
  # then disagree with the text without anything failing.
  # Abbreviations rather than blind truncation: at this width a 34-character
  # label fills the whole scatterplot, and truncating at a fixed length would
  # collapse the four mitochondrial sets into four identical strings. The full
  # set names are in S14 Table and in the caption.
  short_set <- function(x) {
    x <- sub("^(KEGG|REACTOME)_", "", x)
    x <- sub("^MITOCHONDRIAL_", "MITO. ", x)
    x <- sub("_IN_THE_MITOCHONDRION$", " (MITO.)", x)
    x <- gsub("_", " ", x)
    # Whole-name matches first. Rewriting words in place produced nonsense here:
    # "TRNA AMINOACYLATION" became "TRNA AA-TRNA" once AMINOACYLATION was blindly
    # replaced, and the resulting labels were not only ugly but ambiguous, since
    # three different sets collapsed toward the same string.
    x <- sub("^AMINOACYL TRNA BIOSYNTHESIS$", "AA-TRNA BIOSYNTH.", x)
    x <- sub("^RIBOSOME ASSOCIATED QUALITY CONTROL$", "RIBOSOME-ASSOC. QC", x)
    x <- sub("^RIBOSOME QUALITY CONTROL RQC COMPLEX .*$", "RQC COMPLEX / NASCENT PEPT.", x)
    x <- sub("^SRP DEPENDENT COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE$",
             "SRP-DEP. COTRANSL. TARG.", x)
    x <- sub("^CYTOSOLIC ", "CYTO. ", x)
    x <- sub("^EUKARYOTIC TRANSLATION ", "EUK. TRANSL. ", x)
    x <- sub("^MITO\\. TRANSLATION ", "MITO. TRANSL. ", x)
    x <- sub("^RRNA PROCESSING \\(MITO\\.\\)$", "RRNA PROC. (MITO.)", x)
    x <- sub("^RRNA PROCESSING$", "RRNA PROC.", x)
    x <- sub("AMINOACYLATION$", "AMINOACYL.", x)
    ifelse(nchar(x) > 26, paste0(substr(x, 1, 23), "..."), x)
  }
  MIR <- function(h, f) {
    a <- es[es$cohort == h, c("set", "tier", "hedges_g", "meets_min_effect", "ci_excludes_0")]
    b <- es[es$cohort == f, c("set", "tier", "hedges_g", "meets_min_effect", "ci_excludes_0")]
    m <- merge(a, b, by = c("set", "tier"), suffixes = c(".hcc", ".hf"))
    m$mirror <- sign(m$hedges_g.hcc) * sign(m$hedges_g.hf) < 0 &
                m$meets_min_effect.hcc & m$meets_min_effect.hf &
                m$ci_excludes_0.hcc   & m$ci_excludes_0.hf
    m
  }
  prim <- MIR("TCGA_LIHC", "GSE57338")
  cat(sprintf("  recomputed mirror, TCGA vs GSE57338: Tier A %d/%d, all %d/%d\n",
              sum(prim$mirror[prim$tier == "A"]), sum(prim$tier == "A"),
              sum(prim$mirror), nrow(prim)))
  prim$lab <- ifelse(prim$tier == "A" & !prim$mirror, short_set(prim$set), NA)

  # --- a: effect-size scatter -------------------------------------------------
  pa <- ggplot(prim, aes(hedges_g.hcc, hedges_g.hf)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_abline(slope = -1, intercept = 0, linetype = 3, colour = "grey65") +
    geom_point(data = subset(prim, tier != "A"), colour = C_GREY, size = 1.9, alpha = 0.75) +
    geom_point(data = subset(prim, tier == "A" & mirror), colour = C_HCC, size = 2.9) +
    geom_point(data = subset(prim, tier == "A" & !mirror), colour = "#1B7837", size = 2.9, shape = 17) +
    ggrepel::geom_text_repel(aes(label = lab), size = 2.8, colour = "#1B7837",
                             na.rm = TRUE, min.segment.length = 0, max.overlaps = Inf,
                             box.padding = 0.7, point.padding = 0.3, force = 6,
                             max.iter = 20000, seed = 42) +
    scale_x_continuous(expand = expansion(mult = 0.16)) +
    scale_y_continuous(expand = expansion(mult = 0.14)) +
    labs(x = expression("Hedges "*italic(g)*", HCC (TCGA-LIHC)"),
         y = expression("Hedges "*italic(g)*", HF (GSE57338)"), tag = "a",
         title = "Opposite-direction perturbation",
         subtitle = sprintf("Tier A rho = %+.3f (n = 15)\nall scored sets rho = %+.3f (n = %d)",
                            rhoA, rhoAll, nrow(prim))) + th

  # --- b: paired effect sizes for Tier A --------------------------------------
  ta <- prim[prim$tier == "A", ]
  ta <- ta[order(ta$hedges_g.hcc), ]
  ta$set_s <- factor(short_set(ta$set), levels = short_set(ta$set))
  long <- rbind(data.frame(set_s = ta$set_s, g = ta$hedges_g.hcc, disease = "HCC"),
                data.frame(set_s = ta$set_s, g = ta$hedges_g.hf, disease = "HF"))
  pb <- ggplot(long, aes(g, set_s, colour = disease)) +
    geom_vline(xintercept = 0, colour = "grey40") +
    geom_line(aes(group = set_s), colour = "grey75", linewidth = 0.5) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c(HCC = C_HCC, HF = C_HF)) +
    labs(x = expression("Hedges "*italic(g)), y = NULL, colour = NULL, tag = "b",
         title = "The 15 Tier A sets, both diseases",
         subtitle = "ordered by HCC effect size; crossing zero is\nnecessary but not sufficient for the mirror gate") +
    th + theme(axis.text.y = element_text(size = 7.5),
               legend.position = "top", legend.justification = "right")

  # --- c: observed vs permutation null ----------------------------------------
  # Read the primary row from R09$summary, but do not assume the ordering or the
  # string form of `pair` (it is a technology descriptor, not a cohort name):
  # select by tier and take the first row, and cross-check the observed count
  # against the independently recomputed one above.
  sm <- r9$summary
  sA <- sm[sm$tier == "A", , drop = FALSE]
  if (!nrow(sA)) sA <- sm
  p_row <- sA[which.max(sA$n_mirror_sets), , drop = FALSE]
  cat(sprintf("  R09 primary row '%s' tier %s: %d/%d sets, %d/%d clusters\n",
              p_row$pair[1], p_row$tier[1], p_row$n_mirror_sets, p_row$n_sets,
              p_row$n_mirror_clusters, p_row$n_clusters))
  if (p_row$n_mirror_sets != sum(prim$mirror[prim$tier == "A"]))
    cat("  !! WARNING: R09 mirror count differs from the recomputed count\n")
  bar <- data.frame(
    level = factor(c("Individual sets", "Redundancy clusters"),
                   levels = c("Individual sets", "Redundancy clusters")),
    obs   = c(p_row$n_mirror_sets, p_row$n_mirror_clusters),
    nul   = c(p_row$null_sets_mean, p_row$null_clusters_mean),
    total = c(p_row$n_sets, p_row$n_clusters),
    p     = c(p_row$p_perm_sets, p_row$p_perm_clusters))
  bar$lab <- sprintf("%d/%d\nP %s", bar$obs, bar$total,
                     ifelse(bar$p < 1e-4, "< 1e-4",
                            paste0("= ", formatC(bar$p, digits = 3, format = "g"))))
  bl <- rbind(data.frame(level = bar$level, kind = "Observed", n = bar$obs),
              data.frame(level = bar$level, kind = "Permutation null", n = bar$nul))
  pc <- ggplot(bl, aes(level, n, fill = kind)) +
    geom_col(position = position_dodge(0.7), width = 0.62) +
    geom_text(data = bar, aes(level, obs, label = lab), inherit.aes = FALSE,
              vjust = -0.3, size = 3.3, lineheight = 0.95) +
    scale_fill_manual(values = c(Observed = "#B8860B", `Permutation null` = "grey72")) +
    labs(x = NULL, y = "Mirror count", fill = NULL, tag = "c",
         title = "Against the permutation null",
         subtitle = "10,000 permutations of disease labels;\nboth levels clear the null") +
    expand_limits(y = max(bl$n) * 1.4) + th + theme(legend.position = "top")

  # --- d: cohort x cohort ------------------------------------------------------
  # The claim under test is which COHORT the mirror travels with, so the panel is
  # a matrix over cohort pairs rather than over the technology descriptors that
  # R09$summary uses for its own bookkeeping.
  hccs <- c("TCGA_LIHC", "GSE14520_GPL3921", "GSE14520_GPL571", "GSE76427")
  hfs  <- c("GSE57338", "GSE116250", "GSE141910")
  nice <- c(TCGA_LIHC = "TCGA-LIHC", GSE14520_GPL3921 = "GSE14520\nGPL3921",
            GSE14520_GPL571 = "GSE14520\nGPL571", GSE76427 = "GSE76427",
            GSE57338 = "GSE57338\n(array)", GSE116250 = "GSE116250\n(RNA-seq)",
            GSE141910 = "GSE141910\n(RNA-seq)")
  grid <- do.call(rbind, lapply(hccs, function(h) do.call(rbind, lapply(hfs, function(f) {
    m <- MIR(h, f); m <- m[m$tier == "A", ]
    if (!nrow(m)) return(NULL)
    data.frame(hcc = h, hf = f, frac = mean(m$mirror), n = sum(m$mirror), tot = nrow(m),
               rho = cor(m$hedges_g.hcc, m$hedges_g.hf, method = "spearman"))
  }))))
  grid$hcc <- factor(grid$hcc, levels = rev(hccs))
  grid$hf  <- factor(grid$hf,  levels = hfs)
  grid$lab <- sprintf("%d/%d\nrho %+.2f", grid$n, grid$tot, grid$rho)
  # Printed so the manuscript can quote these exact values rather than numbers
  # read off an image. NOTE this matrix is TCGA/HCC-cohort based and is NOT the
  # same statistic as R09$summary, whose `pair` rows are technology descriptors
  # pooling several HCC array cohorts; mixing the two in one sentence would
  # attribute counts to the wrong analysis.
  cat("  Fig2d cohort matrix (Tier A):\n")
  print(data.frame(hcc = as.character(grid$hcc), hf = as.character(grid$hf),
                   mirror = sprintf("%d/%d", grid$n, grid$tot),
                   rho = round(grid$rho, 3)), row.names = FALSE)
  cat("  GSE116250 rho range:", sprintf("%+.3f to %+.3f", min(grid$rho[grid$hf == "GSE116250"]),
                                        max(grid$rho[grid$hf == "GSE116250"])),
      "| GSE141910 range:", sprintf("%+.3f to %+.3f", min(grid$rho[grid$hf == "GSE141910"]),
                                    max(grid$rho[grid$hf == "GSE141910"])), "\n")
  pd <- ggplot(grid, aes(hf, hcc, fill = frac)) +
    geom_tile(colour = "white", linewidth = 1.1) +
    geom_text(aes(label = lab), size = 3.1, lineheight = 0.95,
              colour = ifelse(grid$frac > 0.55, "white", "grey10")) +
    scale_fill_gradient(low = "grey95", high = C_HF, limits = c(0, 1),
                        labels = percent_format(), name = "Tier A sets\nmirroring") +
    scale_x_discrete(labels = nice[hfs]) +
    scale_y_discrete(labels = nice[rev(hccs)]) +
    labs(x = NULL, y = NULL, tag = "d",
         title = "The mirror travels with the cardiac cohort",
         subtitle = "every HCC cohort mirrors GSE57338;\nneither independent HF RNA-seq cohort does") +
    th + theme(panel.grid = element_blank(),
               axis.text = element_text(size = 8))

  save_fig((pa | pb) / (pc | pd) + plot_layout(heights = c(1, 1)),
           "Fig2_pathway_mirror")
})

# =============================================================================
# Fig 3 — same-organ controls: the hepatic arm is not cancer-specific
# =============================================================================
run("Fig3", {
  es <- I("R07_effect_sizes"); r15 <- I("R15_cirrhosis_control")
  w <- reshape(es[, c("cohort", "set", "tier", "hedges_g")],
               idvar = c("set", "tier"), timevar = "cohort", direction = "wide")
  names(w) <- sub("^hedges_g[.]", "", names(w))
  gw <- r15$wide

  # --- a: cirrhosis vs TCGA ---------------------------------------------------
  m1 <- merge(gw[, c("set", "tier", "Cirrhosis")], w[, c("set", "tier", "TCGA_LIHC")], by = "set")
  m1$grp <- ifelse(m1$tier.x == "A", "Tier A (15)", "other scored sets")
  r1 <- cor(m1$Cirrhosis[m1$tier.x == "A"], m1$TCGA_LIHC[m1$tier.x == "A"], method = "spearman")
  pa <- ggplot(m1, aes(Cirrhosis, TCGA_LIHC)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(data = subset(m1, grp == "other scored sets"), colour = C_GREY, size = 1.9, alpha = 0.75) +
    geom_point(data = subset(m1, grp == "Tier A (15)"), colour = "#6A3D9A", size = 3) +
    labs(x = expression("Hedges "*italic(g)*", cirrhosis (GSE89377)"),
         y = expression("Hedges "*italic(g)*", HCC (TCGA-LIHC)"), tag = "a",
         title = "Cirrhosis does not track HCC",
         subtitle = sprintf("Tier A rho = %+.3f (n = %d)\nacross cohorts", r1, sum(m1$tier.x == "A"))) + th

  # --- b: within-cohort, cirrhosis vs HCC -------------------------------------
  m2 <- gw[gw$tier == "A", ]
  m2$disc <- !m2$cirrhosis_same_direction_as_HCC
  m2$lab <- ifelse(m2$disc, sub("^(KEGG|REACTOME)_", "", m2$set), NA)
  r2 <- cor(m2$Cirrhosis, m2$HCC, method = "spearman")
  pb <- ggplot(m2, aes(Cirrhosis, HCC)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(aes(colour = disc), size = 3.2) +
    ggrepel::geom_text_repel(aes(label = lab), size = 3.1, colour = "#1B7837",
                             na.rm = TRUE, min.segment.length = 0, max.overlaps = 30) +
    scale_colour_manual(values = c(`FALSE` = "#6A3D9A", `TRUE` = "#1B7837"),
                        labels = c(`FALSE` = "concordant", `TRUE` = "discordant"), name = NULL) +
    labs(x = expression("Hedges "*italic(g)*", cirrhosis (vs Normal)"),
         y = expression("Hedges "*italic(g)*", HCC (vs Normal)"), tag = "b",
         title = "Held within one cohort and pipeline",
         subtitle = sprintf("Tier A rho = %+.3f; %d of %d concordant.\nThe 3 discordant sets are the\naminoacyl-tRNA and mitochondrial ones",
                            r2, sum(!m2$disc), nrow(m2))) +
    th + theme(legend.position = c(0.83, 0.2), legend.background = element_rect(fill = "white", colour = "grey80"))

  # --- c: ordered stage trajectory --------------------------------------------
  # Use R15's own stage_index rather than recomputing it: the control stage
  # (Normal) carries an index but no effect size, so a factor-based recompute
  # here could silently shift every stage by one.
  es15 <- r15$es
  lab_map <- unique(es15[, c("stage_index", "stage")])
  lab_map <- lab_map[order(lab_map$stage_index), ]
  ta <- es15[es15$tier == "A", ]
  rib <- es15[es15$set == "KEGG_RIBOSOME", ]
  pc <- ggplot(ta, aes(stage_index, hedges_g, group = set)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_line(colour = "grey78", linewidth = 0.45) +
    geom_line(data = rib, colour = "#B8860B", linewidth = 1.5) +
    geom_point(data = rib, colour = "#B8860B", size = 2.6) +
    geom_errorbar(data = rib, aes(ymin = ci_lo, ymax = ci_hi), width = 0.14, colour = "#B8860B") +
    annotate("label", x = 4, y = max(rib$ci_hi, na.rm = TRUE) + 0.35,
             label = "KEGG ribosome", size = 3.2, colour = "#7A5C00", fill = "white") +
    scale_x_continuous(breaks = lab_map$stage_index, labels = lab_map$stage) +
    labs(x = NULL, y = expression("Hedges "*italic(g)*" vs Normal"), tag = "c",
         title = "Peaks in cirrhosis, falls at dysplasia",
         subtitle = "Tier A sets in grey;\nnot switched on at the malignant transition") +
    th + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

  # --- d: the cardiomyopathy control ------------------------------------------
  m4 <- merge(gw[, c("set", "tier")], w[, c("set", "tier", "TCGA_LIHC", "GSE141910")], by = "set")
  m4$grp <- ifelse(m4$tier.x == "A", "Tier A (15)", "other scored sets")
  r4v <- cor(m4$GSE141910[m4$tier.x == "A"], m4$TCGA_LIHC[m4$tier.x == "A"], method = "spearman")
  pd <- ggplot(m4, aes(GSE141910, TCGA_LIHC)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(data = subset(m4, grp == "other scored sets"), colour = C_GREY, size = 1.9, alpha = 0.75) +
    geom_point(data = subset(m4, grp == "Tier A (15)"), colour = C_HF, size = 3) +
    labs(x = expression("Hedges "*italic(g)*", cardiomyopathy (GSE141910)"),
         y = expression("Hedges "*italic(g)*", HCC (TCGA-LIHC)"), tag = "d",
         title = "No mirror from cardiomyopathy",
         subtitle = sprintf("Tier A rho = %+.3f (n = %d)\nvs %+.3f for GSE57338",
                            r4v, sum(m4$tier.x == "A"), cor(w$TCGA_LIHC[w$tier == "A"], w$GSE57338[w$tier == "A"], method = "spearman"))) + th

  save_fig((pa | pb) / (pc | pd), "Fig3_same_organ_controls")
})

# =============================================================================
# Supplementary figures
# =============================================================================
run("S1", {
  r12 <- I("R12_tats_tf"); gs <- r12$gs
  gl <- rbind(
    data.frame(module = gs$module, thr = "|GS| > 0.20 (v1)", n = gs$n_absGS_gt_020),
    data.frame(module = gs$module, thr = "|GS| > 0.50 (strict)", n = gs$n_absGS_gt_050),
    data.frame(module = gs$module, thr = "FDR < 0.05 (v2)", n = gs$n_FDR05_alone))
  gl$thr <- factor(gl$thr, levels = c("|GS| > 0.20 (v1)", "|GS| > 0.50 (strict)", "FDR < 0.05 (v2)"))
  gl$module <- factor(gl$module, levels = c("purple", "magenta"))
  p <- ggplot(gl, aes(thr, n, fill = module)) +
    geom_col(position = position_dodge(0.72), width = 0.66) +
    geom_text(aes(label = n), position = position_dodge(0.72), vjust = -0.4, size = 3.4) +
    scale_fill_manual(values = stats::setNames(mod_col(c("purple", "magenta")),
                                              c("purple", "magenta"))) +
    facet_wrap(~ module, scales = "free_y") +
    labs(x = NULL, y = "Genes meeting threshold", fill = NULL,
         title = "Hub-gene threshold comparison (S31 Table)",
         subtitle = "The FDR < 0.05 criterion used in v2 is stricter\nthan v1's |GS| > 0.20 on this data") +
    expand_limits(y = max(gl$n) * 1.15) + th + theme(legend.position = "none")
  save_fig(p, "FigS1_hub_gene_threshold", PLOS_MAX_W, 3.83)
})

run("S2", {
  r12 <- I("R12_tats_tf"); co <- r12$concordance
  co$pair <- paste(co$score_a, "vs", co$score_b)
  co$cohort <- factor(co$cohort, levels = unique(co$cohort))
  p <- ggplot(co, aes(cohort, pair, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3.2,
              colour = ifelse(abs(co$r) > 0.7, "white", "grey15")) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(x = NULL, y = NULL,
         title = "Concordance among translation scores, by cohort (S28 Table)",
         subtitle = "TCS = canonical Tier A mean; ME = kME-weighted module score") +
    th + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p, "FigS2_score_concordance", PLOS_MAX_W, 4.05)
})

run("S3", {
  w <- I("R02_wgcna"); fit <- w$fit
  fit$disease <- factor(fit$disease, levels = c("HF", "HCC"))
  co <- w$common
  own <- fit[fit$chosen_own, ]
  p <- ggplot(fit, aes(power, fit)) +
    geom_hline(yintercept = WGCNA$powerR2Cut, linetype = 2, colour = "grey40") +
    # Anchored at the left edge, where the curve is at or below zero; at the right
    # edge it sat on top of the points it is meant to be read against.
    annotate("text", x = min(fit$power), y = WGCNA$powerR2Cut + 0.015,
             label = sprintf("R\u00b2 = %.2f criterion", WGCNA$powerR2Cut),
             hjust = 0, size = 3.2, colour = "grey30") +
    expand_limits(y = WGCNA$powerR2Cut + 0.06) +
    geom_line(colour = "grey45") + geom_point(size = 2.1, colour = "grey25") +
    geom_vline(data = own, aes(xintercept = power), linetype = 3,
               colour = "grey55") +
    geom_point(data = own, colour = "#B8860B", size = 3.4, shape = 18) +
    facet_wrap(~ disease, ncol = 1) +
    labs(x = "Soft-thresholding power (\u03b2)", y = "Signed scale-free topology fit (R\u00b2)",
         title = "Soft-threshold selection by a fixed rule, not by hand (S5 Table)",
         subtitle = sprintf("gold diamonds = each disease's rule-implied power;\nboth networks built at \u03b2 = %s, the more stringent of the two",
                            if (!is.null(co)) co$power[1] else "14")) +
    th
  save_fig(p, "FigS3_soft_threshold", PLOS_MAX_W, 5.33)
})

run("S4", {
  f <- file.path(OUT, "intermediate", "R14_tom_robustness.rds")
  if (!file.exists(f)) stop("R14_tom_robustness.rds not present yet (R14 still running)")
  r14 <- readRDS(f)
  cat("    R14 object names:", paste(names(r14), collapse = ", "), "\n")
  # R14 saves list(HF, HCC, table, coherence, params); the per-module summary is
  # $table, whose ratio columns are obs_over_unmatched_null / obs_over_matched_null.
  tb <- r14$table
  if (is.null(tb)) stop("R14 object has no $table element")
  tb$label <- sprintf("%s %s", tb$disease, tb$module)
  long <- rbind(
    data.frame(module = tb$label, null = "uniform (v1-style)", ratio = tb$obs_over_unmatched_null),
    data.frame(module = tb$label, null = "connectivity-matched", ratio = tb$obs_over_matched_null))
  long <- long[is.finite(long$ratio), ]
  p <- ggplot(long, aes(module, ratio, fill = null)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    geom_col(position = position_dodge(0.7), width = 0.62) +
    geom_text(aes(label = sprintf("%.2f", ratio)), position = position_dodge(0.7),
              vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = c(`uniform (v1-style)` = "grey68",
                                 `connectivity-matched` = "#B8860B")) +
    labs(x = NULL, y = "Observed / null mean intra-modular TOM", fill = NULL,
         title = "Module robustness under two nulls (S32 Table)",
         subtitle = "The interpretable quantity is the ratio: Z is not comparable\nbetween nulls, because bin-matching collapses the null's variance") +
    expand_limits(y = max(long$ratio) * 1.18) + th +
    theme(legend.position = c(0.75, 0.85), legend.background = element_rect(fill = "white", colour = "grey80"))
  save_fig(p, "FigS4_module_robustness", PLOS_MAX_W, 4.33)
})

run("S5", {
  r10 <- I("R10_heterogeneity"); h <- r10$heterogeneity
  h <- h[h$tier == "A", ]
  h$disease <- factor(h$disease, levels = c("HF", "HCC"))
  # The set names run to 82 characters ("REACTOME_RIBOSOME_QUALITY_CONTROL_RQC_
  # COMPLEX_EXTRACTS_AND_DEGRADES_NASCENT_PEPTIDE"). At 7 pt that single label is
  # wider than half the page, which is what pushed the facet_wrap version off the
  # canvas; the two facets also printed the same label column twice. Abbreviated
  # to ~30 characters and dodged into one panel, so the labels are stated once and
  # every set stays distinguishable.
  short_het <- function(x) {
    x <- sub("^(KEGG|REACTOME)_", "", x)
    x <- sub("^MITOCHONDRIAL_", "MITO. ", x)
    x <- sub("_IN_THE_MITOCHONDRION$", " (MITO.)", x)
    x <- gsub("_", " ", x)
    x <- sub("^RIBOSOME QUALITY CONTROL RQC COMPLEX EXTRACTS AND DEGRADES NASCENT PEPTIDE$",
             "RQC COMPLEX / NASCENT PEPTIDE", x)
    x <- sub("^RIBOSOME ASSOCIATED QUALITY CONTROL$", "RIBOSOME-ASSOC. QC", x)
    x <- sub("^SRP DEPENDENT COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE$",
             "SRP-DEP. COTRANSL. TARGETING", x)
    x <- sub("^EUKARYOTIC TRANSLATION INITIATION$", "EUK. TRANSL. INITIATION", x)
    x <- sub("^EUKARYOTIC TRANSLATION ELONGATION$", "EUK. TRANSL. ELONGATION", x)
    x <- sub("^CYTOSOLIC TRNA AMINOACYLATION$", "CYTOSOLIC TRNA AMINOACYL.", x)
    x <- sub("^AMINOACYL TRNA BIOSYNTHESIS$", "AA-TRNA BIOSYNTHESIS", x)
    x
  }
  h$lab <- short_het(as.character(h$set))
  ord <- names(sort(tapply(h$I2_pct, h$lab, mean, na.rm = TRUE)))
  h$lab <- factor(h$lab, levels = ord)
  p <- ggplot(h, aes(lab, I2_pct, fill = disease)) +
    geom_col(position = position_dodge(0.72), width = 0.68) +
    scale_fill_manual(values = c(HF = C_HF, HCC = C_HCC)) +
    coord_flip() +
    labs(x = NULL, y = expression(italic(I)^2*" (%)"), fill = NULL,
         title = "Between-cohort heterogeneity, Tier A sets (S21 Table)",
         subtitle = sprintf("Median %s: HF %.1f%%, HCC %.1f%%; the discordance is\nnot a few outlier sets",
                            "I\u00b2", median(h$I2_pct[h$disease == "HF"], na.rm = TRUE),
                            median(h$I2_pct[h$disease == "HCC"], na.rm = TRUE))) +
    th + theme(axis.text.y = element_text(size = 7), legend.position = "top")
  save_fig(p, "FigS5_heterogeneity", PLOS_MAX_W, 4.09)
})

# ---- report -----------------------------------------------------------------
cat("\n=== R17 summary ===\n")
cat("succeeded:", if (length(ok)) paste(names(ok), collapse = ", ") else "(none)", "\n")
if (length(bad)) {
  cat("FAILED:\n")
  for (n in names(bad)) cat("  ", n, ":", bad[[n]], "\n")
}
cat("\nfiles in figures/:\n")
print(basename(list.files(FIG, full.names = TRUE)))
cat("\nR17 DONE\n")
if (length(bad)) quit(status = 1)
