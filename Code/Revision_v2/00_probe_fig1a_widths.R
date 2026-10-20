# Measure the Fig 1a schematic labels instead of guessing at the font size.
#
# Panel a carries three analysis tracks (module preservation; translation module
# and derived scores; pathway effect sizes and mirror), each with three two-line
# boxes. The width of every line is measured on the same cairo device, at the
# same 8 pt floor, and the column extents printed at the end are copied into
# R17_figures.R. The heights are arithmetic and are derived in R17.
#
# WIDTH. strwidth() on the same cairo device, same font family, same point size
# gives the width of each line in inches, which is compared against the physical
# width of a column.
#
# PANEL_W is the full 7.40 in the panel gets after Fig 1 was composed with
# cowplot::plot_grid(align = "none"). It used to be much narrower when panel a
# was stacked with patchwork's `/`, which aligned it to rows b-e and left its
# plotting region at 4.15 in; the columns were sized against that region and
# every box overflowed. A line that fits here fits in the render: the widths
# this prints are ~3-6% larger than the same strings measure on the rendered PNG.
#
# To re-measure after a layout change, edit `cands` and run this file. The
# printed COL_X is the source of truth; the copy in R17 is the thing to check
# against it.
PANEL_W <- 7.40

cairo_pdf(tempfile(fileext = ".pdf"), width = PANEL_W, height = 2)
par(family = "sans", ps = 12)

# Only the labels that are actually in the panel, grouped by the column they
# sit in. Alternatives can be added while choosing; they must be removed before
# the printed COL_X is copied, or the column will be sized to a line that is not
# there.
cands <- list(
  "1" = c("1 · Inputs and modules", "HF and HCC networks",
          "177 failing LV · 148 tumours", "Networks + 15 Tier A sets",
          "→ HF purple; HCC magenta", "15 Tier A sets", "explicit membership"),
  "2" = c("2 · Derived objects", "All named modules", "23 HF; 20 HCC",
          "TCS and module score", "mean z; kME-weighted",
          "ssGSEA per two-arm cohort", "pathway scores per sample"),
  "3" = c("3 · Separate analyses", "Module preservation",
          "both directions; Zsummary", "Score–trait association",
          "TF correlation; survival", "Hedges g → mirror",
          "opp. sign; |g| ≥ 0.20; CI≠ 0")
)

size <- 8 / (72.27 / 25.4)
cex  <- size * (72.27 / 25.4) / 12
cat(sprintf("\ngeom_text size = %.1f  (%.1f pt);  panel %.2f in wide\n",
            size, size * 72.27 / 25.4, PANEL_W))
for (k in names(cands)) {
  cat(sprintf("\n--- column %s candidates ---\n", k))
  for (s in cands[[k]])
    cat(sprintf("  %-34s %5.2f in  (%d chars)\n", s,
                strwidth(s, units = "inches", cex = cex), nchar(s)))
}

# Longest line each column must hold, plus padding, solved against the panel
# width. 0.40 in of padding inside a box (0.20 each side) is the minimum that
# does not read as text touching the border; this is exactly the complaint the
# 4.15 in panel produced, where column 3 was 1.47 in for a 1.49 in line and the
# clearance was 3 px.
LONG_IN <- vapply(names(cands), function(k)
  max(strwidth(cands[[k]], units = "inches", cex = cex)), numeric(1))
PAD_IN  <- 0.40
BOX_IN  <- LONG_IN + PAD_IN
MARG_IN <- 0.15                       # page margin at each end of the panel
GAP_IN  <- (PANEL_W - 2 * MARG_IN - sum(BOX_IN)) / 2
u_per_in <- 100 / PANEL_W
m <- MARG_IN * u_per_in
g <- GAP_IN * u_per_in
b <- BOX_IN * u_per_in
cat(sprintf("\nlongest line per column: %.2f / %.2f / %.2f in\n",
            LONG_IN[1], LONG_IN[2], LONG_IN[3]))
cat(sprintf("box width (line + %.2f pad): %.2f / %.2f / %.2f in\n",
            PAD_IN, BOX_IN[1], BOX_IN[2], BOX_IN[3]))
cat(sprintf("gap between columns: %.2f in each   page margin: %.2f in each\n",
            GAP_IN, MARG_IN))
x <- c(m,
       m + b[1],
       m + b[1] + g,
       m + b[1] + g + b[2],
       m + b[1] + g + b[2] + g,
       m + b[1] + g + b[2] + g + b[3])
cat(sprintf("\nCOL_X <- list(c(%.1f, %.1f), c(%.1f, %.1f), c(%.1f, %.1f))\n",
            x[1], x[2], x[3], x[4], x[5], x[6]))
dev.off()
