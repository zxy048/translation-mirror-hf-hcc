# Measure the Fig 1a schematic labels instead of guessing at the font size.
#
# Two sizing questions govern panel a, and both were answered wrongly the first
# time its labels were raised to the 8 pt floor. This script answers the one that
# can be answered by measurement: the width of a line, against the width of the
# column that has to hold it. The other -- the height of a box against the text
# it holds -- is arithmetic, and is derived in R17_figures.R.
#
# WIDTH. strwidth() on the same cairo device, same font family, same point size
# gives the width of each line in inches, which is compared against the physical
# width of a column.
#
# Box geometry, from R17_figures.R: x runs 0-100 across the panel and the three
# columns span 18, 23 and 26 units. Those extents are themselves derived from
# this script's output rather than the other way round -- the columns are sized
# to the longest line each has to hold, so a new or lengthened label means
# re-running this and copying the COL_X it prints into R17_figures.R.
#
# PANEL_W is NOT the 7.5 in figure width. Panel a is one band of a five-panel
# patchwork and its plotting region is narrower than the page. This constant has
# now been wrong twice -- first 7.5, then 5.07 -- and each time the probe
# reported that labels fitted and the render disagreed, because the columns were
# narrower than the probe believed. It is therefore measured rather than
# reasoned about: the rendered Fig 1 PNG (300 dpi) puts three box fills at known
# x in data units, and all three pairs agree on 12.8676 px per unit (box 1 at
# unit 16 -> px 694.5, box 3 at unit 50 -> px 1132.0, box 5 at unit 84 ->
# px 1569.5). 100 units is thus 1286.8 px = 4.29 in, so equal thirds would leave
# every column 1.29 in and the widest label here (1.48 in) would overflow --
# which is why the columns are content-sized instead.
#
# To re-measure after a layout change: locate the fills of the three boxes by
# their fill colours (#DCE7F2, #EDEDED, #FBEFD0) in the top band of the PNG and
# take the distance between their centres over the distance between their data
# units.
#
# The widths this prints are ~3-6% larger than the same strings measure in the
# render, i.e. err on the side of caution. A line that fits here fits there.
PANEL_W <- 7.40
# The column extents are not equal thirds: each column is sized to the longest
# line it has to hold, which is what this script measures. Equal thirds would
# leave column 3 at 1.29 in, and the eigengene box's second line needs 1.48 in.
# There used to be a COL_U here holding those extents by hand, and it is gone on
# purpose. Nothing read it -- only the COL_X the script *prints* at the end is
# used -- so it was free to drift from R17's COL_X while its comment claimed they
# matched, and it had: hand-editing around the 4.15 in panel left it at
# c(2,21) c(34,59) against R17's c(2,20) c(35,58). The printed COL_X is computed
# from LONG_IN and PAD_IN below, so it cannot drift; that output is the source of
# truth, and the copy in R17 is the thing to check against it.

cairo_pdf(tempfile(fileext = ".pdf"), width = PANEL_W, height = 2)
par(family = "sans", ps = 12)

# Candidates for the shortened labels, grouped by the column they would sit in.
# Column 3 carries three boxes and is the binding constraint on height, so it is
# the column asked to give up a line.
cands <- list(
  "1" = c("HF (GSE57338)", "HCC (GSE141198)", "HF network", "HCC network",
          "177 failing LV", "23 named modules", "148 tumours", "20 named modules"),
  "2" = c("A priori pathway tiers", "15 Tier A sets", "explicit membership",
          "(KEGG + Reactome)", "ssGSEA per network",
          "HF purple (153 genes)", "HCC magenta (146 genes)",
          "→ designated module:"),
  "3" = c("TCS — disease-agnostic", "mean z, 15 Tier A sets",
          "Module eigengene score", "kME-weighted; both directions",
          "kME-weighted, both directions", "both directions",
          "kME-weighted mean z", "Preservation test and",
          "Hedges g comparison")
)

size <- 8 / (72.27 / 25.4)
cex  <- size * (72.27 / 25.4) / 12
cat(sprintf("
geom_text size = %.1f  (%.1f pt);  panel %.2f in wide
",
            size, size * 72.27 / 25.4, PANEL_W))
for (k in names(cands)) {
  cat(sprintf("
--- column %s candidates ---
", k))
  for (s in cands[[k]])
    cat(sprintf("  %-32s %5.2f in  (%d chars)
", s, strwidth(s, units="inches", cex=cex), nchar(s)))
}

# Longest line each column must hold, plus padding, solved against the panel width.
# WID_IN is the widest label per column; the two gaps and the two page margins are
# what is left over. 0.40 in of padding inside a box (0.20 each side) is the
# minimum that does not read as text touching the border -- this is exactly the
# complaint the 4.15 in panel produced, where column 3 was 1.47 in for a 1.49 in
# line and the clearance was 3 px.
LONG_IN <- vapply(names(cands), function(k)
  max(strwidth(cands[[k]], units="inches", cex=cex)), numeric(1))
PAD_IN  <- 0.40
BOX_IN  <- LONG_IN + PAD_IN
MARG_IN <- 0.15                       # page margin at each end of the panel
GAP_IN  <- (PANEL_W - 2*MARG_IN - sum(BOX_IN)) / 2
u <- function(inch) inch / (PANEL_W - 2*MARG_IN) * (100 - 2*MARG_IN/PANEL_W*100) + 0
cat(sprintf("
longest line per column: %.2f / %.2f / %.2f in
", LONG_IN[1], LONG_IN[2], LONG_IN[3]))
cat(sprintf("box width (line + %.2f pad): %.2f / %.2f / %.2f in
", PAD_IN, BOX_IN[1], BOX_IN[2], BOX_IN[3]))
cat(sprintf("gap between columns: %.2f in each   page margin: %.2f in each
", GAP_IN, MARG_IN))
u_per_in <- 100 / PANEL_W
m <- MARG_IN * u_per_in
g <- GAP_IN * u_per_in
b <- BOX_IN * u_per_in
x <- c(m, m + b[1], m + b[1] + g, m + b[1] + g + b[2], m + b[1] + g + b[2] + g, m + b[1] + g + b[2] + g + b[3])
cat(sprintf("
COL_X <- list(c(%.0f, %.0f), c(%.0f, %.0f), c(%.0f, %.0f))
",
            x[1], x[2], x[3], x[4], x[5], x[6]))
dev.off()
