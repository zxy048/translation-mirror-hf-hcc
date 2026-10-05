# =============================================================================
# R00_io.R -- shared GEO / expression-matrix I/O helpers
# -----------------------------------------------------------------------------
# Sourced by R00_config.R, so every Rxx script gets these. Extracted verbatim
# from R01_build_universe.R (where they were first written and validated) so
# R05 and R01 cannot drift apart -- v1's failure mode was three scripts with
# three different regexes and three different filters.
# =============================================================================

#' Parse a GEO series-matrix file into (expr matrix, sample annotation).
#' Annotation rows are the !Sample_characteristics_ch1 rows, keyed by the text
#' before the first colon ("disease", "Sex", "age", ...).
parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt"); lines <- readLines(con, warn = FALSE); close(con)
  tb <- grep("^!series_matrix_table_begin", lines)
  te <- grep("^!series_matrix_table_end", lines)
  if (!length(tb) || !length(te)) stop("no series matrix table in ", path)

  hdr <- strsplit(lines[tb + 1], "\t", fixed = TRUE)[[1]]
  hdr <- gsub('^"|"$', "", hdr)
  body <- lines[(tb + 2):(te - 1)]
  parts <- strsplit(body, "\t", fixed = TRUE)

  # Row IDs are quoted in some series matrices ("1007_s_at") and bare in others.
  # The header is unquoted above; the IDs must be too, or every probe->symbol
  # lookup silently misses and the cohort collapses to zero genes.
  ids <- gsub('^"|"$', "", vapply(parts, `[`, character(1), 1))
  mat <- matrix(NA_real_, nrow = length(parts), ncol = length(hdr) - 1,
                dimnames = list(ids, hdr[-1]))
  for (i in seq_along(parts)) {
    v <- suppressWarnings(as.numeric(parts[[i]][-1]))
    length(v) <- ncol(mat)
    mat[i, ] <- v
  }

  ch <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
  char <- lapply(ch, function(l) {
    v <- strsplit(l, "\t", fixed = TRUE)[[1]][-1]
    gsub('^"|"$', "", v)
  })
  names(char) <- vapply(char, function(v) sub(":.*$", "", v[1]), character(1))
  char <- lapply(char, function(v) sub("^[^:]*:\\s*", "", v))

  gsm <- hdr[-1]
  ann <- data.frame(gsm = gsm, stringsAsFactors = FALSE)
  for (k in names(char)) {
    v <- char[[k]]
    if (length(v) == nrow(ann)) ann[[make.names(k)]] <- v
  }
  list(expr = mat, ann = ann)
}

#' Read a GEO .soft platform table into a data.frame.
parse_soft <- function(path) {
  con <- gzfile(path, "rt"); out <- NULL; inside <- FALSE
  repeat {
    l <- readLines(con, n = 20000, warn = FALSE)
    if (!length(l)) break
    if (!inside) {
      i <- grep("^!platform_table_begin", l)
      if (!length(i)) next
      inside <- TRUE; l <- l[(i[1] + 1):length(l)]
    }
    j <- grep("^!platform_table_end", l)
    if (length(j)) { out <- c(out, l[seq_len(j[1] - 1)]); break }
    out <- c(out, l)
  }
  close(con)
  if (is.null(out) || length(out) < 2) stop("no platform table in ", path)
  hdr <- strsplit(out[1], "\t", fixed = TRUE)[[1]]
  rows <- strsplit(out[-1], "\t", fixed = TRUE)
  rows <- lapply(rows, function(r) { length(r) <- length(hdr); r })
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(df) <- hdr
  df
}

#' Collapse duplicate symbols, keeping the row with the highest mean expression.
#' Deterministic; avoids summing probes to the same gene.
collapse_by_mean <- function(mat, symbols, min_mean = -Inf) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "---"
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  m <- rowMeans(mat, na.rm = TRUE)
  ord <- order(symbols, -m)
  mat <- mat[ord, , drop = FALSE]; symbols <- symbols[ord]
  dup <- duplicated(symbols)
  dropped <- data.frame(symbol = symbols[dup], mean_expr = m[ord][dup],
                        stringsAsFactors = FALSE)
  mat <- mat[!dup, , drop = FALSE]
  rownames(mat) <- symbols[!dup]
  list(expr = mat, dropped = dropped)
}

#' Probe -> SYMBOL lookup from a .soft platform table, using whichever column
#' holds symbols. Column names differ between platforms ("Symbol", "Gene Symbol",
#' "Gene symbol", "GENE_SYMBOL"), so try them in order and fail loudly if none
#' is present -- silently returning all-NA would drop every gene in the cohort.
soft_probe2symbol <- function(soft_df, id_col = "ID") {
  cand <- c("Symbol", "Gene Symbol", "Gene symbol", "GENE_SYMBOL", "GENE")
  hit <- cand[cand %in% names(soft_df)][1]
  if (is.na(hit))
    stop("no symbol column in platform table; columns are: ",
         paste(names(soft_df), collapse = ", "))
  setNames(soft_df[[hit]], soft_df[[id_col]])
}

#' ENSG (with or without version suffix) -> SYMBOL, via org.Hs.eg.db.
ensg2symbol <- function(ensg) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  key <- sub("\\..*$", "", ensg)          # strip Ensembl version suffix
  suppressWarnings(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = unique(key), column = "SYMBOL",
    keytype = "ENSEMBL", multiVals = "first")[key])
}

#' Detection filter, declared per technology. Applied identically to discovery
#' and validation cohorts of the same technology.
filter_array <- function(mat) rowSums(mat > stats::median(mat, na.rm = TRUE),
                                      na.rm = TRUE) >= 0.5 * ncol(mat)
filter_counts <- function(mat, min_count = 10) rowSums(mat >= min_count,
                                                       na.rm = TRUE) >= 0.5 * ncol(mat)
filter_nonzero <- function(mat, frac = 0.5) rowSums(mat > 0, na.rm = TRUE) >= frac * ncol(mat)
