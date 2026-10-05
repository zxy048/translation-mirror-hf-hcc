# 00_probe7.R -- verify the GPL11532 annotation route for GSE57338.
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

cat("=== hugene11sttranscriptcluster.db ===\n")
ok <- requireNamespace("hugene11sttranscriptcluster.db", quietly = TRUE)
cat("available:", ok, "\n")
if (ok) {
  library(hugene11sttranscriptcluster.db)
  k <- AnnotationDbi::keys(hugene11sttranscriptcluster.db, keytype = "PROBEID")
  cat("total PROBEIDs:", length(k), "\n")
  cat("head:", paste(head(k, 5), collapse = ", "), "\n")

  # do the matrix IDs match?
  con <- gzfile(COHORTS$GSE57338$file, "rt"); l <- readLines(con, n = 78); close(con)
  ids <- sub('^"([^"]+)".*', "\\1", l[74:78])
  cat("matrix ID_REF head:", paste(ids, collapse = ", "), "\n")
  cat("matrix IDs present in db:", sum(ids %in% k), "/", length(ids), "\n")

  # symbol coverage on a sample of 2000 probes
  set.seed(STATS$SEED)
  s <- sample(k, 2000)
  sym <- suppressWarnings(AnnotationDbi::mapIds(
    hugene11sttranscriptcluster.db, keys = s, column = "SYMBOL",
    keytype = "PROBEID", multiVals = "first"))
  cat("symbol coverage in random 2000:", sum(!is.na(sym)), "->",
      sprintf("%.1f%%", 100 * mean(!is.na(sym))), "\n")
  cat("examples:", paste(na.omit(sym)[1:10], collapse = ", "), "\n")

  # do the Tier A translation genes survive?
  A <- resolve_sets(SETS_TIER_A)
  target <- unique(unlist(A))
  allsym <- suppressWarnings(AnnotationDbi::mapIds(
    hugene11sttranscriptcluster.db, keys = k, column = "SYMBOL",
    keytype = "PROBEID", multiVals = "first"))
  present <- target[target %in% allsym]
  cat("\nTier A genes on GPL11532:", length(present), "/", length(target),
      sprintf(" (%.1f%%)", 100 * length(present) / length(target)), "\n")
  miss <- setdiff(target, present)
  if (length(miss)) cat("MISSING:", paste(miss, collapse = ", "), "\n")
}

cat("\n=== generic .soft platform-table parser sanity check ===\n")
parse_soft <- function(path) {
  con <- gzfile(path, "rt"); out <- NULL; inside <- FALSE
  repeat {
    l <- readLines(con, n = 5000, warn = FALSE)
    if (!length(l)) break
    if (!inside) {
      i <- grep("^!platform_table_begin", l)
      if (length(i)) { inside <- TRUE; l <- l[(i[1] + 1):length(l)] } else next
    }
    j <- grep("^!platform_table_end", l)
    if (length(j)) { out <- c(out, l[seq_len(j[1] - 1)]); break }
    out <- c(out, l)
  }
  close(con)
  if (is.null(out) || length(out) < 2) return(NULL)
  hdr <- strsplit(out[1], "\t")[[1]]
  body <- out[-1]
  # keep only as many fields as the header
  rows <- strsplit(body, "\t", fixed = TRUE)
  rows <- lapply(rows, function(r) { length(r) <- length(hdr); r })
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(df) <- hdr
  df
}
for (f in c("GPL10558.soft.gz", "GPL3921.soft.gz", "GPL571.soft.gz")) {
  d <- tryCatch(parse_soft(file.path("D:/R_projects", f)),
                error = function(e) { cat(f, "ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(d)) {
    cat(sprintf("  %-20s %d rows x %d cols; has SYMBOL col: %s\n", f, nrow(d), ncol(d),
                any(grepl("symbol", names(d), ignore.case = TRUE))))
    cat("     cols:", paste(head(names(d), 14), collapse = ", "), "\n")
  }
}
