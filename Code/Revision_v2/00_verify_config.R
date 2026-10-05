# 00_verify_config.R -- prove R00_config.R is internally valid BEFORE any
# expensive step runs: every frozen set name must resolve, and the tier
# arithmetic must match what we claim in the manuscript text.
source("D:/R_projects/revision_analysis/Code/Revision_v2/R00_config.R")

fails <- 0
chk <- function(cond, msg) {
  cat(if (isTRUE(cond)) "  PASS  " else "  FAIL  ", msg, "\n", sep = "")
  if (!isTRUE(cond)) fails <<- fails + 1
}

cat("=== 1. Tier A resolves ===\n")
A <- tryCatch(resolve_sets(SETS_TIER_A), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
chk(!is.null(A) && length(A) == length(SETS_TIER_A),
    sprintf("Tier A: %d/%d sets resolved", if (is.null(A)) 0L else length(A), length(SETS_TIER_A)))

cat("\n=== 2. Tier B resolves ===\n")
B <- tryCatch(resolve_sets(c(SETS_TIER_A, SETS_TIER_B_EXTRA)), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
chk(!is.null(B) && length(B) == length(SETS_TIER_A) + length(SETS_TIER_B_EXTRA),
    sprintf("Tier B: %d sets resolved", if (is.null(B)) 0L else length(B)))

cat("\n=== 3. Tier C (v1's 33) resolves ===\n")
C <- tryCatch(resolve_sets(SETS_TIER_C_V1), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
chk(!is.null(C) && length(C) == 33L,
    sprintf("Tier C: %d/33 sets resolved", if (is.null(C)) 0L else length(C)))

cat("\n=== 4. signature sets resolve ===\n")
sig <- unique(unlist(SIGNATURE_SETS))
S <- tryCatch(resolve_sets(sig), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
chk(!is.null(S) && length(S) == length(sig),
    sprintf("signature sets: %d/%d resolved", if (is.null(S)) 0L else length(S), length(sig)))

cat("\n=== 5. the audit claim: how much of v1's 33 is NOT translation? ===\n")
excl <- unique(unlist(SETS_EXCLUDED))
overlap <- intersect(SETS_TIER_C_V1, excl)
cat("  v1 sets total            :", length(SETS_TIER_C_V1), "\n")
cat("  flagged non-translation  :", length(overlap), "\n")
cat("  remaining 'core'         :", length(setdiff(SETS_TIER_C_V1, excl)), "\n")
cat("  by bucket:\n")
for (b in names(SETS_EXCLUDED))
  cat(sprintf("    %-20s %d\n", b, length(intersect(SETS_TIER_C_V1, SETS_EXCLUDED[[b]]))))
chk(length(setdiff(SETS_TIER_C_V1, excl)) == 13L,
    "exactly 13 of v1's 33 are genuine core translation machinery")
chk(length(overlap) == 20L, "exactly 20 of v1's 33 are off-topic by our a priori buckets")
chk(setequal(setdiff(SETS_TIER_C_V1, excl), intersect(SETS_TIER_A, SETS_TIER_C_V1)),
    "core-from-v1 == Tier A intersect Tier C")

cat("\n=== 6. set sizes as they will be used ===\n")
for (nm in SETS_TIER_A) cat(sprintf("  %-78s n=%d\n", substr(nm, 1, 78), length(A[[nm]])))

cat("\n=== 7. stats helpers ===\n")
set.seed(1); x <- rnorm(30, 1); y <- rnorm(25, 0)
cd <- cohens_d(x, y)
chk(is.finite(cd["d"]) && is.finite(cd["g"]) && cd["g"] < cd["d"],
    sprintf("cohens_d: d=%.3f g=%.3f (g shrunk toward 0)", cd["d"], cd["g"]))
chk(abs(n_effective(diag(5)) - 5) < 1e-9, "n_effective(identity(5)) == 5")
chk(abs(n_effective(matrix(1, 5, 5)) - 1) < 1e-9, "n_effective(all-ones(5)) == 1")

cat("\n=== 8. jaccard_dist ===\n")
jd <- jaccard_dist(list(a = c("X","Y","Z"), b = c("X","Y","Z"), c = c("P","Q")))
chk(abs(as.matrix(jd)[1,2] - 0) < 1e-9, "identical sets -> distance 0")
chk(abs(as.matrix(jd)[1,3] - 1) < 1e-9, "disjoint sets -> distance 1")

cat("\n===========================================\n")
cat(if (fails == 0) "ALL CHECKS PASSED\n" else sprintf("%d CHECK(S) FAILED\n", fails))
cat("===========================================\n")
if (fails > 0) quit(status = 1)
