## The acceptance suite.
##
## Usage:  Rscript tests/run_all.R          # everything except ThreeME 8x8
##         Rscript tests/run_all.R full     # everything, ~6 minutes cold
##
## These are scripts rather than testthat files on purpose: what they check is
## whether a real model solves to a real answer in a plausible amount of time,
## and the timings and residual magnitudes are as much the output as the
## pass/fail is.

full <- identical(commandArgs(trailingOnly = TRUE)[1], "full")

scripts <- list(
  c("tests/test_roundtrip.R", ""),
  c("tests/test_opale.R",     "all"),
  c("tests/test_isolation.R", ""),
  c("tests/test_cache.R",     ""),
  c("tests/test_threeme.R",   "4x4 all")
)
if (full) scripts <- c(scripts, list(c("tests/test_threeme.R", "8x8 all")))

results <- character(0)
for (s in scripts) {
  label <- paste(basename(s[1]), s[2])
  cat("\n", strrep("=", 70), "\n", label, "\n", strrep("=", 70), "\n", sep = "")
  t <- system.time(
    code <- system2("Rscript", c(s[1], strsplit(s[2], " ")[[1]]))
  )[["elapsed"]]
  results <- c(results, sprintf("%-28s %s  (%.0f s)", label,
                                if (code == 0) "PASS" else "FAIL", t))
  if (code != 0) {
    cat("\n", paste(results, collapse = "\n"), "\n", sep = "")
    stop(label, " failed", call. = FALSE)
  }
}

cat("\n", strrep("=", 70), "\n", paste(results, collapse = "\n"), "\n", sep = "")
if (!full) cat("\n(ThreeME 8x8 not run; add `full` to include it)\n")
