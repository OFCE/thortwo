## The timings display of thor_model() and thor_solve().
##
## On by default, independent of `verbose` for the build, and for the solve
## shown only when `verbose` is off (the verbose summary already has the time).
## `timings = FALSE` and options(thortwo.timings = FALSE) both switch it off.
##
## Usage:  Rscript tests/test_timings.R

source("tests/helper.R")
options(thortwo.model.cache = TRUE)    # helper.R turns it off; a cache hit is tested here

work <- new_tmpdir("thortwo_timings_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

f <- file.path(work, "tiny.txt")
writeLines(c("endogenous :", "y, c", "exogenous :", "g", "coefficients :", "",
             "equations :", "y = c + g", "c = 0.6*y + 0.2*lag(c,1)"), f)
d <- data.frame(date = 1:6, y = 1, c = 1, g = seq(1, 2, length.out = 6))
line <- function(out) grep("^Timings:", out, value = TRUE)

for (backend in c("sparse", "sparse-r")) {
  cat("\n---- ", backend, " ----\n", sep = "")
  nm <- paste0("tim_", sub("-", "", backend))
  store <- file.path(work, backend)

  out <- capture.output(m <- thor_model(nm, f, backend = backend, cache = store, verbose = FALSE))
  ok("a quiet build prints one timings line",
     length(line(out)) == 1L && grepl("build .* s, compile .* s, total", line(out)), line(out))
  tm <- m@meta$timings
  ok("and keeps them in the model",
     identical(names(tm), c("build", "compile")) && all(tm >= 0) && isFALSE(attr(tm, "from_cache")))

  out <- capture.output(m2 <- suppressMessages(thor_model(nm, f, backend = backend, cache = store, verbose = FALSE)))
  ok("a cache hit says so", grepl("from the cache", line(out)) && isTRUE(attr(m2@meta$timings, "from_cache")),
     line(out))

  out <- capture.output(m3 <- thor_model(nm, f, backend = backend, cache = FALSE, compile = FALSE, verbose = FALSE))
  ok("no compile time when nothing is compiled",
     is.na(m3@meta$timings[["compile"]]) && !grepl("compile", line(out)), line(out))

  out <- capture.output(m4 <- thor_model(nm, f, backend = backend, cache = FALSE, verbose = TRUE))
  ok("a verbose build prints it once", length(line(out)) == 1L)

  out <- capture.output(m4 <- thor_model(nm, f, backend = backend, cache = FALSE, verbose = FALSE, timings = FALSE))
  ok("timings = FALSE prints nothing", length(out) == 0L)

  out <- capture.output(r <- thor_solve(m, 2, 6, d, verbose = FALSE))
  ok("a quiet solve prints its time", length(line(out)) == 1L && grepl("solve .* s \\(5 periods\\)", line(out)), line(out))
  out <- capture.output(r <- thor_solve(m, 2, 6, d, verbose = TRUE))
  ok("a verbose solve does not repeat it", length(line(out)) == 0L && any(grepl("^Solved 5 periods in", out)))
  out <- capture.output(r <- thor_solve(m, 2, 6, d, verbose = FALSE, timings = FALSE))
  ok("timings = FALSE prints nothing", length(out) == 0L)

  old <- options(thortwo.timings = FALSE)
  out <- capture.output({ m4 <- thor_model(nm, f, backend = backend, cache = FALSE, verbose = FALSE)
                          r <- thor_solve(m, 2, 6, d, verbose = FALSE) })
  options(old)
  ok("options(thortwo.timings = FALSE) switches both off", length(out) == 0L)
}

fs <- thortwo:::format_seconds
ok("durations read well", identical(c(fs(0.5239), fs(38.44), fs(59.97), fs(125.4)),
                                    c("0.52 s", "38.4 s", "1 min 00 s", "2 min 05 s")))

cat("\nPASS\n")
