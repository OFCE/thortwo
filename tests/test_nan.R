## A value that is not a number must never come back as a solution.
##
## The compiled backends measured a block's residual with Eigen's
## lpNorm<Infinity>(), which skips a NaN unless it is the first element:
## (1, NaN, 3) gave 3. A block whose residual was NaN in any equation but its
## first was therefore not seen to have failed, and the solver could report
## convergence, and a small residual, while returning NaN for some variables.
## The R backend never had the problem, because R's max() propagates NaN.
##
## Every case here puts the NaN somewhere other than the first equation of
## its block, on all three backends.
##
## Usage:  Rscript tests/test_nan.R

source("tests/helper.R")

work <- new_tmpdir("thortwo_nan_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

write_model <- function(name, endo, exo, eqs) {
  f <- file.path(work, paste0(name, ".txt"))
  writeLines(c("endogenous :", endo, "exogenous :", exo, "coefficients :", "",
               "equations :", eqs), f)
  f
}

## 1. A square root sitting exactly at zero, in a block solved by Newton: the
##    derivative of `s = ss^0.5` with respect to ss is infinite there. On the
##    sparse backend this model used to "solve" with s = NaN.
root <- write_model("root", "a, b, s, ss, y", "x",
                    c("y = 0.5 * y + x", "a = x * y", "b = x * y",
                      "ss = a^2 + b^2", "s = ss^0.5"))
d_root <- data.frame(year = 2000:2004, x = 0, a = 0, b = 0, s = 0, ss = 0, y = 0)

## 2. The log of a negative number, in the second equation of its block.
logm <- write_model("logm", "p, q, r", "x, z",
                    c("p = x + 0.5 * q", "q = log(z) + 0.5 * p", "r = p + q"))
d_log <- data.frame(year = 2000:2004, x = 1, z = -1, p = 1, q = 1, r = 2)

for (backend in c("sparse", "dense-cpp", "sparse-r")) {
  cat("\n=== ", backend, " ===\n", sep = "")

  m <- thor_model(paste0("nan_root_", gsub("-", "", backend)), root, backend = backend,
                  workdir = file.path(work, paste0("root", backend)), verbose = FALSE)
  out <- tryCatch(thor_solve(m, 2001, 2004, d_root, "year", verbose = FALSE),
                  error = function(e) conditionMessage(e))
  ok("a root at zero is solved with real numbers, or refused",
     is.character(out) || !anyNA(out[out$year >= 2001, m@vars$endo]),
     if (is.character(out)) "refused" else "solved")

  m <- thor_model(paste0("nan_log_", gsub("-", "", backend)), logm, backend = backend,
                  workdir = file.path(work, paste0("log", backend)), verbose = FALSE)
  ## (the R backend also warns "NaNs produced", which is the point, not news)
  res <- suppressWarnings(model_residuals(m, d_log, periods = 2001, index_time = "year"))
  ok("a NaN residual is reported as NaN", any(is.nan(res)),
     paste(sprintf("%s %s", colnames(res), format(res[1, ])), collapse = ", "))
  out <- suppressWarnings(tryCatch(thor_solve(m, 2001, 2004, d_log, "year", verbose = FALSE),
                                   error = function(e) conditionMessage(e)))
  ok("a model that evaluates to NaN is refused",
     is.character(out) || !anyNA(out[out$year >= 2001, m@vars$endo]),
     if (is.character(out)) "refused" else "solved")
}

cat("\nPASS\n")
