## thor_solve(reuse_jacobian = TRUE): keep the factorised jacobian from one
## Newton iteration to the next within a period, instead of recomputing and
## refactorising it every time. The default on the compiled backends.
##
## On a large model the factorisation is nearly the whole cost of an
## iteration -- 1.36 s on the heart of ThreeME 29x33, against milliseconds for
## everything else -- so trading it for a few more, cheaper iterations made
## that solve 2.8 times faster (208 s to 74 s). The option changes how the
## solution is reached, not what it is, so that is what is tested:
##
##   1. off, nothing changes: the iteration count is the one it always was;
##   2. on, the solution is Newton's and the residuals are as small. This is
##      not automatic: steps with old factors converge linearly and, stopped
##      at the ordinary tolerance, left residuals of 3e-6 on Opale;
##   3. on, it really does take more iterations -- otherwise it is not on;
##   4. it combines with `sequential = TRUE`;
##   5. the default is on for the compiled backends and off for sparse-r,
##      where more iterations mean more equations evaluated in R: slower.
##
## Usage:  Rscript tests/test_reuse.R

source("tests/helper.R")

work <- new_tmpdir("thortwo_reuse_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

data_opale <- readRDS("inst/Opale/donnees_opale.rds")
coeffs     <- readRDS("inst/Opale/coefficients_opale.rds")
data_opale <- add_coeffs(coeffs, data_opale, pos.coeff.name = 2, pos.coeff.value = 1)
dates <- as.character(data_opale$date); n <- length(dates)
from <- dates[n - 39]; to <- dates[n]; rows <- (n - 39):n

for (backend in c("sparse", "dense-cpp", "sparse-r")) {
  cat("\n=== Opale / ", backend, " ===\n", sep = "")
  m <- thor_model(paste0("ru_", gsub("-", "", backend)), "inst/models/opale.txt", backend = backend,
                  workdir = file.path(work, backend), verbose = FALSE)
  off <- thor_solve(m, from, to, data_opale, "date", verbose = FALSE, diagnostics = TRUE,
                    reuse_jacobian = FALSE)
  on  <- thor_solve(m, from, to, data_opale, "date", verbose = FALSE, diagnostics = TRUE,
                    reuse_jacobian = TRUE)
  dflt <- thor_solve(m, from, to, data_opale, "date", verbose = FALSE, diagnostics = TRUE)
  ok(if (backend == "sparse-r") "the default is off" else "the default is on",
     identical(endo_matrix(dflt, m, rows),
               endo_matrix(if (backend == "sparse-r") off else on, m, rows)))

  ok("off: the iteration count is unchanged", sum(attr(off, "iterations")) == 459L,
     paste(sum(attr(off, "iterations")), "iterations"))
  ok("on: more, cheaper iterations",
     sum(attr(on, "iterations")) > 2 * sum(attr(off, "iterations")),
     paste(sum(attr(on, "iterations")), "iterations"))
  check_converged(on)
  check_residuals(m, on, dates[rows], "date", tol = 1e-8)
  check_agree(endo_matrix(off, m, rows), endo_matrix(on, m, rows), label = "on and off agree")
  worst <- max(rel_diff(endo_matrix(off, m, rows), endo_matrix(on, m, rows))$rel, na.rm = TRUE)
  ok("... on every variable", worst < 1e-8, sprintf("largest relative difference %.2g", worst))
}

cat("\n=== ThreeME 4x4 / sparse, with and without sequential ===\n")
data_3me <- readRDS("inst/ThreeME/data3me_4x4.rds")
rows3 <- match(2016:2050, data_3me$year)
for (sq in c(FALSE, TRUE)) {
  m <- thor_model(paste0("ru_3me_", sq), "inst/ThreeME/threeme_4x4_thor.txt", sequential = sq,
                  workdir = file.path(work, paste0("t", sq)), verbose = FALSE)
  off <- thor_solve(m, 2016, 2050, data_3me, "year", verbose = FALSE, reuse_jacobian = FALSE)
  on  <- thor_solve(m, 2016, 2050, data_3me, "year", verbose = FALSE, diagnostics = TRUE,
                    reuse_jacobian = TRUE)
  check_converged(on, if (sq) "converged, with sequential" else "converged")
  check_residuals(m, on, 2016:2050, "year")
  check_agree(as.matrix(off[rows3, m@vars$endo]), as.matrix(on[rows3, m@vars$endo]),
              label = "on and off agree")
}

cat("\nPASS\n")
