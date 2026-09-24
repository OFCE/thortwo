## Two models in one session (D3).
##
## tresthor's classic backend called `Rcpp::sourceCpp(model@rcpp_source)` on
## *every* solve, which loaded `Rcpp_solver` into the global environment. With
## two classic models live at once, the second silently overwrote the first,
## and the first was then solved with the second's code -- quietly, with no
## error, producing a plausible wrong answer.
##
## Every generated file here exports the same two names
## (`thor_cpp_solve`, `thor_cpp_residuals`), so the failure mode is still
## available; what prevents it is that each model is loaded into its own
## environment, keyed by the hash of its code. This test would fail on
## tresthor.
##
## Usage:  Rscript tests/test_isolation.R

source("tests/helper.R")

work <- new_tmpdir("thortwo_isolation_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

data_opale <- readRDS("inst/Opale/donnees_opale.rds")
coeffs     <- readRDS("inst/Opale/coefficients_opale.rds")
data_opale <- add_coeffs(coeffs, data_opale, pos.coeff.name = 2, pos.coeff.value = 1)
dates <- as.character(data_opale$date); n <- length(dates)
from <- dates[n - 9]; to <- dates[n]

## Two genuinely different models, so that solving one with the other's code
## would give a different answer rather than the same one by luck. The second
## is Opale without its decomposition: same equations, one block, different
## generated code.
a <- thor_model("iso_a", "inst/models/opale.txt", backend = "dense-cpp",
                workdir = file.path(work, "a"), verbose = FALSE)
b <- thor_model("iso_b", "inst/models/opale.txt", backend = "dense-cpp",
                decompose = FALSE, workdir = file.path(work, "b"), verbose = FALSE)

ok("the two models really are different",
   !identical(a@generated$hash, b@generated$hash),
   sprintf("%d blocks vs %d", length(a@blocks), length(b@blocks)))

## Solve a, then b, then a again. If loading b clobbered a, the third solve
## disagrees with the first.
res_a1 <- thor_solve(a, from, to, data_opale, "date", verbose = FALSE)
res_b  <- thor_solve(b, from, to, data_opale, "date", verbose = FALSE)
res_a2 <- thor_solve(a, from, to, data_opale, "date", verbose = FALSE)

rows <- match(dates[(n - 9):n], dates)
ok("solving a second model does not clobber the first",
   identical(endo_matrix(res_a1, a, rows), endo_matrix(res_a2, a, rows)))

## And the two models must agree on the answer, since they are the same
## equations solved two ways: that is what shows the isolation is real rather
## than both of them silently running the same code.
check_agree(endo_matrix(res_a1, a, rows), endo_matrix(res_b, b, rows),
            label = "decomposed and undecomposed agree")

## Residuals are evaluated through the same registry, so they must be isolated
## too.
ra <- model_residuals(a, res_a1, periods = dates[(n - 9):n], index_time = "date")
rb <- model_residuals(b, res_b,  periods = dates[(n - 9):n], index_time = "date")
ok("residuals come from the right model",
   ncol(ra) == length(a@blocks) && ncol(rb) == length(b@blocks) &&
     max(ra) < 1e-6 && max(rb) < 1e-6,
   sprintf("%d blocks / %d blocks", ncol(ra), ncol(rb)))

## The R backend shares the registry, so mixing backends in one session must
## work too.
cat("\n=== mixed backends ===\n")
r <- thor_model("iso_r", "inst/models/opale.txt", backend = "sparse-r",
                workdir = file.path(work, "r"), verbose = FALSE)
res_r  <- thor_solve(r, from, to, data_opale, "date", verbose = FALSE)
res_a3 <- thor_solve(a, from, to, data_opale, "date", verbose = FALSE)
ok("an R-backend model does not clobber a compiled one",
   identical(endo_matrix(res_a1, a, rows), endo_matrix(res_a3, a, rows)))
check_agree(endo_matrix(res_a1, a, rows), endo_matrix(res_r, r, rows),
            label = "dense-cpp and sparse-r agree")

cat("\nPASS\n")
