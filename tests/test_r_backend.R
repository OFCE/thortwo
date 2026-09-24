## The pure-R backend, `sparse-r`.
##
## It replaced the dense R backend in 2026-10: same generated R code for the
## equations, but a sparse jacobian solved with Matrix's sparse LU instead of a
## dense one with base solve(). On ThreeME 8x8 that took the solve from nearly
## 7 minutes to 8 s. This checks the three promises that came with it:
##
##   1. it gives the compiled `sparse` backend's answer, in the same number of
##      Newton iterations -- the iteration count is what caught a stale cached
##      factorisation that still converged, only 3x more slowly;
##   2. `backend = "dense-r"` still works, as the old name of `sparse-r`;
##   3. a model saved with the old dense R backend still loads and solves.
##      `old_dense_r_opale.rds` was saved by thortwo before the change.
##
## Usage:  Rscript tests/test_r_backend.R

source("tests/helper.R")

work <- new_tmpdir("thortwo_rbackend_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

data_opale <- readRDS("inst/Opale/donnees_opale.rds")
coeffs     <- readRDS("inst/Opale/coefficients_opale.rds")
data_opale <- add_coeffs(coeffs, data_opale, pos.coeff.name = 2, pos.coeff.value = 1)
dates <- as.character(data_opale$date); n <- length(dates)
from <- dates[n - 39]; to <- dates[n]; rows <- (n - 39):n

## ---- 1. same answer as compiled sparse ---------------------------------------
s <- thor_model("rb_sparse", "inst/models/opale.txt", backend = "sparse",
                workdir = file.path(work, "s"), verbose = FALSE)
r <- thor_model("rb_r", "inst/models/opale.txt", backend = "sparse-r",
                workdir = file.path(work, "r"), verbose = FALSE)
## plain Newton on both sides: by default the compiled backend reuses its
## factorised jacobian and takes more iterations (test_reuse.R)
res_s <- thor_solve(s, from, to, data_opale, "date", verbose = FALSE, diagnostics = TRUE,
                    reuse_jacobian = FALSE)
res_r <- thor_solve(r, from, to, data_opale, "date", verbose = FALSE, diagnostics = TRUE,
                    reuse_jacobian = FALSE)

ok("same Newton iterations as sparse",
   identical(attr(res_s, "iterations"), attr(res_r, "iterations")),
   sprintf("%d vs %d", sum(attr(res_r, "iterations")), sum(attr(res_s, "iterations"))))
check_agree(endo_matrix(res_s, s, rows), endo_matrix(res_r, r, rows),
            label = "sparse and sparse-r agree")
check_residuals(r, res_r, dates[rows], "date")

## ---- 2. the old name -----------------------------------------------------------
a <- suppressMessages(thor_model("rb_alias", "inst/models/opale.txt", backend = "dense-r",
                                 workdir = file.path(work, "a"), verbose = FALSE))
ok("backend = \"dense-r\" builds sparse-r", identical(a@backend, "sparse-r"))

## ---- 3. a model saved with the old dense R backend ----------------------------
old <- suppressMessages(thor_load("tests/old_dense_r_opale.rds",
                                  workdir = file.path(work, "old")))
res_old <- thor_solve(old, from, to, data_opale, "date", verbose = FALSE, diagnostics = TRUE)
ok("old dense-r model loads as sparse-r",
   identical(old@backend, "sparse-r") && inherits(old@jacobian$heart, "thor_sparse_jacobian"))
## Not bit-for-bit: the saved model carries the derivatives as Deriv wrote
## them, a fresh build has them from stats::D. Same formulas, different
## spelling, so the last digit can differ.
ok("and takes the same Newton iterations as a fresh sparse-r build",
   identical(attr(res_old, "iterations"), attr(res_r, "iterations")))
check_agree(endo_matrix(res_old, old, rows), endo_matrix(res_r, r, rows),
            label = "and gives the same solution")

## ---- 4. a model saved under one sort order, loaded under another ---------------
## "Alphabetical" depends on the session's collation: the underscore sorts
## before the digits under C.UTF-8 and after the capitals under C. A model's
## variable order is fixed when it is built, so loading one on a machine that
## sorts differently must neither be refused (it was: the validity check
## compared the stored order with sort()) nor renumber the columns (the code
## generators sorted the names again, which the dense-r conversion would have
## hit).
cat("\n=== another collation ===\n")
rds <- file.path(work, "saved_r.rds")
thor_save(r, rds)

before <- Sys.getlocale("LC_COLLATE")
other  <- if (identical(sort(c("a_b", "a1")), c("a1", "a_b"))) "en_US.UTF-8" else "C"
changed <- tryCatch({ Sys.setlocale("LC_COLLATE", other); TRUE }, warning = function(w) FALSE)
ok("this test can switch to a collation that sorts differently",
   changed && !identical(r@vars$all, sort(r@vars$all)),
   sprintf("%s -> %s", before, other))

there <- thor_load(rds, workdir = file.path(work, "there"))
res_there <- thor_solve(there, from, to, data_opale, "date", verbose = FALSE)
ok("a saved model loads there and solves identically",
   identical(there@vars$all, r@vars$all) &&
     identical(endo_matrix(res_there, there, rows), endo_matrix(res_r, r, rows)))

old_there <- suppressMessages(thor_load("tests/old_dense_r_opale.rds",
                                        workdir = file.path(work, "old_there")))
res_old_there <- thor_solve(old_there, from, to, data_opale, "date", verbose = FALSE)
check_agree(endo_matrix(res_old_there, old_there, rows), endo_matrix(res_r, r, rows),
            label = "an old dense-r model is converted correctly there")
invisible(Sys.setlocale("LC_COLLATE", before))

cat("\nPASS\n")
