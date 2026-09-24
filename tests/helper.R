## Shared helpers for the acceptance scripts.
##
## Usage:  Rscript tests/test_opale.R
##         Rscript tests/test_threeme.R [4x4|8x8]
##         Rscript tests/test_roundtrip.R
##
## Each script loads the package from source, builds the model on every
## backend it is meant to support, and checks the four invariants of §7 of the
## build spec.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

ok <- function(label, pass, detail = "") {
  cat(if (pass) "  ok   " else "  FAIL ", label,
      if (nzchar(detail)) paste0("  [", detail, "]") else "", "\n", sep = "")
  if (!pass) stop("Failed: ", label, call. = FALSE)
  invisible(TRUE)
}

## --- invariant 1: the Newton step converged -------------------------------
check_converged <- function(res, label = "converged") {
  worst <- max(attr(res, "convergence"))
  ok(label, worst <= 1, sprintf("worst scaled step %.3g (converges at 1)", worst))
}

## --- invariant 2: the equations are actually satisfied at the solution ----
## This is the one that catches code-generation bugs. Newton converges on
## whatever equations it was given; only the residuals notice if those are not
## the equations in the model file.
check_residuals <- function(model, res, periods, index_time, tol = 1e-6,
                            label = "residuals") {
  r <- model_residuals(model, res, periods = periods, index_time = index_time)
  worst <- max(r)
  ok(label, worst < tol, sprintf("max |residual| %.3g over %d periods x %d blocks",
                                 worst, nrow(r), ncol(r)))
  invisible(r)
}

## --- invariant 3: two backends agree on the same model -------------------
## Medians, not maxima: a few variables of every real model are genuinely
## ill-conditioned (in ThreeME, `pds_*`, a price over a near-zero stock
## change), and a maximum reports those rather than the solver.
rel_diff <- function(A, B) {
  sc <- pmax(abs(A), abs(B))
  list(rel = abs(A - B) / sc, sc = sc)
}

check_agree <- function(A, B, tol = 1e-10, label = "backends agree") {
  d <- rel_diff(A, B)
  med <- stats::median(d$rel[d$sc > 0], na.rm = TRUE)
  ok(label, med < tol, sprintf("median relative difference %.3g", med))
  invisible(med)
}

## --- invariant 4: solving over history reproduces history ----------------
## With the historical exogenous path as input, the historical endogenous path
## is the solution. This is the real test that the generated code means what
## the model file says, independently of whether Newton converged.
check_history <- function(A, B, label = "reproduces history",
                          med_tol = 1e-10, max_tol = 1e-8) {
  d <- rel_diff(A, B)
  med <- stats::median(d$rel[d$sc > 0], na.rm = TRUE)
  ## Compared only where the value is of meaningful magnitude. Opale stores a
  ## handful of contribution series (contpib*) with an exact 0.0 in one
  ## quarter, against which any non-zero solution scores a relative difference
  ## of 1; that is an artefact of the stored data, not a solver error.
  worst <- max(d$rel[d$sc > 1], na.rm = TRUE)
  ok(label, med < med_tol && worst < max_tol,
     sprintf("median %.3g, max %.3g for |value| > 1", med, worst))

  big <- d$rel; big[d$sc <= 1e-3] <- NA
  off <- sort(apply(big, 2L, max, na.rm = TRUE), decreasing = TRUE)
  cat("       worst variables (|value| > 1e-3): ",
      paste(sprintf("%s %.1g", names(off)[1:3], off[1:3]), collapse = ", "), "\n", sep = "")
  invisible(med)
}

## A scratch directory per run, cleaned up at the end.
new_tmpdir <- function(prefix) {
  d <- tempfile(prefix); dir.create(d, recursive = TRUE); d
}

endo_matrix <- function(data, model, rows) {
  as.matrix(data[rows, model@vars$endo, drop = FALSE])
}
