## The pure-R solver, for models built with `backend = "dense-r"`.
##
## It is a line-by-line mirror of `newton_block()` in the generated C++
## (codegen-runtime.R): same start value, same convergence measure, same
## backtracking. That is deliberate -- it is what lets the R backend be
## checked against the compiled ones to machine precision rather than merely
## approximately.
##
## It is also, unavoidably, one to two orders of magnitude slower: every
## residual and every jacobian entry is an interpreted R expression. Use it
## when there is no compiler, or on small models.

#' Scale-free convergence measure
#'
#' A step is acceptable when `|dx_i| <= rtol*|x_i| + atol`; this returns that
#' in normalised form, so `<= 1` means converged. See the long comment in
#' `CPP_PREAMBLE` for why it is written this way round.
#'
#' @param dx Newton step
#' @param x current iterate
#' @param rtol,atol tolerances
#' @return numeric scalar
#' @keywords internal
step_measure <- function(dx, x, rtol, atol) {
  max(abs(dx) / (rtol * abs(x) + atol))
}

#' Newton on one block at one period
#'
#' @param B a generated block: `list(n, endo, res, jac, seed)`
#' @param J the block's jacobian, pre-seeded with its constant entries
#' @param M numeric data matrix
#' @param t 1-based row
#' @param rtol,atol,max_iter,damping see [thor_solve()]
#' @return list with `M`, `iter` (-1 if it did not converge), `resid`, `conv`
#' @keywords internal
newton_block_r <- function(B, J, M, t, rtol, atol, max_iter, damping) {

  cols <- B$endo
  x <- M[t, cols]
  f <- B$res(t, M)
  rnorm <- max(abs(f))
  conv <- Inf

  for (it in seq_len(max_iter)) {

    if (is.na(rnorm)) return(list(M = M, iter = -1L, resid = rnorm, conv = conv))

    J <- B$jac(t, M, J)
    dx <- tryCatch(solve(J, f), error = function(e) NULL)
    if (is.null(dx) || !all(is.finite(dx))) {
      return(list(M = M, iter = -1L, resid = rnorm, conv = conv))
    }

    ## Full Newton step, backtracked only if it makes the residual worse.
    lambda <- 1
    ok <- FALSE
    xtry <- x
    new_rnorm <- rnorm
    for (b in seq_len(if (damping) 12L else 1L)) {
      xtry <- x - lambda * dx
      M[t, cols] <- xtry
      f <- B$res(t, M)
      new_rnorm <- max(abs(f))
      if (!damping || (!is.na(new_rnorm) && new_rnorm <= rnorm)) { ok <- TRUE; break }
      lambda <- lambda / 2
    }

    conv <- step_measure(lambda * dx, x, rtol, atol)

    if (!ok) {
      ## No descent direction. Near the solution this just means the residual
      ## has bottomed out in floating point, which is a success.
      if (conv <= 1) return(list(M = M, iter = it, resid = rnorm, conv = conv))
      M[t, cols] <- x
      return(list(M = M, iter = -1L, resid = rnorm, conv = conv))
    }

    x <- xtry
    rnorm <- new_rnorm

    if (conv <= 1) return(list(M = M, iter = it, resid = rnorm, conv = conv))
  }

  list(M = M, iter = -1L, resid = rnorm, conv = conv)
}

#' Run the R solver over a span of periods
#'
#' @param env environment returned by [model_env()]
#' @param M numeric data matrix, variables in alphabetical order
#' @param first,last 1-based rows of the first and last period to solve
#' @param rtol,atol,max_iter,damping,verbose see [thor_solve()]
#' @return list with `data`, `iterations`, `residuals`, `convergence`
#' @keywords internal
solve_r <- function(env, M, first, last, rtol, atol, max_iter, damping, verbose) {

  blocks <- env$blocks
  ## The constant entries of each jacobian are planted once, exactly as the
  ## compiled backends do in Block::setup().
  seeds <- lapply(blocks, function(B) B$seed())

  n_per <- last - first + 1L
  iters <- integer(n_per); resid <- numeric(n_per); conv <- numeric(n_per)

  for (t in seq.int(first, last)) {
    ## carry the previous period forward as the starting point
    for (B in blocks) M[t, B$endo] <- M[t - 1L, B$endo]

    tot <- 0L; worst <- 0; worstc <- 0
    for (k in seq_along(blocks)) {
      B <- blocks[[k]]
      r <- newton_block_r(B, seeds[[k]], M, t, rtol, atol, max_iter, damping)
      M <- r$M
      if (r$iter < 0L) {
        stop(sprintf(paste0("Newton did not converge on block '%s' at row %d after %d ",
                            "iterations: scaled step %.3e (converges at 1.0, rtol %.3e), ",
                            "max |residual| %.3e. Try a looser rtol, a higher max_iter, ",
                            "or check the data at that period."),
                     B$name, t, max_iter, r$conv, rtol, r$resid), call. = FALSE)
      }
      tot <- tot + r$iter
      worst <- max(worst, r$resid); worstc <- max(worstc, r$conv)
    }

    i <- t - first + 1L
    iters[i] <- tot; resid[i] <- worst; conv[i] <- worstc
    if (verbose) cat("  ", t, " (", tot, " it) ", sep = "")
  }
  if (verbose) cat("\n")

  list(data = M, iterations = iters, residuals = resid, convergence = conv)
}

#' Residuals of an R-backend model at one row
#'
#' @param env environment returned by [model_env()]
#' @param M numeric data matrix
#' @param row 1-based row
#' @return named numeric vector, one entry per block
#' @keywords internal
residuals_r <- function(env, M, row) {
  vapply(env$blocks, function(B) max(abs(B$res(row, M))), numeric(1))
}
