## Driver for the two compiled backends.
##
## The R side does nothing per-iteration: it hands the prepared matrix to the
## generated `thor_cpp_solve()` and takes the result back. The sparse and the
## dense-cpp backends generate the same entry point, so there is nothing to
## branch on here.

#' Run a compiled model's solver
#'
#' @param env environment returned by [model_env()]
#' @param M numeric data matrix, variables in alphabetical order
#' @param first,last 1-based rows of the first and last period to solve
#' @param rtol,atol,max_iter,damping,verbose see [thor_solve()]
#' @return list with `data`, `iterations`, `residuals`, `convergence`
#' @keywords internal
solve_cpp <- function(env, M, first, last, rtol, atol, max_iter, damping, verbose) {
  env$thor_cpp_solve(M,
                     as.integer(first - 1L),   # 0-based rows for C++
                     as.integer(last - 1L),
                     rtol, atol, as.integer(max_iter),
                     isTRUE(damping), isTRUE(verbose))
}

#' Residuals of a compiled model at one row
#'
#' @param env environment returned by [model_env()]
#' @param M numeric data matrix
#' @param row 1-based row
#' @return named numeric vector, one entry per block
#' @keywords internal
residuals_cpp <- function(env, M, row) {
  env$thor_cpp_residuals(M, as.integer(row - 1L))
}
