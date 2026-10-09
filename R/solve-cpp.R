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
#' @param labels the periods as the time index names them, one per row of `M`;
#'   used in the progress line and in error messages
#' @return list with `data`, `iterations`, `residuals`, `convergence`
#' @keywords internal
solve_cpp <- function(env, M, first, last, rtol, atol, max_iter, damping, verbose,
                      reuse = FALSE, reuse_asked = FALSE,
                      labels = as.character(seq_len(nrow(M)))) {
  ## A model saved before 2026-10 carries generated code whose entry point has
  ## no `reuse` argument. It still solves; it just cannot reuse.
  if (!"reuse" %in% names(formals(env$thor_cpp_solve))) {
    if (isTRUE(reuse) && isTRUE(reuse_asked)) {    # by default: just solve as it can
      warning("This model was built by an older thortwo, so `reuse_jacobian` is ",
              "ignored. Rebuild it with thor_model() to use it.", call. = FALSE)
    }
    return(env$thor_cpp_solve(M, as.integer(first - 1L), as.integer(last - 1L),
                              rtol, atol, as.integer(max_iter),
                              isTRUE(damping), isTRUE(verbose)))
  }
  ## Likewise one saved before periods were reported by name: it reports rows.
  if (!"labels" %in% names(formals(env$thor_cpp_solve))) {
    return(env$thor_cpp_solve(M, as.integer(first - 1L), as.integer(last - 1L),
                              rtol, atol, as.integer(max_iter),
                              isTRUE(damping), isTRUE(reuse), isTRUE(verbose)))
  }
  env$thor_cpp_solve(M,
                     as.integer(first - 1L),   # 0-based rows for C++
                     as.integer(last - 1L),
                     rtol, atol, as.integer(max_iter),
                     isTRUE(damping), isTRUE(reuse), isTRUE(verbose),
                     as.character(labels))
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

#' Every equation's residual of a compiled model at one row
#'
#' @param env environment returned by [model_env()]
#' @param M numeric data matrix
#' @param row 1-based row
#' @return numeric vector, one entry per equation, blocks in solve order
#' @keywords internal
residuals_cpp_all <- function(env, M, row) {
  env$thor_cpp_residuals_all(M, as.integer(row - 1L))
}
