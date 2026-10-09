## The pure-R solver, for models built with `backend = "sparse-r"`.
##
## It is a line-by-line mirror of `newton_block()` in the generated C++
## (codegen-runtime.R): same start value, same convergence measure, same
## backtracking. That is deliberate -- it is what lets the R backend be
## checked against the compiled ones to machine precision rather than merely
## approximately.
##
## The linear algebra is sparse, as in the compiled `sparse` backend: the
## jacobian is a Matrix::dgCMatrix whose pattern is fixed when the block is
## set up, and each Newton iteration only refills its values and factorises
## it with Matrix's sparse LU. Until 2026-10 this backend used a dense matrix
## and base solve(), whose cost grows with the cube of the block size; on
## ThreeME 4x4 the dense solves were nearly half of a 48 s run, and the rest of
## the gap is below.
##
## R's byte-code compiler is switched off while the generated code runs. The
## generated functions are a few thousand statements each, and compiling them
## on first call cost more than it saved: ~20 s of that 48 s run.
##
## It is still slower than the compiled backends -- every residual and
## jacobian entry is an interpreted R expression -- so it is the choice when
## there is no compiler.

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
#' @param B a generated block: `list(n, endo, res, jac, seed, i, p)`
#' @param J the block's jacobian, a `dgCMatrix` whose values are pre-seeded
#'   with its constant entries
#' @param M numeric data matrix
#' @param t 1-based row
#' @param rtol,atol,max_iter,damping see [thor_solve()]
#' @param reuse logical. Keep the factorised jacobian between iterations; see
#'   `newton_block` in the generated C++, which this mirrors line by line.
#' @return list with `M`, `iter` (-1 if it did not converge), `resid`, `conv`
#' @keywords internal
newton_block_r <- function(B, J, M, t, rtol, atol, max_iter, damping, reuse = FALSE) {

  cols <- B$endo
  x <- M[t, cols]
  f <- B$res(t, M)
  rnorm <- max(abs(f))
  conv <- Inf
  fail <- function() list(M = M, iter = -1L, resid = rnorm, conv = conv)

  lu <- NULL
  stale <- FALSE      # factors from an earlier iteration are worth trying

  for (it in seq_len(max_iter)) {

    if (is.na(rnorm)) return(fail())

    fresh <- !stale
    repeat {                      # old factors first, then fresh ones
      if (fresh) {
        J@x <- B$jac(t, M, J@x)
        ## Matrix caches a matrix's factorisation inside the object, and lu()
        ## returns the cached one. The values were just replaced behind its
        ## back, so the cache must go, or this would silently be the first
        ## iteration's factorisation every time.
        J@factors <- list()
        lu <- tryCatch(suppressWarnings(Matrix::lu(J)), error = function(e) NULL)
        if (is.null(lu)) return(fail())
      }
      ## `tol = 0`: Matrix refuses to solve when the smallest and the largest
      ## pivot are more than 1/eps apart, calling the matrix "computationally
      ## singular". That is a test of scaling, not of singularity, and a
      ## model's jacobian can fail it while being perfectly solvable: the
      ## epilogue of ThreeME 4x4 has pivots from 9e-13 to 1e7, each in an
      ## equation of its own, and the step it gives is exact to 1e-16. The
      ## compiled backends make no such test. A matrix that really is singular
      ## is still caught, by lu() or by the check on `dx` just below.
      dx <- tryCatch(as.numeric(suppressWarnings(Matrix::solve(lu, f, tol = 0))),
                     error = function(e) NULL)
      if (is.null(dx) || !all(is.finite(dx))) {
        if (!fresh) { fresh <- TRUE; next }
        return(fail())
      }

      ## Full Newton step, backtracked only if it makes the residual worse.
      f0 <- f
      lambda <- 1
      ok <- FALSE
      xtry <- x
      new_rnorm <- rnorm
      for (b in seq_len(if (fresh && damping) 12L else 1L)) {
        xtry <- x - lambda * dx
        M[t, cols] <- xtry
        f <- B$res(t, M)
        new_rnorm <- max(abs(f))
        if ((fresh && !damping) || (!is.na(new_rnorm) && new_rnorm <= rnorm)) { ok <- TRUE; break }
        lambda <- lambda / 2
      }

      conv <- step_measure(lambda * dx, x, rtol, atol)

      if (ok || fresh || conv <= 1) break
      ## The old factors gave a step that does not help: take it back.
      M[t, cols] <- x
      f <- f0
      fresh <- TRUE
    }

    if (!ok) {
      ## No descent direction. Near the solution this just means the residual
      ## has bottomed out in floating point, which is a success.
      if (conv <= 1) return(list(M = M, iter = it, resid = rnorm, conv = conv))
      M[t, cols] <- x
      return(fail())
    }

    stale <- reuse && new_rnorm <= 0.5 * rnorm
    x <- xtry
    rnorm <- new_rnorm

    ## a step with old factors must be much smaller to count as the last one:
    ## see REUSE_TIGHT in the generated C++
    if (conv <= (if (fresh) 1 else 1e-3)) {
      return(list(M = M, iter = it, resid = rnorm, conv = conv))
    }
  }

  fail()
}

#' Run the R solver over a span of periods
#'
#' @param env environment returned by [model_env()]
#' @param M numeric data matrix, variables in alphabetical order
#' @param first,last 1-based rows of the first and last period to solve
#' @param rtol,atol,max_iter,damping,verbose see [thor_solve()]
#' @param reuse see `reuse_jacobian` in [thor_solve()]
#' @param labels the periods as the time index names them, one per row of `M`;
#'   used in the progress line and in error messages
#' @return list with `data`, `iterations`, `residuals`, `convergence`
#' @keywords internal
solve_r <- function(env, M, first, last, rtol, atol, max_iter, damping, verbose,
                    reuse = FALSE, labels = as.character(seq_len(nrow(M)))) {

  without_jit()

  blocks <- env$blocks
  ## The sparsity pattern is fixed and the constant entries of each jacobian
  ## are planted once, exactly as the compiled backends do in Block::setup().
  seeds <- lapply(blocks, function(B) if (is.null(B$seq)) block_matrix_r(B))

  n_per <- last - first + 1L
  iters <- integer(n_per); resid <- numeric(n_per); conv <- numeric(n_per)

  for (t in seq.int(first, last)) {
    ## carry the previous period forward as the starting point
    for (B in blocks) M[t, B$endo] <- M[t - 1L, B$endo]

    tot <- 0L; worst <- 0; worstc <- 0
    for (k in seq_along(blocks)) {
      B <- blocks[[k]]
      r <- if (is.null(B$seq)) {
        newton_block_r(B, seeds[[k]], M, t, rtol, atol, max_iter, damping, reuse)
      } else {
        sequential_block_r(B, M, t, rtol, atol, max_iter, damping, labels[t])
      }
      M <- r$M
      if (r$iter < 0L) {
        stop(sprintf(paste0("Newton did not converge on block '%s' at %s (row %d) after %d ",
                            "iterations: scaled step %.3e (converges at 1.0, rtol %.3e), ",
                            "max |residual| %.3e. Try a looser rtol, a higher max_iter, ",
                            "or check the data at that period."),
                     B$name, labels[t], t, max_iter, r$conv, rtol, r$resid), call. = FALSE)
      }
      tot <- tot + r$iter
      worst <- max(worst, r$resid); worstc <- max(worstc, r$conv)
    }

    i <- t - first + 1L
    iters[i] <- tot; resid[i] <- worst; conv[i] <- worstc
    if (verbose) cat("  ", labels[t], " (", tot, " it) ", sep = "")
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
  without_jit()
  vapply(env$blocks, function(B) max(abs(B$res(row, M))), numeric(1))
}

#' Every equation's residual of an R-backend model at one row
#'
#' @param env environment returned by [model_env()]
#' @param M numeric data matrix
#' @param row 1-based row
#' @return numeric vector, one entry per equation, blocks in solve order
#' @keywords internal
residuals_r_all <- function(env, M, row) {
  without_jit()
  unlist(lapply(env$blocks, function(B) B$res(row, M)), use.names = FALSE)
}

#' A block's jacobian as a sparse matrix, constants planted
#'
#' Built straight from the generated CSC pattern, so the order of its `x`
#' slot is the order the generated `<block>_jac()` fills.
#'
#' @param B a generated block
#' @return a `Matrix::dgCMatrix`
#' @keywords internal
block_matrix_r <- function(B) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The sparse-r backend needs the Matrix package, which comes with R; ",
         "reinstall it with install.packages(\"Matrix\").", call. = FALSE)
  }
  methods::new("dgCMatrix", i = as.integer(B$i), p = as.integer(B$p),
               x = B$seed(), Dim = c(B$n, B$n))
}

#' Switch R's byte-code compiler off until the calling function returns
#'
#' See the note at the top of this file.
#'
#' @param envir the frame whose exit restores the previous setting
#' @keywords internal
without_jit <- function(envir = parent.frame()) {
  old <- compiler::enableJIT(0L)
  do.call(on.exit, list(substitute(compiler::enableJIT(old), list(old = old)), add = TRUE),
          envir = envir)
  invisible(old)
}

#' Solve one block sequentially at one period
#'
#' The R counterpart of the sequential branch of `thor_cpp_solve`: runs the
#' generated `<block>_seq()` and reports what [newton_block_r()] reports.
#'
#' @param B a generated block with a `seq` function
#' @param M numeric data matrix
#' @param t 1-based row
#' @param rtol,atol,max_iter,damping see [thor_solve()]
#' @param label the period at row `t`, as the time index names it
#' @return list with `M`, `iter`, `resid`, `conv`
#' @keywords internal
sequential_block_r <- function(B, M, t, rtol, atol, max_iter, damping,
                               label = as.character(t)) {
  r <- B$seq(t, M, rtol, atol, max_iter, damping)
  if (r$fail > 0L) {
    stop(sprintf(paste0("Sequential solve failed on block '%s' at %s (row %d): equation '%s' could ",
                        "not be solved for '%s'. Its derivative with respect to that variable is ",
                        "zero or not finite, or it did not converge in %d iterations. Check the ",
                        "data at that period."),
                 B$name, label, t, B$seq_eq[r$fail], B$seq_var[r$fail], max_iter), call. = FALSE)
  }
  resid <- max(abs(B$res(t, r$M)))
  if (is.na(resid)) {
    stop(sprintf(paste0("Sequential solve of block '%s' at %s (row %d) produced a value that is not ",
                        "a number. Check the data at that period with calibration_check()."),
                 B$name, label, t), call. = FALSE)
  }
  list(M = r$M, iter = r$iters, resid = resid, conv = r$conv)
}
