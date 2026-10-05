## The single solve entry point.
##
## tresthor had `thor_solver()` (308 lines, with an R solver and an
## RcppArmadillo solver interleaved and selected by an `rcpp` flag) and
## `thor_solver_sparse()`. The backend is a property of the model, not of the
## call, so it is stored on the object and dispatched on here.

#' Prepare the data matrix for a solve
#'
#' Checks the time index, resolves the period bounds and builds the numeric
#' matrix whose columns are the model variables in alphabetical order -- the
#' layout every backend indexes into.
#'
#' D7 -- tresthor converted `database`'s time column to character in place and
#' converted it back at the end. Nothing is written to the caller's data here:
#' the time column is read into a separate key vector and left alone.
#'
#' @param model a `thor_model`
#' @param from,to first and last period to solve
#' @param data data.frame holding the data
#' @param index_time name of the time column
#' @return list with `M`, `first`, `last` (1-based rows) and `key`
#' @keywords internal
thor_prepare <- function(model, from, to, data, index_time = "date") {

  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!index_time %in% names(data)) {
    stop("'", index_time, "' is not a column of the data.", call. = FALSE)
  }

  key <- as.character(data[[index_time]])
  if (anyNA(key)) stop("The time index '", index_time, "' has missing values.", call. = FALSE)
  if (anyDuplicated(key)) {
    stop("The time index '", index_time, "' has duplicated values, so a period ",
         "cannot be located unambiguously.", call. = FALSE)
  }

  first <- match(as.character(from), key)
  last  <- match(as.character(to),   key)
  if (is.na(first)) stop("The first period is not found in '", index_time, "'.", call. = FALSE)
  if (is.na(last))  stop("The last period is not found in '", index_time, "'.", call. = FALSE)
  if (last < first) stop("The last period comes before the first period.", call. = FALSE)
  ## The solver needs `max_lag` complete observations before the first solved
  ## one: one to start Newton from, and as many as the deepest lag in the
  ## model reads. Reading past the start of the data gives NA, which would
  ## otherwise surface as a non-convergence with an infinite residual.
  need <- max(1L, as.integer(model@meta$max_lag %||% 1L))
  if (first <= need) {
    stop("The first period to solve is ", key[first], ", row ", first,
         " of the data, but this model reaches ", need,
         " period(s) back, so it needs ", need,
         " complete observation(s) before it.\n",
         "Start at ", key[min(need + 1L, length(key))], " or later",
         if (isTRUE(model@meta$variable_lag))
           ", and note that the model also uses a lag given as a variable, whose depth is not known until the data is read"
         else "", ".", call. = FALSE)
  }

  mv <- model@vars$all
  missing <- setdiff(mv, names(data))
  if (length(missing)) {
    stop("Missing variables in the database: ",
         paste(utils::head(missing, 10), collapse = ", "),
         if (length(missing) > 10) sprintf(" (and %d more)", length(missing) - 10) else "",
         call. = FALSE)
  }

  M <- as.matrix(data[, mv, drop = FALSE])
  storage.mode(M) <- "double"

  ## Those observations also have to be complete: they seed the Newton start
  ## value and supply every lagged term.
  na_prev <- mv[apply(is.na(M[seq(first - need, first - 1L), , drop = FALSE]), 2L, any)]
  if (length(na_prev)) {
    stop("The observation(s) before the first period have missing values, so the ",
         "solver cannot initialise: ",
         paste(utils::head(na_prev, 10), collapse = ", "),
         if (length(na_prev) > 10) sprintf(" (and %d more)", length(na_prev) - 10) else "",
         call. = FALSE)
  }

  ## Exogenous variables and coefficients must be present over the whole span.
  fixed <- c(model@vars$exo, model@vars$coeff)
  if (length(fixed)) {
    span <- M[first:last, fixed, drop = FALSE]
    bad <- fixed[apply(is.na(span), 2L, any)]
    if (length(bad)) {
      stop("Exogenous variables or coefficients are missing over the solved period: ",
           paste(utils::head(bad, 10), collapse = ", "),
           if (length(bad) > 10) sprintf(" (and %d more)", length(bad) - 10) else "",
           call. = FALSE)
    }
  }

  list(M = M, first = first, last = last, key = key)
}

#' Solve a model
#'
#' Solves the model period by period, and within each period block by block,
#' by Newton's method. Which solver runs is decided by the model's backend.
#'
#' @param model a `thor_model`, built by [thor_model()]
#' @param from first period to solve, as found in `index_time`
#' @param to last period to solve
#' @param data data.frame holding the data. Required: unlike tresthor's
#'   `thor_solver()`, whose `database = t_data` default made the solver depend
#'   on a global variable, there is no default.
#' @param index_time name of the time column in `data`
#' @param rtol relative step tolerance. A block has converged when every
#'   endogenous variable satisfies `|dx| <= rtol*|x| + atol`. Relative rather
#'   than absolute because macro models mix variables of very different
#'   magnitude. Default 1e-10.
#' @param atol absolute step tolerance, which governs variables sitting at or
#'   near zero. Default 1e-8.
#' @param max_iter maximum Newton iterations per block and period. Default 100.
#' @param damping logical. Backtrack the Newton step when it would increase the
#'   residual. Default TRUE.
#' @param reuse_jacobian logical, or NULL (the default) to let the backend
#'   decide: TRUE on the compiled backends, FALSE on `sparse-r`.
#'
#'   TRUE keeps the factorised jacobian from one Newton iteration to the next
#'   within a period, and only recomputes it when a step taken with it stops
#'   reducing the residual fast enough. Each iteration is then much cheaper
#'   and more of them are needed. On a large model solved with a compiled
#'   backend, where the factorisation is nearly the whole cost of an
#'   iteration, the solve is about three times faster (ThreeME 29x33: 208 s to
#'   74 s). The solution is the same to within the tolerances, which are
#'   unchanged; the number of iterations reported is higher.
#'
#'   FALSE recomputes and refactorises the jacobian at every iteration: plain
#'   Newton. It is the default on `sparse-r`, where evaluating the equations
#'   in R is the main cost and more iterations make the solve slower.
#' @param cache directory in which to cache the compiled object between
#'   sessions, FALSE to disable, or NULL (the default) for [thor_cache_dir()].
#' @param verbose logical. Print progress. Default TRUE.
#' @param timings logical. TRUE (the default) to print how long the solve
#'   took when `verbose` is off; with `verbose` on, the summary line already
#'   says. `options(thortwo.timings = FALSE)` switches it off everywhere.
#' @param diagnostics logical. Attach per-period iteration counts, residuals
#'   and convergence measures to the result as attributes. Default FALSE.
#'
#' @return `data`, with the endogenous variables solved over `from`..`to`.
#' @export
thor_solve <- function(model, from, to, data,
                       index_time = "date",
                       rtol = 1e-10, atol = 1e-8, max_iter = 100L,
                       damping = TRUE, reuse_jacobian = NULL, cache = NULL,
                       verbose = TRUE, diagnostics = FALSE,
                       timings = getOption("thortwo.timings", TRUE)) {

  if (!methods::is(model, "thor_model")) {
    stop("`model` must be a thor_model, as returned by thor_model().", call. = FALSE)
  }
  stopifnot(rtol > 0, rtol < 0.01, atol >= 0, max_iter >= 1L)
  ## On by default where it pays: see the note on `reuse` in the generated
  ## newton_block. Asked for explicitly, it is honoured on any backend.
  reuse_asked <- !is.null(reuse_jacobian)
  if (reuse_asked) {
    assertthat::assert_that(is.logical(reuse_jacobian), length(reuse_jacobian) == 1L,
                            !is.na(reuse_jacobian))
  }
  reuse <- if (reuse_asked) reuse_jacobian else !is_r_backend(model@backend)

  prep <- thor_prepare(model, from, to, data, index_time)

  env <- model_env(model, cache = cache)

  t_run <- system.time({
    out <- if (is_r_backend(model@backend)) {
      solve_r(env, prep$M, prep$first, prep$last, rtol, atol,
              as.integer(max_iter), isTRUE(damping), isTRUE(verbose),
              reuse = reuse)
    } else {
      solve_cpp(env, prep$M, prep$first, prep$last, rtol, atol,
                as.integer(max_iter), isTRUE(damping), isTRUE(verbose),
                reuse = reuse, reuse_asked = reuse_asked)
    }
  })

  solved <- out$data
  colnames(solved) <- model@vars$all
  data[, model@vars$endo] <- solved[, model@vars$endo, drop = FALSE]

  if (verbose) {
    cat("Solved ", prep$last - prep$first + 1L, " periods in ",
        round(t_run[["elapsed"]], 3), " s (",
        sum(out$iterations), " Newton iterations, worst scaled step ",
        format(max(out$convergence), digits = 3), " (converges at 1), max |residual| ",
        format(max(out$residuals), digits = 3), ")\n", sep = "")
  } else if (isTRUE(timings)) {
    cat("Timings: solve ", format_seconds(t_run[["elapsed"]]), " (",
        prep$last - prep$first + 1L, " periods)\n", sep = "")
  }

  if (diagnostics) {
    span <- prep$key[prep$first:prep$last]
    attr(data, "iterations")  <- stats::setNames(out$iterations, span)
    attr(data, "residuals")   <- stats::setNames(out$residuals, span)
    attr(data, "convergence") <- stats::setNames(out$convergence, span)
    attr(data, "elapsed")     <- unname(t_run[["elapsed"]])
  }
  data
}
