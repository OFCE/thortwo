## Residual check.
##
## In tresthor this lived at the end of 1_11_compile_model.R and worked only
## on sparse models, because only they generated a `sparse_residuals()`. Every
## backend generates the same entry point here, so it works on all three --
## which is what makes it usable as the cross-backend check it was always
## meant to be.

#' Residuals of a model's equations at given periods
#'
#' Evaluates every equation of the model at the given rows of the data and
#' returns the largest absolute residual per block. This is the check that
#' catches code-generation bugs: a Newton iteration can converge on the wrong
#' equations perfectly happily, and only the residuals notice.
#'
#' @param model a `thor_model`
#' @param data data.frame holding the data, typically the result of
#'   [thor_solve()]
#' @param periods the periods to check, as found in `index_time`. Default: all
#'   rows except the first.
#' @param index_time name of the time column in `data`
#' @param cache see [thor_solve()]
#' @return a matrix of maximum absolute residuals, periods by block
#' @export
model_residuals <- function(model, data, periods = NULL, index_time = "date",
                            cache = NULL) {

  if (!methods::is(model, "thor_model")) {
    stop("`model` must be a thor_model, as returned by thor_model().", call. = FALSE)
  }
  if (!index_time %in% names(data)) {
    stop("'", index_time, "' is not a column of the data.", call. = FALSE)
  }

  key <- as.character(data[[index_time]])
  if (is.null(periods)) {
    rows <- seq_along(key)[-1L]
  } else {
    rows <- match(as.character(periods), key)
    if (anyNA(rows)) {
      stop("Periods not found in '", index_time, "': ",
           paste(as.character(periods)[is.na(rows)], collapse = ", "), call. = FALSE)
    }
  }

  mv <- model@vars$all
  missing <- setdiff(mv, names(data))
  if (length(missing)) {
    stop("Missing variables in the database: ",
         paste(utils::head(missing, 10), collapse = ", "), call. = FALSE)
  }
  M <- as.matrix(data[, mv, drop = FALSE])
  storage.mode(M) <- "double"

  env <- model_env(model, cache = cache)
  f <- if (model@backend == "dense-r") residuals_r else residuals_cpp

  nb <- length(model@blocks)
  ## vapply returns one column per row of data; the result is periods by
  ## blocks. Built explicitly rather than with t(), which collapses to the
  ## wrong shape when the model has a single block.
  vals <- vapply(rows, function(r) f(env, M, r), numeric(nb))
  out <- matrix(vals, nrow = length(rows), ncol = nb, byrow = TRUE,
                dimnames = list(key[rows], names(model@blocks)))
  out
}
