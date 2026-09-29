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
#' @param by `"block"` (the default) for the largest absolute residual per
#'   block, or `"equation"` for every equation's residual. Use `"equation"` to
#'   find out *which* equation is off -- when checking a calibration, for
#'   instance, where the answer is a list of equations to go and look at.
#' @param cache see [thor_solve()]
#' @return a matrix of residuals, periods by block (absolute values, `"block"`)
#'   or periods by equation (signed, `"equation"`, with the equation names as
#'   column names)
#' @export
model_residuals <- function(model, data, periods = NULL, index_time = "date",
                            by = c("block", "equation"), cache = NULL) {

  by <- match.arg(by)

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
  r_fun <- if (by == "equation") {
    if (model@backend == "dense-r") residuals_r_all else residuals_cpp_all
  } else {
    if (model@backend == "dense-r") residuals_r else residuals_cpp
  }

  ## The generated code returns a block's equations in the order the model
  ## object lists them, and the blocks in solve order, so the labels come
  ## straight off the object.
  labels <- if (by == "equation") {
    ids <- unlist(lapply(model@blocks, `[[`, "equations"), use.names = FALSE)
    nm <- model@equations$name[match(ids, model@equations$id)]
    ifelse(is.na(nm), ids, nm)
  } else {
    names(model@blocks)
  }
  ncols <- length(labels)

  ## vapply returns one column per row of data; the result is periods by
  ## columns. Built explicitly rather than with t(), which collapses to the
  ## wrong shape when there is only one column.
  vals <- vapply(rows, function(r) r_fun(env, M, r), numeric(ncols))
  matrix(vals, nrow = length(rows), ncol = ncols, byrow = TRUE,
         dimnames = list(key[rows], labels))
}

#' Which equations are not satisfied at a given period
#'
#' A calibration check: evaluate every equation at one period and return the
#' ones whose residual exceeds `tolerance`, worst first. An empty result means
#' the data satisfies the model there.
#'
#' @param model a `thor_model`
#' @param data data.frame holding the data
#' @param period the period to check, as found in `index_time`
#' @param index_time name of the time column in `data`
#' @param tolerance the largest residual to accept
#' @param cache see [thor_solve()]
#' @return a data.frame with `equation`, `part`, `formula` and `residual`, with
#'   no rows when every equation is satisfied
#' @export
calibration_check <- function(model, data, period, index_time = "date",
                              tolerance = 1e-6, cache = NULL) {

  if (length(period) != 1L) {
    stop("`period` must be a single period.", call. = FALSE)
  }
  r <- model_residuals(model, data, periods = period, index_time = index_time,
                       by = "equation", cache = cache)
  res <- r[1L, ]

  off <- which(!is.finite(res) | abs(res) >= tolerance)
  eq <- model@equations
  pos <- match(names(res)[off], eq$name)

  out <- data.frame(equation = names(res)[off],
                    part     = eq$part[pos],
                    formula  = eq$equation[pos],
                    residual = unname(res[off]),
                    stringsAsFactors = FALSE)
  out[order(abs(out$residual), decreasing = TRUE), , drop = FALSE]
}
