## Operators that appear in model formulas.
##
## These are needed at *solve* time, not only at build time: `delta()` occurs
## in the formulas themselves, and a user evaluating a formula by hand, or the
## R backend falling back on an unusual construct, will call it. Easy to leave
## behind when porting, and the resulting failure ("could not find function
## delta") is obscure.

#' Differences that keep the vector's length
#'
#' @param n number of periods to difference over
#' @param x numeric vector
#' @return numeric vector the same length as `x`, with `n` leading NAs
#' @rdname newdiff_delta
#' @export
newdiff <- function(n = 1, x) {
  c(rep(NA_real_, n), diff(x, lag = n))
}

#' @rdname newdiff_delta
#' @export
delta <- function(n = 1, x) {
  c(rep(NA_real_, n), diff(x, lag = n))
}

#' Add coefficients to a database
#'
#' Estimated coefficients usually live outside the database, as a two-column
#' table. The solver needs them as constant columns alongside the data, which
#' is what this does.
#'
#' @param listcoeff data.frame with a column of names and a column of values
#' @param database the database to add them to
#' @param pos.coeff.name which column of `listcoeff` holds the name
#' @param pos.coeff.value which column holds the value
#' @param overwrite logical. TRUE (the default) to replace coefficients
#'   already present in the database.
#' @return `database` with the coefficient columns added
#' @rdname add_coeffs
#' @export
add_coeffs <- function(listcoeff, database, pos.coeff.name = 1, pos.coeff.value = 2,
                       overwrite = TRUE) {

  if (missing(database)) stop("`database` is required.", call. = FALSE)

  names_col <- tolower(as.character(listcoeff[, pos.coeff.name]))
  values    <- as.numeric(listcoeff[, pos.coeff.value])

  lco <- as.data.frame(matrix(rep(values, each = nrow(database)),
                              nrow = nrow(database),
                              dimnames = list(NULL, names_col)),
                       stringsAsFactors = FALSE)

  existing <- intersect(names_col, names(database))
  new      <- setdiff(names_col, names(database))

  if (overwrite) {
    cbind(database[, !names(database) %in% existing, drop = FALSE], lco)
  } else {
    if (length(new) == 0L) stop("No new coefficients to add.", call. = FALSE)
    cbind(database, lco[, new, drop = FALSE])
  }
}
