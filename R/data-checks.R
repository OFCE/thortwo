## Standalone data checks.
##
## `thor_solve()` runs the checks it needs itself, on the matrix it has just
## built, and reports what is missing. These are for looking at a database
## before committing to a solve.

#' Are all the model's variables in the database
#'
#' @param model a `thor_model`
#' @param database the database to check
#' @param quiet logical. FALSE to report success as well as failure.
#' @return TRUE if every model variable is present
#' @export
data_model_checks <- function(model, database, quiet = TRUE) {

  if (!methods::is(model, "thor_model")) {
    stop("`model` must be a thor_model.", call. = FALSE)
  }
  have <- names(database)

  report <- function(vars, label) {
    absent <- setdiff(vars, have)
    if (length(absent)) {
      cat("The following ", label, " are missing from the database:\n", sep = "")
      print(absent)
    }
    length(absent) == 0L
  }

  ok <- all(c(report(model@vars$endo,  "endogenous variables"),
              report(model@vars$exo,   "exogenous variables"),
              report(model@vars$coeff, "coefficients")))

  if (ok && !quiet) cat("All model variables are in the database.\n")
  ok
}

#' Is this column usable as a time index
#'
#' It must have no duplicates, no NAs, and be sorted in ascending order.
#'
#' @param database the database to check
#' @param index_time name of the column
#' @return TRUE if the column can be used as a time index
#' @export
time_model_checks <- function(database, index_time) {

  assertthat::assert_that(assertthat::is.string(index_time))
  if (!index_time %in% names(database)) {
    stop("'", index_time, "' is not a column of the database.", call. = FALSE)
  }

  time_vec <- database[[index_time]]
  check <- TRUE

  if (anyDuplicated(time_vec)) {
    cat("The chosen index_time has duplicated values. It is not a viable time index.\n")
    check <- FALSE
  }
  if (anyNA(time_vec)) {
    cat("The chosen index_time has missing values. It is not a viable time index.\n")
    check <- FALSE
  }
  if (is.unsorted(order(time_vec))) {
    cat("The chosen index_time is not sorted in ascending order. It is not a viable time index.\n")
    check <- FALSE
  }

  check
}

#' Which variables have missing values over given periods
#'
#' @param database the database to check
#' @param times periods to test, as found in `index_time`
#' @param variables variables to check
#' @param index_time name of the time column
#' @return character vector of the variables with an NA over `times`
#' @export
na_report_variables_times <- function(database, times, variables, index_time = "date") {

  if (!time_model_checks(database, index_time)) {
    stop("Problems with the index_time variable specified.", call. = FALSE)
  }
  key <- as.character(database[[index_time]])
  rows <- match(as.character(times), key)
  if (anyNA(rows)) {
    stop("Periods not found in '", index_time, "': ",
         paste(as.character(times)[is.na(rows)], collapse = ", "), call. = FALSE)
  }

  missing <- setdiff(variables, names(database))
  if (length(missing)) {
    stop("Not in the database: ", paste(utils::head(missing, 10), collapse = ", "),
         call. = FALSE)
  }

  sub <- database[rows, variables, drop = FALSE]
  names(sub)[vapply(sub, anyNA, logical(1))]
}

#' Does the solver have enough data to run over these periods
#'
#' Checks that the observation before the first period is complete, and that
#' no exogenous variable or coefficient is missing over the span. Those are the
#' two conditions that make a solve fail for reasons that have nothing to do
#' with the model.
#'
#' @param model a `thor_model`
#' @param database the database to check
#' @param index_time name of the time column
#' @param times periods to test
#' @return named logical vector, one entry per period
#' @export
time_solver_test_run <- function(model, database, index_time = "date", times) {

  if (!methods::is(model, "thor_model")) {
    stop("`model` must be a thor_model.", call. = FALSE)
  }
  if (!time_model_checks(database, index_time)) {
    stop("Problems with the index_time variable specified.", call. = FALSE)
  }
  if (!data_model_checks(model, database)) {
    stop("Missing model variables in the database.", call. = FALSE)
  }

  times <- sort(times)
  key <- as.character(database[[index_time]])
  rows <- match(as.character(times), key)
  if (anyNA(rows)) {
    stop("Periods not found in '", index_time, "': ",
         paste(as.character(times)[is.na(rows)], collapse = ", "), call. = FALSE)
  }

  fixed <- c(model@vars$exo, model@vars$coeff)
  stats::setNames(vapply(rows, function(r) {
    if (r == 1L) {
      cat("The first observation of the database cannot be solved: ",
          "there is no previous observation to initialise from.\n", sep = "")
      return(FALSE)
    }
    prev_ok <- !anyNA(database[r - 1L, model@vars$all])
    if (!prev_ok) {
      cat("The observation before ", key[r], " is incomplete, so the solver ",
          "cannot initialise there.\n", sep = "")
    }
    fixed_ok <- length(fixed) == 0L || !anyNA(database[r, fixed])
    if (!fixed_ok) {
      cat("Exogenous variables or coefficients are missing at ", key[r], ".\n", sep = "")
    }
    prev_ok && fixed_ok
  }, logical(1)), as.character(times))
}
