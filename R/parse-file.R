## Reading and checking a model source file.

#' Read a model source file
#'
#' D6 -- tresthor read the model file by fixed line numbers
#' (`model_input[2]`, `[5]`, `[8]`, `[11:length]`), which is why the shipped
#' model files carry three `####DO NOT ADD BLANK LINES OR REMOVE THIS LINE####`
#' banners: one stray blank line silently shifted every section, and the model
#' was then mis-parsed rather than rejected.
#'
#' Sections are found by their header here (`endo`, `exo`, `coeff`,
#' `equations`, in that order), so blank lines and separator banners are
#' irrelevant. The legacy layout still reads correctly, since its headers are
#' the same. Anything that cannot be resolved is an error.
#'
#' @param path path to the `.txt` model file
#' @return list with `endogenous`, `exogenous`, `coefficients`, `equations`
#' @export
read_model_source <- function(path) {

  if (!file.exists(path)) stop("Model file not found: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)

  ## A header is a line that is only a section keyword, optionally followed by
  ## punctuation ("endo:", "exogenous variables :", "coefficients"). A line
  ## holding an actual variable list or equation never matches, because those
  ## contain commas, equal signs or other names.
  headers <- c(endogenous   = "^\\s*endo(genous)?( +variables?)?\\s*:?\\s*$",
               exogenous    = "^\\s*exo(genous)?( +variables?)?\\s*:?\\s*$",
               coefficients = "^\\s*coeff(icients?)?\\s*:?\\s*$",
               equations    = "^\\s*equations?\\s*:?\\s*$")

  at <- vapply(headers, function(p) {
    hit <- grep(p, lines, ignore.case = TRUE)
    if (length(hit) == 0L) NA_integer_ else hit[1L]
  }, integer(1))

  missing_sections <- names(at)[is.na(at)]
  if (length(missing_sections)) {
    stop("Could not find the ", paste(missing_sections, collapse = " and "),
         " section header in '", basename(path), "'.\n",
         "A model file must contain four section headers, in this order: ",
         "endogenous, exogenous, coefficients, equations. ",
         "See model_source_example().", call. = FALSE)
  }
  if (is.unsorted(at, strictly = TRUE)) {
    stop("The sections of '", basename(path),
         "' are out of order. Expected endogenous, exogenous, coefficients, ",
         "then equations; found them at lines ",
         paste(at[order(at)], collapse = ", "), ".", call. = FALSE)
  }

  ## Everything between one header and the next, minus blank lines and
  ## separator banners. A separator is a line opening on a run of punctuation:
  ## ` ##############`, but also `####DO NOT ADD BLANK LINES...####`, which
  ## carries a message and so is not punctuation-only. Nothing else in a model
  ## file can start that way -- a variable name, and therefore also an
  ## equation, has to start with a letter -- so the rule is safe.
  body <- function(k) {
    from <- at[k] + 1L
    to   <- if (k < length(at)) at[k + 1L] - 1L else length(lines)
    if (from > to) return(character(0))
    x <- lines[from:to]
    x <- x[!grepl("^\\s*$", x)]
    x[!grepl("^\\s*[[:punct:]]{2,}", x)]
  }

  split_names <- function(x) {
    if (length(x) == 0L) return(character(0))
    v <- trimws(unlist(strsplit(paste(x, collapse = ","), ",", fixed = TRUE)))
    tolower(unique(v[nzchar(v)]))
  }

  eqs <- body(4L)
  if (length(eqs) == 0L) {
    stop("The equations section of '", basename(path), "' is empty.", call. = FALSE)
  }

  list(endogenous   = split_names(body(1L)),
       exogenous    = split_names(body(2L)),
       coefficients = split_names(body(3L)),
       equations    = unique(eqs))
}

#' Check that a model file looks like one before parsing it
#'
#' @param path path to the file
#' @return TRUE, invisibly; stops on a malformed file
#' @keywords internal
check_model_file <- function(path) {
  parts <- read_model_source(path)
  if (length(parts$endogenous) == 0L) {
    stop("'", basename(path), "' declares no endogenous variables.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Basic equation checks
#'
#' @param equation_vector character vector of equations, as written
#' @return TRUE, invisibly; stops if any check fails
#' @keywords internal
check_equation_input <- function(equation_vector) {
  ev <- equation_vector
  ok <- TRUE
  wrongcharacters <- "@|\\$|%|#|~|\\||;"

  if (!all(grepl("=", ev, fixed = TRUE))) {
    cat("The following equations do not have the = sign: \n")
    print(ev[!grepl("=", ev, fixed = TRUE)])
    ok <- FALSE
  }
  if (any(stringr::str_count(ev, "=") > 1)) {
    cat("The following equations have too many equal signs: \n")
    print(ev[stringr::str_count(ev, "=") > 1])
    ok <- FALSE
  }
  if (any(stringr::str_count(ev, ":") > 1)) {
    cat("Colons are used to give names to equations. The following equations have too many : signs: \n")
    print(ev[stringr::str_count(ev, ":") > 1])
    ok <- FALSE
  }
  if (any(grepl(wrongcharacters, ev))) {
    cat("The following equations contain some invalid characters: \n")
    print(ev[grepl(wrongcharacters, ev)])
    ok <- FALSE
  }

  if (!ok) stop("Please correct the model's equations.", call. = FALSE)
  invisible(TRUE)
}

#' Check a vector of variable names
#'
#' @param variable_vector character vector of names
#' @param name what is being checked, for the message
#' @return TRUE, invisibly; stops if any check fails
#' @keywords internal
check_var_vector <- function(variable_vector, name = "variables") {
  vv <- variable_vector
  if (length(vv) == 0L) return(invisible(TRUE))
  ok <- TRUE

  if (!all(grepl("^[a-z]", vv))) {
    cat("Something is amiss with ", name, ". Variables must start with a lower-case letter. \n", sep = "")
    print(vv[!grepl("^[a-z]", vv)])
    ok <- FALSE
  }
  if (any(grepl("[^a-zA-Z0-9_]", vv))) {
    cat("Something is amiss with ", name,
        ". Only alphanumeric characters and underscores are allowed in variable names. \n", sep = "")
    print(vv[grepl("[^a-zA-Z0-9_]", vv)])
    ok <- FALSE
  }

  if (!ok) stop("Please correct the ", name, " vector.", call. = FALSE)
  invisible(TRUE)
}

#' Check that two sets of variables do not overlap
#'
#' @param var1,var2 character vectors
#' @param what1,what2 what each vector is, for the message
#' @return TRUE, invisibly; stops on an overlap
#' @keywords internal
check_variable_conflict <- function(var1, var2, what1 = "the first set",
                                    what2 = "the second set") {
  both <- intersect(var1, var2)
  if (length(both)) {
    cat("The following variables appear in both ", what1, " and ", what2, ":\n", sep = "")
    print(both)
    stop("A variable cannot be declared twice.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that the model is square and that the parser found what was declared
#'
#' @param endo declared endogenous variables
#' @param eqns contemporaneous-endogenous table, one row per equation
#' @param names_of_equations equation names, in the order of `eqns`
#' @return TRUE, invisibly; stops if any check fails
#' @keywords internal
check_eq_var_identification <- function(endo, eqns, names_of_equations) {
  ok <- TRUE
  found <- stats::na.omit(unique(as.vector(as.matrix(eqns))))

  undetected <- setdiff(endo, found)
  if (length(undetected)) {
    cat("The following declared endogenous variables were not detected in the equations at t = 0. \n")
    print(undetected)
    ok <- FALSE
  }

  undeclared <- setdiff(found, endo)
  if (length(undeclared)) {
    cat("The following names were found in the equations but are not declared endogenous. ",
        "They may be exogenous variables, or the syntax of the input may be wrong. \n", sep = "")
    print(undeclared)
    ok <- FALSE
  }

  if (length(endo) != nrow(eqns)) {
    cat("The model must be square: as many equations as endogenous variables. Found ",
        nrow(eqns), " equations for ", length(endo), " endogenous variables. \n", sep = "")
    ok <- FALSE
  }

  empty <- apply(eqns, 1L, function(x) all(is.na(x)))
  if (any(empty)) {
    cat("The following equations have no endogenous variable at t = 0: \n")
    print(names_of_equations[empty])
    ok <- FALSE
  }

  if (!ok) stop("Please check the model input and try again.", call. = FALSE)
  invisible(TRUE)
}

#' Every variable name occurring in a set of equations
#'
#' The expensive half of [is_in_formulas()]: on ThreeME 29x33 it takes 2 s.
#' [thor_model()] checks three lists of names against the same equations, so
#' it works this out once and passes it on.
#'
#' @param formulas character vector of equations
#' @return character vector of names, lower case
#' @keywords internal
variables_in_formulas <- function(formulas) {
  ## drop an equation name; `name : eq`, with spaces, is as valid as `name:eq`
  ## (with the spaces left in, "demand : y = ..." read as the token "demand:y"
  ## and `y` was dropped as unused)
  formulas <- gsub("^\\s*\\w+\\s*:", "", formulas)
  get_variables_from_string(tolower(paste(formulas, collapse = "-")))
}

#' Which of these variables actually occur in the equations
#'
#' @param variables_to_test character vector of names
#' @param formulas character vector of equations
#' @param print_type what is being tested, for the message
#' @param verbose report what was found as well as what was dropped
#' @param present the names occurring in `formulas`, as returned by
#'   [variables_in_formulas()], when the caller already has them
#' @return the subset of `variables_to_test` that occurs in `formulas`
#' @keywords internal
is_in_formulas <- function(variables_to_test, formulas, print_type = "", verbose = TRUE,
                           present = variables_in_formulas(formulas)) {
  if (length(variables_to_test) == 0L) return(character(0))
  variables_to_test <- tolower(variables_to_test)

  absent  <- setdiff(variables_to_test, present)

  if (verbose) {
    if (length(absent) == 0L) {
      cat("All ", print_type, " variables are present in the equations list. \n", sep = "")
    } else {
      cat("The following ", print_type,
          " variable(s) do not appear in any equation and will be dropped:\n", sep = "")
      print(absent)
    }
  }
  setdiff(variables_to_test, absent)
}

#' Is this an acceptable variable name
#'
#' @param text character vector to test
#' @return logical scalar
#' @keywords internal
acceptable_var_name <- function(text) {
  check <- TRUE
  if (sum(grepl("^[a-z]", text)) != length(text)) {
    cat("\nNames must start with a lower-case letter.\n")
    print(text[!grepl("^[a-z]", text)])
    check <- FALSE
  }
  if (any(nchar(text) == 0L)) {
    cat("\nNames must have at least one character.\n")
    check <- FALSE
  }
  leftover <- gsub("[a-z_0-9]", "", text)
  if (any(nchar(leftover) > 0L)) {
    cat("\nInvalid characters found. Only lower-case alphanumerics and underscores are allowed.\n")
    print(text[nchar(leftover) > 0L])
    check <- FALSE
  }
  check
}

#' Can the lags and deltas be parsed safely
#'
#' @param formula_list character vector of equations
#' @return logical scalar
#' @keywords internal
parser_lag_delta_check <- function(formula_list) {
  check <- TRUE
  all_in_one <- tolower(paste(formula_list, collapse = "-"))

  n_lags       <- stringr::str_count(all_in_one, "(mylg|lag)\\(")
  correct_lags <- stringr::str_count(all_in_one, "(mylg|lag)\\([a-z](\\w+)?,")
  if (correct_lags != n_lags) {
    cat("\nSome lagged variables are not written in a form the parser can read.",
        "\nOnly one variable may be lagged at a time:",
        "\n  ok      : lag(my_variable42, 1)",
        "\n  invalid : lag(my_variable42 + var2, 1)\n")
    check <- FALSE
  }

  n_deltas       <- stringr::str_count(all_in_one, "delta\\(")
  correct_deltas <- stringr::str_count(all_in_one, "delta\\([0-9]+,")
  if (correct_deltas != n_deltas) {
    cat("\nSome deltas are not written in a form the parser can read.",
        "\nThe first argument of delta() must be a literal integer:",
        "\n  ok      : delta(1, my_variable)",
        "\n  invalid : delta(n, my_variable)\n")
    check <- FALSE
  }

  check
}
