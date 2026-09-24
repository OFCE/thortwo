## Turning raw equation strings into the equations table, and working out
## which endogenous variables occur contemporaneously in each of them.
##
## tresthor used `splitstackshape::cSplit()` for both jobs. It is a heavy
## dependency for what is a `strsplit()` followed by padding, so the splitting
## is done in base R here. The regular expressions that do the actual parsing
## are unchanged.

#' Generate equation ids
#'
#' @param equation_list character vector of equations
#' @return character vector of `eq_###` ids, zero-padded to a fixed width
#' @keywords internal
generate_equations_index <- function(equation_list) {
  n <- length(equation_list)
  paste("eq", formatC(seq_len(n), width = nchar(as.character(n)),
                      format = "d", flag = "0"), sep = "_")
}

#' Split each string on the first occurrence of a separator
#'
#' @param x character vector
#' @param sep single-character separator
#' @return character matrix with two columns; the second is NA when `x` has no
#'   separator
#' @keywords internal
split_first <- function(x, sep) {
  at <- regexpr(sep, x, fixed = TRUE)
  left  <- ifelse(at > 0L, substr(x, 1L, at - 1L), x)
  right <- ifelse(at > 0L, substring(x, at + 1L), NA_character_)
  cbind(left, right, deparse.level = 0L)
}

#' Build the equations table
#'
#' @param eq_list character vector of equations, as written in the model file
#' @param verbose report names that had to be altered
#' @return data.frame with `id`, `name`, `equation`, `LHS`, `RHS`, `formula`
#' @keywords internal
create_equations_list <- function(eq_list, verbose = TRUE) {

  eqlist <- tolower(gsub("\\s+|\\\n", "", eq_list))
  id <- generate_equations_index(eqlist)

  ## "name : lhs = rhs", where the name is optional
  named <- split_first(eqlist, ":")
  has_name <- !is.na(named[, 2L])
  name     <- ifelse(has_name, named[, 1L], id)
  equation <- ifelse(has_name, named[, 2L], named[, 1L])

  if (anyDuplicated(name)) {
    if (verbose) {
      cat("Some equation names are duplicated. They will be made unique. ",
          "To avoid this, change the names in the model input.\n", sep = "")
    }
    name <- make.unique(name, sep = "_")
  }

  ## A user-given name that looks like an id but sits on a different equation
  ## would make `equations[id, ]` and `equations[name, ]` disagree.
  clash <- (name %in% id) & (name != id)
  if (any(clash)) {
    if (verbose) {
      cat("The following equation names clash with generated ids and will revert to the id: \n")
      print(name[clash])
    }
    name[clash] <- id[clash]
  }
  check_var_vector(name, "equation names")

  sides <- split_first(equation, "=")
  out <- data.frame(id = id, name = name, equation = equation,
                    LHS = sides[, 1L], RHS = sides[, 2L],
                    stringsAsFactors = FALSE)
  ## The residual form: every equation becomes `LHS - (RHS)`, whose root is
  ## the solution. Everything downstream differentiates and evaluates this.
  out$formula <- paste0(out$LHS, "-(", out$RHS, ")")
  rownames(out) <- out$id
  out
}

#' Which endogenous variables occur contemporaneously in each equation
#'
#' Lagged occurrences do not count: a variable that only appears under `lag()`
#' is predetermined at solve time, so it belongs to neither the block
#' decomposition nor the jacobian.
#'
#' @param formula_list character vector of residual formulas
#' @param endogenous,exogenous,coefflist declared variables
#' @param equations_index equation ids, used as row names
#' @param functionsthor functions the parser should treat as calls
#' @return data.frame, one row per equation, holding that equation's
#'   contemporaneous endogenous variables padded with NA
#' @keywords internal
table_contemporaneous_endos <- function(formula_list, endogenous, exogenous,
                                        coefflist, equations_index,
                                        functionsthor = thor_functions_supported) {

  fn <- c(functionsthor, toupper(functionsthor), "delta", "newdiff")
  function_pattern <- paste0("-|\\*|/|\\^|\\(|\\)|,|",
                             paste(paste0(fn, "\\("), collapse = "|"))

  s <- formula_list
  ## Blank out whole lagged terms, including any variable used as the lag
  ## amount: `lag(x, trim + 1)` contributes nothing at t = 0.
  s <- gsub("(mylg|lag)\\(\\w+,(-)?([0-9]|\\(\\w+\\+[0-9]\\)|\\(\\w+\\))\\)", "LAGFLAG", s)
  ## delta(n, x) -> deltan( x ) -> + x, so x is kept as contemporaneous
  s <- gsub("delta\\(([0-9]+),", "delta\\1\\(", s)
  s <- gsub("delta[0-9]+", "+", s)
  s <- gsub(function_pattern, "+", s)
  s <- gsub("[0-9]+\\.[0-9]+e", "+", s)   # scientific-notation exponents

  endo_set <- endogenous
  per_equation <- lapply(strsplit(s, "+", fixed = TRUE), function(tok) {
    tok <- tok[nzchar(tok)]
    tok <- tok[!grepl("^[0-9]+\\.?[0-9]*$", tok)]   # numeric literals
    unique(tok[tok %in% endo_set])
  })

  width <- max(1L, max(lengths(per_equation)))
  m <- vapply(per_equation, function(v) c(v, rep(NA_character_, width - length(v))),
              character(width))
  ## vapply gives one column per equation; the table is one row per equation.
  m <- if (width == 1L) matrix(m, ncol = 1L) else t(m)

  out <- as.data.frame(m, stringsAsFactors = FALSE)
  names(out) <- paste0("v", seq_len(ncol(out)))
  rownames(out) <- equations_index
  out
}

#' Rewrite formulas into a form `Deriv` and the code generators can read
#'
#' `lag(x, 2)` becomes the *symbol* `lag.x.2`, so that `Deriv` treats a lagged
#' term as an atom whose derivative with respect to the contemporaneous
#' variable is zero. The code generators decode the symbol again.
#'
#' @param formula_list character vector of residual formulas
#' @return character vector of rewritten formulas
#' @keywords internal
formatting_formulas <- function(formula_list) {
  f <- formula_list

  ## a zero lag is no lag
  f <- gsub("(mylg|lag)\\((\\w+),-?0\\)|(mylg|lag)\\((LOG\\(\\w+\\)),-?0\\)", "\\2", f)
  f <- gsub("(mylg|lag)\\((LOG\\(\\w+\\)),-?0\\)", "\\2", f)
  f <- gsub("(mylg|lag)\\((log\\(\\w+\\)),-?0\\)", "\\2", f)

  f <- gsub("(mylg|lag)\\((\\w+),-?(\\w+)\\)", "lag\\.\\2\\.\\3", f)
  f <- gsub("(mylg|lag)\\((\\w+),-?\\((\\w+)\\)\\)", "lag\\.\\2\\.\\3", f)
  f <- gsub("(mylg|lag)\\((\\w+),-?\\((\\w+)\\+(\\w+)\\)\\)", "lag\\.\\2\\.\\3\\.\\4", f)
  f <- gsub("(mylg|lag)\\((LOG)\\((\\w+)\\),-?(\\w+)\\)", "lag\\.\\2\\.\\3\\.\\4", f)
  f <- gsub("(mylg|lag)\\((LOG)\\((\\w+)\\),-?\\((\\w+)\\)\\)", "lag\\.\\2\\.\\3\\.\\4", f)
  f <- gsub("(mylg|lag)\\((LOG)\\((\\w+)\\),-?\\((\\w+)\\+(\\w+)\\)\\)", "lag\\.\\2\\.\\3\\.\\4.\\4", f)

  f <- gsub("LOG", "log", f)
  f <- gsub("lag\\.log\\.(\\w+)\\.(\\w+)\\.(\\w+)", "log\\(lag\\.\\1\\.\\2\\.\\3\\)", f)
  f <- gsub("lag\\.log\\.(\\w+)\\.(\\w+)", "log\\(lag\\.\\1\\.\\2\\)", f)

  f
}

#' Every variable name occurring in a formula
#'
#' @param string a formula, as a single character string
#' @return character vector of the names found
#' @export
get_variables_from_string <- function(string) {
  string <- gsub("\\s+", "", string)
  fn <- c(thor_functions_supported, toupper(thor_functions_supported),
          "delta", "newdiff", "lag", "mylg")
  function_pattern <- paste(paste0(fn, "\\("), collapse = "|")
  symbol_pattern   <- "(\\+|\\*|,|-|/|\\^|\\(|\\)|=)|\\\\"

  x <- gsub(function_pattern, "@", string)
  x <- gsub(symbol_pattern, "@", x)
  x <- gsub("@+", ",", x)
  x <- gsub(",[0-9]+(\\.[0-9]+)?", ",", x)
  x <- gsub("^[0-9]+(\\.[0-9]+)?", "", x)
  x <- gsub("^,|,$", "", x)

  out <- unique(strsplit(x, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}
