## AST-based translation of thoR formulas into C++ expressions.
##
## tresthor's legacy path wrote R source to disk, read it back as text and
## rewrote it with regular expressions. That is fragile and scales badly. Here
## we parse the formula with R's own parser and walk the resulting call tree,
## which removes the text round-trip entirely.
##
## Every backend goes through this one translator, so the three of them agree
## on what a formula means by construction rather than by testing.
##
## Conventions shared with the generated C++:
##  * `M` is the data matrix, observations in rows, variables in columns,
##    columns sorted alphabetically (same ordering as `t_data` in the solver).
##  * `vidx` maps a variable name to its 0-based column index in `M`.
##  * A time offset is carried as a list(k = integer, vterms = character):
##    the observation read is `t - k - vterms[1] - ...`.

#' The map from variable name to column of the data matrix
#'
#' Two things about it matter.
#'
#' The order is the one it is given, and it is the caller's job to give the
#' order stored on the model (`model@vars$all`), which is also the order
#' thor_solve() lays the data out in. The generators used to sort the names
#' again themselves. "Alphabetical" depends on the session's collation -- the
#' underscore sorts before the digits under C.UTF-8 and after the capitals
#' under C -- so regenerating a model's code on a machine with other settings
#' than the one that built it (which converting an old `dense-r` model does)
#' could number the columns differently from the stored order, and the solver
#' would have read the wrong ones without complaint.
#'
#' Lookups by name go through a hash table, attached as the `hash` attribute.
#' `vidx[[name]]` on a named vector is a linear scan, and the generators ask
#' once per variable reference in every equation and derivative: on ThreeME
#' 29x33 that was 42 s of a 56 s code generation. [var_col()] reads it.
#'
#' @param all_model_vars character vector of every model variable, in the
#'   model's order
#' @param base 0 for C++ column numbers, 1 for R
#' @return named integer vector with a `hash` attribute
#' @keywords internal
var_index <- function(all_model_vars, base = 0L) {
  if (anyDuplicated(all_model_vars)) {
    stop("The model's variable list has duplicated names: ",
         paste(utils::head(unique(all_model_vars[duplicated(all_model_vars)]), 5),
               collapse = ", "), call. = FALSE)
  }
  vidx <- stats::setNames(seq_along(all_model_vars) - 1L + as.integer(base), all_model_vars)
  attr(vidx, "hash") <- list2env(as.list(vidx), hash = TRUE, parent = emptyenv())
  vidx
}

#' Column of one variable
#'
#' @param vidx a map from [var_index()]; a plain named vector also works
#' @param name variable name
#' @return the column number, or NULL if the variable is unknown
#' @keywords internal
var_col <- function(vidx, name) {
  h <- attr(vidx, "hash", exact = TRUE)
  if (!is.null(h)) return(h[[name]])
  if (name %in% names(vidx)) vidx[[name]] else NULL
}

#' Zero time offset
#' @keywords internal
cpp_t0 <- function() list(k = 0L, vterms = character(0))

#' Add a lag to a time offset
#'
#' @param toff time offset, as built by `cpp_t0()`
#' @param k integer number of periods to go back
#' @param vterms character vector of C++ expressions to additionally go back by
#' @keywords internal
cpp_tshift <- function(toff, k = 0L, vterms = character(0)) {
  list(k = toff$k + as.integer(k), vterms = c(toff$vterms, vterms))
}

#' Render a time offset as a C++ row index
#'
#' Variable lag amounts are read out of the data matrix and so arrive as
#' doubles. They must be turned back into an integer explicitly: a row index
#' of type double makes `M(i, j)` resolve to `Eigen::IndexedView` rather than
#' scalar element access, which does not compile. `lround` rather than a plain
#' cast, since a whole number stored as a double may be 2.9999999.
#' @keywords internal
cpp_trender <- function(toff) {
  out <- "t"
  if (toff$k != 0L) out <- paste0(out, " - ", toff$k)
  for (v in toff$vterms) out <- paste0(out, " - (int)std::lround(", v, ")")
  out
}

#' Read a variable at a given time offset
#' @keywords internal
cpp_var <- function(name, vidx, toff) {
  idx <- var_col(vidx, name)
  if (is.null(idx) || is.na(idx)) {
    stop(sprintf("Variable '%s' is used in the model but is not a known model variable.", name))
  }
  paste0("M(", cpp_trender(toff), ",", idx, ")")
}

## Functions that map straight onto a C++ standard-library call.
cpp_unary_map <- c(
  log = "std::log", exp = "std::exp", sqrt = "std::sqrt", abs = "std::fabs",
  log2 = "std::log2", log10 = "std::log10", log1p = "std::log1p",
  expm1 = "std::expm1", sin = "std::sin", cos = "std::cos", tan = "std::tan",
  asin = "std::asin", acos = "std::acos", atan = "std::atan",
  sinh = "std::sinh", cosh = "std::cosh", tanh = "std::tanh",
  asinh = "std::asinh", acosh = "std::acosh", atanh = "std::atanh"
)

#' Decode a `lag.<var>.<n>` / `lag.<var>.<a>.<b>` symbol
#'
#' `formatting_formulas()` rewrites `lag(x, 2)` into the *symbol* `lag.x.2` so
#' that `Deriv` treats a lagged term as atomic (its derivative with respect to
#' the contemporaneous variable is zero). This undoes that encoding.
#'
#' Lag amounts given as variables are read at the current period `t`, matching
#' the behaviour of the legacy generator.
#'
#' @return NULL if `nm` is not a lag symbol, otherwise list(var, k, vterms)
#' @keywords internal
cpp_decode_lag <- function(nm, vidx) {
  if (!startsWith(nm, "lag.")) return(NULL)
  parts <- strsplit(nm, ".", fixed = TRUE)[[1]]
  if (length(parts) < 3L) {
    stop(sprintf("Malformed lag symbol '%s'.", nm))
  }
  ## A variable name never contains a dot, so parts[2] is the variable and
  ## everything after it makes up the lag amount.
  var <- parts[2L]
  amounts <- parts[-c(1L, 2L)]
  k <- 0L
  vterms <- character(0)
  for (a in amounts) {
    if (grepl("^[0-9]+$", a)) {
      k <- k + as.integer(a)
    } else {
      vterms <- c(vterms, cpp_var(a, vidx, cpp_t0()))
    }
  }
  list(var = var, k = k, vterms = vterms)
}

#' Translate one parsed R expression into a C++ expression
#'
#' @param e a language object, symbol or constant
#' @param vidx named integer vector: variable name -> 0-based column index
#' @param toff time offset (see `cpp_t0`)
#' @return a character scalar of C++ code
#' @keywords internal
cpp_expr <- function(e, vidx, toff = cpp_t0()) {

  ## --- constants -----------------------------------------------------------
  if (is.numeric(e)) {
    if (length(e) != 1L) stop("Unexpected vector constant in a model formula.")
    ## always emit a double literal so integer division never happens in C++
    ## (format() alone gives "1" for 1.0, so add the ".0" ourselves)
    s <- sprintf("%.17g", as.double(e))
    if (!grepl("[.eEn]", s)) s <- paste0(s, ".0")
    return(s)
  }

  ## --- symbols -------------------------------------------------------------
  if (is.symbol(e)) {
    nm <- as.character(e)
    lg <- cpp_decode_lag(nm, vidx)
    if (!is.null(lg)) {
      return(cpp_var(lg$var, vidx, cpp_tshift(toff, lg$k, lg$vterms)))
    }
    return(cpp_var(nm, vidx, toff))
  }

  if (!is.call(e)) stop("Unsupported element in a model formula: ", class(e))

  fn <- as.character(e[[1L]])
  args <- as.list(e)[-1L]

  ## --- grouping ------------------------------------------------------------
  if (fn == "(") return(paste0("(", cpp_expr(args[[1L]], vidx, toff), ")"))

  ## --- arithmetic ----------------------------------------------------------
  if (fn %in% c("+", "-", "*", "/")) {
    if (length(args) == 1L) {  # unary + / -
      return(paste0("(", fn, cpp_expr(args[[1L]], vidx, toff), ")"))
    }
    return(paste0("(", cpp_expr(args[[1L]], vidx, toff), " ", fn, " ",
                  cpp_expr(args[[2L]], vidx, toff), ")"))
  }

  if (fn == "^") {
    return(paste0("std::pow(", cpp_expr(args[[1L]], vidx, toff), ", ",
                  cpp_expr(args[[2L]], vidx, toff), ")"))
  }

  ## --- delta(n, x) = x - lag(x, n) ----------------------------------------
  ## Handled here rather than expanded earlier so that the shift applies to
  ## every variable reference inside `x`, however deeply nested.
  if (fn %in% c("delta", "newdiff")) {
    if (length(args) < 2L) stop("delta() needs two arguments: delta(n, x).")
    n <- args[[1L]]
    if (!is.numeric(n)) stop("The first argument of delta() must be a literal integer.")
    n <- as.integer(n)
    now   <- cpp_expr(args[[2L]], vidx, toff)
    before <- cpp_expr(args[[2L]], vidx, cpp_tshift(toff, n))
    return(paste0("(", now, " - ", before, ")"))
  }

  ## --- lag(x, n), should normally already be encoded as a symbol -----------
  if (fn %in% c("lag", "mylg")) {
    if (length(args) < 2L) stop("lag() needs two arguments: lag(x, n).")
    n <- args[[2L]]
    if (is.numeric(n)) {
      return(cpp_expr(args[[1L]], vidx, cpp_tshift(toff, as.integer(n))))
    }
    vt <- cpp_expr(n, vidx, cpp_t0())
    return(cpp_expr(args[[1L]], vidx, cpp_tshift(toff, 0L, vt)))
  }

  ## --- plain maths ---------------------------------------------------------
  if (fn %in% names(cpp_unary_map)) {
    if (fn == "log" && length(args) == 2L) {  # log(x, base)
      return(paste0("(std::log(", cpp_expr(args[[1L]], vidx, toff), ") / std::log(",
                    cpp_expr(args[[2L]], vidx, toff), "))"))
    }
    return(paste0(cpp_unary_map[[fn]], "(", cpp_expr(args[[1L]], vidx, toff), ")"))
  }

  if (fn == "logb") {
    return(paste0("(std::log(", cpp_expr(args[[1L]], vidx, toff), ") / std::log(",
                  cpp_expr(args[[2L]], vidx, toff), "))"))
  }

  if (fn == "sign") {
    a <- cpp_expr(args[[1L]], vidx, toff)
    return(paste0("(((", a, ") > 0.0) - ((", a, ") < 0.0))"))
  }

  stop(sprintf("Function '%s' is not supported by the C++ code generator.", fn))
}

#' Translate a formula string into a C++ expression
#'
#' @param text character scalar, a formula in thoR's `new_formula` syntax
#' @param vidx named integer vector: variable name -> 0-based column index
#' @return character scalar of C++ code
#' @keywords internal
cpp_from_formula <- function(text, vidx) {
  e <- parse(text = text, keep.source = FALSE)
  if (length(e) != 1L) stop("A formula must be a single expression: ", text)
  cpp_expr(e[[1L]], vidx, cpp_t0())
}
