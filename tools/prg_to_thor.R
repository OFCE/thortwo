## Translate a compiled EViews model program into a thortwo model file.
##
## The compiler emits two things: `model.prg`, a list of `<model>.append <eq>`
## lines in EViews syntax, and `calib.csv`, the calibration database with the
## year expressed relative to the base year. This turns them into the four
## sections thortwo's parser reads, plus a database with the coefficients
## already broadcast onto it.
##
## Deliberately dependency-free: base R plus thortwo. Source it, or run
## tools/run_rl_model.R, which uses it.
##
##   source("tools/prg_to_thor.R")
##   m <- prg_to_thor("model.prg", "calib.csv", base_year = 2019,
##                    out_file = "model.txt")
##
## ---------------------------------------------------------------------------
## Relationship to ermeeth2::translate_modelprg()
##
## Same job, same output conventions (`elem_<expr>_<year>` coefficients,
## endogenous = the first variable on the left-hand side), so a model
## translated either way lands in the same place. Four things are done
## differently, each because the ermeeth version has a failure mode on this
## input:
##
##  1. `@elem` is found by scanning balanced parentheses, not by three
##     hand-written regexes for the three shapes seen so far. ermeeth2 matches
##     `@elem(<word>,<year>)`, `@elem(<word>(-<k>),<year>)` and
##     `@elem(<word><op><word>,<year>)`; anything else -- two operators, a
##     function call, a parenthesised sub-expression -- is silently left in
##     place and then reaches the solver as an undefined variable.
##
##  2. The value is obtained by *evaluating* the inner expression against the
##     calibration row, rather than by a switch over the four arithmetic
##     operators. Any expression the compiler can emit therefore works, and
##     there is one code path instead of three.
##
##  3. Occurrences are keyed on (expression, year) and substituted as literal
##     text, longest first. ermeeth2 keys on the variable name alone, so two
##     `@elem` of the same variable at different years collide and both take
##     the value of whichever was seen last.
##
##  4. `rbind(elem_table_1, elem_table_2)` in ermeeth2 errors outright when a
##     model uses lagged `@elem` but no plain ones (or vice versa), because the
##     missing table was never created. There is one table here.
##
## Comparisons are handled differently too. EViews lets an equation contain a
## logical test that evaluates to 0 or 1 (`(pgd_cind/pgm_cind<0.99999)*...`),
## which thortwo cannot parse. ermeeth2 deals with these in two ways: three
## hard-coded ThreeME substitutions that replace a named test with a literal
## 1 or 0, and a generic rewrite that only matches `<word><op><word><cmp><num>`.
##
## Here every comparison goes through one generic rewrite instead, so no test
## is silently pinned to a constant and no shape is missed -- see
## `prg_rewrite_comparisons()` for the form used and its one caveat.
## ---------------------------------------------------------------------------

## ---- small helpers --------------------------------------------------------

## Functions that may appear in an equation and must not be mistaken for
## variables when rewriting reads of the calibration data.
PRG_FUNCTIONS <- c("log", "exp", "sqrt", "abs", "d", "dlog", "lag", "delta",
                   "newdiff", "sin", "cos", "tan", "log10", "log2")

#' Rewrite EViews logical tests as smooth indicators
#'
#' `(A < B)` evaluates to 1 or 0 in EViews. thortwo has no comparison
#' operators -- and a step function has no useful derivative -- so each test
#' becomes the algebraic indicator
#'
#'   A > B   ->   0.5 * (1 + (A-B) / (|A-B| + eps))
#'   A < B   ->   0.5 * (1 - (A-B) / (|A-B| + eps))
#'
#' which is exactly 1 or 0 away from the crossing, and whose derivative is
#' `eps/(|A-B|+eps)^2`, i.e. numerically zero there. `>=` and `<=` are treated
#' as `>` and `<`; they differ only on the measure-zero set A == B, where this
#' form returns 0.5 rather than the 0/0 that the difference-quotient version
#' used by ermeeth2 returns.
#'
#' **Caveat.** Right at the crossing the indicator is not differentiable, so a
#' model that sits exactly on a threshold can stall. Nothing in the translation
#' can fix that; it is a property of writing a discontinuity into a model that
#' is then solved by Newton.
#'
#' Each test must be parenthesised on its own, which is how the compiler emits
#' them; anything else is reported rather than guessed at.
#'
#' @param eqs character vector of equations
#' @return list(eqs, n) -- the rewritten equations and how many tests were found
#' @keywords internal
prg_rewrite_comparisons <- function(eqs) {

  eps <- "1e-30"
  n <- 0L

  for (k in seq_along(eqs)) {
    repeat {
      m <- regexpr("(<=|>=|<|>)", eqs[k])
      if (m == -1L) break
      op  <- regmatches(eqs[k], m)
      pos <- as.integer(m)
      len <- attr(m, "match.length")
      s   <- eqs[k]

      ## innermost enclosing parentheses
      left <- 0L; depth <- 0L
      for (i in seq(pos - 1L, 1L)) {
        ch <- substr(s, i, i)
        if (ch == ")") depth <- depth + 1L
        else if (ch == "(") {
          if (depth == 0L) { left <- i; break }
          depth <- depth - 1L
        }
      }
      right <- 0L; depth <- 0L
      for (i in seq(pos + len, nchar(s))) {
        ch <- substr(s, i, i)
        if (ch == "(") depth <- depth + 1L
        else if (ch == ")") {
          if (depth == 0L) { right <- i; break }
          depth <- depth - 1L
        }
      }
      if (left == 0L || right == 0L) {
        stop("A comparison is not enclosed in its own parentheses, so its ",
             "operands cannot be identified:\n  ", s, call. = FALSE)
      }

      A <- substr(s, left + 1L, pos - 1L)
      B <- substr(s, pos + len, right - 1L)
      sign <- if (substr(op, 1L, 1L) == ">") "+" else "-"
      diff <- paste0("((", A, ")-(", B, "))")
      new  <- sprintf("(0.5*(1.0%s%s/(abs%s+%s)))", sign, diff, diff, eps)

      eqs[k] <- paste0(substr(s, 1L, left - 1L), new, substring(s, right + 1L))
      n <- n + 1L
    }
  }
  list(eqs = eqs, n = n)
}

#' Expand scientific notation into plain decimals
#'
#' `1e-05` survives R's parser but not thortwo's variable extractor, which
#' would read a stray `e` as a variable name.
#'
#' @param eqs character vector of equations
#' @return the equations with every literal written out in full
#' @keywords internal
prg_expand_scientific <- function(eqs) {
  out <- eqs
  for (k in seq_along(out)) {
    m <- gregexpr("[0-9]+(\\.[0-9]+)?e[-+]?[0-9]+", out[k], perl = TRUE)[[1]]
    if (m[1] == -1L) next
    lits <- regmatches(out[k], gregexpr("[0-9]+(\\.[0-9]+)?e[-+]?[0-9]+", out[k], perl = TRUE))[[1]]
    for (lit in unique(lits)) {
      out[k] <- gsub(lit, format(as.numeric(lit), scientific = FALSE, digits = 17),
                     out[k], fixed = TRUE)
    }
  }
  out
}

#' Find every `@elem(...)` occurrence, honouring nested parentheses
#'
#' A regular expression cannot do this: `@elem(pk_scon(-1), 2019)` closes on
#' its second `)`, not its first.
#'
#' @param x character vector of equations
#' @return character vector of the distinct occurrences, as written
#' @keywords internal
prg_find_elem <- function(x) {
  out <- character(0)
  for (line in x) {
    start <- gregexpr("@elem(", line, fixed = TRUE)[[1]]
    if (start[1] == -1L) next
    for (s in start) {
      depth <- 0L
      i <- s + 5L                       # position of the opening "("
      n <- nchar(line)
      repeat {
        ch <- substr(line, i, i)
        if (ch == "(") depth <- depth + 1L
        if (ch == ")") {
          depth <- depth - 1L
          if (depth == 0L) break
        }
        i <- i + 1L
        if (i > n) stop("Unbalanced parentheses in: ", line, call. = FALSE)
      }
      out <- c(out, substr(line, s, i))
    }
  }
  unique(out)
}

#' Split `@elem(<expr>, <year>)` into its expression and its year
#'
#' @param e a single occurrence, as returned by [prg_find_elem()]
#' @return list(expr, year)
#' @keywords internal
prg_parse_elem <- function(e) {
  inner <- substr(e, 7L, nchar(e) - 1L)          # strip "@elem(" and ")"
  at <- max(gregexpr(",", inner, fixed = TRUE)[[1]])
  if (at == -1L) stop("@elem without a year: ", e, call. = FALSE)
  list(expr = trimws(substr(inner, 1L, at - 1L)),
       year = as.integer(trimws(substring(inner, at + 1L))))
}

#' Name for the coefficient an `@elem` becomes
#'
#' Follows the convention already in the translated ThreeME models:
#' `@elem(chd_cind/chm_cind, 2015)` becomes `elem_chd_cind_chm_cind_2015`, and
#' a lag is folded into the year, so `@elem(pk_sind(-1), 2019)` becomes
#' `elem_pk_sind_2018`.
#'
#' @param expr the inner expression
#' @param year the year
#' @return a valid lower-case variable name
#' @keywords internal
prg_elem_name <- function(expr, year) {
  e <- gsub("\\s+", "", expr)
  lag_only <- regmatches(e, regexec("^([a-z][a-z0-9_]*)\\(-([0-9]+)\\)$", e))[[1]]
  if (length(lag_only) == 3L) {
    return(sprintf("elem_%s_%d", lag_only[2], year - as.integer(lag_only[3])))
  }
  slug <- gsub("_+", "_", gsub("[^a-z0-9_]", "_", tolower(e)))
  sprintf("elem_%s_%d", gsub("^_|_$", "", slug), year)
}

#' Evaluate an EViews expression against one row of the calibration data
#'
#' Variable reads become lookups by year, so a lag inside the expression is
#' honoured: in `@elem(a/b(-1), 2019)`, `a` is read at 2019 and `b` at 2018.
#'
#' @param expr the inner expression
#' @param year the year to read at
#' @param calib the calibration data.frame, with an absolute `year` column
#' @return a numeric scalar, or NA if a variable is absent from `calib`
#' @keywords internal
prg_eval_elem <- function(expr, year, calib) {

  e <- gsub("\\s+", "", expr)
  ## `x(-k)` -> `.v("x", year - k)`
  e <- gsub("([a-z][a-z0-9_]*)\\(-([0-9]+)\\)", '.v("\\1",.y-\\2)', e)
  ## remaining bare names -> `.v("x", year)`, skipping function calls and
  ## anything already rewritten
  e <- gsub('(?<![."a-z0-9_])([a-z][a-z0-9_]*)(?!\\s*\\(|["a-z0-9_])',
            '.v("\\1",.y)', e, perl = TRUE)

  .y <- year
  .v <- function(nm, yr) {
    if (!nm %in% names(calib)) return(NA_real_)
    row <- which(calib$year == yr)
    if (length(row) != 1L) return(NA_real_)
    as.numeric(calib[[nm]][row])
  }
  env <- list2env(list(.v = .v, .y = .y), parent = baseenv())

  tryCatch(eval(parse(text = e, keep.source = FALSE), envir = env),
           error = function(err) NA_real_)
}

## ---- the translator -------------------------------------------------------

#' Translate a compiled model program into a thortwo model
#'
#' @param prg_file path to the compiled `.prg`
#' @param calib_file path to the calibration `.csv`, whose `year` column is
#'   relative to `base_year`
#' @param base_year the calendar year that relative year 0 corresponds to
#' @param out_file where to write the thortwo `.txt`. NULL to skip writing.
#' @param first_year,last_year optional bounds to trim the database to
#' @param model_prefix the EViews model object name prefixed to each `.append`
#'   line; the default matches anything of the form `<name>.append`
#' @param verbose report progress
#'
#' @return a list with `equations`, `endo`, `exo`, `coeff`, `data`, `elem`
#'   (the coefficient table) and `warnings`
prg_to_thor <- function(prg_file,
                        calib_file,
                        base_year,
                        out_file = NULL,
                        first_year = NULL,
                        last_year = NULL,
                        model_prefix = "^[a-z_0-9]+\\.append",
                        verbose = TRUE) {

  say <- function(...) if (verbose) cat(..., sep = "")
  warnings_out <- character(0)
  warn <- function(msg) {
    warnings_out <<- c(warnings_out, msg)
    say("   ! ", msg, "\n")
  }

  stopifnot(file.exists(prg_file), file.exists(calib_file))

  ## ---- 1. equations -------------------------------------------------------
  say("1. reading ", basename(prg_file), "\n")
  raw <- readLines(prg_file, warn = FALSE)
  eqs <- tolower(raw)
  eqs <- sub(model_prefix, "", eqs)
  eqs <- gsub("'.*$", "", eqs)                    # EViews end-of-line comments
  eqs <- gsub("\\s+", "", eqs)
  eqs <- eqs[nzchar(eqs)]
  say("   ", length(eqs), " equations\n")

  ## ---- 2. calibration -----------------------------------------------------
  say("2. reading ", basename(calib_file), "\n")
  calib <- utils::read.csv(calib_file, check.names = FALSE)
  names(calib) <- tolower(names(calib))
  calib <- calib[, names(calib) != "baseyear", drop = FALSE]
  if (!"year" %in% names(calib)) {
    stop("The calibration file has no `year` column.", call. = FALSE)
  }
  calib$year <- as.integer(round(calib$year)) + base_year
  if (!is.null(first_year)) calib <- calib[calib$year >= first_year, , drop = FALSE]
  if (!is.null(last_year))  calib <- calib[calib$year <= last_year,  , drop = FALSE]
  rownames(calib) <- NULL
  say("   ", nrow(calib), " periods (", min(calib$year), "-", max(calib$year),
      "), ", ncol(calib) - 1L, " variables\n")

  ## ---- 3. @elem -> coefficients ------------------------------------------
  say("3. resolving @elem\n")
  occ <- prg_find_elem(eqs)
  elem <- NULL

  if (length(occ)) {
    parsed <- lapply(occ, prg_parse_elem)
    elem <- data.frame(
      original = occ,
      expr     = vapply(parsed, `[[`, character(1), "expr"),
      year     = vapply(parsed, `[[`, integer(1),   "year"),
      stringsAsFactors = FALSE
    )
    elem$name  <- mapply(prg_elem_name, elem$expr, elem$year, USE.NAMES = FALSE)
    elem$value <- mapply(prg_eval_elem, elem$expr, elem$year,
                         MoreArgs = list(calib = calib), USE.NAMES = FALSE)

    say("   ", nrow(elem), " distinct @elem -> ",
        length(unique(elem$name)), " coefficients\n")

    bad <- elem[is.na(elem$value), , drop = FALSE]
    if (nrow(bad)) {
      warn(sprintf("%d @elem could not be evaluated against the calibration (e.g. %s)",
                   nrow(bad), paste(utils::head(bad$original, 3), collapse = ", ")))
    }

    ## Two occurrences that reduce to the same name must reduce to the same
    ## number, or the name is ambiguous.
    dup <- tapply(elem$value, elem$name, function(v) length(unique(round(v, 12))) > 1L)
    if (any(dup, na.rm = TRUE)) {
      warn(sprintf("coefficient name(s) map to more than one value: %s",
                   paste(names(dup)[which(dup)], collapse = ", ")))
    }

    ## Substituted as literal text, longest first, so that no occurrence can be
    ## rewritten inside another one.
    ord <- order(nchar(elem$original), decreasing = TRUE)
    for (k in ord) {
      eqs <- gsub(elem$original[k], elem$name[k], eqs, fixed = TRUE)
    }

    left <- sum(grepl("@elem", eqs, fixed = TRUE))
    if (left) warn(sprintf("%d @elem left in the equations after substitution", left))

    ## onto the database, as constant columns
    coeff_tbl <- elem[!duplicated(elem$name), c("name", "value")]
    calib <- calib[, setdiff(names(calib), coeff_tbl$name), drop = FALSE]
    calib <- cbind(calib,
                   as.data.frame(matrix(rep(coeff_tbl$value, each = nrow(calib)),
                                        nrow = nrow(calib),
                                        dimnames = list(NULL, coeff_tbl$name))))
  } else {
    say("   none\n")
  }

  ## ---- 4. EViews -> thortwo syntax ---------------------------------------
  say("4. rewriting lags and differences\n")
  ## `x(-1)` -> `lag(x,1)`. Done before `d(` so that a lag inside a difference
  ## is already in thortwo form.
  eqs <- gsub("([a-z][a-z0-9_]*)\\(-([0-9]+)\\)", "lag(\\1,\\2)", eqs)
  ## `d(...)` -> `delta(1,...)`, but not the `d` of an identifier ending in d
  eqs <- gsub("(?<![a-z0-9_])d\\(", "delta(1,", eqs, perl = TRUE)
  ## `dlog(x)` would have become `delta(1,og(x)` above; it does not occur in
  ## this compiler's output, but catch it rather than emit nonsense
  if (any(grepl("delta(1,og(", eqs, fixed = TRUE))) {
    stop("dlog() is not handled by this translator.", call. = FALSE)
  }
  eqs <- gsub("+-", "-", eqs, fixed = TRUE)

  ## ---- 4b. logical tests and scientific notation --------------------------
  n_cmp <- sum(grepl("[<>]", eqs))
  if (n_cmp) {
    rc <- prg_rewrite_comparisons(eqs)
    eqs <- rc$eqs
    say("   rewrote ", rc$n, " logical test(s) in ", n_cmp, " equation(s) as indicators\n")
  }
  eqs <- prg_expand_scientific(eqs)

  leftover <- grep("@|<|>", eqs)
  if (length(leftover)) {
    warn(sprintf("%d equation(s) still contain @, < or > and will not parse (e.g. %s)",
                 length(leftover), utils::head(eqs[leftover], 1)))
  }

  ## ---- 5. classify the variables -----------------------------------------
  say("5. classifying variables\n")
  ## The endogenous variable of an equation is the first one on its left-hand
  ## side -- the compiler emits equations already normalised that way.
  lhs <- sub("=.*$", "", eqs)
  endo <- vapply(lhs, function(s) {
    v <- thortwo::get_variables_from_string(s)
    if (length(v) == 0L) NA_character_ else v[1L]
  }, character(1), USE.NAMES = FALSE)

  if (anyNA(endo)) {
    warn(sprintf("%d equation(s) have no variable on the left-hand side", sum(is.na(endo))))
    endo <- endo[!is.na(endo)]
  }
  if (anyDuplicated(endo)) {
    d <- unique(endo[duplicated(endo)])
    warn(sprintf("%d variable(s) are the left-hand side of more than one equation (e.g. %s)",
                 length(d), paste(utils::head(d, 3), collapse = ", ")))
  }

  all_vars <- thortwo::get_variables_from_string(paste(eqs, collapse = "+"))
  coeff <- if (is.null(elem)) character(0) else sort(unique(elem$name))
  coeff <- intersect(coeff, all_vars)
  endo  <- sort(unique(endo))
  exo   <- sort(setdiff(all_vars, c(endo, coeff)))

  say("   ", length(endo), " endogenous, ", length(exo), " exogenous, ",
      length(coeff), " coefficients\n")
  if (length(endo) != length(eqs)) {
    warn(sprintf("%d equations for %d endogenous variables: the model is not square",
                 length(eqs), length(endo)))
  }

  missing <- setdiff(c(endo, exo, coeff), names(calib))
  if (length(missing)) {
    warn(sprintf("%d model variable(s) are absent from the calibration (e.g. %s)",
                 length(missing), paste(utils::head(missing, 5), collapse = ", ")))
  }

  ## ---- 6. write the model file -------------------------------------------
  if (!is.null(out_file)) {
    say("6. writing ", basename(out_file), "\n")
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    writeLines(c(
      "endogenous variables :", paste(endo,  collapse = ","), "##############",
      "exogenous variables :",  paste(exo,   collapse = ","), "##############",
      "coefficients :",         paste(coeff, collapse = ","), "##############",
      "equations :",            eqs), out_file)
  }

  invisible(list(equations = eqs, endo = endo, exo = exo, coeff = coeff,
                 data = calib, elem = elem, warnings = warnings_out,
                 file = out_file))
}

#' Print what a translation did and what it could not do
#'
#' @param tr the value of [prg_to_thor()]
#' @return `tr`, invisibly
translate_report <- function(tr) {
  cat("\n--- translation report ---\n")
  cat("equations   :", length(tr$equations), "\n")
  cat("endogenous  :", length(tr$endo), "\n")
  cat("exogenous   :", length(tr$exo), "\n")
  cat("coefficients:", length(tr$coeff), "\n")
  cat("database    :", nrow(tr$data), "periods x", ncol(tr$data) - 1L, "variables\n")
  if (length(tr$warnings) == 0L) {
    cat("warnings    : none\n")
  } else {
    cat("warnings    :\n")
    cat(paste0("  - ", tr$warnings, collapse = "\n"), "\n")
  }
  invisible(tr)
}
