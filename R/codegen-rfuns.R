## The R backend (`sparse-r`): generate R closures that evaluate a block's
## residuals and jacobian values, for users with no compiler.
##
## Two things changed relative to tresthor's 1_4_function_command_writers.R.
##
## D1 -- the build no longer writes to the filesystem at all. tresthor did
## `dir.create("temp_paprfn",)` in the *current working directory* (note the
## stray comma), wrote six .R files into it, `source()`d them and `unlink`ed
## the directory at the end. A failed build left litter behind, two concurrent
## builds in one directory raced, and a read-only working directory broke the
## build outright. The generated text is evaluated with `eval(parse())` into a
## dedicated environment instead.
##
## The commented-out attempt at that in tresthor (1_0_create_model.R:251-256)
## was abandoned because it kept `eval(parse(text = ...))` without
## `keep.source = FALSE`: R then attaches a srcref to every one of the
## thousands of generated statements, which on Opale costs more memory and
## time than the round-trip through disk it was meant to replace. With source
## references switched off the approach is strictly better, and that is what
## is done here.
##
## The second change is that the translation itself is an AST walk over the
## parsed formula, mirroring `cpp_expr()`, rather than a chain of a dozen
## regular expressions rewriting the formula as text. The two generators then
## read the same formula the same way by construction.

#' Translate one parsed R expression into R source against the data matrix
#'
#' The generated code reads `M[t - k, j]`, where `M` is the numeric data
#' matrix whose columns are the model variables in alphabetical order -- the
#' same layout the C++ backends use.
#'
#' @param e a language object, symbol or constant
#' @param vidx named integer vector: variable name -> 1-based column index
#' @param toff time offset, as built by [cpp_t0()]
#' @return character scalar of R code
#' @keywords internal
r_expr <- function(e, vidx, toff = cpp_t0()) {

  render_t <- function(toff) {
    out <- "t"
    if (toff$k != 0L) out <- paste0(out, " - ", toff$k)
    for (v in toff$vterms) out <- paste0(out, " - round(", v, ")")
    out
  }
  read_var <- function(name, toff) {
    idx <- var_col(vidx, name)
    if (is.null(idx) || is.na(idx)) {
      stop(sprintf("Variable '%s' is used in the model but is not a known model variable.",
                   name), call. = FALSE)
    }
    paste0("M[", render_t(toff), ",", idx, "]")
  }

  if (is.numeric(e)) {
    if (length(e) != 1L) stop("Unexpected vector constant in a model formula.", call. = FALSE)
    return(format(as.double(e), digits = 17L, scientific = FALSE, trim = TRUE))
  }

  if (is.symbol(e)) {
    nm <- as.character(e)
    lg <- cpp_decode_lag_r(nm, vidx, read_var)
    if (!is.null(lg)) return(read_var(lg$var, cpp_tshift(toff, lg$k, lg$vterms)))
    return(read_var(nm, toff))
  }

  if (!is.call(e)) stop("Unsupported element in a model formula: ", class(e), call. = FALSE)

  fn <- as.character(e[[1L]])
  args <- as.list(e)[-1L]

  if (fn == "(") return(paste0("(", r_expr(args[[1L]], vidx, toff), ")"))

  ## No parentheses are added around a binary operation. Every formula reaches
  ## here through parse(), so the tree already holds a `(` node wherever the
  ## text needed one, and printing the operands as they are reproduces the
  ## original precedence. Adding them anyway nested a long sum one level per
  ## term -- ThreeME 13x13's `verif_all`, 121 terms, went 121 deep -- past what
  ## R's parser accepts ("contextstack overflow").
  if (fn %in% c("+", "-", "*", "/", "^")) {
    if (length(args) == 1L) return(paste0("(", fn, r_expr(args[[1L]], vidx, toff), ")"))
    return(paste0(r_expr(args[[1L]], vidx, toff), " ", fn, " ",
                  r_expr(args[[2L]], vidx, toff)))
  }

  ## delta(n, x) = x - lag(x, n); the shift applies to every variable
  ## reference inside x, however deeply nested.
  if (fn %in% c("delta", "newdiff")) {
    if (length(args) < 2L) stop("delta() needs two arguments: delta(n, x).", call. = FALSE)
    if (!is.numeric(args[[1L]])) {
      stop("The first argument of delta() must be a literal integer.", call. = FALSE)
    }
    n <- as.integer(args[[1L]])
    ## Both halves are parenthesised: this is the one place where an operand
    ## is built here rather than taken from the parsed text, so nothing else
    ## guarantees `x - (x shifted)` keeps its grouping when x is a sum.
    return(paste0("((", r_expr(args[[2L]], vidx, toff), ") - (",
                  r_expr(args[[2L]], vidx, cpp_tshift(toff, n)), "))"))
  }

  if (fn %in% c("lag", "mylg")) {
    if (length(args) < 2L) stop("lag() needs two arguments: lag(x, n).", call. = FALSE)
    n <- args[[2L]]
    if (is.numeric(n)) return(r_expr(args[[1L]], vidx, cpp_tshift(toff, as.integer(n))))
    return(r_expr(args[[1L]], vidx, cpp_tshift(toff, 0L, r_expr(n, vidx, cpp_t0()))))
  }

  ## Everything else maps onto the R function of the same name, which is what
  ## the formula meant in the first place.
  if (fn %in% thor_functions_supported) {
    return(paste0(fn, "(", paste(vapply(args, r_expr, character(1), vidx = vidx, toff = toff),
                                 collapse = ", "), ")"))
  }

  stop(sprintf("Function '%s' is not supported by the R code generator.", fn), call. = FALSE)
}

#' Decode a `lag.<var>.<n>` symbol for the R generator
#'
#' Same encoding as [cpp_decode_lag()]; only the way a variable lag amount is
#' rendered differs, so the reader is passed in.
#'
#' @param nm symbol name
#' @param vidx variable index
#' @param read_var function(name, toff) rendering a variable read
#' @return NULL if `nm` is not a lag symbol, otherwise list(var, k, vterms)
#' @keywords internal
cpp_decode_lag_r <- function(nm, vidx, read_var) {
  if (!startsWith(nm, "lag.")) return(NULL)
  parts <- strsplit(nm, ".", fixed = TRUE)[[1L]]
  if (length(parts) < 3L) stop(sprintf("Malformed lag symbol '%s'.", nm), call. = FALSE)
  var <- parts[2L]
  k <- 0L
  vterms <- character(0)
  for (a in parts[-c(1L, 2L)]) {
    if (grepl("^[0-9]+$", a)) k <- k + as.integer(a)
    else vterms <- c(vterms, read_var(a, cpp_t0()))
  }
  list(var = var, k = k, vterms = vterms)
}

#' Translate a formula string into R source
#'
#' @param text character scalar, a formula in `new_formula` syntax
#' @param vidx named integer vector: variable name -> 1-based column index
#' @return character scalar of R code
#' @keywords internal
r_from_formula <- function(text, vidx) {
  e <- parse(text = text, keep.source = FALSE)
  if (length(e) != 1L) stop("A formula must be a single expression: ", text, call. = FALSE)
  r_expr(e[[1L]], vidx, cpp_t0())
}

#' Generate the R source of a model's block functions
#'
#' Produces, for each block, `<block>_res(t, M)` returning the residual vector,
#' `<block>_jac(t, M, x)` overwriting the varying entries of the vector of
#' jacobian values, and `<block>_seed()` returning that vector with its
#' constant entries planted; plus a `blocks` list tying them together with
#' the jacobian's sparsity pattern.
#'
#' The values are in the jacobian's CSC order (by column, then row), which is
#' exactly the order of the `x` slot of a `Matrix::dgCMatrix` built from the
#' pattern -- so the solver fills the matrix by assigning one vector, with no
#' sorting or triplet conversion per iteration.
#'
#' @param model_name name of the model
#' @param blocks named list of `list(jac, formulas)`
#' @param all_model_vars character vector of every model variable, in the
#'   model's order (`model@vars$all`)
#' @param verbose print one line per block
#' @return character scalar: the full generated R source
#' @keywords internal
generate_r_source <- function(model_name, blocks, all_model_vars, verbose = TRUE) {

  vidx <- var_index(all_model_vars, base = 1L)

  active <- names(blocks)[vapply(blocks, function(b) !is.null(b) && b$jac$n > 0L,
                                 logical(1))]

  lines <- c(
    paste0("## Generated by thortwo for model '", model_name, "' (R backend). Do not edit."),
    ""
  )

  for (nm in active) {
    b <- blocks[[nm]]
    jac <- b$jac
    if (verbose) {
      cat("   - block '", nm, "': ", jac$n, " equations, ",
          if (is.null(b$seq)) paste(length(jac$i), "jacobian entries")
          else "solved sequentially", "\n", sep = "")
    }

    sp <- cpp_split_constants(jac)

    res_stmts <- vapply(seq_along(b$formulas), function(k) {
      paste0("  f[", k, "] <- ", r_from_formula(b$formulas[k], vidx))
    }, character(1))

    jac_stmts <- vapply(sp$vary, function(k) {
      paste0("  x[", k, "] <- ", r_from_formula(jac$expr[k], vidx))
    }, character(1))

    lines <- c(lines,
      paste0(nm, "_res <- function(t, M) {"),
      paste0("  f <- numeric(", jac$n, ")"),
      res_stmts,
      "  f", "}", "",
      paste0(nm, "_jac <- function(t, M, x) {"),
      jac_stmts,
      "  x", "}", "",
      paste0(nm, "_seed <- function() {"),
      paste0("  x <- numeric(", length(jac$i), ")"),
      if (length(sp$cst_v)) {
        paste0("  x[c(", paste(which(sp$is_const), collapse = ","), ")] <- c(",
               paste(format(sp$cst_v, digits = 17L, scientific = FALSE, trim = TRUE),
                     collapse = ","), ")")
      },
      "  x", "}", ""
    )
    if (!is.null(b$seq)) lines <- c(lines, r_emit_sequential(nm, b$seq, vidx))
  }

  ## column pointers of the CSC pattern, 0-based as dgCMatrix stores them
  colptr <- function(jac) c(0L, cumsum(tabulate(jac$j, nbins = jac$n)))

  lines <- c(lines,
    "blocks <- list(",
    paste0("  ", vapply(active, function(nm) {
      jac <- blocks[[nm]]$jac
      seq <- blocks[[nm]]$seq
      quoted <- function(x) paste0('"', x, '"', collapse = ",")
      sprintf(paste0('%s = list(name = "%s", n = %dL, endo = c(%s), res = %s_res, ',
                     'jac = %s_jac, seed = %s_seed, i = c(%s), p = c(%s)%s)'),
              nm, nm, jac$n,
              paste(as.integer(vidx[jac$endo]), collapse = ","),
              nm, nm, nm,
              if (length(jac$i)) paste0(as.integer(jac$i) - 1L, "L", collapse = ",") else "",
              paste0(colptr(jac), "L", collapse = ","),
              if (is.null(seq)) "" else
                sprintf(", seq = %s_seq, seq_eq = c(%s), seq_var = c(%s)", nm,
                        quoted(seq$order$name), quoted(seq$order$variable)))
    }, character(1)), collapse = ",\n"),
    ")", ""
  )

  paste0(paste(lines, collapse = "\n"), "\n")
}
