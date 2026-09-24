## Sparse symbolic jacobian.
##
## `symbolic_jacobian()` materialises an n x n character matrix and the
## downstream writers then loop over all n^2 cells. On ThreeME 8x8 the heart
## block alone is 2517 x 2517 = 6.3M cells holding 9k real entries; at 25k
## equations that representation cannot be built at all.
##
## Here the jacobian is kept as triplets (i, j, expr) from the start. Nothing
## proportional to n^2 is ever allocated.

#' Compute a block's symbolic jacobian in triplet form
#'
#' Only the partial derivatives that can be non-zero are computed: for each
#' equation, we differentiate with respect to the contemporaneous endogenous
#' variables that actually occur in it.
#'
#' @param equations_list_df data.frame of equations
#' @param eqns_vars_list named list: equation id -> variables occurring in it
#' @param endo_vec character vector of the block's endogenous variables
#' @param equations_subset character vector of the block's equation ids
#' @param id_col column of `equations_list_df` holding the equation id
#' @param formula_col column holding the formula to differentiate
#' @return a `thor_sparse_jacobian`: list with `n`, `equations`, `endo`,
#'   `i`, `j` (1-based, CSC order) and `expr`
#' @keywords internal
symbolic_jacobian_sparse <- function(equations_list_df, eqns_vars_list, endo_vec,
                                     equations_subset, id_col = "id",
                                     formula_col = "new_formula") {

  eqs  <- sort(equations_subset)
  endo <- sort(endo_vec)

  if (length(eqs) != length(endo)) {
    stop(sprintf("Block is not square: %d equations for %d endogenous variables.",
                 length(eqs), length(endo)), call. = FALSE)
  }

  ids <- as.character(equations_list_df[[id_col]])
  pos <- match(eqs, ids)
  if (anyNA(pos)) {
    stop("Some equations of the block were not found in the equation list: ",
         paste(eqs[is.na(pos)], collapse = ", "), call. = FALSE)
  }
  formulas <- as.character(equations_list_df[[formula_col]])[pos]

  ## column lookup: endogenous variable -> index within the block
  endo_col <- stats::setNames(seq_along(endo), endo)

  ## Upper bound on the number of entries, so we fill pre-allocated vectors
  ## instead of growing them equation by equation.
  cand <- lapply(eqs, function(id) {
    intersect(unique(as.character(eqns_vars_list[[id]])), endo)
  })
  cap <- sum(lengths(cand))

  out_i <- integer(cap); out_j <- integer(cap); out_e <- character(cap)
  k <- 0L

  for (r in seq_along(eqs)) {
    vars <- cand[[r]]
    if (length(vars) == 0L) next
    ds <- equation_derivatives(formulas[r], vars)
    ## "0" for variables that only appear inside lag.* atoms
    keep <- !(ds == "0" | ds == "0L")
    nk <- sum(keep)
    if (nk == 0L) next
    at <- k + seq_len(nk)
    out_i[at] <- r
    out_j[at] <- endo_col[vars[keep]]
    out_e[at] <- ds[keep]
    k <- k + nk
  }

  length(out_i) <- k; length(out_j) <- k; length(out_e) <- k

  ## canonical CSC order: by column, then by row. This is the layout Eigen
  ## uses, so the generated code can write straight into the value array.
  o <- order(out_j, out_i)

  structure(list(n         = length(endo),
                 equations = eqs,
                 endo      = endo,
                 i         = out_i[o],
                 j         = out_j[o],
                 expr      = out_e[o]),
            class = "thor_sparse_jacobian")
}

## How an equation is differentiated.
##
## stats::D() first, the Deriv package only when D cannot do it.
##
## Deriv simplifies every derivative it returns, and its simplifier is the
## expensive part: superlinear in the size of the expression, and run once per
## variable. Two measurements made the case for not using it by default:
##
##   * ThreeME's `verif_all`, a root of a sum of squares over every verif_*
##     variable, took 278 s on its own at 13x13 (121 variables) and would have
##     taken about 55 min at 29x33 (286 variables). D: 0.01 s.
##   * Across a whole model the simplifier was 89 of the 115 s the jacobian
##     took on ThreeME 29x33. On 13x13, D for every equation took the R side
##     of the build from 29.8 s to 16.9 s.
##
## D's derivatives are numerically identical (3e-16 on the 13x13 heart), and
## the simplification Deriv adds does not make the generated code smaller:
## the jacobian expressions total 763 KB with D against 762 KB with Deriv on
## 13x13, and the C++ compiler does that tidying anyway.
##
## D knows the arithmetic operators and the usual elementary functions. Two
## things it does not know are dealt with here:
##
##   * `delta(n, x)` / `newdiff(n, x)`. With respect to a current-period
##     variable, d/dv delta(n, x) = dx/dv: the lagged half is made of `lag.*`
##     symbols, which are constants here. So delta(n, x) is replaced by x
##     before differentiating -- the same rule Deriv is given in .onLoad.
##   * anything else (abs, sign, asinh, acosh, atanh, logb, log with a base):
##     the whole equation goes to Deriv, as before.

#' Replace `delta(n, x)` by `x` throughout an expression
#'
#' Only valid for differentiating with respect to a current-period variable;
#' see the note above.
#'
#' @param e a language object
#' @return a language object
#' @keywords internal
strip_delta <- function(e) {
  if (!is.call(e)) return(e)
  fn <- e[[1L]]
  if (is.symbol(fn) && as.character(fn) %in% c("delta", "newdiff") && length(e) >= 3L) {
    return(call("(", strip_delta(e[[3L]])))
  }
  for (k in seq_along(e)[-1L]) e[[k]] <- strip_delta(e[[k]])
  e
}

#' The derivatives of one equation with respect to several variables
#'
#' @param text the equation's residual formula, in `new_formula` syntax
#' @param vars character vector of current-period variable names
#' @return character vector, one derivative per variable ("0" where it vanishes)
#' @keywords internal
equation_derivatives <- function(text, vars) {
  f <- parse(text = text, keep.source = FALSE)
  g <- strip_delta(f[[1L]])
  ## D either handles every function in the equation or none of the calls
  ## below would succeed, so one failure settles it for all the variables.
  out <- tryCatch(
    vapply(vars, function(v) {
      d <- stats::D(g, v)
      ## D leaves constants unevaluated (`1/1000`, `-(1/1000)`). Evaluating
      ## them here lets the generators recognise the entry as constant and
      ## plant it once, instead of recomputing it at every Newton iteration.
      if (!is.numeric(d) && length(all.vars(d)) == 0L) d <- eval(d, baseenv())
      if (is.numeric(d)) sprintf("%.17g", d)
      else deparse1(d, collapse = "", width.cutoff = 500L)
    }, character(1), USE.NAMES = FALSE),
    error = function(e) NULL)
  if (is.null(out)) {
    out <- vapply(vars, function(v) paste(Deriv::Deriv(f, v, cache.exp = FALSE)),
                  character(1), USE.NAMES = FALSE)
  }
  out
}

#' Print a sparse symbolic jacobian
#'
#' @param x a `thor_sparse_jacobian`
#' @param ... ignored
#' @return `x`, invisibly
#' @export
print.thor_sparse_jacobian <- function(x, ...) {
  nnz <- length(x$i)
  cat(sprintf("Sparse symbolic jacobian: %d x %d, %d non-zero (%.4f%% dense), %.2f per row\n",
              x$n, x$n, nnz, 100 * nnz / max(1, x$n^2), nnz / max(1, x$n)))
  invisible(x)
}

#' Convert a dense symbolic jacobian to triplet form
#'
#' Lets the dense backends share the constant/varying entry split and the
#' expression translation with the sparse one.
#'
#' @param m character matrix as returned by [symbolic_jacobian()]
#' @return a `thor_sparse_jacobian`
#' @keywords internal
as_sparse_jacobian <- function(m) {
  nz <- which(m != "0", arr.ind = TRUE)
  o <- order(nz[, "col"], nz[, "row"])
  nz <- nz[o, , drop = FALSE]
  structure(list(n         = ncol(m),
                 equations = rownames(m),
                 endo      = colnames(m),
                 i         = as.integer(nz[, "row"]),
                 j         = as.integer(nz[, "col"]),
                 expr      = m[nz]),
            class = "thor_sparse_jacobian")
}
