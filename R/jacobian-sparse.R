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
    f <- parse(text = formulas[r], keep.source = FALSE)
    for (v in vars) {
      d <- paste(Deriv::Deriv(f, v, cache.exp = FALSE))
      ## Deriv returns "0" for variables that only appear inside lag.* atoms
      if (identical(d, "0") || identical(d, "0L")) next
      k <- k + 1L
      out_i[k] <- r
      out_j[k] <- endo_col[[v]]
      out_e[k] <- d
    }
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
