## Dense symbolic jacobian: an n x n character matrix of derivative
## expressions, "0" where the derivative is structurally zero.
##
## Ported from tresthor's 1_5_symbolic_jacobian.R. Only the purrr calls are
## gone, and the matrix is filled in place rather than rebuilt row by row
## through a tibble.

#' Compute a block's symbolic jacobian as a dense character matrix
#'
#' Only the entries that can be non-zero are differentiated: for each
#' equation, with respect to the contemporaneous endogenous variables that
#' actually occur in it.
#'
#' @param equations_list_df data.frame of equations
#' @param eqns_vars_list named list: equation id -> variables occurring in it
#' @param endo_vec character vector of the block's endogenous variables
#' @param equations_subset character vector of the block's equation ids
#' @param id_col column of `equations_list_df` holding the equation id
#' @param formula_col column holding the formula to differentiate
#' @return character matrix, equations by endogenous variables, with the
#'   equation ids as row names and the variable names as column names
#' @keywords internal
symbolic_jacobian <- function(equations_list_df, eqns_vars_list, endo_vec,
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

  jacobian <- matrix("0", nrow = length(eqs), ncol = length(endo),
                     dimnames = list(eqs, endo))

  for (r in seq_along(eqs)) {
    vars <- intersect(unique(as.character(eqns_vars_list[[eqs[r]]])), endo)
    if (length(vars) == 0L) next
    f <- parse(text = formulas[r], keep.source = FALSE)
    for (v in vars) {
      jacobian[r, v] <- paste(Deriv::Deriv(f, v, cache.exp = FALSE))
    }
  }

  jacobian
}
