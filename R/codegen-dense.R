## Generate the model-specific C++ solver with dense jacobians.
##
## This is the classic backend. tresthor generated it by writing R source to a
## temporary directory, reading it back as text and rewriting it with regular
## expressions (1_7_rcpp_source_builder.R), then solved with RcppArmadillo's
## `inv(J) * f`. Neither survives here:
##
##   * the expressions come from the shared AST translator in codegen-expr.R,
##     so the dense and sparse backends cannot disagree about what a formula
##     means, and nothing is round-tripped through the filesystem;
##   * the linear algebra is Eigen's `PartialPivLU`, which drops the
##     RcppArmadillo dependency entirely and never forms an explicit inverse.
##
## What stays dense is the jacobian itself: one n x n matrix per block,
## allocated once at setup and overwritten in place.

## The dense block type. `newton_block` in CPP_PREAMBLE is templated over it.
CPP_DENSE_BLOCK <- '
struct Block {
  const char* name;
  int n, nvar, ncst;
  const int* var_i; const int* var_j;          // varying entries
  const int* cst_i; const int* cst_j; const double* cst_v;
  const int* endo;
  void (*jac)(const MapMat&, int, double*);
  void (*res)(const MapMat&, int, double*);

  Eigen::MatrixXd A;
  Eigen::PartialPivLU<Eigen::MatrixXd> lu;
  std::vector<double> v;

  Block() {}

  // The zeros and the constant entries are written once; a Newton iteration
  // then only touches the entries that actually vary.
  void setup(){
    A = Eigen::MatrixXd::Zero(n, n);
    for (int c = 0; c < ncst; ++c) A(cst_i[c], cst_j[c]) = cst_v[c];
    v.assign(nvar > 0 ? nvar : 1, 0.0);
  }

  void refresh(const MapMat& M, int t){
    jac(M, t, v.data());
    for (int k = 0; k < nvar; ++k) A(var_i[k], var_j[k]) = v[k];
  }
  bool factorize(){ lu.compute(A); return true; }
  bool solve(const Eigen::VectorXd& f, Eigen::VectorXd& dx){
    dx = lu.solve(f);
    // PartialPivLU never reports failure: a singular jacobian surfaces as
    // inf/nan in the step, so that is what gets checked.
    return dx.allFinite();
  }
};
'

#' Emit the C++ for one dense block
#'
#' @param name block name ("prologue", "heart", "epilogue")
#' @param jac a `thor_sparse_jacobian`, obtained from the dense character
#'   matrix by [as_sparse_jacobian()]. Only the non-zero entries are ever
#'   written out; the matrix they are written into is dense.
#' @param formulas character vector of the block's residual formulas
#' @param vidx named integer vector: variable name -> 0-based column in M
#' @return character vector of C++ lines
#' @keywords internal
cpp_emit_block_dense <- function(name, jac, formulas, vidx) {

  n <- jac$n
  sp <- cpp_split_constants(jac)

  cst_i <- jac$i[sp$is_const] - 1L
  cst_j <- jac$j[sp$is_const] - 1L
  var_i <- jac$i[sp$vary] - 1L
  var_j <- jac$j[sp$vary] - 1L

  jac_stmts <- vapply(seq_along(sp$vary), function(k) {
    paste0("v[", k - 1L, "]=", cpp_from_formula(jac$expr[sp$vary[k]], vidx), ";")
  }, character(1))

  res_stmts <- vapply(seq_along(formulas), function(k) {
    paste0("f[", k - 1L, "]=", cpp_from_formula(formulas[k], vidx), ";")
  }, character(1))

  endo_cols <- vidx[jac$endo]
  if (anyNA(endo_cols)) {
    stop("Block '", name, "': endogenous variables missing from the variable map.",
         call. = FALSE)
  }

  c(
    paste0("// ===================== block: ", name, " ====================="),
    sprintf("static const int %s_n   = %d;", name, n),
    sprintf("static const int %s_nvar = %d;", name, length(var_i)),
    sprintf("static const int %s_ncst = %d;", name, length(cst_i)),
    sprintf("static const int %s_var_i[] = {\n%s};", name, cpp_int_array(var_i)),
    sprintf("static const int %s_var_j[] = {\n%s};", name, cpp_int_array(var_j)),
    sprintf("static const int %s_cst_i[] = {\n%s};", name, cpp_int_array(cst_i)),
    sprintf("static const int %s_cst_j[] = {\n%s};", name, cpp_int_array(cst_j)),
    sprintf("static const double %s_cst_v[] = {\n%s};", name, cpp_dbl_array(sp$cst_v)),
    sprintf("static const int %s_endo[] = {\n%s};", name, cpp_int_array(as.integer(endo_cols))),
    "",
    cpp_chunked_function(jac_stmts, paste0(name, "_jac"),
                         "const MapMat& M, int t, double* v", "M,t,v"),
    cpp_chunked_function(res_stmts, paste0(name, "_res"),
                         "const MapMat& M, int t, double* f", "M,t,f"),
    ""
  )
}

#' Generate the dense C++ solver source for a model
#'
#' @param model_name name of the model
#' @param blocks named list of blocks; each element is
#'   `list(jac = <thor_sparse_jacobian>, formulas = <character>)`
#' @param all_model_vars character vector of every model variable, in the
#'   model's order (`model@vars$all`)
#' @param verbose print one line per block
#' @return character scalar: the full contents of the generated .cpp file
#' @keywords internal
generate_cpp_dense <- function(model_name, blocks, all_model_vars, verbose = TRUE) {

  vidx <- var_index(all_model_vars, base = 0L)

  active <- names(blocks)[vapply(blocks, function(b) !is.null(b) && b$jac$n > 0L,
                                 logical(1))]
  ## blocks solved one equation at a time rather than by Newton (sequential.R)
  seq_blocks <- active[vapply(blocks[active], function(b) !is.null(b$seq), logical(1))]

  lines <- c(
    paste0("// Generated by thortwo for model '", model_name, "' (dense backend). Do not edit."),
    paste0("// ", length(all_model_vars), " variables; blocks: ",
           paste(active, collapse = ", ")),
    CPP_PREAMBLE,
    CPP_DENSE_BLOCK,
    if (length(seq_blocks)) CPP_SEQ_RUNTIME,
    ""
  )

  for (nm in active) {
    b <- blocks[[nm]]
    if (verbose) {
      cat("   - block '", nm, "': ", b$jac$n, " equations, ",
          if (is.null(b$seq)) paste(length(b$jac$i), "jacobian entries")
          else "solved sequentially", "\n", sep = "")
    }
    lines <- c(lines, cpp_emit_block_dense(nm, b$jac, b$formulas, vidx))
    if (!is.null(b$seq)) lines <- c(lines, cpp_emit_sequential(nm, b$seq, vidx))
  }

  k <- stats::setNames(seq_along(active) - 1L, active)
  reg <- unlist(lapply(active, function(nm) c(
    sprintf('  B[%d].name="%s"; B[%d].n=%s_n; B[%d].nvar=%s_nvar; B[%d].ncst=%s_ncst;',
            k[[nm]], nm, k[[nm]], nm, k[[nm]], nm, k[[nm]], nm),
    sprintf('  B[%d].var_i=%s_var_i; B[%d].var_j=%s_var_j;',
            k[[nm]], nm, k[[nm]], nm),
    sprintf('  B[%d].cst_i=%s_cst_i; B[%d].cst_j=%s_cst_j; B[%d].cst_v=%s_cst_v;',
            k[[nm]], nm, k[[nm]], nm, k[[nm]], nm),
    sprintf('  B[%d].endo=%s_endo; B[%d].jac=&%s_jac; B[%d].res=&%s_res;',
            k[[nm]], nm, k[[nm]], nm, k[[nm]], nm)
  )))

  lines <- c(lines,
             "static void thor_register(Block* B){", reg, "}", "",
             cpp_entry_points(active, seq_blocks))

  paste0(paste(lines, collapse = "\n"), "\n")
}
