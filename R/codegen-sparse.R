## Generate the model-specific C++ solver, sparse from end to end.
##
## What this replaces in tresthor's dense generator:
##   * `mat Jacobian_n(n, n, fill::zeros)` allocated on every Newton iteration
##     of every period, then scanned element-by-element to build a sparse copy.
##   * the constant entries of the jacobian passed as a runtime-parsed
##     `arma::uvec("0 5 19 ...")` string literal.
##   * `spsolve()`, which redoes the fill-reducing ordering and the symbolic
##     analysis on every call.
##
## The generated code instead holds one Eigen sparse matrix per block whose
## pattern is fixed at build time. Each Newton iteration only overwrites the
## numerical values of the entries that actually vary, and the fill-reducing
## ordering is computed once per block for the whole simulation.

## The sparse block type. `newton_block` in CPP_PREAMBLE is templated over it.
CPP_SPARSE_BLOCK <- '
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

typedef Eigen::SparseMatrix<double> SpMat;

struct Block {
  const char* name;
  int n, nnz, ncst;
  const int* Ap; const int* Ai;
  const int* cst_k; const double* cst_v;
  const int* endo;
  void (*jac)(const MapMat&, int, double*);
  void (*res)(const MapMat&, int, double*);

  SpMat A;
  Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> > lu;

  Block() {}

  // Build the sparsity pattern once, plant the constant entries, and run the
  // fill-reducing ordering + symbolic factorisation a single time.
  void setup(){
    std::vector<double> z(nnz, 0.0);
    A = Eigen::Map<const SpMat>(n, n, nnz, Ap, Ai, z.data());
    if (A.nonZeros() != nnz)
      Rcpp::stop(std::string("block ") + name + ": sparsity pattern was not preserved.");
    double* v = A.valuePtr();
    for (int c = 0; c < ncst; ++c) v[cst_k[c]] = cst_v[c];
    // NB: Eigen only initialises SparseLU::m_info in factorize(), so info()
    // must not be consulted after analyzePattern() alone -- it would read
    // uninitialised memory. Failures surface on the first factorize().
    lu.analyzePattern(A);
  }

  void refresh(const MapMat& M, int t){ jac(M, t, A.valuePtr()); }
  bool factorize(){ lu.factorize(A); return lu.info() == Eigen::Success; }
  bool solve(const Eigen::VectorXd& f, Eigen::VectorXd& dx){
    dx = lu.solve(f);
    // SparseLU reports success on a jacobian holding NaN or inf; the step
    // then carries them, so that is what gets checked, as in the dense block.
    return lu.info() == Eigen::Success && dx.allFinite();
  }
};
'

#' Emit the C++ for one sparse block
#'
#' @param name block name ("prologue", "heart", "epilogue")
#' @param jac a `thor_sparse_jacobian`
#' @param formulas character vector of the block's residual formulas, in the
#'   same order as `jac$equations`
#' @param vidx named integer vector: variable name -> 0-based column in M
#' @return character vector of C++ lines
#' @keywords internal
cpp_emit_block_sparse <- function(name, jac, formulas, vidx) {

  n   <- jac$n
  nnz <- length(jac$i)

  ## --- CSC pattern ---------------------------------------------------------
  ## jac$i / jac$j are already sorted by (column, row).
  Ai <- jac$i - 1L                                   # 0-based row indices
  Ap <- c(0L, cumsum(tabulate(jac$j, nbins = n)))    # column pointers

  sp <- cpp_split_constants(jac)
  cst_k <- which(sp$is_const) - 1L                   # 0-based slot in valuePtr

  jac_stmts <- vapply(sp$vary, function(k) {
    paste0("v[", k - 1L, "]=", cpp_from_formula(jac$expr[k], vidx), ";")
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
    sprintf("static const int %s_nnz = %d;", name, nnz),
    sprintf("static const int %s_ncst = %d;", name, length(cst_k)),
    sprintf("static const int %s_Ap[] = {\n%s};", name, cpp_int_array(Ap)),
    sprintf("static const int %s_Ai[] = {\n%s};", name, cpp_int_array(Ai)),
    sprintf("static const int %s_cst_k[] = {\n%s};", name, cpp_int_array(cst_k)),
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

#' Generate the sparse C++ solver source for a model
#'
#' @param model_name name of the model
#' @param blocks named list of blocks; each element is
#'   `list(jac = <thor_sparse_jacobian>, formulas = <character>)`, named
#'   "prologue" / "heart" / "epilogue". Blocks that are absent are skipped.
#' @param all_model_vars character vector of every model variable, in the
#'   model's order (`model@vars$all`)
#' @param verbose print one line per block
#' @return character scalar: the full contents of the generated .cpp file
#' @keywords internal
generate_cpp_sparse <- function(model_name, blocks, all_model_vars, verbose = TRUE) {

  vidx <- var_index(all_model_vars, base = 0L)

  active <- names(blocks)[vapply(blocks, function(b) !is.null(b) && b$jac$n > 0L,
                                 logical(1))]
  ## blocks solved one equation at a time rather than by Newton (sequential.R)
  seq_blocks <- active[vapply(blocks[active], function(b) !is.null(b$seq), logical(1))]

  lines <- c(
    paste0("// Generated by thortwo for model '", model_name, "' (sparse backend). Do not edit."),
    paste0("// ", length(all_model_vars), " variables; blocks: ",
           paste(active, collapse = ", ")),
    CPP_PREAMBLE,
    CPP_SPARSE_BLOCK,
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
    lines <- c(lines, cpp_emit_block_sparse(nm, b$jac, b$formulas, vidx))
    if (!is.null(b$seq)) lines <- c(lines, cpp_emit_sequential(nm, b$seq, vidx))
  }

  k <- stats::setNames(seq_along(active) - 1L, active)
  reg <- unlist(lapply(active, function(nm) c(
    sprintf('  B[%d].name="%s"; B[%d].n=%s_n; B[%d].nnz=%s_nnz; B[%d].ncst=%s_ncst;',
            k[[nm]], nm, k[[nm]], nm, k[[nm]], nm, k[[nm]], nm),
    sprintf('  B[%d].Ap=%s_Ap; B[%d].Ai=%s_Ai; B[%d].cst_k=%s_cst_k; B[%d].cst_v=%s_cst_v;',
            k[[nm]], nm, k[[nm]], nm, k[[nm]], nm, k[[nm]], nm),
    sprintf('  B[%d].endo=%s_endo; B[%d].jac=&%s_jac; B[%d].res=&%s_res;',
            k[[nm]], nm, k[[nm]], nm, k[[nm]], nm)
  )))

  lines <- c(lines,
             "static void thor_register(Block* B){", reg, "}", "",
             cpp_entry_points(active, seq_blocks))

  paste0(paste(lines, collapse = "\n"), "\n")
}
