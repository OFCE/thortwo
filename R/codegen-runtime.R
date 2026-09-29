## The part of the generated C++ that does not depend on the model, and is
## shared by the sparse and the dense-cpp backends.
##
## `newton_block` is a template over the block type, so the two backends
## differ only in how a block stores, refreshes, factorises and solves its
## jacobian -- not in the Newton iteration itself. That is what makes the two
## backends agree to machine precision rather than merely approximately.

## Statements per generated function. Compilers are superlinear in the size of
## a single function body, so large blocks are split across several.
CPP_CHUNK <- 1500L

#' Format an integer vector as a C array initialiser
#' @keywords internal
cpp_int_array <- function(x, per_line = 20L) {
  if (length(x) == 0L) return("0")
  grp <- split(x, ceiling(seq_along(x) / per_line))
  paste(vapply(grp, function(g) paste(g, collapse = ","), character(1)),
        collapse = ",\n")
}

#' Format a double vector as a C array initialiser
#' @keywords internal
cpp_dbl_array <- function(x, per_line = 12L) {
  if (length(x) == 0L) return("0.0")
  s <- format(x, digits = 17L, scientific = FALSE, trim = TRUE)
  grp <- split(s, ceiling(seq_along(s) / per_line))
  paste(vapply(grp, function(g) paste(g, collapse = ","), character(1)),
        collapse = ",\n")
}

#' Split a vector of C++ statements into chunked functions
#'
#' @param stmts character vector of statement lines
#' @param fname base name of the generated function
#' @param sig argument list of the generated function
#' @param callargs arguments forwarded to each chunk
#' @return character vector of C++ lines defining `fname` and its chunks
#' @keywords internal
cpp_chunked_function <- function(stmts, fname, sig, callargs) {
  if (length(stmts) == 0L) {
    return(c(paste0("static void ", fname, "(", sig, "){ (void)M; (void)t; }"), ""))
  }
  idx <- split(seq_along(stmts), ceiling(seq_along(stmts) / CPP_CHUNK))
  out <- character(0)
  for (c_i in seq_along(idx)) {
    out <- c(out,
             paste0("static void ", fname, "_", c_i - 1L, "(", sig, "){"),
             stmts[idx[[c_i]]],
             "}", "")
  }
  c(out,
    paste0("static void ", fname, "(", sig, "){"),
    paste0("  ", fname, "_", seq_along(idx) - 1L, "(", callargs, ");"),
    "}", "")
}

#' Split a block's jacobian entries into constant and varying
#'
#' The constant entries are planted once when the block is set up; only the
#' varying ones are recomputed on each Newton iteration. On Opale that is
#' roughly half the entries.
#'
#' @param jac a `thor_sparse_jacobian`
#' @return list with `is_const`, `cst_v` (values) and `vary` (indices into
#'   `jac$expr`)
#' @keywords internal
cpp_split_constants <- function(jac) {
  num <- suppressWarnings(as.numeric(jac$expr))
  is_const <- !is.na(num)
  list(is_const = is_const, cst_v = num[is_const], vary = which(!is_const))
}

## Includes, the convergence measure and the Newton iteration.
CPP_PREAMBLE <- '
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]
//
// C++17 is pinned deliberately. R >= 4.5 compiles Rcpp sources as gnu++20 by
// default, and on clang that makes this file take ~25x longer to compile
// (394 s versus 16 s for ThreeME 4x4) for no benefit here.
//
// Only the Eigen headers actually needed are included: pulling in the whole
// RcppEigen.h umbrella drags in the unsupported modules and costs seconds per
// translation unit.
#include <Rcpp.h>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include <cstdio>

typedef Eigen::Map<Eigen::MatrixXd> MapMat;

// Scale-free convergence measure. Macro models mix variables of very
// different magnitude -- in ThreeME 4x4 the heart block spans 0 to 1.2e7, and
// 338 of its variables are identically zero in a given scenario -- so neither
// a purely absolute nor a purely relative step test works.
//
// This is the standard mixed criterion: a step is acceptable when
//     |dx_i| <= rtol*|x_i| + atol
// which is returned in normalised form, so <= 1 means converged. Writing it
// instead as |dx_i|/(atol+|x_i|) <= rtol would be 1/rtol times stricter on
// variables sitting at zero, where it would measure nothing but rounding
// noise (dx ~ 5e-16 against a 1e-8 floor reads as 5e-8).
static double step_measure(const Eigen::VectorXd& dx, const Eigen::VectorXd& x,
                           double rtol, double atol)
{
  double c = 0.0;
  for (int i = 0; i < dx.size(); ++i) {
    double d = std::fabs(dx[i]) / (rtol * std::fabs(x[i]) + atol);
    if (d > c) c = d;
  }
  return c;
}

// Newton on one block at one period. Returns the number of iterations used,
// or -1 if it did not converge.
//
// Templated over the block type so that the sparse and the dense backends run
// the identical iteration, and differ only in refresh/factorize/solve.
template <class BlockT>
static int newton_block(BlockT& B, MapMat& M, int t,
                        double rtol, double atol, int max_iter, bool damping,
                        double& final_resid, double& final_conv)
{
  const int n = B.n;
  Eigen::VectorXd x(n), f(n), dx(n), xtry(n);

  for (int i = 0; i < n; ++i) x[i] = M(t, B.endo[i]);

  B.res(M, t, f.data());
  double rnorm = f.lpNorm<Eigen::Infinity>();
  double conv = std::numeric_limits<double>::infinity();

  for (int it = 0; it < max_iter; ++it) {

    if (!(rnorm == rnorm))                      // NaN in the residual
      { final_resid = rnorm; final_conv = conv; return -1; }

    B.refresh(M, t);
    if (!B.factorize()) { final_resid = rnorm; final_conv = conv; return -1; }
    if (!B.solve(f, dx)) { final_resid = rnorm; final_conv = conv; return -1; }

    // Full Newton step, backtracked only if it makes the residual worse.
    double lambda = 1.0;
    double new_rnorm = rnorm;
    int tries = damping ? 12 : 1;
    bool ok = false;
    for (int b = 0; b < tries; ++b) {
      xtry = x - lambda * dx;
      for (int i = 0; i < n; ++i) M(t, B.endo[i]) = xtry[i];
      B.res(M, t, f.data());
      new_rnorm = f.lpNorm<Eigen::Infinity>();
      if (!damping || (new_rnorm == new_rnorm && new_rnorm <= rnorm)) { ok = true; break; }
      lambda *= 0.5;
    }

    conv = step_measure(lambda * dx, x, rtol, atol);

    if (!ok) {
      // No descent direction. Near the solution this just means the residual
      // has bottomed out in floating point, which is a success, not a failure.
      if (conv <= 1.0) { final_resid = rnorm; final_conv = conv; return it + 1; }
      for (int i = 0; i < n; ++i) M(t, B.endo[i]) = x[i];
      final_resid = rnorm; final_conv = conv;
      return -1;
    }

    x = xtry;
    rnorm = new_rnorm;

    if (conv <= 1.0) { final_resid = rnorm; final_conv = conv; return it + 1; }
  }

  final_resid = rnorm; final_conv = conv;
  return -1;
}
'

#' The exported entry points of a generated file
#'
#' The sparse and dense generators emit the same two functions with the same
#' signatures, so [compile_model_cpp()], [thor_solve()] and
#' [model_residuals()] do not need to know which backend they are holding.
#'
#' @param active character vector of the block names, in solve order
#' @return character vector of C++ lines
#' @keywords internal
cpp_entry_points <- function(active) {
  c(
    sprintf("static const int NB = %d;", length(active)),
    sprintf("static const char* BLOCK_NAMES[%d] = {%s};", length(active),
            paste0('"', active, '"', collapse = ",")),
    "",
    "// Maximum absolute residual of each block at one observation. Lets the",
    "// caller check that a solution really does satisfy the equations, and",
    "// compare backends on the same footing.",
    "// [[Rcpp::export]]",
    "Rcpp::NumericVector thor_cpp_residuals(Rcpp::NumericMatrix data, int row)",
    "{",
    "  MapMat M(data.begin(), data.nrow(), data.ncol());",
    "  Block B[NB];",
    "  thor_register(B);",
    "  Rcpp::NumericVector out(NB);",
    "  for (int b = 0; b < NB; ++b) {",
    "    std::vector<double> f(B[b].n);",
    "    B[b].res(M, row, f.data());",
    "    double mx = 0.0;",
    "    for (int i = 0; i < B[b].n; ++i) { double a = std::fabs(f[i]); if (a > mx) mx = a; }",
    "    out[b] = mx;",
    "  }",
    "  out.attr(\"names\") = Rcpp::CharacterVector(BLOCK_NAMES, BLOCK_NAMES + NB);",
    "  return out;",
    "}",
    "",
    "// Every equation's residual at one observation, blocks in solve order and",
    "// equations within a block in the order the model object lists them. Lets",
    "// the caller see *which* equation is off, not just that one is.",
    "// [[Rcpp::export]]",
    "Rcpp::NumericVector thor_cpp_residuals_all(Rcpp::NumericMatrix data, int row)",
    "{",
    "  MapMat M(data.begin(), data.nrow(), data.ncol());",
    "  Block B[NB];",
    "  thor_register(B);",
    "  int total = 0;",
    "  for (int b = 0; b < NB; ++b) total += B[b].n;",
    "  Rcpp::NumericVector out(total);",
    "  int k = 0;",
    "  for (int b = 0; b < NB; ++b) {",
    "    std::vector<double> f(B[b].n);",
    "    B[b].res(M, row, f.data());",
    "    for (int i = 0; i < B[b].n; ++i) out[k++] = f[i];",
    "  }",
    "  return out;",
    "}",
    "",
    "// [[Rcpp::export]]",
    "Rcpp::List thor_cpp_solve(Rcpp::NumericMatrix data,",
    "                          int first_date, int last_date,",
    "                          double rtol, double atol, int max_iter,",
    "                          bool damping, bool verbose)",
    "{",
    "  MapMat M(data.begin(), data.nrow(), data.ncol());",
    "  // Block holds an Eigen factorisation object, which is not copyable, so",
    "  // this has to be a fixed array rather than a std::vector.",
    "  Block B[NB];",
    "  thor_register(B);",
    "  for (int b = 0; b < NB; ++b) B[b].setup();",
    "",
    "  Rcpp::IntegerVector iters(last_date - first_date + 1);",
    "  Rcpp::NumericVector resid(last_date - first_date + 1);",
    "  Rcpp::NumericVector conv(last_date - first_date + 1);",
    "",
    "  for (int t = first_date; t <= last_date; ++t) {",
    "    // carry the previous period forward as the starting point",
    "    for (int b = 0; b < NB; ++b)",
    "      for (int i = 0; i < B[b].n; ++i)",
    "        M(t, B[b].endo[i]) = M(t - 1, B[b].endo[i]);",
    "",
    "    int tot = 0; double worst = 0.0; double worstc = 0.0;",
    "    for (int b = 0; b < NB; ++b) {",
    "      double r = 0.0, c = 0.0;",
    "      int it = newton_block(B[b], M, t, rtol, atol, max_iter, damping, r, c);",
    "      if (it < 0) {",
    "        char buf[512];",
    "        std::snprintf(buf, sizeof(buf),",
    "          \"Newton did not converge on block '%s' at row %d after %d iterations: \"",
    "          \"scaled step %.3e (converges at 1.0, rtol %.3e), max |residual| %.3e. \"",
    "          \"Try a looser rtol, a higher max_iter, or check the data at that period.\",",
    "          B[b].name, t + 1, max_iter, c, rtol, r);",
    "        Rcpp::stop(buf);",
    "      }",
    "      tot += it; if (r > worst) worst = r; if (c > worstc) worstc = c;",
    "    }",
    "    iters[t - first_date] = tot;",
    "    resid[t - first_date] = worst;",
    "    conv[t - first_date] = worstc;",
    "    if (verbose) Rcpp::Rcout << \"  \" << (t + 1) << \" (\" << tot << \" it) \";",
    "  }",
    "  if (verbose) Rcpp::Rcout << std::endl;",
    "",
    "  return Rcpp::List::create(Rcpp::_[\"data\"] = data,",
    "                            Rcpp::_[\"iterations\"] = iters,",
    "                            Rcpp::_[\"residuals\"] = resid,",
    "                            Rcpp::_[\"convergence\"] = conv);",
    "}"
  )
}

#' Write generated code to a file, but only if it has changed
#'
#' Rcpp's compile cache is invalidated by the source file's modification time,
#' not only by its contents, so rewriting an identical file would force a full
#' recompile. Code generation is deterministic, so a model whose equations have
#' not changed keeps its file, and its cached object, untouched.
#'
#' @param code character scalar, the full file contents
#' @param path where to write it
#' @return `path`, with an `unchanged` attribute
#' @keywords internal
write_if_changed <- function(code, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  unchanged <- file.exists(path) &&
    identical(readChar(path, file.size(path), useBytes = TRUE), code)
  if (!unchanged) writeChar(code, path, eos = NULL)
  attr(path, "unchanged") <- unchanged
  path
}
