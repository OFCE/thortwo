## The part of the generated C++ that does not depend on the model, and is
## shared by the sparse and the dense-cpp backends.
##
## `newton_block` is a template over the block type, so the two backends
## differ only in how a block stores, refreshes, factorises and solves its
## jacobian -- not in the Newton iteration itself. That is what makes the two
## backends agree to machine precision rather than merely approximately.

## Statements per generated function. Compilers are superlinear in the size of
## a single function body, so large blocks are split across several.
##
## The chunks must also be kept from being inlined: each is static and called
## exactly once, so at -O2 clang inlines them straight back into the caller
## and the split is undone -- smaller chunks then made things *worse* (150
## statements: 156 s). With `noinline` they stay separate and small is fast.
## Measured on ThreeME 13x13 (2 MB of generated C++, -O2):
##
##   statements per chunk     1500    500    150     50     20
##   as before (inlined)      83 s    67 s  156 s
##   noinline                 76 s    62 s   38 s   20 s   15 s
##
## The fixed part (Rcpp, Eigen, the Newton template) is ~4 s of that.
CPP_CHUNK <- 25L

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
             paste0("THOR_NOINLINE static void ", fname, "_", c_i - 1L, "(", sig, "){"),
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

// Keeps the generated chunks out of their caller: inlined, they would be
// compiled as one giant function again (see CPP_CHUNK in codegen-runtime.R).
#if defined(__GNUC__) || defined(__clang__)
#define THOR_NOINLINE __attribute__((noinline))
#else
#define THOR_NOINLINE
#endif

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
//
// A NaN anywhere makes the result NaN. `if (d > c)` alone would skip it, since
// every comparison with NaN is false, and a step that is not a number would
// then be measured by its other components and could pass as converged.
static double step_measure(const Eigen::VectorXd& dx, const Eigen::VectorXd& x,
                           double rtol, double atol)
{
  double c = 0.0;
  for (int i = 0; i < dx.size(); ++i) {
    double d = std::fabs(dx[i]) / (rtol * std::fabs(x[i]) + atol);
    if (!(d == d)) return d;
    if (d > c) c = d;
  }
  return c;
}

// Largest absolute value of n numbers, NaN if any of them is NaN.
//
// This is what Newton tests to see whether a residual has gone bad. It used
// to be Eigen`s lpNorm<Infinity>(), which skips a NaN unless it is the first
// element: (1, NaN, 3) gave 3. A block whose residual was NaN in any equation
// but its first was therefore not seen to have failed, and the solver could
// return NaN for some variables while reporting convergence.
static double max_abs(const double* f, int n)
{
  double m = 0.0;
  for (int i = 0; i < n; ++i) {
    double a = std::fabs(f[i]);
    if (!(a == a)) return a;
    if (a > m) m = a;
  }
  return m;
}

// Newton on one block at one period. Returns the number of iterations used,
// or -1 if it did not converge.
//
// Templated over the block type so that the sparse and the dense backends run
// the identical iteration, and differ only in refresh/factorize/solve.
//
// `reuse` keeps the factorised jacobian from one iteration to the next within
// the period, instead of recomputing and refactorising it every time. On a
// large block the factorisation is nearly the whole cost of an iteration: on
// the heart of ThreeME 29x33 (18,380 equations) it is 1.36 s, against 5 ms to
// use the factors and a few ms to evaluate the equations, because the factors
// hold 93 times more entries than the jacobian. A step taken with factors
// that are slightly out of date is a little less accurate, so more steps are
// needed, but each costs next to nothing.
//
// It is safeguarded so that the worst case is the ordinary iteration:
//   * old factors are tried first, with a full step only. If that step does
//     not reduce the residual it is discarded, and the iteration is redone
//     with a fresh jacobian, exactly as without `reuse`;
//   * factors are kept for the next iteration only while they at least halve
//     the residual.
//
// One more thing differs: when to stop. A step is accepted as the last one
// when it is smaller than the tolerance. A fresh Newton step converges
// quadratically, so the first step under the tolerance leaves an error far
// below it. A step with old factors converges only linearly and would stop
// just under the tolerance: on Opale that left residuals of 3e-6 where Newton
// leaves 2e-10. So a step with old factors has to be REUSE_TIGHT times
// smaller than the tolerance to count as the last one; such steps are cheap,
// and the few extra ones buy back the accuracy.
//
// Every period starts with a fresh factorisation. Carrying the factors over
// from the previous period was tried and dropped: on ThreeME 29x33 it was no
// faster (76 s against 74 s), because the old factors rarely halve the
// residual from a new period`s starting point, and it left residuals ten
// times larger.
//
// With `reuse` false every iteration is fresh and this is the iteration it
// always was.
static const double REUSE_TIGHT = 1e-3;

template <class BlockT>
static int newton_block(BlockT& B, MapMat& M, int t,
                        double rtol, double atol, int max_iter, bool damping, bool reuse,
                        double& final_resid, double& final_conv)
{
  const int n = B.n;
  Eigen::VectorXd x(n), f(n), f0(n), dx(n), xtry(n);

  for (int i = 0; i < n; ++i) x[i] = M(t, B.endo[i]);

  B.res(M, t, f.data());
  double rnorm = max_abs(f.data(), n);
  double conv = std::numeric_limits<double>::infinity();
  bool stale = false;       // factors from an earlier iteration are worth trying

  for (int it = 0; it < max_iter; ++it) {

    if (!(rnorm == rnorm))                      // NaN in the residual
      { final_resid = rnorm; final_conv = conv; return -1; }

    bool fresh = !stale;
    double lambda = 1.0;
    double new_rnorm = rnorm;
    bool ok = false;

    for (;;) {                                  // old factors first, then fresh ones
      if (fresh) {
        B.refresh(M, t);
        if (!B.factorize()) { final_resid = rnorm; final_conv = conv; return -1; }
      }
      if (!B.solve(f, dx)) {
        if (!fresh) { fresh = true; continue; }
        final_resid = rnorm; final_conv = conv; return -1;
      }

      // Full Newton step, backtracked only if it makes the residual worse.
      f0 = f;
      lambda = 1.0; new_rnorm = rnorm; ok = false;
      int tries = (fresh && damping) ? 12 : 1;
      for (int b = 0; b < tries; ++b) {
        xtry = x - lambda * dx;
        for (int i = 0; i < n; ++i) M(t, B.endo[i]) = xtry[i];
        B.res(M, t, f.data());
        new_rnorm = max_abs(f.data(), n);
        if ((fresh && !damping) || (new_rnorm == new_rnorm && new_rnorm <= rnorm)) { ok = true; break; }
        lambda *= 0.5;
      }

      conv = step_measure(lambda * dx, x, rtol, atol);

      if (ok || fresh || conv <= 1.0) break;
      // The old factors gave a step that does not help: take it back.
      for (int i = 0; i < n; ++i) M(t, B.endo[i]) = x[i];
      f = f0;
      fresh = true;
    }

    if (!ok) {
      // No descent direction. Near the solution this just means the residual
      // has bottomed out in floating point, which is a success, not a failure.
      if (conv <= 1.0) { final_resid = rnorm; final_conv = conv; return it + 1; }
      for (int i = 0; i < n; ++i) M(t, B.endo[i]) = x[i];
      final_resid = rnorm; final_conv = conv;
      return -1;
    }

    stale = reuse && (new_rnorm <= 0.5 * rnorm);
    x = xtry;
    rnorm = new_rnorm;

    if (conv <= (fresh ? 1.0 : REUSE_TIGHT)) { final_resid = rnorm; final_conv = conv; return it + 1; }
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
#' @param seq character vector of the blocks solved sequentially (see
#'   sequential.R). When empty, the output is exactly what it was before
#'   sequential blocks existed.
#' @return character vector of C++ lines
#' @keywords internal
cpp_entry_points <- function(active, seq = character(0)) {

  has_seq <- length(seq) > 0L
  per_block <- function(fmt) paste(ifelse(active %in% seq, sprintf(fmt, active), "0"),
                                   collapse = ", ")
  seq_tables <- if (has_seq) c(
    "// Which blocks are solved sequentially, and how (0: by Newton).",
    sprintf("static int (*const THOR_SEQ[NB])(MapMat&, int, SeqCtl&) = {%s};", per_block("&%s_seq")),
    sprintf("static const char* const* const THOR_SEQ_EQ[NB] = {%s};", per_block("%s_seq_eq")),
    sprintf("static const char* const* const THOR_SEQ_VAR[NB] = {%s};", per_block("%s_seq_var")),
    "")

  setup_line <- if (has_seq) "  for (int b = 0; b < NB; ++b) if (!THOR_SEQ[b]) B[b].setup();"
                else         "  for (int b = 0; b < NB; ++b) B[b].setup();"

  newton_call <- "      int it = newton_block(B[b], M, t, rtol, atol, max_iter, damping, reuse, r, c);"
  solve_call <- if (!has_seq) newton_call else c(
    "      int it;",
    "      if (THOR_SEQ[b]) {",
    "        SeqCtl S; S.rtol = rtol; S.atol = atol; S.max_iter = max_iter;",
    "        S.damping = damping; S.conv = 0.0; S.iters = 1;",
    "        int bad = THOR_SEQ[b](M, t, S);",
    "        if (bad >= 0) {",
    "          char buf[640];",
    "          std::snprintf(buf, sizeof(buf),",
    "            \"Sequential solve failed on block '%s' at row %d: equation '%s' could \"",
    "            \"not be solved for '%s'. Its derivative with respect to that variable is \"",
    "            \"zero or not finite, or it did not converge in %d iterations. Check the \"",
    "            \"data at that period.\",",
    "            B[b].name, t + 1, THOR_SEQ_EQ[b][bad], THOR_SEQ_VAR[b][bad], max_iter);",
    "          Rcpp::stop(buf);",
    "        }",
    "        // the residuals of the block, as Newton would report them",
    "        std::vector<double> f(B[b].n);",
    "        B[b].res(M, t, f.data());",
    "        r = max_abs(f.data(), B[b].n);",
    "        if (!(r == r)) {",
    "          char buf[512];",
    "          std::snprintf(buf, sizeof(buf),",
    "            \"Sequential solve of block '%s' at row %d produced a value that is not \"",
    "            \"a number. Check the data at that period with calibration_check().\",",
    "            B[b].name, t + 1);",
    "          Rcpp::stop(buf);",
    "        }",
    "        c = S.conv; it = S.iters;",
    "      } else {",
    "        it = newton_block(B[b], M, t, rtol, atol, max_iter, damping, reuse, r, c);",
    "      }")

  c(
    sprintf("static const int NB = %d;", length(active)),
    sprintf("static const char* BLOCK_NAMES[%d] = {%s};", length(active),
            paste0('"', active, '"', collapse = ",")),
    "",
    seq_tables,
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
    "    out[b] = max_abs(f.data(), B[b].n);",
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
    "                          bool damping, bool reuse, bool verbose)",
    "{",
    "  MapMat M(data.begin(), data.nrow(), data.ncol());",
    "  // Block holds an Eigen factorisation object, which is not copyable, so",
    "  // this has to be a fixed array rather than a std::vector.",
    "  Block B[NB];",
    "  thor_register(B);",
    setup_line,
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
    solve_call,
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
