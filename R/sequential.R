## Sequential solution of the prologue and the epilogue.
##
## Both blocks are recursive by construction (decompose.R): in the right order
## every equation contains exactly one variable that is not yet known. By
## default they are nevertheless handed to Newton as a simultaneous system,
## like the heart. With `thor_model(sequential = TRUE)` they are solved the way
## their structure allows: one equation at a time, in dependency order, in a
## single pass.
##
## Each equation is one of two kinds.
##
##   direct   -- written `x = f(...)`, with x the variable the equation
##               determines and x absent from the right-hand side. The
##               generated code is the assignment itself. No derivative.
##   scalar   -- anything else (`log(y) = ...`, `pk * f_k = ...`,
##               `delta(1, log(x)) = ...`). Solved for its one unknown by a
##               one-variable Newton iteration, which needs one derivative:
##               the equation with respect to that variable.
##
## What that buys is in the build. No jacobian is derived or compiled for
## these blocks, only the scalar equations' own derivatives; on ThreeME 13x13
## the two blocks hold 3,768 of the 20,281 jacobian entries. And an equation
## like ThreeME's `verif_all`, a root of a sum of squares over hundreds of
## variables, costs nothing here: it is direct.
##
## The existing path is untouched. A sequential block is passed to the block
## generators with an empty jacobian, so its residual code is still emitted
## (model_residuals() and calibration_check() need it) and nothing about the
## Block type changes; the sequential code is emitted next to it, and the
## solve loop calls it instead of Newton for those blocks only. With
## `sequential = FALSE` the generated code is byte-for-byte what it was.

#' The order in which a recursive block can be solved
#'
#' Repeatedly takes an equation with exactly one undetermined variable of the
#' block, and marks that variable as determined.
#'
#' @param eqs character vector of the block's equation ids
#' @param endo character vector of the block's endogenous variables
#' @param eqns_vars_list named list: equation id -> contemporaneous variables
#' @return data.frame with `equation` and `variable`, in solve order; or NULL
#'   if the block is not recursive
#' @keywords internal
sequential_order <- function(eqs, endo, eqns_vars_list) {

  n <- length(eqs)
  if (n != length(endo)) return(NULL)

  vars <- lapply(eqs, function(id) {
    intersect(unique(as.character(eqns_vars_list[[id]])), endo)
  })
  ## variable -> the equations it occurs in
  where <- split(rep(seq_len(n), lengths(vars)), unlist(vars))

  open <- lengths(vars)                 # undetermined block variables per equation
  done <- logical(n)
  known <- new.env(hash = TRUE, parent = emptyenv())
  out_eq <- integer(n); out_var <- character(n); k <- 0L

  queue <- which(open == 1L)
  while (length(queue)) {
    e <- queue[1L]; queue <- queue[-1L]
    if (done[e]) next
    v <- vars[[e]]
    v <- v[!vapply(v, exists, logical(1), envir = known, inherits = FALSE)]
    if (length(v) != 1L) next           # its variable was claimed by another equation
    k <- k + 1L
    out_eq[k] <- e; out_var[k] <- v
    done[e] <- TRUE
    assign(v, TRUE, envir = known)
    for (e2 in where[[v]]) {
      if (done[e2]) next
      open[e2] <- open[e2] - 1L
      if (open[e2] == 1L) queue <- c(queue, e2)
    }
  }

  if (k < n) return(NULL)
  data.frame(equation = eqs[out_eq], variable = out_var, stringsAsFactors = FALSE)
}

#' Plan the sequential solution of one block
#'
#' @param block a block: `list(name, endo, equations, ...)`
#' @param equations_list the model's equations table, with `new_formula`
#' @param eqns_vars_list named list: equation id -> contemporaneous variables
#' @return NULL if the block is not recursive; otherwise a list with `order`
#'   (data.frame: `equation`, `name`, `variable`, `direct`), `rhs` (the
#'   right-hand sides, in formula syntax), `formula` (the residual formulas)
#'   and `deriv` (the scalar equations' derivatives, NA for direct ones)
#' @keywords internal
sequential_plan <- function(block, equations_list, eqns_vars_list) {

  ord <- sequential_order(block$equations, block$endo, eqns_vars_list)
  if (is.null(ord)) return(NULL)

  pos <- match(ord$equation, as.character(equations_list$id))
  lhs <- equations_list$LHS[pos]
  rhs <- formatting_formulas(equations_list$RHS[pos])
  formula <- equations_list$new_formula[pos]

  ## Direct: the left-hand side is the variable itself and the right-hand side
  ## does not use it in the current period. A lagged occurrence is a different
  ## symbol (`lag.x.1`) by then, so all.vars() only sees current-period ones.
  on_rhs <- mapply(function(v, r) v %in% all.vars(parse(text = r, keep.source = FALSE)),
                   ord$variable, rhs, USE.NAMES = FALSE)
  direct <- lhs == ord$variable & !on_rhs

  deriv <- rep(NA_character_, nrow(ord))
  for (k in which(!direct)) {
    d <- equation_derivatives(formula[k], ord$variable[k])
    if (identical(d, "0") || identical(d, "0L")) {
      stop("Equation '", equations_list$name[pos[k]], "' has to be solved for '",
           ord$variable[k], "', but its derivative with respect to that variable is zero.",
           call. = FALSE)
    }
    deriv[k] <- d
  }

  list(order = data.frame(equation = ord$equation, name = equations_list$name[pos],
                          variable = ord$variable, direct = direct,
                          stringsAsFactors = FALSE),
       rhs = rhs, formula = formula, deriv = deriv)
}

#' A jacobian with no entries, for a block that is solved sequentially
#'
#' @param block a block: `list(endo, equations, ...)`
#' @return a `thor_sparse_jacobian`
#' @keywords internal
empty_jacobian <- function(block) {
  structure(list(n = length(block$endo), equations = sort(block$equations),
                 endo = sort(block$endo), i = integer(0), j = integer(0),
                 expr = character(0)),
            class = "thor_sparse_jacobian")
}

## The part of the generated C++ that sequential blocks share. Emitted only
## when a model has one.
##
## thor_scalar is newton_block for a single equation and a single unknown: the
## same start value (whatever is in M, i.e. the previous period), the same
## backtracking, the same convergence measure.
CPP_SEQ_RUNTIME <- '
struct SeqCtl {
  double rtol, atol; int max_iter; bool damping;
  double conv;     // worst scaled step over the block, as reported by Newton
  int iters;       // most iterations any one equation needed
};
typedef double (*ThorFn)(const MapMat&, int);

static bool thor_scalar(ThorFn fr, ThorFn fd, MapMat& M, int t, int c, SeqCtl& S)
{
  double x = M(t, c);
  double r = fr(M, t);
  double s = std::numeric_limits<double>::infinity();

  for (int it = 0; it < S.max_iter; ++it) {
    if (!(r == r)) return false;
    double dx = r / fd(M, t);
    if (!std::isfinite(dx)) return false;

    double lambda = 1.0, rn = r;
    int tries = S.damping ? 12 : 1;
    bool ok = false;
    for (int b = 0; b < tries; ++b) {
      M(t, c) = x - lambda * dx;
      rn = fr(M, t);
      if (!S.damping || (rn == rn && std::fabs(rn) <= std::fabs(r))) { ok = true; break; }
      lambda *= 0.5;
    }

    s = std::fabs(lambda * dx) / (S.rtol * std::fabs(x) + S.atol);
    bool converged = s <= 1.0;
    if (!ok && !converged) { M(t, c) = x; return false; }
    if (ok) { x = M(t, c); r = rn; }
    if (converged) {
      if (s > S.conv) S.conv = s;
      if (it + 1 > S.iters) S.iters = it + 1;
      return true;
    }
  }
  return false;
}
'

#' Emit the C++ that solves one block sequentially
#'
#' @param name block name
#' @param seq a plan from [sequential_plan()]
#' @param vidx named integer vector: variable name -> 0-based column in M
#' @return character vector of C++ lines
#' @keywords internal
cpp_emit_sequential <- function(name, seq, vidx) {

  n <- nrow(seq$order)
  helpers <- character(0)
  stmts <- character(n)

  for (k in seq_len(n)) {
    col <- var_col(vidx, seq$order$variable[k])
    if (seq$order$direct[k]) {
      stmts[k] <- sprintf("M(t,%d)=%s;", col, cpp_from_formula(seq$rhs[k], vidx))
    } else {
      helpers <- c(helpers,
        sprintf("static double %s_r%d(const MapMat& M, int t){ return %s; }",
                name, k - 1L, cpp_from_formula(seq$formula[k], vidx)),
        sprintf("static double %s_d%d(const MapMat& M, int t){ return %s; }",
                name, k - 1L, cpp_from_formula(seq$deriv[k], vidx)))
      stmts[k] <- sprintf("if(!thor_scalar(&%s_r%d,&%s_d%d,M,t,%d,S)) return %d;",
                          name, k - 1L, name, k - 1L, col, k - 1L)
    }
  }

  ## Chunked like the residual and jacobian code, and for the same reason; a
  ## chunk returns the position of the equation that failed, or -1.
  idx <- split(seq_len(n), ceiling(seq_len(n) / CPP_CHUNK))
  chunks <- unlist(lapply(seq_along(idx), function(ci) c(
    sprintf("THOR_NOINLINE static int %s_seq_%d(MapMat& M, int t, SeqCtl& S){", name, ci - 1L),
    stmts[idx[[ci]]],
    "return -1;", "}", "")))

  quoted <- function(x) paste0('"', x, '"', collapse = ",")

  c(sprintf("// --------------------- sequential: %s ---------------------", name),
    helpers, "",
    chunks,
    sprintf("static int %s_seq(MapMat& M, int t, SeqCtl& S){", name),
    "  int k;",
    sprintf("  if((k=%s_seq_%d(M,t,S))>=0) return k;", name, seq_along(idx) - 1L),
    "  return -1;", "}",
    sprintf("static const char* %s_seq_eq[] = {\n%s};", name, quoted(seq$order$name)),
    sprintf("static const char* %s_seq_var[] = {\n%s};", name, quoted(seq$order$variable)),
    "")
}

#' Emit the R code that solves one block sequentially
#'
#' Produces `<block>_seq(t, M, rtol, atol, max_iter, damping)`, returning the
#' updated `M`, the position of the equation that failed (0 if none), the
#' worst scaled step and the most iterations any equation needed. The scalar
#' Newton is the same iteration as `thor_scalar` in the generated C++.
#'
#' It is written out in line for each scalar equation rather than called as a
#' function: a callee cannot update `M` in place, and returning a modified
#' copy of the whole data matrix once per equation would cost far more than
#' the solve.
#'
#' @param name block name
#' @param seq a plan from [sequential_plan()]
#' @param vidx named integer vector: variable name -> 1-based column in M
#' @return character vector of R lines
#' @keywords internal
r_emit_sequential <- function(name, seq, vidx) {

  n <- nrow(seq$order)
  helpers <- character(0)
  body <- character(0)

  for (k in seq_len(n)) {
    col <- var_col(vidx, seq$order$variable[k])
    at <- sprintf("M[t,%d]", col)
    if (seq$order$direct[k]) {
      body <- c(body, sprintf("  %s <- %s", at, r_from_formula(seq$rhs[k], vidx)))
      next
    }
    fr <- sprintf("%s_r%d", name, k); fd <- sprintf("%s_d%d", name, k)
    helpers <- c(helpers,
      sprintf("%s <- function(t, M) %s", fr, r_from_formula(seq$formula[k], vidx)),
      sprintf("%s <- function(t, M) %s", fd, r_from_formula(seq$deriv[k], vidx)))
    body <- c(body,
      sprintf("  .x <- %s; .r <- %s(t, M); .ok <- FALSE; .s <- Inf", at, fr),
      "  for (.i in seq_len(max_iter)) {",
      "    if (is.na(.r)) break",
      sprintf("    .dx <- .r / %s(t, M)", fd),
      "    if (!is.finite(.dx)) break",
      "    .l <- 1; .g <- FALSE",
      "    for (.b in seq_len(if (damping) 12L else 1L)) {",
      sprintf("      %s <- .x - .l * .dx; .rn <- %s(t, M)", at, fr),
      "      if (!damping || (!is.na(.rn) && abs(.rn) <= abs(.r))) { .g <- TRUE; break }",
      "      .l <- .l / 2",
      "    }",
      "    .s <- abs(.l * .dx) / (rtol * abs(.x) + atol)",
      sprintf("    if (!.g && .s > 1) { %s <- .x; break }", at),
      sprintf("    if (.g) { .x <- %s; .r <- .rn }", at),
      "    if (.s <= 1) { .ok <- TRUE; break }",
      "  }",
      sprintf("  if (!.ok) return(list(M = M, fail = %dL, conv = .s, iters = iters))", k),
      "  if (.s > conv) conv <- .s",
      "  if (.i > iters) iters <- .i")
  }

  c(helpers, "",
    sprintf("%s_seq <- function(t, M, rtol, atol, max_iter, damping) {", name),
    "  conv <- 0; iters <- 1L",
    body,
    "  list(M = M, fail = 0L, conv = conv, iters = iters)",
    "}", "")
}
