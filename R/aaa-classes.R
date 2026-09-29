## The model object.
##
## tresthor's `thoR.model` carried 24 slots, six of which held `print("none")`
## stubs whenever the backend was not the R one, three of which held an empty
## `matrix()` whenever the backend was sparse, and none of which could hold a
## sparse jacobian -- so the sparse path smuggled its real jacobians through
## `attr(model, "sparse_jacobians")`. Here there is one class with a
## backend-specific payload in one slot, and nothing is smuggled.

#' Functions a model formula may use
#'
#' The code generators translate exactly these, plus the arithmetic operators,
#' `lag()`/`mylg()` and `delta()`/`newdiff()`.
#'
#' @export
thor_functions_supported <- c(
  "abs", "acos", "acosh", "asin", "asinh", "atan", "atanh", "cos", "cosh",
  "exp", "expm1", "log", "log10", "log1p", "log2", "logb", "sign", "sin",
  "sinh", "sqrt", "tan", "tanh"
)

#' The backends a model can be built with
#'
#' * `sparse`    -- triplet jacobians, generated Eigen `SparseLU` C++.
#' * `dense-cpp` -- dense jacobians, generated Eigen `PartialPivLU` C++.
#' * `dense-r`   -- dense jacobians, generated R closures. No compiler needed.
#'
#' @export
thor_backends <- c("sparse", "dense-cpp", "dense-r")

#' A thortwo model
#'
#' Built by [thor_model()], solved by [thor_solve()].
#'
#' @slot name character. Name of the model.
#' @slot backend character. One of [thor_backends].
#' @slot equations data.frame. One row per equation: `id`, `name`, `equation`,
#'   `LHS`, `RHS`, `formula`, `new_formula`, `part`.
#' @slot vars list. `endo`, `exo`, `coeff`: lower-case, sorted character
#'   vectors, and `all`, their sorted union. `all` is the column ordering of
#'   the data matrix handed to every backend.
#' @slot blocks list. One element per block (`prologue`, `heart`, `epilogue`),
#'   each `list(name, present, endo, equations)`. Absent blocks are dropped.
#' @slot jacobian list. Backend-specific, one element per present block:
#'   a `thor_sparse_jacobian` for `sparse`, a character matrix otherwise.
#' @slot generated list. `path`, `code`, `hash`, and for `dense-r` the
#'   generated R source. `code` is authoritative; `path` is a hint.
#' @slot meta list. `version`, `built_at`, `source`, `source_hash`, `decompose`.
#'
#' @export
setClass("thor_model", slots = c(
  name      = "character",
  backend   = "character",
  equations = "data.frame",
  vars      = "list",
  blocks    = "list",
  jacobian  = "list",
  generated = "list",
  meta      = "list"
))

## A malformed model used to surface as an error three steps later, usually
## inside generated code. Catch it at construction instead.
setValidity("thor_model", function(object) {
  problems <- character(0)

  if (length(object@name) != 1L || !nzchar(object@name)) {
    problems <- c(problems, "name must be a single non-empty string")
  }
  if (length(object@backend) != 1L || !object@backend %in% thor_backends) {
    problems <- c(problems, paste0("backend must be one of ",
                                   paste(thor_backends, collapse = ", ")))
  }

  needed <- c("id", "name", "equation", "formula", "new_formula", "part")
  missing_cols <- setdiff(needed, names(object@equations))
  if (length(missing_cols)) {
    problems <- c(problems, paste0("equations is missing the column(s) ",
                                   paste(missing_cols, collapse = ", ")))
  }

  v <- object@vars
  if (!all(c("endo", "exo", "coeff", "all") %in% names(v))) {
    problems <- c(problems, "vars must have the elements endo, exo, coeff, all")
  } else {
    if (!identical(v$all, sort(unique(c(v$endo, v$exo, v$coeff))))) {
      problems <- c(problems, "vars$all must be the sorted union of endo, exo and coeff")
    }
    for (nm in c("endo", "exo", "coeff")) {
      if (length(v[[nm]]) && !identical(v[[nm]], tolower(v[[nm]]))) {
        problems <- c(problems, paste0("vars$", nm, " must be lower-case"))
      }
    }
    if (length(intersect(v$endo, v$exo)) || length(intersect(v$endo, v$coeff)) ||
        length(intersect(v$exo, v$coeff))) {
      problems <- c(problems, "endo, exo and coeff must not overlap")
    }
  }

  if (length(object@blocks)) {
    ## Every block is square, and together the blocks partition the
    ## endogenous variables and the equations exactly once.
    for (b in object@blocks) {
      if (length(b$endo) != length(b$equations)) {
        problems <- c(problems, sprintf(
          "block '%s' is not square: %d endogenous variables for %d equations",
          b$name, length(b$endo), length(b$equations)))
      }
    }
    all_endo <- unlist(lapply(object@blocks, `[[`, "endo"), use.names = FALSE)
    all_eqs  <- unlist(lapply(object@blocks, `[[`, "equations"), use.names = FALSE)
    if (anyDuplicated(all_endo)) {
      problems <- c(problems, "an endogenous variable appears in more than one block")
    }
    if (anyDuplicated(all_eqs)) {
      problems <- c(problems, "an equation appears in more than one block")
    }
    if (!is.null(v$endo) && !setequal(all_endo, v$endo)) {
      problems <- c(problems, "the blocks do not cover every endogenous variable exactly once")
    }
    if (length(object@jacobian) &&
        !setequal(names(object@jacobian), names(object@blocks))) {
      problems <- c(problems, "jacobian and blocks do not describe the same blocks")
    }
  }

  if (length(problems)) problems else TRUE
})

#' @param object a `thor_model`
#' @rdname thor_model-class
#' @export
setMethod("show", "thor_model", function(object) {
  cat("<thor_model> '", object@name, "'  [", object@backend, "]\n", sep = "")
  cat("  ", nrow(object@equations), " equations, ",
      length(object@vars$endo), " endogenous, ",
      length(object@vars$exo), " exogenous, ",
      length(object@vars$coeff), " coefficients\n", sep = "")

  if (length(object@blocks)) {
    for (b in object@blocks) {
      line <- sprintf("  %-9s %5d equations", b$name, length(b$equations))
      j <- object@jacobian[[b$name]]
      if (inherits(j, "thor_sparse_jacobian")) {
        line <- paste0(line, sprintf("   jacobian %d non-zero (%.2f%% dense)",
                                     length(j$i), 100 * length(j$i) / max(1, j$n^2)))
      } else if (is.matrix(j)) {
        nz <- sum(j != "0")
        line <- paste0(line, sprintf("   jacobian %d non-zero (%.2f%% dense)",
                                     nz, 100 * nz / max(1, length(j))))
      }
      cat(line, "\n", sep = "")
    }
  } else {
    cat("  not decomposed\n")
  }

  g <- object@generated
  if (length(g$code)) {
    cat("  generated ", if (object@backend == "dense-r") "R" else "C++",
        ": ", format(nchar(g$code) / 1024, digits = 3), " KB",
        if (length(g$hash)) paste0(", hash ", substr(g$hash, 1, 10)) else "",
        "\n", sep = "")
  }
  if (length(object@meta$built_at)) {
    cat("  built ", format(object@meta$built_at, "%Y-%m-%d %H:%M:%S"),
        " with thortwo ", as.character(object@meta$version), "\n", sep = "")
  }
  invisible(object)
})

## D2 -- tresthor wrote into Deriv's rule table from inside both builders and
## never restored it, and also set `options(stringsAsFactors = FALSE)`, a
## no-op since R 4.0. The rules are genuinely needed, so they are registered
## once here instead of on every build.
##
## `delta(n, x)` is linear in x, so d/dx delta(n, x) = 1: the lagged part is
## carried by `lag.*` symbols, which Deriv treats as unrelated atoms.
## `Deriv::drule` is an environment, so the rules are assigned into it rather
## than through `Deriv::drule[[...]] <- `, which R would read as an assignment
## to `::` itself.
.onLoad <- function(libname, pkgname) {
  rules <- Deriv::drule
  assign("delta",   alist(x = 1, y = NULL), envir = rules)  # y is just a flag
  assign("newdiff", alist(x = 1, y = NULL), envir = rules)
  ## d|x|/dx is sign(x). tresthor registered this as
  ## `ifelse(x == 0, 0, sign(x))`, which is the same function -- sign(0) is
  ## already 0 -- but it puts an `ifelse` into every jacobian entry derived
  ## from an abs(), and no code generator can translate that. The bug was
  ## latent there only because none of the shipped models used abs().
  assign("abs", list(x = quote(sign(x))), envir = rules)
  invisible()
}

## Null-coalescing, for reading optional metadata off a model saved before a
## field existed. Base R gained `%||%` in 4.4; thortwo supports 4.1.
## @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
