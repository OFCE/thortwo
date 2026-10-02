## Diagnosing a broken C++ toolchain.
##
## Two of the three backends compile generated C++, and when that fails the
## compiler error is long, points into system headers, and looks like a
## thortwo bug. Usually it is not: the machine cannot compile *any* Rcpp code.
##
## The usual culprit is a personal Makevars written to work around one
## compiler or SDK update and then forgotten. Seen in practice:
##
##   CXXFLAGS=-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1
##
## added to fix a removed MacOSX14.4 SDK. After the SDK moved to 26.5 and
## clang to 21 it broke libc++'s own header ordering (`FP_NAN` undeclared,
## `std::memcpy` unresolved) -- for every Rcpp compile, not just ours. And
## compile_model_cpp() starts from `R CMD config CXXFLAGS`, which includes the
## user's Makevars, so the flag was carried into every model build too.
##
## So: when a compile fails, say which Makevars are in effect and which of
## their flags look wrong, and offer thor_check_toolchain() to settle whether
## the machine can compile anything at all.

#' Makevars files R applies to compilations, in the order R reads them
#'
#' @param makevars_user value of `R_MAKEVARS_USER` to assume, NA for unset
#' @return a named character vector of existing files (`site`, `user`)
#' @keywords internal
makevars_in_effect <- function(makevars_user = Sys.getenv("R_MAKEVARS_USER", unset = NA)) {

  site <- Sys.getenv("R_MAKEVARS_SITE", unset = "")
  if (!nzchar(site)) {
    site <- file.path(R.home("etc"), Sys.getenv("R_ARCH"), "Makevars.site")
  }

  if (!is.na(makevars_user) && nzchar(makevars_user)) {
    user <- makevars_user
  } else {
    home <- path.expand("~/.R")
    cand <- if (.Platform$OS.type == "windows") {
      file.path(home, c(if (.Machine$sizeof.pointer == 8L) "Makevars.win64", "Makevars.win"))
    } else {
      file.path(home, c(paste0("Makevars-", R.version$platform), "Makevars"))
    }
    user <- cand[file.exists(cand)][1L]
  }

  out <- c(site = site, user = user)
  out[!is.na(out) & file.exists(out)]
}

#' Flag include and sysroot options that are likely to break compilation
#'
#' Two patterns: a hand-added path into a C++ standard library (`.../c++/v1`,
#' `.../include/c++/...`), which shadows the compiler's own and breaks on the
#' next compiler or SDK update; and a path that does not exist, typically an
#' SDK that has since been removed.
#'
#' @param text character vector of compiler flags or Makevars lines
#' @return character vector of problems, one per offending flag
#' @keywords internal
suspicious_flags <- function(text) {
  text <- sub("#.*$", "", text)
  m <- regmatches(text, gregexpr("(-isystem|-isysroot|-I|--sysroot=)\\s*\"?[^\"[:space:]]+", text))
  flags <- unique(trimws(unlist(m)))
  out <- character(0)
  for (f in flags) {
    path <- sub("^(-isystem|-isysroot|-I|--sysroot=)\\s*\"?", "", f)
    if (grepl("$(", path, fixed = TRUE) || grepl("${", path, fixed = TRUE)) next
    if (grepl("c\\+\\+/v1/?$|/include/c\\+\\+(/|$)", path)) {
      out <- c(out, paste0(f, " puts the C++ standard library headers on the include ",
                           "path by hand; this breaks after a compiler or SDK update"))
    } else if (!dir.exists(path.expand(path))) {
      out <- c(out, paste0(f, " points to a directory that does not exist"))
    }
  }
  out
}

#' What can be said about the toolchain without compiling anything
#'
#' @param makevars_user value of `R_MAKEVARS_USER` to assume
#' @return a list: `makevars` (files in effect), `problems` (character)
#' @keywords internal
toolchain_findings <- function(makevars_user = Sys.getenv("R_MAKEVARS_USER", unset = NA)) {
  mv <- makevars_in_effect(makevars_user)
  problems <- character(0)
  for (i in seq_along(mv)) {
    p <- suspicious_flags(readLines(mv[[i]], warn = FALSE))
    if (length(p)) problems <- c(problems, paste0(mv[[i]], ": ", p))
  }
  list(makevars = mv, problems = problems)
}

#' Error message for a failed model compilation
#'
#' @param e the condition raised by `Rcpp::sourceCpp()`
#' @param makevars_user the user's own `R_MAKEVARS_USER`, NA if unset
#' @return a character scalar
#' @keywords internal
compile_failure_message <- function(e, makevars_user) {
  f <- toolchain_findings(makevars_user)
  lines <- c(
    paste0("Compiling the model failed: ", conditionMessage(e)),
    "",
    "If the compiler errors above point into system headers (math.h, cstdint, ...)",
    "rather than into the model's .cpp, the machine's C++ setup is at fault, not",
    "the model. Run thortwo::thor_check_toolchain() to find out.",
    if (length(f$makevars)) c("", "Makevars in effect:", paste0("  ", names(f$makevars), ": ", f$makevars)),
    if (length(f$problems)) c("", "Likely cause:", paste0("  ", f$problems))
  )
  paste(lines, collapse = "\n")
}

#' Check that this machine can compile thortwo models
#'
#' Compiles a few lines of C++ that include Rcpp and RcppEigen, exactly as a
#' model build would, and reports what it found: the compiler, the Makevars
#' files R is applying, any flags in them that are known to break
#' compilation, and, on failure, the first compiler errors.
#'
#' Run it when [thor_model()] fails to compile on the `sparse` or `dense-cpp`
#' backend. If this check fails too, the problem is the machine's setup and
#' not the model; the `dense-r` backend needs no compiler and works meanwhile.
#'
#' @param quiet logical. TRUE to print nothing and only return the result.
#' @return invisibly, a list: `ok` (logical), `problems` (character),
#'   `makevars` (files in effect), `output` (compiler output)
#' @export
thor_check_toolchain <- function(quiet = FALSE) {

  say <- function(...) if (!quiet) cat(..., "\n", sep = "")
  rbin <- file.path(R.home("bin"), "R")
  rconf <- function(v) {
    x <- tryCatch(suppressWarnings(system2(rbin, c("CMD", "config", v),
                                           stdout = TRUE, stderr = TRUE)),
                  error = function(e) character(0))
    trimws(paste(x, collapse = " "))
  }

  cxx <- rconf("CXX")
  cc_version <- if (nzchar(cxx)) {
    tryCatch(suppressWarnings(system2(strsplit(cxx, "\\s+")[[1L]][1L], "--version",
                                      stdout = TRUE, stderr = TRUE))[1L],
             error = function(e) NA_character_)
  } else NA_character_

  say("thortwo toolchain check")
  say("  R          ", as.character(getRversion()), "  (", R.version$platform, ")")
  say("  compiler   ", if (nzchar(cxx)) cxx else "(none configured)",
      if (!is.na(cc_version)) paste0("  [", cc_version, "]") else "")

  has <- vapply(c("Rcpp", "RcppEigen"), requireNamespace, logical(1), quietly = TRUE)
  for (p in names(has)) {
    say("  ", formatC(p, width = -10), " ",
        if (has[[p]]) as.character(utils::packageVersion(p)) else "NOT INSTALLED")
  }

  f <- toolchain_findings()
  ## the flags R actually ends up with, as well as the files they come from
  conf_problems <- suspicious_flags(c(rconf("CPPFLAGS"), rconf("CXXFLAGS")))
  conf_problems <- setdiff(conf_problems,
                           sub("^.*?: ", "", f$problems))
  problems <- c(f$problems, if (length(conf_problems)) paste0("R CMD config: ", conf_problems))

  if (length(f$makevars)) {
    for (i in seq_along(f$makevars)) {
      say("  Makevars   ", names(f$makevars)[i], ": ", f$makevars[[i]],
          if (names(f$makevars)[i] == "user" && !is.na(Sys.getenv("R_MAKEVARS_USER", unset = NA)))
            "  (from R_MAKEVARS_USER)" else "")
      body <- trimws(readLines(f$makevars[[i]], warn = FALSE))
      body <- body[nzchar(body) & !startsWith(body, "#")]
      for (l in body) say("               ", l)
    }
  } else {
    say("  Makevars   none")
  }

  output <- character(0)
  ok <- FALSE
  if (all(has) && nzchar(cxx)) {
    dir <- tempfile("thortwo_toolchain_")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    writeLines(c(
      "#include <cmath>",
      "#include <cstring>",
      "#include <RcppEigen.h>",
      "// [[Rcpp::depends(RcppEigen)]]",
      "// [[Rcpp::export]]",
      "double thortwo_probe(double x) {",
      "  Eigen::VectorXd v = Eigen::VectorXd::Constant(2, x);",
      "  return std::sqrt(v.sum());",
      "}"), file.path(dir, "probe.cpp"))

    incl <- paste0("-I\"", c(system.file("include", package = "Rcpp"),
                             system.file("include", package = "RcppEigen")), "\"",
                   collapse = " ")
    old <- Sys.getenv("PKG_CPPFLAGS", unset = NA)
    Sys.setenv(PKG_CPPFLAGS = incl)
    output <- withr_dir(dir, suppressWarnings(
      system2(rbin, c("CMD", "SHLIB", "probe.cpp"), stdout = TRUE, stderr = TRUE)))
    if (is.na(old)) Sys.unsetenv("PKG_CPPFLAGS") else Sys.setenv(PKG_CPPFLAGS = old)
    ok <- is.null(attr(output, "status")) && file.exists(file.path(dir, paste0("probe", .Platform$dynlib.ext)))
  }

  if (ok) {
    say("  compile    ok")
  } else if (!all(has)) {
    say("  compile    skipped: install the missing packages first")
  } else {
    say("  compile    FAILED -- this machine cannot compile Rcpp code, so the")
    say("             problem is its setup, not the model. First errors:")
    err <- grep("error", output, value = TRUE, ignore.case = TRUE)
    for (l in utils::head(if (length(err)) err else output, 5L)) say("               ", l)
  }

  if (length(problems)) {
    say("")
    say("  Likely cause:")
    for (p in problems) say("    ", p)
    say("  Remove or fix that line (or unset R_MAKEVARS_USER if it points there),")
    say("  then restart R.")
  } else if (!ok && all(has)) {
    say("")
    say("  No known-bad flags found. Check that the compiler itself works")
    say(switch(Sys.info()[["sysname"]],
      Darwin  = "  (`xcode-select --install`, or reinstall the Command Line Tools).",
      Windows = "  (that the Rtools matching your R version is installed).",
                "  (that g++ or clang++ and R's development headers are installed)."))
  }

  invisible(list(ok = ok, problems = problems, makevars = f$makevars, output = output))
}

#' Evaluate an expression with a different working directory
#'
#' @param dir directory
#' @param code expression
#' @keywords internal
withr_dir <- function(dir, code) {
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  force(code)
}
