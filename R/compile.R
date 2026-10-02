## Compilation, loading and caching of generated model code.
##
## Generated model code is unlike hand-written code: a few thousand very large
## arithmetic expressions with no control flow. One of R's default compiler
## settings behaves pathologically on it.
##
## Measured on ThreeME 4x4 (1729 equations, 370 KB of generated C++):
##
##   R's default flags (-g -O2)      299 s
##   same without -g                  12 s
##
## `-g` makes clang emit debug metadata for every subexpression, which is
## worthless here -- nobody steps through generated code in a debugger -- and
## costs 25x the compile time. The effect grows with model size, so on a 25k
## equation model it is the difference between minutes and an hour.

## D3 -- Loaded models, keyed by the hash of their generated code.
##
## Every generated file exports functions with the same names
## (`thor_cpp_solve`, `thor_cpp_residuals`), so each model is loaded into its
## own environment rather than into the global one. tresthor's classic backend
## called `Rcpp::sourceCpp(model@rcpp_source)` on *every solve*, loading
## `Rcpp_solver` into the global environment: with two classic models in one
## session the second silently overwrote the first, and the first was then
## solved with the wrong code. Both compiled backends go through this registry.
##
## Keying on the code hash rather than on the file path is what makes a model
## loaded from a .rds and a model just built indistinguishable: same
## equations, same environment, no second compilation.
.thor_loaded <- new.env(parent = emptyenv())

#' Directory used to cache compiled models between sessions
#'
#' Without a persistent cache, `Rcpp::sourceCpp()` builds into `tempdir()` and
#' so recompiles every model once per R session -- around 7 s for Opale, 14 s
#' for ThreeME 4x4 and 30 s for 8x8, and considerably more for larger models.
#' With one, the compilation happens once per machine and per version of the
#' model.
#'
#' Rcpp keys its cache on the hash of the source, the platform and its own
#' version, so an edited model or a changed toolchain rebuilds by itself. It
#' does not key on the R version, and R binaries are not compatible across
#' minor releases, so that is added here.
#'
#' Set `options(thortwo.cache.dir = "...")` to move it, or
#' `options(thortwo.cache.dir = FALSE)` to disable caching entirely.
#'
#' @return the cache directory, created if needed, or NULL if caching is off
#' @export
thor_cache_dir <- function() {

  d <- getOption("thortwo.cache.dir", NULL)
  if (isFALSE(d)) return(NULL)

  if (is.null(d)) {
    d <- tryCatch(tools::R_user_dir("thortwo", "cache"),
                  error = function(e) file.path(tempdir(), "thortwo-cache"))
  }

  d <- file.path(path.expand(d), paste0("R-", getRversion()))
  if (!dir.exists(d)) {
    ok <- dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(d)) {
      warning("Could not create the thortwo cache directory '", d,
              "'. Models will be recompiled each session.", call. = FALSE)
      return(NULL)
    }
  }
  d
}

#' Default directory for a model's generated source
#'
#' D5 -- Rcpp keys its compile cache on the source *path*, and tresthor's
#' builders wrote to whatever `rcpp_path` they were given: every caller,
#' including both of its own test scripts, passed a fresh `tempfile()`. The
#' result was a full recompile on every run, plus a ~1.3 MB cache entry left
#' behind each time. Measured: the same model built twice with a stable path
#' costs 7.0 s then 0.3 s; with tempdir paths it is 7 s every time.
#'
#' So the default is a stable, model-derived location. A temporary directory
#' is now something you opt into by passing `workdir` explicitly.
#'
#' @param name model name
#' @return a directory path, created if needed
#' @export
thor_workdir <- function(name) {
  base <- tryCatch(tools::R_user_dir("thortwo", "cache"),
                   error = function(e) file.path(tempdir(), "thortwo-cache"))
  d <- file.path(path.expand(base), "src", name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Remove cached compiled models
#'
#' @param confirm logical. FALSE to skip the interactive confirmation.
#' @return the number of files removed, invisibly
#' @export
clear_model_cache <- function(confirm = interactive()) {
  d <- thor_cache_dir()
  if (is.null(d) || !dir.exists(d)) {
    cat("Nothing cached.\n")
    return(invisible(0L))
  }
  files <- list.files(d, recursive = TRUE, full.names = TRUE)
  size <- sum(file.info(files)$size, na.rm = TRUE)
  cat("Cache: ", d, "\n", length(files), " files, ",
      round(size / 1024^2, 1), " MB\n", sep = "")
  if (confirm && !isTRUE(utils::askYesNo("Delete?"))) return(invisible(0L))
  unlink(d, recursive = TRUE)
  ## also drop anything loaded in this session, so the next solve rebuilds
  rm(list = ls(.thor_loaded, all.names = TRUE), envir = .thor_loaded)
  invisible(length(files))
}

#' Hash a string
#'
#' Used to key the in-session registry and to tell whether a saved model's
#' code is the one already compiled. `tools::md5sum()` only hashes files, so
#' the string goes through a temporary one.
#'
#' @param x character scalar
#' @return a 32-character hex string
#' @keywords internal
thor_hash <- function(x) {
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  writeChar(x, f, eos = NULL)
  unname(tools::md5sum(f))
}

#' Compile a generated model source file and return its functions
#'
#' Compiles with the platform's usual flags minus `-g`, by way of a temporary
#' `R_MAKEVARS_USER` file. The user's own `~/.R/Makevars` is left untouched.
#'
#' @param path path to the generated .cpp file
#' @param key registry key; defaults to the hash of the file's contents
#' @param rebuild logical. TRUE to compile from scratch, ignoring both the
#'   in-session and the on-disk cache.
#' @param cache directory in which to cache the compiled object between
#'   sessions, FALSE to disable, or NULL (the default) for [thor_cache_dir()].
#' @param debug logical. TRUE to keep `-g` (much slower; only useful when
#'   debugging the code generator itself). Default FALSE.
#' @param quiet logical. TRUE to suppress compiler output. Default TRUE.
#' @return an environment holding the model's compiled functions
#' @keywords internal
compile_model_cpp <- function(path, key = NULL, rebuild = FALSE, cache = NULL,
                              debug = FALSE, quiet = TRUE) {

  stopifnot(file.exists(path))

  ## The generated source carries `// [[Rcpp::depends(RcppEigen)]]`, which
  ## Rcpp resolves against the installed package at compile time. Check it
  ## here so the failure is a clear message rather than a compiler error
  ## about a missing Eigen header.
  if (!requireNamespace("RcppEigen", quietly = TRUE)) {
    stop("Package 'RcppEigen' is required to compile a model. ",
         "Install it with install.packages(\"RcppEigen\").", call. = FALSE)
  }

  if (is.null(key)) key <- thor_hash(readChar(path, file.size(path), useBytes = TRUE))

  if (!rebuild && !is.null(.thor_loaded[[key]])) return(.thor_loaded[[key]])

  env <- new.env(parent = globalenv())

  ## the user's own setting, for the diagnosis should the compile fail
  user_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA)

  if (!debug) {
    ## Start from the platform's own flags so we keep -arch, -falign-functions
    ## and anything else the build was configured with, and only drop -g.
    cxxflags <- tryCatch(system2("R", c("CMD", "config", "CXXFLAGS"),
                                 stdout = TRUE, stderr = FALSE),
                         error = function(e) character(0))
    cxxflags <- paste(cxxflags, collapse = " ")
    if (!nzchar(trimws(cxxflags))) cxxflags <- "-O2"
    flags <- strsplit(trimws(cxxflags), "\\s+")[[1L]]
    flags <- flags[!grepl("^-g([0-9]|gdb.*)?$", flags)]   # keep an explicit -g0
    cxxflags <- paste(c(flags, "-g0"), collapse = " ")

    mk <- tempfile(pattern = "thortwo_makevars_")
    writeLines(c(paste0("CXXFLAGS = ", cxxflags),
                 paste0("CXX17FLAGS = ", cxxflags),
                 paste0("CXX20FLAGS = ", cxxflags)), mk)

    Sys.setenv(R_MAKEVARS_USER = mk)
    on.exit({
      if (is.na(user_makevars)) Sys.unsetenv("R_MAKEVARS_USER")
      else Sys.setenv(R_MAKEVARS_USER = user_makevars)
      unlink(mk)
    }, add = TRUE)
  }

  cache_dir <- if (is.null(cache)) thor_cache_dir() else if (isFALSE(cache)) NULL else cache
  if (is.null(cache_dir)) cache_dir <- tempdir()

  ## A failure here is far more often the machine's toolchain than the
  ## generated code (see toolchain.R), so say so instead of leaving the user
  ## with a page of errors from system headers.
  tryCatch(
    Rcpp::sourceCpp(path, env = env, rebuild = rebuild,
                    cacheDir = cache_dir, verbose = !quiet),
    error = function(e) stop(compile_failure_message(e, user_makevars), call. = FALSE)
  )

  if (!exists("thor_cpp_solve", envir = env, inherits = FALSE)) {
    stop("'", basename(path), "' does not define thor_cpp_solve(). ",
         "Was it generated by thor_model()?", call. = FALSE)
  }

  .thor_loaded[[key]] <- env
  env
}

#' Evaluate a model's generated R source into its own environment
#'
#' D1 -- nothing is written to disk and nothing is `source()`d from the
#' working directory. `keep.source = FALSE` matters: with source references on,
#' parsing the thousands of generated statements of a real model costs more
#' than the filesystem round-trip it replaces.
#'
#' @param code character scalar of generated R source
#' @param key registry key; defaults to the hash of `code`
#' @return an environment holding the model's `blocks` list
#' @keywords internal
load_model_r <- function(code, key = NULL) {
  if (is.null(key)) key <- thor_hash(code)
  if (!is.null(.thor_loaded[[key]])) return(.thor_loaded[[key]])

  env <- new.env(parent = getNamespace("thortwo"))
  eval(parse(text = code, keep.source = FALSE), envir = env)

  if (!exists("blocks", envir = env, inherits = FALSE)) {
    stop("The generated R source does not define the model's blocks.", call. = FALSE)
  }
  .thor_loaded[[key]] <- env
  env
}

#' Make a built or loaded model ready to solve
#'
#' Compiles or evaluates the model's generated code if this session has not
#' seen it yet, and returns the environment holding its functions. A model
#' just built and the same model read back from a `.rds` reach this function
#' with the same code and the same hash, so only the first of them pays.
#'
#' @param model a `thor_model`
#' @param cache see [compile_model_cpp()]
#' @param rebuild logical. TRUE to compile from scratch.
#' @return an environment
#' @keywords internal
model_env <- function(model, cache = NULL, rebuild = FALSE) {

  g <- model@generated
  if (!length(g$code)) {
    stop("Model '", model@name, "' carries no generated code. Rebuild it with thor_model().",
         call. = FALSE)
  }
  key <- g$hash

  if (!rebuild && !is.null(.thor_loaded[[key]])) return(.thor_loaded[[key]])

  if (model@backend == "dense-r") return(load_model_r(g$code, key = key))

  ## `generated$path` is only a hint: the file may be gone, or may be somebody
  ## else's. `generated$code` is what the model actually is, so the file is
  ## rewritten from it whenever the two disagree. Writing only on a difference
  ## keeps the modification time, and so Rcpp's cache entry, intact.
  path <- g$path
  if (!length(path) || is.na(path)) path <- file.path(thor_workdir(model@name),
                                                      paste0(model@name, ".cpp"))
  path <- write_if_changed(g$code, path)

  compile_model_cpp(as.character(path), key = key, cache = cache, rebuild = rebuild)
}
