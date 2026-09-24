## Saving and loading a model.
##
## tresthor's `save_model()` wrote the object with `saveRDS`, and
## `load_model()` read it back and `sourceCpp`d `model@rcpp_source` -- an
## *absolute path* recorded at build time, which pointed into whatever
## temporary directory the model happened to be built in. A saved model was
## therefore machine-local and usually session-local: it could not be mailed
## to a colleague, committed, or reused after a reboot cleared /tmp.
##
## Here the generated code travels inside the .rds, as a string. `generated$path`
## is a hint; `generated$code` and `generated$hash` are authoritative. That makes
## "saving or not saving" genuinely free: a model just built and one read back
## from disk are the same object, and whether compilation actually happens is
## decided by the cache, not by which of the two you are holding.

#' Save a model
#'
#' The result is self-contained: it carries the generated source, so it can be
#' copied to another machine and loaded there.
#'
#' @param model a `thor_model`
#' @param path file to write. A directory is also accepted, in which case the
#'   file is named after the model.
#' @return `path`, invisibly
#' @export
thor_save <- function(model, path) {
  if (!methods::is(model, "thor_model")) {
    stop("`model` must be a thor_model.", call. = FALSE)
  }
  if (!length(model@generated$code)) {
    stop("Model '", model@name, "' carries no generated code, so it cannot be ",
         "saved in a self-contained form. Rebuild it with thor_model().", call. = FALSE)
  }

  if (dir.exists(path)) path <- file.path(path, paste0(model@name, ".rds"))
  if (!grepl("\\.rds$", path, ignore.case = TRUE)) path <- paste0(path, ".rds")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  saveRDS(model, file = path)
  invisible(path)
}

#' Load a saved model
#'
#' @param path the `.rds` written by [thor_save()]
#' @param compile logical. TRUE (the default) to make the model ready to solve
#'   before returning. The generated code is written to `workdir` and compiled,
#'   which hits the compile cache whenever the code has not changed -- so
#'   loading a model that this machine has compiled before costs a fraction of
#'   a second, not a full build.
#' @param workdir directory for the generated source. Defaults to the same
#'   stable, model-derived location a build would use, which is what lets the
#'   cache recognise it.
#' @param cache see [thor_solve()]
#' @return a `thor_model`
#' @export
thor_load <- function(path, compile = TRUE, workdir = NULL, cache = NULL) {

  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  model <- readRDS(path)
  if (!methods::is(model, "thor_model")) {
    stop("'", basename(path), "' does not hold a thor_model.", call. = FALSE)
  }
  methods::validObject(model)

  ## Point the model at a work directory on *this* machine before anything
  ## tries to use the path recorded at build time.
  if (is.null(workdir)) workdir <- thor_workdir(model@name)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  model@generated$path <- file.path(
    normalizePath(workdir, mustWork = TRUE),
    paste0(model@name, if (model@backend == "dense-r") ".R" else ".cpp"))

  if (compile) model_env(model, cache = cache)
  model
}

#' Export a model back to a `.txt` model file
#'
#' Round-trips with [read_model_source()], so
#' `thor_model(export_model(m))` reproduces `m`'s equations table. That makes
#' it the cheapest available regression test on the parser.
#'
#' @param model a `thor_model`
#' @param filename file to write
#' @return `filename`, invisibly
#' @export
export_model <- function(model, filename = "model.txt") {

  if (!methods::is(model, "thor_model")) {
    stop("`model` must be a thor_model.", call. = FALSE)
  }
  if (!grepl("\\.[[:alnum:]]+$", filename)) filename <- paste0(filename, ".txt")
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  eq <- model@equations
  ## A generated id carries no information, so it is not written back out; a
  ## user-given name is.
  written <- ifelse(eq$name == eq$id, eq$equation, paste0(eq$name, ":", eq$equation))

  writeLines(c(
    "endogenous variables :",
    paste(model@vars$endo, collapse = ","),
    "##############",
    "exogenous variables :",
    paste(model@vars$exo, collapse = ","),
    "##############",
    "coefficients :",
    paste(model@vars$coeff, collapse = ","),
    "##############",
    "equations :",
    written), filename)

  invisible(filename)
}

#' Path to a sample model file
#'
#' @param path_only TRUE to return the path instead of opening the file
#' @return the path to the example model
#' @examples
#' model_source_example(TRUE)
#' @export
model_source_example <- function(path_only = FALSE) {
  p <- system.file("models", "model_type.txt", package = "thortwo")
  if (path_only) return(p)
  utils::browseURL(p)
  invisible(p)
}
