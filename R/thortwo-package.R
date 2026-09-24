#' thortwo: build and solve large macroeconomic models
#'
#' `thortwo` is the model solver extracted from `tresthor`: it parses a model
#' written as a text file, decomposes it into blocks, differentiates it
#' symbolically, generates and compiles a model-specific Newton solver, and
#' solves it against a database.
#'
#' The whole package is four functions:
#'
#' * [thor_model()] builds a model, on one of three interchangeable backends.
#' * [thor_solve()] solves it.
#' * [thor_save()] / [thor_load()] move a built model between sessions and
#'   machines, generated source and all.
#' * [model_residuals()] checks that a solution really does satisfy the
#'   equations, per block or per equation; [calibration_check()] names the
#'   equations a dataset fails at a given period.
#'
#' @keywords internal
#' @importFrom methods new validObject is setClass setValidity setMethod show
#' @importFrom stats na.omit setNames median
#' @importFrom utils head packageVersion askYesNo browseURL
#' @importFrom tools md5sum R_user_dir
#' @importFrom Deriv Deriv
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
"_PACKAGE"
