## The build pipeline: text file in, ready-to-solve model out.
##
## tresthor had two builders, `create_model()` (357 lines) and
## `create_model_sparse()` (216), which were near-identical up to step 4 and
## diverged only in how the jacobian was represented and which generator was
## called. They drifted apart, and a fix to one never reached the other. Here
## there is one function and a `backend` argument: steps 1-3 run once, and
## only steps 4 and 5 branch.

#' Build a model
#'
#' Reads a model, checks it, decomposes it into blocks, differentiates it
#' symbolically, generates a model-specific solver and compiles it.
#'
#' @param name character. Name of the model.
#' @param source path to the `.txt` model file. If NULL, the model is taken
#'   from `endogenous` / `exogenous` / `coefficients` / `equations`.
#' @param endogenous,exogenous,coefficients,equations character vectors, used
#'   only when `source` is NULL.
#' @param backend one of `"sparse"` (the default), `"dense-cpp"` or
#'   `"dense-r"`. See [thor_backends].
#' @param decompose logical. TRUE (the default) to split the model into
#'   prologue, heart and epilogue; FALSE to solve it as one block.
#' @param workdir directory for the generated source. Defaults to a stable,
#'   model-derived location under [thor_workdir()], which is what makes the
#'   compile cache work; pass a `tempfile()` to opt out of it.
#' @param compile logical. TRUE (the default) to compile or evaluate the
#'   generated code immediately, so that the returned model is ready to solve.
#' @param cache directory in which to cache the compiled object between
#'   sessions, FALSE to disable caching, or NULL (the default) for
#'   [thor_cache_dir()].
#' @param verbose logical. TRUE (the default) to report progress.
#'
#' @return a [thor_model-class] object
#' @examples
#' \dontrun{
#' m <- thor_model("opale", system.file("models", "opale.txt", package = "thortwo"))
#' m
#' }
#' @export
thor_model <- function(name,
                       source = NULL,
                       endogenous = NULL,
                       exogenous = NULL,
                       coefficients = NULL,
                       equations = NULL,
                       backend = c("sparse", "dense-cpp", "dense-r"),
                       decompose = TRUE,
                       workdir = NULL,
                       compile = TRUE,
                       cache = NULL,
                       verbose = TRUE) {

  t_start <- Sys.time()
  backend <- match.arg(backend)
  assertthat::assert_that(assertthat::is.string(name), nzchar(name))
  assertthat::assert_that(is.logical(decompose), length(decompose) == 1L)

  say <- function(...) if (verbose) cat(..., sep = "")
  say("Building model '", name, "' (", backend, ")\n\n")

  ################################
  #### 1. Read and check
  ################################
  say("Step 1: reading and checking the model...\n")

  if (!is.null(source)) {
    parts <- read_model_source(source)
    eqlist <- parts$equations
    endo   <- parts$endogenous
    exo    <- parts$exogenous
    coeff  <- parts$coefficients
    source_hash <- unname(tools::md5sum(source))
  } else {
    if (is.null(endogenous) || is.null(equations)) {
      stop("No `source` was given and no `equations` were specified. Provide either a ",
           "model file, or the equations together with their endogenous variables.",
           call. = FALSE)
    }
    eqlist <- unique(equations)
    endo   <- tolower(unique(endogenous))
    exo    <- tolower(unique(exogenous))
    coeff  <- tolower(unique(coefficients))
    source_hash <- NA_character_
  }

  endo <- sort(endo); exo <- sort(exo); coeff <- sort(coeff)
  eqlist <- gsub("mylg\\(", "lag(", eqlist)

  check_equation_input(eqlist)
  check_var_vector(endo,  "Endogenous variables")
  check_var_vector(exo,   "Exogenous variables")
  check_var_vector(coeff, "Coefficient variables")
  check_variable_conflict(endo, exo,   "the endogenous", "the exogenous variables")
  check_variable_conflict(endo, coeff, "the endogenous variables", "the coefficients")
  check_variable_conflict(exo,  coeff, "the exogenous variables", "the coefficients")

  exo   <- is_in_formulas(exo,   eqlist, "exogenous",   verbose)
  endo  <- is_in_formulas(endo,  eqlist, "endogenous",  verbose)
  coeff <- is_in_formulas(coeff, eqlist, "coefficient", verbose)
  all_model_variables <- sort(c(endo, exo, coeff))

  equations_list <- create_equations_list(eqlist, verbose)
  if (!parser_lag_delta_check(equations_list$equation)) {
    stop("Parser error on the lags and/or deltas. Please check the model's formulas.",
         call. = FALSE)
  }

  ################################
  #### 2. Endogenous variables per equation
  ################################
  say("Step 2: identifying the endogenous variables in the equations...\n")
  eqns <- table_contemporaneous_endos(formula_list = equations_list$formula,
                                      endogenous = endo, exogenous = exo,
                                      coefflist = coeff,
                                      equations_index = equations_list$id)
  check_eq_var_identification(endo = endo, eqns = eqns,
                              names_of_equations = equations_list$name)

  ################################
  #### 3. Decomposition
  ################################
  say("Step 3: ", if (decompose) "decomposing the model into blocks..." else
                  "creating the single block of the model...", "\n")
  d <- decomposing_model(endogenous_variables = endo, eq_var_matrix = eqns,
                         decomposition = decompose, verbose = verbose)

  equations_list$part <- "tbd"
  equations_list$part[equations_list$id %in% d$prologue_equations] <- "prologue"
  equations_list$part[equations_list$id %in% d$heart_equations]    <- "heart"
  equations_list$part[equations_list$id %in% d$epilogue_equations] <- "epilogue"
  if ("tbd" %in% equations_list$part) {
    stop("Some equations were not assigned to a block by the decomposition: ",
         paste(utils::head(equations_list$id[equations_list$part == "tbd"], 5),
               collapse = ", "), call. = FALSE)
  }
  if (verbose) print(table(equations_list$part))

  equations_list$new_formula <- formatting_formulas(equations_list$formula)

  ## Blocks are kept in solve order; an empty block is simply absent.
  block_spec <- list(
    prologue = list(present = d$prologue, endo = sort(d$prologue_endo),
                    equations = d$prologue_equations),
    heart    = list(present = d$heart,    endo = sort(d$heart_endo),
                    equations = d$heart_equations),
    epilogue = list(present = d$epilogue, endo = sort(d$epilogue_endo),
                    equations = d$epilogue_equations))
  block_spec <- block_spec[vapply(block_spec,
                                  function(b) isTRUE(b$present) && length(b$equations) > 0L,
                                  logical(1))]
  blocks <- lapply(names(block_spec), function(nm) {
    c(list(name = nm), block_spec[[nm]][c("present", "endo", "equations")])
  })
  names(blocks) <- names(block_spec)

  ################################
  #### 4. Symbolic jacobians
  ################################
  say("\nStep 4: computing the symbolic jacobians...\n")
  eqns_as_list <- lapply(split(as.matrix(eqns), row(eqns)),
                         function(x) unique(stats::na.omit(x)))
  names(eqns_as_list) <- rownames(eqns)

  sparse_path <- identical(backend, "sparse")

  jacobian <- lapply(blocks, function(b) {
    if (sparse_path) {
      symbolic_jacobian_sparse(equations_list_df = equations_list,
                               eqns_vars_list = eqns_as_list,
                               endo_vec = b$endo, equations_subset = b$equations)
    } else {
      symbolic_jacobian(equations_list_df = equations_list,
                        eqns_vars_list = eqns_as_list,
                        endo_vec = b$endo, equations_subset = b$equations)
    }
  })
  names(jacobian) <- names(blocks)

  ## The generators all consume triplets; the dense backends keep the matrix
  ## on the object, because that is what "dense jacobian" means to a caller.
  gen_blocks <- lapply(names(blocks), function(nm) {
    j <- jacobian[[nm]]
    tj <- if (sparse_path) j else as_sparse_jacobian(j)
    pos <- match(tj$equations, as.character(equations_list$id))
    list(jac = tj, formulas = equations_list$new_formula[pos])
  })
  names(gen_blocks) <- names(blocks)

  if (verbose) {
    for (nm in names(gen_blocks)) {
      cat("   ", nm, ": ", sep = ""); print(gen_blocks[[nm]]$jac)
    }
  }

  ################################
  #### 5. Code generation
  ################################
  say("\nStep 5: generating the ", if (backend == "dense-r") "R" else "C++",
      " solver...\n")

  code <- switch(backend,
    "sparse"    = generate_cpp_sparse(name, gen_blocks, all_model_variables, verbose),
    "dense-cpp" = generate_cpp_dense(name, gen_blocks, all_model_variables, verbose),
    "dense-r"   = generate_r_source(name, gen_blocks, all_model_variables, verbose))

  if (is.null(workdir)) workdir <- thor_workdir(name)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(normalizePath(workdir, mustWork = TRUE),
                    paste0(name, if (backend == "dense-r") ".R" else ".cpp"))
  path <- write_if_changed(code, path)
  say("   ", if (isTRUE(attr(path, "unchanged"))) "unchanged: " else "written to ",
      as.character(path), " (", round(nchar(code) / 1024), " KB)\n")

  ################################
  #### 6. The model object
  ################################
  model <- methods::new("thor_model",
    name      = name,
    backend   = backend,
    equations = equations_list,
    vars      = list(endo = endo, exo = exo, coeff = coeff,
                     all = all_model_variables),
    blocks    = blocks,
    jacobian  = jacobian,
    generated = list(path = as.character(path), code = code, hash = thor_hash(code)),
    meta      = list(version = utils::packageVersion("thortwo"),
                     built_at = Sys.time(),
                     source = if (is.null(source)) NA_character_ else normalizePath(source),
                     source_hash = source_hash,
                     decompose = decompose))
  methods::validObject(model)

  ################################
  #### 7. Compile
  ################################
  if (compile) {
    say("\nStep 6: ", if (backend == "dense-r") "evaluating" else "compiling", "...\n")
    el <- system.time(model_env(model, cache = cache))[["elapsed"]]
    say("   ready in ", round(el, 1), " s", if (el < 1) "  (from cache)" else "", "\n")
  }

  say("\nModel built in ",
      round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), " s\n")
  model
}
