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
#'   `"sparse-r"` (no compiler needed). `"dense-r"` is accepted as the old
#'   name of `"sparse-r"`. See [thor_backends].
#' @param decompose logical. TRUE (the default) to split the model into
#'   prologue, heart and epilogue; FALSE to solve it as one block.
#' @param sequential logical. FALSE (the default) to solve every block by
#'   Newton's method. TRUE to solve the prologue and the epilogue one equation
#'   at a time, in dependency order, which they allow because they are
#'   recursive; the heart is still solved by Newton. See the section below.
#' @param workdir directory for the generated source. Defaults to a stable,
#'   model-derived location under [thor_workdir()], which is what makes the
#'   compile cache work; pass a `tempfile()` to opt out of it.
#' @param compile logical. TRUE (the default) to compile or evaluate the
#'   generated code immediately, so that the returned model is ready to solve.
#' @param cache directory in which to cache the built model and its compiled
#'   code between sessions, FALSE to disable caching, or NULL (the default) for
#'   [thor_cache_dir()].
#' @param recompile logical. FALSE (the default) to reuse a cached build of
#'   this exact model if there is one; TRUE to build and compile from scratch,
#'   replacing it. See the section on the cache.
#' @param verbose logical. TRUE (the default) to report progress.
#'
#' @section Sequential prologue and epilogue:
#' An advanced option. With `sequential = TRUE`, each equation of the prologue
#' and of the epilogue is solved for the one variable it determines: by
#' evaluating it directly when it is written `x = f(...)`, and by a
#' one-variable Newton iteration otherwise. The solution is the same, to
#' within the solver's tolerance. What changes is the build: no jacobian is
#' derived or compiled for these blocks, which shortens it, most of all for
#' models with very long equations in the epilogue. It has no effect when
#' `decompose = FALSE`.
#'
#' @section The cache:
#' A model is identified by a hash of everything that determines the result:
#' its name, its variables and equations, `backend`, `decompose`, and
#' thortwo's own code. If a model with the same hash has been built before, on
#' this computer and with this version of R, it is loaded from the cache in a
#' fraction of a second instead of being built again, and a message says so.
#' Changing an equation, a variable list or an option changes the hash, so a
#' stale model is never returned; changing only spacing, or upper and lower
#' case, does not.
#'
#' `recompile = TRUE` ignores the cache for this call. `cache = FALSE`, or
#' `options(thortwo.model.cache = FALSE)`, switch the lookup off, and
#' [clear_model_cache()] empties the cache.
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
                       backend = c("sparse", "dense-cpp", "sparse-r"),
                       decompose = TRUE,
                       sequential = FALSE,
                       workdir = NULL,
                       compile = TRUE,
                       cache = NULL,
                       recompile = FALSE,
                       verbose = TRUE) {

  t_start <- Sys.time()
  backend <- match.arg(normalise_backend(backend[1L]), thor_backends)
  assertthat::assert_that(assertthat::is.string(name), nzchar(name))
  assertthat::assert_that(is.logical(decompose), length(decompose) == 1L)
  assertthat::assert_that(is.logical(sequential), length(sequential) == 1L)

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

  ## Has this exact model been built before? Looked up before any of the
  ## expensive steps: on ThreeME 13x13 a hit costs 0.55 s where the build costs
  ## 30 s plus the compile.
  store <- model_cache_file(cache, model_cache_key(name, backend, decompose, sequential,
                                                   endo, exo, coeff, eqlist))
  if (!is.null(store) && !isTRUE(recompile)) {
    cached <- model_cache_read(store)
    if (!is.null(cached)) {
      message("This version of the model has already been compiled before; ",
              "using the existing cache.\n",
              "To recompile, use recompile = TRUE and rerun the build.")
      if (!is.null(source)) cached@meta$source <- normalizePath(source)
      return(attach_model(cached, compile = compile, workdir = workdir, cache = cache))
    }
  }

  check_equation_input(eqlist)
  check_var_vector(endo,  "Endogenous variables")
  check_var_vector(exo,   "Exogenous variables")
  check_var_vector(coeff, "Coefficient variables")
  check_variable_conflict(endo, exo,   "the endogenous", "the exogenous variables")
  check_variable_conflict(endo, coeff, "the endogenous variables", "the coefficients")
  check_variable_conflict(exo,  coeff, "the exogenous variables", "the coefficients")

  used  <- variables_in_formulas(eqlist)
  exo   <- is_in_formulas(exo,   eqlist, "exogenous",   verbose, present = used)
  endo  <- is_in_formulas(endo,  eqlist, "endogenous",  verbose, present = used)
  coeff <- is_in_formulas(coeff, eqlist, "coefficient", verbose, present = used)
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

  ## How many complete observations the solver will need before the first
  ## solved period. Cheap to work out here, and it turns an "inf residual"
  ## at solve time into a sentence the caller can act on.
  maxlag <- model_max_lag(equations_list$new_formula)
  say("   the model reaches ", maxlag$k, " period(s) back",
      if (maxlag$variable) ", plus at least one lag given as a variable" else "", "\n")

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
  eqns_as_list <- lapply(split(as.matrix(eqns), row(eqns)),
                         function(x) unique(stats::na.omit(x)))
  names(eqns_as_list) <- rownames(eqns)

  ## Advanced: the prologue and the epilogue solved one equation at a time
  ## (sequential.R). Planned here, before the jacobians, because a block that
  ## is solved sequentially does not need one.
  seq_plans <- list()
  if (sequential) {
    if (!decompose) {
      say("\n`sequential` has no effect without the decomposition: there is no ",
          "prologue or epilogue.\n")
    }
    for (nm in intersect(c("prologue", "epilogue"), names(blocks))) {
      p <- sequential_plan(blocks[[nm]], equations_list, eqns_as_list)
      if (is.null(p)) {
        warning("The ", nm, " could not be put in a solving order, so it will be ",
                "solved by Newton as usual.", call. = FALSE)
        next
      }
      seq_plans[[nm]] <- p
      say("   ", nm, " solved sequentially: ", sum(p$order$direct), " equations directly, ",
          sum(!p$order$direct), " by one-variable Newton\n")
    }
  }

  say("\nStep 4: computing the symbolic jacobians...\n")

  sparse_path <- backend %in% c("sparse", "sparse-r")

  jacobian <- lapply(blocks, function(b) {
    if (!is.null(seq_plans[[b$name]])) return(empty_jacobian(b))
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

  gen_blocks <- codegen_blocks(jacobian, equations_list, seq_plans)

  if (verbose) {
    for (nm in names(gen_blocks)) {
      cat("   ", nm, ": ", sep = ""); print(gen_blocks[[nm]]$jac)
    }
  }

  ################################
  #### 5. Code generation
  ################################
  say("\nStep 5: generating the ", if (is_r_backend(backend)) "R" else "C++",
      " solver...\n")

  code <- switch(backend,
    "sparse"    = generate_cpp_sparse(name, gen_blocks, all_model_variables, verbose),
    "dense-cpp" = generate_cpp_dense(name, gen_blocks, all_model_variables, verbose),
    "sparse-r"  = generate_r_source(name, gen_blocks, all_model_variables, verbose))

  if (is.null(workdir)) workdir <- thor_workdir(name)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(normalizePath(workdir, mustWork = TRUE),
                    paste0(name, if (is_r_backend(backend)) ".R" else ".cpp"))
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
                     decompose = decompose,
                     sequential = lapply(seq_plans, `[[`, "order"),
                     max_lag = maxlag$k,
                     variable_lag = maxlag$variable))
  methods::validObject(model)

  ################################
  #### 7. Compile
  ################################
  if (compile) {
    say("\nStep 6: ", if (is_r_backend(backend)) "evaluating" else "compiling", "...\n")
    el <- system.time(model_env(model, cache = cache,
                                rebuild = isTRUE(recompile)))[["elapsed"]]
    say("   ready in ", round(el, 1), " s", if (el < 1) "  (from cache)" else "", "\n")
  }

  model_cache_write(model, store)

  say("\nModel built in ",
      round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), " s\n")
  model
}

#' Pair each block's jacobian with its formulas, for the code generators
#'
#' The generators all consume triplets; the dense backend keeps the matrix on
#' the object, because that is what "dense jacobian" means to a caller.
#'
#' @param jacobian named list, one jacobian per block (triplets or matrix)
#' @param equations_list the model's equations table
#' @param seq_plans named list of [sequential_plan()]s, for the blocks that
#'   are solved sequentially
#' @return named list of `list(jac, formulas)`, with a `seq` element for the
#'   sequential blocks
#' @keywords internal
codegen_blocks <- function(jacobian, equations_list, seq_plans = list()) {
  out <- lapply(names(jacobian), function(nm) {
    j <- jacobian[[nm]]
    tj <- if (inherits(j, "thor_sparse_jacobian")) j else as_sparse_jacobian(j)
    pos <- match(tj$equations, as.character(equations_list$id))
    b <- list(jac = tj, formulas = equations_list$new_formula[pos])
    if (!is.null(seq_plans[[nm]])) b$seq <- seq_plans[[nm]]
    b
  })
  names(out) <- names(jacobian)
  out
}
