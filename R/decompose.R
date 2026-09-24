## Block decomposition.
##
## The algorithm is tresthor's (1_3_decomposition_algo.R), which works. It is
## implemented with running counts rather than by recounting on every pass;
## see the note in decomposing_model().
##
## The model is split into three blocks that can be solved in sequence rather
## than all at once:
##
##   prologue -- equations that, once the previously determined ones are
##               substituted in, contain a single undetermined endogenous
##               variable. Solved first, recursively.
##   epilogue -- equations holding an endogenous variable that occurs nowhere
##               else. Those variables are determined by the rest of the model,
##               so they can be computed last. Also recursive.
##   heart    -- whatever is left: a genuinely simultaneous system.
##
## Newton cost is cubic in block size, so on Opale this turns one 496 x 496
## solve into 201 + 98 + 197.

#' Decompose a model into prologue, heart and epilogue
#'
#' @param endogenous_variables character vector of endogenous variables
#' @param eq_var_matrix contemporaneous-endogenous table, one row per equation
#'   with the equation ids as row names
#' @param decomposition logical. FALSE to put everything in a single block.
#' @param verbose report empty blocks
#' @return list with the three block flags, their endogenous variables and
#'   their equation ids
#' @keywords internal
decomposing_model <- function(endogenous_variables, eq_var_matrix,
                              decomposition = FALSE, verbose = TRUE) {

  note <- function(...) if (verbose) cat(...)

  prologue <- FALSE; heart <- FALSE; epilogue <- FALSE
  prologue_endo <- character(0); heart_endo <- character(0); epilogue_endo <- character(0)
  prologue_equations <- character(0); heart_equations <- character(0)
  epilogue_equations <- character(0)

  remaining_endo   <- endogenous_variables
  equations_vector <- rownames(eq_var_matrix)

  if (!decomposition) {
    heart <- TRUE
    heart_endo <- remaining_endo
    heart_equations <- equations_vector

  } else {

    ## The table is turned once into two index structures -- each equation's
    ## variables, and each variable's equations -- and the peeling below keeps
    ## running counts up to date instead of recounting. tresthor's version,
    ## kept here until 2026-10, rebuilt `table(unlist(eqns_mat))` and re-ran
    ## `apply()` over the whole table on every pass: 22 s on ThreeME 29x33,
    ## for a result that takes a fraction of a second this way. The passes are
    ## the same, in the same order, so the blocks are the same.
    m <- as.matrix(eq_var_matrix)
    n_eq <- nrow(m)
    eq_vars <- lapply(seq_len(n_eq), function(i) { v <- m[i, ]; unname(v[!is.na(v)]) })
    lens <- lengths(eq_vars)
    var_names <- unique(unlist(eq_vars))
    eq_vid <- relist(match(unlist(eq_vars), var_names), eq_vars)   # variables as integers
    var_eqs <- split(rep(seq_len(n_eq), lens), unlist(eq_vid))     # variable -> equations
    var_eqs <- var_eqs[as.character(seq_along(var_names))]

    alive  <- lens > 0L                 # equations still in play
    is_det <- logical(length(var_names))  # variables determined by the prologue
    open   <- lens                      # undetermined variables per equation

    ## ---- prologue -------------------------------------------------------
    ## Repeatedly take the equations with exactly one undetermined endogenous
    ## variable, mark that variable as determined, and blank it out of the
    ## remaining equations.
    repeat {
      take <- which(alive & open == 1L)
      if (length(take) == 0L) break

      newly <- unique(vapply(take, function(e) {
        v <- eq_vid[[e]]; v[!is_det[v]][1L]
      }, integer(1)))

      prologue_equations <- c(prologue_equations, equations_vector[take])
      prologue_endo <- c(prologue_endo, var_names[newly])

      is_det[newly] <- TRUE
      for (v in newly) {
        e <- var_eqs[[v]]
        open[e] <- open[e] - 1L
      }
      ## an equation with nothing left undetermined is out, whether it was
      ## taken above or merely had all its variables determined by others
      alive <- alive & open > 0L
      if (!any(alive)) break
    }
    if (length(prologue_endo) > 0L) prologue <- TRUE else note("Prologue block is empty.\n")

    remaining_endo <- setdiff(remaining_endo, prologue_endo)

    ## ---- epilogue -------------------------------------------------------
    ## An endogenous variable occurring in only one remaining equation is
    ## determined by that equation alone once everything else is known, so
    ## both can be peeled off the end. Repeat until no such variable is left.
    if (length(remaining_endo) > 0L && any(alive)) {
      ## occurrences of each undetermined variable among the remaining equations
      occ <- integer(length(var_names))
      for (e in which(alive)) {
        v <- eq_vid[[e]]; v <- v[!is_det[v]]
        occ[v] <- occ[v] + 1L
      }

      repeat {
        once <- which(occ == 1L)
        if (length(once) == 0L) break
        epilogue_endo <- c(epilogue_endo, var_names[once])

        holds_one <- unique(unlist(lapply(once, function(v) {
          e <- var_eqs[[v]]; e[alive[e]]
        })))
        epilogue_equations <- c(epilogue_equations, equations_vector[holds_one])

        alive[holds_one] <- FALSE
        for (e in holds_one) {
          v <- eq_vid[[e]]; v <- v[!is_det[v]]
          occ[v] <- occ[v] - 1L
        }
        if (!any(alive)) break
      }
    }

    epilogue_endo <- intersect(epilogue_endo, endogenous_variables)
    if (length(epilogue_endo) > 0L) epilogue <- TRUE else note("Epilogue block is empty.\n")

    remaining_endo <- setdiff(remaining_endo, epilogue_endo)

    ## ---- heart ----------------------------------------------------------
    if (length(remaining_endo) > 0L) {
      heart <- TRUE
      heart_endo <- remaining_endo
      heart_equations <- equations_vector[alive]
    } else {
      note("Heart block is empty.\n")
    }
  }

  list(prologue = prologue, heart = heart, epilogue = epilogue,
       prologue_endo = prologue_endo,
       heart_endo    = heart_endo,
       epilogue_endo = epilogue_endo,
       prologue_equations = sort(prologue_equations),
       heart_equations    = sort(heart_equations),
       epilogue_equations = sort(epilogue_equations))
}
