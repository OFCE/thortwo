## Block decomposition.
##
## Ported from tresthor's 1_3_decomposition_algo.R, which works; the only
## changes are the dplyr verbs, replaced by base subsetting.
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
  eqns_mat <- as.data.frame(eq_var_matrix, stringsAsFactors = FALSE)

  if (!decomposition) {
    heart <- TRUE
    heart_endo <- remaining_endo
    heart_equations <- equations_vector

  } else {

    ## ---- prologue -------------------------------------------------------
    ## Repeatedly take the equations with exactly one undetermined endogenous
    ## variable, mark that variable as determined, and blank it out of the
    ## remaining equations.
    n_endo <- function(df) apply(df, 1L, function(y) sum(!is.na(y)))

    counts <- n_endo(eqns_mat)
    while (length(counts) && min(counts) == 1L) {
      unique_endo <- eqns_mat[counts == 1L, , drop = FALSE]

      prologue_equations <- c(prologue_equations, rownames(unique_endo))
      prologue_endo <- c(prologue_endo,
                         stats::na.omit(unique(as.vector(as.matrix(unique_endo)))))

      m <- as.matrix(eqns_mat)
      m[m %in% prologue_endo] <- NA
      eqns_mat <- as.data.frame(m, stringsAsFactors = FALSE)

      counts <- n_endo(eqns_mat)
      eqns_mat <- eqns_mat[counts > 0L, , drop = FALSE]
      counts <- counts[counts > 0L]
      if (length(counts) == 0L) break
    }
    if (length(prologue_endo) > 0L) prologue <- TRUE else note("Prologue block is empty.\n")

    remaining_endo <- setdiff(remaining_endo, prologue_endo)

    ## ---- epilogue -------------------------------------------------------
    ## An endogenous variable occurring in only one remaining equation is
    ## determined by that equation alone once everything else is known, so
    ## both can be peeled off the end. Repeat until no such variable is left.
    if (length(remaining_endo) > 0L && nrow(eqns_mat) > 0L) {
      occurrences <- table(unlist(eqns_mat))

      while (length(occurrences) && min(occurrences) == 1L) {
        epilogue_endo <- c(epilogue_endo, names(occurrences)[occurrences == 1L])

        m <- as.matrix(eqns_mat)
        m[!m %in% epilogue_endo] <- NA
        holds_one <- rowSums(!is.na(m)) > 0L
        epilogue_equations <- c(epilogue_equations, rownames(eqns_mat)[holds_one])

        eqns_mat <- eqns_mat[!holds_one, , drop = FALSE]
        if (nrow(eqns_mat) == 0L) break
        occurrences <- table(unlist(eqns_mat))
      }
    }

    epilogue_endo <- intersect(epilogue_endo, endogenous_variables)
    if (length(epilogue_endo) > 0L) epilogue <- TRUE else note("Epilogue block is empty.\n")

    remaining_endo <- setdiff(remaining_endo, epilogue_endo)

    ## ---- heart ----------------------------------------------------------
    if (length(remaining_endo) > 0L) {
      heart <- TRUE
      heart_endo <- remaining_endo
      heart_equations <- rownames(eqns_mat)
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
