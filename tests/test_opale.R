
## Opale, on all three backends, over 40 quarters.
##
## Opale is the quarterly French macro model shipped with the package: 496
## equations, decomposing into 201 / 98 / 197. It is much smaller than ThreeME
## but structurally harder -- quarterly data, heavy use of lags and of the
## `trim` seasonal index, and a large epilogue of accounting identities -- so
## it exercises the jacobians and the code generators on a shape of model

## ThreeME does not reach.
##
## Usage:  Rscript tests/test_opale.R [sparse|dense-cpp|sparse-r|all]

source("tests/helper.R")

which_backends <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(which_backends) || which_backends == "all") {
  which_backends <- c("sparse", "dense-cpp", "sparse-r")
}

model_file <- "inst/models/opale.txt"
stopifnot(file.exists(model_file))
work <- new_tmpdir("thortwo_opale_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

## ---- data ---------------------------------------------------------------
## The estimated coefficients live outside the database and have to be
## broadcast onto it as constant columns before solving.
data_opale <- readRDS("inst/Opale/donnees_opale.rds")
coeffs     <- readRDS("inst/Opale/coefficients_opale.rds")
data_opale <- add_coeffs(coeffs, data_opale, pos.coeff.name = 2, pos.coeff.value = 1)

dates <- as.character(data_opale$date)
n <- length(dates)
periods <- dates[seq(n - 39, n)]
first_period <- periods[1]; last_period <- periods[length(periods)]
rows <- match(periods, dates)

solutions <- list()

for (backend in which_backends) {

  cat("\n=== ", backend, " ===\n", sep = "")

  t_build <- system.time(
    m <- thor_model("opale", model_file, backend = backend,
                    workdir = file.path(work, backend), verbose = FALSE)
  )[["elapsed"]]

  ## block sizes are a property of the model, not of the backend
  ok("decomposition 201/98/197",
     identical(unname(vapply(m@blocks, function(b) length(b$equations), integer(1))),
               c(201L, 98L, 197L)))

  res <- thor_solve(m, first_period, last_period, data_opale, "date",
                    verbose = FALSE, diagnostics = TRUE)

  cat("       build ", round(t_build, 1), " s   solve ",
      round(attr(res, "elapsed"), 3), " s for 40 periods, ",
      sum(attr(res, "iterations")), " Newton iterations\n", sep = "")

  check_converged(res)
  check_residuals(m, res, periods, "date")
  check_history(endo_matrix(data_opale, m, rows), endo_matrix(res, m, rows))

  solutions[[backend]] <- endo_matrix(res, m, rows)

  ## ---- save -> load -> solve ---------------------------------------------
  ## A saved model carries its generated source, so loading it needs neither
  ## the original .txt nor the directory it was built in.
  rds <- file.path(work, paste0(backend, ".rds"))
  thor_save(m, rds)
  unlink(file.path(work, backend), recursive = TRUE)   # the build dir is gone
  m2 <- thor_load(rds, workdir = file.path(work, "reloaded", backend))
  res2 <- thor_solve(m2, first_period, last_period, data_opale, "date", verbose = FALSE)
  ok("save -> load -> solve",
     identical(endo_matrix(res, m, rows), endo_matrix(res2, m2, rows)),
     "identical to the pre-save solution")
}

## ---- backends agree -----------------------------------------------------
if (length(solutions) > 1L) {
  cat("\n=== agreement ===\n")
  nms <- names(solutions)
  for (i in seq_along(nms)[-1]) {
    check_agree(solutions[[nms[1]]], solutions[[nms[i]]],
                label = paste(nms[1], "vs", nms[i]))
  }
}

cat("\nPASS\n")
