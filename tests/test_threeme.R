## ThreeME, on the compiled backends.
##
## ThreeME is the large end of the range: 1729 equations at the 4x4
## classification, 3628 at 8x8. It is what the sparse path exists for -- at
## 8x8 the heart block alone is 2517 x 2517, whose dense character jacobian is
## 6.3M cells holding 9k real entries.
##
## The R backend is not run here. It is not broken; it is simply the wrong
## tool at this size, and timing it out is not a useful test. Use
## `Rscript tests/test_opale.R dense-r` for its coverage.
##
## Usage:  Rscript tests/test_threeme.R [4x4|8x8] [sparse|dense-cpp|all]

source("tests/helper.R")

args <- commandArgs(trailingOnly = TRUE)
classification <- if (is.na(args[1])) "4x4" else args[1]
stopifnot(classification %in% c("4x4", "8x8"))

which_backends <- if (is.na(args[2]) || args[2] == "all") c("sparse", "dense-cpp") else args[2]

model_file <- file.path("tests", paste0("threeme_", classification, "_thor.txt"))
data_file  <- file.path("tests", paste0("data3me_", classification, ".rds"))
stopifnot(file.exists(model_file), file.exists(data_file))

work <- new_tmpdir("thortwo_threeme_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

data_3me <- readRDS(data_file)
first_period <- 2016
last_period  <- 2050
periods <- first_period:last_period
rows <- match(as.character(periods), as.character(data_3me$year))

solutions <- list()

for (backend in which_backends) {

  cat("\n=== ", classification, " / ", backend, " ===\n", sep = "")

  t_build <- system.time(
    m <- thor_model(paste0("threeme", classification), model_file, backend = backend,
                    workdir = file.path(work, backend), verbose = FALSE)
  )[["elapsed"]]

  cat("       ", nrow(m@equations), " equations, blocks ",
      paste(vapply(m@blocks, function(b) length(b$equations), integer(1)), collapse = "/"),
      "\n", sep = "")

  res <- thor_solve(m, first_period, last_period, data_3me, "year",
                    verbose = FALSE, diagnostics = TRUE)

  cat("       build ", round(t_build, 1), " s   solve ",
      round(attr(res, "elapsed"), 3), " s for ", length(periods), " periods, ",
      sum(attr(res, "iterations")), " Newton iterations\n", sep = "")

  check_converged(res)
  check_residuals(m, res, periods, "year")

  solutions[[backend]] <- as.matrix(res[rows, m@vars$endo, drop = FALSE])

  rds <- file.path(work, paste0(backend, ".rds"))
  thor_save(m, rds)
  m2 <- thor_load(rds, workdir = file.path(work, "reloaded", backend))
  res2 <- thor_solve(m2, first_period, last_period, data_3me, "year", verbose = FALSE)
  ok("save -> load -> solve",
     identical(as.matrix(res[rows, m@vars$endo]), as.matrix(res2[rows, m2@vars$endo])),
     "identical to the pre-save solution")
}

if (length(solutions) > 1L) {
  cat("\n=== agreement ===\n")
  ## Medians, not maxima: a handful of ThreeME variables (`pds_*`, a price
  ## over a near-zero stock change) are genuinely ill-conditioned, and a
  ## maximum reports those rather than the solvers.
  check_agree(solutions[[1]], solutions[[2]],
              label = paste(names(solutions)[1], "vs", names(solutions)[2]))
}

cat("\nPASS\n")
