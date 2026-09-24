## ThreeME FRA 29x33, end to end from what the EViews compiler emits.
##
##   model_29_33.prg + calib_29_33.csv  --translate-->  .txt model + database
##                                      --thor_model-->  compiled sparse solver
##                                      --thor_solve-->  solution
##
## This is the largest model in the suite: 26,439 equations, about four times
## ThreeME 13x13. Only the sparse backend is run; a dense jacobian at this size
## is not an option. Every step is timed, and the timings are summarised at
## the end.
##
## The first run compiles the generated C++, which is the slowest step at this
## size (13x13 compiles in ~20 s). Later runs of an unchanged model hit thortwo's
## compile cache: pass `rebuild` to compile from scratch and time it again
## (the cache itself is left alone, so other models are not affected).
##
## The translation step is ermeeth2::prg_to_thor(), so ermeeth2 has to be
## installed; thortwo itself does not depend on it.
##
## Usage:  Rscript tests/test_threeme_29x33.R [n_periods|all] [rebuild]
##   e.g.  Rscript tests/test_threeme_29x33.R            # all periods
##         Rscript tests/test_threeme_29x33.R 5          # first 5 periods
##         Rscript tests/test_threeme_29x33.R all rebuild
##
## Or source it from the repository root in an interactive session; it then
## runs with the defaults.

source("tests/helper.R")

if (!requireNamespace("ermeeth2", quietly = TRUE)) {
  stop("This test needs ermeeth2 for prg_to_thor().", call. = FALSE)
}

args      <- commandArgs(trailingOnly = TRUE)
n_periods <- suppressWarnings(as.integer(args[1]))    # NA = all periods
rebuild   <- "rebuild" %in% args

PRG       <- "tests/model_29_33.prg"
CALIB     <- "tests/calib_29_33.csv"
BASE_YEAR <- 2019          # relative year 0 in the calibration csv
NAME      <- "threeme29x33"
stopifnot(file.exists(PRG), file.exists(CALIB))

work <- new_tmpdir("thortwo_threeme29x33_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

## ---- time keeping --------------------------------------------------------
## timed("label", expr) evaluates expr, prints how long it took, records it for
## the summary at the end, and returns the value of expr.
timings <- data.frame(step = character(0), seconds = numeric(0))
timed <- function(label, expr) {
  cat("\n---- ", label, " ----\n", sep = "")
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  dt <- proc.time()[["elapsed"]] - t0
  timings[nrow(timings) + 1L, ] <<- list(label, dt)
  cat(sprintf("     %s: %s\n", label, fmt_time(dt)))
  invisible(value)
}
fmt_time <- function(s) {
  if (s < 60) sprintf("%.2f s", s) else sprintf("%d min %04.1f s", s %/% 60, s %% 60)
}
t_total <- proc.time()[["elapsed"]]

cat("ThreeME 29x33  |  ", format(Sys.time(), "%Y-%m-%d %H:%M"), "  |  thortwo ",
    as.character(utils::packageVersion("thortwo")), ", ermeeth2 ",
    as.character(utils::packageVersion("ermeeth2")), "\n", sep = "")

## ---------------------------------------------------------------------------
## 1. Translate the .prg and the calibration
## ---------------------------------------------------------------------------
model_txt <- file.path(work, paste0(NAME, ".txt"))

tr <- timed("1. translate (prg_to_thor)",
  ermeeth2::prg_to_thor(PRG, CALIB, base.year = BASE_YEAR,
                        out_file = model_txt, verbose = FALSE)
)
ermeeth2::translate_report(tr)
if (length(tr$warnings)) {
  stop("The translation reported problems; fix them before solving.", call. = FALSE)
}
data_3me <- tr$data
names(data_3me) <- tolower(names(data_3me))

## A sector with no capital at all leaves its price of capital undetermined
## (prices.mdl: PK[s] * F[K, s] = ...): this is what stopped ThreeME 13x13,
## see ThreeME_V4/FROM_THORTWO.md. Flag it up front rather than after the build.
fk <- grep("^f_k_s[a-z0-9]+$", names(data_3me), value = TRUE)
no_capital <- fk[vapply(fk, function(v) all(data_3me[[v]] == 0, na.rm = TRUE), logical(1))]
if (length(no_capital)) {
  cat("\n  WARNING  zero capital in every period: ", paste(no_capital, collapse = ", "),
      "\n           the matching pk_* is undetermined and the solve will likely fail.\n",
      sep = "")
}

## ---------------------------------------------------------------------------
## 2. Build and compile
## ---------------------------------------------------------------------------
## No `workdir`: the default is stable across runs, which is what lets the
## compile cache recognise an unchanged model. `rebuild` compiles without the
## cache (cache = FALSE), so the full compilation is timed.
m <- timed(if (rebuild) "2. build + compile, no cache" else "2. build + compile (thor_model)",
  thor_model(NAME, model_txt, backend = "sparse",
             cache = if (rebuild) FALSE else NULL, verbose = TRUE)
)

blocks <- vapply(m@blocks, function(b) length(b$equations), integer(1))
cat("\n     ", nrow(m@equations), " equations, blocks ",
    paste(names(blocks), blocks, sep = " ", collapse = " / "),
    ", max lag ", m@meta$max_lag, "\n", sep = "")

## ---------------------------------------------------------------------------
## 3. Solve
## ---------------------------------------------------------------------------
## The model reaches `max_lag` periods back, so that many complete
## observations have to sit before the first solved one.
years <- data_3me$year
first_period <- years[1L + max(1L, as.integer(m@meta$max_lag))]
last_period  <- if (is.na(n_periods)) max(years) else
  min(max(years), first_period + n_periods - 1L)
periods <- first_period:last_period
rows    <- match(periods, years)
cat("\n     solving ", first_period, " to ", last_period,
    " (", length(periods), " periods)\n", sep = "")

res <- timed("3. solve (thor_solve)",
  tryCatch(
    thor_solve(m, first_period, last_period, data_3me, "year",
               verbose = TRUE, diagnostics = TRUE),
    error = function(e) {
      ## Name the equations the data already misses at the first period: a
      ## non-convergence there is usually the calibration, not the solver.
      cat("\n  Solve failed: ", conditionMessage(e), "\n\n", sep = "")
      cc <- calibration_check(m, data_3me, first_period, "year")
      cat("  Equations the calibration misses at ", first_period, ": ", NROW(cc), "\n", sep = "")
      if (NROW(cc)) {
        cc$formula <- substr(cc$formula, 1, 90)
        print(utils::head(cc, 10), right = FALSE)
      }
      stop("Solve failed; see above.", call. = FALSE)
    })
)
cat(sprintf("     %d Newton iterations, %.3f s per period\n",
            sum(attr(res, "iterations")), attr(res, "elapsed") / length(periods)))

## ---------------------------------------------------------------------------
## 4. Check the answer
## ---------------------------------------------------------------------------
check_converged(res)

## Convergence is not correctness: the residuals are what say the generated
## code means what the .prg said.
timed("4. residuals (model_residuals)",
  check_residuals(m, res, periods, "year")
)

## Solving over the calibration with its own exogenous path should return it.
## Only the median is tested: the ThreeME calibrations do not satisfy every
## equation exactly (13x13 missed some by up to 2.8), so the worst variable
## measures the calibration rather than the solver. It is reported, not tested.
A <- as.matrix(data_3me[rows, m@vars$endo]); B <- as.matrix(res[rows, m@vars$endo])
d <- rel_diff(A, B)
med   <- stats::median(d$rel[d$sc > 0], na.rm = TRUE)
worst <- max(d$rel[d$sc > 1], na.rm = TRUE)
ok("reproduces the calibration (median)", med < 1e-8,
   sprintf("median %.3g; max %.3g for |value| > 1, not tested", med, worst))
big <- d$rel; big[d$sc <= 1e-3] <- NA
off <- suppressWarnings(sort(apply(big, 2L, max, na.rm = TRUE), decreasing = TRUE))
off <- off[is.finite(off)]
cat("       worst variables (|value| > 1e-3): ",
    paste(sprintf("%s %.1g", names(off)[1:3], off[1:3]), collapse = ", "), "\n", sep = "")

## ---------------------------------------------------------------------------
## 5. Save -> load -> solve
## ---------------------------------------------------------------------------
rds <- file.path(work, paste0(NAME, ".rds"))
timed("5a. save (thor_save)", thor_save(m, rds))
cat(sprintf("     %.1f MB on disk\n", file.size(rds) / 1024^2))
## Near zero here: the same code is already loaded in this session. In a new
## session it costs a cache lookup, or a full compile on a new machine.
m2 <- timed("5b. load (thor_load)", thor_load(rds))
res2 <- timed("5c. solve the loaded model",
  thor_solve(m2, first_period, last_period, data_3me, "year", verbose = FALSE)
)
ok("save -> load -> solve", identical(B, as.matrix(res2[rows, m2@vars$endo])),
   "identical to the pre-save solution")

## ---------------------------------------------------------------------------
## Timings
## ---------------------------------------------------------------------------
timings[nrow(timings) + 1L, ] <- list("total", proc.time()[["elapsed"]] - t_total)
cat("\n==== timings ====\n")
for (i in seq_len(nrow(timings))) {
  cat(sprintf("  %-34s %14s\n", timings$step[i], fmt_time(timings$seconds[i])))
}

cat("\nPASS\n")
