## Translate the compiled model, build it, and solve it over 30 periods.
##
## End to end from what the compiler emits:
##
##   rl_model.prg + rl_calib.csv  --translate-->  rl_model.txt + database
##                                --thor_model-->  compiled sparse solver
##                                --thor_solve-->  30 periods
##
## Usage:  Rscript tools/run_rl_model.R [n_periods]
##
## The first run compiles the generated C++ (a minute or two at this size);
## later runs hit thortwo's compile cache and start solving almost immediately.

suppressMessages(library(thortwo))
source("tools/prg_to_thor.R")

n_periods <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
if (is.na(n_periods)) n_periods <- 30L

PRG   <- "tests/rl_model.prg"
CALIB <- "tests/rl_calib.csv"
BASE_YEAR <- 2019          # relative year 0 in the calibration csv
WORK  <- file.path(tempdir(), "rl")
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(PRG), file.exists(CALIB))

## ---------------------------------------------------------------------------
## 1. Translate
## ---------------------------------------------------------------------------
cat("\n================ 1. translate ================\n")
model_txt <- file.path(WORK, "rl_model.txt")

t_tr <- system.time(
  tr <- prg_to_thor(PRG, CALIB, base_year = BASE_YEAR, out_file = model_txt)
)[["elapsed"]]

translate_report(tr)
cat("translated in", round(t_tr, 1), "s ->", model_txt, "\n")

if (length(tr$warnings)) {
  stop("The translation reported problems; fix them before solving.", call. = FALSE)
}

data_rl <- tr$data

## ---------------------------------------------------------------------------
## 2. Build and compile
## ---------------------------------------------------------------------------
cat("\n================ 2. build ================\n")
t_build <- system.time(
  m <- thor_model("rl_model", model_txt, backend = "sparse",
                  workdir = file.path(WORK, "src"), verbose = FALSE)
)[["elapsed"]]

m
cat("built and compiled in", round(t_build, 1), "s\n")

## ---------------------------------------------------------------------------
## 3. Solve
## ---------------------------------------------------------------------------
cat("\n================ 3. solve ================\n")

## The model reaches `max_lag` periods back -- several equations here are of
## the form delta(1, log(lag(x, 1))), which reads x at t-2 -- so that many
## complete observations have to sit before the first solved one.
years  <- data_rl$year
need   <- max(1L, as.integer(m@meta$max_lag))
first  <- need + 1L
cat("the model reaches", need, "period(s) back, so solving starts at row", first, "\n")

first_period <- years[first]
last_period  <- years[min(length(years), first + n_periods - 1L)]
cat("solving", first_period, "to", last_period,
    sprintf("(%d periods)\n", last_period - first_period + 1L))

res <- thor_solve(m, from = first_period, to = last_period,
                  data = data_rl, index_time = "year",
                  verbose = FALSE, diagnostics = TRUE)

cat("solved in", round(attr(res, "elapsed"), 3), "s,",
    sum(attr(res, "iterations")), "Newton iterations\n")
cat("worst scaled step:", signif(max(attr(res, "convergence")), 3),
    "(converges at 1)\n")

## ---------------------------------------------------------------------------
## 4. Check the answer
## ---------------------------------------------------------------------------
cat("\n================ 4. check ================\n")

## Convergence is not correctness: Newton converges on whatever equations it
## was handed. The residuals are what tell us the generated code means what
## the .prg said.
periods <- first_period:last_period
resid <- model_residuals(m, res, periods = periods, index_time = "year")
cat("max |residual| per block:\n"); print(signif(apply(resid, 2, max), 3))

worst <- max(resid)
cat("\nmax |residual| overall:", signif(worst, 3), "\n")
stopifnot(max(attr(res, "convergence")) <= 1)
stopifnot(worst < 1e-6)

## The strongest check available here: the calibration file already holds a
## consistent baseline, so re-solving over it with its own exogenous path must
## return it. This is what verifies the *translation* -- that the .txt means
## what the .prg said -- independently of whether Newton converged.
rows <- match(periods, res$year)
endo <- m@vars$endo
A <- as.matrix(data_rl[rows, endo]); B <- as.matrix(res[rows, endo])
sc <- pmax(abs(A), abs(B)); rel <- abs(A - B) / sc

med <- stats::median(rel[sc > 0], na.rm = TRUE)
cat("\nagainst the calibration baseline:\n")
cat("  median relative difference          :", signif(med, 3), "\n")
cat("  max relative difference, |value| > 1:", signif(max(rel[sc > 1], na.rm = TRUE), 3), "\n")
stopifnot(med < 1e-10)

## Reported on variables of meaningful magnitude only. A few prices sit over a
## numerically-zero quantity -- pdsd_ctrd is pds over dsd_ctrd ~ -2e-11 -- and
## are genuinely ill-conditioned, so they are compared by median, not maximum.
big <- rel; big[sc <= 1e-3] <- NA
off <- suppressWarnings(sort(apply(big, 2L, max, na.rm = TRUE), decreasing = TRUE))
off <- off[is.finite(off)]
cat("  worst variables (|value| > 1e-3)    :",
    paste(sprintf("%s %.1g", names(off)[1:3], off[1:3]), collapse = ", "), "\n")

## A few headline series, to show the solution is not merely convergent but
## populated with plausible numbers.
show_vars <- intersect(c("gdp", "p", "u", "unr", "ems_co2", "nos"), names(res))
if (length(show_vars)) {
  cat("\nheadline series:\n")
  idx <- match(seq(first_period, last_period, length.out = min(6, length(periods))) |>
                 round(), res$year)
  print(res[unique(stats::na.omit(idx)), c("year", show_vars)], row.names = FALSE)
}

cat("\nPASS\n")
