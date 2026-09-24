## thor_model(sequential = TRUE): the prologue and the epilogue solved one
## equation at a time instead of by Newton.
##
## The option changes how two blocks are solved, not what the model is, so
## the test is that nothing else changes:
##
##   1. the solution is Newton's, on every backend, on Opale and on ThreeME;
##   2. each equation is classified correctly -- evaluated directly when it is
##      written `x = f(...)`, solved by a one-variable Newton otherwise -- and
##      the equations are taken in an order that works;
##   3. a failure names the equation and the variable, which is the one thing
##      this mode can say that Newton on a whole block cannot;
##   4. with the option off, the generated code contains none of it.
##
## Usage:  Rscript tests/test_sequential.R

source("tests/helper.R")

work <- new_tmpdir("thortwo_sequential_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)
wd <- function(x) file.path(work, x)

## ---- 1. a small model with every case in it ------------------------------------
## Prologue: a chain, written out of order, with one direct equation, one that
##           needs the scalar Newton because of the log, and one because the
##           variable is on both sides.
## Heart:    y and c determine each other.
## Epilogue: direct, and implicit through a product.
cat("=== a hand-made model ===\n")
f <- file.path(work, "small.txt")
writeLines(c(
  "endogenous :", "a, b, d, y, c, s, pk", "exogenous :", "x, fk", "coefficients :", "",
  "equations :",
  "d = b + a",                    # prologue, needs a and b first
  "log(b) = 0.5 * log(a) + 0.1",  # prologue, scalar: left-hand side is not the variable
  "a = 2 * x + 0.5 * lag(a, 1)",  # prologue, direct
  "y = c + d",                    # heart
  "c = 0.6 * y + 1",              # heart
  "s = y - c",                    # epilogue, direct
  "capital : pk * fk = 0.9 * lag(pk, 1) * lag(fk, 1) + s"), f)   # epilogue, scalar

years <- 2000:2010
d <- data.frame(year = years, x = seq(1, 2, length.out = 11), fk = 10,
                a = 3, b = 2, d = 5, y = 15, c = 10, s = 5, pk = 1)

n0 <- thor_model("sq_small", f, backend = "sparse-r", workdir = wd("n"), verbose = FALSE)
s1 <- thor_model("sq_small", f, backend = "sparse-r", sequential = TRUE,
                 workdir = wd("s"), verbose = FALSE)

plan <- s1@meta$sequential
ok("the prologue is taken in dependency order",
   identical(plan$prologue$variable, c("a", "b", "d")))
ok("direct and scalar equations are told apart",
   identical(plan$prologue$direct, c(TRUE, FALSE, TRUE)) &&
     identical(plan$epilogue$direct[match(c("s", "pk"), plan$epilogue$variable)], c(TRUE, FALSE)))
ok("the heart is left to Newton", is.null(plan$heart) && length(s1@jacobian$heart$i) > 0)
ok("no jacobian is built for the sequential blocks",
   length(s1@jacobian$prologue$i) == 0 && length(s1@jacobian$epilogue$i) == 0)

rn <- thor_solve(n0, 2001, 2010, d, "year", verbose = FALSE)
rs <- thor_solve(s1, 2001, 2010, d, "year", verbose = FALSE, diagnostics = TRUE)
endo <- n0@vars$endo
ok("same solution as Newton", max(abs(as.matrix(rn[, endo]) - as.matrix(rs[, endo]))) < 1e-8,
   sprintf("largest difference %.2g", max(abs(as.matrix(rn[, endo]) - as.matrix(rs[, endo])))))
check_converged(rs)
check_residuals(s1, rs, 2001:2010, "year")

## ---- 2. a failure names the equation and the variable ----------------------------
d0 <- d; d0$fk <- 0          # pk * 0 = ... : pk cannot be determined
msg <- tryCatch({ thor_solve(s1, 2001, 2010, d0, "year", verbose = FALSE); "" },
                error = conditionMessage)
ok("an undetermined variable is reported by name",
   grepl("capital", msg, fixed = TRUE) && grepl("'pk'", msg, fixed = TRUE), msg)

## ---- 3. off means off -------------------------------------------------------------
cat("\n=== the option is opt-in ===\n")
ok("without it, the generated code has no sequential part",
   !grepl("_seq", n0@generated$code, fixed = TRUE))
c0 <- thor_model("sq_small", f, backend = "sparse", compile = FALSE, workdir = wd("c"), verbose = FALSE)
ok("... in C++ either", !grepl("SeqCtl|THOR_SEQ|_seq", c0@generated$code))
u <- thor_model("sq_small", f, backend = "sparse-r", sequential = TRUE, decompose = FALSE,
                workdir = wd("u"), verbose = FALSE)
ok("without the decomposition there is nothing to solve sequentially",
   length(u@meta$sequential) == 0 && !grepl("_seq", u@generated$code, fixed = TRUE))

## ---- 4. Opale, every backend ---------------------------------------------------------
data_opale <- readRDS("inst/Opale/donnees_opale.rds")
coeffs     <- readRDS("inst/Opale/coefficients_opale.rds")
data_opale <- add_coeffs(coeffs, data_opale, pos.coeff.name = 2, pos.coeff.value = 1)
dates <- as.character(data_opale$date); n <- length(dates)
from <- dates[n - 39]; to <- dates[n]; rows <- (n - 39):n

for (backend in c("sparse", "dense-cpp", "sparse-r")) {
  cat("\n=== Opale / ", backend, " ===\n", sep = "")
  mn <- thor_model("sq_opale_n", "inst/models/opale.txt", backend = backend,
                   workdir = wd(paste0("on", backend)), verbose = FALSE)
  ms <- thor_model("sq_opale_s", "inst/models/opale.txt", backend = backend, sequential = TRUE,
                   workdir = wd(paste0("os", backend)), verbose = FALSE)
  res_n <- thor_solve(mn, from, to, data_opale, "date", verbose = FALSE)
  res_s <- thor_solve(ms, from, to, data_opale, "date", verbose = FALSE, diagnostics = TRUE)
  cat("       ", paste(sprintf("%s: %d direct, %d scalar", names(ms@meta$sequential),
                               vapply(ms@meta$sequential, function(p) sum(p$direct), 1L),
                               vapply(ms@meta$sequential, function(p) sum(!p$direct), 1L)),
                       collapse = "; "), "\n", sep = "")
  check_converged(res_s)
  check_residuals(ms, res_s, dates[rows], "date")
  check_agree(endo_matrix(res_n, mn, rows), endo_matrix(res_s, ms, rows),
              label = "sequential and Newton agree")

  if (backend == "sparse") {
    rds <- file.path(work, "seq.rds")
    thor_save(ms, rds)
    ml <- thor_load(rds, workdir = wd("loaded"))
    res_l <- thor_solve(ml, from, to, data_opale, "date", verbose = FALSE)
    ok("a saved sequential model reloads and solves identically",
       identical(endo_matrix(res_s, ms, rows), endo_matrix(res_l, ml, rows)))
  }
}

## ---- 5. ThreeME 4x4 -------------------------------------------------------------------
cat("\n=== ThreeME 4x4 / sparse ===\n")
data_3me <- readRDS("inst/ThreeME/data3me_4x4.rds")
rows3 <- match(2016:2050, data_3me$year)
tn <- thor_model("sq_3me_n", "inst/ThreeME/threeme_4x4_thor.txt", workdir = wd("tn"), verbose = FALSE)
ts <- thor_model("sq_3me_s", "inst/ThreeME/threeme_4x4_thor.txt", sequential = TRUE,
                 workdir = wd("ts"), verbose = FALSE)
r3n <- thor_solve(tn, 2016, 2050, data_3me, "year", verbose = FALSE)
r3s <- thor_solve(ts, 2016, 2050, data_3me, "year", verbose = FALSE, diagnostics = TRUE)
check_converged(r3s)
check_residuals(ts, r3s, 2016:2050, "year")
check_agree(as.matrix(r3n[rows3, tn@vars$endo]), as.matrix(r3s[rows3, ts@vars$endo]),
            label = "sequential and Newton agree")

cat("\nPASS\n")
