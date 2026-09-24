## The compile cache, and the two axes the build spec asks to be independently
## switchable: caching on/off, and persistence.
##
## The point of D5 is that the cache must work *by default*. tresthor's cache
## was defeated by the usual call pattern: Rcpp keys its cache on the source
## path, and every caller passed a fresh tempfile(), so the same model
## recompiled from scratch every run. This checks the fix from the outside --
## a second build of the same model must be fast, without anyone having to
## know about `workdir`.
##
## Usage:  Rscript tests/test_cache.R

source("tests/helper.R")

model_file <- "inst/models/opale.txt"
cache <- new_tmpdir("thortwo_cache_")
work  <- new_tmpdir("thortwo_cachework_")
on.exit(unlink(c(cache, work), recursive = TRUE), add = TRUE)

options(thortwo.cache.dir = cache)

## A fresh R session is the case that matters: within one session the model is
## already loaded, so nothing would be recompiled either way.
build_in_new_session <- function(extra = "") {
  script <- sprintf('
    suppressMessages(pkgload::load_all(".", quiet = TRUE))
    options(thortwo.cache.dir = "%s")
    t <- system.time(m <- thor_model("cachetest", "%s", backend = "sparse",
                                     workdir = "%s", verbose = FALSE))[["elapsed"]]
    %s
    cat(round(t, 2), "\n")', cache, model_file, file.path(work, "src"), extra)
  as.numeric(system2("Rscript", c("-e", shQuote(script)), stdout = TRUE))
}

cat("=== cache on (the default) ===\n")
cold <- build_in_new_session()
cat("       first build  ", cold, " s\n", sep = "")
warm <- build_in_new_session()
cat("       second build ", warm, " s\n", sep = "")

## The spec's target: a second build of the same model compiles in under a
## second. Generation itself costs a few tenths of a second, so the whole
## build is allowed a little more than the compile step alone.
ok("second build hits the cache", warm < 3 && warm < cold / 2,
   sprintf("%.2f s versus %.2f s cold", warm, cold))

cat("\n=== cache off ===\n")
nocache <- function() {
  script <- sprintf('
    suppressMessages(pkgload::load_all(".", quiet = TRUE))
    t <- system.time(m <- thor_model("cacheoff", "%s", backend = "sparse",
                                     workdir = "%s", cache = FALSE,
                                     verbose = FALSE))[["elapsed"]]
    cat(round(t, 2), "\n")', model_file, file.path(work, "off"))
  as.numeric(system2("Rscript", c("-e", shQuote(script)), stdout = TRUE))
}
off1 <- nocache(); off2 <- nocache()
cat("       first build  ", off1, " s\n       second build ", off2, " s\n", sep = "")
ok("cache = FALSE really disables it", off2 > 3,
   sprintf("second build still %.2f s", off2))

cat("\n=== the cache directory is where it says it is ===\n")
files <- list.files(cache, recursive = TRUE)
ok("objects were cached", length(files) > 0, sprintf("%d files under %s", length(files), cache))
ok("clear_model_cache empties it",
   { clear_model_cache(confirm = FALSE); length(list.files(cache, recursive = TRUE)) } == 0)

cat("\n=== build without compiling, then solve ===\n")
## `compile = FALSE` must produce a usable model: the code is there, and the
## first solve compiles it. That is what makes "saving or not saving" free.
m <- thor_model("lazy", model_file, backend = "sparse",
                workdir = file.path(work, "lazy"), compile = FALSE, verbose = FALSE)
ok("uncompiled model carries its code", nchar(m@generated$code) > 1000)

data_opale <- readRDS("inst/Opale/donnees_opale.rds")
coeffs     <- readRDS("inst/Opale/coefficients_opale.rds")
data_opale <- add_coeffs(coeffs, data_opale, pos.coeff.name = 2, pos.coeff.value = 1)
dates <- as.character(data_opale$date); n <- length(dates)
res <- thor_solve(m, dates[n - 9], dates[n], data_opale, "date", verbose = FALSE,
                  diagnostics = TRUE)
check_converged(res, "solving an uncompiled model compiles it")

cat("\n=== a loaded model does not rebuild ===\n")
rds <- file.path(work, "saved.rds")
thor_save(m, rds)
t_load <- system.time(m2 <- thor_load(rds, workdir = file.path(work, "loaded")))[["elapsed"]]
ok("thor_load is cheap when the code is unchanged", t_load < 3,
   sprintf("%.2f s", t_load))
ok("the loaded model has the same generated code", identical(m@generated$hash, m2@generated$hash))

cat("\nPASS\n")
