## The built-model cache.
##
## thor_model() keeps every model it builds under a hash of what determines
## it, and hands the stored one back when asked for the same model again --
## in this session or a later one -- instead of parsing, differentiating and
## generating code a second time. On ThreeME 13x13 that is 0.4 s instead of
## 30 s plus the compile.
##
## A cache like this is only worth having if it is never wrong, so most of
## what is checked here is when it must NOT answer: an equation edited, a
## different backend, `recompile = TRUE`, `cache = FALSE`, a damaged file.
##
## Usage:  Rscript tests/test_model_cache.R

source("tests/helper.R")

cache <- new_tmpdir("thortwo_modelcache_")
work  <- new_tmpdir("thortwo_modelcachework_")
on.exit(unlink(c(cache, work), recursive = TRUE), add = TRUE)
options(thortwo.cache.dir = cache, thortwo.model.cache = TRUE)   # helper.R turns it off

## Runs a build and reports whether it came from the cache.
build <- function(...) {
  msgs <- character(0)
  m <- withCallingHandlers(thor_model(..., verbose = FALSE),
    message = function(c) { msgs <<- c(msgs, conditionMessage(c)); invokeRestart("muffleMessage") })
  list(model = m, hit = any(grepl("already been compiled before", msgs)),
       says_how = any(grepl("recompile = TRUE", msgs)))
}

write_model <- function(file, mpc = "0.6", lhs = "y = ch + i + g") {
  writeLines(c("endogenous :", "ch, i, y", "", "exogenous :", "g", "", "coefficients :", "",
               "equations :", lhs, paste0("consumption : ch = ", mpc, " * lag(y, 1)"),
               "i = 0.25 * y"), file)
  file
}
f <- write_model(file.path(work, "mini.txt"))
wd <- function(x) file.path(work, x)

cat("=== hits ===\n")
b1 <- build("mc", f, backend = "sparse-r", workdir = wd("a"))
ok("first build is a real build", !b1$hit)
b2 <- build("mc", f, backend = "sparse-r", workdir = wd("a"))
ok("second build comes from the cache, and says how to recompile", b2$hit && b2$says_how)
ok("the cached model is the same model",
   identical(b1$model@generated$hash, b2$model@generated$hash) &&
     identical(b1$model@equations, b2$model@equations))

## the same model, written differently
g <- file.path(work, "mini_reformatted.txt")
writeLines(c("ENDOGENOUS :", "y,   ch,i", "exogenous :", "G", "coefficients :", "", "equations :",
             "Y   =   CH + I + G", "consumption:ch=0.6*lag(y,1)", "i = 0.25 * y"), g)
ok("spacing, case and the order of the variable lists do not matter",
   build("mc", g, backend = "sparse-r", workdir = wd("a"))$hit)

d <- data.frame(year = 2015:2020, g = 20, y = 100, ch = 55, i = 25)
r1 <- thor_solve(b1$model, 2016, 2020, d, "year", verbose = FALSE)
r2 <- thor_solve(b2$model, 2016, 2020, d, "year", verbose = FALSE)
ok("a cached model solves like a built one", identical(r1, r2))

cat("\n=== misses ===\n")
ok("an edited equation",
   !build("mc", write_model(file.path(work, "edited.txt"), mpc = "0.61"),
          backend = "sparse-r", workdir = wd("b"))$hit)
ok("another backend", !build("mc", f, backend = "sparse", workdir = wd("c"))$hit)
ok("another decomposition",
   !build("mc", f, backend = "sparse-r", decompose = FALSE, workdir = wd("d"))$hit)
ok("another name", !build("mc2", f, backend = "sparse-r", workdir = wd("e"))$hit)

cat("\n=== switches ===\n")
ok("recompile = TRUE rebuilds",
   !build("mc", f, backend = "sparse-r", workdir = wd("a"), recompile = TRUE)$hit)
ok("... and the result is cached again", build("mc", f, backend = "sparse-r", workdir = wd("a"))$hit)
ok("cache = FALSE does not look",
   !build("mc", f, backend = "sparse-r", workdir = wd("a"), cache = FALSE)$hit)
options(thortwo.model.cache = FALSE)
ok("options(thortwo.model.cache = FALSE) does not look",
   !build("mc", f, backend = "sparse-r", workdir = wd("a"))$hit)
options(thortwo.model.cache = TRUE)

cat("\n=== a damaged cache file is ignored ===\n")
stored <- list.files(file.path(list.files(cache, full.names = TRUE)[1], "models"), full.names = TRUE)
ok("models are stored where the cache says", length(stored) >= 4, sprintf("%d files", length(stored)))
for (s in stored) writeLines("not an rds file", s)
b3 <- build("mc", f, backend = "sparse-r", workdir = wd("a"))
ok("the model is rebuilt instead", !b3$hit && identical(b3$model@generated$hash, b1$model@generated$hash))

cat("\n=== across sessions, with compiled code ===\n")
in_new_session <- function() {
  script <- sprintf('
    suppressMessages(pkgload::load_all(".", quiet = TRUE))
    options(thortwo.timings = FALSE)    # stdout is parsed below
    options(thortwo.cache.dir = "%s")
    hit <- FALSE
    t <- system.time(m <- withCallingHandlers(
      thor_model("mc_opale", "inst/models/opale.txt", workdir = "%s", verbose = FALSE),
      message = function(c) { hit <<- TRUE; invokeRestart("muffleMessage") }))[["elapsed"]]
    cat(round(t, 2), hit, "\n")', cache, wd("opale"))
  out <- strsplit(trimws(system2("Rscript", c("-e", shQuote(script)), stdout = TRUE)), " ")[[1]]
  list(seconds = as.numeric(out[1]), hit = as.logical(out[2]))
}
s1 <- in_new_session(); s2 <- in_new_session()
cat("       first session  ", s1$seconds, " s\n       second session ", s2$seconds, " s\n", sep = "")
ok("a new session finds the model built by an earlier one", !s1$hit && s2$hit)
ok("and gets it quickly", s2$seconds < 2 && s2$seconds < s1$seconds / 3,
   sprintf("%.2f s versus %.2f s", s2$seconds, s1$seconds))

ok("clear_model_cache empties it",
   { clear_model_cache(confirm = FALSE); length(list.files(cache, recursive = TRUE)) } == 0)

cat("\nPASS\n")
