## Parser regression tests.
##
## `export_model(thor_model(f)) |> thor_model()` must produce an identical
## equations table. Cheap, and it pins the parser: any change to the reader,
## the namer, the LHS/RHS split or the residual form shows up here.
##
## Also checks D6 directly -- the layout of a model file must not depend on
## line numbers.
##
## Usage:  Rscript tests/test_roundtrip.R

source("tests/helper.R")

work <- new_tmpdir("thortwo_roundtrip_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

build <- function(f, ...) {
  thor_model("rt", f, compile = FALSE, workdir = file.path(work, "src"),
             verbose = FALSE, ...)
}

for (f in c("inst/models/opale.txt", "tests/threeme_4x4_thor.txt")) {

  cat("\n=== ", basename(f), " ===\n", sep = "")

  m1 <- build(f)
  out <- file.path(work, "exported.txt")
  export_model(m1, out)
  m2 <- build(out)

  ok("equations table round-trips",
     identical(m1@equations[, c("id", "name", "equation", "LHS", "RHS",
                                "formula", "new_formula", "part")],
               m2@equations[, c("id", "name", "equation", "LHS", "RHS",
                                "formula", "new_formula", "part")]),
     sprintf("%d equations", nrow(m1@equations)))
  ok("variables round-trip", identical(m1@vars, m2@vars))
  ok("blocks round-trip", identical(m1@blocks, m2@blocks))
}

## ---- D6: the layout must not matter -------------------------------------
cat("\n=== layout independence (D6) ===\n")
src <- readLines("inst/models/model_type.txt")
base <- build("inst/models/model_type.txt")

## Same model, but with the separator banners removed and blank lines added
## where tresthor's fixed line numbers would have shifted every section.
shuffled <- file.path(work, "shuffled.txt")
writeLines(c("", "endogenous variables:", "", "endovar1,endovar2", "endovar3,endovar4",
             "", "", "exo :", "exovar1,exovar2,exovar3", "exovar4,exovar5",
             "#########", "", "coefficients", "cf1,cf2,cf3", "", "",
             "equations :", "",
             src[grep("^(delta|equation_var3|equilibrium)", src)]), shuffled)
alt <- build(shuffled)

ok("blank lines and missing banners do not shift the sections",
   identical(base@vars, alt@vars) &&
   identical(base@equations$formula, alt@equations$formula))

## A file whose sections cannot be found must fail loudly rather than
## mis-parse: that is the whole point of the change.
broken <- file.path(work, "broken.txt")
writeLines(c("endogenous variables:", "a,b", "equations:", "a=b", "b=1"), broken)
ok("a missing section is an error",
   inherits(try(read_model_source(broken), silent = TRUE), "try-error"))

## ---- manual input matches file input ------------------------------------
cat("\n=== manual input ===\n")
p <- read_model_source("inst/models/model_type.txt")
man <- thor_model("rt", endogenous = p$endogenous, exogenous = p$exogenous,
                  coefficients = p$coefficients, equations = p$equations,
                  compile = FALSE, workdir = file.path(work, "src2"), verbose = FALSE)
ok("manual input gives the same model",
   identical(base@equations$new_formula, man@equations$new_formula) &&
   identical(base@vars, man@vars))

cat("\nPASS\n")
