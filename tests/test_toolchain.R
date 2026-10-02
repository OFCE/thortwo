## A broken C++ toolchain is diagnosed, not reported as a thortwo bug.
##
## A forgotten personal Makevars once added
##   CXXFLAGS=-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1
## which, after an SDK update, broke every Rcpp compile on the machine with a
## page of errors from system headers. This checks that the flag is named as
## the likely cause, both by thor_check_toolchain() and in the error from a
## failed model build. The broken Makevars used here fails on every platform,
## not only on the macOS setup that started it.
##
## Usage:  Rscript tests/test_toolchain.R

source("tests/helper.R")

work <- new_tmpdir("thortwo_toolchain_")
on.exit(unlink(work, recursive = TRUE), add = TRUE)

## ---- this machine ------------------------------------------------------
res <- thor_check_toolchain(quiet = TRUE)
ok("this machine compiles Rcpp + RcppEigen", res$ok,
   if (!res$ok) paste(res$problems, collapse = "; ") else "")

## ---- known-bad flags are recognised -------------------------------------
p <- suspicious_flags(c(
  "CXXFLAGS=-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1",
  "CPPFLAGS += -isysroot /no/such/MacOSX14.4.sdk",
  "PKG_CPPFLAGS = -I$(R_HOME)/include",        # make variables are left alone
  sprintf("CPPFLAGS = -I%s", R.home("include")) # an existing ordinary path is fine
))
ok("libc++ include and missing sysroot flagged, nothing else", length(p) == 2L &&
   grepl("standard library", p[1]) && grepl("does not exist", p[2]),
   paste(length(p), "flagged"))

## ---- a broken Makevars, end to end --------------------------------------
## `-include` of a missing header fails on any compiler.
bad <- file.path(work, "Makevars")
writeLines(c("CXXFLAGS = -I/no/such/sdk/usr/include/c++/v1 -include /no/such/header.h"), bad)
old <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
Sys.setenv(R_MAKEVARS_USER = bad)

res <- thor_check_toolchain(quiet = TRUE)
ok("check fails under the broken Makevars", !res$ok)
ok("check names the Makevars and the flag",
   unname(res$makevars["user"]) == bad && any(grepl("c++/v1", res$problems, fixed = TRUE)))

msg <- tryCatch({
  thor_model("toolchain_probe", "inst/models/model_type.txt", backend = "sparse",
             workdir = file.path(work, "build"), verbose = FALSE)
  ""
}, error = conditionMessage)
ok("model build fails with a diagnosis",
   grepl("thor_check_toolchain", msg) && grepl(bad, msg, fixed = TRUE) &&
   grepl("Likely cause", msg))

if (is.na(old)) Sys.unsetenv("R_MAKEVARS_USER") else Sys.setenv(R_MAKEVARS_USER = old)

cat("\nPASS\n")
