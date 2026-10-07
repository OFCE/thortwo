# thortwo (development version)

Changes made in October 2026, most recent first. Timings are on an Apple M1
with 16 GB.

## Where ThreeME 29x33 stands

26,439 equations; the heart block alone has 18,380.

| | early October | now |
|---|---|---|
| translating the `.prg` (ermeeth2) | 192 s | 19 s |
| building: work in R | about 54 min | 38 s |
| building: C++ compile | 24 min | 44 s |
| building the same model again | the above, minus the compile | 0.5 s |
| solving 38 years | 196 s | 74 s |

The compile and build figures are with `Verif_ALL` written without its square
root (see below); with it, the compile is 112 s.

## Timings

* **`thor_model()` and `thor_solve()` report how long they took**, by default
  and whether or not `verbose` is on:

      Timings: build 38.4 s, compile 44.1 s, total 1 min 23 s
      Timings: solve 1 min 14 s (38 periods)

  `build` is the work done in R, `compile` the C++ compilation (or, from the
  cache, the time to load it). The build's figures are kept in
  `model@meta$timings`. `timings = FALSE` on either function, or
  `options(thortwo.timings = FALSE)`, switches the display off. A script that
  reads what a quiet build prints will now see this extra line.

## Solving

* **`thor_solve()` reuses the factorised jacobian within a period**
  (`reuse_jacobian`). On a large model, factorising the jacobian is nearly
  the whole cost of a Newton iteration: 1.36 s on the heart of ThreeME 29x33,
  against milliseconds for everything else. The solver now factorises once
  per period and reuses the result for that period's other iterations,
  redoing it only when a step stops reducing the residual. ThreeME 29x33
  solves in 74 s instead of 208 s, 13x13 in 1.2 s instead of 3.8 s. Same
  solution to 1e-10, same residuals; the number of iterations reported is
  about twice as high.

  **This is the default on the compiled backends** (`sparse`, `dense-cpp`).
  It is off by default on `sparse-r`, where evaluating the equations is the
  main cost. `reuse_jacobian = FALSE` gives plain Newton on any backend.

* **A value that is not a number is never returned as a solution.** The
  compiled backends measured a block's residual with a function that skips a
  `NaN` unless it is the first element, so a block whose residual was `NaN`
  in any other equation was not seen to have failed, and the solver could
  report convergence while returning `NaN` for some variables. The solve now
  stops with an error. `model_residuals()` had the same blind spot and now
  reports `NaN`.

* **New advanced option `thor_model(sequential = TRUE)`.** The prologue and
  the epilogue are solved one equation at a time, in dependency order,
  instead of by Newton: an equation written `x = f(...)` is evaluated, any
  other is solved for its one unknown by a one-variable Newton. Same solution.
  Its main use is the error it gives when a variable cannot be determined,
  which names the equation and the variable. Off by default; with it off the
  generated code is unchanged.

## The pure-R backend

* **`"sparse-r"` replaces `"dense-r"`.** The backend that needs no compiler
  now solves with a sparse jacobian and the Matrix package's sparse LU
  instead of a dense matrix and base `solve()`. ThreeME 4x4 solves in 2.3 s
  instead of 48 s, 8x8 in 8 s instead of nearly 7 minutes, and 13x13, which
  could not be built at all, in 19 s. `"dense-r"` is still accepted as the old
  name, and models saved with it are converted when loaded.

## Building

* **Built models are cached.** `thor_model()` stores every model it builds
  under a hash of its name, variables, equations, options and thortwo's own
  code. Asked for the same model again, in this session or a later one, it
  loads the stored one and says so. `recompile = TRUE` forces a fresh build;
  `cache = FALSE` or `options(thortwo.model.cache = FALSE)` switch the lookup
  off.

* **Differentiation uses `stats::D` first, and Deriv only as a fallback.**
  Deriv simplifies every derivative, and its simplifier was 89 of the 115 s
  the jacobian took on ThreeME 29x33; one equation, `verif_all`, would have
  taken about 55 minutes on its own. `stats::D` gives numerically identical
  derivatives. `delta()` is handled directly. An equation containing a
  function `stats::D` does not know (`abs`, `sign`, `asinh`, `acosh`, `atanh`,
  `logb`) still goes to Deriv.

* **The C++ compiles about five times faster.** The generated code was meant
  to be split into small functions, but the compiler was merging them back
  into a few very large ones. They are now kept separate and much smaller (25
  statements instead of 1,500).

* **The decomposition keeps running counts** instead of recounting every
  variable's occurrences on each pass: 22 s to 0.3 s on ThreeME 29x33, with
  identical blocks.

* **Reading the model is quicker**: 20 s to 3 s on ThreeME 29x33. The
  endogenous variables of every equation are now identified in one pass
  (12.6 s to 1.2 s), and the names used in the equations are extracted once
  instead of once for each of the three lists of variables (6 s to 2 s).
  The generated code is unchanged.

* **Variable positions are looked up through a hash table** during code
  generation: 56 s to 9 s on ThreeME 29x33. The generated code is unchanged.

* **Saved models are portable across language settings.** "Alphabetical" is
  not the same on every computer (the underscore sorts differently), and a
  model saved on one could be refused on another, because the validity check
  compared the stored variable order with a fresh sort. The order is now fixed
  when the model is built and never sorted again.

* **`thor_check_toolchain()`**, and a diagnosis appended to compile errors.
  When a model fails to compile with errors pointing into system headers, the
  cause is usually the computer's C++ setup, typically a forgotten personal
  `Makevars`. The check compiles a few lines of Rcpp, lists the `Makevars`
  files in effect and flags the lines known to break compilation.

## Fixes

* **The pure-R backend no longer rejects a badly scaled block as singular.**
  The Matrix package refuses to solve when a jacobian's smallest and largest
  pivots are more than about 1e16 apart. That is a test of scaling, not of
  singularity, and the compiled backends make no such test. It stopped
  ThreeME 4x4 on `"sparse-r"` at the first period ("Newton did not converge
  on block 'epilogue' ... scaled step Inf"), where `"sparse"` solved it. The
  test is now switched off; a matrix that really is singular is still caught.
  Needs Matrix 1.6-0 or later.

* A named equation written with a space before the colon
  (`demand : y = ...`) made the first variable after the colon look unused,
  and it was dropped.

## Documentation

* Three guides in `doc/`: building and solving a model, writing equations,
  and choosing a backend, with benchmarks from 3 to 6,791 equations.
* `PLAN.md` records what was tried and measured for the alternative solving
  strategies, including what was dropped.

## Test fixtures

* `Verif_ALL` is written without its square root in `tests/model_13_13.prg`,
  `tests/model_29_33.prg` and `tests/threeme_13x13_thor.txt`, following the
  change made in ThreeME. With the root, each of its derivatives repeated the
  whole sum (2.3 MB of generated code and about 90 s of the 29x33 compile),
  and its derivative did not exist at zero, where the calibration puts it.
