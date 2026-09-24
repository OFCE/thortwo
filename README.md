# thortwo

The model solver from [tresthor](https://github.com/OFCE/tresthor), extracted
into a package that does one thing: parse a model, build it, solve it, store
the result.

```r
m   <- thor_model("opale", "opale.txt")          # parse, differentiate, compile
res <- thor_solve(m, "2015Q1", "2024Q4", data)   # Newton, per block, per period
```

## Documentation

In `doc/`, as Quarto sources with rendered `.html` next to them:

- [Building and solving a model](doc/getting-started.qmd): a first script,
  what the data frame must look like, and every option of `thor_model()` and
  `thor_solve()`.
- [Writing model equations](doc/equations.qmd): the model file, the equation
  syntax, supported functions, conditionals, and how differentiation works.
- [Choosing a backend](doc/backends.qmd): dense or sparse, with benchmarks
  from 3 to 6,791 equations.
- [Solving ThreeME with thortwo](doc/threeme.qmd): a worked example.

## What it does

```
model .txt  ──parse──>  equations table + endo/exo/coeff vectors
            ──decompose──>  prologue / heart / epilogue blocks
            ──differentiate──>  symbolic jacobian
            ──generate──>  C++ source, or R closures
            ──compile──>  .so, cached between sessions
            ──solve──>  Newton per block per period
```

Everything else tresthor does — model modification, estimation, single-equation
solving, plots, contributions, the variable dictionary — is deliberately out of
scope. The boundary is the point.

## Backends

Three, chosen at build time with `backend =` and stored on the model, so
`thor_solve()` has one entry point:

| backend | jacobian | linear algebra | when |
|---|---|---|---|
| `sparse` (default) | triplets | Eigen `SparseLU`, compiled | normal work: the fastest solve, and cached |
| `sparse-r` | triplets | `Matrix` sparse LU, pure R | no compiler, or one quick run of an edited model |
| `dense-cpp` | dense matrix | Eigen `PartialPivLU`, compiled | an independent cross-check of `sparse` |

`dense-r`, the former pure-R backend, is still accepted and now means
`sparse-r`; models saved with it are converted when loaded.

All three read the formulas the same way and run the same Newton iteration,
so they agree to machine precision. Measured on Opale (496
equations, 40 quarters), all three take the same 459 Newton iterations and
agree to a median relative difference of ~1e-15.

| | build | solve, 40 periods |
|---|---|---|
| sparse | 7.2 s | 0.013 s |
| sparse-r | 1.4 s | 0.65 s |
| dense-cpp | 5.8 s | 0.10 s |

On ThreeME 8x8 (3628 equations, a 2517×2517 heart block) the solve takes
0.80 s on sparse, 8.1 s on sparse-r and 56 s on dense-cpp, which is what the
sparse path exists for. Full benchmarks: [doc/backends.qmd](doc/backends.qmd).

## How the solve is kept fast on large models

On a large model, factorising the jacobian is nearly the whole cost of a
Newton iteration: 1.36 s on the heart of ThreeME 29x33, against milliseconds
for everything else. So on the compiled backends `thor_solve()` factorises
once per period and reuses the result for that period's other iterations,
redoing it only when a step stops reducing the residual. ThreeME 29x33 solves
in 74 s instead of 208 s, with the same solution to 1e-10 and about twice the
iteration count. `reuse_jacobian = FALSE` gives plain Newton; it is the
default on `sparse-r`, where evaluating the equations is the main cost.

What changed, and when, is in [NEWS.md](NEWS.md).

## Advanced: sequential prologue and epilogue

`thor_model(..., sequential = TRUE)` solves the two recursive blocks one
equation at a time, in dependency order, instead of handing them to Newton as
systems: an equation written `x = f(...)` is evaluated, any other is solved for
its one unknown by a one-variable Newton. The heart is unchanged. The solution
is Newton's to within tolerance, on every backend. Off by default, and with it
off the generated code is byte-for-byte what it was.

On the models measured the effect on build and solve time is small in either
direction (ThreeME 13x13: build 29 → 26 s in R, 18 → 17 s compiling). What it
adds is an error that names the equation and the variable when one cannot be
determined. Details in [doc/getting-started.qmd](doc/getting-started.qmd);
what is planned next in [PLAN.md](PLAN.md).

## Caching and persistence

All three work by default.

**Built-model cache.** `thor_model()` stores every model it builds under a
hash of what determines it: name, variables, equations, backend,
decomposition, and thortwo's own code. Asked for the same model again, in this
session or a later one, it loads the stored one and says so, instead of
parsing, differentiating and generating code again: 0.55 s instead of 30 s plus
the compile on ThreeME 13x13. An edited equation changes the hash;
reformatting the file does not. `recompile = TRUE` forces a fresh build,
`cache = FALSE` or `options(thortwo.model.cache = FALSE)` switch the lookup
off.

**Compile cache.** Generated source goes to a stable, model-derived location
under `thor_workdir()`, so Rcpp's cache recognises it across sessions. A second
build of an unchanged model costs 1.7 s instead of 8.3 s. Turn it off with
`cache = FALSE`, move it with `options(thortwo.cache.dir = "...")`, empty it
with `clear_model_cache()`.

**Persistence.** `thor_save()` writes one `.rds` holding the model *and* its
generated source as a string; `thor_load()` writes that source out on the
machine reading it and compiles it, hitting the cache when the code has not
changed. A saved model is therefore portable, including to a computer whose
language settings sort names differently: it can be mailed, committed, or
reused after a reboot. Whether compilation actually happens is decided by the
cache, not by whether you saved.

## When a model fails to compile

If `thor_model()` fails on `sparse` or `dense-cpp` with compiler errors that
point into system headers (`math.h`, `cstdint`, ...) rather than into the
model's `.cpp`, the machine cannot compile Rcpp code at all. Run:

```r
thortwo::thor_check_toolchain()
```

It compiles a few lines of Rcpp + RcppEigen the way a model build does, lists
the Makevars files R is applying, and flags the lines in them known to break
compilation. The usual culprit is a personal `~/.R/Makevars` (or a file named
by `R_MAKEVARS_USER`) written to work around one compiler or SDK update, such as

```
CXXFLAGS=-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1
```

which then breaks on the next one. Remove the line, restart R. In the
meantime `backend = "sparse-r"` needs no compiler.

## Tests

```sh
Rscript tests/run_all.R        # ~2 min
Rscript tests/run_all.R full   # ~6 min, adds ThreeME 8x8
```

Four invariants, in increasing order of strength:

1. **Converged** — worst scaled Newton step ≤ 1.
2. **Residuals** — `model_residuals()` over every block and period < 1e-6.
   This is the one that catches code-generation bugs; convergence alone does
   not.
3. **Backends agree** — median relative difference < 1e-10.
4. **History reproduces** — solving over history with historical exogenous
   inputs returns the historical endogenous path. On Opale: median 1.4e-15,
   max 1.8e-11 for variables with |value| > 1.

Plus `export_model() |> thor_model()` round-tripping the equations table,
which pins the parser.
