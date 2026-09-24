# thortwo

The model solver from [tresthor](https://github.com/OFCE/tresthor), extracted
into a package that does one thing: parse a model, build it, solve it, store
the result.

```r
m   <- thor_model("opale", "opale.txt")          # parse, differentiate, compile
res <- thor_solve(m, "2015Q1", "2024Q4", data)   # Newton, per block, per period
```

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
| `sparse` (default) | triplets | Eigen `SparseLU` | anything, and the only option above a few thousand equations |
| `dense-cpp` | dense matrix | Eigen `PartialPivLU` | small models; the classic path |
| `dense-r` | dense matrix | R `solve()` | no compiler available |

All three go through the same expression translator and the same Newton
iteration, so they agree to machine precision. Measured on Opale (496
equations, 40 quarters), all three take the same 459 Newton iterations and
agree to a median relative difference of ~1e-15.

| | build | solve, 40 periods |
|---|---|---|
| sparse | 7.2 s | 0.010 s |
| dense-cpp | 6.3 s | 0.097 s |
| dense-r | 1.3 s | 4.3 s |

On ThreeME 8x8 (3628 equations, a 2517×2517 heart block) the gap is 0.72 s
against 53 s, which is what the sparse path exists for.

## Caching and persistence

Both are independently switchable, and both work by default.

**Compile cache.** Generated source goes to a stable, model-derived location
under `thor_workdir()`, so Rcpp's cache recognises it across sessions. A second
build of an unchanged model costs 1.7 s instead of 8.3 s. Turn it off with
`cache = FALSE`, move it with `options(thortwo.cache.dir = "...")`, empty it
with `clear_model_cache()`.

**Persistence.** `thor_save()` writes one `.rds` holding the model *and* its
generated source as a string; `thor_load()` writes that source out on the
machine reading it and compiles it, hitting the cache when the code has not
changed. A saved model is therefore portable: it can be mailed, committed, or
reused after a reboot. Whether compilation actually happens is decided by the
cache, not by whether you saved.

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
