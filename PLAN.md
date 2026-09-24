# Plan: an "advanced mode" with selectable solving strategies

Status: **option 1 is done** (2026-10-04, as `thor_model(sequential = TRUE)`),
and so is **factorisation reuse**, the cheap precursor of option 3
(2026-10-05, as `thor_solve(reuse_jacobian = )`, on by default on the
compiled backends); results below. Broyden
itself and option 2 (Gauss–Seidel) are not started. Written 2026-10-03 to record
what was decided and measured, so the work can be picked up later.

## The idea

Three alternative solving strategies, to be added and tested one at a time:

1. **Sequential prologue and epilogue**: solve these blocks one equation at a
   time, in dependency order, instead of as a Newton system.
2. **Gauss–Seidel** on the heart.
3. **Broyden** on the heart.

## Ground rules

- **Not a rewrite.** The current behaviour (Newton on every block) stays
  the default and stays as it is. With advanced mode off, the generated code
  must be byte-for-byte what it is today, so the compile cache, saved models
  and every existing test are unaffected.
- **Opt-in, in one place.** One options object selects the strategy; nothing
  else in the interface changes.
- **Each option is judged against Newton** on the existing fixtures before it
  is documented as usable: same solution, and a measured benefit.

## Interface

Option 1 was built as a plain argument, as asked:

```r
model <- thor_model("m", "m.txt", sequential = TRUE)   # default FALSE
result <- thor_solve(model, from, to, data)            # unchanged
```

The choice is made at **build** time, because it needs different generated
code, and is recorded in `model@meta$sequential`.

For options 2 and 3, a second argument in the same style would do, e.g.
`heart = c("newton", "broyden", "gauss-seidel")`. Gauss–Seidel is a build-time
choice for the same reason as option 1; Broyden could be chosen at solve time,
since it uses the code Newton already has.

## Shared groundwork (needed by 1 and 2)

**A. Which variable each equation determines.**

- Prologue and epilogue: `decomposing_model()` already finds this, and the
  order to solve in, as it peels equations off. `thor_model()` then sorts the
  equations and variables alphabetically and the order is lost
  (`R/build.R`, `block_spec`). Keep it: return the (equation, variable) pairs
  in peel order from `decomposing_model()`.
- Heart: needs a matching of equations to variables (maximum bipartite
  matching). A scratch version using igraph exists from the component
  analysis; thortwo should not depend on igraph for this, so write a small
  augmenting-path matcher, or start from "the variable on the left-hand side"
  and repair what is left.

**B. Turning an equation into an assignment.** Two cases:

- *Explicit*: the equation is `x = f(...)` with `x` its matched variable and
  `x` not on the right. Generate `x <- f(...)` directly. Measured share on
  ThreeME 13x13: prologue 95%, heart 68%, epilogue 81%.
- *Implicit*: anything else (`log(y) = ...`, `pk * f_k = ...`). Solve for the
  one variable with a scalar Newton step, which needs one derivative: the
  equation with respect to its own variable.

**C. Code generation** for the assignment form, in R (`codegen-rfuns.R`) and
C++ (`codegen-sparse.R`), next to the existing residual/jacobian generators,
not in place of them.

Start with the R backend (`sparse-r`): no compile step, so it is quick to
iterate on, and since 2026-10 it is fast enough to test at ThreeME size.

## Option 1: sequential prologue and epilogue (done)

**Result (2026-10-04).** Implemented in `R/sequential.R`, for all three
backends; tested in `tests/test_sequential.R`.

- Correct: same solution as Newton to within tolerance on Opale, ThreeME 4x4,
  8x8 and 13x13 (median relative difference at most 1e-14).
- With `sequential = FALSE` the generated code is byte-for-byte unchanged
  (checked by hash against the previous commit, three backends, two models).
- **The time benefit is small, less than expected below.** ThreeME 13x13,
  sparse: R part 29 → 26 s, compile 18 → 17 s, code 2,064 → 1,740 KB, solve
  unchanged. ThreeME 4x4: slightly *slower* to build (code 396 → 475 KB).
  Reasons: a sequential block still needs its residual code (for
  `model_residuals()` and the convergence report), so its equations are
  emitted twice; and the equation that made these blocks expensive,
  `verif_all`, had already been dealt with by `stats::D`.
- What it does give: an error naming the equation and the variable that
  cannot be determined, and a faster solve where the blocks' linear algebra
  was the cost (dense-cpp on Opale: 0.095 → 0.009 s; sparse-r on 13x13:
  19 → 16 s).
- Share of the two blocks' equations solved directly: Opale 74%, ThreeME 4x4 77%, 13x13 86%.
- Possible follow-up: emit each equation once and build both the assignment
  and the residual from it, which would remove the duplication.
- **Correction, 2026-10-05: on 29x33 the benefit was large**, as long as
  `verif_all` was a root of a sum of squares: compile 112 → 56 s, code
  9.8 → 7.8 MB. The 257 derivatives of that one epilogue equation, each
  repeating its 286-term sum, were 2.3 MB of code and about 90 s of the
  compile, and sequential mode needs none of them. With the root removed from
  the model (compile 44 s in the default mode), the two modes are close again.

The text below is the plan as written before it was built.

**What it is.** These blocks are recursive by construction: in the right order
each equation only needs variables already computed. So one pass in that
order is exact. No iteration, and no "repeat until convergence" across blocks
either: the prologue does not depend on the heart, the heart only on the
prologue, the epilogue on both.

**Expected benefit.** Mostly in the build:

- Derivatives: only each implicit equation's derivative with respect to its
  own variable. On 13x13 that removes about 3,800 of 20,300 jacobian entries
  (19%) and their C++.
- It removes the `verif_all` problem at the root. That equation is in the
  epilogue and is explicit, so it needs no derivative at all. (Today it is
  handled by differentiating it with `stats::D`; see `equation_derivatives()`.)
- Solve: little. These blocks are already cheap (their jacobians are
  triangular); the heart dominates. Measure per-block time first.

**Risks.** An implicit equation whose derivative with respect to its own
variable is zero fails exactly as today (`pk_sgzx` in 13x13), only with a
clearer message. Otherwise low: the result is deterministic and must equal
Newton's to rounding.

**Size.** Groundwork A (prologue/epilogue part) and B, C. The smallest of the
three once the groundwork exists, and the groundwork is reused by option 2.

## Option 2: Gauss–Seidel on the heart

**What it is.** Sweep the heart's equations in a fixed order, each updating its
own variable from the latest values of the others; repeat until the sweep stops
changing anything.

**Expected benefit.** The build: no jacobian for the heart at all. The heart
holds 81% of the jacobian entries on 13x13, and at 29x33 the jacobian is half
the R-side build time and most of the generated C++.

**Risks: high.** Gauss–Seidel converges only if within-period feedback is weak
enough, and how fast depends on the equation order. ThreeME's heart is a
**single strongly connected component** (4,663 of 4,663 equations at 13x13,
18,380 of 18,380 at 29x33, measured), so there is no ordering that makes it
recursive. Expect hundreds of sweeps where Newton takes about 4 steps, if it
converges at all.

**Do this first:** a throwaway prototype in R, using the existing residual
functions and a matching, on ThreeME 4x4 and 13x13. Does it converge? In how
many sweeps? With what damping? Only build the generator if the answer is
encouraging.

**Size.** The largest: groundwork A (heart matching), B, C, plus the sweep
loop, convergence test and damping in both R and C++.

## Option 3: Broyden on the heart

**Factorisation reuse is done (2026-10-05): `thor_solve(reuse_jacobian = )`.**
Implemented in `newton_block` (generated C++) and `newton_block_r`; tested in
`tests/test_reuse.R`. **It is the default on the compiled backends**, by
decision of 2026-10-05, and off by default on `sparse-r`. This is the one
exception to the rule above that the default behaviour stays as it was: the
solution is unchanged to 1e-10, but iteration counts are not.

- Where the solve time goes on ThreeME 29x33 (5.0 s per period, 12 Newton
  iterations per period): one factorisation of the heart's jacobian takes
  1.36 s and there are 4 per period. Using the factors takes 5 ms, evaluating
  equations and derivatives a few ms. The factors hold 6.3 M entries for a
  jacobian of 67,610 (fill 93x), with COLAMD, by far the best ordering Eigen
  offers (AMD and natural ordering: 870x, 24 s).
- Result: 29x33 208 → 74 s (2.8x), 13x13 3.8 → 1.2 s (3.2x). Iterations
  456 → about 1,000. Same solution to 1e-10, same residuals.
- Two safeguards: a step with old factors that does not reduce the residual
  is discarded and redone fresh; factors are kept only while they halve the
  residual.
- Steps with old factors must be 1,000 times smaller than the tolerance to
  count as converged (`REUSE_TIGHT`). Without that, residuals were 3e-6 on
  Opale instead of 2e-10: such steps converge linearly, not quadratically.
- **Tried and dropped: carrying the factors across periods.** No faster on
  29x33 (76 s against 74 s) and ten times less accurate. A new period needs
  a fresh factorisation anyway.
- No gain on `sparse-r` (evaluation in R dominates); slower on small models.
- **So the remaining floor is one factorisation per period** (1.36 s on
  29x33, i.e. about 52 of the 74 s). Broyden would not lower it: it also
  starts from one factorisation. What would is a cheaper factorisation: the
  fill comes largely from a few variables that appear in a great many
  equations (`pch` in 177, `p` in 136, the `pyq_*` in 114 each). Setting those
  aside (a bordered or Schur-complement solve) is untested: a first probe was
  badly designed and proved nothing.

The text below is the plan as written before this was built.

**What it is.** Newton recomputes and refactorises the jacobian at every
iteration. Broyden factorises once and then corrects the solution using
residuals only. A dense Broyden update is impossible at this size (the 29x33
heart would need 2.7 GB), so it has to be the limited-memory form: keep the
sparse LU of the starting jacobian and a short list of correction vectors;
restart from a fresh jacobian if convergence slows.

**Expected benefit.** Solve time only, perhaps 2x: more iterations, each much
cheaper. Nothing for the build, since the starting jacobian still needs the
derivatives.

**Try this first:** reuse the factorisation for several iterations, or across
periods, without any Broyden correction (about 20 lines in each solver). This
was observed by accident on 2026-10-03, through a stale cached factorisation
in the R backend: Opale converged in 1,514 iterations instead of 459, each
cheaper. If that already gives most of the gain, Broyden is not worth adding.
Measure what share of the solve is factorisation before either.

**Risks.** Low: it can always fall back to a Newton step.

**Size.** Medium: a second loop next to `newton_block` in the C++ runtime
(`codegen-runtime.R`) and next to `newton_block_r` (`solve-r.R`). No
groundwork from A–C needed.

## Suggested order

1. ~~Groundwork A–C for the prologue and epilogue, then **option 1**.~~ Done.
   The ordering (`sequential_order()`), the direct/scalar classification and
   the assignment-form generators are there for option 2 to reuse; what
   option 2 still needs is the matching of heart equations to variables.
2. ~~**Factorisation reuse** (the cheap precursor of option 3), measured on
   13x13 and 29x33. Decide on Broyden from the result.~~ Done; Broyden is not
   worth adding (see option 3).
3. **Gauss–Seidel prototype** in R. Decide on option 2 from the result.

## How each option is tested

For every option, on Opale, ThreeME 4x4, 8x8 and 13x13:

- the four existing invariants (`tests/helper.R`): converged, residuals below
  1e-6, backends agree, history reproduced;
- the solution agrees with Newton's (median relative difference below 1e-10);
- with advanced mode off, the generated code hash is unchanged from before
  the feature existed;
- build time, solve time and iteration counts recorded next to Newton's, and
  added to `doc/backends.qmd` once an option is worth recommending.

## Reference measurements (2026-10-03, Apple M1)

| | 13x13 | 29x33 |
|---|---|---|
| equations (prologue / heart / epilogue) | 792 / 4,663 / 1,336 | heart 18,380 of 26,439 |
| jacobian entries (prologue / heart / epilogue) | 1,095 / 16,513 / 2,673 | |
| Newton iterations per period, all blocks | 12 | |
| build, R part / compile (sparse) | 29 s / 18 s (2026-10-05: 10 s / 18 s) | 2026-10-05: 55 s / 112 s |
| solve, sparse / sparse-r | 3.8 s / 19 s | about 3 min 15 s / not measured |
