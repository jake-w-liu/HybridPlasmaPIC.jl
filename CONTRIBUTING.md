# Contributing to HybridPlasmaPIC.jl

Thanks for your interest in improving HybridPlasmaPIC.jl. This document explains how
to run the test suite, add a benchmark, and the correctness standard every
contribution is held to.

## Running the tests

The full suite lives in `test/runtests.jl` (one `include` per test file). From
the repository root:

```bash
# Instantiate the package environment once
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run the whole suite (this is exactly what CI runs)
julia --project=. -e 'using Pkg; Pkg.test()'
```

`Pkg.test()` activates the `[targets].test` environment, which adds the
test-only dependencies (`Test`, `Statistics`, `Aqua`, `JET`, `HDF5`). To iterate on a
single test file while developing, run it against the package project:

```bash
julia --project=. test/test_hybrid.jl
```

Note that test files relying on a test-only dependency must be run through
`Pkg.test()` so those packages are on the load path. In particular,
`test_quality.jl` requires Aqua and `test_extensions.jl` requires HDF5. The
`test_jet.jl` file is the exception: it skips when JET is absent.

### Static analysis (JET) and quality (Aqua)

- `test/test_quality.jl` runs `Aqua.test_all` — no method ambiguities, no
  piracy, no stale or under-constrained dependencies, no undefined exports.
- `test/test_jet.jl` runs `JET.report_package` as a **report-only** smoke test:
  it prints the findings and asserts the analysis completed, but does not gate
  CI on JET reports that originate in upstream packages (FFTW, Base). When you
  add code, skim the JET output for new instabilities in *our* methods.

## Formatting

The repository ships a `.JuliaFormatter.toml` (default style, 100-column
margin). Before opening a pull request, format your changes:

```bash
julia -e 'using Pkg; Pkg.add("JuliaFormatter")'   # once
julia -e 'using JuliaFormatter; format(".")'
```

CI runs a formatting check and will flag a diff if the tree is not formatted.

## Adding a benchmark

Physics benchmarks are verification tests that compare a solver output against
an independent analytic oracle. To add one:

1. Add a `test/test_<name>.jl` file. State the analytic oracle in a comment at
   the top (the closed-form result, its source, and the physics-justified
   tolerance).
2. Inside an `@testset`, run the solver and `@test` the measured quantity
   against the oracle within tolerance. Print the measured value and the error
   so the log is self-documenting (e.g.
   `@info "whistler dispersion" measured expected rel_err`).
3. Register the file by adding `include("test_<name>.jl")` to
   `test/runtests.jl`.
4. Keep grids small enough to run in CI (seconds, not minutes) while still
   resolving the physics you are validating.

## The CRC / verification standard

Every contribution must satisfy, in priority order:

1. **Correct** — verified against the stated analytic oracle (or, for
   infrastructure, against observed behavior). Trace edge cases: empty inputs,
   boundaries, off-by-one, error paths. Do not present code as correct unless
   you have run it.
2. **Robust** — handles the full realistic input range, not just the happy
   path. No hard-coded magic numbers unless inherent to the problem (justify
   them in a comment). No workarounds that only make code *appear* to work
   (faked results, swallowed errors, hard-coded expected outputs).
3. **Complete** — production-grade: real error handling, sound resource and
   performance behavior, no stubs / `TODO`s / silently skipped cases unless
   explicitly called out.

Verification is by *doing*, not recalling: read the source, run the test, check
with a tool. A correct answer late beats a wrong answer fast.

## Design rules specific to this package

- **Design for 3D, verify in 1D.** New kernels should be written dimension-
  parametrically (the spatial dimension `D` is a type parameter); add a small
  1D verification test.
- **No global mutable state.** Configuration is passed explicitly through
  `HybridStepper` / `ParticleSet` / `HybridFields`.
- Particle positions carry exactly `D` coordinates; velocities always carry 3.
- Units are Ω_ci-normalized (see the README normalization table). Use
  `set_density_weight!` so a uniform load deposits `n = n0 = 1`.

## Pull requests

- Branch off `main`; do not commit directly to it.
- Ensure `Pkg.test()` is green and the tree is formatted before requesting
  review.
- Update `CHANGELOG.md` under `[Unreleased]` describing your change.
