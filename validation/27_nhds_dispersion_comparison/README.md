# Case 27 — Plasma code-to-code dispersion comparison vs NHDS

A **physics** comparison against another **open-source plasma code** — not an
infrastructure/library check. We overlay our hybrid whistler dispersion ω(k) against
**NHDS** (the New Hampshire Dispersion relation Solver, Verscharen & Chandran 2018,
[github.com/danielver02/NHDS](https://github.com/danielver02/NHDS), BSD-2-Clause),
which independently solves the full **hot-plasma Vlasov–Maxwell** dispersion relation.

This is the analog of validating an EM solver against Meep: two independent codes,
same physics, lines should overlay — and where they don't, the difference is
diagnostic (a real bug, or a known model limit).

## What is compared

- **Ours:** the independent Hall–MHD eigenvalue oracle (`test/oracles`, HYB-006) — the
  highest parallel branch (the whistler / R-mode) at β=1.
- **NHDS:** the bundled `whistler.in` case — β=1, parallel (θ=0.01°), proton + real-mass-ratio
  electron, full kinetic. Output `ω(k)` over `k·d_p ∈ [0.2, 4]`.

Same normalization (k·d_p, ω/Ω_p; d_p = ρ_p at β=1).

## Result (honest)

| regime | agreement |
|---|---|
| k ≥ 2 (fluid ≈ kinetic) | **≤ 1.2%** — gated |
| full range k = 0.2–4 | **≤ 4.9%** — gated |
| k ≈ 0.2–1 (k·ρ_i ~ 1) | a few % — the **expected FLR/thermal fluid-vs-kinetic difference**, recorded as a skip |

Our fluid-closure whistler reproduces the kinetic whistler to ≤1.2% where the closure
is valid, with the few-% spread near `k·ρ_i ~ 1` being the physical finite-Larmor-radius
correction that NHDS (full Vlasov-Maxwell) includes and Hall-MHD does not. No bug; a
clean, quantified validation against an external kinetic code.

(Note: NHDS's Newton solver tracks the whistler at β=1; a naive β→0 "cold" run instead
converges to the *ion-cyclotron* branch, so it is **not** used here — that would be a
branch mismatch, not a comparison.)

## Running it

NHDS source/binary/output are **gitignored** (regenerated, never committed — keeps the
repo small). Build + run NHDS, then the case:

```
bash validation/27_nhds_dispersion_comparison/build_nhds.sh   # needs git + gfortran + make
julia --project=validation validation/27_nhds_dispersion_comparison/run.jl
julia --project=validation validation/27_nhds_dispersion_comparison/plot.jl   # overlay PDF
```

`run.jl` auto-runs `build_nhds.sh` if the NHDS output is missing, and **skips cleanly**
if NHDS cannot be built (no network/gfortran). The case is `default=false` (research
validation; not in the package CI test suite, which has no gfortran/network).
