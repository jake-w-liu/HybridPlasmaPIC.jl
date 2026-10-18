# Hybrid energy accounting and pressure extension

Status: the diagnostic documentation is corrected to match the implemented
force and reduced inertia model. A coupled evolving-pressure closure is in
isolated implementation for the ideal massless case. Dissipative recoil/heating and complete
electron inertia remain separate model extensions.

## Verified balance and scope

For smooth periodic scalar-pressure states with inactive density floors, write
`J_i = n u_i`, `u_e = u_i − J/n`, `R = ηJ − ηH∇²J` and `a = de2`. The code uses
the full E in the ion force (`extrapolated_leapfrog.jl::_push_hybrid_population!`
and `particles/boris.jl`) and forms `E₀ = −u_e×B − ∇p/n + R`
(`electrons/ohms_law.jl::_ohm_Efield!`). The inertia filter obeys
`(I + a curl curl)E = E₀`, hence Faraday gives `E = E₀ + a∂tJ`.

With electron internal-energy source Q, the continuum balance is

```
d(K_i + W_B + U_e + a∫|J|²/2)/dt
    = ∫(J_i − J)·R + Q + a J_i·∂tJ.
```

The scalar algebraic closures have no independently evolving heat source.
Conventional Ohmic/hyperresistive heating cancels the integrated `−J·R`
contribution but leaves `J_i·R` unless the ion force also includes recoil.
Choosing `Q = −n u_e·R` would cancel the balance algebraically but is
sign-indefinite; it is not a justified irreversible heating law.

**VERIFIED:** 17 exact symbolic assertions check the identities and two smooth
periodic witnesses in `/tmp/hybrid-pressure-energy-contract-20261016/witness.py`
(SymPy 1.14.0, exit 0; exact values in `witness.json`):

- With `n=1`, constant pressure, `B=(0,0,−sin x)` and `u_i=2J`, the adiabatic
  total rate is `+π(η+ηH)`. Conventional heating makes it `+2π(η+ηH)`; paired
  recoil removes the remaining source.
- With the same B/n/pressure, `R=0` and `u_i=(1,sin x,0)`, the reduced inertia
  model gives total rate `πa/(1+a)`, despite including its J² reservoir.

These are continuum/model findings, not measured timestep errors. The existing
`test/test_budgets.jl` inertia conservation test explicitly uses fixed ions and
uniform density; it does not establish the moving-ion claim. Scalar ideal
massless continuum conservation also does not establish exact conservation by
the current explicit particle/field integrator.

The paired-force requirement agrees with the ion force and pressure source in
[Stanier, Chacón and Chen, equations 1–4 and Appendix A](https://arxiv.org/html/1803.07158).
The terms omitted by reduced inertia schemes are distinguished from the full
electron momentum model in [Muñoz et al., section 2](https://arxiv.org/html/1612.03818).
The balance above is derived from this repository's equations.

## Implementation boundaries

- Correct the `energy_budget` and `resistive_dissipation` documentation without
  changing their returned values, the ion force, or existing validation targets.
- Add persistent, stage-consistent pressure evolution first for
  `η=ηH=de2=0`. Define positivity, rollback, MPI, restart and explicit pressure
  input to diagnostics before implementing the new closure.
- Design dissipative recoil and heat partition together before enabling
  resistive/hyperresistive evolution with a closed-energy claim.
- Treat full finite-electron-inertia dynamics as a separate extension;
  relabeling the quadratic diagnostic cannot supply omitted momentum terms.

## Evolving-pressure implementation tracker

The reviewed design is `/tmp/hybrid-evolving-pressure-design-v1/implementation-plan.md`
(SHA256 `628edafd3d54b15c1bb42576a6d428468bfff850ea26563f5219d5c81005a244`).
Its baseline probes pass 31 assertions on each Julia 1.12.7 and 1.13.2.
Implementation acceptance is in progress. The baseline confirms that current Ohm
preparation overwrites `pe`, restart captures no pressure, budgets reconstruct
pressure from density, and the standalone Euler routine is first order.

The isolated candidate's 185 source/test files remain frozen in
`/tmp/hybrid-evolving-pressure-v1-freeze1-hashes.json`. Its final test overlay
passes 118 assertions on each runtime, including direct production RHS versus
an independent DFT, periodic pressure work, coupled convergence, and genuine
3D translation on a 7×9×11 grid. The exact tested accuracy file has SHA256
`0f0dbee909fa4f7a7e778bd65f14e1b52acca4912b918a7c4d6fda97681f9b63`;
commands and terminal statuses are in
`/tmp/hybrid-evolving-pressure-final-test-overlay/results.json`.

Real MPI checks pass at 1/2/4 ranks on 1.13.2 and 1/2 ranks on 1.12.7,
including restart, migration, rank-local rejection and rollback. Commands,
logs and frozen-file checks are recorded in
`/tmp/hybrid-evolving-pressure-io-mpi-v1/freeze1-results.json`.
Review reproduced direct Ohm-law mutation when pressure shared a gradient
buffer. The shared storage preflight now rejects it before mutation; the
independent recheck exits 0 on 1.13.2
(`/tmp/hybrid-pressure-direct-ohm-alias-recheck.log`).

The first package-quality run on 1.12.7 adds one JET report (23 versus 22),
identifying a captured `count` variable in `_save_pressure_transaction!`.
Concrete calls pass, but the report remains a promotion gate. The isolated v2
candidate gives the total and per-population counts distinct bindings; final
full-report comparisons and affected native gates are running. No ceiling is raised.
The original source and report are retained; the live implementation is
unchanged. Public documentation is prepared separately, with its example
verified on both native runtimes.

- [x] Define the first slice: opt-in `EvolvingElectrons(gamma>1)`, periodic
  `HybridStepper`, `η=ηH=de2=0`; reject unsupported CAMCL/Hall-MHD use explicitly.
- [x] Define persistent state: users initialize `fields.pe` before `init!`;
  finite nonnegative pressure, including exact zero, is accepted. No clipping
  or implicit density-derived reset.
- [ ] Share the pressure RHS with the standalone equation; advance B and pe
  together through each RK4 stage, using stage current and frozen midpoint
  ion moments to form electron velocity.
- [ ] Integrate serial single/collection populations and replicated MPI with
  one field-advance owner. Keep existing staggering and variable-step behavior.
- [ ] Reject invalid pressure/stages with contextual errors. For new numerical
  rejection, restore public fields, particle positions/velocities and clocks
  exactly; coordinate MPI failure before subsequent collective operations.
- [ ] Preserve pe in checkpoint/archive/distributed restart, validate before
  restore mutation, and require explicit `pe=` in evolving-pressure budgets.
- [ ] Retain old algebraic outputs and steady-state allocations. Allocate
  pressure/rollback buffers only for the new closure and measure their cost.
- [ ] Verify analytic advection, compression, Hall transport and battery
  coupling; joint-stage and coupled convergence; zero/variable timesteps;
  rollback/retry, restart and MPI failure; compiler/quality gates on both Julia
  versions. Then push the complete slice and require exact-revision green CI.

The proposed explicit spectral pressure update is not positivity preserving;
smaller timesteps or better-resolved positive data may be necessary. The
existing magnetic timestep recommendation does not bound pressure advection or
compression. These constraints must be documented with the implementation.
