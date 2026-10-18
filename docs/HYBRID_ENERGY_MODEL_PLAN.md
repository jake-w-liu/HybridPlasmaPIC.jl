# Hybrid energy accounting and pressure extension

Status: the diagnostic documentation is corrected to match the implemented
force and reduced inertia model. The coupled evolving-pressure closure is
implemented for the ideal massless case and passes both native Julia acceptance
gates. Exact-revision GitHub CI remains pending. Dissipative recoil/heating and
complete electron inertia remain separate model extensions.

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
Local implementation acceptance is complete. The pre-implementation baseline
confirmed that Ohm preparation overwrote `pe`, restart captured no pressure,
budgets reconstructed pressure from density, and the standalone Euler routine
was first order. The new closure uses persistent pressure and joint RK4 stages.

The final isolated v4 contains 186 frozen source/test files, recorded in
`/tmp/hybrid-evolving-pressure-v4-hashes.json`. Compared with v2, only the stored
pressure energy owner and its regression registration change; 183 files remain
byte-identical. Removing the new energy helpers and reverting their single call
recovers the entire v2 integrator exactly (`v4-verification/scope.json`).

The physical accuracy overlay passes 118 assertions on each Julia 1.12.7 and
1.13.2: direct production RHS versus an independent DFT, periodic pressure work,
coupled convergence, and genuine 3D translation on a 7×9×11 grid. Its SHA256 is
`0f0dbee909fa4f7a7e778bd65f14e1b52acca4912b918a7c4d6fda97681f9b63`;
commands and terminal statuses are in
`/tmp/hybrid-evolving-pressure-final-test-overlay/results.json`. Four deliberately
incorrect implementations fail the unchanged tests on both runtimes: replacing
electron velocity by ion velocity, freezing the pressure gradient, changing the
compression coefficient, and omitting z transport.

Real MPI checks pass at 1/2/4 ranks on 1.13.2 and 1/2 ranks on 1.12.7,
including restart, migration, rank-local rejection and rollback. Commands and
frozen-source checks are in
`/tmp/hybrid-evolving-pressure-io-mpi-v1/freeze1-results.json`. The v2 two-rank
reruns also pass on both runtimes. Each runtime's 24 concrete I/O/MPI analyses
match v1's complete report trees; all eight Serialization diagnostics match
the baseline Serialization subtrees
(`v2-verification/iojet{112,1132}-comparison.json`). This is a
no-added-diagnostics result, not a zero-report claim.

Review repaired direct Ohm-law mutation through shared pressure/gradient storage
and an additional package-analysis report caused by a captured particle count.
The latter has zero `Core.Box` occurrences on both runtimes. The complete v3
package report on 1.12.7 matches the baseline: 22 reports in 14 full trees.
Only source paths, offsets of unchanged source lines and generated identifiers
are normalized (`v3-verification/canonical-delta112.log`). A separate same-object
rendering experiment explains the Serialization display difference in the older
v1 comparison; the original added-count report is retained. No report ceiling
is raised. The full v2 gates exit 0 on both native runtimes: 2,290 assertions
each across 66 testsets, including all 48 concrete signatures, existing
integration regressions, pressure physics and I/O, JET and Aqua. Their complete
package comparisons match the baseline with independently verified same-object
Serialization display translations (`v2-verification/v2-report-delta112.log`
and `v2-report-delta-proved1132.log`). Final v4 complete report comparisons
also match the baseline on both
runtimes: 22 reports/14 trees on 1.12.7 and 21 reports/13 trees on 1.13.2
(`v4-verification/canonical-delta{112,1132}.log`).

The frozen v2 resource rerun passes 336 assertions on each Julia 1.12.7 and
1.13.2 (both exit 0),
including exact algebraic outputs and allocations, unchanged standalone Euler,
adjacent CAM-CL/Hall behavior and grid-size-independent field-stage allocations.
Two subcycles allocate 2,624 bytes in 1D and 2,880 in 2D/3D at 4, 8 and 16 cells
per axis, for Float32 and Float64. These are measured allocations, not a portable
runtime claim. Commands and frozen-source checks are in
`/tmp/hybrid-evolving-pressure-v2-verification/support-results.json`.

**VERIFIED candidate defects and repairs:** the new stored-pressure integral
initially loses finite results through sum, quotient and volume intermediates.
The independent rational witness passes 32 assertions on each runtime. The v3
repair passes 875 integrated assertions and 260 four-thread checks on 1.12.7,
but a subsequent eight-assertion witness on each version shows that rounded
finite sums can conceal true output overflow. The v4 owner adds a conservative
power-of-two bound and retains a streamed exact fallback through `_from_exact`.
The final energy, concrete-compiler and Aqua gate passes 996 assertions on each
runtime with exit 0; a separate four-thread gate passes 650 assertions per
runtime, including both sides of the overflow boundary. The complete package
report trees match the baseline on both versions. Commands, elapsed times and
frozen hashes are recorded in
`/tmp/hybrid-evolving-pressure-v4-verification/results.json`.
The broader v2 regressions and repeated MPI/I/O compiler gates are complete
on both versions (`v2-verification/{minimum,latest,support}-results.json`).
All failed/predecessor evidence remains in
`v2-verification` and `v3-verification`; the pre-edit contracts are
`v2-verification/pressure-energy-brief.md` and `v4-verification/brief.md`.

The public guide example runs on both native versions. The staged documentation
and formatting preserve executable ASTs across all 186 candidate Julia files;
all 20 changed Julia files pass JuliaFormatter 1.0.62 idempotence. This does not
replace exact-revision GitHub CI. An environment
review found stale v2 `[sources]` entries alongside v4 manifest paths; direct
and Aqua-style child probes confirmed the v4 source actually loaded. Corrected
private environment copies pass Aqua on both versions. Latest-Julia resolution
updates six bundled runtime version metadata entries; both environment probes
load identical bundled package versions and roots. The original failed graph
identity assertion and its exact differences are retained in
`v4-verification/environment-source-audit.json` and its linked evidence.

- [x] Define the first slice: opt-in `EvolvingElectrons(gamma>1)`, periodic
  `HybridStepper`, `η=ηH=de2=0`; reject unsupported CAMCL/Hall-MHD use explicitly.
- [x] Define persistent state: users initialize `fields.pe` before `init!`;
  finite nonnegative pressure, including exact zero, is accepted. No clipping
  or implicit density-derived reset.
- [x] Share the pressure RHS with the standalone equation; advance B and pe
  together through each RK4 stage, using stage current and frozen midpoint
  ion moments to form electron velocity.
- [x] Integrate serial single/collection populations and replicated MPI with
  one field-advance owner. Keep existing staggering and variable-step behavior.
- [x] Reject invalid pressure/stages with contextual errors. For new numerical
  rejection, restore public fields, particle positions/velocities and clocks
  exactly; coordinate MPI failure before subsequent collective operations.
- [x] Preserve pe in checkpoint/archive/distributed restart, validate before
  restore mutation, and require explicit `pe=` in evolving-pressure budgets.
- [x] Preserve finite stored-pressure integrals across sum, quotient and volume
  range failures; check exact rounding boundaries and unchanged ordinary costs.
- [x] Retain old algebraic outputs and steady-state allocations. Allocate
  pressure/rollback buffers only for the new closure and measure their cost.
- [x] Verify analytic advection, compression, Hall transport and battery
  coupling; joint-stage and coupled convergence; zero/variable timesteps;
  rollback/retry, restart and MPI failure; compiler/quality gates on both Julia
  versions.
- [ ] Complete exact-revision GitHub CI for the pushed implementation.

The explicit spectral pressure update is not positivity preserving;
smaller timesteps or better-resolved positive data may be necessary. The
existing magnetic timestep recommendation does not bound pressure advection or
compression. These constraints must be documented with the implementation.
