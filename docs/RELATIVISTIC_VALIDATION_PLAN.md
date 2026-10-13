# Relativistic velocity validation

Status: shared boundary repair verified on Julia 1.12.7 and 1.13.1 and pushed
as `43881de`.

## Confirmed defect and ownership

VERIFIED: exact rational sums of the stored velocity components exceed `c²`
for Float16/32/64 examples that the divided/squared/summed guard accepts.
The shared `_require_subluminal_velocities` in `src/particles/boris.jl` serves
radiation, both relativistic electromagnetic PIC APIs, and electromagnetic
checkpoint validation for electrons and mobile ions.

The public baseline probe passes 28 witness assertions on each Julia version:
radiation accepts the three precisions in 1D–3D, and both PIC APIs accept the
Float32/64 examples. Logs:
`/tmp/hybrid-relativistic-public-baseline112.log` and
`/tmp/hybrid-relativistic-public-baseline113.log`.
The first probe stopped at a fixture error: default particle charge was +1,
while electromagnetic electrons require -1. The corrected probe sets -1 and
completed on both versions; its eight PIC results are explicitly `:accepted`.

## Repair and compatibility

Retain the finite-component and rounded-norm preflight, then compare the exact
stored squared components with exact `c²` when the accepted rounded norm lies
within eight machine epsilons of one. This repair rejects a truly invalid state
before the first particle or field is changed. It changes no public signature,
state representation, or checkpoint schema.

For hardware floats, the normal arithmetic has at most five accumulated
relative-rounding factors: the quotient contributes twice when squared, the
multiplication once, and the two positive additions once each. Bounding each
by one machine epsilon gives a lower bound of `(1-eps)^5 >= 1-5eps` for the
retained contribution. The eight-epsilon guard also leaves room for tiny
underflowed contributions; even Float16's minimum normal squared is less than
one epsilon by over five orders of magnitude. Six exact arithmetic checks of
the normal-error and conservative tiny-contribution bounds pass for the three
hardware types (`/tmp/hybrid-relativistic-roundoff-bound112.log`). This is a boundary predicate,
not a replacement numerical integrator.

Inputs whose rounded norm is already at least one remain rejected. Admitting
additional exactly subluminal states requires a separate repair of downstream
Lorentz-factor arithmetic, which still uses rounded squared norms. This slice
does not claim to solve those downstream range limitations or the open
[radiation arithmetic defect](RADIATION_REACTION_PLAN.md).

## Acceptance checks

- Exact-rational boundary witnesses, both sides of the boundary, all axis/sign
  permutations, three power-of-two light-speed scales and 1D–3D particles.
- Seeded nearby random states classified by an independent exact norm.
- Radiation rejection before cooling the earlier valid particle, including
  zero-force calls; PIC rejection for electron and mobile-ion inputs.
- Checkpoint rejection preserving prior file bytes and destination state;
  continued acceptance of a valid checkpoint.
- Both native Julia versions, neighboring forcing/particle tests, concrete
  compiler calls, package analysis against the unchanged gate, ordinary and
  exceptional allocation measurements, and final formatter/diff review.

Wrong patches distinguished by these checks: rejecting every near-boundary
state loses valid neighbors; checking the rounded `v/c` tuple exactly still
loses division bits at nonunit scales; checking only the first particle or
checking after mutation breaks public rollback assertions. Non-power-of-two
scale witnesses exercise division rounding independently and reject that second
wrong patch on all three hardware precisions.

## Observed verification

The frozen v1 source passes 12,706 focused assertions on each Julia version:
12,486 exact-boundary/permutation/random checks, 36 radiation rollback checks,
and 184 PIC/checkpoint transaction checks. The final v2 adds only two testsets:
15 non-power-of-two division checks and 80 BigFloat stored-precision/rounding
checks, all passing on each version. Final source and retained test text match
v1 exactly; a byte comparison confirms the two added testsets are the entire
v2 test delta. The complete new regression file therefore has 12,801 verified
assertions. Logs:
`/tmp/hybrid-relativistic-v1-focused112.log`,
`/tmp/hybrid-relativistic-v1-focused113.log`,
`/tmp/hybrid-relativistic-v2-added-regressions112.log`,
`/tmp/hybrid-relativistic-v2-added-regressions113.log`.
The last two logs also contain all 1,810 neighboring forcing, particle and
long-orbit checks; their total is 1,905 passing assertions each.

All 18 concrete compiler calls per version have zero reports (36 analysis and
runtime assertions). Package JET retains 26/25 reports against the unchanged
ceiling of 30, matching the baseline; see the existing
[MPI report classification](MPI_CHECKPOINT_PLAN.md). Logs:
`/tmp/hybrid-relativistic-v1-jet112.log` and
`/tmp/hybrid-relativistic-v1-jet113.log`.

All 27 warmed ordinary Float16/32/64, 1D–3D, 0/1/256-particle guard cases remain
unchanged and allocate 0 B on Julia 1.12.7. The three valid boundary-neighbor
cases allocate 3,168/3,168/3,248 B respectively, compared with the old guard's
0 B; exact boundary discrimination has an explicit resource cost. Each resource
run passes 57 assertions. Logs:
`/tmp/hybrid-relativistic-baseline-resources112.log` and
`/tmp/hybrid-relativistic-v1-resources112.log`.
No optimization or whole-project completion claim follows from these checks.
