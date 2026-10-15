# Complete electromagnetic charge combination

Status: confirmed follow-up to the isolated MOMENT-001 per-species repair.
Implementation and final acceptance remain open. This plan covers charge density;
combined current, hybrid bulk-velocity ratios, pressure and MPI reductions need
separate repairs and verification.

The isolated `/tmp/hybrid-empic-charge-combination-v1` prototype retains private
normalization status and routes both model families through one charge owner.
Its first Julia 1.13.1 focused run passes **56 assertions**, exit 0, covering
overflow/cancellation, final scratch contents and both subnormal boundary cases.
Log: `/tmp/hybrid-empic-charge-combination-v1-focused113.log`. It remains outside
the shared source; broader compatibility, integration and resource gates are
still required.

## Verified boundaries

With finite inputs, both `EMPIC1D` and `EMPIC` still combine already-rounded
species densities. On the final MOMENT-001 candidate, the retained 22-assertion
probe reproduces opposite overflowing species densities/currents becoming NaN
instead of exact zero, and background subtraction becoming −Inf instead of
finite `-floatmax(T)/2`. The Julia 1.12.7 replay exits 0:
`/tmp/hybrid-moment-combination-baseline112.log`. These are passing defect
witnesses, not successful repaired outputs.

Two independent boundary controls pass **28 witness assertions on each native
Julia version**, both exit 0. Let `u` be the smallest positive subnormal and
the cell width be 8:

- An ion of weight `3u` and an electron of weight `-3u` each contributes `3u/8`
  charge density. Separately rounded densities vanish. Their complete sum is
  `3u/4`, which rounds to `u`; the implementation returns zero.
- A background of `u` and an electron of weight `4u` have complete charge
  density `u/2`, a tie that rounds to even zero. Rounding electron density first
  gives zero, so the implementation incorrectly returns `u`.

Both controls cover Float32/64 and both EMPIC interfaces. Rational identities
verify the exact fractions independently of production rounding helpers.
Script: `/tmp/hybrid-charge-combination-boundaries.jl`; logs:
`/tmp/hybrid-charge-combination-boundaries{112,113}.log`.
They show why detecting only nonfinite combined outputs is insufficient.

## Pre-edit brief

- Mode: repair complete charge reduction at the electromagnetic model boundary.
- Acceptance: when finite represented inputs lose range during deposition,
  normalization, species scaling or combination, round the complete cell charge
  expression once. Preserve ordinary order, visible nonfinite propagation,
  particle state and final electron-density scratch contents.
- Baseline: the overflow/cancellation and underflow/midpoint witnesses above run
  against the final isolated MOMENT-001 source, so correcting per-species
  normalization alone does not close this defect.
- Cause: **VERIFIED**, both `_charge_density!` implementations in
  `src/models/full_pic.jl` call `density!` before combining species or subtracting
  from `n0`. The exact contributions in the witnesses are lost at that boundary.
- Invariants: mobile charge is `q_i*n_i + q_e*n_e`; immobile charge is `n0-n_e`.
  Public EMPIC validators require electron/ion charges −1/+1. Retain ion-then-
  electron ordinary order and the existing immobile expression. Empty species,
  signed weights, deterministic shapes and writable particle views remain valid.
- Owners and surfaces: both model-family `_charge_density!` entry points and
  their species-add helpers; the shared private density/normalization owner
  carries range/eligibility status without changing `density!`'s public return.
  Initialization and each step already call these charge owners. Reuse the
  per-species exact moment accumulator and common final-rounding owner.
- Compatibility: no new public diagnostic, field, checkpoint schema or species
  model. Preserve ordinary arithmetic and nonfinite/unsupported-coefficient
  behavior. Finite exceptional combined outputs intentionally change.
- Strategy: retain normalization status through the private density path; one
  shared model charge owner combines those statuses with coefficient/addition
  checks. Eligible range loss evaluates all included species and the background
  in rational cell expressions before one conversion. Retain the normal scratch
  writes, and never infer eligibility from rounded output alone.
- Verification: exact rational rounding cells, overflow and subnormal/midpoint
  controls, nonfinite/unsupported coefficients, zero/unequal populations, all
  shapes and dimensions, ordinary bitwise/scratch comparisons, initialization
  and neighboring step tests, concrete/package JET, Aqua, allocations/timing,
  both native versions, formatter, exact diff and GitHub CI.

The cold path is bounded to one complete combination replay per call and reuses
the deposition arithmetic owner. Its work is proportional to particles/stencils
and grid cells, with operand-dependent BigInt storage. Measure its actual costs;
do not claim a fixed memory or latency bound.

Three incomplete repairs must be rejected: testing only Inf/NaN misses the
subnormal controls; adding rounded exact per-species values still loses the
midpoint case; making all calls exact changes ordinary rounding and cost.

- [x] Reproduce post-MOMENT-001 overflow/cancellation witnesses.
- [x] Verify independent underflow/midpoint witnesses on both Julia versions.
- [x] Trace both canonical charge owners and their initialization/step callers.
- [ ] Implement the shared private status and complete-combination path.
- [ ] Pass focused, neighboring, compiler, quality and resource gates.
- [ ] Commit/push the verified repair and observe GitHub acceptance.
