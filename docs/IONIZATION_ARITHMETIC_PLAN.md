# Ionization energy and optical-depth arithmetic

Status: COLL-001 confirmed on Julia 1.12.7 and 1.13.2; production repair is
not yet implemented. This is separate from the GROWTH-001 storage preflight
and GROWTH-002 migration repairs.

## Verified baseline

The frozen reference source at
`/tmp/hybrid-rank-allocation-v2-20261020/source` matches pushed `b3e2eb6`
across all 92 source/extension files. `src/particles/collisions.jl:1029–1042`
specifies strict eligibility `m|v|²/2 > E_iz`, probability
`1-exp(-nσ_iz|v|dt)`, and primary cooling by `E_iz`.

**VERIFIED:** `/tmp/hybrid-ionization-energy-range-review-20261020/probe.jl`
passes 132 assertions and fails 30 on each runtime, with no fixture errors.
The matrix covers Float32/64 and 1D–3D; every velocity component, mass and
energy input is finite. Exact rational energy ratios and a serial 512-bit
square-root oracle establish eligibility and finite cooled components.
Ordinary controls pass. When both represented speed and threshold are
infinite, six eligible cases return zero births. In six further cases, the
birth occurs but the primary retains every original velocity component.
Complete records, failure-line counts and all log/source hashes are retained
in `energy-summary.json`.

**VERIFIED cause:** the owner first collapses the norm to `T` at line 1101
and compares it with the collapsed threshold from line 1097. The recorded
`Inf > Inf` comparison skips eligible particles. At lines 1168–1176, a
finite threshold divided by an infinite norm produces a zero loss ratio and
unit cooling factor. The recorded unchanged primary values reproduce this
path.

**VERIFIED:** the independent `probability.jl` probe passes 144 assertions
and fails 18 per runtime. With three `floatmax(T)` components and rate
`floatmin(T)/128`, exact arithmetic bounds the finite optical depth below
`1/16`. A uniform sample of `1/2` must reject the event, but all six extreme
cases create a pair and advance both counters. Ordinary rare-event and
saturated controls pass. `probability-summary.json` retains all records and
hashes. The source passes the already infinite norm to
`_collision_probability` at line 1103; its finite-factor reconstruction at
lines 67–90 cannot recover the discarded speed.

The first additional threshold probe used an ambiguous parametric RNG
overload. Its failures are retained and are not accepted evidence. The revised
probe at `/tmp/hybrid-ionization-boundary-review-v2-20261021` uses explicit
Float32/64 methods with native RNG controls, preserving all behavior assertions.
**VERIFIED:** both runtimes retain 77 passes and 19 behavior failures, with
every error result `nothing`; `summary.json` checks all hashes, failure lines
and twelve records per runtime. Float32 creates a birth at exact equality
`v=(1,0,0), m=1, E_iz=0.5`; Float64 skips a birth for `nextfloat(1)`. Tiny
orthogonal components expose additional lost eligibility/cooling, including
cases where a completed primary component is representable while a standalone
cooling scale would underflow. Below-threshold controls pass.

## Repair scope and acceptance

Repair the existing ionization owner, preserving the public signature,
validation order, newborn sampling order, ID semantics and storage guards.
Retain complete finite products and compare the strict energy threshold
before narrowing range. Cooling must preserve representable completed
components even when an intermediate speed or scale is unrepresentable.
Keep ordinary cases allocation-free at the scalar arithmetic layer and
measure full-call cost against the accepted parent.

The neutral-MCC caller checks that its relative speed is finite before using
the shared probability helper (`collisions.jl:976–980`). Do not silently
change its documented range-error behavior as part of this repair.

Before production edits, complete a source-backed design brief. Check existing exact-arithmetic helpers before
adding a primitive. Avoid BigFloat transcendental calls in concurrent
production paths while NUM-001 remains open.

Required gates:

- Preserve all original failed probes and verify the repaired public call.
- Cover equality, immediately adjacent energies, tiny orthogonal components,
  subnormal/overflow intermediates, sign/permutation symmetry and empty/no-event
  cases. Verify complete populations and counters.
- Check error paths, selected-batch validation and storage rejection together
  with the revised GROWTH-001 guard.
- Verify any new arithmetic primitive independently with exact midpoint and
  range boundaries, then threaded controls on both Julia versions.
- Run collision, particle-growth and integration regressions; concrete and
  complete package compiler comparisons; formatting; full-output and resource
  comparisons. Record any behavior or allocation change explicitly.
- Promote this repair separately after its parent and require all six GitHub
  jobs green for the exact commit.
