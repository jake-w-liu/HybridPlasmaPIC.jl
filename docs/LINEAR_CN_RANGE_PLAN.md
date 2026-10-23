# Linear Crank–Nicolson range repair

Status: confirmed baseline defect; no implementation accepted. This is the
remaining linear-reference part of CAP-005A. The accepted
[shock field solve](SHOCK_CN_SOLVE_PLAN.md) is a separate owner. Nonlinear
particle–field coupling and general implicit solvers remain CAP-005B/C.

## Checked behavior and ownership

`cn_multiplier` in [semi_implicit.jl](../src/integrators/semi_implicit.jl)
implements `(1 - im*ω*dt/2)/(1 + im*ω*dt/2)`. `run_whistler` applies it to each
Fourier amplitude; `compare_integrators_whistler` calls that driver. The
existing [phase-final tests](../test/test_phase_final.jl) cover ordinary phase,
stiff energy conservation and invalid driver configuration.

**VERIFIED:** both Julia 1.12.7 and 1.13.2 retain 18 control passes and 35 failed
assertions against the same source as `d6bce16`. The witness computes real and
imaginary components independently from `Rational{BigInt}` inputs:
`h=ω*dt/2`, `Re=(1-h²)/(1+h²)`, `Im=-2h/(1+h²)`.

| Input | Observed | Independent reference |
| --- | --- | --- |
| `floatmax(Float64), 4.0` | NaN components | `-1 - 5.562684646268003e-309im` |
| `floatmax(Float32), 4f0` | NaN components | `-1 - 2.938736e-39im` in Float32 |
| `typemax(Int), 2` | `0 + 1im` | Approximately `-1 - 2.168404344971009e-19im` |
| Four-mode packet, finite maximum Float64 step | NaN energy ratio | Unit modulus for each CN multiplier |

The first three failures follow the intermediate product in the canonical
helper at line 29: finite floating products overflow and fixed-width integer
multiplication wraps before division. Source hashes, executed commands and all
20 recorded profiles per runtime are retained under the persistent
`HybridPlasmaPIC-audit-evidence/2026-10-23/linear-cn-range/` directory in
`baseline-commands.json`, `baseline-*.toml` and their logs.

## Expanded baseline

Both Julia versions retain **2,571 passes and 237 failures** across 686 recorded
profiles in `extended-baseline-commands.json` and `extended-baseline-*.toml`.
The binary/integer matrix has 2,316 passes / 180 failures, the isolated compact
BigFloat exponent matrix 242 / 40, and complete frequency construction 13 / 17.
There are no test errors. The BigFloat fixture restores exponent limits,
precision and rounding settings; its exact oracle reads stored input bits.

**VERIFIED:** frequency construction also loses finite complete frequencies.
The six recorded cases cover ordinary construction, coefficient overflow and
underflow, squared-wavenumber overflow and underflow, and zero guide field.
The oracle combines stored `Float64(2π)`, mode index, guide field, density and
length as one rational expression. Repairing `cn_multiplier` alone cannot close
these recorded frequency failures.

The subsequent rational/large-exponent matrix records 132 passes / 16 failures
and 79 profiles per runtime. It retains 11 `OverflowError` and two `InexactError`
observations for rational inputs rather than counting them as successful phase
results. Large canceling BigFloat exponents pass; overflowing products lose
16 checked components. The compatibility choice for narrow rational types
remains open.

An independent near-quarter-turn witness records another 12 passes / 12 failures
per runtime. `nextfloat(one(T))*prevfloat(T(2))`, with both argument orders and
signs for Float16/32/64, loses a nonzero real component when the intermediate
product rounds to two. An overflow-only repair would miss this defect.

## Compatibility decisions to verify before implementation

The completed V2 characterization records Float16/32/64 outputs, promoted
mixed floating inputs, fixed-width integer outputs as ComplexF64, BigInt outputs
as Complex{BigFloat}, and exact rational outputs with their rational element
type. Positive and negative zero preserve the phase sign. Direct nonfinite
inputs currently propagate NaN; the driver already rejects nonfinite inputs.
Changing that direct-helper policy requires an explicit documented decision.
The first characterization fixture had a metadata dictionary-type error; its
records are retained as invalid evidence, separately from the corrected V2.

The smallest positive Float64 phase currently becomes signed zero. The repair
must retain representable phase components as well as unit modulus. Checking
only energy would also miss the integer witness, whose wrong result has norm
one. Ordinary complete outputs and return types must be compared to the parent;
any intentional correction beyond the original failures needs its own oracle.

The existing exact conversion owner is
[exact_arithmetic.jl](../src/utils/exact_arithmetic.jl); scaled product helpers
also exist in [collisions.jl](../src/particles/collisions.jl). Reuse these where
their contracts apply. Do not add unconditional multiprecision arithmetic to
ordinary floating calls. Before selecting a fallback, check mixed signed and
unsigned integers, rational overflow behavior, BigFloat stored precision and
compact/large MPFR exponent ranges. A finite quotient must not require an
unbounded rational expansion of canceling input exponents.

## Implementation and acceptance tracker

- [x] Map canonical helper, frequency builder, both public callers and tests.
- [x] Reproduce the scalar and complete-packet failures on both Julia versions.
- [x] Characterize ordinary public return types, zero signs and NaN propagation.
- [x] Extend baseline probes to promoted fixed-width integers, stored BigFloat
  precision/compact exponent boundaries and complete frequency construction.
- [x] Characterize rational inputs and very large BigFloat exponents, including
  explicit exception records. Reproduce near-quarter-turn component loss.
- [ ] Fix the rational compatibility contract before editing. Complete the
  [shared BigFloat materializer repair](EXACT_BIGFLOAT_CONVERSION_PLAN.md)
  through its remaining exact-SHA CI gate before accepting a multiplier that
  relies on it. Its bounded-rational contract does not yet cover arbitrary
  separate binary exponents.
- [ ] Implement the repair at the canonical multiplier. Preserve ordinary
  arithmetic, finite phase signs and representable subnormal components.
- [ ] Audit `_whistler_omega` separately: coefficient division, `k²`, zero mode
  and finite output representability. A repaired multiplier alone does not
  establish range-correct frequency construction.
- [ ] Add exact-rational and independent high-precision regression oracles,
  exponent boundaries, swapped inputs, sign/zero controls, complete trajectories
  and both public callers. Run original witnesses unchanged on both runtimes.
- [ ] Freeze scalar and whole-packet cost profiles before measuring. Require
  ordinary complete outputs/types and warmed allocations to match; retain
  all cold-path costs. Define slowdown controls and fixed latency gates before
  measurement, then run once per frozen candidate and preserve failures.
- [ ] Review concrete/full compiler trees, neighboring integration, resource
  lifetimes and a fresh final source pass. No suppressed diagnostics.
- [ ] Correct the mathematical description: the numerator and denominator are
  degree one, so this is the [1/1] Padé approximant; its expansion agrees with
  `exp(-z)` through `z²`, with local error `O(z³)`.
- [ ] Synchronize README, changelog and capability/audit trackers; commit and
  push this repair separately to main, then require all six exact-SHA CI jobs.

No correctness, performance or full-simulator completion claim follows from
this plan. The unresolved capabilities remain in the
[implementation roadmap](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md).
