# Linear Crank–Nicolson range repair

Status: V5c landed on `main` after `4fafe4a` (2026-10-26). In its landing verification,
`test_linear_cn_range.jl` (16,762 assertions), `test_phase_final.jl`,
`test_exact_arithmetic.jl`, `test_shock_cn.jl`, `test_dispersion.jl`, Aqua, JET
(24 / 21, parent signatures apart from line shifts) and validation case 19 pass
on Julia 1.12.7 and 1.13.2. Ordinary Float64 `cn_multiplier` costs 15.9 / 13.9
ns per call versus 12.9 / 11.9 ns on the parent (1.23× / 1.17×, measured on one
shared Apple M5 core); this reference helper is not a solver hot path.
`03e4eb3` passes all six [exact-SHA CI jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40275287900); complete compiler review is recorded in the [audit tracker](DEEP_DEBUG_TRACKER.md#current-progress).
This is the linear-reference
part of CAP-005A. The accepted
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
baseline helper at line 29: finite floating products overflow and fixed-width integer
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
16 checked components. The candidate computes rational-only results exactly,
retaining the promoted legacy integer type when the result fits and otherwise
widening to `BigInt`.

An independent near-quarter-turn witness records another 12 passes / 12 failures
per runtime. `nextfloat(one(T))*prevfloat(T(2))`, with both argument orders and
signs for Float16/32/64, loses a nonzero real component when the intermediate
product rounds to two. An overflow-only repair would miss this defect.

## Implemented compatibility contract

The completed V2 characterization records Float16/32/64 outputs, promoted
mixed floating inputs, fixed-width integer outputs as ComplexF64, BigInt outputs
as Complex{BigFloat}, and exact rational outputs with their rational element
type. Positive and negative zero preserve the phase sign. Direct nonfinite
inputs currently propagate NaN; the driver already rejects nonfinite inputs.
The candidate preserves that direct-helper policy.
The first characterization fixture had a metadata dictionary-type error; its
records are retained as invalid evidence, separately from the corrected V2.

The baseline loses the smallest positive Float64 phase to signed zero. The
candidate retains representable phase components as well as unit modulus. Checking
only energy would also miss the integer witness, whose wrong result has norm
one. Ordinary complete outputs and return types must be compared to the parent;
any intentional correction beyond the original failures needs its own oracle.

The existing exact conversion owner is
[exact_arithmetic.jl](../src/utils/exact_arithmetic.jl); scaled product helpers
also exist in [collisions.jl](../src/particles/collisions.jl). New scaled exact
comparisons and materialization extend the shared conversion owner. Hardware
hot paths retain native arithmetic; generic exact evaluation covers mixed
signed/unsigned integers, rational overflow and stored BigFloat precision.
Finite quotients do not expand canceling absolute binary exponents into integers.

## Candidate implementation

The canonical helper retains ordinary hardware arithmetic outside cancellation
and range-loss regions. A product residual corrects the near-quadrature real
component. Exact exceptional evaluation stores bounded significands and a
separate binary exponent; shared exact comparisons round the complete component
without allocating integers proportional to the absolute exponent. Production
code changes no global MPFR precision, rounding mode or exponent limit.

Ordinary BigFloat inputs form `x=ω*dt`, `x²`, and `4±x²` exactly in private
MPFR values whose individual precision follows the stored input widths and
product exponent. Only the two final component divisions round. A conservative
guard proves that every intermediate fits the current exponent interval and
bounds temporary precision by input width; remote exponents or insufficient
range use the separate-significand path. Input values are never mutated.

Pure rational calls retain exact results, widening only when needed. Mixed
operands enter the hardware path only when conversion preserves both values.
Direct nonfinite inputs retain their prior expression behavior and signed-zero
phases retain their sign. The complete frequency builder now receives `B0` and
`n0` separately; an unrepresentable finite Float64 frequency raises
`ArgumentError` before random packet construction. Float16/32 coefficients
are divided at the Float64 output precision. Validated grid counts use signed
mode indexing so unsigned counts preserve negative band-limited modes. Integer
and mixed physical parameters use the hardware frequency path when their
conversion to Float64 is exact; this includes the three-argument frequency
wrapper's unit density.

The unchanged retained matrices pass on Julia 1.12.7 and 1.13.2: 2,808 expanded
assertions, 24 near-quadrature assertions, 148 rational/huge-exponent assertions
and the original 53 package assertions per runtime. The registered suite passes
16,762 parent assertions on each version, including 550 complete-frequency and
trajectory checks; its isolated child passes 115,300 compact rounding-cell and
300 huge-exponent assertions, plus 15,570 admission-boundary, stored-input and
rounding-cell checks, with settings restored. Existing exact-conversion
checks pass 49,083 assertions and the original CN phase/caller testset passes 20.

Both versions pass 20 concrete JET checks and nine Aqua checks. Complete compiler
trees match the clean parent after checked unchanged-line mappings (29 reports
on 1.12.7; 27 on 1.13.2), with no unmapped frames or suppressed findings.
Twelve configured formatter/syntax checks pass. The external evidence
directory `linear-cn-range/repair-v1` preserves the V4 source freeze, commands,
logs, compiler comparisons and failed fixtures. V5 revalidation is retained
separately in `linear-cn-range/repair-v5`.

**VERIFIED resource failure and repair:** the unchanged ordinary allocation
preflight failed 21 of 24 V4 profiles on each runtime. At precision 256,
`cn_multiplier(BigFloat(3),BigFloat(2))` allocated 590,528 bytes versus 3,360
parent bytes. Integer and mixed coefficient packets allocated about 2.5 MB
versus 23,984/24,112 parent bytes. Allocation stacks identify per-bit exact
comparison loops and per-mode rational frequency construction as the owners.
V5 passes all 24 unchanged profiles on both versions: ordinary BigFloat
precision 24/53/128/256 uses 464/528/672/928 bytes, and every ordinary packet
matches the parent exactly. Hardware scalar calls remain allocation-free.
The fixed latency contract and exceptional-path cost reports remain open;
functional test durations do not establish latency acceptance.

## Implementation and acceptance tracker

- [x] Map canonical helper, frequency builder, both public callers and tests.
- [x] Reproduce the scalar and complete-packet failures on both Julia versions.
- [x] Characterize ordinary public return types, zero signs and NaN propagation.
- [x] Extend baseline probes to promoted fixed-width integers, stored BigFloat
  precision/compact exponent boundaries and complete frequency construction.
- [x] Characterize rational inputs and very large BigFloat exponents, including
  explicit exception records. Reproduce near-quarter-turn component loss.
- [x] Fix the rational compatibility contract before editing. The
  [shared BigFloat materializer repair](EXACT_BIGFLOAT_CONVERSION_PLAN.md)
  now passes all six exact-SHA CI jobs and full compiler review. Its
  bounded-rational owner now also supports separately scaled conversion and
  exact component comparisons for the multiplier.
- [x] Implement the repair at the canonical multiplier. Preserve ordinary
  arithmetic, finite phase signs and representable subnormal components.
- [x] Audit `_whistler_omega` separately: coefficient division, `k²`, zero mode
  and finite output representability. A repaired multiplier alone does not
  establish range-correct frequency construction.
- [x] Add exact-rational and independent high-precision regression oracles,
  exponent boundaries, swapped inputs, sign/zero controls, complete trajectories
  and both public callers. Run original witnesses unchanged on both runtimes.
- [ ] Run the frozen scalar and whole-packet cost contract. Require
  ordinary complete outputs/types and warmed allocations to match; retain
  all cold-path costs. Define slowdown controls and fixed latency gates before
  measurement, then run once per frozen candidate and preserve failures.
- [x] Preserve the V4 allocation-gate failure, profile its canonical owners,
  implement the bounded V5 repair, and pass the unchanged 24-profile allocation
  preflight on both runtimes. Keep the original latency ceilings unchanged.
- [x] Review concrete/full compiler trees, neighboring integration, resource
  lifetimes and a fresh final source pass. No suppressed diagnostics.
- [x] Correct the mathematical description: the numerator and denominator are
  degree one, so this is the [1/1] Padé approximant; its expansion agrees with
  `exp(-z)` through `z²`, with local error `O(z³)`.
- [ ] Synchronize README, changelog and capability/audit trackers; commit and
  push this repair separately to main, then require all six exact-SHA CI jobs.

No correctness, performance or full-simulator completion claim follows from
this plan. The unresolved capabilities remain in the
[implementation roadmap](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md).

## Current implementation and allocation finding (2026-10-25)

The isolated V4 implementation passes 16,534 registered parent assertions and
115,600 isolated BigFloat child assertions per runtime, plus the retained
witnesses and neighboring exact-arithmetic checks. Concrete compiler/Aqua gates
pass; complete compiler trees match the parent after checked unchanged-line
mapping (29/27 reports). These checks do not establish performance acceptance.

Root's next ordinary allocation preflight fails 21 of 24 profiles on **each**
Julia 1.12.7/1.13.2. Ten alternating same-input measurements have identical byte
counts per profile. At precision 256, `cn_multiplier(BigFloat(3), BigFloat(2))`
allocates 590,528 bytes versus 3,360 in the clean parent. Integer/mixed coefficient
packets at n = 128 allocate about 2.5–2.62 MB versus roughly 24 KB. Ordinary Float64
packets retain exact byte equality. The underlying allocation paths require
profiling; these measurements alone do not assign every allocation to an owner.

V4 and its failures remain frozen. Canonical resource optimization is active
before any timing campaign or promotion. Exact component/rounding-cell tests,
compact and huge exponent behavior, return types and unchanged resource limits
remain required. Evidence is in
`linear-cn-range/repair-v1/root-resource-preflight-v1-{112,1132}.{log,toml}` and
`root-resource-preflight-v1-commands.json` under the retained evidence directory.
