# Reference-unit normalization range repair

Status: **NORM-001 V7 is pushed as `1097076`; local acceptance and all six
exact-revision CI jobs pass on Julia 1.12.7 and 1.13.2.**
This is a prerequisite for accepting the new local collision normalization
helper, and also affects existing public unit conversion and RayCon adapters.

## Evidence and canonical ownership

The keyword `PlasmaUnits` constructor accepts finite positive density, field
and mass references. [normalization.jl](../src/verification/normalization.jl)
owns `alfven_speed`, `gyrofrequency`, `inertial_length` and `_scale`, used by
`to_SI` / `to_normalized`. The RayCon normalized equilibrium and ray adapters
call these accessors in [raycon_normalized.jl](../src/raycon/raycon_normalized.jl).
The isolated collision candidate also derives its setup scales from them.

**VERIFIED:** the unchanged source at `6e7051b` produces 45 passing assertions
and 25 failures per runtime in a four-profile witness. Its independent 512-bit
reference combines the stored physical constants and inputs before conversion
to Float64. All three accessors, eight conversion kinds and representable
inverse conversions are checked. Ordinary solar-wind values are controls.

| References `(n0, B0, mi)` | Public result | Independent reference |
| --- | --- | --- |
| `(1e300, 1e200, 1e100)` | Alfvén speed `0.0` | `892.0620578335771` |
| `(1e-300, 1e-200, 1e-100)` | Alfvén speed `Inf` | `892.0620578335771` |
| `(1e300, 1e-300, 1e100)` | Inertial length `NaN` | `5.567813428950414e-79` |
| `(1e-300, 1e-200, 1e-100)` | Current scale `0.0` | Approximately `1.429241e-316` |

The checked arithmetic forms `mu0*n0*mi` before its square root (line 61),
then divides independently rounded Alfvén and gyro scales (line 75). `_scale`
composes those results further. These operations are the verified failing
owners; the complete-scale oracle and recorded outputs establish the finite
results that they lose. True output overflow is recorded separately.

Evidence is retained under
`HybridPlasmaPIC-audit-evidence/2026-10-23/normalization-range/`:
`source.json`, `witness.jl`, `baseline-commands.json`, both full logs/TOMLs and
`baseline-review.json`. The source hash is checked by each executed witness.

## Complete conversion evidence and implementation contract

The expanded baseline has **416 passes / 304 failures on each runtime**.
All 64 Float32/Float64 profiles have finite complete references but lose them
when pressure, electric, current or time scales are materialized first. Eight
BigFloat profiles additionally cancel enormous exponents analytically; both
conversion directions fail. Original logs, 72 records per runtime and source
hashes remain in `complete-baseline-*` under the evidence directory above.

A separate characterization passes 54 checks per runtime: 50 ordinary
input/direction records cover integer and rational promotion, complex numbers,
nonfinite values, signed zero, missing values, arrays, views, tuples, empty
containers and zero-dimensional broadcasting. Four checks establish that raw
constructors retain stored 512-bit BigFloat values at 24-bit working precision.

**VERIFIED fixture correction:** the first candidate passes 69 of the original
70 assertions. The remaining assertion incorrectly requires the inverse of a
rounded subnormal reference to equal one within normal-value epsilon. An
independent 200-digit Decimal computation gives `1.00000001413406` for that
complete inverse, agreeing with the candidate. `witness-v2.jl` uses an
independent complete inverse oracle with the original tolerance; the original
failed run and `subnormal-inverse-fixture-review.json` remain preserved.

The pre-edit contract and live execution checklist are retained in
`normalization-range/pre-edit-v1.md`. The canonical implementation is being
built in a separate checkout based on `f1410b7`:

- Generalize the existing exact square-root owner to a bounded rational with
  a separate Int128 exponent, explicit rounding and sign. Its existing
  two-argument collision entry retains nearest rounding. The first candidate
  passes **429,664 independent signed rounding-cell checks and all 86,988
  unchanged root checks per runtime**. This is numerical evidence, not final
  acceptance; compiler, resource and cost gates remain open.
- Express each complete scale's square as a rational monomial of the stored
  references. Combine the caller's value and exponent before final rounding
  on exceptional paths. Preserve the ordinary formulas and their output bits.
- Enforce finite positive values for all six unit fields in the inner
  constructor. Preserve valid same-type positional construction, typed
  conversion and stored BigFloat precision; retain the mixed-positional
  `MethodError`. Rejecting previously accepted invalid raw constants is an
  intentional compatibility change, pending verification and documentation.
- Treat signed zero and nonfinite values as multiplication/division by a
  positive physical scale. Preserve ordinary broadcast behavior. Verify the
  exceptional container/device path explicitly before acceptance.

The first normalization candidate passes the 720-check complete-conversion
witness, 54 compatibility checks, 143 existing normalization checks and Aqua's
nine checks on Julia 1.12. Its broader BigFloat run fails because `issubnormal`
has no BigFloat method. The failed source/logs are preserved. A second candidate
uses a separate MPFR endpoint guard. It passes all 720 complete-conversion,
54 compatibility, 143 existing normalization, 36 huge-exponent and nine Aqua
checks on each runtime. The corrected 70-check witness passes both; the
unchanged parent still fails 25 checks. Three directed-rounding cases fail
on the parent and pass on this second candidate.

Further independent checks keep acceptance open:

| Surface | Retained result | Required correction |
| --- | --- | --- |
| Broad CPU formulas and constructors | 28,019 passes / four errors per runtime | Preserve Julia's floating promotion for `π` and `ℯ`; they do not provide exact rational decomposition |
| Metal ordinary conversions | 48 checks pass on each candidate/runtime | Retain native device output type, input immutability and ordinary values |
| Metal complete exceptional conversions | Parent: 56 passes / four numerical failures; candidate: 48 passes / four GPU compilation errors | Keep Rational/BigInt work on the host and return corrected values to the original backend |
| Final value operation | 157 binary16 failures on Julia 1.12, despite an ordinary scale | Guard final multiplication/division too; an exactly representable Float64 oracle holds all three binary16 significands |

The initial Metal harness used a nonexistent status field; its failure logs
are retained separately. The table uses the corrected harness and actual native
runs. V3 addresses the verified gaps through the canonical conversion owner
and shared device staging, with thin CUDA/Metal adapters. V3 passes 28,027
broader CPU assertions, 3,800 final-operation range checks and 60 native Metal
conversion checks per runtime. Device layout, compiler, resource and cost
acceptance remain open.
Device inspection adds synchronization; its complete-call resource and timing
costs must be measured before acceptance. Native CUDA execution is not verified.
No candidate is promoted on the strength of partial results.

## Current candidate: V7

The V7 permanent run passes 99,481 assertions per runtime, plus the isolated
child suites: 429,664 signed square-root cells,
40,500 endpoint-guard checks, 1,100 working-precision checks, three directed
endpoint cases and 16 huge-exponent conversions. The child also checks restoration
of numerical settings. Native Metal coverage includes 624 singleton boundary
assertions and 288 storage/layout assertions; empty, zero-dimensional, strided
and transposed arrays are included. A further 784 flag boundary/lifetime checks
and 96 complex-device checks pass per runtime. The final Metal-only rerun passes
all 1,792 assertions after adding the complex cases to the permanent suite.
The original normalization checks remain.

| Additional validation | Verified result on Julia 1.12.7 and 1.13.2 |
| --- | --- |
| Neighboring arithmetic, collisions, ionization and RayCon | 141,936 passes each on V5; V6/V7 changes only normalization guards, special outputs and device inspection |
| Four-thread numerical-context preservation | 6,525 passes each, with each worker's MPFR exponent limits set and restored explicitly |
| Exceptional resource profiles | 720 complete conversions each, 2,880 checks, zero measured compilation; maximum 59,200 bytes against the frozen significand/precision bound |
| Concrete compiler analysis | All 60 signatures have zero reports each |
| Package quality/compiler | Aqua and JET gates pass; full 29/27 report trees equal the reviewed parent after exact unchanged-line mapping |
| CUDA 6.4.1 interface and joint extension loading | Ten routing checks and one ambiguity check each; native NVIDIA execution unavailable |
| Final ordinary timing/allocation profiles | V2 passes all 306 profiles and six controls per runtime; all measured compilation is zero and all predeclared resource limits pass. V1 failures remain retained |

Earlier failed candidates remain evidence. V3's expanded layout checks exposed
wrapper dispatch and zero-dimensional return changes; V4 fixes these. V4's
256-bit, 128-element BigFloat conversion allocated 84,096 bytes versus the
parent's 16,688. V5 preserves the exact endpoint predicate using an exponent
check, reducing that profile to 15,680 bytes. All 30 ordinary allocation profiles
per runtime retained complete parent outputs and zero measured compilation.
These are allocation results, not latency measurements.

V5 singleton device checks then retained 608 passes/16 failures per runtime.
A native predicate probe verifies that Float32 subnormals survive transfer but
numeric zero comparison classifies them as zero. V6 inspects stored IEEE bits
before deciding whether a lost finite input needs complete conversion. A second
V5 witness retained 400 passes/200 failures: exceptional BigFloat special values
kept stored 512-bit precision. V6 constructs arithmetic outputs at working
precision while leaving constructor reference precision intact. Both original
regressions now pass, alongside negative, nonfinite and complex controls.

Full command records, source snapshots, complete compiler trees and results are
in `normalization-range/*-v7-*` and `compiler-v7-review.json`. The resource
protocol is `resources-v6.md`; the original failures are never overwritten.
The failed V1 fixed run retains all samples: its six ratio controls pass, but
12 zero-compilation assertions fail before any candidate profiles run.
**VERIFIED:** native traces identify Metal's shape-specialized broadcast
transition after ten calls, beyond the original three warmups. Separate
component probes reproduce recurring inference in `Metal.return_type` used by
`mapreduce`; changing warmup alone cannot remove that per-call query. Dependency
sources and full probes are retained as `gpu-compile-counter-v2-*` and `-v3-*`.

V7 replaces that Metal inspection reduction with an initialized, per-call
atomic flag using the same numerical predicate. It covers empty inputs,
thread-group boundaries, concurrent flag writes and lifetime through
synchronization. All 60 warmed compilation observations per runtime are zero.
CUDA retains the generic reduction. This change does not promote the separate
held gathered-push implementation.

The frozen V2 run passes all 210 CPU and 96 Metal profiles on each runtime,
plus six identical/doubled-work controls. Every calibration and paired sample
has zero measured compilation. Machine CPU calls add no bytes; BigFloat public
conversions stay at or below parent allocation, while the three accessors use
at most the predeclared one-temporary allowance (80/96/112 bytes at 24/128/256
bits). Metal adds at most 1,716 / 1,682 host bytes. No failed sample was removed,
no unchanged run retried, and no threshold relaxed.

Measured ratios of median complete-call times against the unchanged parent:

| Profile group | Julia 1.12.7 ratio range | Julia 1.13.2 ratio range |
| --- | --- | --- |
| CPU Float16 | 1.1398–4.0971 | 1.1085–4.3082 |
| CPU Float32 | 1.2252–3.6217 | 1.1846–3.4115 |
| CPU Float64 | 1.2967–2.5934 | 1.1177–2.6213 |
| CPU BigFloat p=24 | 0.1235–1.4393 | 0.1230–1.6082 |
| CPU BigFloat p=128 | 0.0886–1.3157 | 0.1061–1.3199 |
| CPU BigFloat p=256 | 0.0854–1.3394 | 0.0818–1.2406 |
| Native Metal, all 96 profiles | 1.4830–2.6198 | 1.7813–2.4828 |

These ratios report local costs, including synchronization for Metal. The
predeclared correctness-repair protocol has validity and allocation limits,
with no latency-improvement claim. The measured CPU and device regressions are
open optimization work. All samples and complete outputs remain in
`costs-v2-*.toml` and `.outputs`; `cost-v2-review.json` checks every frozen gate.
A separate output review compares values, types and stored precision.

The final source/caller review finds no additional confirmed defect in this
bounded conversion surface. API docstring edits preserve the executable AST;
formatting and relative documentation links pass.
[All six exact-revision CI jobs pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40042288497).
Complete logs verify 723,136 / 723,137 assertions on Julia 1.12.7 / 1.13.2.
Full compiler trees match reviewed `ba3694c` after checked unchanged-line maps
(24 / 21 reports), with no suppressed diagnostics. CAP-001 and broader GPU
evolution remain open.

RayCon has additional compound temperature, frequency and amplitude conversions
in its adapter. Repairing these shared accessors alone does not establish the
range correctness of every RayCon expression; that integration remains open.

## Repair and acceptance tracker

- [x] Trace accessors, public scalar/array conversions, RayCon callers and tests.
- [x] Reproduce complete-scale failures against independent references on both
  supported Julia versions; preserve the ordinary controls.
- [x] Characterize raw positional constructors, stored floating types and invalid
  constants: 112 records per runtime retain the mixed-positional `MethodError`
  and show that raw typed constructors accept zero, negative and nonfinite
  constants. Eight huge-exponent BigFloat profiles retain 20 passes / 16 failures.
  Precision, rounding and exponent limits are restored after characterization.
- [x] Record constructor/zero-sign compatibility before the isolated edit.
  Local implementation, validation and exact-revision CI pass. The
  [shared BigFloat materializer repair](EXACT_BIGFLOAT_CONVERSION_PLAN.md)
  now passes local acceptance and all six exact-SHA CI jobs with full compiler
  review. Both repairs retain their separate evidence and scope.
- [x] Implement the canonical complete scale in the isolated candidate, retaining ordinary outputs and types
  where the existing evaluation is safe. Reuse exact conversion/square-root
  ownership where applicable. Bound exceptional work; avoid expanding enormous
  BigFloat exponents into unbounded rational integers or changing global MPFR
  state. Do not patch only the collision normalization helper.
- [x] Reproduce full-value conversion failures separately from scale failures.
- [x] Complete scalar/container/device compatibility checks for
  the public conversion implementation, including irrational promotion, final
  output endpoints, GPU staging, empty arrays, views and zero-dimensional arrays.
- [x] Add independent boundary and exponent matrices, subnormal outputs,
  forward/inverse scalar and array checks, and retained-precision inputs.
  Reject zero-clamping, independent-factor rounding and round-trip-only patches
  with explicit negative controls.
- [x] Run the corrected independent witness, existing normalization/RayCon tests,
  concrete/full compiler comparisons and threaded state checks on both versions.
  The original failed inverse-subnormal fixture remains recorded; its oracle
  correction and independent decimal check are documented above.
- [ ] Run collision unit-equivalence tests on the final collision source; this
  remains part of CAP-001 acceptance, not evidence supplied by NORM-001.
- [x] Freeze and pass exceptional resource profiles; retain every failed ordinary
  measurement and keep the original validity/resource thresholds.
- [x] Execute and review all final ordinary complete-call cost profiles.
- [x] Update API docs, README and changelog; review final source, callers and diff.
- [x] Integrate this repair as a separate change on main.
- [x] Require all six exact-SHA CI jobs for this repair; full compiler reports reviewed.

## Measured optimization follow-up

- [ ] Profile CPU scalar/vector paths against this accepted repair. The worst
  observed Float16 profile is over four times the old parent. Test whether
  separating array conversion and range inspection reduces complete-call cost;
  this is an untested design hypothesis, not an attributed cause.
- [ ] Profile native Metal conversion and inspection. Test a single fused device
  kernel against the two-stage implementation, preserving the exact predicate,
  signed/subnormal behavior, storage/layout, per-call ownership and host fallback.
- [ ] Before each implementation's measurements, freeze all 306 profiles, both
  runtime lanes, controls and allocation limits. Require no profile above 1.25
  times this repair's cost, and at least a 20% reduction in the specifically
  targeted worst profiles. Retain every failed candidate. These are future
  optimization targets; they do not retroactively change V2's declared gates.

This plan does not establish correctness of the unaccepted collision candidate
or close CAP-001. Its corrected unit-equivalence and coupled-channel tests still
require execution on the final source.
