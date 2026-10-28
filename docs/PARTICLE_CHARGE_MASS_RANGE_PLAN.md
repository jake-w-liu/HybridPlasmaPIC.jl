# Particle charge/mass factor range

Status: confirmed current-main defect. The isolated P1 CPU adapter preserves
original factors only for nonnormal quotients and passes a 964-assertion
source/AST preflight. Numerical, compiler, resource and integration acceptance
remain open. A new isolated production candidate is being implemented; its
trace/carrier/certificate source has no numerical or performance acceptance yet.
Bounded scalar and whole-batch experiments are recorded separately.
This work is separate from C1's supported relativistic output closure.
The numerical contract uses the original stored q, m, dt and field operands;
a rounded derived ratio must not silently erase a representable result.

## Verified reproductions

VERIFIED on current main `cd8b94d`, actual Julia 1.12.7 and 1.13.3:
public `push_uniform!` and `push_gathered!` each return normally with zero
velocity and coordinates for the following nonzero-charge inputs. Each uses
four particles at rest in 3D, B=0, E=(0,0,dt), both charge signs, and independent
gathered midpoint storage. All exact kick, position and midpoint results are
representable dyadics.

| Type | Charge magnitude | Mass | dt and electric amplitude | Exact kick magnitude |
| --- | --- | --- | --- | --- |
| Float16 | 2^-24 | 2^10 | 2^10 | 2^-14 |
| Float32 | 2^-149 | 2^80 | 2^60 | 2^-109 |
| Float64 | 2^-1074 | 2^600 | 2^400 | 2^-874 |

On each runtime, all 12 nonzero-charge cases return the wrong zero result;
all six exact-zero-charge controls pass. The independent reference evaluates
q*dt*E/m as an exact rational and converts its representable dyadic through
8,192-bit arithmetic. Every actual output bit and complete command/log hash is
retained in dated takeover evidence `particle-ratio-underflow/`, specifically
`probe-v1.jl`, `observations-112.toml`, `observations-11333.toml` and
`root-review-v1.json`. A zero command exit means observations were collected;
it does not turn these physical failures into passed tests.

Source owner: `src/particles/boris.jl::_validated_qm` rounds nonzero q/m to
signed zero; the later product guard sees an apparently exact zero factor.
That entire function is byte-identical in current main and frozen C1.

A related public EMPIC witness uses q=1 and minimum-subnormal positive mass.
On actual Julia 1.13.3, both EMPIC owners at fresh/established clocks throw after
NaN mutation in frozen `fd0a168`. C1 restores complete physical state in all
eight cases, but still rejects a finite representable complete kick. This is
additional rollback evidence and an explicit coefficient-range limitation;
minimum-runtime and repaired-result acceptance have not been established.
See dated takeover `empic-mass-ratio-range/FINDING_V1.md`.

## Canonical implementation stages

| Stage | Owners and contract | Required acceptance |
| --- | --- | --- |
| P1 CPU uniform/gathered | Preserve validated q and positive m through every IEEE species path, including normal quotients. Extend the existing kick/repair owner; keep ordinary native operation order and direct explicit-q/m semantics. Define and verify the original-factor result-range/publication policy before writes; retain direct-scalar semantics separately. | All original witnesses and controls, both signs, normal/subnormal/overflow ratio boundaries, empty/heterogeneous particles, aliases, nonzero B, genuine output overflow/error state and both actual runtimes. Ordinary bit/zero-allocation replay and a separately declared cold full-call resource budget. |
| P2 shared ES/hybrid owners | Inventory `_validated_qm` consumers and direct species divisions in `electrostatic.jl` and `extrapolated_leapfrog.jl`. Preserve loaded phase, predictor and variable-step semantics through the same factor owner. | Independent complete-stage references, public init/prime/step/restart, neighboring conservation and rollback, exact-source/compiler/resource comparisons. |
| P3 EMPIC priming/push | Preserve original ion q/m factors in both public owners and their canonical priming/Boris arithmetic. Keep C1 admission and whole-step undo; preserve all ordinary valid historical bits. | Both owners/dimensions, fresh/established/variable/subcycled clocks, electric/magnetic/zero-source controls, unsupported stored outputs, checkpoint/archive continuation and unchanged physical state on error. |
| P4 device/backend | Audit the Metal host preflight and gathered kernels plus every other affected caller. Use the documented native domain and reject unsupported original factors before launch when a backend lacks a cold owner. | Actual device execution or explicit preflight rejection on both supported runtimes; no silent CPU fallback, dropped charge or unsupported resident-solver claim. |

The complete caller map and cold operation/storage budget must be frozen before
production edits. Ray-dispersion charge/mass formulas have a separate model
contract and are not automatically covered by particle fixes.

Normal-ratio native replay is a compatibility requirement, not an accuracy
certificate. Include the original q/m division perturbation in the unchanged
four-ulp physical oracle and the full local rounding budget before edits. Four
normal-ratio Float16 witnesses on both actual runtimes exceed the unchanged
four-ulp bound against both original-factor and explicit-qm exact references.
Native/scalar/uniform/gathered bits match; an independent rational 3×3 solve
confirms all 32 outputs in takeover `p1-normal-native-accuracy-root-review-v1.json`.
The external scalar 48-node experiment passes 13,074 assertions per runtime on
actual Julia 1.12.7 and 1.13.3. Independent root replay checks every center/radius
bit and exact target: 192 coefficient-meaning contexts and 8,208 finite node
enclosures per runtime, plus 594 outward-primitive checks. Its sufficient test
accepts 66 contexts after the unchanged range guard. This bounded experiment
does not establish a universal or production certificate, batch hit rate or cost.
The separate whole-batch experiment passes 144,590 assertions per runtime.
Root complete source/derivation/replay review checks 570 batch contexts,
66,664 magnitude and 66,656 error enclosures, and 222 sufficient certificates
per runtime. All original inputs, 95 independently invalid native contexts and
38 invalid range-only batch bypasses remain recorded. This is bounded evidence;
production and cost gates remain unaccepted. The quotient-range adapter does
not repair this gap. Preserve
explicit scalar result-range semantics separately from public original-factor
publication; define any changed error/staging policy before implementation.

The planned CPU carrier retains concrete q,m,qm once per public push, including
normal quotients. Native arithmetic uses the already computed qm in the one
existing expression; the species accuracy certificate and exact cold reference
retain original q,m. Direct scalar calls keep their stored explicit-qm meaning.
No per-particle ratio division is introduced. Every unresolved scalar or batch
case must retain the correct factor meaning through exact replay; a range-only
batch flag cannot authorize the native accuracy bypass. Independently valid
finite parent bits are retained under the unchanged physical tolerance and
nearest-overflow rule. Required measurements include complete public uniform
and gathered calls, Array/strided midpoint storage, zero warmed ordinary heap
allocation and zero added retained storage. The new isolated candidate
implements this flow in source, with source/formatter proofs passing. Its first
minimum-runtime bounded run retains 615 passes, 43 failed position/midpoint
assertions and one gathered-field container error. The log classifies 32
position and 11 midpoint failures; it lacks initial-position/route labels, so
it does not establish the cause of every failure. The test used undefined
position storage while expecting zero origins, and a vector where the gathered
API requires a three-tuple. An isolated V2 changes exactly these fixtures:
three permanent-test edits and two external-harness zero initializations.
Complete inverses preserve production code, assertions, references and
tolerances. All 489 V1/V2 content pins and 468 Julia files/modes were checked.
The V1 test was edited and restored from immutable formatter bytes; its
timestamp changed and is not claimed untouched. The four-call V2 bounded run
remains active. A separately isolated V3 restores the historical generic
three-argument square sum: both runtimes pass 45 checks, and root replays all
228 stored values per runtime plus 48 trace nodes. Default generic dispatch
returns to 13 binary / one ternary addition; all three warmed BigFloat samples
match the parent in each of 12 contexts. Only this bounded generic compatibility
is accepted. A new four-call V3 physical harness adds external phase records
around unchanged tests; its source inverse and frozen inputs are reviewed and
released. Physical, compiler, resource, backend and integration gates remain open.

A separate unchanged-converter hardware-mode
probe first fails on an unavailable public method, then terminates with SIGSEGV
without result data. Both failures are retained; no mode/converter cause is
established for that crash. The revised typed probe completes 2,976 conversions
per runtime and restores raw rounding/flush state after every call. Independent
root replay finds 76 converter mismatches per runtime under changed hardware
modes; every nearest/gradual row is correct. IEEE packing of the already-rounded
integer is a separate canonical repair in progress. A second typed probe
checks stored-input extraction: 1,680 rows per runtime retain the input bits
and restore mode state; independent integer-field replay finds 144 wrong
rationals per runtime, 48 per IEEE type, all with flush-to-zero enabled. All
non-flush rows are correct. The canonical candidate extends the existing BigFloat
`Base.decompose` extraction body to IEEE types and packs already-rounded output
integers without floating scaling. Both bounded stdlib children pass 87,507
assertions per runtime; root replays every 14,240 rows/runtime, including all
original inputs and output witnesses. All 693 frozen inputs remain exact.
Broader package/caller, generic/nonfinite/BigFloat, compiler, strict ordinary
resource and actual Metal gates remain open. Neither this probe nor the output probe
accepts the complete cold Boris pipeline under changed hardware modes.
This is separate from the normal-mode production campaign.

The separately frozen 96-row native-node observation corpus adds three
exact-overflow failures on both actual runtimes, one each for Float16/32/64.
Let M=floatmax(T), U=ulp(T,M), q=m=1, dt=2, v=(M,0,0),
E=(U/4,0,0), B=0. The exact Boris target is (M+U/2,0,0), whose first
component nearest-rounds to infinity. The actual raw/native output remains
(M,0,0), and the old scalar range guard returns true. This new three-row
observation exercised the raw/native owner and actual guard only; public
scalar/uniform/gathered calls were not run for these three inputs. The earlier
four Float16 witnesses retain their separately verified public-call coverage.
Complete source/48-node/compiler captures and the independent exact references
are retained under dated `original-factor-p1-20261028/normal-accuracy-prototype-v1/`:
`prototype-nodes-112.toml`, `prototype-nodes-11333.toml`,
`prototype-observations-terminal-review-v1.json` and
`prototype-observations-terminal-receipt-v1.json`. Independent root review is
`takeover/p1-prototype-observations-root-review-v1.json`. All six observation
commands exit zero; that records collection integrity, not repaired numerical
acceptance. Both normal-accuracy and exact-output range closure remain open.

The original hardware widening proof does not extend automatically: Float16
q/m can reach about 2^40 from finite q,m, so its downstream squared products can
exceed the old Float32 repair envelope. Derive the complete factor expression
range or use the existing exact owner; merely changing qm's type is insufficient.
Keep local numerical errors at the declared operation scale and record all
new storage, copy and arithmetic costs. No full-optimization claim follows.

Wrong fixes excluded by the gates include clamping or replacing q/m by zero,
using a wider ratio that still overflows later products, duplicating the Boris
formula, validating only the first particle, and repairing velocity while
losing the coordinate or midpoint contribution.

## Progress

- [x] Record min/latest original-factor underflow reproductions and zero controls.
- [x] Related actual latest-Julia EMPIC overflow and C1 whole-state rollback.
- [x] Freeze the original P1 adapter source/caller contract and cold budget;
  preserve its seven source/test snapshots and 475 pins.
- [x] Pass the bounded 964-assertion source/AST preflight.
- [x] Retain the bounded three-type raw/native exact-overflow observations
  and independent root review; no public-call coverage is inferred for them.
- [x] Pass and independently replay the bounded actual Julia scalar-radius
  prototype on minimum/latest Julia, preserving all 96 inputs and all 48 nodes.
- [x] Complete bounded whole-batch source/derivation/evidence review on both
  runtimes, preserving the original inputs and every rejected numerical context.
- [ ] Complete production source proofs and normal-accuracy verification,
  preserving all four confirmed witnesses, the three raw overflow rows and
  the original four-ulp/nearest-infinity contract.
- [ ] Implement and CRC-verify P1, then P2–P4 as separate coherent slices.
- [ ] Current-main integration, docs, separate pushes and exact-SHA CI.

The [audit tracker](DEEP_DEBUG_TRACKER.md) owns live acceptance and the
[completeness plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md) retains the wider
simulator goal. No P1 numerical or publication acceptance is claimed.
