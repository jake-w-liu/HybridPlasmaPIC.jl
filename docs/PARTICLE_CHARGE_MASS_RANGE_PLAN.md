# Particle charge/mass factor range

Status: confirmed current-main defect. The isolated P1 CPU adapter is written
and passes a 964-assertion source/AST preflight; numerical, compiler, resource
and integration acceptance remain open. This work is separate from C1's supported relativistic output closure.
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
| P1 CPU uniform/gathered | Preserve validated q and positive m through the exceptional branch. Extend the existing kick/repair owner; keep ordinary native operation order and direct explicit-q/m semantics. Define and verify the original-factor result-range/publication policy before writes; retain direct-scalar semantics separately. | All original witnesses and controls, both signs, normal/subnormal/overflow ratio boundaries, empty/heterogeneous particles, aliases, nonzero B, genuine output overflow/error state and both actual runtimes. Ordinary bit/zero-allocation replay and a separately declared cold full-call resource budget. |
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
A 48-node outward-error certificate is proposed; its proof/prototype and cost
gates are unaccepted. The quotient-range adapter does not repair this gap. Preserve
explicit scalar result-range semantics separately from public original-factor
publication; define any changed error/staging policy before implementation.

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

- [x] Fixed min/latest original-factor underflow reproductions and zero controls.
- [x] Related actual latest-Julia EMPIC overflow and C1 whole-state rollback.
- [x] Freeze the original P1 adapter source/caller contract and cold budget;
  preserve its seven source/test snapshots and 475 pins.
- [x] Pass the bounded 964-assertion source/AST preflight.
- [ ] Prove and prototype the canonical normal-accuracy certificate; retain
  all four confirmed witnesses and the original four-ulp contract.
- [ ] Implement and CRC-verify P1, then P2–P4 as separate coherent slices.
- [ ] Current-main integration, docs, separate pushes and exact-SHA CI.

The [audit tracker](DEEP_DEBUG_TRACKER.md) owns live acceptance and the
[completeness plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md) retains the wider
simulator goal. No P1 numerical or publication acceptance is claimed.
