# Particle charge/mass factor range

P1 CPU kick code is implemented; final GitHub CI is pending. Its current
contract and tracker live in [CPU Boris accuracy](BORIS_ACCURACY_IMPLEMENTATION.md).
Nonrelativistic priming/predictor factor routing, sufficient term-scale accuracy
admission and exact endpoint replay are implemented; GitHub acceptance remains
open. Both relativistic EMPIC owners also retain original factors through priming,
kick repair and output publication; GitHub verification is pending. Metal's
native pushes now validate original factors and reject unsupported quotients
before launch; verification is pending. Device factor recovery and full normal
relativistic accuracy remain open. This plan does not broaden C1's domain.
CPU hybrid/CAM-CL and ElectrostaticPIC collection validation admits finite
original charge and positive mass through the canonical factor owner, including
derived quotient overflow and underflow. Sign, storage, neutrality and actual
arithmetic limits remain; essential full-model and restart checks await GitHub.

The physical expression uses the original stored q, positive m, dt and field
operands. A rounded derived ratio must not silently erase a representable result.
The standalone `boris_kick(..., qm, dt)` still means the supplied stored `qm`;
it cannot recover factors that the caller already divided.

## Retained failure evidence

VERIFIED on historical main `cd8b94d`, actual Julia 1.12.7 and 1.13.3:
public `push_uniform!` and `push_gathered!` returned zero velocity and coordinates
for nonzero-charge inputs below. Each used particles at rest, B=0, E=(0,0,dt),
both charge signs and independent gathered midpoint storage. All exact kick,
position and midpoint results are representable dyadics.

| Type | Charge magnitude | Mass | dt and electric amplitude | Exact kick magnitude |
| --- | --- | --- | --- | --- |
| Float16 | 2^-24 | 2^10 | 2^10 | 2^-14 |
| Float32 | 2^-149 | 2^80 | 2^60 | 2^-109 |
| Float64 | 2^-1074 | 2^600 | 2^400 | 2^-874 |

The historical `_validated_qm` rounded the quotient to signed zero; a later
product guard treated it as an exact zero factor. Original observations,
independent rational references, failed prototypes and their raw outputs remain
in the dated audit evidence. A successful observation command did not turn its
wrong outputs into passes.

Separate retained normal-ratio Float16 witnesses exceed the four-ulp contract,
and raw Float16/32/64 boundary witnesses retain finite floatmax where exact
nearest rounding selects infinity. These motivate the accuracy guard as well
as factor preservation. Bounded experiments and physical replay passed before
integration; they do not establish whole-project or performance acceptance.

The shared conversion repair is integrated in `7300c4f`. Its bit packing and
stored-input extraction address the earlier converter-mode failures. That does
not by itself verify every cold caller under all hardware contexts. The earlier
probe SIGSEGV remains unexplained and is not attributed to a converter instruction.

## Implementation stages

| Stage | Canonical work | Required verification |
| --- | --- | --- |
| P1 CPU kicks | Implemented: immutable q,m,qm carrier, one native expression, scalar/batch accuracy admission, exact replay and CPU model kick threading | Essential regressions and complete package/compiler, formatter and MPI checks in GitHub; performance and backend guarantees remain separately unverified |
| P2 CPU priming and prediction | Implemented factor routing, range/accuracy admission and exact endpoint replay through hybrid, electrostatic and EMPIC priming and the hybrid endpoint predictor | Essential factor, overflow endpoint, normal accuracy, phase and checkpoint regressions await GitHub CI; complete public resource acceptance remains open |
| P3 relativistic factors | Implemented original-factor range replay through priming and kicks in both EMPIC owners, retaining C1 admission and whole-step undo | Essential axial-reference, variable/subcycled phase, restart and rollback checks await GitHub CI; full normal-range relativistic accuracy remains open |
| P4a native Metal admission | Implemented shared host validation of original q,m and rejection of charged nonfinite/zero/subnormal quotients before either push | Essential portable admission and conditional device no-mutation controls await GitHub/suitable remote hardware; no device accuracy or cost acceptance |
| P4b device recovery and priming | Preserve original factors through a device numerical owner, with explicit supported arithmetic and state ownership | Scaled/exact device arithmetic design, suitable remote hardware and backend equivalence; no silent CPU fallback |

Before each remaining stage, inspect its canonical callers, staggered time-level
contract, error/publication boundary and checkpoint ownership. Extend shared
numerical owners. Keep native operation order for admitted inputs and preserve
valid native bits under the declared normwise tolerance. Range-only batch
admission cannot substitute for an accuracy check.

Record real allocation or timing regressions against existing gates; do not
relax a limit to obtain a pass. Add focused physical and boundary regressions,
then run them and the existing suites in GitHub Actions. Large local protocol,
provenance and replay campaigns are not prerequisites for code review.

Relativistic native range admission requires a normal derived quotient (or an
original zero charge) in the existing nearest/gradual context. Kicks additionally
use the bounded native/DD residual and joint batch accuracy owners in the
[normal-accuracy plan](RELATIVISTIC_NORMAL_ACCURACY_PLAN.md). A certified DD
attempt precedes the inherited exact-gap 256-bit replay; both retain original
q,m rather than a lossy derived quotient. This does not prove arbitrary-extreme
accuracy of that replay.
Priming reuses the scalar DD certificate only when that quotient is exact and
its EFT context is admitted; the decision is cached on the first cold particle
per species loop. Cold allocations and non-exact-quotient output-bit changes are
not covered by a new resource guarantee. Explicit-qm scalar semantics remain.
This range slice does not prove full normal relativistic accuracy, every hardware
mode pipeline or trajectory-wide correctness.

Metal uniform and gathered pushes remain native Float32 operations. Both validate
finite original q and positive finite m, allowing true signed-zero charge or a
finite normal quotient. A charged quotient rounded to zero, subnormal or infinity
raises `ArgumentError` before launch/mutation, including empty populations and
dt=0. Original q,m may be subnormal when their host ratio is normal. No device
kernel, ordinary arithmetic or backend dispatch changes were needed for that
host-factor slice. Uniform coefficient range admission remains. The separate
[gathered preflight candidate](METAL_GATHERED_VALIDATION_PLAN.md) now adds a
finite-field scan, guarded publication scheduling and explicit native field
storage admission; remote correctness, compiler and unchanged cost gates remain
open. These checks do not recover unsupported factors. Device priming and
GPU-resident coupled stepping remain unimplemented; CUDA still supplies storage/telemetry only.

## Stored scalar admission and checkpoint compatibility

Original IEEE factor admission now reads stored sign/zero for Float16/32/64
after the existing declared-type conversion and finite checks. The shared
`_validated_charge_mass` owner accepts finite signed charge and positive finite
mass; `ParticleSet` and checkpoint scalar checks delegate to it. Hybrid
collections retain positive ion charge. Hybrid/ES/EM collection weights retain
finite nonnegative admission, including signed zero, with negative subnormals
classified by storage. Generic number types retain their comparison semantics;
the global positive/nonnegative parameter validators are unchanged.

Scalar Hybrid/CAM-CL save and restore now enforce that existing stepping factor
domain. Previously type-correct nonfinite charge or nonpositive/nonfinite mass
could pass this scalar checkpoint gate. This is an intentional compatibility
tightening, with no new schema: valid signed charge and finite extreme ratios
still restore verbatim. Collection positive-ion policy, normalized ES/EM,
transport q=1, layout/phase rules and other legacy snapshot field/weight value
policies are unchanged. Rejection precedes file or destination publication.

Native zero charge returns the original signed zero before division, including
with a positive subnormal mass. Nonzero native quotients still require finite
normal storage; this adds no device recovery. Fixed-row tests record requested
and active FTZ mode and original words, restore the caller's mode in `finally`,
and cover public CPU pushes, collection rejection and checkpoint preservation.
These controls are unrun on the candidate; all Julia checks remain GitHub-only.
Source-level exposure to permitted FTZ comparison behavior is distinct from a
measured hardware failure or support for every model pipeline under FTZ.

## Progress

- [x] Preserve the original-factor, normal-accuracy and overflow witnesses.
- [x] Implement P1 canonical arithmetic, admission and exact replay.
- [x] Thread original factors through CPU public and nonrelativistic model kicks.
- [x] Compose the shared IEEE conversion repair and focused regressions.
- [ ] Complete exact-revision GitHub package/compiler, formatter and MPI checks.
- [x] Implement P2 original-factor nonrelativistic priming/predictor range replay.
- [x] Implement sufficient Euler term-scale admission and exact overflow classification.
- [ ] Verify P2 accuracy, endpoints and phase behavior in GitHub.
- [x] Align CPU hybrid/CAM-CL and electrostatic collection factor admission.
- [x] Share stored IEEE sign/zero factor admission and scalar Hybrid/CAM-CL checkpoint validation.
- [ ] Verify actual-mode factor/sign/dt and checkpoint preservation essentials in GitHub.
- [ ] Verify full-model factor range and checkpoint continuation in GitHub.
- [x] Implement P3 relativistic original-factor range routing and publication.
- [ ] Verify P3 phase/restart/rollback in GitHub and audit normal accuracy.
- [x] Implement P4a shared native Metal original-factor admission.
- [ ] Verify P4a portable and remote-device rejection/admission checks.
- [ ] Complete P4b device factor recovery, priming and hardware checks.
- [ ] Establish full public resource and comparative performance results.

The [audit tracker](DEEP_DEBUG_TRACKER.md) owns current acceptance; the
[completeness plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md) retains the wider
simulator goal. No full original-factor trajectory or global optimization claim
is established.

The [normal relativistic accuracy work order](RELATIVISTIC_NORMAL_ACCURACY_PLAN.md)
records the parent's six production failures and the bounded correction's
minimum/latest observation with zero production violations. Focused particle
essentials also pass on both runtimes; full exact-revision CI, retained legacy
gates and resource acceptance remain pending.
The inherited extreme 256-bit replay does not gain a general accuracy proof.
