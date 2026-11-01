# Particle charge/mass factor range

P1 CPU kick code is implemented; final GitHub CI is pending. Its current
contract and tracker live in [CPU Boris accuracy](BORIS_ACCURACY_IMPLEMENTATION.md).
Nonrelativistic priming/predictor factor routing, sufficient term-scale accuracy
admission and exact endpoint replay are implemented; GitHub acceptance remains
open. Both relativistic EMPIC owners also retain original factors through priming,
kick repair and output publication; GitHub verification is pending. Device
support and full normal-range relativistic accuracy remain open. This plan does
not broaden C1's supported relativistic domain. CPU hybrid/CAM-CL and
ElectrostaticPIC collection validation now admits finite original charge and
positive mass through the canonical factor owner, including derived quotient
overflow and underflow. Sign, storage, neutrality and actual arithmetic limits
remain enforced; the essential full-model and restart cases await GitHub.

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
| P4 devices and parallel owners | Preserve original factors or reject unsupported factors explicitly before launch through the backend's existing owner | Suitable remote hardware, backend equivalence, launch/error behavior and documented native domain; no silent CPU fallback |

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
original zero charge) in the existing nearest/gradual context. Cold evaluation
uses the original rational q/m in the existing 256-bit BigFloat expression.
Priming reuses the scalar DD certificate only when that quotient is exact and
its EFT context is admitted; the decision is cached on the first cold particle
per species loop. Cold allocations and non-exact-quotient output-bit changes are
not covered by a new resource guarantee. Explicit-qm scalar semantics remain.
This range slice does not prove full normal relativistic accuracy, every hardware
mode pipeline or trajectory-wide correctness.

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
- [ ] Verify full-model factor range and checkpoint continuation in GitHub.
- [x] Implement P3 relativistic original-factor range routing and publication.
- [ ] Verify P3 phase/restart/rollback in GitHub and audit normal accuracy.
- [ ] Complete P4 device factors and their hardware checks.
- [ ] Establish full public resource and comparative performance results.

The [audit tracker](DEEP_DEBUG_TRACKER.md) owns current acceptance; the
[completeness plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md) retains the wider
simulator goal. No full original-factor trajectory or global optimization claim
is established.
