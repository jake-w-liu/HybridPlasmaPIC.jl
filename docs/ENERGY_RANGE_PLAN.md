# Energy diagnostic range repairs

Status: field-energy repair verified on Julia 1.12.7 and 1.13.1. Relativistic
particle energy remains separate planned work; ordinary diagnostic performance
remains an optimization target.

## Verified behavior

`em_field_energy` returns NaN for finite light speeds and B=1/c although the
complete energy is finite. Both EMPIC interfaces are affected. Squaring fields
or materializing cell volume also loses finite integrated energies in 1D–3D.
An exact-rational oracle verifies 40 witness assertions per Julia version:
`/tmp/hybrid-field-energy-baseline112.log` and `113.log`.
The first probe included four non-failing controls as expected failures; their
36-pass/4-fail log is retained as `...baseline-v1-112.log`. The corrected probe
uses field amplitudes whose squares actually exceed the hardware range.

`kinetic_energy_relativistic` rejects valid v=c/2 when c² underflows or overflows.
It also returns NaN for an empty population at large c, NaN instead of 0.5 for
v=1 at large c, and zero when a tiny weight has a finite mass-weighted energy.
An independent 4096-bit stable expression verifies 18 witness assertions per
version: `/tmp/hybrid-relativistic-energy-baseline112.log` and `113.log`.
These runs use frozen `4f4fec4` production bytes for the affected functions.

## Field-energy pre-edit brief

- Mode: repair numerical reduction at the public diagnostic owner.
- Requirement: finite represented fields, c and grid spacings must retain a
  representable complete integral, including zero fields at extreme c. Preserve
  ordinary arithmetic order, return type, state and nonfinite-input visibility.
- Baseline: the exact-rational witness above isolates c², field squares, products,
  positive accumulation and cell-volume materialization from stepping/FFT errors.
- Cause: both overloads in `src/models/full_pic.jl` form these intermediates
  before the final ½-volume factor can restore range.
- Contract: this is the energy of the stored fields, using the existing EM
  normalization. A genuinely unrepresentable final hardware result remains
  infinite. No field mutation, extra state or alternate public diagnostic.
- Owner and surfaces: both `em_field_energy` overloads, one shared private exact
  complete-integral helper, existing exact conversion utilities, and independent
  diagnostic tests. Current correction and stepping are outside this slice.
- Compatibility: keep the specialized 1D and general component/summation order
  on ordinary inputs. Preserve invalid/nonfinite field behavior. Other energy,
  momentum and work diagnostics require separate assessment.
- Strategy: guard all range-sensitive products and accumulation, including
  spacings before their product. Eligible exceptional finite inputs use a
  streaming rational sum and the existing final-rounding owner. No retained
  full-field copy, transcendental operation, retry or hidden clamp.
- Verification: independent rational rounding cells; both interfaces, precisions,
  dimensions and components; large/small c, fields and volumes; zero fields,
  true output overflow, signed/subnormal values and nonfinite controls. Compare
  ordinary results and allocations against frozen baseline, check concrete and
  package compiler reports, run EM/diagnostic neighbors on both Julia versions,
  then inspect the exact diff and formatter output.

Changing c²B² to (cB)² alone fails field/volume witnesses. Scaling by one global
field maximum can erase small contributions later recovered by volume. Computing
an exact sum only after a floating volume product preserves the volume defect.
The regression set must distinguish all three incomplete repairs.

## Final field-energy verification

The final v3 retains the original positive summation order and tracks the minimum
square of nonzero fields. The minimum detects square/product underflow; overflow
remains visible in the positive total. The original guard and minimum selector
agree in 800 seeded comparisons on each Julia version (the 1.12 prototype differs
from final code only in local names and explanatory comments).

The final immutable source passes **27,896 ES/EM assertions on each version**,
including 462 field-energy assertions. Independent rational rounding cells cover
a tiny square beyond an exact midpoint, positive-sum overflow before a small
volume factor, all spacing associations, true final overflow, nonfinite controls
and bitwise ordinary results. Logs: `/tmp/hybrid-field-energy-v3-full112.log` and
`113.log`; focused results are in `...v3-focused112.log` and `113.log`.

All 24 concrete compiler signatures per version are clean (40 runtime/analysis
assertions each). Package report messages match the baseline exactly: 26 on 1.12
and 25 on 1.13, with no suppressed report or raised ceiling. Logs:
`...v3-jet112.log`, `...v3-jet113.log`. JuliaFormatter 1.0.62 is idempotent at the
repository's configured margin 100 (`...v3-format-configured.log`). An initial
read-only check used an incorrect margin 92 and failed; it changed no source.

Final source SHA256, recorded with the immutable snapshot in
`/tmp/hybrid-field-energy-v3-manifest.json`:
`3dfda8965b222a34254c23d609dafc779e2841f87137ce6fa716542ea5c18fb1`.
The snapshot excludes the concurrent continuity, publication and radiation work.

## Resource evidence and open optimization

Both versions pass 152 ordinary comparison/allocation assertions and 48 cold
resource assertions in `/tmp/hybrid-field-energy-v3-resources-warmed112.log` and
`113.log`. Eight ordinary configurations retain exact baseline values and 0 B
allocation. The 24 cold configurations measure 24,120–3,543,640 B per call across
Float32/64, 1D–3D, two sizes and extreme light speeds. Nine balanced timing pairs
per configuration give candidate/baseline median ratios 1.13–2.40 on 1.12 and
1.53–3.20 on 1.13. These shared-machine measurements establish an optimization
target, not an isolated performance ceiling.

The first guard candidate measured 1.88–3.78 times baseline latency; the final
minimum-square selector removes its repeated per-product guard work. The local
timing series are noisy, so they do not establish a controlled speedup between
candidates. Preserve the complete-integral contract during further optimization.

Earlier resource runs fail exact equality between repeated cold allocation
measurements: one minimum-field prototype case, two minimum-square cases, and
two final-v3 cases. Numerical and ordinary allocation checks pass in those runs.
The final-v3 raw samples include 1,065,408 versus 1,885,824 B per three calls.
The failed logs remain `...minima-resources112.log`, `...squares-resources112.log`
and `...v3-resources112.log`; successful uninstrumented replays do not explain them.

VERIFIED: explicit first-time FFT-plan finalization allocates 410,208–411,104 B
per plan and records compiler time; repeated finalization of the same plan type
allocates zero. The standalone probe passes 30 assertions:
`/tmp/hybrid-field-energy-finalizer-probe112.log`. Compiler traces from the resource
run also show those finalizers compiling. This proves that unrelated plan cleanup
can enter process-wide allocation measurements, but does not individually trace
every allocation in the earlier failed samples.

The final resource harness warms destruction of disposable plans for both
precisions and all dimensions before measurement. It keeps the original numerical,
zero-allocation and exact repeated-cold-allocation assertions unchanged. No
production code, assertion, GC policy or arithmetic was changed to pass that gate.
The first tracing attempt used `@timed` and measured 16 B extra on ordinary calls;
its failed log is retained as `...v3-resources-traced112.log`. The corrected trace
uses the original `@allocated` expression (`...v3-resources-traced-v2-112.log`).

## Relativistic-energy follow-up

Trace the diagnostic's independent speed check, stable gamma-minus-one formula,
weight reduction and mass/light-speed factors. Define the accepted boundary
against exact stored components and the shared pusher's narrower supported
range before editing. Preserve the Newtonian limit, valid empty populations,
error visibility, writable views and ordinary cost. This plan does not claim
that the radiation candidate or field-energy repair fixes particle energy.
