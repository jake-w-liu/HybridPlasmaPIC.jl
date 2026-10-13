# Energy diagnostic range repairs

Status: two public defects reproduced on Julia 1.12.7 and 1.13.1. Field-energy
repair is the next implementation slice; relativistic particle energy remains
separate planned work.

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

## Relativistic-energy follow-up

Trace the diagnostic's independent speed check, stable gamma-minus-one formula,
weight reduction and mass/light-speed factors. Define the accepted boundary
against exact stored components and the shared pusher's narrower supported
range before editing. Preserve the Newtonian limit, valid empty populations,
error visibility, writable views and ordinary cost. This plan does not claim
that the radiation candidate or field-energy repair fixes particle energy.
