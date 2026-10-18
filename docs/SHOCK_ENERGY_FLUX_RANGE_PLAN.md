# Shock boundary energy-flux range

Status: ENERGY-008 is confirmed on Julia 1.12.7 and 1.13.2. Repair is planned
as a separate change after the current polytropic-energy integration.

## Evidence and owning path

**VERIFIED:** both runtimes pass 34 witness assertions against the shock-flux,
shock-constructor and particle-storage sources byte-identical to `15662e1`.
`boundary_energy_flux` returns infinity or zero for finite complete bulk and
particle kinetic fluxes. Its electron enthalpy flux also returns infinity for
a finite complete product. Evidence, exact rational fixtures, source hashes,
commands and terminal results are in
`/tmp/hybrid-shock-energy-flux-range-20261018` (`source-proof.json`,
`witness.jl`, `results.json`, `witness112.log`, `witness1132.log`).

For example, Float64 density/weight `2^-1000` and x velocity `2^600` give
bulk flux `2^799`, while the current owner returns `Inf`. Velocity `2^-600`
and density/weight `2^1000` give `2^-801`, while the owner returns zero.
The particle witness includes the stored SBP boundary norm before rounding.
With gamma 1.5, pressure `floatmax(T)/2` and x velocity 0.25, complete enthalpy
is finite while the intermediate pressure product overflows.

The canonical owners are `boundary_energy_flux` and `_ion_energy_flux` in
`src/diagnostics/shock_surface.jl` (lines 42–103 at the checked revision).
The existing `test/test_shock_diag.jl` already checks constructed boundary
states, H-weighted kinetic-flux conservation and invalid particle velocities.
This is separate from volume-integrated `electron_internal_energy`.

A further 18 nodal witness assertions pass per runtime (`nodal-results.json`).
They establish lost subnormal contributions before accumulation, overflowing
opposite contributions whose complete wall flux is one, NaN in a zero-weight
neighbor, and finite flux lost before division by a large SBP norm. These cases
require retaining the complete nodal expression, not only repairing a scalar
per-particle value.

## Repair contract and plan

- Retain each finite stored factor through bulk kinetic and electron enthalpy
  flux evaluation when native intermediate arithmetic loses range. Reuse the
  shared exact conversion owner; preserve ordinary arithmetic where certified.
- For particle flux, include every deposited term, CIC weight and final SBP norm
  in the complete nodal expression. Correcting an individual particle's flux
  after materialization is insufficient when deposition, accumulation or
  normalization changes the representable range.
- Preserve generic/nonfinite behavior and the existing singular gamma-one
  enthalpy convention. Keep particle validation and exception order, inputs,
  deposited node locations, boundary ordering and the current mass convention.
  This repair changes arithmetic, not the shock model or its evolution.
- Audit scalar total-flux combination separately before extending the contract;
  do not silently claim a complete total when only its components are repaired.
- Keep ordinary allocation equal to the existing output-vector requirement.
  Measure complete cold-call cost, temporary storage and grid/particle scaling.

Acceptance must cover the original witnesses, signed and zero factors,
subnormal and overflow boundaries, exact/neighboring ties, zero CIC weights,
multiple particles whose individually rounded contributions lose the nodal
sum, views, empty populations, validation failures and unchanged inputs.
Reject patches that only check `isfinite`, only repair squared velocities,
or round per-particle flux before its stencil and H normalization.

Run independent rational nodal/flux oracles, all shock diagnostic and solver
regressions, ordinary bitwise/allocation comparisons, concrete and complete
compiler checks, Aqua, four-thread shared-input checks and formatting on both
supported runtimes. Commit and push the verified repair separately; require
green CI for that exact SHA. This item remains open until those gates pass.
