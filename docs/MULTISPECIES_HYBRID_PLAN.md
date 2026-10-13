# Multispecies hybrid evolution

Status: design; implementation and acceptance checks pending.

## Verified gap

`compute_moments_multi!` in `src/models/hybrid_pic.jl` already deposits charge
density and charge-weighted velocity from several ion populations. The public
`HybridStepper` and `CAMCLStepper` initialization and stepping methods accept
one `ParticleSet`. A capability check at `c76e39c` on Julia 1.13.1 confirmed that
the collection calls have no applicable method; log:
`/tmp/hybrid-multispecies-capability-baseline.log`.

## Contract and ownership

- Add collection overloads to `init!`, `step!`, `init_camcl!`, and `step_camcl!`
  for `AbstractVector{<:ParticleSet{D,T}}`. Each population retains its own
  charge, mass, weights, identity, and storage. All use the stepper's shared
  grid, shape, closure, electromagnetic fields, and timestep history.
- Use the existing quasineutral moment definition:
  `n_e = sum(q_s*n_s)` and `u_i = sum(q_s*n_s*u_s)/max(n_e,nfloor)`.
  Apply the density floor after summing populations. Push each population with
  its own `q_s/m_s`.
- The collection interface represents positive ions; validate positive finite
  charge and mass and a finite charge-to-mass ratio before mutation. Validate
  every population's array lengths, axes, finite phase space, nonnegative
  weights, and independent storage before advancing any population. Keep
  existing single-population behavior compatible.
- Retain the documented staggered velocity and integer position/field levels.
  Midpoint moments must include every population before magnetic advancement.
  The final velocity predictor must use the same carried fields for every
  population before their combined moments produce the next carried field.
- Empty populations and a typed empty collection deposit zero moments and
  retain the existing density-floor treatment. A zero timestep remains a
  no-op. Populations must share a time level; call initialization again when
  introducing newly loaded physical velocities into an existing run.
- Preserve CAM-CL's persistent magnetic copies and existing closure restriction.
  Do not label the existing predictor algorithm an exact implementation of
  every step of the Matthews paper.

## Implementation plan

1. Extract the shared raw charge-moment accumulation and normalization from
   `compute_moments_multi!`, keeping its current public scratch arguments.
   Reuse those operations in collection stepping.
2. Extract the existing particle gather/prime/kick/drift/midpoint operation so
   scalar and collection paths use the same kernel. Reuse the stepper's
   particle buffers at the total particle count; contiguous views select each
   population's workspace. Do not concatenate or copy population state.
3. Reuse existing grid buffers according to their lifetime: midpoint totals
   reside in `fn/fui`; per-population deposits use public moment arrays. After
   the magnetic subcycle, `fn/fui` are available as deposition scratch for
   integer and predictor moments. Keep final field integration canonical.
4. Extend checkpoint/archive capture, validation and restore to preserve the
   ordered collection, its per-population parameters and all persistent
   integrator state. Validate every destination before changing any one of
   them. Reuse single-population serialization helpers.
5. Add an executable two-ion example, public API documentation, independent
   physics tests, and a reproducible allocation/runtime benchmark. Update the
   main audit tracker and push each verified implementation slice to `main`.

The MPI collection interface is a separate follow-on: its deposition must
reduce the total charge moments across ranks before normalization, and migrate
each population without losing species identity. GPU and nonperiodic paths
remain separately tracked capabilities.

## Acceptance checks

- One proton population supplied as a collection agrees with the scalar path
  for both integrators, all shapes, Float32/64 and 1D–3D. Compare the old scalar
  implementation against the extracted kernel for fixed and changing steps.
- Uniform two-ion relaxation in a constant magnetic field has an independent
  analytic oracle. With `a_s=q_s*n_s/n_e`, the relative velocity obeys
  `d(v_1-v_2)/dt = (a_2*q_1/m_1+a_1*q_2/m_2)*(v_1-v_2) × B`;
  mass-weighted momentum is constant. Check the staggered velocity time,
  integrated positions, timestep refinement and unequal charge/mass/density.
- Build the linear transverse multispecies wave matrix directly from the
  Lorentz force, Ampere, Ohm and Faraday equations. Compare small-amplitude
  evolution and temporal convergence with its matrix exponential, including
  field subcycles and changing steps.
- Check stationary pressure-balanced/drifting cases and species-order
  equivalence within roundoff. Exercise empty populations, different counts,
  writable views and population-count changes.
- Reject malformed later populations, invalid charge/mass, aliasing and clock
  overflow before changing earlier populations, fields or clock. Check zero
  timestep and existing CAM-CL closure restrictions.
- Compare uninterrupted and restored multi-ion runs exactly, including changed
  timesteps/subcycles, externally changed magnetic fields, malformed archives,
  independent views and incompatible destinations.
- Run focused regression suites, JET/Aqua/format checks, and full package suites
  on Julia 1.12.7 and latest stable Julia. Measure warmed allocations and cost
  for realistic unequal population sizes; do not claim a speedup without data.

Wrong implementations these checks must reject include using number-weighted
velocity, flooring each population before summation, reusing the first species'
charge-to-mass ratio, advancing fields separately for each population, and
recomputing the predictor field partway through the population loop.
