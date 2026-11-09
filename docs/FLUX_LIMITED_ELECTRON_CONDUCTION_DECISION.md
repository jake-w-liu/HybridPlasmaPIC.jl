# Flux-limited electron conduction decision record

Status: **design decisions settled for implementation; no production source has
been added by this record.** This closes the pre-edit decision step requested by
P40 and the frozen read-only proposal. It does not close any execution,
scientific or resource gate. The original proposal remains
`handoff/2026-11-09/flux-proposal/PRE_EDIT_CONTRACT.md`; this record is the
operator disposition of its explicitly unapproved choices.

## Admitted model and API

`FluxLimitedConductingElectrons(gamma=5/3; kappa, electron_mass,
alpha=0.4*sqrt(2/pi), conductivity=:constant, reference_temperature=1)` is a new
stored-pressure closure:

- `gamma > 1`, `kappa >= 0`, `electron_mass > 0`, `alpha > 0`, and
  `reference_temperature > 0` are validated in the original precision and after
  grid conversion. A nonzero original coefficient that underflows to zero at
  grid precision is rejected.
- `conductivity=:constant` uses exact stored `kappa`.
- `conductivity=:spitzer` uses `kappa(Theta) = kappa_ref*(Theta/Theta_ref)^(5/2)`.
  The normalized reference value is supplied by the caller; no SI numerical
  coefficient is silently assigned.
- The host is the serial periodic CPU `HybridStepper`, D=1–3, Float32/Float64,
  NGP/CIC/TSC, scalar unit-charge positive-ion collections, `eta>=0`,
  `etaH>=0`, `de2=0`, with the existing pressure/RK4/rollback owner. Other
  backends, arbitrary coefficient arrays, anisotropic limiting, wall fluxes,
  implicit evolution, drag/inertia and evolving tensors are separate slices.

The flux is `q = -kappa(Theta) grad(Theta)/(1 +
kappa(Theta)|grad(Theta)|/q_sat)` with `Theta=pe/ne`,
`q_sat = alpha*pe*sqrt(Theta/mu_e)`, `mu_e=me/mi`. Active nonconstant `Theta`
requires strictly positive `pe` and `ne` at every actual RK trial. A globally
constant `Theta`, including the all-zero pressure extension, returns exact zero
heat before any reciprocal-temperature calculation. **Mixed zero/positive active
cells are rejected transactionally; clipping or an invented cold-front rate is
not admitted.**

## Discretization decisions

- Use the canonical stored odd derivative `D_d` and the positive Nyquist
  multiplier `N_d`; form the enlarged gradient `G=(D_1 Theta,...,N_1 Theta,...)`.
- One shared mobility `a = sqrt(K)/(1+sqrt(K*g2/C))` limits the entire vector;
  per-axis independent limiting is not admitted because it does not prove a
  vector cap.
- The physical flux uses the first D components. The Nyquist components are an
  explicitly declared diagonal completion, not extra spatial directions.
- Prepare the original quotient and the generally nondyadic flux with the
  existing two-term owner; project the computed source to exact zero mean on
  each axis line before publication.
- Certify the **actual prepared flux inputs** (with removed anchors restored)
  against `sum(q_tilde_c^2)<=C` at every cell. A failure rejects; there is no
  clipping, silent re-preparation or epsilon relaxation.
- Root evaluation uses the existing scaled positive-rational square-root owner
  with directed enclosures; accept only `0<aL<=aH` and `aH-aL<=16u*aL`. The
  lower enclosure must certify the cap. `kappa=0` must reduce to the existing
  TransportElectrons arithmetic with no quotient/transform/root/entropy work.
- `alpha` is stored as the actual converted `T` value of
  `0.4*sqrt(2/pi)`; the evaluation/conversion order is documented and saved.

## Energy, entropy and persistence

- The complete signed pressure RHS is the only rounded publication owner.
  Conductive input is `R_cond=(gamma-1) H`, with `H=-div(q)` in the declared
  discrete map.
- On the new active family, compute the projected-source entropy rate
  `sigma=V*(gamma-1)*sum(H/Theta)` from actual original positive quotients and
  the full prepared/projected source. Reject a strictly negative `sigma`.
  This is a projected conductive-source gate, not an exact RK/full-step entropy
  theorem.
- Persist an inner `conduction_policy_version=4` with exact tags for the law,
  original factors, two preparations and projected admission. Keep outer
  schemas 1–3 and CGL 4 literal. Persist the closed scalar factors, `pe`, carried
  force, fields, clocks and particle phase; never persist quotient/mobility/FFT
  scratch.
- Restore validates the incoming closure/policy/pressure/ownership graph before
  detaching or publishing. Stale destination values never select incoming
  physics. Same-concrete-type zero-to-active restore is supported because the
  derived scratch always exists.
- Timeline rollback reuses the existing HybridStepper transaction for pressure,
  fields, particles, carried force, clock and RNG. Collisional composition leaves
  pressure fixed during ion–ion scattering. Optional thermal exchange requires
  the new closure electron mass to equal the exchange-bound mass and its
  existing grid gamma policy.

## Implementation scope

Planned owners: `src/electrons/flux_limited_conduction.jl` (new),
`src/electrons/closures.jl`, `src/electrons/energy_equation.jl`,
`src/integrators/extrapolated_leapfrog.jl`, `src/diagnostics/conservation.jl`,
`src/io/checkpoint.jl`, `src/integrators/collisional_hybrid.jl`,
`src/electrons/thermal_exchange.jl`, `src/utils/exact_arithmetic.jl` (only the
existing normalized-root extraction/forwarder if needed), plus one focused test
and the feature guide. Old policies 1/2/3 retain literal arithmetic and
persistence meanings.

## Gates still open

Independent direct-DFT source/cap/Nyquist/entropy oracles, constant/zero/classical
limits, Spitzer/coefficient range rescue, positivity and late rollback, raw/archive
continuation, exact-head minimum/latest CI/JET/format, and complete-step
allocation/peak-memory/FFT/time resource evidence. No implementation is accepted
by this decision record.
