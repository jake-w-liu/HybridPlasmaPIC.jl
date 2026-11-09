# Field-aligned scalar electron conduction

`FieldAlignedConductingElectrons(γ=5/3; κ_parallel, κ_perpendicular=0)` adds
constant field-aligned heat transport to the stored scalar pressure of a
serial periodic CPU `HybridStepper`. It uses the existing coupled magnetic and
pressure RK4 stages, carried ideal ion force and Ohmic/hyperresistive heating.
The admitted grid precision is Float32/64, D=1–3, with NGP/CIC/TSC, a single
unit-charge ion population or an ordered positive-ion collection. Fixed-step
ion-ion composition uses `CollisionalHybridStepper`.

Seed finite nonnegative `st.fields.pe` and finite `st.fields.B` before `init!`.
Both conductivities must be finite, nonnegative and survive conversion to grid
precision independently; checking only their maximum would lose a smaller
positive coefficient. Gamma must remain greater than one at grid precision.
The model admits constant `η, ηH >= 0` and requires `de2=0`, as for scalar
[transport](EVOLVING_ELECTRONS.md#constant-hyperresistive-transport).

## Tensor and represented channels

At each actual joint RK stage, let `Θ=Q(pe)/Q(n)` use the original stored pressure
and frozen midpoint density. Q denotes exact rational conversion. The effective
grid coefficients define

```text
A = κ_perpendicular I
    + (κ_parallel − κ_perpendicular) B Bᵀ / (B·B).
```

All three original B components contribute to `B·B`, including in 1D/2D. Unequal
effective coefficients require nonzero B at every participating cell, stage and
completed endpoint. No direction floor or isotropic fallback is inserted.
Exact rational projector arithmetic avoids an overflowing or lost native B²
sum. Equal coefficients permit zero B and delegate to the existing isotropic
operator. Two zero coefficients retain the transport reduction without a
pressure/density quotient evaluation.

The real odd Fourier derivative D has the canonical zero-Nyquist policy. N_d
is the positive real Nyquist projector times its wavenumber on each even axis;
it vanishes on odd axes. The declared heat operator is

```text
L_A Θ = Σ_d D_d (Σ_j A_dj D_j Θ) − Σ_d N_d (A_dd N_d Θ).
Fp_heat = s L_A Θ,      s = Q(RN_T(γ_grid − 1)).
```

Only diagonal Nyquist complements are present. Off-diagonal Nyquist products
would need an extra sign convention for real modes; this discretization does
not introduce one. It is a specified discrete operator, not a claim that the
continuum tensor has a unique Nyquist extension.

For exact channel actions on a periodic grid, the tensor is positive
semidefinite, the source has zero mean and

```text
<Θ,L_A Θ> = −Σ_cells [(DΘ)ᵀ A (DΘ) + Σ_d A_dd (N_d Θ)²] <= 0.
```

The implementation uses FFTW with two bounded preparations. It removes each
channel's exact axis-null anchor, prepares two normalized IEEE terms and uses
the existing dyadic FFT owner. Tensor flux is generally nondyadic, so it
undergoes its own anchor removal and two-term preparation before the outgoing
derivative/projector. Final exact line-mean projection removes only a channel's
legitimate null component. This does not assert exact FFTs or exact computed
self-adjointness.

Signed advection, compression, Ohmic heat, hyperresistive heat and conduction
accumulate before the complete pressure rate is rounded once and range-checked.
A partial gradient, tensor flux or heat term need not fit the grid precision;
the complete required rate and every published RK pressure must fit. A nonzero
complete rate that rounds entirely to zero rejects. All pressure stages and
endpoints must be nonnegative; there is no clipping or positivity theorem.

## Timestep and diagnostics

`recommended_dt(st; NB, safety, d_i, Ω_ci)` also uses the maximum effective
conductivity and `nfloor` in the frozen pure-conduction RK4 cap, including full
Nyquist modes. The stage guard uses the actual frozen midpoint density. This
does not bound the combined varying-B/advection/particle/RK map, establish
positivity, or estimate physical resolution error. Refine the actual timestep
and grid for the intended simulation.

`electron_conduction_diagnostics(st, species)` privately deposits fresh density
at the current integer-time positions. It therefore does not depend on the
derived density omitted by plain checkpoint restore. It returns exact rational
integrals of the computed unrounded heat source: internal-energy rate,
`variance_rate=<Θ,Fp_heat>`, nonnegative prepared-channel quadratic dissipation,
and the measured `adjoint_defect=variance_rate+s*channel_dissipation`.
The internal-energy denominator is the canonical `Q(γ_grid)−1`.

These are instantaneous observations with explicit time/phase metadata. They
are not the prior RK stage's work or the rounded complete source. FFT and
preparation error can give a nonzero adjoint defect. Public physical state and
RNG are preserved; the usual nonconcurrent Fourier scratch contract applies.

## Persistence and composition

Policy 3 retains the tensor, diagonal-Nyquist and two-preparation conventions,
original closure coefficients, stored pe, carried E/E_force/B and timestep
phase. Plain raw/archive restore uses the canonical checkpoint owner, with
incoming direction admission against the saved model and B. Same-concrete-type
zero-to-active parameter restore reuses allocated derived workspace. Norm and
gradient grids are scratch and are not serialized. Isotropic policy 2 and
transport policy 1 retain their existing meanings.

Fixed-step ion-ion composition leaves pe fixed during scattering and retains
the outer physical-state, carried-force and RNG rollback/restore contract. It
does not add electron-ion heat exchange or a Coulomb-derived conductivity.

Source/caller review, frozen adoption and compact independent operator/lifecycle/IO
registrations are complete; all new remote execution remains pending. No scientific or
performance acceptance is inferred. State-dependent coefficients, flux limiting,
electron-ion exchange, CGL pressure, full electron inertia, physical boundaries,
implicit transport, MPI and device execution remain open in
[CAP-011](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-011--electron-transport-and-complete-energy-accounting).
