# Flux-limited scalar electron conduction: read-only contract

Status: proposed implementation; no production approval, edits or runtime acceptance.
Source base is immutable `5de6522bbc58a7a8ecd3988b5c0b5dfa587e48e3` plus the released CGL-v1 overlay explicitly recorded in `SOURCE_PINS.json`. This composition is not represented as an integrated Git head or remote main. Root owns the CGL integration and two pending doc precision refinements.

## Verified source boundary

1. `closures.jl:112–262` owns existing constant isotropic policy 2 and field-aligned policy 3; retain their types, arithmetic, zero reductions and meanings.
2. `energy_equation.jl:288–465` owns original `Q(pe)/Q(n)`, two-term dyadic preparation, axis-null removal, complete signed pressure accumulation and one final rate publication. An independently overflowing heat component is not a valid rejection of a representable complete signed rate.
3. `field_aligned_conduction.jl:34–84,93–162` owns reusable real odd-gradient/positive-Nyquist channels, exact computed line-mean projection, two preparations and measured adjoint qualification. No new FFT or derivative implementation is proposed.
4. `extrapolated_leapfrog.jl:939–1023,1074–1203` supplies actual `ptrial`, `Btrial`, frozen actual midpoint density/current and common RK4 stages; `:1208–1241` supplies substep h. The existing carried ideal ion force, finite/nonnegative stages, public rollback and variable-step half-time particle phase remain canonical.
5. `conservation.jl:378–427` uses private fresh charge-density deposition and an instantaneous unrounded-source budget. Do not reinterpret this as accepted-step history.
6. Plain IO `checkpoint.jl:268–335,785–855,1064–1157` validates incoming model/policy and restores carried force without recentering. Collisional IO `collisional_checkpoint.jl:49–183,267–339` uses one prepared restore and actual lock/RNG transaction. `restart.jl:156–168` supplies locked archive capture; both archive-bang delegates already exist.
7. `exact_arithmetic.jl:315–399` owns directed scaled square-root rounding. `thermal_exchange.jl:182–193` already normalizes a positive rational before that root, returning a rational scale beyond grid floating-point range; reuse this algorithm rather than materializing an unnecessary native speed or classical flux.

## Physical law and normalization: settled derivation

Cowie–McKee (1977), Eq. 7, gives the approximate Maxwellian electron saturation model
`q_sat = 0.4 ne kB Te sqrt(2 kB Te/(pi me))`; its 0.4 current correction is assumed, not a universal kinetic theorem. The actual primary PDF's Eq. 7 discussion was checked through the web tool. Orlando et al. (2005), p. 507, Eqs. 2–4, directly uses a harmonic classical/saturated transition. Its hydrodynamic equal-temperature setup is not transferred to this kinetic-ion model.

Propose the isotropic vector model `q = -kappa grad(Theta)/(1 + kappa |grad(Theta)|/q_sat)`, with `q_sat = alpha pe sqrt(Theta/mu_e)` and `Theta=pe/ne` in this package's units. `mu_e=me/mi`; `alpha` is a stored positive model parameter, default expression `0.4*sqrt(2/pi)` evaluated once, validated/converted during closure construction and persisted as its actual stored value (the eventual implementation must document this evaluation/conversion order). This makes the empirical saturation coefficient explicit instead of claiming exact pi-dependent physical data.

From `normalization.jl:5–8,16–22,125–135,160–185,354–368`, pressure scale is `p0=B0^2/mu0=n0 mi vA^2`, velocity scale vA, mass scale mi. Dividing the SI heat flux by `p0*vA` gives exactly the displayed normalized cap. No `weight_to_number` multiplier enters this flux after density/pressure normalization; use the same reference plasma. This is heat redistribution within the existing electron reservoir, with no bath, extra energy reservoir, drag or electron inertia.

## Proposed bounded family and API

Proposed `FluxLimitedConductingElectrons(gamma=5/3; kappa, electron_mass, alpha=0.4*sqrt(2/pi), conductivity=:constant, reference_temperature=1)`.

- CPU periodic FourierGrid D=1,2,3; Float32/64 grids and closure scalar storage first. Original stored closure factors may have a different supported IEEE precision from the grid; keep them as exact Q in the nonlinear flux. Finite original `gamma>1`, `kappa>=0`, `mu_e>0`, `alpha>0`, reference temperature >0; grid-converted gamma must remain >1. Reject a nonzero requested parameter disappearing during closure construction. No standalone conversion of kappa, q_sat, speed, gradient or classical product to grid T is required.
- `conductivity=:constant` means exact stored kappa; `:spitzer` means `kappa(Theta)=kappa_ref*(Theta/Theta_ref)^(5/2)`. This second option is proposed, not approved. Orlando Eq. 3 independently verifies the physical temperature power; no SI numerical coefficient is silently assigned to normalized units. A caller supplies its normalized reference value.
- Nonlinear effective conductivity varies spatially even in constant-reference mode through saturation. Evaluate it inside each divergence, never as pointwise `a*laplacian(Theta)`. Prescribed mutable coefficient arrays/callbacks, anisotropic flux limiting and wall fluxes are separate remaining slices, not implicitly supported.
- Existing scalar transport force/heating, eta/etaH>=0, de2=0, actual density >=nfloor, scalar unit-charge / positive-ion collections, NGP/CIC/TSC remain the supported host. Existing unsupported integrators/backends continue rejection through pressure/transport traits.
- `kappa=0`: exact TransportElectrons arithmetic with no quotient/transform/root/entropy calculation, including mixed cold pressure. Allocate the new closure's derived scratch even at zero for same-concrete-type zero-to-active restart.
- Active nonconstant Theta requires strictly positive pressure at every input and RK trial; a mixed-zero state rejects transactionally. Globally constant Theta, including all-zero pressure, gives exact zero heat before any reciprocal-temperature calculation. A zero-source extension at uniform cold temperature does not claim finite thermodynamic entropy there.

## Discrete full-Nyquist map: settled identities, proposed convention

For each spatial axis let D_d be the stored odd Fourier derivative (Nyquist zero) and N_d the nonnegative multiplier |k_N| on that axis's even Nyquist plane, zero otherwise. Reuse the canonical channel owner. Define the enlarged gradient `G=(D_1 Theta,...,D_D Theta,N_1 Theta,...,N_D Theta)`; omit nonexistent odd-grid N channels.

At each cell set `K=kappa(Theta)^2`, `C=alpha^2*Q(pe)^2*Theta/Q(mu_e)`, `g2=sum(G_c^2)`, and
`a = sqrt(K)/(1+sqrt(K*g2/C))`, `q_c=-a*G_c`.
The first D components are the physical odd heat flux. N components are an explicitly declared diagonal Nyquist completion, not additional resolved spatial directions. Use one shared a for all channels: then `sum(q_c^2)<=C`, so the physical D-component norm also obeys the cap. Limiting each axis independently would not prove a vector cap and is excluded.

The target conductive heat source is `H=-sum(D_d q_Dd)+sum(N_d q_Nd)` and pressure term `R_cond=s H`, where `s=Q(RN_T(gamma_grid-1))`, preserving the existing coefficient owner.
Since `D*=-D`, `N*=N`, and both kill constants, in exact channel arithmetic
`sum H=0`, and `<Theta,H>=-sum a*|G|^2<=0` for spatially varying a. These identities include every Nyquist plane. No off-diagonal Nyquist products are introduced.
As alpha increases without bound at fixed positive stored state, `a -> kappa(Theta)`. For constant kappa, `H -> kappa*sum(D_d^2-N_d^2)Theta = kappa*L_full Theta`, because modewise `kvec_d^2 + k_Nd^2 = kfull_d^2`. This is an ideal operator identity, not bitwise equality between two FFT preparations and old policy-2's direct Laplacian. For Spitzer coefficients the limit is the variable divergence-form operator, not `kappa*L_full Theta`.

Actual map: prepare original quotient channels using existing two-term normalization; form original-factor rational fluxes; prepare the generally nondyadic flux channels a second time; apply canonical channel transforms and subtract computed exact per-line means before addition to the complete pressure equation. Exact `sum H_computed=0` follows from that projection. The constitutive cap alone does NOT certify the second prepared flux: for channel c and its axis-line anchor A, the actual transform input corresponds to `q_tilde_c=q_c(A)+dyadic2(q_c-q_c(A))`. Its error is bounded by `u^2*abs(q_c-q_c(A))`, but unequal local caps prevent a universal relative cap guarantee. Propose a strict exact joint check `sum(q_tilde_c^2)<=C` at every cell, on the actual prepared inputs with removed anchors restored, before source admission. A failed check rejects without clipping/repeated preparation; it is an explicit bounded-domain limit. This certifies the real-space flux input, not FFTW derivative accuracy. Expose the actual preparation through a bounded optional observation seam in the existing channel owner, preserving its old default rows. FFTW/metadata errors and the source mean correction remain qualified. Keep old policies 2/3 literally; this is a new policy.

## Range and coefficient accuracy

Never first form native `kappa*|G|`, `sqrt(Theta/mu_e)`, q_sat or Theta. Their overflow/underflow can be irrelevant to the finite saturated flux.
Use exact positive Q K,C,g2; zero kappa and zero g2 branch before division. Normalize roots through the existing scaled-root algorithm. For Spitzer, `K=kappa_ref^2*(Theta/Theta_ref)^5` is exact Q; constant K is exact without a root approximation.
With directed root enclosures `[kL,kH]` for sqrt(K), `[rL,rH]` for sqrt(K*g2/C), use `aL=kL/(1+rH)`, `aH=kH/(1+rL)` and select aL. Then `aL<=a_exact<=aH`; exact `aL^2*g2<=C` certifies the cap independently of native range. Require `0<aL<=aH` and `aH-aL<=16u*aL`, `u=2^-precision(T)`. This is a proposed deterministic accuracy criterion, not a measured error budget. Unresolved certification rejects before public stores; no epsilon relaxation or arbitrary iteration.
Only the complete signed pressure RHS is rounded/published by `_publish_complete_rate!`; preserve its finite-rate and entirely-lost-nonzero rejection. No standalone conduction-rate fit requirement or pressure clipping.

## Energy and entropy: exact scope of the proposed admission

Use exact volume `V=prod Q(dx)`, `g=Q(gamma_grid)-1>0`. The existing electron energy is `Ue=V*sum(pe)/g`. Conductive input gives `dUe/dt=V*s*sum(H_computed)/g=0` before complete-rate rounding. This is term-specific semidiscrete balance, not exact native RK or total plasma energy conservation.
For fixed actual stage density, scalar entropy can be normalized as `Se=V/g*sum(ne*log(pe/ne^gamma_grid))` up to reference constants. Its conductive directional rate is therefore
`sigma_computed = V*s/g * sum(H_computed/Theta)`.
This Q sum uses actual original positive quotients and the full prepared/projected source; it needs no logarithm approximation. Propose rejecting strictly negative sigma_computed on this new active family only. Do not gate the ideal quadratic dissipation and call it entropy; do not add advection/compression/heating to this conductive gate.
A general collocation reciprocal chain rule is not available: rowwise `sum_j D_ij/Theta_j` is not identically `-sum_j D_ij*Theta_j/Theta_i^2` when an off-diagonal D_ij is nonzero. This algebra explains why positive mobility alone does not establish the required discrete entropy sign. No actual negative-rate runtime witness has been executed.
For a diagnostic use prepared reciprocal channels Z=(D(1/Theta),N(1/Theta)). Define `Pi=V*s/g*sum(a*|G|^2/Theta^2)>=0`, chain defect `-V*s/g*sum(a*G.Z)-Pi`, and measured adjoint defect `sigma_computed+V*s/g*sum(a*G.Z)`. Their exact sum reconstructs sigma_computed, including the actual projection/preparations. This separates physical target, discrete chain defect and numerical adjoint defect.
At final publication let exact complete rate R contain R_other+sH, and delta=Q(RN_T(R))-R. Actual rate-entropy decomposition is `V/g*sum(R_other/Theta) + sigma_computed + V/g*sum(delta/Theta)`. The source gate proves no sign for the final two-term sum after publication, the RK combination, moving density or full field/particle step. Record/check delta at the actual complete-RHS boundary when a budget is requested; never attribute all delta silently to conduction or round conduction separately for admission. An instantaneous conduction-only diagnostic must state that it has no last-stage delta history.

## Canonical integration, timestep, diagnostics, transaction

- Extend the existing pressure trait and complete conduction adapter, passing actual ptrial for the cap. Both ordinary and hyperheat paths already supply actual B/ptrial; preserve their order, current and Ohm owners. No new integrator or source scheduler.
- Keep original scalar gamma/heating and all ion-force gather, priming/predictor and endpoint refresh semantics. Heat flux changes pressure; canonical pressure force and Faraday then see that changed pressure through the existing lifecycle.
- Refactor routing to accept an exact positive reference coefficient for the new closure instead of requiring native `kappa::T`; old methods/default values/expressions retain literal inverse checks. Do not return an invented one-valued conductivity sentinel.
- For constant reference kappa, retain a conservative frozen-mobility pure-diffusion cap with a<=kappa. For optional Spitzer use a directed exact upper bound on max kappa(Theta) at EACH actual RK trial, with actual frozen density; the recommendation may use pe/nfloor. Pass h through a new-closure stage admission seam; old stage bodies remain unchanged. A frozen-coefficient negative-real RK4 cap is not a nonlinear/Jacobian or positivity theorem; retain existing stage rejection.
- New scratch: original quotient, all D/N channels and rational mobility/source/budget roles; enumerate every role for shape and mutual aliases against fields, FFT/metadata, RK, predictor, particles, carried/saved/outer force, rollback and RNG. Derived scratch is never persistence carry and must be overwritten on retry.
- Reuse `electron_conduction_diagnostics`: private fresh deposited density, instantaneous full computed source/mean correction, constitutive and actual prepared cap bounds, variance/dissipation/adjoint error, entropy target/chain/adjoint decomposition and explicit entropy domain. No RNG draws, public mutation or false accepted-step history. Existing field-aligned return remains unchanged.
- Plain public-step rollback already snapshots pressure/fields/particles/clock/force; inspect every new failure after staging and valid retry. Collisional composition reuses existing outer transaction/locks and leaves heat pressure fixed during ion-ion scattering. With optional thermal exchange require exact equality between the new original conduction electron mass and the bound exchange mass, and its existing grid gamma=5/3 requirement. Do not alter other closures' binding.

## Versioned raw/archive restart

Propose new inner `conduction_policy_version=4`, with exact tags for joint odd+diagonal-Nyquist harmonic CM77, original factors, two preparations and projected-source entropy admission. Persist actual closure scalars/law, stored pe, carried Eforce, B/E and existing time/last_dt/particle phase. No quotient, mobility, entropy scratch or FFT plans are saved.
Validate the incoming closure/policy/pressure and ownership graph before detach/publication; no stale destination physical value may select incoming physics. Use canonical model type validation and restored parameter values. Same-concrete-type zero-to-active restore is supported because scratch always exists. Invalid policies/tags/scalars/aliases reject unchanged. Preserve old policy 1/2/3 meanings and explicit reserved-new-field rejection; no generic graph or codec suppression.
Existing collisional outer schemas 1/2/3 and upcoming CGL4 retain literal routes. The new scalar model is distinguished by its inner policy; no gratuitous outer version or alternate restore. Bind every new parameter/law in `_collisional_closure_values`; raw and archive paths delegate to the existing checksum, lock, validation, detach, publication and rollback owners.

## Proposed exclusive implementation scope, not approval

Existing production: `src/electrons/closures.jl`; `src/electrons/energy_equation.jl`; `src/electrons/field_aligned_conduction.jl` (adapter signature plus optional actual-preparation observation seam; old default rows preserved); `src/integrators/extrapolated_leapfrog.jl`; `src/diagnostics/conservation.jl`; `src/io/checkpoint.jl`; `src/integrators/collisional_hybrid.jl`; `src/electrons/thermal_exchange.jl`; `src/utils/exact_arithmetic.jl`.
New production: `src/electrons/flux_limited_conduction.jl`. New feature guide: `docs/FLUX_LIMITED_ELECTRON_CONDUCTION.md`.
The two numerical-owner edits are a literal extraction of the existing positive-Q normalized-root helper into exact_arithmetic and an old-name forwarding wrapper, with a whole-owner inverse; no root algorithm change. Thermal owner additionally has only new-closure mass-match admission. If keeping the existing helper name/caller ownership avoids this move, root must settle that explicitly before implementation. No pi helper, stencil, field solver, deposition, sorting or shared nonlinear changes are proposed.
Root/OpenCode integration owns module include/export, one compact essential test, registration/workflows/global docs/tracker. No collisional IO algorithm edit is currently required: its inner canonical validation/restore already carries new inner policy and bound model. Recheck after actual CGL integration before approving any scope.

## Essential independent remote acceptance

1. Independent direct-DFT small D1/2/3 source oracle, odd/even/mixed grids: actual original quotient, full Nyquist, shared constitutive AND actual prepared vector caps, varying density/pressure, original mass/alpha and constant/Spitzer mobility. Verify divergence form, source mean, target/actual defects and exact entropy decomposition; no operator-only substitute for lifecycle.
2. Zero kappa bitwise Transport reduction, exact constant positive Theta with varying density, all-cold constant extension, mixed-zero active rejection, positive-temperature entropy rejection/acceptance with independent Q dot product. A negative-entropy fixture must first be constructed/proved remotely; it is not currently an executed witness.
3. Weak-limiter convergence toward full stored-k constant-kappa operator including pure Nyquist; strongly saturated cap and isotropy across D components. Check the ideal identity separately from FFT arithmetic qualification.
4. Extreme finite original-factor witnesses where native classical flux/speed overflows or underflows but saturated complete source fits; exact cap and deterministic root interval criterion. Also final complete-rate overflow/loss rejection and real signed-source cancellation rescue, never relaxed tolerances.
5. Charged nonconstant p/B Hybrid init+step with eta and etaH paths, actual changing RK pressure, variable dt/NB, force/predictor phase and fresh diagnostic density; invalid later stage/entropy/alias then bytewise public-state and RNG rollback plus valid retry.
6. Fixed-dt ion-ion wrapper with real RNG and reduced thermal exchange (where mass/gamma/domain match), source/field pressure-force refresh, fixed-config mismatch rejection, outer rollback and continuation.
7. Raw save/load and archive load! of plain and wrapper states, same-type zero-to-active restore, all carried fields/force/particle phase/clock/RNG equality and multiple continuation steps; malformed policy/tag/model/mass/alias rejection unchanged. Retain old policy/schema fixtures.
8. Minimum/latest Julia focused/full unsuppressed tests and JET, formatter and complete-step memory/allocation/FFT counts on GitHub at one immutable adopted head. No local workload has been run and no acceptance is inferred from old green primitive jobs.

## Decisions still requiring approval / end state

Unapproved choices are the new API/type/storage domain, optional Spitzer branch in this first family, conservative lower mobility plus 16u certificate, mixed-zero rejection, strict actual-prepared-flux cap check, strict projected-source entropy gate, complete-RHS budget hook shape and normalized-root extraction. Their equations/proofs above are design evidence, not execution.
No confirmed current-source defect was found in this read-only conduction pass. New flux-limited functionality is absent at the pinned base. Full anisotropic limiting, positive cold-front scheme, arbitrary coefficient fields, implicit/MPI/device/boundary conduction, full drag/inertia/evolving pressure tensors and exact-head science/resource acceptance remain separate full-goal work. Stop here: no implementation permission is inferred from this packet.
