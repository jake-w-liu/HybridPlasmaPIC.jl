# CAP-011 S2: complete represented transport sources

Updated: 2026-11-01. The serial CPU implementation is composed with the shared
IEEE conversion and particle-factor changes. GitHub validation is pending;
CAP-011 remains incomplete. This page owns the implemented source domain and
remaining work. The [transport guide](EVOLVING_ELECTRONS.md#constant-hyperresistive-transport)
provides setup; the [model contract](TRANSPORT_HYPERRESISTIVITY_PLAN.md) derives
the positive heat and periodic discrete balance.

## Implemented scope

`TransportElectrons` accepts constant finite `η ≥ 0`, `ηH ≥ 0`, `de2 = 0`
on serial periodic CPU Float32/64 grids in one to three dimensions. A nonzero
coefficient must stay finite and nonzero at grid precision. Density must remain
at least `nfloor`; existing scalar/collection charge restrictions still apply.
MPI, CAM-CL, Hall-MHD and device transport remain unsupported.

Each RK4 stage computes its own current. Faraday uses
`E = E_force + ηJ − ηH∇²J`; ions gather the carried ideal `E_force`.
Pressure receives `(γ−1)` times positive Ohmic and hyperresistive heat. The
latter includes centered derivative squares and even-axis Nyquist projector
squares, pairing the full Laplacian's periodic global SBP identity. This is a
semidiscrete energy identity, with RK4 field/pressure and second-order coupled
particle time error; it is not exact fully discrete conservation.

## Canonical owners and range handling

- `src/electrons/ohms_law.jl` owns the complete ideal/Ohmic/hyperresistive field
  and source-aware Faraday composition. Cold evaluation combines signed terms
  before conversion, including cancellation of separately overflowing terms.
- `src/electrons/energy_equation.jl` owns the complete signed pressure rate:
  advection, compression, Ohmic heat and Nyquist-completed hyperresistive heat.
  The zero-`ηH` stage shares this owner with ideal evolving pressure; its
  [range scan and essential regressions](SCALAR_PRESSURE_RATE_RANGE_PLAN.md)
  retain ordinary ordering and skip inactive heat and Laplacian channels.
- `src/integrators/extrapolated_leapfrog.jl` selects ordinary or cold stage
  evaluation and retains the existing RK arrays, stage order and transaction.
  Activity checks select a path; they do not prove a nonzero complete source.
- SpectralOperators owns normalized Fourier channels for Laplacian,
  derivative, Nyquist projection and derivative of Laplacian. HPIC consumes each
  borrowed result before calling another operator on that grid. Its physical
  value is `real(values) * input_scale * 2^exponent`; normalized scratch is never
  published as a physical field.
  The registered 0.1.8 consumer floor adds certified axis-null subtraction to
  directional cold preparation. Its exact binary scaling/subtraction gates leave
  the full Laplacian and uncertified fallback unchanged. The dependency passes
  all five checks; exact-revision HPIC GitHub acceptance remains pending.
- `src/utils/exact_arithmetic.jl` owns exact binary rescaling and final IEEE
  conversion. The transport addition retains the shared converter implementation.

The ordinary branch retains the existing arithmetic order; `ηH = 0` takes the
S1 branch, whose force-plus-ηJ publication also completes range-flagged products
before final conversion. Its force remains independently carried. Its signed
pressure rate also requests the shared complete owner on derivative/product
range risk; a separately overflowing heat is not rejected before signed
completion. Risk routing
considers source amplitudes and active curl axes rather
than the norm of an already rounded total E. Cold completion differentiates the
represented ideal source and composes transport terms before final publication.
It uses transient exact accumulators, with no new persistent arrays or schema.
Its FFT coefficients remain numerical represented data, not symbolic exact DFTs.

Required physical `lapJ` and E must be finite, even if later coefficient scaling
could rescue an otherwise overflowing scratch value. Complete required physical
pressure and magnetic rates must fit their storage type. Exact zero and complete
signed cancellation are valid. An exactly nonzero complete rate norm whose every
component rounds to zero is rejected. A finite update rescued only by multiplying
an unrepresentable rate by `dt` requires a different RK representation and remains
unsupported.

### Original-duration RK4 affine combinations

The `_RK4Affine` source candidate in
`src/integrators/extrapolated_leapfrog.jl` owns the scalar stage and final
combinations shared by `_rk4_B!` and `_rk4_Bp!`. For finite represented rates,
an unresolved coefficient, product or sum completes the original `h` expression
before one final rounding through the shared IEEE converter. This includes a
subnormal `h/2` or `h/6` and an overflowing weighted sum.
Ordinary nearest/gradual arithmetic keeps the existing coefficient and operation
order. Finite IEEE inputs outside that arithmetic context use exact completion;
nonfinite operands retain IEEE propagation. Other scalar types retain their
native arithmetic. Existing rate, pressure-candidate and rollback policies
remain in force; no RK rate normalization or persistent array is added.

The `rk4-affine` GitHub essentials suite includes independent rational affine
controls, ordinary scalar/array bit comparisons, and public Ohmic heating at
`dt=2*floatmin(T)`. The heating fixture has finite physical rates but an
overflowing native RK weighted sum; its complete pressure update is finite.
All new runtime, compiler, formatting and resource results remain pending.

The separate `_validated_step_dt` preflight at
[`7254c7e`](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/commit/7254c7e638df5ee6ca62e777397d4070f0cf3766)
rejects a nonzero full timestep when `dt/NB` rounds to zero before evolution;
its runtime acceptance remains pending. This affine candidate leaves that
preflight and CAM-CL unchanged and supports represented subnormal field steps.

Raw `ohms_law!` may leave partial work after a later error. Public `init!` and
`step!` retain their existing rollback, carried force and clock ownership.
Checkpoints and archives retain physical carried force and transport policy 1;
the zero-`ηH` payload meaning is unchanged.

`CollisionalHybridStepper` admits this existing S2 model in its fixed-step,
first-order serial CPU composition. It reuses canonical field/pressure stages
and carried force, holds pressure fixed during ion-ion scattering, and snapshots
physical `lapJ` with all public fields. Schema 1 and transport policy 1 are
unchanged; no additional persistent state is introduced. The
[composition tracker](COLLISION_COMPOSITION_RESTART_PLAN.md#s2-hyperresistive-collision-composition-2026-11-01)
records essential regressions and pending exact-revision GitHub acceptance.

## Verification and remaining implementation

- [x] Compose the full hyperresistive prerequisite and complete-source changes
  with current IEEE and CPU particle owners; inspect their numerical owners,
  scratch lifetimes, stage ordering and publication boundaries directly.
- [x] Include focused permanent independent-DFT/SBP, damping/energy, Nyquist,
  cancellation/range, phase, rollback and restart regressions.
- [x] Historical prerequisite: registered SpectralOperators 0.1.7 passed its four
  dependency CI jobs. [Dependency evidence](https://github.com/jake-w-liu/SpectralOperators.jl/blob/main/docs/SCALED_FOURIER_VERIFICATION.md)
  does not establish HPIC integration acceptance.
- [x] Adopt registered 0.1.8 after all five dependency checks and General publication.
- [ ] Pass HPIC's signed mixed-axis pressure and represented complete-Faraday
  exact-matrix controls on GitHub. The previous 0.1.7 evidence does not certify this update; [bounded consumer scope](SCALAR_PRESSURE_RATE_RANGE_PLAN.md#spectraloperators-018-consumer-adoption).
- [x] Extend complete signed pressure routing to the zero-hyperresistivity
  transport and ideal evolving-pressure stage, with essential cancellation,
  underflow and derivative-rescue regressions prepared for GitHub.
- [ ] Pass the existing minimum/latest Julia package/compiler, formatter and
  MPI checks on the exact composed revision in GitHub. Record failures and fix
  their causes without weakening existing gates.
- [ ] Establish complete public resource costs and suitable hardware results
  before making allocation, peak-memory or comparative performance claims.
- [ ] Complete [original-factor ideal fields and Faraday](IDEAL_FORCE_SOURCE_RANGE_WORK_ORDER.md).
  The bounded scalar `ηH=0` Hybrid route now captures original represented
  constitutive factors on finite product-range events and differentiates their
  exact dyadic source through the existing FFT owner. GitHub verification is
  pending. Hyperresistive completion still starts from represented ideal data;
  normal swallowed variation and general FFT conditioning remain open.
- [x] Integrate original-duration affine RK completion for represented finite
  rates; its [GitHub and remaining normalized-storage gates](SCALAR_PRESSURE_RATE_RANGE_PLAN.md)
  remain separate from source publication.
- [x] Implement the bounded [represented-temperature conduction source](EVOLVING_ELECTRONS.md#constant-isotropic-conduction-source-candidate);
  exact-revision GitHub acceptance and original-quotient completion remain open.
- [ ] Implement broader conduction, equilibration, anisotropy, inertia and backend stages in the
  [capability plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-011--electron-transport-and-complete-energy-accounting).

Earlier bounded candidates, failed harnesses and their raw outputs remain in
the dated audit evidence. Their 504 numerical assertions and 48 ordinary-bit
profiles per runtime belong to those frozen sources. The old 36/46-command local
campaigns are retained history, not prerequisites for this source review.
All new runtime validation runs remotely under the current working policy.
The [audit tracker](DEEP_DEBUG_TRACKER.md) records publication and CI status.
