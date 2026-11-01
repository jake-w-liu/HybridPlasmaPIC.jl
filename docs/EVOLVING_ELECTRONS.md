# Evolving electron pressure

`EvolvingElectrons(γ = 5/3)` gives periodic `HybridStepper` runs an independent
scalar electron pressure. It supports one proton population, ordered ion
collections, and replicated MPI. Set `stepper.fields.pe` before `init!`;
initialization preserves it. The default zero pressure is a valid cold state.
For a single species with a different charge, use a one-element ion collection
to obtain charge-weighted moments; the scalar interface retains proton normalization.

The closure requires finite `γ > 1` and `η = ηH = de2 = 0`, including after
conversion to the grid precision. CAM-CL and Hall-MHD reject this closure.
Constant Ohmic heating is provided by
[`TransportElectrons`](#constant-ohmic-transport), with its hyperresistive
[extension](#constant-hyperresistive-transport). The bounded isotropic
[conduction source candidate](#constant-isotropic-conduction-source-candidate)
adds heat flux; anisotropic pressure evolution remains separate.

## Initialize and advance

```julia
using HybridPlasmaPIC

g = FourierGrid((32,), (2pi,))
ions = ParticleSet{1,Float64}(128)
load_lattice_1d!(ions, 0.0, g.L[1])
set_density_weight!(ions, 1.0, g)

closure = EvolvingElectrons(5 / 3)
st = HybridStepper(g, HybridModel(closure), CIC(), nparticles(ions))
x = (0:g.n[1]-1) .* g.dx[1]
st.fields.pe .= 0.7 .+ 0.04 .* cos.(x)
init!(st, ions)
step!(st, ions, 0.001; NB = 2)

budget = energy_budget(ions, st.fields.B, st.fields.n, closure, g;
                       pe = st.fields.pe)
```

The explicit timestep in this example is for this initial condition. Choose
and refine the timestep for the actual pressure, flow, and magnetic scales.
`recommended_dt` estimates magnetic and gyro limits; it does not bound
electron-pressure advection, compression, or ion-acoustic accuracy.

The equation is

```
∂p_e/∂t = −u_e·∇p_e − γ p_e ∇·u_e,
u_e = u_i − J/max(n, nfloor).
```

Each magnetic subcycle advances B and pressure together through RK4. Every
stage uses its own B, current, and pressure, with charge-weighted ion moments
held at the particle-step midpoint. The outer particle/field coupling is
second order; RK4's fourth order applies to the subsystem with fixed moments.
Changing `dt` or `NB` preserves the pressure state.

After a step, B and pressure are at integer time; particle velocities retain
their half-time level. Recomputing carried E reads the advanced pressure
without evolving it again. Direct `ohms_law!` also reads pressure without
changing it. `electron_pressure!(pe, n, closure)` rejects this closure because
density alone cannot reconstruct the stored pressure.

For fixed-step CPU ion-ion collisions, use
[`CollisionalHybridStepper`](COLLISIONAL_HYBRID.md). Its first-order composition
leaves pressure fixed during the collision substep and retains pressure/RNG
rollback and restart. The evolving extension awaits exact-revision GitHub CI.

## Rejected states and retries

Pressure must be finite and nonnegative on entry, at every RK4 stage, and in
each completed subcycle. Zero is accepted. Any negative value is rejected,
including a small negative roundoff value. There is no clipping or pressure
floor. The explicit spectral method does not preserve positivity for all
data and timesteps; a smaller timestep or better-resolved positive profile
may be necessary.

Coupled scalar pressure rates complete signed advection, compression and
applicable Ohmic heat when derivative/product range checks request cold replay.
Exact cancellation is valid; the complete required rate must fit its storage.
A nonzero complete norm that rounds wholly to zero is rejected. See the
[range contract and pending GitHub verification](SCALAR_PRESSURE_RATE_RANGE_PLAN.md).
The standalone Euler pressure building block keeps its separate contract.

A rejected serial `init!` or `step!` restores public field arrays, particle
positions and velocities, the density-floor counter, and simulation clocks.
Array identities are preserved, so existing views remain valid. Private
workspaces remain scratch. A smaller-step retry can use the restored state.
A zero timestep returns without validating or changing pressure.

Pressure must have independent storage. Public Ohm-law, stepping, and
persistence entry points reject pressure that overlaps their other state or
work arrays before mutation.

MPI ranks must start with identical replicated fields, including pressure,
and matching grid/model configuration. Pressure validation failures are
coordinated across ranks before further moment reductions or migration; the
numerical transaction restores each rank's local public state. MPI transport,
process, and migration failures are outside this guarantee.

## Restart and energy diagnostics

`save_checkpoint`, `load_checkpoint!`, `archive_run`, and distributed
checkpoints preserve pressure along with the existing field and particle
state. A saved evolving-pressure state must contain valid pressure; restore
does not reconstruct it from density or call `init!`. Use the same numerical
environment and stepping inputs for bitwise continuation.

Pass `pe = st.fields.pe` to `electron_internal_energy` or `energy_budget`.
The electron reservoir is `∫ p_e/(γ−1) dV`. Omitting `pe` for this closure is
an error; algebraic closures reject an alternate `pe` argument. Collections
share one electron reservoir.
The stored-pressure integral preserves the ordinary floating-point evaluation
order and uses exact arithmetic when intermediate range loss or possible final
overflow requires it. This protection applies to the electron-pressure term;
[other diagnostic range defects](HYBRID_DIAGNOSTIC_RANGE_PLAN.md) remain tracked.

Energy diagnostics use the supplied particle velocities without a hidden
time-centering kick. Compare reservoirs at a common time when studying energy
convergence. Ideal, periodic, floor-inactive continuum energy conservation
does not imply exact conservation by the explicit integrator. See the
[energy-model balance](HYBRID_ENERGY_MODEL_PLAN.md) for the scope of that
identity and the separate dissipative and inertial model limitations.

## Workspace cost

`EvolvingElectrons` and `TransportElectrons` share pressure stages and transaction
buffers; algebraic closures use none. For D
dimensions, Ncell grid cells and Nparticle local particles, the shared numeric
storage is `(27 + D)*Ncell + (D + 3)*Nparticle` grid-precision values, plus
array headers and bookkeeping. Particle buffers resize with local particle
counts. Transport adds the six force arrays described below. Successful steps
copy public state into the saved buffers; rejected steps also copy it back.

## Constant Ohmic transport

`TransportElectrons(γ = 5/3)` (CAP-011 stage S1) keeps this stored pressure
and adds constant Ohmic transport for serial periodic CPU `HybridStepper` runs,
Float32/64, D = 1–3, with one unit-charge population or an ordered positive-ion
collection. Use `HybridModel(TransportElectrons(γ); η)` with constant finite
`η ≥ 0` and `ηH = de2 = 0`; a nonzero `η` that overflows or becomes zero at
grid precision is rejected. Ohm's law is split as

```
E_force = −u_i×B + (J×B)/n − ∇p_e/n,       E = E_force + ηJ,
∂B/∂t = −∇×E,
∂p_e/∂t = −u_e·∇p_e − γ p_e ∇·u_e + (γ−1) η|J|²,     u_e = u_i − J/n.
```

Faraday's law and `fields.E` use the full field. Ions are pushed by `E_force`
in the first-step priming, every Boris kick and the endpoint predictor, so each
species receives the recoil `−q_s n_s ηJ` while the electrons gain the heat
`η|J|²`. This is a prescribed reduced friction model, not a Coulomb collision
or electron–ion equilibration model. `E_force` is computed by the ideal Ohm
kernel before `ηJ` is added; with `η = 0`, on states whose density stays at or
above `nfloor`, every field, particle and clock value equals
`EvolvingElectrons` bitwise. Each RK4 stage heats with its own current. The
heat uses power-of-two scaling when its ordinary product loses range. The
coupled stage completes the signed ideal/Ohmic pressure equation before final
range classification, so separately overflowing terms may cancel. A genuinely
unrepresentable complete rate rejects the step; see the
[shared scalar pressure plan](SCALAR_PRESSURE_RATE_RANGE_PLAN.md).

The full S1 field combines represented `E_force + ηJ` exactly when the
product overflows or enters the subnormal range, or the native sum is nonfinite.
Ordinary operation order and the `η=0` copy remain. True rounded overflow and raw
nonfinite operands retain their existing caller validation policy. The force is
never reconstructed by subtracting rounded fields. GitHub verification of this
pointwise range update is pending.

```julia
using HybridPlasmaPIC

g = FourierGrid((32,), (2pi,))
ions = ParticleSet{1,Float64}(128)
load_lattice_1d!(ions, 0.0, g.L[1])
set_density_weight!(ions, 1.0, g)

closure = TransportElectrons(5 / 3)
st = HybridStepper(g, HybridModel(closure; η = 0.01), CIC(), nparticles(ions))
x = (0:g.n[1]-1) .* g.dx[1]
st.fields.B[2] .= 0.1 .* cos.(x)
st.fields.B[3] .= 0.1 .* sin.(x)
st.fields.pe .= 0.5
init!(st, ions)
step!(st, ions, 0.001; NB = 2)

budget = energy_budget(ions, st.fields.B, st.fields.n, closure, g;
                       pe = st.fields.pe)
```

Pass the same `η` to `recommended_dt` for its resistive limit; it still does
not bound pressure advection, compression or heating.

**Density.** The deposited density must be finite and at least `nfloor` in
every cell at initialization and at each midpoint, endpoint and predictor
deposition. Floor activation, empty populations and missing particles are
rejected with rollback; the floor never supplies an untracked source. A scalar
`ParticleSet` must have `q = 1`; use a one-element collection otherwise.

**Energy.** Kinetic, magnetic and `∫p_e/(γ−1) dV` energies form the only
reservoirs; with periodic boundaries the continuum total is conserved. The
semidiscrete Ohmic exchange matches to roundoff. With fixed ions the RK4
field/pressure update closes this total at fourth order, and coupled particle
runs close the common-time total at second order (both measured by
`test/test_transport_electrons.jl`). This is not exact discrete conservation.

**State, rollback and restart.** The stepper owns the carried force and its
rollback copy: six additional `Ncell` arrays and no particle buffers.
Rejected `init!`/`step!` calls also restore the carried force; zero `dt` is a
true no-op. Checkpoints and archives require `transport_policy_version = 1` and
the force `Eforce`, validated and restored verbatim without `init!` or
reconstruction from half-time moments. The archive schema remains 4; its
checksum covers the new fields. Earlier checkpoints contain no carried force.

Fixed-step CPU ion-ion collision composition is available through
[`CollisionalHybridStepper`](COLLISIONAL_HYBRID.md) for S1 and its S2 hyperresistive
extension (`eta>=0`, `etaH>=0`, `de2=0`). The field stages supply the applicable
Ohmic/hyperresistive heat; ion collisions hold pressure fixed, then both
carried force and full E are refreshed. The wrapper adds its own force rollback
copy and retains pressure/force/RNG continuation. Exact-revision GitHub validation
is pending; electron-ion thermal exchange, higher-order splitting and backend
composition remain separate work.

MPI, CAM-CL and Hall-MHD reject this closure, and the hybrid stepper has no
device path. MPI ranks agree on the rejection collectively before mutation when
every rank uses an evolving-pressure model; other model mismatches remain
outside the existing MPI contract.

## Constant hyperresistive transport

`HybridModel(TransportElectrons(γ); η, ηH)` extends the same serial periodic CPU
contract to constant finite `η ≥ 0`, `ηH ≥ 0` and `de2 = 0`. Each nonzero
coefficient must remain finite and nonzero at grid precision. The equations are

```
R = ηJ − ηH∇²J,                         E = E_force + R,
∂B/∂t = −∇×E,
∂p_e/∂t = −u_e·∇p_e − γ p_e ∇·u_e + (γ−1)(η|J|² + H_H).
```

For the full Fourier Laplacian, positive `H_H` includes both ordinary
current-gradient squares and the even-axis Nyquist projector squares. Their
periodic sum pairs the magnetic loss. A pure curl-null mode has `J = 0` and no
heat. The [model and discrete balance](TRANSPORT_HYPERRESISTIVITY_PLAN.md)
explain this collocation convention and its boundaries.

Ions retain the carried ideal force. The zero-`ηH` branch keeps S1's arithmetic,
and zero `η` and `ηH` reduce to the ideal evolving-pressure closure. Existing
pressure admissibility, density rejection, transaction, carried-force and
policy-1 restart rules apply. No persistent field or particle arrays are added.
Pass both coefficients to `recommended_dt`; its magnetic stability bound still
does not bound pressure transport or heating.

The [complete-source contract and verification tracker](TRANSPORT_SOURCE_RANGE_PLAN.md)
defines cold completion before separately rounded E, pressure and Faraday
sources lose information. Public `lapJ` and E must remain finite physical
arrays; complete required rates must also be representable. Exact-zero rates
are valid. A nonzero complete required norm that rounds wholly to zero is
rejected. A finite update rescued only by multiplying an unrepresentable rate
by `dt` is outside this representation. Raw `ohms_law!` can retain partial work
on a later failure; `init!` and `step!` restore their public transaction state.

The composed implementation awaits GitHub validation in the linked tracker.
Hyperresistive completion starts from an already represented ideal field. The
separate [original-factor work order](IDEAL_FORCE_SOURCE_RANGE_WORK_ORDER.md)
has implemented scalar/CGL pointwise range completion and a bounded scalar
Hybrid Faraday route for `de2=ηH=0` finite product-range events. Its verification,
normal swallowed-variation admission, broader operators and general FFT
conditioning remain open. Broader conduction, equilibration, anisotropic transport,
complete electron inertia, MPI and device integration remain separate capabilities.

## Constant isotropic conduction source candidate

`ConductingElectrons(γ = 5/3; κ)` adds a declared constant isotropic heat flux
`q_e = −κ∇Te` to the stored scalar-pressure transport closure. Use
`HybridModel(ConductingElectrons(γ; κ); η, ηH, de2=0)` with serial periodic CPU
Float32/64 `HybridStepper`, one unit-charge scalar population or an ordered
positive-ion collection. The current implementation is a source candidate;
its GitHub runtime/compiler and integration acceptance remains pending. This
does not add magnetic anisotropy, flux limiting, equilibration, MPI/device,
CAM-CL or Hall-MHD support. Fixed-step ion-ion collision composition uses the
[existing wrapper](COLLISIONAL_HYBRID.md), with its source and acceptance scope.

S3b prepares the exact quotient `Θ = Q(pe)/Q(n)` from the original stored pressure
and actual deposited density, rather than reciprocal multiplication or a
grid-precision temperature. Active κ requires nearest hardware rounding and
gradual underflow. κ must be finite and nonnegative, and a nonzero κ must resolve
finite and nonzero at grid precision. The derived quotient grid is exact rational
storage; `TransportElectrons` keeps its original layout and semantics. A quotient
need not fit the grid's floating-point type.

For each axis `d`, prepare a fresh residual
`r_d[I] = Θ[I] − Θ[I with index d replaced by 1]`. The removed background is
constant along that axis and its ideal axis Laplacian is zero. Exact constant
and axis-null inputs, including nondyadic or unrepresentable temperatures,
therefore require no transform after input and metadata validation. The original
quotient grid is retained unchanged between axes.

Each nonzero residual is approximated by two normalized nearest-even IEEE
mantissas with their powers of two retained exactly. If their dyadic sum is `a_d`,
the retained preparation error satisfies `|r_d−a_d| ≤ u²|r_d|`, where
`u=2^(-precision(T))`. This is a bounded approximation of a general quotient,
not an exact quotient representation. The canonical dyadic FFT owner then
applies each axis's full stored-wavenumber multiplier `−Q(kfull[d])²`, including
even-grid Nyquist; composing first derivatives would lose that mode.

After the inverse transform and real part, each axis result is projected by
subtracting its exact rational mean along that axis's lines. This projection
has norm at most one and preserves the ideal axis Laplacian. Its rational line
sums, and hence the complete conductive contribution's global sum, are exactly
zero before final rate rounding. Only conduction is projected: the complete
pressure rate and pressure state are not projected. Compatible generated Fourier
metadata is assumed, with `kfull`, `ik` and `kvec` treated as read-only; existing
full-gain metadata admission remains in force even for a null input.

The target pressure equation is

```
∂pe/∂t = −ue·∇pe − γ pe ∇·ue
          + (γ−1) [η|J|² + H_H + κ L(Θ)].
```

`L` here is the ideal Fourier operator using the stored full wavenumbers, not
infinitely precise continuum wavenumbers. The computed projected axis results
enter the existing complete signed advection/compression/heating accumulator.
The existing grid-precision `γ−1` subtraction is retained; this does not complete
an original γ subtraction. Quotient, Laplacian, coefficient product and partial
conductive source need not fit T. Only the complete computed rate is converted:
it must be finite, and a nonzero rate that rounds wholly to zero rejects. Stored
rates and RK stages must fit T even if `dt*rate` could fit. Field, heating and
original-duration RK owners remain unchanged. `κ=0` bypasses quotient preparation,
transforms and projection and retains transport operation ordering.

With `Λ_d=max(kfull[d]²)` and `S_cond=Σ_d Λ_d*||r_d||₂`, preparation alone adds at
most `u²*S_cond` to the operator's l2 error. Exact multipliers and signed sums
add no rational arithmetic error; final rate conversion adds its IEEE rounding.
FFTW is still a finite-precision transform, with no verified universal all-plan
error constant. Independent exact-four-point and general-DFT controls remain
required GitHub qualifications. Severe signed cancellation can magnify relative
final-rate error, so this preparation bound does not certify arbitrary componentwise
accuracy or exact-DFT range classification near overflow/zero boundaries.

For the ideal periodic semidiscrete operator, conduction has zero mean, so it
redistributes the single electron internal energy `∫pe/(γ−1) dV`. At frozen
positive density, the unrounded-temperature equation dissipates
`V = ∫pe²/(2n) dV`. Those are ideal target properties, not a sign guarantee for
the prepared/FFTW/RK map. Final stored rates may have nonzero mean from their
individual rounding. The full spectral Laplacian is not positivity preserving.
Every explicit pressure stage and final candidate must remain nonnegative;
failure rolls back public fields, force, particles and clocks without clipping.
No exact fully discrete energy, variance or bitwise mean conservation is claimed.

Each nonzero resolved substep `h = dt/NB` must satisfy the strict scoped bound
`h*(γ−1)*κ*Kmax²/min(n) < 2.785`, where `Kmax²` sums the largest squared full
mode on each axis. This is a frozen-density pure-conduction RK4 guard, not a
combined stability or positivity theorem. `recommended_dt(stepper; NB, safety,
d_i, Ω_ci)` adds a conservative cap using `nfloor` and rounds downward to grid
precision, respecting the strict bound after resolved substep division even at
`safety=1`; a cap below the public `2*floatmin(T)` full-step resolution or an
unresolved zero substep returns zero. The existing grid-based overload is
unchanged. Advection/compression and nonlinear coupled accuracy still require
timestep refinement.

Checkpoint/archive state saves κ in the concrete closure, pressure and carried
force verbatim, existing `transport_policy_version=1`, and required
integer `conduction_policy_version=2`, including when κ=0, under unchanged archive
container 4. Missing, policy-1, boolean or differently typed versions reject
before mutation; there is no hidden legacy mode or migration. Quotient and
transform scratch is not persisted and is recomputed at each active stage; the
quotient scratch array identity survives restore. A matching concrete conducting
destination may replace its κ, including zero with positive κ. Legacy transport
payloads retain their meanings. Fixed collisional wrappers retain their stricter
saved-model match and outer schema 1.

Essential candidate tests cover independent DFT stages, variable-density
constant-temperature/axis-null and Nyquist controls, original quotient variation
and range rescue, signed source range cancellation,
RK4 decay refinement, scoped timestep/positivity rejection, full rollback,
zero-κ legacy reduction, malformed checkpoint policies and unchanged resume.
They await GitHub acceptance on the composed source. Complete-step BigInt/FFT
cost and broader conditioning are unmeasured; no optimized-cost claim is made.

Historical S3a used `round_T(Q(pe)/Q(n))`, required a representable nonzero
physical temperature and saved conduction policy 1. Its dated source/CI evidence
remains history; S3b replaces that numerical policy rather than silently
reinterpreting its checkpoints.
