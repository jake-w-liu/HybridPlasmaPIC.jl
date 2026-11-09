# CGL electrons with kinetic ion–ion collisions

`CollisionalHybridStepper` composes the existing algebraic `CGLElectrons` field
model with the real weighted kinetic ion–ion Coulomb operator:

```julia
model = HybridModel(CGLElectrons(0.4, 0.6, 1.0, 1.0))
inner = HybridStepper(g, model, CIC(), sum(nparticles, ions))
inner.fields.B[3] .= 1
stepper = CollisionalHybridStepper(
    inner, ions;
    dt, rng, coulomb_log, epsilon0, weight_to_number,
)
init!(stepper, ions)
step!(stepper, ions, dt)
budget = energy_budget(ions, inner.fields.B, inner.fields.n, model.closure, g)
```

Here `g` is a periodic Fourier grid and `ions` is an ordered vector of positive
ion populations. Seed the magnetic field before `init!`. Use the existing
`coulomb_normalization` conventions for the collision coefficients and number
normalization. The wrapper owns a copy of the explicit Xoshiro or MersenneTwister
stream; independent runs require independent inner steppers and populations.

The admitted composition uses Float32/Float64 grid and standard `Vector`
particle storage, one through three spatial dimensions and three velocity
components, NGP/CIC/TSC coupling, fixed counts/charges/masses, and fixed positive
timestep. Set `η=ηH=de2=0`. A zero timestep is a no-op; zero Coulomb logarithms
preserve the plain inner physical evolution and consume no collision draws.
The scalar `electron_ion` option must remain `nothing`, including when its
requested logarithms would be zero.

## The pressure model and its limits

The closure retains two global coefficients:

```text
p_perp = Cperp n |B|
p_parallel = Cpar n³/|B|²
P_e = p_perp I + (p_parallel-p_perp) b b
```

It does not store or advect independent perpendicular/parallel pressure fields.
Equal reference pressures specify an isotropic reference in the ideal
parameterization; they do not make the whole law isotropic. For positive stored
coefficients and positive density/field, constitutive isotropy requires
`B³/n²=Cpar/Cperp`. The ideal equal-reference ratio is `B0³/n0²`; actual stored
coefficients retain constructor rounding. A tensor field reduces to scalar
pressure only when its two pressures agree throughout.

Original invariant coefficients must be finite and nonnegative. Their standard
Float16/32/64 or BigFloat precision remains independent of grid precision; no
standalone conversion to the grid type narrows the model. The existing pressure
owner first converts density and field magnitude to the coefficient type and
then computes the pressure before field conversion. Its range/rounding behavior
is retained. Mixed grid/coefficient precision does not imply a universal accuracy
or range guarantee. BigFloat calculations require the same active arithmetic
environment for reproducible continuation; the wrapper does not reconfigure
global MPFR precision, rounding or exponent limits.

At every magnetic RK stage the existing field owner recomputes the full
three-component `div(P_e)` from actual trial B. It does not freeze the field
direction for a whole step. The force uses its existing regularized magnitude
`max(hypot(B), T(1e-12))`, with direction B divided by that magnitude. This
composition does not change that floor or add a pressure limiter.

Pressure anisotropy can drive firehose or mirror instabilities. CGL and kinetic
mirror thresholds differ; an ordering such as `p_parallel>p_perp` alone does not
establish stability. These MHD results do not supply a universal threshold for
this kinetic-ion/Hall implementation. The electron closure adds no FLR/Landau
stabilization or numerical stability guarantee. See
[Snyder, Hammett & Dorland (1997)](https://w3.pppl.gov/~hammett/gyrofluid/papers/snyder_pop97.pdf)
and [Parker (1958)](https://journals.aps.org/pr/abstract/10.1103/PhysRev.109.1874).

## Phase, scattering and energy

The existing cycle is preserved:

1. Push ions and deposit midpoint charge moments.
2. Advance B with its trial-field CGL force at each RK stage.
3. Scatter the actual half-level ion velocities at the new integer positions.
4. Recompute integer-position moments, CGL Ohm and predictor/recentered carried E.

The collision interval is the full fixed field timestep. `NB` refines field
stages only. Positions and B remain at integer time; velocity remains at the
half time. This is a first-order Lie composition. It adds neither a second
field integrator nor a new source scheduler.

Ion–ion scattering holds positions, weights, charges, masses and B fixed.
Fresh deposition at those positions therefore leaves n and the constitutive
electron tensor unchanged during that subflow. Ion velocity anisotropy can
relax through the existing three-component rotations. Equal-weight pair
conservation and native roundoff retain the collision operator's contract;
unequal-weight conservation remains statistical. No ion energy drift is paid
from an electron bath, and no electron isotropization or thermal exchange occurs.
Changed ion moments still affect the carried E and subsequent field evolution.

The electron energy is the gyrotropic internal energy
`integral(p_perp+p_parallel/2) dV`. Use the B-taking diagnostics shown above.
`HybridFields.pe` is CGL pressure/divergence scratch, not a scalar heat reservoir;
passing it as an alternate `pe` to the CGL budget is rejected. Float32/Float64
coefficient variants retain the existing checked native energy/replay path;
Float16/BigFloat variants retain the generic dispatch and its context/range
qualifications. The energy diagnostic uses raw B and keeps the parallel-pressure
singularity at zero field visible, even though the force uses a floor.
Neither adding this reservoir nor composing ion collisions establishes exact
full-step energy conservation, a kinetic entropy theorem or stable evolution.

## Failure and complete restart

The outer transaction owns full-sized anisotropic `pforce` and `cgl_dp` buffers,
even when both coefficients are zero. It snapshots every public field/cache,
particle position/velocity, clock and complete advanced RNG state. A rejected
initialization, field stage, collision or carried-force calculation restores
those values in place. Private bins/RK/FFT scratch may change. Independent
rollback roles and RNG recurrence/cache buffers are checked; stock read-only
RNG seed/advance metadata may be shared by its canonical copy operation.

CGL composition uses collisional schema 4. Required policy metadata is
`cgl_policy_version=1`, the grid-type `cgl_force_floor=T(1e-12)` and
`cgl_force_operator=:global_invariants_trial_B_tensor`. The payload keeps the
original model/coefficient type and values, shape, fixed collision coefficients
and timestep, full ordered particles/provenance, public fields and caches,
source clock/phase and advanced RNG. Caches include n/ui/B/E/J/pe/gradp/ninv/
lapJ/pforce/cgl_dp and the floor counter; scalar-pressure positivity rules are
not imposed on CGL scratch.

`save_checkpoint`/`load_checkpoint!` and `archive_run`/`load_archive!` use the
existing locked, prepared restore; archive loading checks the checksum first.
Restore copies prepared values without initialization, field evolution, FFT or
collision replay. Keep the same supported runtime/platform and arithmetic
environment. Reinitialization instead resets the clock while preserving the
current RNG stream. A seed alone is not an advanced-stream checkpoint.

Schema 4 cannot be interchanged with scalar schema 1, NGP thermal schema 2 or
CIC/TSC thermal schema 3. Model, shape, normalization, logarithms, counts/charges/
masses and policy must match; malformed metadata or anisotropic storage rejects
before publication. Existing scalar and scalar-thermal formats retain their
own behavior.

Remaining stages include spatially evolving two-pressure electron invariants,
electron collisional anisotropy relaxation/heat/drag, consistent dissipative or
inertial CGL energy accounting, kinetic electron closures, variable/higher-order
composition and distributed/device storage/RNG. The source implements the
existing ideal algebraic CGL plus ion–ion composition. Exact-head numerical,
compiler/format/CI, scientific and measured resource acceptance remains pending.
