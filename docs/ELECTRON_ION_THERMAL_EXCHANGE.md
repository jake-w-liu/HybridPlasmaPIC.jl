# Electron–ion thermal moment exchange

`ElectronIonThermalExchange` adds two-way exchange between stored scalar electron
pressure and ion thermal motion to `CollisionalHybridStepper`:

```julia
exchange = ElectronIonThermalExchange(
    electron_mass = me_over_mref,
    coulomb_log = [log_ei_species1, log_ei_species2],
)
stepper = CollisionalHybridStepper(
    inner, ions;
    dt, rng, coulomb_log = ion_ion_logs,
    epsilon0, weight_to_number,
    electron_ion = exchange,
)
init!(stepper, ions)
step!(stepper, ions, dt)
receipt = electron_ion_exchange_diagnostics(stepper, ions)
```

The scalar logarithm form expands to one value per ion species. Electron mass,
per-species logarithms and `log_bits` bind at construction to the grid's numeric
type; those stored values define this reduced model. The request and its arrays
are copied. Nonzero logarithms cannot disappear during conversion. The wrapper
owns the bound logarithms and rejects later changes. `log_bits` defaults to
`precision(T)+16`, with an admitted range from that value through 4096.

This is a Maxwellian thermal moment approximation. It does not add a kinetic
electron population, drift friction, pitch-angle scattering, electron recoil
or a finite electron momentum reservoir. Existing kinetic `collide_coulomb!`
supports charged kinetic electron and ion species under its separate contracts.
Prescribed resistivity, hyperresistivity and conduction retain their existing
meanings; this option does not derive their coefficients from a shared Coulomb
law or add another Joule heating term.

## Supported composition and phase

The thermal option admits periodic Cartesian CPU Float32/Float64, one through
three spatial dimensions and all three velocity components, positive-ion
standard `Vector` storage and NGP coupling. Use stored `EvolvingElectrons`,
`TransportElectrons`, `ConductingElectrons` or
`FieldAlignedConductingElectrons` with grid `gamma == T(5/3)`, `de2 == 0` and
actual deposited electron density at least `nfloor`. The remaining collisional
wrapper guards still apply. Initialize stored pressure before `init!`.

It does not admit isothermal or CGL electrons, implicit stepping, non-Cartesian
or open boundaries, CIC/TSC, distributed/device storage, or variable timesteps.
The future CIC/TSC extension needs a positive adjoint overlap construction;
using floor-cell bins would change the present physical discretization.

One fixed step advances the existing field/pressure subsystem, applies existing
ion–ion collisions, applies this thermal source, then recomputes the carried
moments, pressure-dependent force/E and predictor. NGP nearest-node bins match
actual deposition, including periodic wrapping and upward ties. Ion–ion
collisions keep their existing floor-cell bins. `NB` refines field stages only.
Pressure is at the new integer time and ions remain at the velocity half time:
this is a first-order Lie composition, not a synchronous energy-conserving
integrator for the full field–particle system.

## Reservoirs, rate and target

At one grid node, let `W_s` be the sum of ion weights, `V` the cell volume and
`u_s` their exact weighted mean velocity. The normalized reservoirs are

```text
Ui_s = (m_s/2) sum w |v-u_s|²       C_s = 3 W_s/2
Ne   = ne V                        Ce  = Ne/(gamma-1)
Ue   = pe V/(gamma-1)               theta_a = U_a/C_a
```

`theta` means `k_B T` in normalized energy units. It differs from a velocity
variance by the species mass. Use `coulomb_normalization` for `epsilon0` and
`weight_to_number = A`; in fewer than three dimensions the omitted normalized
extents are one. Reservoir energy has no extra factor `A`. The physical rate
contains the real ion density `A W_s/V`:

```text
nu_es = 8 sqrt(2π) A (W_s/V) (qe qs)² log_es /
        [3 (4π epsilon0)² me ms (theta_e/me + theta_s/ms)^(3/2)]
```

Here `qe = -1`. This is the weak-coupling, nondegenerate, nonrelativistic
Maxwellian rate in the SI-permittivity convention; see
[Benedict et al., PRE 86, 046406, equations 1–2](https://link.aps.org/accepted/10.1103/PhysRevE.86.046406).
It is not a dense-plasma or non-Maxwellian rate model.

Rates are frozen from the source's incoming reservoirs. Rational Machin bounds
for π and directed, normalized exact square roots enclose
`beta_s = lambda_s/(1+lambda_s)`, `lambda_s = h Ce nu_es/C_s`.
The selected rational midpoint has relative error at most `4*2^-precision(T)`;
an unresolved bound rejects the step. No standalone native collision frequency
or scale factor is required to fit the grid type. The shared positive target is

```text
alpha_s = C_s beta_s
theta_e* = (Ce theta_e + sum alpha_s theta_s)/(Ce + sum alpha_s)
theta_s* = (1-beta_s) theta_s + beta_s theta_e*
```

This target conserves its exact reservoir energy and does not decrease its
positive Maxwellian moment entropy. That statement applies to the target map. Actual
stored velocities and pressure pass the separate gates below.

## Publication, rejection and diagnostics

Warm ion deviations scale about their exact incoming mean. Cold ions receiving
heat require at least two positive markers and use independent Gaussian seed
components, centered by their exact weighted mean. Every final velocity component
is rounded once from the complete rational mean-plus-scaled-deviation expression.
Zero-weight tracers are untouched. Warm, equilibrium and zero-transfer branches
consume no thermal random draws.

The source bounds each stored momentum change and kinetic-energy rounding
defect using IEEE spacing. Electron pressure pays the **actual** ion kinetic
change, including bulk rounding. The remaining combined-energy defect is bounded
by the final pressure publication's spacing and recorded explicitly. Positive
active target reservoirs must retain positive actual energy. A fixed-precision
signed logarithm enclosure must certify nonnegative actual moment entropy;
negative or unresolved signs reject. This is not a finite-marker kinetic
H theorem or an unconditional stability claim.

An entirely cold active component is an identity. A mixed component containing
a singular zero-temperature electron/ion pair rejects; no artificial cold floor
or infinite rate is substituted. One positive marker cannot receive thermal
energy while retaining momentum. Rate and sparse-cell checks precede every
thermal draw. The wrapper restores public fields, pressure, force, particles,
clock, source report and RNG when a source or later carried-E calculation fails.
Private bins, candidates and other scratch may change. Independent wrappers
require independent inner steppers and particle storage.

`electron_ion_exchange_diagnostics` rechecks the completed numeric receipt against
actual stored endpoint particles and pressure without drawing RNG values. It
returns a detached copy, or `nothing` after initialization. The receipt includes
the physical coefficient enclosure, shared target, scale/seed energy, actual
reservoirs, momentum/kinetic bounds, combined-energy defect, charge-deposition
defect and actual moment-entropy enclosure. It describes the last source split;
it does not report total simulation conservation.

## Restart and remaining work

Enabled thermal exchange uses collisional schema 2 with its complete fixed
configuration, numeric source receipt and source clock, alongside the existing
full fields/pressure/force, particles, fixed timestep and advanced RNG state.
Raw checkpoints and checksum-verified archives use the canonical prepared
restore without replaying exchange. A seed alone cannot reproduce an advanced
RNG stream. The existing same-Julia-runtime/platform requirement remains.
Enabled all-zero logarithms preserve the explicit schema-2 configuration while
skipping thermal cell work and draws. `electron_ion=nothing` keeps the previous
schema-1 and ion–ion-only path. `dt=0` remains a no-op. Initialization resets the
report and source clock while retaining the current RNG stream.

Remaining stages include CIC/TSC overlap exchange; Coulomb-consistent multi-ion
drag, Ohm/Faraday coupling, recoil and drift heating; finite electron momentum
and inertial bulk energy; CGL and variable/higher-order composition; distributed
and device RNG ownership; and non-Maxwellian/dense rates. This reduced option
does not close the complete electron–ion or collisional simulator scope.
