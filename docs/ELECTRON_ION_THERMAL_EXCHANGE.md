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
standard `Vector` storage and NGP, CIC or TSC coupling. Use stored `EvolvingElectrons`,
`TransportElectrons`, `ConductingElectrons` or
`FieldAlignedConductingElectrons` with grid `gamma == T(5/3)`, `de2 == 0` and
actual deposited electron density at least `nfloor`. The remaining collisional
wrapper guards still apply. Initialize stored pressure before `init!`.

It does not admit isothermal or CGL electrons, implicit stepping, non-Cartesian
or open boundaries, distributed/device storage, or variable timesteps.

One fixed step advances the existing field/pressure subsystem, applies existing
ion–ion collisions, applies this thermal source, then recomputes the carried
moments, pressure-dependent force/E and predictor. NGP nearest-node bins match
actual deposition, including periodic wrapping and upward ties. CIC/TSC use
the distinct shape-overlap frame construction below. Ion–ion collisions keep
their existing floor-cell bins. `NB` refines field stages only.
Pressure is at the new integer time and ions remain at the velocity half time:
this is a first-order Lie composition, not a synchronous energy-conserving
integrator for the full field–particle system.

## Reservoirs, rate and target

For the NGP policy at one node, let `W_s` be the sum of ion weights, `V` the cell volume and
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

## CIC/TSC overlap and edge ordering

CIC/TSC select thermal policy 2. Ion thermal frames are disjoint periodic floor
cells, with physical `W_f=sum w`, `C_f=3W_f/2` and one Maxwellian moment
temperature per species/frame. They are not point deposits at the frame center.
The canonical stored shape factors give positive incidence
`H_fj=sum_(p in f) w_p S_pj`; periodic duplicate images are combined before
making an edge. Products are retained exactly, including small positive support.
The same incidence drives the node-to-frame thermal gather and its transpose.

Stored floating-point shape weights need not sum to exactly one as rational
values. The source preserves those actual weights and reports
`F_f=sum_j H_fj-W_f`, keeping the physical ion capacity unchanged. It also reports
`D_j=ne_j V-sum_f q_f H_fj`. Their global balance is
`sum_j ne_j V-sum_f q_f W_f = sum_j D_j+sum_f q_f F_f`.
The source does not normalize these measured defects away or replace actual
charge/current/pressure-force deposition and gathering with frame-center values.

Each edge uses the physical rate above with real ion density `A H_fj/V`, nodal
electron temperature and frame ion temperature. All rates are frozen from the
source's incoming reservoirs. Their common conductance is `G_fj=Ce_j*nu_fj`:

```text
dUe_j/dt = sum_f G_fj (theta_f-theta_ej)
dUi_f/dt = sum_j G_fj (theta_ej-theta_f)
```

Edges run in fixed `(species order, frame linear index, node linear index)`
order. Each positive conservative backward-Euler pair map uses the current
staged temperatures and the frozen certified beta:

```text
alpha = C_f beta
theta_e* = (Ce theta_e + alpha theta_i)/(Ce+alpha)
theta_i* = (1-beta) theta_i + beta theta_e*
```

After each edge the entire disjoint ion frame is materialized about its current
exact mean. Its **actual** kinetic change is paid at that edge's electron node,
and actual reservoirs feed the next edge. Per-edge pressure rounding defects
sum to the global source-energy defect. Tight component momentum/frame-energy
and actual moment-entropy gates apply to every edge. Connected entirely cold
components are identities; a mixed component with an incoming singular cold
edge rejects before any thermal draw. A cold frame is seeded only at its first
receiving edge, then follows the warm path.

This is a fixed-order first-order split of the frozen overlap network. It is
not a global backward-Euler solve, a symmetric integrator, or a finite-step
equivalent of NGP's simultaneous multi-ion star. Physical coefficient and native
rounding errors remain explicit. The frame-Maxwellian approximation does not
resolve an arbitrary particle distribution's exact local collision operator.

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
by its pressure publications' spacing and recorded explicitly. Positive
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

NGP retains its policy-1 per-node `cells` report. Policy-2 CIC/TSC reports contain
`initial_frames`, exact `incidence`, ordered `edges`, `density`,
`initial_pressure`, `frame_partition`, `charge_defect` and source-wide energy/
momentum budgets. Edge `before`/`after` hold exact statistics of actual stored
components; `target` is ion thermal energy and `electron_target` is electron
temperature. `cold` marks an entirely cold connected component. The beta triple
is `(lower, selected, upper)` frozen from incoming reservoirs. Trace validation
checks node/frame chains and final endpoint state without rerunning collisions.
Historical intermediate spacing bounds use derived IEEE/Cauchy caps based on
that edge's weight and kinetic energy; actual evolution retains its tighter
component bounds. This is numeric history consistency, not proof of unsaved
intermediate particle trajectories.

## Restart and remaining work

NGP thermal exchange retains collisional schema 2; CIC/TSC overlap uses schema 3
with fixed shape/frame/incidence/order/publication tags. Each carries its complete fixed
configuration, numeric source receipt and source clock, alongside the existing
full fields/pressure/force, particles, fixed timestep and advanced RNG state.
Raw checkpoints and checksum-verified archives use the canonical prepared
restore without replaying exchange. A seed alone cannot reproduce an advanced
RNG stream. The existing same-Julia-runtime/platform requirement remains.
Enabled all-zero logarithms preserve the explicit policy/schema configuration while
skipping thermal cell work and draws. `electron_ion=nothing` keeps the previous
schema-1 and ion–ion-only path. `dt=0` remains a no-op. Initialization resets the
report and source clock while retaining the current RNG stream. Policies/shapes
cannot be interchanged or silently migrated on an advanced restart.

Remaining stages include Coulomb-consistent multi-ion
drag, Ohm/Faraday coupling, recoil and drift heating; finite electron momentum
and inertial bulk energy; CGL and variable/higher-order composition; distributed
and device RNG ownership; and non-Maxwellian/dense rates. This reduced option
does not close the complete electron–ion or collisional simulator scope.
