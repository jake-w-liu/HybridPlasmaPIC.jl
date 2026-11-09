# Coupled toroidal annulus hybrid

`ToroidalAnnulusHybridPIC` is a bounded circular-annulus family with full-f
kinetic ions, primary magnetic flux and evolving scalar electron pressure.
Every classical RK4 stage deposits all populations and recalculates their common
fields. The radial boundary represents stationary equipotential electrical
contacts with prescribed incoming electron-pressure reservoirs.

This family excludes the coordinate axis. It does not complete CAP006A's
solid-axis topology, material ion-wall events, closed thermal wall, toroidal
equilibria/imports, general maps or distributed/device backends. Compiler,
formatting, runtime, numerical, boundary-stability, convergence and resource
acceptance are **unverified** until the frozen implementation passes its
independent remote gates. Source implementation and exact work accounting alone
do not establish scientific acceptance.

## Setup and synchronous phase

```julia
using HybridPlasmaPIC

T = Float64
n = (5, 5, 5)
M = (2, 1, 1)
mesh = ToroidalAnnulusGrid(n; R0=3, r_inner=0.5, r_outer=1, T=T)
N = (n[1] - 1) * n[2] * n[3] * prod(M)
ions = ParticleSet{3,T}(N; q=1, m=1)
load_toroidal_annulus_cold!(ions, mesh;
    number_density=1, velocity=(0,0,0), particles_per_cell=M)
species = [ions]
beta0 = toroidal_annulus_vacuum_flux(mesh; strength=1)
simulation = ToroidalAnnulusHybridPIC(mesh, species;
    closure=EvolvingElectrons{T}(T(5)/T(3)), beta0=beta0, pe0=T(0.1),
    pressure_reservoir=(T(0.1), T(0.1)))
init!(simulation, species)
step!(simulation, species, T(1e-4))
report = toroidal_annulus_diagnostics(simulation, species)
save_checkpoint("annulus.bin", simulation, species)
load_checkpoint!(simulation, species, "annulus.bin")
archive_run("annulus-run.bin", simulation, species; rng_seed=0)
load_archive!(simulation, species, "annulus-run.bin")
```

This is an unexecuted API illustration. Its parameters and timestep are not a
validated equilibrium or timestep recommendation. Measure source, force,
boundary and refinement errors for the physical problem being studied.

The admitted storage is ordinary CPU Float32/Float64 arrays, with nearest
rounding and gradual underflow. At least three radial nodes and odd periodic
sizes of at least three are required. The logical radial interval is closed;
the two angular coordinates use normalized turns in `[0,1)`.

Bound `ParticleSet.x` holds logical `(xi,u,v)` positions. Its three velocity
components are physical Cartesian components. Positions, velocities, primary
`beta` and electron pressure share the same integer global time. There is no
half-step priming or field subcycling. Four coupled RK4 stages evolve all these
primaries; a fifth endpoint refresh synchronizes derived fields and admission
receipts without entering the four-stage work quadrature. Positive macro
durations may vary. Zero duration validates the phase and is a no-op; negative,
nonfinite, disappearing nonzero or clock-nonadvancing durations reject.

Call the cold loader before binding a simulation. Ordered species owners,
component identities, counts, positive q/m and nonnegative physical number
weights are fixed after binding. Empty individual populations and zero-weight
tracers are allowed, but the complete collection must deposit positive actual
density at every node. No density floor or vacuum fill is supplied. Explicit
`init!` admits changed evolving primary data and resets the clock and work ledger
within the same fixed configuration. Changing geometry, source inventory,
closure, boundary data or tolerances requires a matching new model.

## Coordinates, metric and cold loading

The physical map, with `r=r_inner+xi`, is

```text
theta = 2*pi*u; phi = 2*pi*v
x = ((R0+r*cos(theta))*cos(phi),
     (R0+r*cos(theta))*sin(phi), r*sin(theta))
```

The geometry requires `0 < r_inner < r_outer < R0`, positive distinct represented
radial nodes and finite admitted transforms. `to_cartesian(mesh, logical_tuple)`
and `to_logical(mesh, Cartesian_tuple)` use this annular coordinate contract.
There is no exact floating-point map/inverse identity or tolerance that moves an
out-of-domain particle into the radial interval.

The canonical `MixedShockOperators` owns SBP-(2,1) radial nodes and derivatives,
with periodic transverse derivatives and logical lengths `(ell,1,1)`, where
`ell=r_outer-r_inner`. Bounded-turn trigonometry and the represented angular
factor `tau=T(2*pi)` define the actual metric. With
`R=R0+r*cos(theta)`, its scale factors are `h=(1,tau*r,tau*R)`,
positive volume Jacobian `J=prod(h)` and diagonal metric `Gii=hi^2`.
Products and divisions retain the actual represented factors through exact
rational completion before final materialization.

For coordinate order `(r,theta,phi)`, the basis is left-handed, so the curl
orientation is `sigma=-1`. Positive physical volume is separate from this
orientation. Primary `beta^i=J*B^i`; its covariant magnetic component is
`B_i=Gii*beta^i/J`. Logical vector components must not be interpreted as
Cartesian components. The public B/current/electric/velocity fields are
Cartesian. `simulation.fields.Eohm` retains the raw Ohm result;
`simulation.fields.E` is the controlled Cartesian field called Eused below, and
`simulation.fields.Ecov` is its separately retained covector called eused below.
Faraday consumes that covector.

Nodal logical quadrature is `Hradial[i]*du*dv`; physical volume is that measure
times the same actual J. Radial face area is `h2*h3*du*dv`. Metric primary
magnetic energy and rounded Cartesian snapshot energy can differ; the
difference is observed rather than assumed to cancel. Collocated flux arrays
have no independent staggered face degrees of freedom or asserted exact local
finite-volume Gauss law.

For `particles_per_cell=M`, the cold loader requires
`(Nr-1)*Nu*Nv*prod(M)` existing particles. It uses adjacent actual logical node
coordinates, with the last angular cell ending at one turn. Markers occupy a
quiet lattice of strict subcell interiors; a coordinate that rounds onto a
subcell edge rejects. At each actual prepared marker it assigns physical number
weight `number_density*J(marker)*logical_cell_volume/prod(M)` and the supplied
constant Cartesian drift. Each complete weight is materialized once. All input,
storage, alias, position and weight checks precede public stores. IDs and tags
are preserved, and retained charge and mass must be finite positive values.
Old positions, velocities and weights may be stale because the complete load
replaces them. The loader neither assigns nor promises uniform deposited density.

`toroidal_annulus_vacuum_flux` returns independent primary arrays for the analytic
profile `B=(strength/R)*e_phi`. Only component three is seeded, from the same
nodal factors as `beta3=J*strength/(h3*R)`. The phi-independent analytic profile
provides a baseline seed; its actual differentiated current, force and discrete
equilibrium residual must be measured. Cold stationary ions and constant
pressure do not constitute a solved toroidal equilibrium.

## Coupling, pressure and radial contacts

The units are coherent Omega_ci hybrid normalization, with elementary charge and
mu0 equal to one. Weights are physical particle numbers. Electron density is
`n=sum(q*w*S)/V`, Ji is the charge-weighted Cartesian ion current, `ui=Ji/n` and
`ue=ui-Jcurrent/n`. Only the ideal nonrelativistic, massless-electron closure
`EvolvingElectrons` with gamma greater than one is admitted. Pressure is stored
and evolved; it is not reconstructed from density. There are no resistive,
hyperresistive, electron-inertia, conductive, collisional or active RNG channels.

Canonical CIC coefficients use folded radial endpoint nodes and periodic angular
indices. Exact whole-collection deposition and complete Cartesian vector gather
use the same stored tensor coefficients. Their finite-sum particle/grid work
identity uses actual physical volumes. Materialized moments/gathers have
measured signed work and charge defects; no deposit renormalization or charge
repair is applied.

Each private stage refreshes actual density/current, B and its covariant form,
canonical mixed curl current, pressure gradient, ideal Ohm field, electron
velocity and mapped pressure divergence. Canonical differentiated-source-risk
rejection applies before any boundary override, including initialization,
endpoint and restore. A finite local electric field does not certify its
differentiated source. Fresh beta must pass the actual physical divergence gate;
no projection, filter or clipped tolerance repairs a rejected field.

Each radial face is a stationary equipotential electrical contact. The complete
Ohm result is retained as Cartesian `Eohm`. At a face, consumed covector `eused`
keeps its normal component and sets both tangential components to stored zero.
The reconstructed Cartesian `Eused` drives every gathered particle force;
Faraday consumes that same intended field through `eused`. Initial normal-beta
face data are frozen and checked exactly; a changed trace or nonzero consumed
normal-beta rate rejects rather than being reset. These traces may be nonzero.
The contact does not impose a Meissner condition.

`pressure_reservoir=(p_inner,p_outer)` supplies fixed nonnegative incoming
pressure data. With outward sign `s=-1` at the inner face and `s=+1` at the outer,
actual electron contraflux `U^r=J*(ue dot er)/h1`, and radial norm entry Hf,

```text
a_in = max(-s*U^r, 0)
nu = a_in/(J*Hf)
S_p = -nu*(pe-p_res)
```

Outflow receives no pressure penalty; zero radial electron flux activates
neither reservoir. The canonical complete primitive pressure rate and exact
incoming SAT are added before one final rounding. Pressure is never overwritten
with p_res or clipped. Every stage and endpoint must retain finite nonnegative
pressure. Exact `h*nu <= 1` admission limits only the local SAT relaxation
stiffness; it is not a full coupled stability bound or an RK4 positivity proof.

The contacts may exchange electrical and thermal energy with the domain; they
are not an isolated adiabatic material wall. Every actual logical radial
particle stage must lie in `[0,ell]`. A stage outside rejects the whole step.
Folded CIC is a shape rule, not specular reflection, a curved-wall event solver
or a guarantee that no crossing occurred between evaluated stages. No ion
injection, removal or reflection is provided. The deposited normal-current trace
is observed separately from any future physical particle-wall event count.

## Work observations, rollback and complete restart

Primary reservoirs use exact rational expressions of admitted stored operands:
particle kinetic K, represented-metric magnetic UB and thermal
`Ue=sum(V*pe)/(gamma-1)`, with `H=K+UB+Ue`. Every actual RK stage independently
observes particle force work, grid-ion/current work, magnetic/thermal power,
pressure work and physical outward Poynting/enthalpy powers. Cartesian rotation,
metric reconstruction, canonical curl, shape/gather and pressure assembly defects
remain signed observations. Actual face Poynting power is retained even when
the prescribed tangential covector is exactly zero.

Electrical controller power is the face sum of
`V*(Ji-Jcurrent) dot (Eused-Eohm)`. Thermal controller power is the face sum of
`V*S_p/(gamma-1)`. Positive power supplies domain energy; negative power removes
it. Physical enthalpy flux uses the evolved face pressure, not p_res. Per-face
controller powers, physical fluxes, incoming relaxation coefficients and complete
pressure-rate assembly defects remain visible.

With the four classical RK4 weights, accumulated outward work Wout, controller
work Wcontroller, stage residual Dstage and finite-step residual DRK satisfy

```text
Hnew-Hold + Wout-Wcontroller = Dstage+DRK
Hcurrent-H0 + Wout_cumulative-Wcontroller_cumulative
    = Dstage_cumulative+DRK_cumulative
```

These exact ledger relations account for the actual stored phase and consumed
stages. They are not zero energy-error, momentum, angular-momentum, stability or
convergence theorems. The fifth endpoint receipt is separate from stage work.

All trial stages and fields are private. Complete expected validation, report
and ledger checks precede fixed-length ordinary publication. A separate committed
primary shadow detects external evolving-state changes; failed scratch work does
not replace it. Rejected init/step/load calls retain the preceding public phase,
fields, provenance, clock, report and ledger. Exclusive operation ownership
prevents reentry; callers must not mutate live state concurrently.

The distinct closed schema `:toroidal_annulus_hybrid` records fixed geometry,
operator/shape/component/source policies, closure, tolerances, pressure
reservoirs, frozen normal-beta data, every ordered particle leaf and provenance,
primary beta/pe, complete derived fields including Eohm/E/Ecov, synchronous
clock, endpoint receipt, complete last four-stage report and cumulative work.
Private FFT plans and warm scratch are rebuilt, not serialized.

Restore checks original incoming aliases before detachment, exact fixed
configuration, fresh source/divergence/boundary receipts and the controller-aware
ledger before publication. Saved fields restore verbatim after validation;
restore does not call init!, project beta, overwrite pressure, recenter particles
or silently change reservoirs. Old destination evolving state may be stale;
fixed owners, components and configuration must still match. Raw checkpoint and
checksummed archive restoration share the same physical restore. The archive's
rng_seed is provenance for this deterministic family.

Fixed configuration includes immutable tuples of the actual represented metric,
trigonometric, SBP-norm and transverse-symbol factors. Current endpoint contact
powers and incoming relaxation coefficients are checked exactly against saved
fields. Refreshed endpoint powers and defects use the configured relative
allowance scaled by both receipts' complete raw powers and physical face powers;
native derivative refresh does not require bit-identical cancelled residuals.
Historical RK stages validate their own saved algebra and intervals.

## Required independent acceptance

The following gates remain unverified until exercised remotely on the exact
frozen revision, in both admitted precisions:

- Map/inverse, orientation and nonconstant-J measures; independent marker-weight
  and nodal-flux references, including range and unresolved-coordinate rejection.
- Nonuniform-weight, multispecies particle/grid work; both-face physical/logical
  curl and pressure work; independently assembled incoming/outgoing SAT signs,
  zero flux and outflow-reservoir independence.
- Genuine nonconstant particles+beta+pe evolution compared with an independent
  direct small-reference RK4 calculation; both Eused consumers and nonzero
  electrical/thermal controller work; observed normal-flux preservation.
- Zero-step, pressure/radial/source/divergence/stiffness rejection and complete
  rollback; exact cumulative ledger; raw/archive continuation into a stale
  destination; original alias, fixed-boundary, schema and checksum rejection.
- Independent mesh/time refinement of a compatible coupled annular wave or
  equilibrium perturbation, boundary-mode growth/reflection and incoming-pressure
  response, plus decay of the electrical contact reaction in its smooth
  conductor limit. An unexplained growing boundary mode or nonconvergent
  reaction fails boundary acceptance even when the ledger closes.
- Minimum/latest supported Julia, package/compiler/format checks and measured
  memory, mapping, gather, FFT and complete-timestep cost. Source inspection
  supplies no speed or allocation claim.

Full CAP006A still requires a solid-axis owner with regularity/parity and
axis-crossing particles; actual curved-wall events and closed field/thermal-wall
physics; finite-temperature/flow equilibria and imported references; general
nonorthogonal maps/topologies; and distributed ownership and device operators.
An annular cutout is not axis support and does not close those requirements.

Source entry points are the [mesh and metric owner](../src/meshes/toroidal_annulus_mesh.jl),
[cold load and seed](../src/initial_conditions/toroidal_annulus_hybrid.jl),
[coupling owner](../src/coupling/toroidal_annulus_hybrid_coupling.jl),
[field/boundary owner](../src/fields/toroidal_annulus_hybrid_operators.jl),
[coupled lifecycle](../src/models/toroidal_annulus_hybrid_pic.jl),
[diagnostics](../src/diagnostics/toroidal_annulus_hybrid.jl) and
[physical restore](../src/io/toroidal_annulus_hybrid_checkpoint.jl).
