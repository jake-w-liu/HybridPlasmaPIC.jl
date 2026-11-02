# Coupled periodic shear hybrid

`ShearedPeriodicHybridPIC` evolves actual full-f kinetic ions and their common
magnetic/electric state on a nonorthogonal periodic mesh. This is the first
bounded curvilinear family. It does not complete CAP-006A's solid-torus axis,
wall, toroidal equilibrium, variable-J, imported-map, AMR, MPI or device work.
Runtime, refinement, independent-solver and resource acceptance remain separate
from source implementation.

```julia
logical = FourierGrid((5, 5, 5), (4.0, 4.0, 4.0))
mesh = PeriodicShearGrid(logical; amplitude=0.15)
ions = ParticleSet{3,Float64}(5^3; q=1, m=1)
load_sheared_uniform!(ions, mesh; number_density=1, velocity=(0,0,0))
species = [ions]
simulation = ShearedPeriodicHybridPIC(mesh, species; Te=0.1, B0=(0,0,1))
init!(simulation, species)
step!(simulation, species, 0.001)
report = sheared_hybrid_diagnostics(simulation)
save_checkpoint("shear.bin", simulation, species)
load_checkpoint!(simulation, species, "shear.bin")
archive_run("shear-run.bin", simulation, species; rng_seed=0)
load_archive!(simulation, species, "shear-run.bin")
```

This is an API illustration, not a validated equilibrium or numerical result.
The loader prepares a real deterministic cold lattice and number weights; it
never assigns an expected deposited density or measured equilibrium residual.
Use physical sampling and convergence studies to choose timestep, marker count
and mesh resolution. The Cartesian whistler recommendation does not certify
this mapped coupled RK4 method.

## Coordinates, components and units

For logical xi on the periodic box, the fixed map is
`x=(xi1+A*sin(2pi*xi2/L2),xi2,xi3)`. Its inverse subtracts the same sinusoidal
shift. `to_cartesian(mesh, logical_tuple)` and `to_logical(mesh, physical_tuple)`
return half-open periodic representatives using the canonical periodic owner.
Trig uses bounded fractions; map and metric rounding remain numerical geometry
error, without an exact floating-point forward/inverse identity.

For `s=(2pi*A/L2)*cos(2pi*xi2/L2)`, the map derivative is
`F=[1 s 0; 0 1 0; 0 0 1]`, with J=1. The metric slope is cached along axis2.
Map phases use sincospi; the slope coefficient uses the recorded rounded tau.
Both retain original factors through range completion. Only CPU Float32/Float64, nearest rounding, gradual underflow,
and odd sizes at least three are admitted. A=0 is the Cartesian limit. Nonzero
amplitude must survive conversion. The family supports this fixed map, not an
arbitrary mapping callback or a variable-J metric.

A bound `ParticleSet.x` contains **logical positions**. Its three `v` components
are **physical Cartesian velocities**. Positive q/m, nonnegative physical number
weights, counts, order, owner and leaf identities are fixed. Empty individual
populations and zero-weight tracers are allowed; the collection must deposit
positive actual density at every grid point. Call the cold loader before binding.
With M particles per axis per logical cell, it requires Ng*M^3 particles and
assigns each number weight n0*dV/M^3. IDs/tags are preserved by loading and stepping.

Units are the package's coherent Omega_ci hybrid normalization, with elementary
charge and mu0 equal to one. Te is normalized electron energy; q/m remain the
actual normalized species values. SI conversion requires consistent external
reference scales. There is no automatic Gaussian-c or independent e/mu0 change.
The fixed model is ideal, collisionless, nonrelativistic and massless-electron,
with isothermal pe=Te*n and eta=etaH=de2=0. There is no density floor, vacuum fill,
kinetic electron, source/controller or active RNG state.

The primary `simulation.beta` is contravariant magnetic **flux density** F^-1 B.
The public `fields` tuple holds n, Ji, ui, B, E, J, pe, gradp and ninv. Vector
entries are physical Cartesian components. n is the electron density, Ji the
charge-weighted ion current and ui=Ji/n. Edit beta for initial magnetic data;
B is a derived snapshot. Positions, velocities or beta edited between calls
require init!; fixed geometry/Te/q/m/weights/counts require a matching new object.
Ready fields are owned by the family between calls.

Quadrature uses the exact product of stored logical spacings. Ng*dV can differ
from the nominal product of lengths by rounding; diagnostics expose that defect.
Oriented area covectors per logical face area are (1,-s,0), (0,1,0), (0,0,1).
The fields are collocated flux densities, with no independent staggered face
values or exact local finite-volume Gauss theorem. Cross-section flux diagnostics
use the recorded logical quadrature and these component conventions.

## Coupling, field calculation and time phase

Canonical TSC values and the existing complete contribution iterator deposit
n and all Cartesian Ji components across the whole species collection before
one materialization per grid entry. The existing complete gather interpolates
Cartesian E/B with those same stored per-axis weights. Thus the exact arithmetic
pair has the particle/grid ion-work identity; materialization introduces rounding.
Logical field components are never silently gathered as Cartesian vectors.

For primary beta, the field calculation is

```
Bphys = F*beta
Bcov  = transpose(F)*Bphys
Jphys = F*curl_logical(Bcov)
gradp = (D1pe, D2pe-s*D1pe, D3pe)
Ephys = -ui cross Bphys + (Jphys cross Bphys-gradp)/n
Ecov  = transpose(F)*Ephys
beta_dot = -curl_logical(Ecov)
```

The canonical logical Fourier operators and canonical pointwise Ohm rows supply
these operations. The periodic derivative algebra makes div(curl)=0 and gives
the magnetic curl work identity in exact arithmetic. The implementation reports
actual defects; it does not transfer those identities into an exact FFT/RK4
trajectory conservation claim.

Every prepared initialization, RK stage, endpoint and incoming restart rejects
the canonical Ohm `source_risk` flag. A finite completed local E does not settle
its differentiated source; original-product mapped Faraday recovery is outside
this family. There is no rounded-E fallback. Actual density, reciprocal,
pressure, transformed fields and rates must materialize within their stated
range. Te>0 with positive pressure disappearing to zero rejects. Fresh magnetic
Fourier data pass a normalized modal component max divergence gate before
publication, with density-independent B/length tolerance units. Comparisons
retain exact scales before narrowing receipt scalars. No projection or filter
repairs failed divergence admission.

The coupled characteristics are `xidot=(vx-s(xi2)*vy,vy,vz)` and
`vdot=(q/m)*(E+v cross B)`. The Lorentz evaluator retains original q/m and all
E/v/B factors and materializes the complete continuous derivative once; it
contains no timestep or finite-step Boris/Euler map. Nonfinite rates or nonzero
exact rates disappearing on materialization reject. Representable subnormal
rates are admitted.

All four classical RK4 stages re-deposit every population and recalculate the
common fields. A fifth endpoint refresh synchronizes beta, particles and every
public derived field. Original-duration RK affine arithmetic is reused. Velocities
and positions share the integer time, without leapfrog priming or a hidden half
interval. Actual zero dt is a validated no-op; disappearing nonzero dt, negative
or nonfinite dt and a nonadvancing stored clock reject.

Every stage and endpoint is private until all gates succeed. Rejection leaves
all public physical leaves, clock and preceding receipt unchanged. Publication
uses fixed-length ordinary stores only. Work and FFT plans are exclusive during
a call. The busy guard prevents reentry; it does not synchronize external mutation.

## Diagnostics and complete restart

Diagnostics return exact number/charge inventories, Cartesian momentum,
kinetic and physical magnetic energy, canonical isothermal free energy
Te*dV*sum(n*log(n)) and their complete total before reservoir narrowing. The
reference density for that logarithm is the normalized unit density. Signed
particle/grid ion-work and magnetic curl-work defects, actual density minimum,
shear, divergence receipts, volume/deposition defects and cross-section fluxes
remain visible. Saturated displayed outputs carry a range flag. Fresh work
observations use private prepared fields and expose their own receipt.

No empirical marker entropy, exact pressure/free-energy balance, exact total
energy/momentum theorem, toroidal angular momentum or unmeasured order/cost is
claimed. Ordinary shape deposition has no exact Fourier-grid continuity law;
this explicit method is distinct from conservative implicit midpoint schemes.
The model/component distinction follows
[Stanier and Chacon, sections 2 and 3.2–3.5](https://arxiv.org/html/2110.03886v1),
whose implicit conservation results do not transfer to this family.

Distinct schema1 `:sheared_periodic_hybrid` includes the actual metric cache and
all fixed coordinate/component/source-risk/normalization/quadrature policies,
all ordered ion leaves and provenance, primary beta, every exposed field, clock,
completed interval and scalar receipt. Exact fixed configuration and inventory
must match the destination. Original aliases are checked before whole-snapshot
detachment and revalidation. Saved density/current are privately re-deposited;
other stored fields pass their relative refresh gates and saved beta its fresh
divergence/source-risk gates. This validation has a cost.

Stored physical fields still restore verbatim, without init!, projection,
recentering or replacement solve. Stale old evolving destination fields, x/v,
provenance and clock can be replaced. Failure preserves the complete public
phase. `load_archive!` first verifies the canonical run schema/checksum, then
uses the same physical restore; `load_archive(path)` remains a read-only reader.
Owners-default checkpoint/archive restore overloads are available. Archive
rng_seed is provenance for this deterministic family, not omitted RNG history.
