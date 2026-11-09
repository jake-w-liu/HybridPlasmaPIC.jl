# Adaptive Cartesian composite electrostatic PIC

This CPU family implements coupled 2D3V and 3D3V electrostatic particles and a
multilevel Cartesian Q1 field. Source review, remote essential runtime checks,
scientific convergence/self-force evidence and measured cost are separate
acceptance stages. No runtime result is asserted by this new guide. CAP006C
still requires periodic/neutral-Neumann/material boundary stages, hybrid/EM
synchronization, subcycling and distributed/device evolution.

## Construct and evolve

`CartesianAMRHierarchy(nbase,L,refinements=())` owns a dyadic forest. Each
refinement `(level,index_tuple)` replaces that logical cell by 2 cells per
axis; base level is0, indices are one-based on `nbase .* 2^level`. Every ancestor
must be refined. Input order is canonicalized; duplicates and invalid ancestors
reject. Valid leaves touching at a face, edge or corner must differ by at most
one level. Original finite positive Float32/64 lengths define exact rational
geometry. All counts/indices must fit Int. No finest-grid array is allocated.

`AdaptiveCompositeElectrostaticPIC(mesh,species;clouds,background_charge=0,
rtol=sqrt(eps(T)),atol=0,maxiter=0)` binds an ordered collection of ordinary
Vector-backed `ParticleSet{D,T}`. Each axis has an existing
`YeeTriangularCloud(a)`; original width type/bits are retained independently,
including mixed Float32/64 widths. Hold these physical widths fixed when
refining. q may be signed/zero, m finite positive and number weights finite
nonnegative. Empty rows and zero-weight charged tracers work. All3 velocities
are retained; D2's third component is unforced. Units are epsilon0=1 and, in
D2, per unit omitted length. Original q/m, row count, cloud/background bits and
solver policies are constructor-bound. No hidden neutralizing source is added.

The following example is an unexecuted API illustration:

```julia
using HybridPlasmaPIC
T = Float64
mesh = CartesianAMRHierarchy((3,3),(T(6),T(6)))
p = ParticleSet{2,T}(1)
p.x[1][1],p.x[2][1] = T(2.75),T(3.25)
p.v[1][1],p.v[2][1] = T(0.125),T(0.125)
species = [p,ParticleSet{2,T}(0)]
st = AdaptiveCompositeElectrostaticPIC(mesh,species;
    clouds=(YeeTriangularCloud(T(0.5)),YeeTriangularCloud(T(0.5))))
init!(st,species)
step!(st,species,T(0.01))
adapt!(st,species;refine_threshold=T(0.02),coarsen_threshold=T(0.005),
    max_level=2,buffer=1)
budget = composite_electrostatic_diagnostics(st,species)
```

Choose duration for the actual physical problem and verify convergence. This
family supplies no universal nonlinear stability certificate and does not copy
the unrelated D1 rest-mode formula. Positive dt is fixed after the first success;
converted nonzero dt disappearing to zero rejects. True signed-zero steps
validate ready state and advance nothing. Global time adds the actual converted
h and must advance strictly; it is not required to equal step*h. Exact Q(h)/2
retains subnormal half durations. Nearest rounding/gradual underflow are required.

`init!` explicitly starts a new run and can rebind matching-factor rows with
changed counts/weights/provenance. Ordinary steps preserve particle objects,
component identities, counts, q/m, weights/IDs/tags. External value/source/field
mutation is checked against independent committed state and requires init!.
Successful mutation can replace `state.mesh` and `state.phi`; refetch them.
There is one primary free-master potential, with a derived full physical field.

## Grounded boundary and complete fixed cloud

All exterior faces are grounded: phi=0, including edges/corners. Particle
centers remain strictly inside the box. Their physical cloud can cross a wall
or span multiple patches: its outside part is bound to the closest face/edge/
corner using the canonical closed-cloud endpoint-tail policy. No cloud charge
is dropped or renormalized. This is a declared soft-cloud boundary model;
material contact, absorption, reflection, injection and moving voltages are
separate unimplemented stages. Each actual KDK drift's exact straight endpoint
and represented endpoint must remain strict-interior. The convex box then
excludes wall contact on that stage. Crossing refinement interfaces is allowed.

For interval [l,u], h=u-l, reuse the canonical triangular primitive P and A:
`m_left=P(x-l)-A`, `m_right=A-P(x-u)`, `m0=P(x-l)-P(x-u)`.
Add `1-P(x)` to the left moment only at the physical left wall and
`P(x-L)` to the right moment only at the physical right wall. Tensor products,
condensed through the same hanging constraints as the field, deposit the full
original q*w. Background integrates the actual unique leaf volume only.
The differentiated axis uses interior `m0/h`; other axes retain tail-inclusive
moments. This derivative gathers the actual physical-potential gradient and
satisfies the source/force virtual-work identity, including near-wall clouds.

Diagnostics separate full particle charge, inside charge, every tensor wall
stratum, background, numerical free-conductor reaction and net wall charge.
`stratum_axes` labels `inventory.strata` with -1/0/+1 meaning left wall/inside/
right wall per axis; the first axis varies fastest. `boundary_faces` and
`interior_faces` label the flux/jump vectors with indices into returned mesh.faces.
The grounded controller pairing `sum(phi_boundary*Delta reaction)` is zero
because its specified voltage is zero. Particle work and mesh exchange are
measured separately and never reset under that boundary condition.

## Composite field and residual meaning

Every unique leaf has a Q1 multilinear basis. Exact hanging rows P express raw
vertices through independent masters, are nonnegative and sum to one. Every
incident coarse trace and every interior face-fragment corner is checked for
exact equality, including D3 edge/corner contacts. Covered parents add no source,
energy or independent field. This serial global solve needs no ghost cells.

The global operator is `K=P' sum(K_leaf) P`, with exact local tensor mass and
gradient integrals. Zero boundary masters give a positive-definite free block A.
A cold call-local Krylov CG solve uses one exact global diagonal scale and
Jacobi preconditioning. Its stopping threshold is derived from the requested
physical tolerance; solver success alone is insufficient. A fresh original-Q
residual `A*phi-rhs` must satisfy `maxabs(r)<=atol+rtol*maxabs(rhs)` in physical
nodal-charge units. No fallback/direct solver, source adjustment, shift or
silent tolerance relaxation occurs. Nonconvergence or required nonfinite
materialization rejects before publication. Extreme-range finite mathematical
problems may therefore be inadmissible. Solver scratch is not restart state.

The primary phi values are T, but hanging values are exact P*Q(phi), without
independent rounding. Force uses the full piecewise gradient, whose transverse
variation matters. `composite_electric_field(st,x::NTuple{D,T})` returns an exact-Q
3-tuple at a closed-box point; outer faces use the one-sided interior trace and
internal ties use the right leaf. `patch_electric_field(st,patch)` returns a
fresh 3-tuple of T arrays of volume means, including covered cells. Those means
are diagnostic/tagging views and are not the particle force.

The residual is weak Galerkin Gauss. Q1 normal E can jump across interfaces;
physical face-flux integrals, normal jumps and boundary flux/reaction discrepancy
are reported separately. No strong finite-volume Gauss, absent interface
self-force or EM-compatible conservative current is claimed. Field energy is
`phi' A phi/2`, the actual gradient integral, not energy of leaf mean E.

## Steps, adaptation and actual work

A physical step performs source-compatible half-kick, full D-coordinate drift,
one composite endpoint solve and final half-kick. The original q/m and exact
full-duration half use the canonical3V Boris owner with B=0. CSR ownership and
grouped straight-segment face/edge/corner crossings update without particle
splitting or sorting. Endpoint ties match bins; tangencies do not invent events.

For unchanged mesh, with actual endpoint sources b and potentials u:
`W_rho=-(u1+u0)'*(b1-b0)/2`. `W_F` is actual half-kick force power using
pre/post mean velocities for both half-kicks. Signed receipts retain
`D_field=DeltaU+W_rho`, `D_kick=DeltaK-W_F`,
`D_trajectory=W_F-W_rho`, and their sum `D_total=DeltaU+DeltaK`.
The field term also equals `(u1+u0)'*(r1-r0)/2`. Actual rounding and trajectory
mismatch remain visible; endpoint nodal work is not an implemented face current.

`regrid!` accepts a canonical balanced mesh of the same original base/domain.
It preserves every particle value/identity, deposits the same cloud on the new
basis and reconstructs the new field. Actual `U_new-U_old` is signed mesh
energy exchange, not a claim of energy-conserving remapping. Equal validated
topology is a true no-op with no reconstruction, clock/epoch/history/counter or
ledger reset. There is no unused E transfer presented as evolved state.

`adapt!` requires finite T refine/coarsen thresholds, explicit max_level and
nonnegative Int-sized buffer. Coarsen<refine, except both may be0. It averages
the actual old piecewise field on virtual same-level neighbor boxes, then uses
the existing undivided-curvature predicate in each direction/component. It
coarsens complete low families bottom-up, refines/buffers high cells and closes
2:1 balance across faces/edges/corners before one reconstruction. The criterion
is a multilevel heuristic, not a physical error estimator. max_level cannot be
below the existing deepest leaf. No hidden callback adapts partway through KDK.

Global time/count and fixed h continue across regrids; epoch origins are stored
verbatim and epoch step count is derived. The preceding physical step remains
on its actual old mesh. A complete detached ending payload is retained only
while the new epoch has no physical step; repeated no-step regrids preserve it.
The next step replaces it. Never interpret old history on current leaves.
Exact cumulative accounting is `H-H0=physical_defect+mesh_work`, with
`physical_defect=field_defect+kick_defect+trajectory_defect`. Init starts a new
origin; zero-step/save/load/regrid do not reset physical accumulators.

## Failure and persistence

One exclusive guard owns each operation. All candidates, constraints, charge,
solve, path/range, report, history and alias checks precede ordinary public
stores. Failure preserves public state and permits a subsequent valid call.
Ordinary storage with reliable Base.mightalias is required; unsafe pointer
aliasing and concurrent external mutation are unsupported. Regrid/step are
separate transactions, so a failed regrid preserves a preceding accepted step.

Distinct schema1 `:adaptive_composite_es` raw/archive records save immutable
source/solver policies, canonical current topology, all particles and phi,
solve/step receipts, global/epoch phase, cumulative ledgers, last regrid and the
complete old-mesh payload when required. Derived geometry/constraints/operator,
charge/means/bins rebuild from saved primary state. No solver workspace or
invented half-step state is serialized. Restore checks original input aliases
before copying, admits changed topology/counts under matching fixed policies,
and ignores stale destination numerical values. Saved primary fields/history/
phase restore verbatim, without solving/reinitializing. Invalid keys/schema,
config, source/residual, history or accounting reject atomically. Raw uses the
shared atomic serializer; archive uses the shared checksum container. Same-
runtime continuation requires remote acceptance; cross-Julia serialization
bitwise portability is not established.

Periodic seam constraints and periodized clouds, neutral-Neumann face charges/
gauge/full-row residual, material/electrode events/work, hybrid/EM compatible
current/flux/divergence synchronization, local subcycling, partition/rebalance
and device storage remain concrete later implementation stages. Scientific
fixed-width convergence, h convergence, wall-image/interface self-force,
physical flux convergence, long-run repeated adaptation and independent-code
comparison remain open. The implementation rebuilds exact derived caches for
admission and allocates private rational candidates and retained history.
Storage includes particle/mesh state, actual crossing-event receipts and the
necessary prior ending payload; no path-independent linear storage or speed
claim is made. Tree/ray queries avoid unconditional particle-by-all-leaves and all-leaf-pairs scans;
physical wide-cloud support may itself touch many leaves. Measure cost and peak
storage at matched physical error, including regrids, before claiming advantage.
