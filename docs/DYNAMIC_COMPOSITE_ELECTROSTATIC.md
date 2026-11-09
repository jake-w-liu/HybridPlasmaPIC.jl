# Dynamic composite electrostatic PIC

This source candidate adds completed-boundary mesh adaptation to the existing
CPU 1D physical-cloud electrostatic solver. It has not been run through Julia,
compiler, formatter or numerical tests locally; exact-revision GitHub acceptance
and scientific/resource validation remain open. It does not complete D2/3,
hybrid/EM, distributed/device or subcycled AMR.

## Construct and evolve

`DynamicCompositeElectrostaticPIC` uses the same original-source Neumann solve,
physical triangular-cloud coupling, original-factor KDK mover, isolated free
plates, finite-range and event-free wall guards as
[the static family](STATIC_COMPOSITE_ELECTROSTATIC.md). There is one physical
field on unique leaves and one potential on shared vertices; covered cells add
neither charge nor energy. All q/m/source types and bits remain constructor-bound.

The following API example is unexecuted; it asserts no numerical result:

```julia
using HybridPlasmaPIC

T = Float64
mesh = StaticAMRHierarchy1D(4, T(8), ((1, 2, 2),))
ions = ParticleSet{1,T}(1; q=1, m=1)
ions.x[1][1], ions.v[1][1] = T(3.5), T(0.1)
species = [ions, ParticleSet{1,T}(0)]
st = DynamicCompositeElectrostaticPIC(mesh, species;
    shape=YeeTriangularCloud(T(0.25)))
init!(st, species)
h = T(0.01)
step!(st, species, h)
adapt!(st, species; threshold=T(0.1), buffer=1)
step!(st, species, h)
regrid!(st, species, StaticAMRHierarchy1D(4, T(8))) # Coarsen at the ready boundary.
budget = composite_electrostatic_diagnostics(st, species)
```

Select h for the actual physical problem and the existing strict admission
bounds. The first positive converted duration fixes all subsequent duration
bits. Choose it for the finest planned hierarchy: refinement can tighten the
canonical next-step rest-mode bound, and regridding does not relax that guard or
silently change the fixed duration. Signed-zero steps validate the ready state
and advance nothing. A nonzero
original duration disappearing during conversion rejects. Nearest rounding and
gradual underflow are required. Species use ordinary Float32/64 `ParticleSet`
vectors; signed q, positive m, nonnegative number weights, empty rows and tracers
retain the static model's meaning and restrictions.

`state.mesh`, `state.E` and `state.phi` are the current physical storage.
`state.time`, `state.step` and `state.last_dt` are plain global scalar fields.
Refetch mesh/E/phi after successful mutation: these references can be replaced.
Physical steps preserve public particle component identities; regrids preserve
every particle value and identity. Restore replaces particle leaves. External
field/particle/source mutation requires explicit reinitialization; ready
validation checks actual public values against the bound private physical state.
`patch_electric_field(state,patch_id)` delegates the same width-restricted
covered/leaf query as the static family and returns detached values.

## Adaptation and field reconstruction

`regrid!(state,species,new_mesh)` accepts any canonical nested balanced 2:1
interval hierarchy with the same nbase, precision and original stored L. It
validates and owns a fresh mesh copy. A matching refinement tuple is a validated
no-op: no solver clone, new-mesh field reconstruction, epoch/counter/history
reset or physical store. Read-only ready validation still checks the existing
source map; it can allocate.

`adapt!(state,species;threshold,buffer=1)` uses the existing width-restricted
base-patch E and `refine_flags` undivided second-difference rule. Threshold must
be finite, nonnegative and in mesh precision; buffer is a nonnegative Int-sized
cell distance. The curvature owner does not tag the two boundary cells; buffer
expansion can include them. Contiguous tagged base runs become disjoint 2x
parent-1 patches, so this automatic policy is always two-level and 2:1 balanced.
No tags means the base mesh. The criterion is a declared heuristic, not a
trajectory error estimate. Deeper nested regrids are available explicitly;
multilevel automatic tagging/hysteresis remain separate work.

No particle splitting, merging, reordering, weight change or momentum kick occurs
on regrid. The same finite physical cloud is deposited onto the new nodal basis
by the canonical exact source owner. Its partition of unity and actual dual
volumes conserve macrocharge and background inventory; the original free plates
stay unchanged. The new E/phi are the unique canonical new-mesh constrained
solution. Old E is not an independently evolved unknown, and no unused E
interpolation is presented as an evolution remap.

The actual field-energy difference is reported as signed **mesh energy
exchange**. This is a discretization change, not physical force work and not a
claim of energy-conserving remapping or absent interface self-force.

Each `step!` and `adapt!`/`regrid!` call is its own atomic operation. In the loop
above, failed adaptation leaves the preceding successful physical step intact
and changes nothing itself. No hidden adaptation callback partially succeeds
inside a physical step.

## Global phase and histories

Global physical time/count and fixed duration continue across regrids. A new
topology begins a local canonical static epoch with local time/count/history
zero, using unchanged synchronous velocities. KDK has no hidden half-time
priming here. Global epoch-start time/count are saved verbatim. Local count plus
epoch-start count equals global count; no float-sum or time=count*h identity is
imposed, since repeated additions can round differently.
An epoch-start count of zero requires stored-zero origin time; a positive origin
count requires positive origin time at least as large as the fixed duration.

The previous physical step belongs to its original mesh. Before a regrid replaces
that epoch, the complete canonical old-mesh physical checkpoint is detached with
its own epoch/global time/count metadata. Repeated regrids without a step retain
that same physical report/current and particles. It is validated only against
its own old mesh and saved particles. The next real step supersedes it with the
current epoch's actual report. A second persistent physical snapshot is kept
only while the current epoch has no step after a preceding global step.

The exact cumulative accounting identity is
`H(current)-initial_energy == physical_defect + mesh_work`.
`physical_defect` sums actual canonical `DeltaK+DeltaU` from completed steps;
`mesh_work` sums actual `DeltaU` from completed regrids. These are measured
arithmetic defects/exchanges, not zeroed conservation claims. Before any global
step, the physical-step defect is zero even if time-zero regrids have occurred.

Diagnostics return detached canonical current-mesh charge/energy/residual data,
current mesh/E/phi, global and local phase, cumulative accounting and the full
last regrid report. The inherited `last_report`/`last_impulse` describe the current
static epoch; `last_physical_step.state` carries the actual preceding physical
step on its own mesh after regridding. At initialization it is absent.
No saved finite history proves an arbitrary past trajectory was actually run;
admission checks the declared present source map, phase and accounting identities.

## Failure and continuation

An exclusive guard rejects overlapping public operations. Every physical step
is cloned and executed by the canonical static solver privately. Regrids prepare
a fresh private canonical state. All configuration, source, alias, range,
history, clock and exact-accounting checks and allocations precede public
ordinary leaf/pointer/scalar stores. Expected rejection retains public physical
state, histories, mesh and clocks. Private candidates may be discarded.

`init!` may rebind owner objects/counts/weights/provenance with the original
ordered q/m configuration; it resets global/epoch/regrid/work history at the
current mesh. Failed init leaves the old usable state. Same-size external
component replacement is rejected until explicit init. This is ordinary host
storage only; aliases must be reported truthfully by the canonical overlap API.

```julia
save_checkpoint("adaptive.chk", st, species)
archive_run("adaptive.run", st, species; rng_seed=0)
load_checkpoint!(st, species, "adaptive.chk")
load_archive!(st, species, "adaptive.run")
```

Distinct schema1 `:dynamic_composite_es` saves complete current topology and
static physical state, global/epoch phase, cumulative energy/work and the optional
detached prior physical mesh/state. Restore accepts changed saved topology/counts
and stale old public values/fields/phase, with current ordinary owner/component
identities and immutable constructor configuration retained. Original payload
links are checked before deepcopy; all replacement graphs are prepared privately
and saved fields/current/report are restored verbatim. Canonical source replay
validates them; it does not replace them with reconstructed values. Raw files use
the shared atomic serializer; archives use the existing checksum container and
the same restore. RNG seed is provenance only. Static/other-family migration and
cross-version/rank/device reproducibility are not promised.

This correctness-first implementation copies O(Np+Nleaf) storage per trial and
retains complete physical history only where required. Geometry validation can
rebuild/sort leaves, alias admission compares array roles pairwise, and exact
rational costs depend on operand sizes. No zero-allocation, speed, resource,
self-force, refinement-order or independent-code acceptance follows from source
review. Those gates and the broader adaptive hierarchy remain open.
