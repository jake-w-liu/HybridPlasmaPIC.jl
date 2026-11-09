# Periodic adaptive Cartesian electrostatic PIC

The periodic boundary mode shares the Cartesian hierarchy, original physical
clouds, Q1 field, KDK, adaptation, accounting and restart owners with
[the grounded family](ADAPTIVE_COMPOSITE_ELECTROSTATIC.md). It is a coupled CPU
2D3V/3D3V electrostatic model on an all-axis torus. Source review, remote runtime
acceptance, scientific convergence and resource acceptance are separate. This
guide asserts no executed runtime result.

## Use and physical source

Construct `CartesianAMRHierarchy(nbase,L,refinements=();periodic=true)`, then
`AdaptiveCompositeElectrostaticPIC(mesh,species;clouds,background_charge=0,
periodic_neutrality=:require,rtol,atol,maxiter)`. Existing numerical-control
defaults remain available. Refine addresses, original Float32/64 lengths,
ordinary particle storage and fixed original q/m/cloud-width policies retain
the Cartesian family's meaning. Positions are half-open `[0,L)` in every axis;
a seam is part of the domain, not a reflecting or absorbing wall.

The following example is an unexecuted API illustration:

```julia
using HybridPlasmaPIC
T = Float64
mesh = CartesianAMRHierarchy((2,2),(T(4),T(4));periodic=true)
a = ParticleSet{2,T}(1;q=T(1),m=T(2))
b = ParticleSet{2,T}(1;q=T(-0.5),m=T(3))
a.x[1][1],a.x[2][1],a.weight[1] = T(0.25),T(0.375),T(1.25)
b.x[1][1],b.x[2][1],b.weight[1] = T(2.75),T(3.25),T(2.5)
a.v[1][1],b.v[1][1] = T(0.125),T(-0.0625)
species = [a,b,ParticleSet{2,T}(0)]
st = AdaptiveCompositeElectrostaticPIC(mesh,species;
    clouds=(YeeTriangularCloud(Float64(0.625)),
            YeeTriangularCloud(Float32(0.5))))
init!(st,species)
step!(st,species,T(1/32))
adapt!(st,species;refine_threshold=T(0.02),
    coarsen_threshold=T(0.005),max_level=2,buffer=1)
diag = composite_electrostatic_diagnostics(st,species)
```

`:require` admits only exact original particle plus declared uniform-background
charge zero. It does not tolerate a nonzero null mode and then remove it.
`:neutralize` explicitly adds uniform immobile charge `-declared_net_charge`,
with exact rational density, to the declared background. That source is deposited
through the same actual leaf volumes. `inventory.declared_background`,
`inventory.neutralizing_charge`, `inventory.background` and `background_density`
separate the declared, added and actual combined source. Weights are unchanged.
Successful `init!` starts a new run; sources remain fixed between initializations.
A nondefault periodic-neutrality policy rejects on a grounded mesh.

The field is `E=-grad(phi)` for a periodic potential. Its harmonic circulation
and original-Q volume mean are zero. There is no prescribed DC electric field,
conductor reaction, wall tail or electrode controller in this mode. A physical
harmonic forcing requires its own source/work model. Units are epsilon0=1 and,
in D2, per unit omitted length; all 3 velocities contribute kinetic energy and
D2's third component is unforced. Empty rows and zero-weight charged tracers
retain their ordinary meanings.

## Seams, fixed clouds and the composite field

The same owned dyadic forest has unique valid leaves and no covered-cell double
counting. Face/edge/corner 2:1 balance includes opposite periodic images. Raw
seam vertices are identified before hanging-node elimination, so a fine cut
vertex can be constrained by its opposite coarse trace. Exact nonnegative unit
rows P must match every incident and translated face-fragment trace. No
post-solve or post-hanging master quotient substitutes for that construction.

The triangular cloud keeps each axis's original positive finite width and may
cross several seams or have width larger than the period. The existing scalar
R/P/A owner supplies interior physical interval moments; their periodic image
sums are formed before tensor products. There are no clamped conducting-wall
tails. Integer image bands have polynomial degree at most 3 and are summed
exactly from at most 4 primitive evaluations per band, without a width-proportional
image loop or fitted polynomial. Wide physical support can still touch many
actual leaves. This is not a measured speed or peak-memory claim.

Charge and force use that same constrained periodic basis. Force uses its full
represented-potential gradient, not a cell-mean field. The same local Q1 tensors
assemble `K=P' sum(K_leaf) P`. Its constant nullspace is fixed by setting the
origin master potential to 0 as a gauge; this point is not an electrode.
The projected free block is positive definite. The cold scaled/Jacobi CG owner
is reused, with fresh original-Q residual admission at ALL masters, including
the omitted gauge row. The physical target is
`atol+rtol*maxabs(full_nodal_source)`. No RHS mean subtraction, physical ridge,
ignored gauge residual, hidden direct solve or relaxed tolerance is used.

`full_gauss_residual`, `gauge_master` and `gauge_residual` describe the complete
weak Galerkin constraint. `harmonic_mean_field` is the actual original-Q gradient
mean. Interior and `seam_normal_flux_jumps` are separate physical-gradient trace
observations with face indices. Q1 normal E can jump across a seam or refinement
interface; weak Gauss is not a strong finite-volume Gauss or an EM current law.
Field energy is the actual gradient integral `phi' A phi/2`.
`composite_electric_field(st,x)` wraps finite T query points to the torus;
`patch_electric_field` remains a detached volume-mean view.

## Actual lifted motion and work

Synchronous KDK retains the original-factor3V half-kicks. Periodic drift first
forms `ideal=Q(old_x)+Q(h)*Q(vhalf)`. Canonical scalar wrapping, nearest-T endpoint
rounding and its rounded-to-L correction reconstruct the actual represented
lift `Q(new_wrapped_x)+winding*Q(L)`. Its difference from ideal is recorded.
No unwrapped T overflow is required merely to reach a finite wrapped endpoint.

Finite multi-period transport is allowed subject to checked event-count/range
and ordinary allocation admission. The actual represented lifted segment is
split at tile boundaries and routed through the tree, with simultaneous seam/
edge/corner events grouped and endpoint ownership explicit. Identical wrapped
endpoints do not erase completed periods; this is not a minimum-image path.
`transport` contains independent per-population `old_positions`,
`half_velocity`, BigInt `winding`, exact `defect`, and `tiles_before`/`tiles_after`
aligned with `last_step_report.crossings`. These are the last physical interval's
provenance, not an all-time unwrapped trajectory origin. A rounded-to-L endpoint
folds to 0 with the corresponding winding and visible rounding defect.

The shared signed receipts keep source work `W_rho=-ubar'*(b1-b0)`, actual
half-kick power `W_F`, `D_field=DeltaU+W_rho`, `D_kick=DeltaK-W_F`,
`D_trajectory=W_F-W_rho`, and `D_total=DeltaU+DeltaK`.
`D_field` also equals `ubar'*(r1-r0)` for the full original residual. Seam
crossing creates no fictitious electrode work. Large or underresolved steps can
have large trajectory defects; neither finite stages nor wrap admission prove
nonlinear stability, conservation, low self-force or trajectory accuracy.
The first positive converted h remains fixed across steps and regrids.

## Adaptation, failure and restart

Virtual same-level neighbor boxes wrap at the cut. Multilevel curvature,
hysteresis, complete-family coarsening, unique wrapped buffers and seam-aware
balance closure precede one fresh constrained field reconstruction. Regrid
preserves all particle values/identities and actual source inventory; signed
`U_new-U_old` enters mesh work. An equal validated topology is a true no-op.
The shared exact identity is `H-H0=physical_defect+mesh_work`, with physical
defect equal to the sum of the 3 work defects. These are measured receipts, not
a claim of energy-conserving remapping.

Global time/count/fixed h continue through regrids. A retained preceding physical
payload owns its original periodic mesh and complete transport. Repeated no-step
regrids retain it; the next real step replaces it. It is never reinterpreted on
new leaves. Every fallible source, solve, path, role, phase and report check
precedes ordinary public stores, and failure preserves the previous state.
Concurrency and unsafe pointer aliases are outside the ordinary-storage contract.

Periodic raw tag `:periodic_adaptive_composite_es`, schema 1, is distinct from
the unchanged grounded format. Current/history payloads additionally retain
actual Q `background_density` and complete `transport` (or nothing before a
real local step). Original incoming roles are checked before copying. Restore
supports changed counts/topology and stale old numerical values under matching
fixed configuration; primary phi, phase, ledgers and old-mesh transport restore
verbatim without solving or reinitializing. The shared atomic raw serializer
and checksummed archive container remain the publication owners. Cross-boundary
migration and cross-Julia bitwise serialization portability are not promised.

Mixed periodic/physical axes, neutral-Neumann/material/electrode boundaries,
hybrid/EM compatible currents and flux/divergence synchronization, subcycling,
MPI/device hierarchy, scientific convergence/self-force/long-run comparisons
and measured resource acceptance remain later work. All runtime/formatter checks
for this stage run on GitHub minimum/latest Julia. No local numerical workload
or unmeasured optimization claim accompanies this guide.
