# Static composite electrostatic PIC

Source implementation candidate. Exact-revision GitHub numerical, compiler,
restart and resource acceptance is pending. No local Julia or tests were run
for this implementation. CAP-006C remains open beyond this bounded first family.

`StaticAMRHierarchy1D(nbase,L,refinements=())` owns nested static 2:1 interval
patches. Each refinement is `(parent_patch,first_parent_cell,last_parent_cell)`;
parents precede children, siblings cannot overlap and adjacent valid leaves must
satisfy 2:1 balance. Original finite positive Float32/64 L defines exact rational
geometry on `[0,L]`. Integer indices/counts must fit Int. There is no finest-grid
field allocation or independent per-level solve.

`StaticCompositeElectrostaticPIC(mesh,species;shape=YeeTriangularCloud(a),
background_charge=0,left_plate_charge=0)` is CPU 1D longitudinal, nonrelativistic,
all-mobile, synchronous electrostatic PIC. Species are an ordered vector of
ordinary Vector-backed `ParticleSet{1,T}` at mesh precision; q may have either
sign, m is positive, weights nonnegative, velocities finite and transverse
velocities stored zero. Centers stay strictly inside the physical interval.
Zero-weight particles are tracers: they feel q/m force without depositing charge
or kinetic energy. The positive physical cloud width is retained at its original
Float32/64 precision and held fixed under refinement.

Original ordered q/m type/bits are constructor-bound. Changing either requires a
new model. `init!` may replace owner objects with the same ordered factor/type
configuration and species count, and replace particle counts/weights/provenance.
It prepares the entire new binding before publication; failed init retains the
old usable owner list. Ordinary steps preserve membership, component identities,
source/provenance arrays and counts. External phase-space mutation is unsupported;
same-sized replacement cannot generally be detected. Call `init!` explicitly to
reset the numerical phase and bindings under the fixed configuration.

Only valid leaves store physical E; unique composite vertices store phi. Covered
coarse fields are derived, not extra degrees of freedom. `patch_electric_field`
returns a fresh vector with width-restricted covered values. Shared vertex maps
provide patch potential samples. This global 1D solve requires zero ghost width;
there are no unused halo arrays. Particle bins contain each original particle
index once in its actual leaf, with a unique right-cell choice on internal nodes.
Crossings update those bins without splitting/reordering particles or weights.

For actual leaf width h and the canonical physical cloud cumulative A/derivative
g, nodal particle coefficients are `1-A[1]`, adjacent A differences, and `A[end]`.
They are nonnegative and sum exactly to one, including the canonical closed-cloud
endpoint tails. Background charge uses the actual nodal dual volumes: endpoint
half widths and adjacent half-width sums. Both leaf and dual volumes total L;
covered cells contribute no additional charge or energy. Units are epsilon0=1
and unit transverse area; weights carry macro-particle number, not density.

The fixed left free-plate charge and derived
`right=-left-background*L-sum(q*w)` balance total charge without a hidden volume
neutralizer. They are isolated free charges, not fixed electrode voltages. With
leaf-node incidence R and H=diag(h), one composite operator is
`K=R' inv(H) R`, `E=-inv(H)R phi`, `-R'E=rhs`. Its constant nullspace is fixed by
phi[1]=0, a gauge only. Exact rational prefix elimination and potential integration
solve this single Neumann system. Each E and phi is rounded once to T; a required
nonfinite materialized field or potential rejects before physical publication.
A finite E with unrepresentable required phi is an explicit admitted-range limit.

Force gathers `sum(h*E*g)` from stored E, using the derivative of the same cloud
coefficients. Original q, m, gathered force and converted full-duration half are
retained through the canonical rational Boris map; the new family promises no
ordinary-force/native-bit equivalence to the rounded uniform Yee current map.
Stored E and phi are separately rounded: all-node Gauss and constitutive residuals
are reported exactly, with no claim of exact floating-point Gauss or relative
accuracy in a cancelling value. Positive-step/diagnostic/save ready validation
replays the declared source solve and requires matching rounded physical fields;
manual field edits cannot silently introduce another forcing policy.

`step!` uses half-kick, full drift, one composite endpoint solve, half-kick.
Positive converted T duration, including subnormals, uses exact Q(dt)/2; a nonzero
original duration that converts to zero rejects. The first accepted positive dt
fixes later dt bits. A true signed-zero duration validates storage/source/phase
and returns without source solve, bin rebuilding or numerical scratch writes.
Clock advance must remain finite and strict. At an advanced boundary the last
duration cannot exceed time; no `time=step*dt` identity is imposed. Nearest rounding
and gradual underflow are required. There is no particle/field time subcycling.

A strict zero-field rest-mode bound uses
`dt^2*min(1/hmin,1/a)*sum(q^2*w/m)<4`. It is not a general nonlinear or finite-T
stability theorem. A conservative physical-wall exclusion tube additionally
rejects possible physical contacts; the staged drift is independently required
strict-interior. Coarse/fine interfaces impose no contact barrier. There is no
reflection, absorption, injection, collision or arbitrary forcing in this family.

`composite_electrostatic_diagnostics` returns exact charges, current, unique-leaf
field/kinetic energy, all-node residuals and a copied full last report. Nodal
source includes plates; volume nodal charge and its exact dual-volume density
exclude free plates. Let
`I=sum(q*w*(Anew-Aold))`, `W=sum(h*(Eold+Enew)*I)/2`. The signed report retains
DeltaU, DeltaK, `field_solve_defect=DeltaU+W`,
`particle_work_defect=DeltaK-W`, `total_energy_defect=DeltaU+DeltaK`, actual residual
magnitudes and per-interface crossing counts. None is forced zero. Required
exact histories are finite Q scalars; large physical energy does not overflow a
required T report. Refinement loses translation symmetry: interface self-force
and fixed-width convergence must be measured, not claimed absent. No universal
nonlinear physical-trajectory error certificate follows from the wall tube.

All validation, field/particle/report/bin preparation and range checks precede
ordinary public stores. Failure preserves public particles, fields, histories,
binding and clock; private candidates may change. One call at a time, with
ordinary storage whose overlap Base.mightalias faithfully reports. No unsafe
pointer-overlap, concurrent mutation or custom store semantics are supported.

Complete distinct schema1 `:static_composite_es` checkpoints save matching mesh/
source/factor configuration, all ordered particles, E, phi, exact derived right
plate, clock/fixed duration, last integrated current and full report. Only ready
initialized boundaries are saved. Derived bins and solve storage are rebuilt,
never serialized as physical state. Restore requires the same immutable config
and currently bound owner objects, but ignores stale destination particle values,
factors/counts, fields, phase and report. All incoming Q leaves, arrays, source
map, phase and report identities are checked before prepared leaf replacement;
saved fields/history are restored verbatim. No silent reinitialization/projection.

`save_checkpoint`/`load_checkpoint!` use the shared atomic raw serializer;
`archive_run`/`load_archive!` use the existing checksummed byte container and the
same state restore. `rng_seed` is archive provenance only; the model has no RNG.
There is no migration from another solver or different hierarchy. Cross-version
serialization/bitwise continuation is not established by source review.

Dynamic tagging/regrid/remap, D2/3 face/edge/corner and halo algorithms, periodic
composite solves, MPI/devices, hybrid/EM synchronization and performance acceptance
remain later work. The current correctness-first admission rebuilds/sorts the
static layout to validate its mutable caches, and field/report work allocates
exact-rational candidates; no zero-allocation or speed claim is made. Cloud and
path routing use binary-searched actual intervals, CSR construction avoids a
species-by-leaf empty-vector table, and crossing counts visit only crossed
interfaces. The existing scalar AMR transfer functions and periodic/closed
uniform solvers retain their separate contracts and tests.
