# Compatible periodic Yee implicit electromagnetic PIC

`ImplicitYeeEMPICStepper` is a distinct source implementation of synchronous
nonrelativistic CPU1D3V midpoint PIC. Its orbit-averaged current and electric
gather consume the same original rational cloud coefficients. Numerical,
compiler and restart regressions registered at main `9b4d363` passed the
[26-job CI run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41813949753).
The independent transverse references subsequently passed on minimum/latest
Julia at `62f985f` in [run 42552199520](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42552199520),
with 366 assertions per job, including 38 new transverse checks. The source implementation is
not a floating-point conservation or nonlinear stability theorem.

## Supported model and use

```julia
g = YeeGrid((32,), (1.0,); boundaries = (:periodic,))
state = ImplicitYeeEMPICStepper(g, species;
    shape = YeeTriangularCloud(physical_half_width), c = light_speed)
# Set finite state.fields.E/B and synchronous particles, then:
init!(state, species)
step!(state, species, dt)
diag = implicit_yee_em_diagnostics(state, species)
save_checkpoint("state.chk", state, species)
load_checkpoint!(state, species, "state.chk")
archive_run("state.run", state, species; rng_seed = provenance_seed)
load_archive!(state, species, "state.run")
```

The first family requires Float32/64 nearest rounding with gradual underflow,
ordinary Vector-backed ordered mobile populations, finite signed q, positive m,
nonnegative weights and exact physical neutrality. Empty populations and typed
empty collections support vacuum. Normalization is epsilon0=1, mu0=1/c².
Both transverse polarizations evolve. Bx must be exactly constant and may be
nonzero; its stored entries, including signed zeros, copy unchanged.

The physical period must satisfy `Q(L)==n*Q(dx)` exactly, and the original cloud
half-width must satisfy `0<a<=L/2`. These are explicit first-tier restrictions.
The model periodizes the existing free-space triangular density with three
images. It does not change closed PEC shape/coupling APIs or their plate tails.
Even grids are supported: the staggered derivative has only a DC nullspace.
There is no spectral Nyquist suppression or current projection.

Population identity/order and original q/m bind at construction; init/restore
bind counts, weights and provenance. External source changes require a
consistent initialized field/source boundary; no injection-current or external
work owner is supplied. One instance/population/grid permits one nonreentrant
caller. Custom setter failures, OS protection changes and concurrent mutation
are outside the ordinary-storage publication contract.

All physical x/v/E/B are synchronized at integer time. Original zero dt is a
validated no-op; positive dt uses the canonical resolved floor, must survive
conversion and must advance the stored clock. Different later positive steps
need no leapfrog recentering, without implying arbitrary variable-step accuracy
or stability. Actual current is interval data with a nominal half-step offset.
No timestamp is reconstructed by subtracting dt from rounded time.

## Actual coupled equations

Every nonlinear trial rebuilds exact endpoint-average velocity, lifted/wrapped
endpoint, periodized charge and orbit-averaged current/gather from the same
velocity/E/By/Bz unknown. The trial's ideal and represented transport must each
be shorter than one period. Winding is retained across the seam; this is not a
minimum-image trajectory. A zero represented path uses the analytic stationary
basis, while any unresolved displacement remains in the recorded path residual.
No curved-orbit integration, subcycling or physical hidden-event theorem is
supplied.

The longitudinal electric basis is the path average of `dx*g_e`; transverse
components use the nodal cloud average. Those identical coefficients deposit
current with original charge/weight and endpoint-average velocity. The complete
signed current rounds once per public entry. Charge snapshots are exact
`Rational{BigInt}` extensive nodal charges. Required final current must fit grid
precision; a more general cancellation with a field update is outside this
stored-current contract. Required normalized residuals must also round finite;
a different algebraic solver reformulation could admit a wider range.

The full Lorentz midpoint residual reuses canonical original-Q Lorentz increment
parts; this is not an Euler advance. Mean fields, factors and canonical forward/
backward Yee curl expressions remain rational through residual completion. Bx
has no1D Maxwell update. No separately rounded gather/mean/curl or frozen field
predictor supplies the force, and no post-solve correction projects fields.

Inside every residual the published endpoint Ex is rebuilt from the exact
endpoint charge and the trial's solved mean by the same cyclic prefix used at
initialization. Its actual midpoint follows from the old and published values.
The electric Newton rows still use the unknown field. The report's
`exact_residuals` gate those solver rows; the stored per-node
`exact_electric_residual` and report `field_residual_work` use the published field.
Reconstruction separately enforces endpoint Gauss up to Ex-storage rounding.

The published longitudinal Ampere residual need not fit the nonlinear tolerance.
In exact stored-operand arithmetic it is `r = P*a + S*(C-G0) + eta`: `a` is the
unknown-field Ampere row, `P` its constant/DC projection, `S` the zero-mean cyclic
primitive, `C` the continuity defect, `G0` the old Gauss defect, and `eta` the final
Ex-storage error. After fixed electric normalization, solver-row tolerance controls only the first
term; the separately admitted source and storage terms also contribute. The final store
need not preserve the intended Ex mean exactly. Physical constraint gates and
actual residual-work receipts remain separate; no stronger published-row norm
or exact physical-equation guarantee is implied by solver convergence.

The shared bounded Newton/forward-Jv/restarted-GMRES/Armijo owner retains its
actual action and fresh-final gates. Defaults are nonlinear_rtol=sqrt(eps(T)),
nonlinear_atol=64eps(T), max_newton=20, max_linear=200, memory=20,
max_backtrack=12 and max_jv_shrink=8. Three scales stay fixed during the solve.
The family also gates the original exact solver-row norms; a residual rounded
to zero cannot certify an out-of-budget equation. Only owned trial-domain errors
allow existing bounded shrink; other errors rethrow. Exhaustion rejects the
whole operation, without a smaller-step retry.

## Charge, work and diagnostics

Initialization deposits exact neutral charge and initializes Ex by a cyclic
prefix while retaining its seeded mean. It preserves transverse seeds and
rejects inadmissible Bx. There are no free plates or imposed longitudinal voltage.
For the actual represented path, continuity error is the sum of the explicit
coordinate-rounding contribution and final-current rounding. Actual continuity
and Gauss are separately gated; defaults are constraint_rtol=256eps(T),
constraint_atol=0. The derived bound does not relax this physical constraint gate.

`CompatibleYeeSolveReport` retains exact residuals, constraints, path/current
bounds, counts and signed work. Its four work terms sum to the complete physical
endpoint energy change:

- DeltaK minus particle electric work, equal to actual Lorentz residual work.
- Particle minus unrounded grid-current work, computed from the paired owners.
- Unrounded minus stored grid-current work, exposing current rounding.
- DeltaU plus stored grid-current work, exposing field residual work and the
  computed periodic curl-adjoint contribution.

No term is filled with a fake zero. Exact midpoint/compatible-basis algebra has
an energy/charge law when its equations and kinematics are solved exactly;
actual stored residuals remain visible. This family removes the spectral
current-correction work mismatch of `ImplicitEMPICStepper`; it does not relabel
that earlier method or inherit its payload meaning. The design follows the
same-current/electric-shape requirement in
[Angus et al. (2023), equations9–11/18–24](https://debog.github.io/Files/2023_Angus_EtAl_JCP.pdf),
with a separately specified physical-width cloud rather than that paper's CC1.

Diagnostics expose physical K/U/H and particle momentum, exact source/Gauss and
the detached last report. Complete signed energies use original stored operands;
true displayed overflow is Inf. There is no exact momentum theorem, universal
large-step stability, time-accuracy or measured efficiency claim.

## Publication and persistence

Every operation validates storage/source/boundary, resets private scratch,
solves/rebuilds/gates all physical data and prepares report, replacements and
capacity before one public commit. Expected rejection preserves physical fields,
particles, source snapshots/bindings, clock and report. Field/current/charge
array and ParticleSet identities remain; particle leaves may be replaced.

Distinct schema1 `:implicit_periodic_yee_cloud_midpoint` stores all6 fields,
three interval currents, exact old/new charges, ordered particles/provenance,
original width TYPE/value, fixed geometry/physics/solver policies, actual
synchronous clock/history and report. Count0 requires zero current/equal sources,
zero history and no report. Advanced history requires `0<last_dt<=time`, not
`time==step*last_dt`. There is no RNG; archive seed is provenance. Raw checkpoints
are atomic without a separate checksum; run-container4 verifies archive checksum.

Restore validates original links before whole-graph detachment and revalidates
all incoming source/field/phase/report data. It prepares resized replacements,
bindings and solver capacity before publication. Old destination physical,
count/qm/source/history/report values may be stale; storage, owners and fixed
configuration still match. Fields/current/exact sources copy verbatim, with no
public redeposition, projection, priming or nonlinear solve. Past report checks
are not a replay proof of old particle trajectories. Other implicit/explicit/
closed-wall schemas reject without silent migration.

D2/3 compatible paths/forms, broader spacing/width/transport, bounded walls,
collisions, relativity, schedules and distributed/device paths remain open.

## Acceptance scope

The registered `test/test_implicit_yee_em.jl` regressions cover an analytic
longitudinal DC response, nonuniform published Ex on odd/even Float32/64 grids,
constraint/work receipts, persistence, failure preservation, constructor guards
and diagnostics. The green run above does not establish transverse physics:
those cases seed longitudinal electric fields and velocities.

`test/test_implicit_yee_transverse.jl` adds two independent Float32/64 references:
a four-cell vacuum Maxwell Cayley map assembled from cyclic Yee incidence,
including both polarizations and the true Yee Nyquist mode; and a magnetized,
exact-neutral two-population cold-plasma DC midpoint map with unequal original
q/m and stationary x. Their exact rational matrices use no production
residual/curl/current helper. They compare every physical endpoint, preserve
constant Bx and absent channels, and check the actual constraints/work receipts.
Both references are registered in full test group1 and the `implicit-yee-em`
focused GitHub selector. The recorded two-job run above passed both references
on Julia 1.12.7 and 1.13.3; every raw job log was checked against that revision.
Full exact-head CI remains required before later integrated main publication.

Actual order/dispersion, nonuniform transverse transport, nonlinear conditioning
and complete-step resource acceptance remain open. Compiler/format and all
registered regression acceptance must be repeated on each integrated head;
past green CI is evidence only for its recorded revision and scope.
