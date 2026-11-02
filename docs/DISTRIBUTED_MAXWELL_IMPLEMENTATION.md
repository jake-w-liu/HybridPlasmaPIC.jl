# Rank-owned periodic vacuum Maxwell

`distributed_maxwell_vacuum(grid; c=1)` is a source implementation of source-free
Maxwell evolution on the existing periodic 3D CPU Float32/Float64 physical
pencils. Import `PencilArrays` and `PencilFFTs` alongside HybridPlasmaPIC to load
the optional extension. It supplies a field-only CAP-003 prerequisite; particle
ownership, deposition/gather halos and coupled hybrid/ES/EM evolution remain open.

The current implementation-first pass adds no tests or runtime execution.
Numerical/convergence, actual-rank failures/liveness, compiler/JET, formatting/CI
and complete-step resource acceptance are deferred to a later validation pass
on the exact integrated revision. Source review does not establish those results.

## Fields and collective lifecycle

Construct a `distributed_fourier_grid(n, L; comm)` first, then pass that grid to
`distributed_maxwell_vacuum`. The model allocates initially zero synchronized
`fields.E` and `fields.B`, each a three-tuple of real rank-owned PencilArrays.
Populate owned cells through `PencilArrays.global_view(model.fields.E[c])` and
the corresponding B view. Model time and step start at zero. `c` converts once
to the grid precision and must be finite, positive and equal on every rank.

Call `step_maxwell!(model, dt)` on every rank in the same order. Precision, grid
dimensions/spacing, converted dt, speed and current clocks must agree. A step
is positive, must advance the finite clock, and must be strictly below
`maxwell_cfl(model)`. Nearest rounding and
gradual underflow are required for construction, stepping and diagnostics.
The grid's metadata remain read-only; the grid/model belong to one caller at a
time, with no concurrent mutation. The caller owns the communicator; this model
neither frees it nor finalizes MPI. Plans and transposes retain the existing
[distributed operator lifecycle](DISTRIBUTED_OPERATOR_IMPLEMENTATION_PLAN.md).

All physical/private fields use the exact grid pencil, precision, axes and
ordinary CPU Array backing. Their storage must be pairwise independent and
disjoint from grid workspaces/metadata. Checks use `Base.mightalias`; callers
must provide storage whose overlap it reports. Finite physical values validate
separately from storage, before any staging write or transform.

## Synchronized spectral map

With `h` and `c` denoting the validated converted scalars, one step is

```text
Bhalf = round_T(Bold - (h/2) * curl(Eold))
Enew  = round_T(Eold + h*c² * curl(Bhalf))
Bnew  = round_T(Bhalf - (h/2) * curl(Enew))
```

Every curl calls the existing collective operator, with shared SpectralOperators
mode arithmetic and actual PencilFFTs transposes. Every active local affine row
uses the existing complete Maxwell expression, retaining exact converted-scalar
coefficients through one final rounding. Half durations and `h*c²` are not first
rounded/materialized in the field precision. A stored-zero curl copies the old
entry, including its signed zero. This is a distinct synchronized spectral map,
not a claim of ordinary EMPIC or Yee field-bit identity.

The derivative source is the represented, finite materialized curl. A curl
outside that operator's supported range rejects even if a later timestep could
make the final affine result finite. Original-field derivative underflow,
overflow recovery, general FFT conditioning and cancellation-relative accuracy
are not closed by the scalar completion. No Gauss projection removes an initial
divergence residual. The usual compatible-curl identity is subject to the actual
FFT/operator roundoff; no whole-trajectory conservation/stability bound is added.

The strict spectral KDK constant-step bound is
`2/(c*sqrt(sum(gain[d]^2)))`, where each gain is the existing axis maximum of
the derivative-resolved, odd-Nyquist-zeroed wavenumbers. The bound rounds down
through the canonical exact square-root owner. Equality rejects the repeated
amplification-eigenvalue endpoint. A completely derivative-null grid returns Inf.
This differs from the Yee stencil bound and does not guarantee stability under
arbitrary variable-step sequences or prescribe an accurate timestep.

All three curls and scalar stages write model-private fields. Only after the
last collective success are all six physical arrays and time/step committed;
their identities stay fixed. Recoverable local validation/numerical failure
preserves physical state on every rank. Private work may change and is overwritten
on reuse. MPI process/transport loss and fatal resource failures retain the
operator contract's explicit exclusions.

## Diagnostics, cost and remaining implementation

`em_field_energy(model)` is collective. It reuses canonical exact owned-cell
energy and sums serialized scalar rank contributions before rounding once to
the field precision; an unrepresentable total returns Inf. No field is gathered.
`distributed_maxwell_diagnostics(model)` adds synchronized time/step and the
global maximum absolute E/B divergence from the existing divergence operator.
Diagnostics do not project or modify physical fields.

The model owns six physical and nine private real arrays on owned cells, in
addition to the existing grid's FFT/operator work and replicated axis metadata.
One step invokes three curls: nine forward and nine inverse transforms. Exact
active-row arithmetic has allocation and computational cost; only coefficient
and array reuse are implemented here. Diagnostics exchange scalar energy payloads
in rank order and are outside the step hot path. No timing, scaling, allocation
ceiling or optimized/zero-allocation claim has been measured.

This initial field slice supplies no checkpoint/archive adapter. A physical-pencil
rank-generation checkpoint with matching layout is the next source implementation;
plans/communicators must not be serialized or passed through the hybrid schema.
Changed-layout restart, coupled particles, mixed/nonperiodic boundaries, backends
and multi-node scaling remain separate CAP-003 work. This source implementation
does not complete that capability.
