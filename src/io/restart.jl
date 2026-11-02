# repro.jl — reproducibility archive (§21.3), operator migration check (Phase-1),
# and sampled particle dumps (§6).
#
# Three additive utilities, all built on existing infrastructure (no new deps):
#
#   * archive_run / load_archive — a provenance-stamped, checksum-verified snapshot
#     of a supported simulation state. This layers the run-metadata
#     wrapper (save_run/load_run + capture_metadata/RunMetadata from metadata.jl)
#     over the same deterministic state that checkpoint.jl restarts from, so an
#     archive is self-describing (who/when/which commit/which seed) AND integrity
#     checked (hash of the serialized state).
#
#   * sample_particles — every stride-th particle's positions and velocities, for
#     SAMPLED phase-space dumps (so a dump is a deliberate decimation, not a
#     full indiscriminate copy of millions of particles).
#
#   * operators_match — verifies the extracted SpectralOperators.deriv! produces
#     results IDENTICAL (to roundoff) to HybridPlasmaPIC.deriv! for a random field, i.e.
#     the Phase-1 operator extraction did not change behaviour. HybridPlasmaPIC
#     re-exports SpectralOperators.deriv!, so this guards against future import
#     drift.
#
# Serialization is already imported module-wide (checkpoint.jl); no new dependency.

# ---------------------------------------------------------------- archive

"""
    archive_run(path, st, ps::ParticleSet; rng_seed,
                normalization="Omega_ci", filter_desc="", boundary_desc="periodic",
                diagnostic_desc="", rank_layout="") -> path
    archive_run(path, empic, electrons, ions; rng_seed, kwargs...) -> path
    archive_run(path, hall::HallMHDState; rng_seed, kwargs...) -> path
    archive_run(path, model::MaxwellVacuum; rng_seed, kwargs...) -> path
    archive_run(path, stepper, species::AbstractVector; rng_seed, kwargs...) -> path

Write a provenance-stamped, checksum-verified archive of the run state at `path`.
Supports [`HybridStepper`](@ref), [`CAMCLStepper`](@ref),
[`Electrostatic1D`](@ref), [`ElectrostaticPIC`](@ref), [`EMPIC1D`](@ref), [`EMPIC`](@ref),
[`HallMHDState`](@ref), [`MaxwellVacuum`](@ref), [`CollisionalHybridStepper`](@ref), [`CollisionalEMPICStepper`](@ref),
[`CollisionalElectrostaticPIC`](@ref),
[`YeeConductiveLayer`](@ref),
[`YeePML`](@ref),
[`YeeLongitudinalPIC`](@ref),
[`YeeElectromagneticPIC`](@ref), [`YeeImpactLongitudinalPIC`](@ref),
and [`YeeAbsorbingLongitudinalPIC`](@ref).

The archive captures the deterministic simulation state (model and shape
configuration, particle phase space + provenance, the canonical magnetic field
`st.fields.B`, the carried electric field `st.fields.E`, the grid identity
`st.g.n`/`st.g.L`, the `time`/`step` counters, and `last_dt` for the particle
velocity time level, plus the persistent magnetic copies and subcycle phase for
CAM-CL — the same state [`save_checkpoint`](@ref) restarts from) together with a
[`RunMetadata`](@ref) record built via [`capture_metadata`](@ref). It is written
through [`save_run`](@ref), so the file carries the schema tag and an integrity
checksum over the state.
Hybrid species collections preserve the ordered per-population state alongside
the common fields, configuration and integrator phase. Evolving electron
pressure is included; scalar transport closures also store their carried ion
force and `transport_policy_version`. `ConductingElectrons` additionally requires
integer `conduction_policy_version=2`, including for κ=0; γ/κ remain in its
saved model and quotient/transform scratch is not persisted. Policy-1 conducting
state cannot be restored under policy 2. All are covered by the unchanged
schema-4 checksum; `load_archive` checks integrity, not destination compatibility.
Electrostatic archives instead contain their electric field, density, immutable
background density, and priming/timestep state alongside the particles and grid.
Electrostatic collections retain ordered species, neutrality policy and charge
tolerance through the same collection checkpoint schema.
Electromagnetic archives include both species for mobile runs, their physics
configuration, electric/magnetic fields, continuity diagnostic snapshots, and
full/fast timestep history. Supply ions exactly when `mobile=true` for the legacy
electron/ion API. Synchronous `EMPIC` collections preserve their ordered arbitrary
populations, mobility, neutrality/background policy, relativistic mode and c,
and common fast velocity half level through collection checkpoint schema 1.
Opted-in generic dividing schedules use schema2/policy1 and per-species strides
and completed intervals. Their destinations must also be opted in; validated
generic schema1 without schedule keys migrates explicitly to all-fast rows there.
`EMPIC1D` collections use distinct schema 2 `:empic1d_collection`, required
`polarization=:ex_ey_bz` and integer `schedule_policy=1`, ordered species strides
and completed species intervals, and their seven specialized physical field arrays.
Validated specialized schema-1 state migrates explicitly to all-fast strides;
schema-1 payloads carrying schema-2 schedule fields are rejected.
Scalar, generic collection and specialized collection payloads cannot be
interchanged; the surrounding archive container remains schema 4.
Hall-MHD archives contain the model, fluid and public derived fields, retained
stage `lapJ`, and clock; no particle argument or timestep history is needed.
Collisional archives add the complete advanced RNG state and completed source
ledger through the same versioned state owner as stochastic checkpoints.
Collisional EMPIC archives use their distinct schema-1 field-then-collision policy
and nested generic EMPIC physical state, including the preceding drift current.
They retain the full advanced RNG, fixed configuration and completed source
ledger; restored fields and charge/current snapshots are copied verbatim.
The explicit positive collision permittivity and count conversion must be equal
at the stored precision, as required by the Maxwell source normalization.
Collisional ES archives use distinct `:collisional_espic` schema 1 with the nested
ES E/ne/particle/phase state, full RNG, fixed equal positive collision permittivity
and count conversion, and sole wrapper source clock/count. They restore physical
fields verbatim and retain the new-position/carried-half collision policy.
Field-only `MaxwellVacuum` archives retain synchronized fields, clock and actual
Yee grid/boundary/speed configuration in a distinct schema, with no particle or
waveform payload. Their defaults describe vacuum units and actual wall axes.
`YeeConductiveLayer` archives additionally retain fixed scalar-conductivity
profile/timestep policies and every exact loss/work/defect ledger in a distinct
schema. Derived caches are rebuilt; no PML auxiliary state or forcing waveform
is invented. The existing container remains version 4.
`YeePML` archives retain every actual directional ADE memory and signed exact
exchange/field-defect history under their distinct policy. No forcing sequence,
scalar-layer substitution or positive auxiliary-energy law is inferred.
Closed-Yee archives retain their distinct completed-state schemas and restore
physical fields verbatim. Impact longitudinal state additionally retains its
root-refinement budget, exact wall/electric impulses, counts and carried error
bounds, with required elastic zero charge/work records and policies. Its schema 1
cannot be interchanged with event-free schema 2; the archive container remains 4.
Absorbing longitudinal archives use a distinct schema-1 tag and retain initial/live
source records, both dynamic free plates, local capture records and all signed
capture/energy ledgers. Their particle destinations require ordinary Vectors;
reflection/event-free payloads and matched-view destinations cannot be substituted.

`rng_seed` records the initialization seed as provenance in the metadata. A seed
does not replace the advanced RNG state needed for stochastic continuation.
The remaining keyword arguments describe the units convention and the
filter/boundary/diagnostic configuration for the audit trail. Returns `path`.
`rank_layout` records the serial/MPI rank topology in the same metadata record;
leave it empty for the serial default.
"""
function archive_run(path::AbstractString, st::_CheckpointStepper, ps::ParticleSet; kwargs...)
    return _archive_run(path, st, _checkpoint_state(st, ps); kwargs...)
end

function archive_run(path::AbstractString, st::HallMHDState; kwargs...)
    return _archive_run(path, st, _checkpoint_state(st); kwargs...)
end

function archive_run(
    path::AbstractString,
    st::Union{_HybridCheckpointStepper,ElectrostaticPIC,EMPIC,EMPIC1D},
    species::AbstractVector{<:ParticleSet};
    kwargs...,
)
    return _archive_run(path, st, _checkpoint_state(st, species); kwargs...)
end

function archive_run(
    path::AbstractString,
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet};
    kwargs...,
)
    return _archive_run(path, st, _checkpoint_state(st, electrons, ions); kwargs...)
end

function archive_run(
    path::AbstractString,
    st::CollisionalHybridStepper,
    species::AbstractVector{<:ParticleSet};
    kwargs...,
)
    _acquire_collisional!(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, species); kwargs...)
    finally
        st.busy[] = false
    end
end

function archive_run(
    path::AbstractString,
    st::CollisionalEMPICStepper,
    species::AbstractVector{<:ParticleSet};
    kwargs...,
)
    _acquire_collisional_empic!(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, species); kwargs...)
    finally
        st.busy[] = false
    end
end

function archive_run(
    path::AbstractString,
    st::CollisionalElectrostaticPIC,
    species::AbstractVector{<:ParticleSet};
    kwargs...,
)
    _acquire_collisional_espic!(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, species); kwargs...)
    finally
        st.busy[] = false
    end
end

function _archive_run(
    path::AbstractString,
    st::Union{
        _CheckpointStepper,
        HallMHDState,
        MaxwellVacuum,
        YeeConductiveLayer,
        YeePML,
        CollisionalHybridStepper,
        CollisionalEMPICStepper,
        CollisionalElectrostaticPIC,
        ImplicitHybridStepper,
        ImplicitEMPICStepper,
        ImplicitYeeEMPICStepper,
        SlabGyrokineticPIC,
        OpenHallHybrid,
        ShearedPeriodicHybridPIC,
        StaticCompositeElectrostaticPIC,
        _ClosedYeeCheckpointStepper,
    },
    state;
    rng_seed::Integer,
    normalization::AbstractString = "Omega_ci",
    filter_desc::AbstractString = "",
    boundary_desc::AbstractString = "periodic",
    diagnostic_desc::AbstractString = "",
    rank_layout::AbstractString = "",
)
    meta = capture_metadata(;
        rng_seed = rng_seed,
        normalization = normalization,
        filter_desc = filter_desc,
        boundary_desc = boundary_desc,
        diagnostic_desc = diagnostic_desc,
        rank_layout = rank_layout,
    )
    _validate_checkpoint_state(state, st)
    return save_run(path, state, meta)
end

"""
    load_archive(path) -> (; meta, state)

Read an archive written by [`archive_run`](@ref). Delegates to [`load_run`](@ref),
which verifies the schema version and recomputes the stored checksum (throwing on
version skew or corruption). Returns the [`RunMetadata`](@ref) `meta` and the
deserialized run `state` NamedTuple.
"""
function load_archive(path::AbstractString)
    loaded = load_run(path)
    return (meta = loaded.meta, state = loaded.state)
end

# ---------------------------------------------------------------- sampled dumps

"""
    sample_particles(ps::ParticleSet{D,T}, stride::Int)
        -> (; index, x::NTuple{D,Vector{T}}, v::NTuple{3,Vector{T}})

Return every `stride`-th particle's positions and velocities — a decimated copy
for SAMPLED phase-space dumps (so a dump is a deliberate sub-sampling rather than
an indiscriminate full copy). With `N` particles this yields `cld(N, stride)` ≈
`N/stride` entries: particles `1, 1+stride, 1+2·stride, …`.

The returned `index` holds the source particle indices; `x[d]`/`v[c]` are fresh
`Vector`s (copies, not views) so the dump is decoupled from later mutation of
`ps`. `stride` must be ≥ 1.
"""
function sample_particles(ps::ParticleSet{D,T}, stride::Int) where {D,T}
    stride >= 1 || throw(ArgumentError("stride must be ≥ 1, got $stride"))
    N = nparticles(ps)
    idx = collect(1:stride:N)                      # 1, 1+stride, …  (cld(N,stride) entries)
    x = ntuple(d -> ps.x[d][idx], D)               # indexing a Vector by an index Vector copies
    v = ntuple(c -> ps.v[c][idx], 3)
    return (index = idx, x = x, v = v)
end

# ---------------------------------------------------------------- migration check

"""
    operators_match(g::FourierGrid) -> Bool

Verify that the extracted [`SpectralOperators`](@ref) derivative agrees with the
package's [`deriv!`](@ref) — the Phase-1 operator-migration guard. For a random
real field on `g`, computes the first derivative along every axis with both
`SpectralOperators.deriv!` and `HybridPlasmaPIC.deriv!` and returns `true` iff they are
identical to roundoff (max absolute difference ≤ `8·eps(T)·(1 + max|result|)`,
i.e. a few ULP scaled by the field magnitude).

Since `HybridPlasmaPIC.deriv!` is re-exported from `SpectralOperators.deriv!`,
agreement is exact in practice; the tolerance only guards against future import
drift. Uses a separate `FourierGrid` for each operator so their in-place scratch
buffers cannot interfere.
"""
function operators_match(g::FourierGrid{D,T}) where {D,T}
    # Independent grids so the two in-place operators do not share scratch.
    g1 = FourierGrid(g.n, g.L)
    g2 = FourierGrid(g.n, g.L)
    rng = Random.MersenneTwister(0x5eed)           # deterministic test field
    f = rand(rng, T, g.n...)
    a = similar(f)
    b = similar(f)
    tol = T(8) * eps(T)
    @inbounds for j = 1:D
        HybridPlasmaPIC.deriv!(a, f, g1, j)
        SpectralOperators.deriv!(b, f, g2, j)
        scale = one(T) + max(maximum(abs, a), maximum(abs, b))
        maxdiff = maximum(abs.(a .- b))
        maxdiff <= tol * scale || return false
    end
    return true
end
