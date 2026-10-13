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

Write a provenance-stamped, checksum-verified archive of the run state at `path`.
Supports [`HybridStepper`](@ref), [`CAMCLStepper`](@ref),
[`Electrostatic1D`](@ref), [`ElectrostaticPIC`](@ref), [`EMPIC1D`](@ref), [`EMPIC`](@ref),
and [`HallMHDState`](@ref).

The archive captures the deterministic simulation state (model and shape
configuration, particle phase space + provenance, the canonical magnetic field
`st.fields.B`, the carried electric field `st.fields.E`, the grid identity
`st.g.n`/`st.g.L`, the `time`/`step` counters, and `last_dt` for the particle
velocity time level, plus the persistent magnetic copies and subcycle phase for
CAM-CL — the same state [`save_checkpoint`](@ref) restarts from) together with a
[`RunMetadata`](@ref) record built via [`capture_metadata`](@ref). It is written
through [`save_run`](@ref), so the file carries the schema tag and an integrity
checksum over the state.
Electrostatic archives instead contain their electric field, density, immutable
background density, and priming/timestep state alongside the particles and grid.
Electromagnetic archives include both species for mobile runs, their physics
configuration, electric/magnetic fields, continuity diagnostic snapshots, and
full/fast timestep history. Supply ions exactly when `mobile=true`.
Hall-MHD archives contain the model, fluid and public derived fields, retained
stage `lapJ`, and clock; no particle argument or timestep history is needed.

`rng_seed` is the seed that makes the run reproducible (stored in the metadata).
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
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet};
    kwargs...,
)
    return _archive_run(path, st, _checkpoint_state(st, electrons, ions); kwargs...)
end

function _archive_run(
    path::AbstractString,
    st::Union{_CheckpointStepper,HallMHDState},
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
