# checkpoint.jl — checkpoint/restart of the deterministic simulation state.
#
# Stepping (step!) uses no RNG, so the state needed for a bitwise-identical
# restart is: the particle phase space + provenance, the canonical magnetic
# field and carried electric field, and the time/step counters. The frozen
# n+1/2 moments are scratch recomputed every step, so they are not stored.
# ponytail: stdlib Serialization (no extra dep). Swap for HDF5 when cross-version
# or cross-language portability is needed.

using Serialization

function _validate_checkpoint_particle_state(s, ::Val{D}) where {D}
    length(s.x) == D || throw(
        ArgumentError("checkpoint has $(length(s.x)) position arrays, expected $D"),
    )
    length(s.v) == 3 || throw(
        ArgumentError("checkpoint has $(length(s.v)) velocity arrays, expected 3"),
    )
    N = length(s.weight)
    for d = 1:D
        length(s.x[d]) == N || throw(
            ArgumentError("checkpoint x[$d] length $(length(s.x[d])) ≠ particle count $N"),
        )
        _check_particle_vector_axes(Symbol(:x, d), s.x[d], N)
    end
    for c = 1:3
        length(s.v[c]) == N || throw(
            ArgumentError("checkpoint v[$c] length $(length(s.v[c])) ≠ particle count $N"),
        )
        _check_particle_vector_axes(Symbol(:v, c), s.v[c], N)
    end
    length(s.id) == N || throw(
        ArgumentError("checkpoint id length $(length(s.id)) ≠ particle count $N"),
    )
    length(s.tag) == N || throw(
        ArgumentError("checkpoint tag length $(length(s.tag)) ≠ particle count $N"),
    )
    _check_particle_vector_axes(:weight, s.weight, N)
    _check_particle_vector_axes(:id, s.id, N)
    _check_particle_vector_axes(:tag, s.tag, N)
    return nothing
end

function _validate_checkpoint_field_state(s, st)
    length(s.B) == 3 || throw(ArgumentError("checkpoint has $(length(s.B)) B arrays, expected 3"))
    length(s.E) == 3 || throw(ArgumentError("checkpoint has $(length(s.E)) E arrays, expected 3"))
    for c = 1:3
        axes(s.B[c]) == axes(st.fields.B[c]) || throw(
            ArgumentError("checkpoint B[$c] axes $(axes(s.B[c])) ≠ $(axes(st.fields.B[c]))"),
        )
        axes(s.E[c]) == axes(st.fields.E[c]) || throw(
            ArgumentError("checkpoint E[$c] axes $(axes(s.E[c])) ≠ $(axes(st.fields.E[c]))"),
        )
    end
    return nothing
end

"""
    save_checkpoint(path, stepper, ps)

Serialize the restartable state (particles, B, carried E, time, step) to `path`.
"""
function save_checkpoint(
    path::AbstractString,
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
) where {D,T}
    state = (
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        x = ps.x,
        v = ps.v,
        weight = ps.weight,
        id = ps.id,
        tag = ps.tag,
        q = ps.q,
        m = ps.m,
        B = st.fields.B,
        E = st.fields.E,
        time = st.time[],
        step = st.step[],
    )
    serialize(path, state)
    return path
end

"""
    load_checkpoint!(stepper, ps, path)

Restore a checkpoint into a `stepper`/`ps` of matching dimension and grid size
(particle arrays are resized to match). After loading, `step!` continues
bitwise-identically to the original run.
"""
function load_checkpoint!(
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
    path::AbstractString,
) where {D,T}
    s = deserialize(path)
    s.D == D || throw(ArgumentError("checkpoint dimension $(s.D) ≠ $D"))
    s.T == T || throw(
        ArgumentError(
            "checkpoint eltype $(s.T) ≠ $T (would silently convert and break bitwise restart)",
        ),
    )
    s.ncell == st.g.n || throw(ArgumentError("checkpoint grid $(s.ncell) ≠ $(st.g.n)"))
    hasproperty(s, :L) || throw(
        ArgumentError("checkpoint is missing grid lengths L and cannot guarantee a bitwise-identical restart"),
    )
    s.L == st.g.L || throw(ArgumentError("checkpoint box lengths $(s.L) ≠ $(st.g.L)"))
    _validate_checkpoint_particle_state(s, Val(D))
    _validate_checkpoint_field_state(s, st)
    for d = 1:D
        resize!(ps.x[d], length(s.x[d]))
        copyto!(ps.x[d], s.x[d])
    end
    for c = 1:3
        resize!(ps.v[c], length(s.v[c]))
        copyto!(ps.v[c], s.v[c])
    end
    resize!(ps.weight, length(s.weight))
    copyto!(ps.weight, s.weight)
    resize!(ps.id, length(s.id))
    copyto!(ps.id, s.id)
    resize!(ps.tag, length(s.tag))
    copyto!(ps.tag, s.tag)
    ps.q = s.q
    ps.m = s.m
    for c = 1:3
        copyto!(st.fields.B[c], s.B[c])
        copyto!(st.fields.E[c], s.E[c])
    end
    st.time[] = s.time
    st.step[] = s.step
    return st
end
