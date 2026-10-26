# checkpoint.jl — checkpoint/restart of the deterministic simulation state.
#
# Stepping (step!) uses no RNG, so the state needed for a bitwise-identical
# restart is: model and shape, particle phase space + provenance, canonical B,
# carried E, evolving electron pressure when present, time/step counters, and
# the previous particle timestep. The frozen n+1/2 moments are scratch
# recomputed every step, so they are not stored.
# CAM-CL also retains its two magnetic copies, prior output, and subcycle phase.
# Electrostatic PIC retains E, density, background density, and priming/timestep
# state; the per-particle field gather buffers are scratch.
# Electromagnetic PIC also retains each species separately, its fast/full-step
# phase, and the charge snapshots/current used by continuity diagnostics.
# Hall-MHD retains the fluid fields, their public derived values, model, and
# clock. Inactive stage lapJ can be carried into the public fields; the other RK
# stage values are scratch. Hall-MHD carries no staggered timestep history.
# ponytail: stdlib Serialization (no extra dep). Swap for HDF5 when cross-version
# or cross-language portability is needed.

using Serialization

const _HybridCheckpointStepper{D,T} = Union{HybridStepper{D,T},CAMCLStepper{D,T}}
const _ElectrostaticCheckpointStepper = Union{Electrostatic1D,ElectrostaticPIC}
const _ElectromagneticCheckpointStepper = Union{EMPIC1D,EMPIC}
const _CheckpointStepper = Union{
    _HybridCheckpointStepper,
    _ElectrostaticCheckpointStepper,
    _ElectromagneticCheckpointStepper,
}

const _CHECKPOINT_REQUIRED_FIELDS = (
    :integrator,
    :D,
    :T,
    :ncell,
    :L,
    :model,
    :shape,
    :x,
    :v,
    :weight,
    :id,
    :tag,
    :q,
    :m,
    :B,
    :E,
    :time,
    :step,
    :last_dt,
)

const _ESPIC_CHECKPOINT_REQUIRED_FIELDS = (
    :integrator,
    :D,
    :T,
    :ncell,
    :L,
    :n0,
    :x,
    :v,
    :weight,
    :id,
    :tag,
    :q,
    :m,
    :E,
    :ne,
    :primed,
    :last_dt,
)

const _EMPIC_CHECKPOINT_REQUIRED_FIELDS = (
    :integrator,
    :D,
    :T,
    :ncell,
    :L,
    :n0,
    :c,
    :shape,
    :relativistic,
    :mobile,
    :mi,
    :n_sub,
    :electrons,
    :ions,
    :fields,
    :time,
    :step,
    :last_dt,
    :last_sub_dt,
)

const _CHECKPOINT_PARTICLE_FIELDS = (:x, :v, :weight, :id, :tag, :q, :m)

const _COLLECTION_CHECKPOINT_REQUIRED_FIELDS = (
    filter(k -> !(k in _CHECKPOINT_PARTICLE_FIELDS), _CHECKPOINT_REQUIRED_FIELDS)...,
    :particle_layout,
    :species,
)

const _HALL_CHECKPOINT_REQUIRED_FIELDS =
    (:integrator, :D, :T, :ncell, :L, :model, :fields, :stage_lapJ, :floor_count, :time, :step)

function _validate_checkpoint_container(s, required = _CHECKPOINT_REQUIRED_FIELDS)
    s isa NamedTuple ||
        throw(ArgumentError("checkpoint file does not contain a valid checkpoint container"))
    if !all(map(k -> haskey(s, k), required))
        missing = Symbol[k for k in required if !haskey(s, k)]
        missing_list = join(string.(missing), ", ")
        throw(ArgumentError("checkpoint is missing required fields: $missing_list"))
    end
    return nothing
end

function _validate_checkpoint_vector(name::Symbol, a, N::Int, ::Type{T}) where {T}
    a isa AbstractVector ||
        throw(ArgumentError("checkpoint $(name) is $(typeof(a)), expected AbstractVector{$T}"))
    eltype(a) === T || throw(
        ArgumentError(
            "checkpoint $(name) eltype $(eltype(a)) ≠ $T (would silently convert on restore)",
        ),
    )
    length(a) == N ||
        throw(ArgumentError("checkpoint $(name) length $(length(a)) ≠ particle count $N"))
    _check_particle_vector_axes(name, a, N)
    return nothing
end

function _validate_checkpoint_array(name::Symbol, a, reference::AbstractArray{T}) where {T}
    a isa AbstractArray ||
        throw(ArgumentError("checkpoint $(name) is $(typeof(a)), expected AbstractArray{$T}"))
    eltype(a) === T || throw(
        ArgumentError(
            "checkpoint $(name) eltype $(eltype(a)) ≠ $T (would silently convert on restore)",
        ),
    )
    axes(a) == axes(reference) ||
        throw(ArgumentError("checkpoint $(name) axes $(axes(a)) ≠ $(axes(reference))"))
    return nothing
end

function _validate_checkpoint_particle_scalars(s, ::Type{T}) where {T}
    s.q isa T || throw(
        ArgumentError("checkpoint q type $(typeof(s.q)) ≠ $T (would silently convert on restore)"),
    )
    s.m isa T || throw(
        ArgumentError("checkpoint m type $(typeof(s.m)) ≠ $T (would silently convert on restore)"),
    )
    return nothing
end

function _validate_checkpoint_scalar_state(s, ::Type{T}) where {T}
    _validate_checkpoint_particle_scalars(s, T)
    _validate_checkpoint_time_state(s, T)
    return nothing
end

function _validate_checkpoint_clock_state(s, ::Type{T}) where {T}
    s.time isa T || throw(
        ArgumentError(
            "checkpoint time type $(typeof(s.time)) ≠ $T (would silently convert on restore)",
        ),
    )
    _require_finite_nonnegative_real("checkpoint time", s.time, T)
    s.step isa Int || throw(
        ArgumentError(
            "checkpoint step type $(typeof(s.step)) ≠ Int (would fail or convert on restore)",
        ),
    )
    s.step >= 0 || throw(ArgumentError("checkpoint step must be nonnegative"))
    return nothing
end

function _validate_checkpoint_time_state(s, ::Type{T}) where {T}
    _validate_checkpoint_clock_state(s, T)
    s.last_dt isa T || throw(ArgumentError("checkpoint last_dt type $(typeof(s.last_dt)) ≠ $T"))
    isfinite(s.last_dt) && s.last_dt >= zero(T) ||
        throw(ArgumentError("checkpoint last_dt must be finite and nonnegative"))
    (s.step == 0 ? iszero(s.last_dt) : s.last_dt > zero(T)) || throw(
        ArgumentError("checkpoint last_dt must be zero before stepping and positive afterward"),
    )
    return nothing
end

_validate_checkpoint_particle_state(s, dims::Val) =
    _validate_checkpoint_particle_state(s, dims, s.T)

function _validate_checkpoint_particle_state(s, ::Val{D}, ::Type{T}) where {D,T}
    s.x isa Tuple || throw(ArgumentError("checkpoint positions must be stored as a tuple"))
    s.v isa Tuple || throw(ArgumentError("checkpoint velocities must be stored as a tuple"))
    length(s.x) == D ||
        throw(ArgumentError("checkpoint has $(length(s.x)) position arrays, expected $D"))
    length(s.v) == 3 ||
        throw(ArgumentError("checkpoint has $(length(s.v)) velocity arrays, expected 3"))
    s.weight isa AbstractVector ||
        throw(ArgumentError("checkpoint weight is $(typeof(s.weight)), expected AbstractVector"))
    N = length(s.weight)
    _validate_checkpoint_vector(:weight, s.weight, N, T)
    for d = 1:D
        _validate_checkpoint_vector(Symbol(:x, d), s.x[d], N, T)
    end
    for c = 1:3
        _validate_checkpoint_vector(Symbol(:v, c), s.v[c], N, T)
    end
    _validate_checkpoint_vector(:id, s.id, N, UInt64)
    _validate_checkpoint_vector(:tag, s.tag, N, UInt32)
    return nothing
end

function _checkpoint_species(s, dims::Val{D}, ::Type{T}) where {D,T}
    _validate_checkpoint_container(s, _CHECKPOINT_PARTICLE_FIELDS)
    _validate_checkpoint_particle_state(s, dims, T)
    _validate_checkpoint_particle_scalars(s, T)
    # A validated wrapper reuses the mover's species and relativistic checks
    # without copying particle arrays or converting saved scalar types.
    return ParticleSet{D,T}(s.x, s.v, s.weight, s.id, s.tag, s.q, s.m)
end

function _validate_checkpoint_field_components(name::Symbol, components, reference)
    components isa Tuple ||
        throw(ArgumentError("checkpoint $name fields must be stored as a tuple"))
    length(components) == length(reference) || throw(
        ArgumentError(
            "checkpoint has $(length(components)) $name arrays, expected $(length(reference))",
        ),
    )
    for c in eachindex(reference)
        _validate_checkpoint_array(Symbol(name, c), components[c], reference[c])
    end
    return nothing
end

function _validate_checkpoint_fields(saved, fields::NamedTuple)
    _validate_checkpoint_container(saved, keys(fields))
    for (name, reference) in pairs(fields)
        if reference isa Tuple
            _validate_checkpoint_field_components(name, getproperty(saved, name), reference)
        else
            _validate_checkpoint_array(name, getproperty(saved, name), reference)
        end
    end
    return nothing
end

function _restore_checkpoint_fields!(fields::NamedTuple, saved)
    saved_fields = NamedTuple{keys(fields)}(saved)
    map(values(fields), values(saved_fields)) do destination, source
        if destination isa Tuple
            map(copyto!, destination, source)
        else
            copyto!(destination, source)
        end
    end
    return nothing
end

function _validate_checkpoint_field_state(s, st::_HybridCheckpointStepper{D,T}) where {D,T}
    _validate_checkpoint_field_components(:B, s.B, st.fields.B)
    _validate_checkpoint_field_components(:E, s.E, st.fields.E)
    if _evolves_pressure(s.model.closure)
        _validate_checkpoint_container(s, (:pe,))
        s.pe isa Array{T,D} || throw(ArgumentError("checkpoint pe must be stored as Array{$T,$D}"))
        _validate_evolving_pressure(s.pe, st.g; context = "checkpoint electron pressure")
    end
    return nothing
end

function _validate_checkpoint_integrator_state(s, st::_HybridCheckpointStepper{D,T}) where {D,T}
    expected = st isa CAMCLStepper ? :camcl : :rk4
    s.integrator === expected ||
        throw(ArgumentError("checkpoint integrator $(s.integrator) ≠ $expected"))
    if st isa CAMCLStepper
        _validate_checkpoint_container(s, (:Blead, :Blag, :Bout, :cl_h, :cl_live))
        _validate_checkpoint_fields(s, (Blead = st.Blead, Blag = st.Blag, Bout = st.Bout))
        s.cl_h isa T || throw(ArgumentError("checkpoint cl_h type $(typeof(s.cl_h)) ≠ $T"))
        s.cl_live isa Bool || throw(ArgumentError("checkpoint cl_live must be a Bool"))
        isfinite(s.cl_h) && s.cl_h >= zero(T) ||
            throw(ArgumentError("checkpoint cl_h must be finite and nonnegative"))
        s.cl_live &&
            iszero(s.cl_h) &&
            throw(ArgumentError("checkpoint live cyclic-leapfrog state requires positive cl_h"))
    end
    return nothing
end

function _validate_checkpoint_grid(s, g::FourierGrid{D,T}) where {D,T}
    s.D == D || throw(ArgumentError("checkpoint dimension $(s.D) ≠ $D"))
    s.T == T || throw(
        ArgumentError(
            "checkpoint eltype $(s.T) ≠ $T (would silently convert and break bitwise restart)",
        ),
    )
    s.ncell == g.n || throw(ArgumentError("checkpoint grid $(s.ncell) ≠ $(g.n)"))
    s.L == g.L || throw(ArgumentError("checkpoint box lengths $(s.L) ≠ $(g.L)"))
    return nothing
end

# A raw pointer is process-local even when it is nested in an isbits parameter.
# Ordinary isbits configuration still avoids graph allocation entirely.
function _checkpoint_pointer_free_bits(::Type{T}) where {T}
    isbitstype(T) || return false
    T <: Ptr && return false
    for field in fieldtypes(T)
        _checkpoint_pointer_free_bits(field) || return false
    end
    return true
end

_checkpoint_backing_arrays(value::BitArray) = (value.chunks,)
_checkpoint_backing_arrays(value::Dict) = (value.slots, value.keys, value.vals)
_checkpoint_backing_arrays(value::IdDict) = (value.ht,)
_checkpoint_backing_arrays(value::Expr) = (value.args,)

function _queue_checkpoint_object!(pending, edges, item, logical, stored = true)
    push!(pending, (item, logical))
    edges === nothing || push!(edges, (item, stored))
    return nothing
end

# Walk logical storage and stored metadata. Array parents do not expand a view's
# logical extent, and standard container capacity stays opaque. A later logical
# visit upgrades a metadata visit; each object is expanded at most once per mode.
# Optional edge records distinguish actual serialized fields from synthetic
# custom-container iteration, without invoking a serializer or constructor.
function _visit_checkpoint_objects(f, roots; backing = nothing, record = nothing)
    pending = Tuple{Any,Bool}[(roots, true)]
    seen = IdDict{Any,Bool}()
    while !isempty(pending)
        value, logical = pop!(pending)
        T = typeof(value)
        (
            (isbitstype(T) && (record === nothing || _checkpoint_pointer_free_bits(T))) ||
            value isa Union{String,Symbol,Module,Type}
        ) && continue
        haskey(seen, value) && (seen[value] || !logical) && continue
        seen[value] = logical
        logical && f(value)
        edges = record === nothing ? nothing : Tuple{Any,Bool}[]
        copied_view = value isa SubArray && parent(value) isa Array
        if value isa BitArray
            for buffer in _checkpoint_backing_arrays(value)
                backing === nothing || backing(buffer)
                record === nothing || _queue_checkpoint_object!(pending, edges, buffer, false)
            end
        elseif value isa AbstractArray && !(value isa AbstractRange)
            if logical && (
                !isbitstype(eltype(value)) ||
                (record !== nothing && !_checkpoint_pointer_free_bits(eltype(value)))
            )
                stored = value isa Union{Array,Memory} || copied_view
                for i in eachindex(value)
                    isassigned(value, i) &&
                        _queue_checkpoint_object!(pending, edges, value[i], true, stored)
                end
            end
            if !(value isa Union{Array,Memory})
                array_parent = parent(value)
                array_parent === value ||
                    _queue_checkpoint_object!(pending, edges, array_parent, false, !copied_view)
                for i = 1:fieldcount(T)
                    isdefined(value, i) || continue
                    item = getfield(value, i)
                    array_parent !== value && item === array_parent && continue
                    # Array-backed views serialize fresh parents and indices.
                    _queue_checkpoint_object!(pending, edges, item, true, !copied_view)
                end
            end
        elseif value isa AbstractDict
            stored = value isa Union{Dict,IdDict}
            for (key, item) in value
                _queue_checkpoint_object!(pending, edges, key, true, stored)
                _queue_checkpoint_object!(pending, edges, item, true, stored)
            end
            if stored
                if backing !== nothing
                    foreach(backing, _checkpoint_backing_arrays(value))
                end
            else
                for i = 1:fieldcount(T)
                    isdefined(value, i) &&
                        _queue_checkpoint_object!(pending, edges, getfield(value, i), true)
                end
            end
        elseif value isa Expr && record !== nothing
            # Expr serializes each argument, rebuilding its args vector.
            _queue_checkpoint_object!(pending, edges, value.head, true)
            _queue_checkpoint_object!(pending, edges, value.args, true, false)
            for item in value.args
                _queue_checkpoint_object!(pending, edges, item, true)
            end
        elseif value isa Set
            _queue_checkpoint_object!(pending, edges, value.dict, true)
        elseif value isa Core.SimpleVector
            for item in value
                _queue_checkpoint_object!(pending, edges, item, true)
            end
        else
            if value isa AbstractSet
                for item in value
                    _queue_checkpoint_object!(pending, edges, item, true, false)
                end
            end
            for i = 1:fieldcount(T)
                isdefined(value, i) &&
                    _queue_checkpoint_object!(pending, edges, getfield(value, i), true)
            end
        end
        record === nothing || record(value, logical, edges)
    end
    return nothing
end

# The owner is the object whose serialized data allocation will hold the logical
# array after loading. Generic wrappers preserve their declared parent; Arrays,
# Memory and each serialized Array-backed SubArray allocate their own data.
function _checkpoint_serialization_owner(array::AbstractArray)
    seen = IdSet{Any}()
    while true
        array in seen && throw(ArgumentError("checkpoint array parent cycle"))
        push!(seen, array)
        array isa Union{Array,Memory} && return array
        array isa BitArray && return array.chunks
        array isa SubArray && parent(array) isa Array && return array
        next = parent(array)
        next === array && return array
        array = next
    end
end

function _validate_checkpoint_configuration_graph(configuration, graph, require_independent)
    # Logical iteration can expose a raw pointer that lives in a wrapper's data.
    # Built-in value codecs reconstruct their private native resources instead.
    reached = IdSet{Any}()
    pending = Any[configuration]
    while !isempty(pending)
        value = pop!(pending)
        value in reached && continue
        push!(reached, value)
        haskey(graph, value) || continue
        if value isa FFTW.FFTWPlan || (value isa Ptr && value != C_NULL)
            throw(
                ArgumentError("checkpoint configuration contains a process-local native resource"),
            )
        end
        for (item, _) in graph[value].edges
            push!(pending, item)
        end
    end

    # Serialization tracks most mutable objects once. Immutable containers and
    # the BigInt/Regex value codecs are emitted per incoming edge. Counts capped
    # at two expose copied aliases while keeping fan-out and cycles bounded.
    counts = IdDict{Any,Int}()
    pending = Any[configuration]
    while !isempty(pending)
        value = pop!(pending)
        haskey(graph, value) || continue
        count = get(counts, value, 0)
        cycle_tracked = ismutable(value) && !(value isa Union{BigInt,Regex})
        count == (cycle_tracked ? 1 : 2) && continue
        counts[value] = count + 1
        for (item, stored) in graph[value].edges
            stored && push!(pending, item)
        end
    end

    claims = Tuple{AbstractArray,Any}[]
    copied = IdSet{Any}()
    tables = Union{Dict,IdDict,Expr}[]
    for (value, node) in graph
        get(counts, value, 0) > 0 || continue
        value isa Union{BigInt,Regex} &&
            get(counts, value, 0) > 1 &&
            throw(
                ArgumentError("checkpoint configuration repeats a value that serialization copies"),
            )
        value isa Union{Dict,IdDict,Expr} && push!(tables, value)
        node.logical && value isa AbstractArray || continue
        owner = _checkpoint_serialization_owner(value)
        push!(claims, (value, owner))
        if owner isa SubArray && parent(owner) isa Array
            get(counts, owner, 0) > 1 && throw(
                ArgumentError("checkpoint configuration repeats a view that serialization copies"),
            )
            if !(owner in copied)
                push!(copied, owner)
                indices = parentindices(owner)
                for (i, index) in enumerate(indices)
                    index isa AbstractArray && !(index isa AbstractRange) || continue
                    for j = 1:(i-1)
                        other = indices[j]
                        other isa AbstractArray &&
                            !(other isa AbstractRange) &&
                            _particle_arrays_mightalias(index, other) &&
                            throw(
                                ArgumentError(
                                    "checkpoint view indices alias storage that serialization copies independently",
                                ),
                            )
                    end
                    push!(claims, (index, owner))
                end
            end
        end
    end
    # Each reconstructed buffer is a separate allocation, including buffers
    # shared by distinct containers before serialization. Keep capacity opaque.
    for table in tables
        for (i, buffer) in enumerate(_checkpoint_backing_arrays(table))
            require_independent(buffer)
            push!(claims, (buffer, (table, i)))
        end
    end
    for i in eachindex(claims)
        array, owner = claims[i]
        for j = 1:(i-1)
            other, other_owner = claims[j]
            owner === other_owner && continue
            _particle_arrays_mightalias(array, other) && throw(
                ArgumentError(
                    "checkpoint configuration aliases storage that serialization copies independently",
                ),
            )
        end
    end
    return nothing
end

function _validate_checkpoint_configuration(s, st, configuration)
    _checkpoint_pointer_free_bits(typeof(configuration)) && return nothing
    owned = IdDict{Any,Nothing}()
    arrays = AbstractArray[]
    function remember(value)
        # Equal immutable metadata can be shared by value without owning state.
        # Its mutable children are still visited and registered independently.
        (ismutable(value) || value isa AbstractArray) || return nothing
        haskey(owned, value) && return nothing
        owned[value] = nothing
        value isa AbstractArray && push!(arrays, value)
        return nothing
    end
    # Configuration is replaced on restore, whereas these arrays are copied into
    # destination storage. Include the payload to reject retained incoming links.
    for name in keys(s)
        name in (:model, :shape) && continue
        _visit_checkpoint_objects(remember, getproperty(s, name); backing = remember)
    end
    # Capture must also reject links to scratch and Refs omitted from the payload.
    # Read-only grid geometry is independent configuration; FFT plans and scratch
    # belong to the solver and cannot be restored as detached configuration.
    owned[st] = nothing
    for name in fieldnames(typeof(st))
        name in (:g, :model, :shape) && continue
        _visit_checkpoint_objects(remember, getfield(st, name); backing = remember)
    end
    g = st.g
    _visit_checkpoint_objects(
        remember,
        (g.plan, g.iplan, g.cbuf, g.tbuf, g.abuf, g.sbuf);
        backing = remember,
    )
    graph = IdDict{Any,NamedTuple{(:logical, :edges),Tuple{Bool,Vector{Tuple{Any,Bool}}}}}()
    function record(value, logical, edges)
        # Serialization has dedicated value codecs for these types; their private
        # pointers are reconstructed rather than emitted as stored fields.
        value isa Union{Regex,BigInt} && empty!(edges)
        graph[value] = (logical = logical, edges = edges)
        return nothing
    end
    function require_independent(value)
        haskey(owned, value) && throw(
            ArgumentError(
                "checkpoint configuration must not reference evolving simulation storage",
            ),
        )
        if value isa AbstractArray
            for array in arrays
                _particle_arrays_mightalias(value, array) && throw(
                    ArgumentError(
                        "checkpoint configuration must not alias evolving simulation storage",
                    ),
                )
            end
        end
    end
    _visit_checkpoint_objects(require_independent, configuration; record)
    _validate_checkpoint_configuration_graph(configuration, graph, require_independent)
    return nothing
end

function _validate_checkpoint_state(s, st::_HybridCheckpointStepper{D,T}) where {D,T}
    collection = s isa NamedTuple && haskey(s, :particle_layout)
    _validate_checkpoint_container(
        s,
        collection ? _COLLECTION_CHECKPOINT_REQUIRED_FIELDS : _CHECKPOINT_REQUIRED_FIELDS,
    )
    _validate_checkpoint_grid(s, st.g)
    typeof(s.model) === typeof(st.model) ||
        throw(ArgumentError("checkpoint model type $(typeof(s.model)) ≠ $(typeof(st.model))"))
    typeof(s.shape) === typeof(st.shape) ||
        throw(ArgumentError("checkpoint shape type $(typeof(s.shape)) ≠ $(typeof(st.shape))"))
    _validate_evolving_model(s.model, T)
    if collection
        s.particle_layout === :collection ||
            throw(ArgumentError("unsupported checkpoint particle layout $(s.particle_layout)"))
        _validate_checkpoint_time_state(s, T)
        s.species isa AbstractVector ||
            throw(ArgumentError("checkpoint species must be an ordered vector"))
        populations =
            ParticleSet{D,T}[_checkpoint_species(record, Val(D), T) for record in s.species]
        _validate_hybrid_species(populations)
    else
        _validate_checkpoint_scalar_state(s, T)
        _validate_checkpoint_particle_state(s, Val(D))
    end
    _validate_checkpoint_field_state(s, st)
    _validate_checkpoint_integrator_state(s, st)
    _validate_checkpoint_configuration(s, st, (s.model, s.shape))
    return nothing
end

_espic_checkpoint_integrator(::Electrostatic1D) = :espic1d
_espic_checkpoint_integrator(st::ElectrostaticPIC) =
    st.collection === nothing ? :espic : :espic_collection

function _validate_checkpoint_state(s, st::_ElectrostaticCheckpointStepper)
    if s isa NamedTuple && haskey(s, :particle_layout)
        st isa ElectrostaticPIC && st.collection !== nothing || throw(
            ArgumentError(
                "collection checkpoint requires a collection ElectrostaticPIC destination",
            ),
        )
        return _validate_espic_collection_checkpoint(s, st)
    end
    st isa ElectrostaticPIC &&
        st.collection !== nothing &&
        throw(ArgumentError("collection ElectrostaticPIC requires a collection checkpoint"))
    _validate_checkpoint_container(s, _ESPIC_CHECKPOINT_REQUIRED_FIELDS)
    _validate_checkpoint_grid(s, st.g)
    expected = _espic_checkpoint_integrator(st)
    s.integrator === expected ||
        throw(ArgumentError("checkpoint integrator $(s.integrator) ≠ $expected"))
    T = eltype(st.ne)
    _validate_checkpoint_particle_scalars(s, T)
    s.q == -one(T) && s.m == one(T) ||
        throw(ArgumentError("electrostatic checkpoint requires electrons with q = -1 and m = 1"))
    s.n0 isa T || throw(ArgumentError("checkpoint n0 type $(typeof(s.n0)) ≠ $T"))
    isfinite(s.n0) && s.n0 >= zero(T) ||
        throw(ArgumentError("checkpoint n0 must be finite and nonnegative"))
    isequal(s.n0, st.n0) || throw(ArgumentError("checkpoint background density $(s.n0) ≠ $(st.n0)"))
    s.primed isa Bool || throw(ArgumentError("checkpoint primed must be a Bool"))
    s.last_dt isa T || throw(ArgumentError("checkpoint last_dt type $(typeof(s.last_dt)) ≠ $T"))
    isfinite(s.last_dt) && s.last_dt >= zero(T) ||
        throw(ArgumentError("checkpoint last_dt must be finite and nonnegative"))
    (s.primed ? s.last_dt > zero(T) : iszero(s.last_dt)) ||
        throw(ArgumentError("checkpoint last_dt must be positive exactly when primed"))
    _validate_checkpoint_particle_state(s, Val(length(st.g.n)))
    _validate_checkpoint_array(:ne, s.ne, st.ne)
    if st isa Electrostatic1D
        _validate_checkpoint_array(:E, s.E, st.E)
    else
        _validate_checkpoint_field_components(:E, s.E, st.E)
    end
    return nothing
end

function _checkpoint_particle_state(ps::ParticleSet)
    return (x = ps.x, v = ps.v, weight = ps.weight, id = ps.id, tag = ps.tag, q = ps.q, m = ps.m)
end

function _checkpoint_state(st::_HybridCheckpointStepper{D,T}, ps::ParticleSet{D,T}) where {D,T}
    _evolves_pressure(st.model.closure) && _validate_pressure_storage(st, ps)
    return _hybrid_checkpoint_state(st, _checkpoint_particle_state(ps))
end

function _checkpoint_state(
    st::_HybridCheckpointStepper{D,T},
    species::AbstractVector{<:ParticleSet},
) where {D,T}
    species isa AbstractVector{<:ParticleSet{D,T}} ||
        throw(ArgumentError("checkpoint species must match the stepper dimension and precision"))
    _evolves_pressure(st.model.closure) && _validate_pressure_storage(st, species)
    particles = (
        particle_layout = :collection,
        species = [_checkpoint_particle_state(ps) for ps in species],
    )
    return _hybrid_checkpoint_state(st, particles)
end

function _hybrid_checkpoint_state(
    st::_HybridCheckpointStepper{D,T},
    particles::NamedTuple,
) where {D,T}
    common = (;
        integrator = st isa CAMCLStepper ? :camcl : :rk4,
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        model = st.model,
        shape = st.shape,
        particles...,
        B = st.fields.B,
        E = st.fields.E,
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
    )
    if _evolves_pressure(st.model.closure)
        return merge(common, (pe = st.fields.pe,))
    end
    if st isa CAMCLStepper
        return merge(
            common,
            (
                Blead = st.Blead,
                Blag = st.Blag,
                Bout = st.Bout,
                cl_h = st.cl_h[],
                cl_live = st.cl_live[],
            ),
        )
    end
    return common
end

function _checkpoint_state(st::_ElectrostaticCheckpointStepper, ps::ParticleSet{D,T}) where {D,T}
    st isa ElectrostaticPIC && _espic_require_scalar(st)
    return (;
        integrator = _espic_checkpoint_integrator(st),
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        n0 = st.n0,
        _checkpoint_particle_state(ps)...,
        E = st.E,
        ne = st.ne,
        primed = st.primed[],
        last_dt = st.last_dt[],
    )
end

_empic_checkpoint_integrator(::EMPIC1D) = :empic1d
_empic_checkpoint_integrator(::EMPIC) = :empic
_empic_checkpoint_fields(st::EMPIC1D) = (
    Ex = st.Ex,
    Ey = st.Ey,
    Bz = st.Bz,
    Jx = st.Jx,
    Jy = st.Jy,
    rho_n = st.rho_n,
    rho_np1 = st.rho_np1,
)
_empic_checkpoint_fields(st::EMPIC) =
    (E = st.E, B = st.B, J = st.J, rho_n = st.rho_n, rho_np1 = st.rho_np1)

function _checkpoint_state(
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet} = nothing,
)
    return (
        integrator = _empic_checkpoint_integrator(st),
        D = length(st.g.n),
        T = typeof(st.n0),
        ncell = st.g.n,
        L = st.g.L,
        n0 = st.n0,
        c = st.c,
        shape = st.shape,
        relativistic = st.relativistic,
        mobile = st.mobile,
        mi = st.mi,
        n_sub = st.n_sub,
        electrons = _checkpoint_particle_state(electrons),
        ions = ions === nothing ? nothing : _checkpoint_particle_state(ions),
        fields = _empic_checkpoint_fields(st),
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
        last_sub_dt = st.last_sub_dt[],
    )
end

function _validate_checkpoint_state(s, st::_ElectromagneticCheckpointStepper)
    _validate_checkpoint_container(s, _EMPIC_CHECKPOINT_REQUIRED_FIELDS)
    _validate_checkpoint_grid(s, st.g)
    expected = _empic_checkpoint_integrator(st)
    s.integrator === expected ||
        throw(ArgumentError("checkpoint integrator $(s.integrator) ≠ $expected"))
    typeof(s.shape) === typeof(st.shape) ||
        throw(ArgumentError("checkpoint shape type $(typeof(s.shape)) ≠ $(typeof(st.shape))"))
    T = s.T
    for name in (:n0, :c, :mi)
        value = getproperty(s, name)
        value isa T || throw(ArgumentError("checkpoint $name type $(typeof(value)) ≠ $T"))
    end
    _require_finite_nonnegative_real("checkpoint n0", s.n0, T)
    _require_finite_positive_real("checkpoint c", s.c, T)
    _require_finite_positive_real("checkpoint mi", s.mi, T)
    s.relativistic isa Bool || throw(ArgumentError("checkpoint relativistic must be a Bool"))
    s.mobile isa Bool || throw(ArgumentError("checkpoint mobile must be a Bool"))
    s.mobile == st.mobile ||
        throw(ArgumentError("checkpoint ion mobility does not match destination"))
    s.n_sub isa Int && s.n_sub >= 1 ||
        throw(ArgumentError("checkpoint n_sub must be a positive Int"))
    _validate_checkpoint_time_state(s, T)
    s.last_sub_dt isa T ||
        throw(ArgumentError("checkpoint last_sub_dt type $(typeof(s.last_sub_dt)) ≠ $T"))
    isfinite(s.last_sub_dt) && zero(T) <= s.last_sub_dt <= s.last_dt ||
        throw(ArgumentError("checkpoint last_sub_dt must lie between zero and last_dt"))
    (s.step == 0 ? iszero(s.last_sub_dt) : s.last_sub_dt > zero(T)) || throw(
        ArgumentError("checkpoint last_sub_dt must be zero before stepping and positive afterward"),
    )
    dims = Val(s.D)
    electrons = _checkpoint_species(s.electrons, dims, T)
    _require_empic_electrons(electrons)
    s.relativistic && _require_subluminal_velocities(electrons, s.c, "checkpoint electrons")
    if s.mobile
        s.ions === nothing && throw(ArgumentError("mobile checkpoint requires ion particles"))
        ions = _checkpoint_species(s.ions, dims, T)
        _require_empic_ions(ions, s.mi)
        s.relativistic && _require_subluminal_velocities(ions, s.c, "checkpoint ions")
        _particle_sets_mightalias(electrons, ions) &&
            throw(ArgumentError("checkpoint electron and ion storage must be independent"))
    else
        s.ions === nothing ||
            throw(ArgumentError("immobile checkpoint cannot contain ion particles"))
    end
    _validate_checkpoint_fields(s.fields, _empic_checkpoint_fields(st))
    _validate_checkpoint_configuration(s, st, s.shape)
    return nothing
end

function _hall_checkpoint_fields(st::HallMHDState)
    f = st.fields
    return (
        n = f.n,
        ui = f.ui,
        B = f.B,
        E = f.E,
        J = f.J,
        pe = f.pe,
        gradp = f.gradp,
        ninv = f.ninv,
        lapJ = f.lapJ,
        pforce = f.pforce,
        cgl_dp = f.cgl_dp,
    )
end

function _checkpoint_state(st::HallMHDState{D,T}) where {D,T}
    return (
        integrator = :hall_mhd,
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        model = st.model,
        fields = _hall_checkpoint_fields(st),
        # Ohm's law leaves inactive lapJ components intact when ηH is zero;
        # a later Hall step copies them from stage into the public fields.
        stage_lapJ = st.stage.lapJ,
        floor_count = st.fields.floor_count[],
        time = st.time[],
        step = st.step[],
    )
end

function _validate_checkpoint_state(s, st::HallMHDState{D,T}) where {D,T}
    _validate_checkpoint_container(s, _HALL_CHECKPOINT_REQUIRED_FIELDS)
    _validate_checkpoint_grid(s, st.g)
    s.integrator === :hall_mhd ||
        throw(ArgumentError("checkpoint integrator $(s.integrator) ≠ hall_mhd"))
    typeof(s.model) === typeof(st.model) ||
        throw(ArgumentError("checkpoint model type $(typeof(s.model)) ≠ $(typeof(st.model))"))
    # Reuse the model's numeric contract without replacing its closure object.
    HallMHDModel(
        s.model.closure;
        Ti = s.model.Ti,
        η = s.model.η,
        ηH = s.model.ηH,
        nfloor = s.model.nfloor,
    )
    _validate_checkpoint_clock_state(s, T)
    s.floor_count isa Int && 0 <= s.floor_count <= length(st.fields.n) ||
        throw(ArgumentError("checkpoint floor_count must be an Int within the grid cell count"))
    _validate_checkpoint_fields(s.fields, _hall_checkpoint_fields(st))
    _validate_checkpoint_field_components(:stage_lapJ, s.stage_lapJ, st.stage.lapJ)
    _validate_checkpoint_configuration(s, st, s.model)
    return nothing
end

function _restore_checkpoint_state!(st::HallMHDState, s)
    _validate_checkpoint_state(s, st)
    st.model = s.model
    _restore_checkpoint_fields!(_hall_checkpoint_fields(st), s.fields)
    foreach(copyto!, st.stage.lapJ, s.stage_lapJ)
    st.fields.floor_count[] = s.floor_count
    st.time[] = s.time
    st.step[] = s.step
    return st
end

function _validate_checkpoint_destination(ps::ParticleSet, N::Int)
    _particle_components_mightalias(ps) &&
        throw(ArgumentError("checkpoint destination particle components must not alias"))
    return _require_particle_resize_support(ps, N, "checkpoint")
end

function _restore_checkpoint_particles!(ps::ParticleSet{D}, s) where {D}
    N = length(s.weight)
    for d = 1:D
        length(ps.x[d]) == N || resize!(ps.x[d], N)
        copyto!(ps.x[d], s.x[d])
    end
    for c = 1:3
        length(ps.v[c]) == N || resize!(ps.v[c], N)
        copyto!(ps.v[c], s.v[c])
    end
    length(ps.weight) == N || resize!(ps.weight, N)
    copyto!(ps.weight, s.weight)
    length(ps.id) == N || resize!(ps.id, N)
    copyto!(ps.id, s.id)
    length(ps.tag) == N || resize!(ps.tag, N)
    copyto!(ps.tag, s.tag)
    ps.q = s.q
    ps.m = s.m
    return ps
end

function _validate_hybrid_checkpoint_restore(
    st::_HybridCheckpointStepper{D,T},
    ps::ParticleSet{D,T},
    s,
) where {D,T}
    _validate_checkpoint_state(s, st)
    haskey(s, :particle_layout) &&
        throw(ArgumentError("collection checkpoint requires a species collection destination"))
    N = length(s.weight)
    _validate_checkpoint_destination(ps, N)
    _evolves_pressure(st.model.closure) && _validate_pressure_storage(st, ps)
    return N
end

function _validate_hybrid_checkpoint_restore(
    st::_HybridCheckpointStepper{D,T},
    species::AbstractVector{<:ParticleSet},
    s,
) where {D,T}
    _validate_checkpoint_state(s, st)
    haskey(s, :particle_layout) ||
        throw(ArgumentError("species collection destination requires a collection checkpoint"))
    species isa AbstractVector{<:ParticleSet{D,T}} ||
        throw(ArgumentError("destination species must match the stepper dimension and precision"))
    length(species) == length(s.species) ||
        throw(ArgumentError("checkpoint population count does not match destination"))
    _require_independent_hybrid_species(species)
    total = 0
    for (ps, record) in zip(species, s.species)
        count = length(record.weight)
        _validate_checkpoint_destination(ps, count)
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            _check_particle_vector_axes(:checkpoint_destination, a, length(a))
        end
        total += count # Saved populations have already passed the total-count overflow check.
    end
    _evolves_pressure(st.model.closure) && _validate_pressure_storage(st, species)
    return total
end

function _restore_validated_hybrid_checkpoint!(st, ps::ParticleSet, s, count::Int)
    _resize_hybrid_particle_workspaces!(st, count)
    _restore_checkpoint_particles!(ps, s)
    return _restore_hybrid_checkpoint_fields!(st, s)
end

function _restore_validated_hybrid_checkpoint!(st, species::AbstractVector, s, count::Int)
    _resize_hybrid_particle_workspaces!(st, count)
    for (ps, record) in zip(species, s.species)
        _restore_checkpoint_particles!(ps, record)
    end
    return _restore_hybrid_checkpoint_fields!(st, s)
end

function _restore_checkpoint_state!(
    st::_HybridCheckpointStepper{D,T},
    ps::ParticleSet{D,T},
    s,
) where {D,T}
    count = _validate_hybrid_checkpoint_restore(st, ps, s)
    return _restore_validated_hybrid_checkpoint!(st, ps, s, count)
end

function _restore_checkpoint_state!(
    st::_HybridCheckpointStepper,
    species::AbstractVector{<:ParticleSet},
    s,
)
    count = _validate_hybrid_checkpoint_restore(st, species, s)
    return _restore_validated_hybrid_checkpoint!(st, species, s, count)
end

function _restore_hybrid_checkpoint_fields!(st::_HybridCheckpointStepper, s)
    st.model = s.model
    st.shape = s.shape
    for c = 1:3
        copyto!(st.fields.B[c], s.B[c])
        copyto!(st.fields.E[c], s.E[c])
    end
    _evolves_pressure(st.model.closure) && copyto!(st.fields.pe, s.pe)
    st.time[] = s.time
    st.step[] = s.step
    st.last_dt[] = s.last_dt
    if st isa CAMCLStepper
        for (destination, saved) in ((st.Blead, s.Blead), (st.Blag, s.Blag), (st.Bout, s.Bout)),
            c = 1:3

            copyto!(destination[c], saved[c])
        end
        st.cl_h[] = s.cl_h
        st.cl_live[] = s.cl_live
    end
    return st
end

function _resize_empic_checkpoint_workspaces!(st::EMPIC1D, ne, ni)
    _ensure_electron_buffers!(st, ne)
    _ensure_ion_buffers!(st, ni)
    return nothing
end

function _resize_empic_checkpoint_workspaces!(st::EMPIC, ne, ni)
    _ensure_empic_electron_buffers!(st, ne)
    _ensure_empic_ion_buffers!(st, ni)
    return nothing
end

_restore_checkpoint_state!(st::_ElectromagneticCheckpointStepper, electrons::ParticleSet, s) =
    _restore_checkpoint_state!(st, electrons, nothing, s)

function _restore_checkpoint_state!(
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet},
    s,
)
    _validate_checkpoint_state(s, st)
    D, T = s.D, s.T
    electrons isa ParticleSet{D,T} || throw(
        ArgumentError("checkpoint electrons do not match destination dimension and precision"),
    )
    if ions === nothing
        s.mobile && throw(
            ArgumentError(
                "mobile checkpoint requires destination ions of matching dimension and precision",
            ),
        )
    else
        s.mobile || throw(ArgumentError("immobile checkpoint does not restore ion particles"))
        ions isa ParticleSet{D,T} || throw(
            ArgumentError(
                "mobile checkpoint requires destination ions of matching dimension and precision",
            ),
        )
        _particle_sets_mightalias(electrons, ions) &&
            throw(ArgumentError("destination electron and ion storage must be independent"))
    end
    ne = length(s.electrons.weight)
    ni = ions === nothing ? 0 : length(s.ions.weight)
    _validate_checkpoint_destination(electrons, ne)
    ions === nothing || _validate_checkpoint_destination(ions, ni)
    _resize_empic_checkpoint_workspaces!(st, ne, ni)
    for name in (:n0, :c, :shape, :relativistic, :mi, :n_sub)
        setproperty!(st, name, getproperty(s, name))
    end
    _restore_checkpoint_particles!(electrons, s.electrons)
    ions === nothing || _restore_checkpoint_particles!(ions, s.ions)
    _restore_checkpoint_fields!(_empic_checkpoint_fields(st), s.fields)
    st.time[] = s.time
    st.step[] = s.step
    st.last_dt[] = s.last_dt
    st.last_sub_dt[] = s.last_sub_dt
    return st
end

function _restore_checkpoint_state!(
    st::_ElectrostaticCheckpointStepper,
    ps::ParticleSet{D,T},
    s,
) where {D,T}
    _validate_checkpoint_state(s, st)
    s.D == D && s.T === T || throw(
        ArgumentError("checkpoint particles do not match destination dimension and precision"),
    )
    N = length(s.weight)
    _validate_checkpoint_destination(ps, N)
    _resize_espic_particle_workspaces!(st, N)
    _restore_checkpoint_particles!(ps, s)
    copyto!(st.ne, s.ne)
    if st isa Electrostatic1D
        copyto!(st.E, s.E)
    else
        for c = 1:3
            copyto!(st.E[c], s.E[c])
        end
    end
    st.primed[] = s.primed
    st.last_dt[] = s.last_dt
    return st
end

"""
    save_checkpoint(path, stepper, ps)
    save_checkpoint(path, stepper, species::AbstractVector)
    save_checkpoint(path, empic, electrons, ions)
    save_checkpoint(path, hall::HallMHDState)

Serialize the restartable state (model, shape, particles, B, carried E, time,
step, and previous timestep for the staggered particle velocities) to `path`.
Supports [`HybridStepper`](@ref) and [`CAMCLStepper`](@ref); CAM-CL checkpoints
also preserve the persistent cyclic-leapfrog field copies and subcycle phase.
Hybrid species collections are stored in order, with each population's charge,
mass, phase space, and provenance. Populations require independent storage and
the positive-ion validation used by collection stepping.
For [`Electrostatic1D`](@ref) and [`ElectrostaticPIC`](@ref), save the electron
particles, electric field, number density, background density, and priming/timestep
state. Collection `ElectrostaticPIC` checkpoints also store ordered arbitrary
species, neutrality policy and charge tolerance. Scalar payloads retain their
original schema; collection schema 1 is distinct. ES collections use CIC and one
shared velocity half level. Background density, policy and tolerance must match
the destination; saved physical species values replace the destination values.
[`EMPIC1D`](@ref) and [`EMPIC`](@ref) preserve their physics configuration,
electric/magnetic fields, current and charge snapshots, and full/fast timestep
history. Supply the ion `ParticleSet` exactly when `mobile=true`; electron and
ion storage must be independent. Particle workspaces are resized on restore;
scratch values are recomputed during stepping.
Mutable model/shape parameters must be independent of particle and evolving
solver storage, including FFT plans, scratch arrays, and clock/phase references. Links to
evolving storage are rejected before serialization can detach them. Array checks
use Base's alias and parent interfaces; custom arrays must implement those
interfaces consistently with their logical storage. Access outside a declared
view is not inspected.
[`HallMHDState`](@ref) stores its model, fluid fields, public derived fields,
and time/step counters without particle arguments. The retained stage `lapJ`
values are saved because inactive hyperresistivity can carry them into the public
fields. Other RK workspaces are rebuilt during stepping; saved derived fields
are restored without reevaluating Ohm's law.

Publish by same-directory rename after serialization, flush and close succeed;
earlier failures preserve the previous file. Cleanup failures accompany the
original error. Existing file rwx permissions are retained; new files use 0600.
Symlinks are followed to their referent, whose inode is replaced; other hardlinks
keep the old contents. Requires directory write access and no concurrent directory
mutation. This does not guarantee durability after an operating-system crash.
"""
function save_checkpoint(path::AbstractString, st::_CheckpointStepper, ps::ParticleSet)
    return _save_checkpoint(path, st, _checkpoint_state(st, ps))
end

function save_checkpoint(path::AbstractString, st::HallMHDState)
    return _save_checkpoint(path, st, _checkpoint_state(st))
end

function save_checkpoint(
    path::AbstractString,
    st::Union{_HybridCheckpointStepper,ElectrostaticPIC},
    species::AbstractVector{<:ParticleSet},
)
    return _save_checkpoint(path, st, _checkpoint_state(st, species))
end

function save_checkpoint(
    path::AbstractString,
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet},
)
    return _save_checkpoint(path, st, _checkpoint_state(st, electrons, ions))
end

function _save_checkpoint(path::AbstractString, st::Union{_CheckpointStepper,HallMHDState}, state)
    _validate_checkpoint_state(state, st)
    _atomic_serialize(path, state)
    return path
end

"""
    load_checkpoint!(stepper, ps, path)
    load_checkpoint!(stepper, species::AbstractVector, path)
    load_checkpoint!(empic, electrons, ions, path)
    load_checkpoint!(hall::HallMHDState, path)

Restore a checkpoint into a `stepper`/`ps` of matching dimension, precision,
grid, and integrator (particle arrays are resized to match).
Destination arrays that do not support resizing must already have the saved
particle count; writable views are restored in place. Particle components must
have independent storage so restoring one cannot overwrite another.
Hybrid destinations require matching concrete model/shape types; the checkpoint
replaces their parameters, including custom mutable parameters, so restored
evolution uses the saved physics.
For `EvolvingElectrons`, stored pressure is required and restored verbatim. Its
shape, precision, finite nonnegative values, and ideal model parameters are
validated before mutation. Algebraic closures retain pressure-free payloads.
Collection destinations require the same number of ordered populations and
matching dimension/precision. Each population keeps its existing storage while
saved charge, mass, provenance and phase space replace its values. Every
population is validated before any is changed; scalar and collection layouts
are distinct and cannot be interchanged on restore.
Electrostatic destinations must have the same background density `n0`, which
is immutable; their saved fields, density, and leapfrog phase are restored directly.
Electromagnetic destinations require the same ion mobility and concrete shape
type; saved physics parameters replace the destination values. Supply independent
electron and ion destinations for mobile runs. Both species are validated before
either is changed.
Hall-MHD destinations require the same concrete model/closure type; saved model
parameters replace destination values while existing field arrays remain in place.
Checkpoints missing fields required by their model family are rejected.
After loading, `step!`, `step_camcl!`, `step_espic!`, or `step_empic!` continues bitwise-identically to the
original run with the same stepping inputs. Hall-MHD continues with
`step_hall_mhd!`, including the same `project_b` inputs.
"""
function load_checkpoint!(st::_CheckpointStepper, ps::ParticleSet, path::AbstractString)
    return _restore_checkpoint_state!(st, ps, deserialize(path))
end

function load_checkpoint!(st::HallMHDState, path::AbstractString)
    return _restore_checkpoint_state!(st, deserialize(path))
end

function load_checkpoint!(
    st::Union{_HybridCheckpointStepper,ElectrostaticPIC},
    species::AbstractVector{<:ParticleSet},
    path::AbstractString,
)
    return _restore_checkpoint_state!(st, species, deserialize(path))
end

function load_checkpoint!(
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet},
    path::AbstractString,
)
    return _restore_checkpoint_state!(st, electrons, ions, deserialize(path))
end


const _ESPIC_COLLECTION_CHECKPOINT_FIELDS = (
    :integrator,
    :D,
    :T,
    :ncell,
    :L,
    :n0,
    :particle_layout,
    :schema_version,
    :neutrality,
    :charge_rtol,
    :species,
    :E,
    :ne,
    :primed,
    :last_dt,
    :initialized,
)

function _validate_espic_collection_checkpoint(s, st::ElectrostaticPIC{D,T}) where {D,T}
    _validate_checkpoint_container(s, _ESPIC_COLLECTION_CHECKPOINT_FIELDS)
    _validate_checkpoint_grid(s, st.g)
    s.integrator === :espic_collection &&
        s.particle_layout === :collection &&
        s.schema_version === 1 ||
        throw(ArgumentError("unsupported electrostatic collection checkpoint schema"))
    work = _espic_collection(st)
    s.n0 isa T && isequal(s.n0, st.n0) ||
        throw(ArgumentError("checkpoint background_charge does not match destination"))
    s.neutrality === work.neutrality ||
        throw(ArgumentError("checkpoint neutrality policy does not match destination"))
    s.charge_rtol isa T && isequal(s.charge_rtol, work.charge_rtol) ||
        throw(ArgumentError("checkpoint charge_rtol does not match destination"))
    s.primed isa Bool && s.initialized isa Bool ||
        throw(ArgumentError("checkpoint phase flags must be Bool"))
    s.last_dt isa T &&
        isfinite(s.last_dt) &&
        (s.primed ? s.last_dt > zero(T) : iszero(s.last_dt)) ||
        throw(ArgumentError("checkpoint last_dt must be positive exactly when primed"))
    s.primed && !s.initialized && throw(ArgumentError("primed checkpoint must be initialized"))
    s.species isa AbstractVector ||
        throw(ArgumentError("checkpoint species must be an ordered vector"))
    populations = ParticleSet{D,T}[_checkpoint_species(record, Val(D), T) for record in s.species]
    _espic_validate_species(st, populations)
    _validate_checkpoint_array(:ne, s.ne, st.ne)
    _validate_checkpoint_field_components(:E, s.E, st.E)
    for a in (s.ne, s.E...)
        all(isfinite, a) || throw(ArgumentError("checkpoint electrostatic fields must be finite"))
        for ps in populations
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("checkpoint particles must not alias fields"))
        end
    end
    return nothing
end

function _checkpoint_state(
    st::ElectrostaticPIC{D,T},
    species::AbstractVector{<:ParticleSet},
) where {D,T}
    work = _espic_collection(st)
    _espic_validate_species(st, species)
    if work.initialized[]
        _espic_require_bound_species(work, species)
    end
    return (
        integrator = :espic_collection,
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        n0 = st.n0,
        particle_layout = :collection,
        schema_version = 1,
        neutrality = work.neutrality,
        charge_rtol = work.charge_rtol,
        species = [_checkpoint_particle_state(ps) for ps in species],
        E = st.E,
        ne = st.ne,
        primed = st.primed[],
        last_dt = st.last_dt[],
        initialized = work.initialized[],
    )
end

function _restore_checkpoint_state!(
    st::ElectrostaticPIC{D,T},
    species::AbstractVector{<:ParticleSet},
    s,
) where {D,T}
    work = _espic_collection(st)
    _validate_checkpoint_state(s, st)
    species isa AbstractVector{<:ParticleSet{D,T}} ||
        throw(ArgumentError("destination species must match grid dimension and precision"))
    length(species) == length(s.species) ||
        throw(ArgumentError("checkpoint population count does not match destination"))
    # Destination particle values are replaced; validate storage without requiring
    # the old destination charge or phase to match the saved physical state.
    _espic_validate_owner_types(work, species)
    _espic_validate_storage(st, species)
    for (ps, record) in zip(species, s.species)
        _validate_checkpoint_destination(ps, length(record.weight))
    end
    for (ps, record) in zip(species, s.species)
        _restore_checkpoint_particles!(ps, record)
    end
    copyto!(st.ne, s.ne)
    foreach(copyto!, st.E, s.E)
    empty!(work.owners)
    append!(work.owners, species)
    work.initialized[] = s.initialized
    st.primed[] = s.primed
    st.last_dt[] = s.last_dt
    return st
end
