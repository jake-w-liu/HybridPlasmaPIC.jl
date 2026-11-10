# checkpoint.jl — checkpoint/restart of the deterministic simulation state.
#
# The deterministic steppers use no RNG, so their restart state contains model
# and shape, particle phase space + provenance, canonical B,
# carried E, evolving electron pressure when present, the scalar transport
# carried ion force with its policy version, time/step counters, and
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
# CollisionalHybridStepper extends these owners in collisional_checkpoint.jl
# with full RNG state, a fixed-step source ledger and every public hybrid field.
# ponytail: stdlib Serialization (no extra dep). Swap for HDF5 when cross-version
# or cross-language portability is needed.

using Serialization

const _HybridCheckpointStepper{D,T} = Union{HybridStepper{D,T},CAMCLStepper{D,T}}
const _ElectrostaticCheckpointStepper = Union{Electrostatic1D,ElectrostaticPIC}
const _ElectromagneticCheckpointStepper = Union{EMPIC1D,EMPIC}
const _ClosedYeeCheckpointStepper = Union{
    YeeLongitudinalPIC,
    YeeElectromagneticPIC,
    YeeImpactLongitudinalPIC,
    YeeAbsorbingLongitudinalPIC,
}
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
    if !all(k -> haskey(s, k), required)
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
    _validated_charge_mass(s.q, s.m, T)
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

function _restore_checkpoint_fields!(
    fields::NamedTuple{names},
    saved::NamedTuple{names},
) where {names}
    for (destination, source) in zip(values(fields), values(saved))
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
    _transports_electrons(s.model.closure) && _validate_checkpoint_transport_state(s, st)
    return nothing
end

# Version of the TransportElectrons payload meaning (carried ion force E_force
# with constant dissipative transport). It is required; old payloads carry no force.
const _TRANSPORT_POLICY_VERSION = 1
# Original pressure/density quotient preparation and projected conduction.
const _CONDUCTION_POLICY_VERSION = 2
const _FIELD_ALIGNED_CONDUCTION_POLICY_VERSION = 3

_validate_checkpoint_conduction_state(s, ::ElectronClosure) = nothing
function _validate_checkpoint_conduction_state(s, ::ConductingElectrons)
    _validate_checkpoint_container(s, (:conduction_policy_version,))
    s.conduction_policy_version === _CONDUCTION_POLICY_VERSION ||
        throw(ArgumentError("unsupported electron conduction checkpoint policy"))
    return nothing
end

function _validate_checkpoint_conduction_state(s, closure::FieldAlignedConductingElectrons)
    _validate_checkpoint_container(
        s,
        (
            :conduction_policy_version,
            :conduction_operator,
            :conduction_direction,
            :conduction_preparation,
        ),
    )
    s.conduction_policy_version === _FIELD_ALIGNED_CONDUCTION_POLICY_VERSION &&
        s.conduction_operator === :odd_tensor_diagonal_nyquist &&
        s.conduction_direction === :original_Q_nonzero_projector &&
        s.conduction_preparation === :two_term_axis_null_gradient_flux ||
        throw(ArgumentError("unsupported field-aligned conduction checkpoint policy"))
    # Admission follows the incoming model, including zero-to-active restore.
    # The derived direction/gradient workspace is recomputed, never loaded.
    T = eltype(s.pe)
    _electron_conductivity(closure, T)
    for array in (s.B..., s.E...)
        _validate_pressure_finite(array, "field-aligned checkpoint field")
    end
    _validate_field_aligned_direction(closure, s.B, T)
    return nothing
end

function _validate_checkpoint_transport_state(s, st::_HybridCheckpointStepper{D,T}) where {D,T}
    _validate_checkpoint_container(s, (:transport_policy_version, :Eforce))
    _validate_checkpoint_conduction_state(s, s.model.closure)
    s.transport_policy_version === _TRANSPORT_POLICY_VERSION || throw(
        ArgumentError(
            "checkpoint transport_policy_version $(repr(s.transport_policy_version)) ≠ $_TRANSPORT_POLICY_VERSION",
        ),
    )
    s.Eforce isa NTuple{3,Array{T,D}} ||
        throw(ArgumentError("checkpoint Eforce must be stored as three Array{$T,$D}"))
    for c = 1:3
        force = s.Eforce[c]
        _require_grid_array(:Eforce, c, force, st.g)
        _validate_pressure_finite(force, "checkpoint transport ion force")
        # A valid writer stores independent arrays; shared payload storage
        # indicates a malformed or hand-edited state.
        aliased = any(k -> _particle_arrays_mightalias(force, s.Eforce[k]), 1:(c-1))
        for other in (s.E..., s.B..., s.pe)
            aliased |= _particle_arrays_mightalias(force, other)
        end
        aliased && throw(ArgumentError("checkpoint Eforce must not alias other saved fields"))
    end
    if !haskey(s, :particle_layout)
        isone(s.q) || throw(ArgumentError("transport scalar checkpoint populations require q = 1"))
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

# The packed chunks holding the elements a parent chain reads, if a BitArray
# stores them. A parent cycle reports none, keeping the generic alias check.
function _checkpoint_packed_chunks(array::AbstractArray)
    seen = IdSet{Any}()
    while !(array in seen)
        array isa BitArray && return array.chunks
        push!(seen, array)
        next = parent(array)
        next === array && return nothing
        array = next
    end
    return nothing
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

# FFT plans and scratch belong to the spectral solver and cannot be restored as
# detached configuration; a plain Yee geometry grid owns no solver scratch.
_checkpoint_owned_scratch(g::FourierGrid) = (g.plan, g.iplan, g.cbuf, g.tbuf, g.abuf, g.sbuf)
_checkpoint_owned_scratch(::YeeGrid) = ()
function _register_checkpoint_configuration_owned!(remember, remember_backing, owned, s, st)
    # Configuration is replaced on restore, whereas these arrays are copied into
    # destination storage. Include the payload to reject retained incoming links.
    for name in keys(s)
        name in (:model, :shape) && continue
        _visit_checkpoint_objects(remember, getproperty(s, name); backing = remember_backing)
    end
    # Capture must also reject links to scratch and Refs omitted from the payload.
    # Read-only grid geometry is independent configuration; FFT plans and scratch
    # belong to the solver and cannot be restored as detached configuration.
    owned[st] = true
    for name in fieldnames(typeof(st))
        name in (:g, :model, :shape) && continue
        _visit_checkpoint_objects(remember, getfield(st, name); backing = remember_backing)
    end
    _visit_checkpoint_objects(remember, _checkpoint_owned_scratch(st.g); backing = remember_backing)
    return nothing
end

function _validate_checkpoint_configuration(s, st, configuration)
    return _validate_checkpoint_configuration_registered(
        configuration,
    ) do remember, remember_backing, owned
        _register_checkpoint_configuration_owned!(remember, remember_backing, owned, s, st)
    end
end

function _validate_checkpoint_configuration_registered(register_owned!, configuration)
    _checkpoint_pointer_free_bits(typeof(configuration)) && return nothing
    # Whether any visit reached an object logically, not only as backing storage.
    owned = IdDict{Any,Bool}()
    arrays = AbstractArray[]
    function remember(value, logical = true)
        # Equal immutable metadata can be shared by value without owning state.
        # Its mutable children are still visited and registered independently.
        (ismutable(value) || value isa AbstractArray) || return nothing
        if haskey(owned, value)
            owned[value] |= logical
        else
            owned[value] = logical
            value isa AbstractArray && push!(arrays, value)
        end
        return nothing
    end
    remember_backing(buffer) = remember(buffer, false)
    register_owned!(remember, remember_backing, owned)
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
            # A value read through a BitArray shares that array's chunks. They link
            # it to state only if state reached them logically; packed backing
            # alone defers to the owned views' logical extents (Julia >= 1.13.3
            # reports every view of a BitArray as sharing its chunks).
            chunks = _checkpoint_packed_chunks(value)
            for array in arrays
                linked = array === chunks ? owned[array] : _particle_arrays_mightalias(value, array)
                linked && throw(
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
    _validate_hybrid_clock_phase(s.time, s.step, s.last_dt)
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

# Stored scalar pressure; transport also stores its policy
# version and the carried ion force, which cannot be rebuilt from half-time
# particle moments.
# Only HybridStepper owns the force; other steppers reject the closure.
_pressure_checkpoint_state(st, ::ElectronClosure) = (pe = st.fields.pe,)
_conduction_checkpoint_state(::ElectronClosure) = NamedTuple()
_conduction_checkpoint_state(::ConductingElectrons) =
    (conduction_policy_version = _CONDUCTION_POLICY_VERSION,)
_conduction_checkpoint_state(::FieldAlignedConductingElectrons) = (
    conduction_policy_version = _FIELD_ALIGNED_CONDUCTION_POLICY_VERSION,
    conduction_operator = :odd_tensor_diagonal_nyquist,
    conduction_direction = :original_Q_nonzero_projector,
    conduction_preparation = :two_term_axis_null_gradient_flux,
)
_pressure_checkpoint_state(st::HybridStepper, closure::_ScalarTransportClosure) = merge(
    (
        pe = st.fields.pe,
        transport_policy_version = _TRANSPORT_POLICY_VERSION,
        Eforce = st.pressure_work.Eforce,
    ),
    _conduction_checkpoint_state(closure),
)

_restore_transport_checkpoint_force!(st, s, ::ElectronClosure) = nothing
function _restore_transport_checkpoint_force!(st::HybridStepper, s, ::_ScalarTransportClosure)
    foreach(copyto!, st.pressure_work.Eforce, s.Eforce)
    return nothing
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
        return merge(common, _pressure_checkpoint_state(st, st.model.closure))
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
_empic_checkpoint_fields(st::EMPIC1D) = _empic_physical_fields(st)
_empic_checkpoint_fields(st::EMPIC) = _empic_physical_fields(st)

function _checkpoint_state(
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet} = nothing,
)
    _require_scalar_empic(st)
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
    if s isa NamedTuple && haskey(s, :particle_layout)
        st.collection !== nothing || throw(
            ArgumentError(
                "collection checkpoint requires an electromagnetic collection destination",
            ),
        )
        return _validate_empic_collection_checkpoint(s, st)
    end
    st.collection !== nothing &&
        throw(ArgumentError("electromagnetic collection requires a collection checkpoint"))
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
    _validate_empic_clock_phase(s.time, s.step, s.last_dt, s.last_sub_dt)
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

_hall_checkpoint_fields(st::HallMHDState) = _hybrid_field_state(st.fields)

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
    hall_fields = _hall_checkpoint_fields(st)
    # Saved payloads may carry extra externally edited keys; project first.
    _restore_checkpoint_fields!(
        hall_fields,
        NamedTuple{keys(hall_fields)}(s.fields),
    )
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
    _restore_transport_checkpoint_force!(st, s, st.model.closure)
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
    _require_scalar_empic(st)
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
    empic_fields = _empic_checkpoint_fields(st)
    # Saved payloads may carry extra externally edited keys; project first.
    _restore_checkpoint_fields!(
        empic_fields,
        NamedTuple{keys(empic_fields)}(s.fields),
    )
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
    save_checkpoint(path, model::MaxwellVacuum)

Serialize the restartable state (model, shape, particles, B, carried E, time,
step, and previous timestep for the staggered particle velocities) to `path`.
Supports [`HybridStepper`](@ref) and [`CAMCLStepper`](@ref); CAM-CL checkpoints
also preserve the persistent cyclic-leapfrog field copies and subcycle phase.
Scalar Hybrid/CAM-CL save and restore require finite signed charge and positive
finite mass through the stepping factor owner, without a quotient restriction.
Other legacy snapshot field/weight value policies are unchanged.
Hybrid species collections are stored in order, with each population's charge,
mass, phase space, and provenance. Populations require independent storage and
the positive-ion validation used by collection stepping.
[`CollisionalHybridStepper`](@ref) saves only initialized completed boundaries,
including its full advanced RNG, source ledger, fixed configuration and all
public hybrid fields. Its versioned payload requires the same Julia runtime and
platform; the RNG seed alone does not describe a stochastic continuation.
[`CollisionalEMPICStepper`](@ref) uses a distinct schema-1 payload containing its
generic EMPIC collection state, fixed field-then-collision policy and timestep,
full RNG and completed source ledger. It saves only initialized ready boundaries;
its collision positions are at the new integer time and velocities remain at the
carried fast half level. No second field or current snapshot is constructed.
Both collision `epsilon0` and `weight_to_number` are explicit, positive and equal
at the stored precision to match the unit Maxwell source coefficient.
[`CollisionalElectrostaticPIC`](@ref) retains its distinct schema-1 ES collection
state, field-then-collision policy, fixed normalization, full RNG and sole wrapper
source clock/count. It saves initialized ready boundaries; count zero is unprimed,
and advanced state retains the common carried half-time velocities. Collision
`epsilon0` and `weight_to_number` are explicit, positive and equal at stored precision.
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
scratch values are recomputed during stepping. `EMPIC` collection schema 1 stores
ordered arbitrary species, per-population mobility, neutrality/background policy
and the common fast velocity phase alongside these fields. Opted-in generic
species schedules instead use schema2/policy1 with `species_stride` and
`last_species_dt`; loading schema2 requires an opted-in destination. Fully
validated generic schema1 without reserved schedule keys migrates only for an
opted-in destination, using all-fast rows and saved last_sub_dt history.
`EMPIC1D` collections use distinct schema 2 `:empic1d_collection` with required
`polarization=:ex_ey_bz`, integer `schedule_policy=1`, ordered `species_stride`
and `last_species_dt`, and physical `Ex/Ey/Bz/Jx/Jy/rho_n/rho_np1` fields.
Only a fully validated specialized schema-1 payload without schedule fields
migrates to all-fast strides; advanced mobile history uses its saved fast interval.
Scalar, generic collection and specialized collection payloads cannot be
interchanged. Collection policy, background, tolerance, mobility and shape type
must match the destination. Saved relativistic mode, c, n_sub, species strides and
physical species values replace destination values. The saved mode and c govern
particle-domain validation, including the original planar EMPIC1D velocities,
before any restore publication. Every incoming leaf is detached and revalidated,
and destination capacity is reserved before physical values are copied in place.
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

`MaxwellVacuum` uses a distinct field-only schema with synchronized E/B, exact
Yee grid/boundary/speed configuration and completed clock. No particles, source
waveform or preceding timestep is invented. Its matching destination arrays and
Refs restore in place after detached incoming validation; no field projection.
`YeeConductiveLayer` has a separate scalar-conductivity schema containing that
canonical vacuum payload, fixed staggered profile/timestep and exact loss/work/
defect ledgers. Derived caches/staging are omitted. Incoming budget identity and
fields detach/validate before publication; PML memory and particles are not supplied.
`YeePML` retains its separate directional ADE policy, every causal memory and
signed exchange/field-defect ledger. Completed count-zero memories must all be
stored zero; advanced memories restore verbatim without reconstructing history.

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
    st::Union{_HybridCheckpointStepper,ElectrostaticPIC,EMPIC,EMPIC1D},
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

function _save_checkpoint(
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
        CUDAHybridStepper,
        SlabGyrokineticPIC,
        OpenHallHybrid,
        ShearedPeriodicHybridPIC,
        ToroidalAnnulusHybridPIC,
        StaticCompositeElectrostaticPIC,
        DynamicCompositeElectrostaticPIC,
        AdaptiveCompositeElectrostaticPIC,
        _ClosedYeeCheckpointStepper,
    },
    state,
)
    _validate_checkpoint_state(state, st)
    _atomic_serialize(path, state)
    return path
end

"""
    load_checkpoint!(stepper, ps, path)
    load_checkpoint!(stepper, species::AbstractVector, path)
    load_checkpoint!(empic, electrons, ions, path)
    load_checkpoint!(hall::HallMHDState, path)
    load_checkpoint!(model::MaxwellVacuum, path)

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
`TransportElectrons` payloads also require `transport_policy_version = 1` and
the carried ion force `Eforce` (three independent finite grid `Array`s); the
force is restored verbatim, never rebuilt from half-time particle moments.
`ConductingElectrons` additionally requires integer `conduction_policy_version = 2`,
including when κ is zero; policy-1 payloads are rejected before mutation. The
saved model replaces γ/κ and derived quotient/transform scratch is not persisted.
`FieldAlignedConductingElectrons` requires integer policy 3 and its exact tensor,
Nyquist, direction and preparation tags. Incoming B is validated against the
saved coefficients, including zero-to-active restore. Both original coefficients
and gamma replace matching-concrete-type model values; derived norm and gradient
scratch retain their identity and are not loaded.
Collection destinations require the same number of ordered populations and
matching dimension/precision. Each population keeps its existing storage while
saved charge, mass, provenance and phase space replace its values. Every
population is validated before any is changed; scalar and collection layouts
are distinct and cannot be interchanged on restore.
Collisional destinations additionally require identical fixed timestep, model,
shape, collision normalization, RNG algorithm, per-population counts, charges
and masses. They restore the saved full RNG and public state in place without
resizing particles. Construct a fresh matching wrapper before loading; do not
call `init!` after loading, since that resets the physical clock.
`CollisionalEMPICStepper` additionally requires matching fixed electromagnetic
configuration and split policy. Its ordinary Vector destinations keep fixed
counts and storage while saved weights, provenance, phase space and nested
physical fields replace their values. The complete incoming RNG and inner state
are detached and validated before publication; the source clock and both fast/full
timestep histories restore with the RNG. No field projection or current refresh
is performed on restore.
`CollisionalElectrostaticPIC` similarly requires matching fixed populations,
configuration and full RNG algorithm/runtime. Its nested E/ne and velocity phase
restore verbatim, with the sole wrapper source clock/count; ES has no inner
time/step counters. Destination old physical, clock and inner-phase values may be
replaced after storage/configuration preflight. No Poisson solve or time
reconstruction from count times timestep is performed.
`YeeConductiveLayer` requires matching fixed grid/speed/profile/timestep policy.
Its nested fields/clock and exact budget restore verbatim; derived axis caches
prepare before any publication. Old destination fields, clock and ledger values
may be stale. No field projection or PML state substitution occurs; runtime
continuation acceptance remains deferred.
`YeePML` additionally validates/detaches all incoming directional memories and
their step-zero/PEC/zero-rate constraints. Old destination memory values are
replaceable; fixed profile/timestep policy and storage must match. No scalar
absorber or vacuum payload is silently upgraded to PML history.
Electrostatic destinations must have the same background density `n0`, which
is immutable; their saved fields, density, and leapfrog phase are restored directly.
Electromagnetic destinations require the same ion mobility and concrete shape
type; saved physics parameters replace the destination values. Supply independent
electron and ion destinations for mobile runs. Both species are validated before
either is changed. Ordered EMPIC1D collections instead require their distinct
collection tag and `:ex_ey_bz` polarization; all populations are validated using
the saved c and relativistic mode, with exactly planar stored velocities.
Their saved n_sub and species strides replace destination configuration together;
each stored species interval is restored independently of the current schedule.
Specialized schema 1 migrates only from its validated original envelope, with
all-fast strides and mobile history from its saved last_sub_dt. A schema-1 payload
carrying any schema-2 schedule fields is rejected. Destination configuration and
history values being replaced do not constrain an otherwise valid restore.
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

"Verify the run-container checksum, then use the canonical hybrid checkpoint restore."
function load_archive!(
    st::HybridStepper{D,T},
    species::_HybridPopulations{D,T},
    path::AbstractString,
) where {D,T}
    loaded = load_run(path)
    return _restore_checkpoint_state!(st, species, loaded.state)
end

function load_checkpoint!(st::HallMHDState, path::AbstractString)
    return _restore_checkpoint_state!(st, deserialize(path))
end

function load_checkpoint!(
    st::Union{_HybridCheckpointStepper,ElectrostaticPIC,EMPIC,EMPIC1D},
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
    # In-memory inputs may reference live destination fields or particles.
    # Detach all leaves together before any copy can overwrite a later source.
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, st)
    for (ps, record) in zip(species, saved.species)
        count = length(record.weight)
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            a isa Vector && length(a) != count && sizehint!(a, count)
        end
    end
    sizehint!(work.owners, length(species))
    for (ps, record) in zip(species, saved.species)
        _restore_checkpoint_particles!(ps, record)
    end
    copyto!(st.ne, saved.ne)
    foreach(copyto!, st.E, saved.E)
    empty!(work.owners)
    append!(work.owners, species)
    work.initialized[] = saved.initialized
    st.primed[] = saved.primed
    st.last_dt[] = saved.last_dt
    return st
end
