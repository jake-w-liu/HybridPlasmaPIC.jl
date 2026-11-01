"""
    YeeChargeWorkspace(g::YeeGrid{1})

Reusable CPU storage for [`yee_particle_charge!`](@ref) and
[`yee_longitudinal_impulse!`](@ref) on a closed 1D PEC grid. The payload is one
Float64 sum vector and one grid-precision candidate, each of length `n+1`.
A workspace belongs to one call at a time; its arrays must remain independent
of outputs and particle/endpoint storage.
"""
struct YeeChargeWorkspace{T<:Union{Float32,Float64}}
    n::Int
    L::T
    sum::Vector{Float64}
    candidate::Vector{T}
end

function YeeChargeWorkspace(g::YeeGrid{1,T}) where {T}
    g.periodic[1] && throw(ArgumentError("Yee charge coupling requires a closed PEC grid"))
    nodes = Base.checked_add(g.n[1], 1)
    return YeeChargeWorkspace{T}(g.n[1], g.L[1], zeros(Float64, nodes), zeros(T, nodes))
end

YeeChargeWorkspace(g::YeeGrid) =
    throw(ArgumentError("Yee charge coupling requires one spatial dimension"))

# Only totals at or below this limit enter the native error bound. Larger
# totals keep one replay sentinel, without narrowing array length to Int.
const _YEE_CHARGE_NATIVE_COUNT_LIMIT = Int128(1) << 51

@inline function _yee_charge_native_count(count::Int128, particles::Integer)
    cap = _YEE_CHARGE_NATIVE_COUNT_LIMIT + one(Int128)
    particles >= cap - count && return cap
    # Here 0 <= particles < cap-count, so conversion is exact and in range.
    return count + convert(Int128, particles)
end

_yee_charge_species(ps::ParticleSet) = (ps,)
_yee_charge_species(species::Tuple) = species
function _yee_charge_species(species::AbstractVector)
    Base.require_one_based_indexing(species)
    return species
end
_yee_charge_species(species) =
    throw(ArgumentError("Yee charge particles must be a ParticleSet, tuple or vector of species"))

function _yee_charge_endpoints(ps::ParticleSet, old, new)
    return (old,), (new,)
end
function _yee_charge_endpoints(species, old, new)
    old isa Union{Tuple,AbstractVector} && new isa Union{Tuple,AbstractVector} ||
        throw(ArgumentError("Yee endpoint collections must be tuples or one-based vectors"))
    old isa AbstractVector && Base.require_one_based_indexing(old)
    new isa AbstractVector && Base.require_one_based_indexing(new)
    length(old) == length(species) && length(new) == length(species) ||
        throw(DimensionMismatch("one old/new position vector is required per Yee charge species"))
    return old, new
end

function _yee_charge_positions(positions, count, ::Type{T}) where {T}
    positions isa AbstractVector{T} ||
        throw(ArgumentError("Yee charge positions must be vectors of grid precision $T"))
    length(positions) == count ||
        throw(DimensionMismatch("Yee endpoint count must equal its species particle count"))
    _check_particle_vector_axes(:yee_position, positions, count)
    return nothing
end

function _yee_charge_aliases(buffers, ps)
    for buffer in buffers
        _mightalias_particle_storage(buffer, ps) &&
            throw(ArgumentError("Yee charge output and workspace must not alias particle storage"))
    end
    return nothing
end

function _yee_endpoint_aliases(buffers, positions)
    for buffer in buffers
        _particle_arrays_mightalias(buffer, positions) &&
            throw(ArgumentError("Yee charge output and workspace must not alias endpoint storage"))
    end
    return nothing
end

function _yee_validate_charge_position(
    x::Union{Float32,Float64},
    length,
    species,
    particle,
    role,
    native_mode,
)
    inside = if !isfinite(x)
        false
    elseif native_mode
        zero(x) <= x <= length
    else
        0 <= _exact_rational(x) <= _exact_rational(length)
    end
    inside || throw(
        ArgumentError("Yee species $species $role position[$particle] must be finite in [0,L]"),
    )
    return nothing
end

@inline function _yee_charge_nonnegative(weight::Union{Float32,Float64}, native_mode)
    return isfinite(weight) && (native_mode ? weight >= zero(weight) : _exact_rational(weight) >= 0)
end

# Every value/shape/alias is checked across all populations before either scratch
# array is written. Unused mass/velocity/provenance values are not operands.
function _yee_validate_charge(
    out,
    species,
    old,
    new,
    g::YeeGrid{1,T},
    work,
    mode,
    native_mode,
) where {T}
    g.periodic[1] && throw(ArgumentError("Yee charge coupling requires a closed PEC grid"))
    work isa YeeChargeWorkspace{T} ||
        throw(ArgumentError("Yee charge workspace must match grid precision $T"))
    work.n == g.n[1] && work.L === g.L[1] ||
        throw(ArgumentError("Yee charge workspace must match the grid cell count and length"))
    nodes = Base.checked_add(g.n[1], 1)
    out isa Vector{T} || throw(ArgumentError("Yee charge output must be an ordinary Vector{$T}"))
    expected = mode isa Val{:charge} ? nodes : g.n[1]
    length(out) == expected || throw(DimensionMismatch("Yee charge output has the wrong length"))
    length(work.sum) == nodes && length(work.candidate) == nodes ||
        throw(DimensionMismatch("Yee charge workspace arrays must have n+1 entries"))
    buffers = (out, work.sum, work.candidate)
    (
        _particle_arrays_mightalias(work.sum, out) ||
        _particle_arrays_mightalias(work.candidate, out) ||
        _particle_arrays_mightalias(work.candidate, work.sum)
    ) && throw(ArgumentError("Yee charge output and workspace arrays must be independent"))
    count = Int128(0)
    magnitude = 0.0
    for k in eachindex(species)
        ps = species[k]
        ps isa ParticleSet{1,T} ||
            throw(ArgumentError("Yee charge species must be one-dimensional at grid precision $T"))
        count, magnitude = _yee_validate_charge_population(
            buffers,
            ps,
            old,
            new,
            g,
            k,
            mode,
            count,
            magnitude,
            native_mode,
        )
    end
    return count, magnitude
end

@noinline function _yee_validate_charge_population(
    buffers,
    ps::ParticleSet{1,T},
    old,
    new,
    g,
    k,
    mode,
    count::Int128,
    magnitude::Float64,
    native_mode,
) where {T}
    particles = _check_particle_storage_axes(ps)
    count = _yee_charge_native_count(count, particles)
    _yee_charge_aliases(buffers, ps)
    if mode isa Val{:impulse}
        before = old[k]
        _yee_charge_positions(before, particles, T)
        after = new[k]
        _yee_charge_positions(after, particles, T)
        _yee_endpoint_aliases(buffers, before)
        _yee_endpoint_aliases(buffers, after)
        before isa StridedVector{T} && after isa StridedVector{T} ||
            throw(ArgumentError("Yee charge endpoints must use CPU strided vectors"))
    else
        before, after = nothing, nothing
    end
    all(a -> a isa StridedVector, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)) ||
        throw(ArgumentError("Yee charge particle components must use CPU strided vectors"))
    isfinite(ps.q) || throw(ArgumentError("Yee species $k charge q must be finite"))
    magnitude = _yee_validate_charge_values(ps, before, after, g, k, mode, magnitude, native_mode)
    return count, magnitude
end

@noinline function _yee_validate_charge_values(
    ps::ParticleSet{1,T},
    old,
    new,
    g,
    k,
    mode,
    magnitude::Float64,
    native_mode,
) where {T}
    for p in eachindex(ps.weight)
        weight = ps.weight[p]
        _yee_charge_nonnegative(weight, native_mode) ||
            throw(ArgumentError("Yee species $k weight[$p] must be finite and nonnegative"))
        if mode isa Val{:charge}
            _yee_validate_charge_position(ps.x[1][p], g.L[1], k, p, :current, native_mode)
        else
            _yee_validate_charge_position(old[p], g.L[1], k, p, :old, native_mode)
            _yee_validate_charge_position(new[p], g.L[1], k, p, :new, native_mode)
        end
        if native_mode
            magnitude =
                _eft_add_upper(magnitude, _eft_mul_upper(abs(Float64(ps.q)), Float64(weight)))
        end
    end
    return magnitude
end

@inline _yee_macrocharge(ps, p, ::Val{false}) = Float64(ps.q) * Float64(ps.weight[p])
@inline _yee_macrocharge(ps, p, ::Val{true}) = _exact_rational(ps.q) * _exact_rational(ps.weight[p])

@inline _yee_macrocharge_admitted(value, ::Val{false}) =
    isfinite(value) && abs(value) >= floatmin(Float64)
@inline _yee_macrocharge_admitted(value, ::Val{true}) = true

@inline function _yee_add_charge!(sums, i, charge, coefficient, ::Val{false})
    iszero(coefficient) && return true
    value = charge * coefficient
    isfinite(value) && abs(value) >= floatmin(Float64) || return false
    updated = sums[i] + value
    isfinite(updated) || return false
    sums[i] = updated
    return true
end
@inline function _yee_add_charge!(sums, i, charge, coefficient, ::Val{true})
    iszero(coefficient) || (sums[i] += charge * coefficient)
    return true
end

function _yee_charge_kernel!(sums, species, old, new, g, exact, mode::Val{:charge})
    for ps in species
        _yee_charge_population!(sums, ps, nothing, nothing, g, exact, mode) || return false
    end
    return true
end

function _yee_charge_kernel!(sums, species, old, new, g, exact, mode::Val{:impulse})
    for k in eachindex(species)
        _yee_charge_population!(sums, species[k], old[k], new[k], g, exact, mode) || return false
    end
    return true
end

@noinline function _yee_charge_population!(sums, ps, old, new, g, exact, ::Val{:charge})
    for p in eachindex(ps.weight)
        if exact isa Val{false}
            (iszero(ps.q) || iszero(ps.weight[p])) && continue
        end
        charge = _yee_macrocharge(ps, p, exact)
        exact isa Val{true} && iszero(charge) && continue
        _yee_macrocharge_admitted(charge, exact) || return false
        cell, weights, admitted = _yee_cic_stencil(ps.x[1][p], g, exact)
        admitted || return false
        _yee_add_charge!(sums, cell + 1, charge, weights[1], exact) || return false
        _yee_add_charge!(sums, cell + 2, charge, weights[2], exact) || return false
    end
    return true
end

@noinline function _yee_charge_population!(sums, ps, old, new, g, exact, ::Val{:impulse})
    for p in eachindex(ps.weight)
        if exact isa Val{false}
            (iszero(ps.q) || iszero(ps.weight[p])) && continue
            first, last = old[p], new[p]
            first == last && continue
            a, b = minmax(first, last)
            charge = _yee_macrocharge(ps, p, exact)
            _yee_macrocharge_admitted(charge, exact) || return false
            first > last && (charge = -charge)
        else
            charge = _yee_macrocharge(ps, p, exact)
            iszero(charge) && continue
            first, last = old[p], new[p]
            direction = cmp(_exact_rational(first), _exact_rational(last))
            iszero(direction) && continue
            a, b = direction < 0 ? (first, last) : (last, first)
            direction > 0 && (charge = -charge)
        end
        lower, upper, weights, admitted = _yee_cic_path_stencil(a, b, g, exact)
        admitted || return false
        if lower == upper
            _yee_add_charge!(sums, lower + 1, charge, weights[1], exact) || return false
        else
            _yee_add_charge!(sums, lower + 1, charge, weights[1], exact) || return false
            for cell = (lower+1):(upper-1)
                _yee_add_charge!(sums, cell + 1, charge, one(charge), exact) || return false
            end
            _yee_add_charge!(sums, upper + 1, charge, weights[2], exact) || return false
        end
    end
    return true
end

function _yee_charge_native_candidate!(work, count, magnitude, entries)
    # Coefficient error <=8u, two products and P additions: gamma_(P+10).
    eta = (Float64(count) + 10.0) * 0x1p-53
    gamma = _eft_div_upper(eta, prevfloat(1.0 - eta))
    error = _eft_mul_upper(gamma, magnitude)
    T = eltype(work.candidate)
    maximum = Float64(floatmax(T))
    for i = 1:entries
        value = work.sum[i]
        isfinite(value) && !issubnormal(value) || return false
        _eft_add_upper(abs(value), error) <= maximum || return false
        candidate = T(value)
        isfinite(candidate) && !issubnormal(candidate) || return false
        work.candidate[i] = candidate
    end
    return true
end

_yee_validate_charge_shape(::CIC) = nothing
_yee_validate_charge_shape(::YeeTriangularCloud) = nothing
_yee_validate_charge_shape(shape) =
    throw(ArgumentError("Yee charge shape must be CIC() or YeeTriangularCloud"))

_yee_charge_kernel!(sums, species, old, new, g, exact, mode, ::CIC) =
    _yee_charge_kernel!(sums, species, old, new, g, exact, mode)

function _yee_charge_kernel!(sums, species, old, new, g, ::Val{true}, mode, shape::YeeTriangularCloud)
    for k in eachindex(species)
        before = mode isa Val{:charge} ? nothing : old[k]
        after = mode isa Val{:charge} ? nothing : new[k]
        _yee_charge_population!(sums, species[k], before, after, g, Val(true), mode, shape)
    end
    return true
end

@noinline function _yee_charge_population!(sums, ps::ParticleSet{1,T}, old, new, g, ::Val{true}, mode, shape::YeeTriangularCloud) where {T}
    h = _exact_rational(g.L[1]) / g.n[1]
    a = _exact_rational(shape.half_width)
    for p in eachindex(ps.weight)
        charge = _yee_macrocharge(ps, p, Val(true))
        iszero(charge) && continue
        if mode isa Val{:charge}
            x = _exact_rational(ps.x[1][p])
            previous = one(h)
            edges = _yee_cloud_edge_span(x, x, h, a, g.n[1])
            for e in edges
                cumulative, _ = _yee_cloud_edge(x, e * h, h, a)
                sums[e+1] += charge * (previous - cumulative)
                previous = cumulative
            end
            sums[last(edges)+2] += charge * previous
        else
            before = _exact_rational(old[p])
            after = _exact_rational(new[p])
            for e in _yee_cloud_edge_span(before, after, h, a, g.n[1])
                left, _ = _yee_cloud_edge(before, e * h, h, a)
                right, _ = _yee_cloud_edge(after, e * h, h, a)
                sums[e+1] += charge * (right - left)
            end
        end
    end
    return true
end

@noinline function _yee_charge_exact_candidate!(work, species, old, new, g, mode, entries, shape)
    sums = fill(zero(Rational{BigInt}), Base.checked_add(g.n[1], 1))
    _yee_charge_kernel!(sums, species, old, new, g, Val(true), mode, shape)
    T = eltype(work.candidate)
    for i = 1:entries
        candidate = _from_exact(T, sums[i])
        isfinite(candidate) ||
            throw(ArgumentError("Yee charge/impulse entry $i does not round to a finite $T value"))
        work.candidate[i] = candidate
    end
    return nothing
end

function _yee_charge_apply!(out, species, old, new, g::YeeGrid{1}, work, mode, shape)
    _yee_validate_charge_shape(shape)
    # The native EFT bounds require nearest rounding and gradual underflow.
    # Other contexts keep all zero/order/geometry decisions on stored rationals.
    native_mode = shape isa CIC && _boris_accuracy_mode()
    count, magnitude = _yee_validate_charge(out, species, old, new, g, work, mode, native_mode)
    fill!(work.sum, 0.0)
    fill!(work.candidate, zero(eltype(work.candidate)))
    entries = length(out)
    # The count bound makes eta<1/2. It selects replay, never rejects a population.
    native = native_mode && count <= _YEE_CHARGE_NATIVE_COUNT_LIMIT && isfinite(magnitude)
    native = native && _yee_charge_kernel!(work.sum, species, old, new, g, Val(false), mode)
    native = native && _yee_charge_native_candidate!(work, count, magnitude, entries)
    native || _yee_charge_exact_candidate!(work, species, old, new, g, mode, entries, shape)
    for i in eachindex(out)
        out[i] = work.candidate[i]
    end
    return out
end

"""
    yee_particle_charge!(out, particles, g::YeeGrid{1}, workspace; shape=CIC())

Overwrite `out` with extensive CIC particle charge on every node of a closed PEC
interval `[0,L]`, including both wall nodes. Accept one ParticleSet or a tuple/
one-based vector of species, with finite q and nonnegative finite weights. This
is charge, without division by dx. Geometry is `s=n*x/L` from stored inputs.
Only CPU strided particle components and ordinary grid-precision output vectors
are supported. Inputs remain unchanged; output/workspace cannot alias them.
Mass/velocity/provenance values are not operands. Validation precedes scratch
writes; arithmetic failure preserves output. Complete nonfinite rounded output
rejects before publication; intermediate range loss replays exact original inputs.
Nondefault hardware rounding or subnormal flushing uses exact replay, including
stored-input admission and geometry; final exact rounding remains nearest-even.
See `docs/BOUNDED_PARTICLE_COUPLING_PLAN.md` for the normwise rounding boundary.
Select `shape=YeeTriangularCloud(a)` for the optional finite physical-width
cloud, retaining its original width precision. Its tails are particle-bound
endpoint support. That shape uses exact whole-entry evaluation and one final
rounding, with the same preflight and publication contracts; the point/CIC native
bound does not apply to cloud polynomials.
"""
function yee_particle_charge!(out, particles, g::YeeGrid{1}, workspace; shape = CIC())
    return _yee_charge_apply!(
        out,
        _yee_charge_species(particles),
        nothing,
        nothing,
        g,
        workspace,
        Val(:charge),
        shape,
    )
end

"""
    yee_longitudinal_impulse!(out, particles, xold, xnew, g::YeeGrid{1}, workspace; shape=CIC())

Overwrite `out` with the longitudinal integrated current compatible with CIC
endpoint charge in a closed PEC interval, epsilon0=1 and unit transverse area.
A single species takes two position vectors; a collection takes matching tuples/
one-based vectors of position vectors. q/weight/identity belong to the same
unchanged population. All consumed positions lie in `[0,L]`, with exact wall
endpoints allowed. The particle set's current position values are not endpoints.

No dt argument or additional dx normalization is needed. In exact arithmetic,
`Qnew-Qold+[out[1],diff(out)...,-out[end]]=0`; independently rounded arrays have
a normwise residual, not bitwise continuity. Reversed/retraced normal paths cancel
and multiple-cell flux is included. No transverse current, force, reflection
mover, absorbed charge, voltage constraint or coupled solver is supplied.
Storage, range and atomic-publication rules follow [`yee_particle_charge!`](@ref).
Use the same optional `YeeTriangularCloud` shape as the endpoint charges and
path voltage. Cloud impulse is the exact cumulative difference, with all
populations accumulated before one final entry rounding.
"""
function yee_longitudinal_impulse!(out, particles, xold, xnew, g::YeeGrid{1}, workspace; shape = CIC())
    species = _yee_charge_species(particles)
    old, new = _yee_charge_endpoints(particles, xold, xnew)
    return _yee_charge_apply!(out, species, old, new, g, workspace, Val(:impulse), shape)
end

yee_particle_charge!(out, particles, g::YeeGrid, workspace; shape = CIC()) =
    throw(ArgumentError("Yee charge coupling requires one spatial dimension"))
yee_longitudinal_impulse!(out, particles, old, new, g::YeeGrid, workspace; shape = CIC()) =
    throw(ArgumentError("Yee charge coupling requires one spatial dimension"))

function _yee_validate_longitudinal_storage(Ex, g::YeeGrid{1,T}) where {T}
    g.periodic[1] && throw(ArgumentError("Yee longitudinal coupling requires a closed PEC grid"))
    Ex isa StridedVector{T} || throw(
        ArgumentError("Yee longitudinal field must be a CPU strided vector of grid precision $T"),
    )
    length(Ex) == g.n[1] || throw(DimensionMismatch("Yee longitudinal field must have n entries"))
    Base.require_one_based_indexing(Ex)
    return nothing
end

"""
    yee_longitudinal_voltage(Ex, xold, xnew, g::YeeGrid{1}; shape=CIC())

Return the oriented longitudinal path voltage paired with
[`yee_longitudinal_impulse!`](@ref) on a closed 1D PEC interval. For the exact
stored-coordinate cell overlaps `C[e]` in `s=n*x/L`, this is
`g.dx[1] * sum(Ex[e]*C[e])`. It includes the longitudinal DC field without
imposing a plate-voltage constraint. Reversing endpoints negates the exact
voltage; a zero path returns zero after validation.

`Ex` is a read-only CPU strided vector of grid precision with n finite entries;
both endpoints have that precision and lie in `[0,L]`. The complete expression
uses exact stored operands and rounds once to the grid precision, independently
of hardware rounding or subnormal flushing. A nonfinite final rounding throws
`ArgumentError`; inputs are unchanged on success or failure. Validation scans
all n field entries; contraction costs O(crossed cells) rational operations.
No allocation or latency bound is supplied.

Before independent rounding, particle work is `q*weight*voltage`, the adjoint
of the existing impulse under the stored-dx field-work metric. For a nonzero
physical displacement the associated path average is `voltage/(xnew-xold)`;
it includes `n*g.dx[1]/g.L[1]`, which need not equal one exactly. Supply a temporal
mean field to pair with source work, allowing for that mean's own rounding.
This operator selects no stationary-node force, particle mover or wall event,
and does not establish coupled energy conservation or force accuracy.
The optional `shape=YeeTriangularCloud(a)` instead contracts that cloud's exact
cumulative differences over all n cells. Its stationary derivative is available
separately through [`yee_longitudinal_field`](@ref), requiring explicit cloud
selection. Point/CIC geometry and default behavior remain unchanged.
"""
function yee_longitudinal_voltage(Ex, xold, xnew, g::YeeGrid{1,T}; shape = CIC()) where {T}
    _yee_validate_charge_shape(shape)
    _yee_validate_longitudinal_storage(Ex, g)
    xold isa T && xnew isa T ||
        throw(ArgumentError("Yee path voltage endpoints must have grid precision $T"))
    _yee_validate_charge_position(xold, g.L[1], 1, 1, :old, false)
    _yee_validate_charge_position(xnew, g.L[1], 1, 1, :new, false)
    all(isfinite, Ex) || throw(ArgumentError("Yee path voltage field must be finite"))
    direction = cmp(_exact_rational(xold), _exact_rational(xnew))
    iszero(direction) && return zero(T)
    shape isa YeeTriangularCloud && return _yee_cloud_field_contraction(Ex, xold, xnew, g, shape)
    a, b = direction < 0 ? (xold, xnew) : (xnew, xold)
    lower, upper, weights, _ = _yee_cic_path_stencil(a, b, g, Val(true))
    total = _exact_rational(Ex[lower+1]) * weights[1]
    if lower != upper
        for cell = (lower+1):(upper-1)
            total += _exact_rational(Ex[cell+1])
        end
        total += _exact_rational(Ex[upper+1]) * weights[2]
    end
    direction > 0 && (total = -total)
    value = _from_exact(T, _exact_rational(g.dx[1]) * total)
    isfinite(value) || throw(ArgumentError("Yee path voltage does not round to a finite $T value"))
    return value
end

yee_longitudinal_voltage(Ex, xold, xnew, g::YeeGrid; shape = CIC()) =
    throw(ArgumentError("Yee path voltage requires one spatial dimension"))

# A nonzero path contracts A(new)-A(old); a stationary field contracts g(x).
# Both share original geometry/width, the stored-dx metric and one final round.
function _yee_cloud_field_contraction(Ex, old, new, g::YeeGrid{1,T}, shape::YeeTriangularCloud) where {T}
    h = _exact_rational(g.L[1]) / g.n[1]
    a = _exact_rational(shape.half_width)
    x = _exact_rational(old)
    y = new === nothing ? nothing : _exact_rational(new)
    total = zero(h)
    for e in _yee_cloud_edge_span(x, y === nothing ? x : y, h, a, g.n[1])
        cumulative, derivative = _yee_cloud_edge(x, e * h, h, a)
        coefficient = if y === nothing
            derivative
        else
            following, _ = _yee_cloud_edge(y, e * h, h, a)
            following - cumulative
        end
        total += _exact_rational(Ex[e+1]) * coefficient
    end
    value = _from_exact(T, _exact_rational(g.dx[1]) * total)
    isfinite(value) || throw(ArgumentError("Yee cloud field contraction does not round to a finite $T value"))
    return value
end

"""
    yee_longitudinal_field(Ex, x, g::YeeGrid{1}; shape::YeeTriangularCloud)

Stationary work-adjoint longitudinal field for the explicitly selected finite
physical-width cloud. Returns `dx*sum(Ex[e]*g_e(x))`, where `g_e` is the same
cloud cumulative's derivative used by charge/current/path voltage. This is a
field, without charge/mass, acceleration or a particle mover. Mesh-node and
overlapping-cloud states have a continuous force basis. No point/CIC stationary
force is selected.

The width retains its original IEEE precision; field and coordinate have grid
precision. Inputs are read-only and validated before calculation. Exact stored
operands round once to grid precision, nearest-even independently of arithmetic
context; nonfinite final rounding rejects. Interior DC includes `n*dx/L`.
Particle-bound endpoint tails give half that DC response at either wall when
the half-width is at most L. This is the declared cloud/surface shape, not an
absorbed/free plate charge or a point-wall force. Coupled PIC remains incomplete.
"""
function yee_longitudinal_field(Ex, x, g::YeeGrid{1,T}; shape::YeeTriangularCloud) where {T}
    _yee_validate_longitudinal_storage(Ex, g)
    x isa T || throw(ArgumentError("Yee cloud coordinate must have grid precision $T"))
    _yee_validate_charge_position(x, g.L[1], 1, 1, :current, false)
    all(isfinite, Ex) || throw(ArgumentError("Yee cloud field must be finite"))
    return _yee_cloud_field_contraction(Ex, x, nothing, g, shape)
end

yee_longitudinal_field(Ex, x, g::YeeGrid; shape::YeeTriangularCloud) =
    throw(ArgumentError("Yee cloud field requires one spatial dimension"))
