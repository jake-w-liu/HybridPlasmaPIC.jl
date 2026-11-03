# Full-f electrostatic gyrocentres in a constant-B periodic Cartesian slab.
# Physical mu=m*vperp^2/(2B), coherent SI-form units, Te in energy units.

"""
    GyrocentreSet{T}(N; q=1, m=1)

CPU full-f gyrocentre population with three position vectors `X`, parallel
velocity `vpar`, physical magnetic moment `mu`, number `weight`, IDs and tags.
Populate initial data before constructing `SlabGyrokineticPIC`. No artificial
perpendicular particle velocities or delta-f weights are stored.
"""
struct GyrocentreSet{T<:Union{Float32,Float64}}
    X::NTuple{3,Vector{T}}
    vpar::Vector{T}
    mu::Vector{T}
    weight::Vector{T}
    id::Vector{UInt64}
    tag::Vector{UInt32}
    q::T
    m::T
end
function GyrocentreSet{T}(N::Integer; q = one(T), m = one(T)) where {T<:Union{Float32,Float64}}
    _boris_accuracy_mode() ||
        throw(ArgumentError("gyrocentres require nearest rounding and gradual underflow"))
    count = _particle_length(N)
    charge = _require_finite_positive_real("gyrocentre q", q, T)
    mass = _require_finite_positive_real("gyrocentre m", m, T)
    return GyrocentreSet{T}(
        ntuple(_ -> zeros(T, count), 3),
        zeros(T, count),
        zeros(T, count),
        ones(T, count),
        _particle_ids(count),
        zeros(UInt32, count),
        charge,
        mass,
    )
end
nparticles(p::GyrocentreSet) = length(p.weight)
Base.eltype(::GyrocentreSet{T}) where {T} = T
_sg_arrays(p::GyrocentreSet) = (p.X..., p.vpar, p.mu, p.weight, p.id, p.tag)

struct _SlabGyroStage{T}
    X::NTuple{3,Vector{T}}
    vpar::Vector{T}
end
_SlabGyroStage(::Type{T}, N::Int) where {T} =
    _SlabGyroStage{T}(ntuple(_ -> zeros(T, N), 3), zeros(T, N))
_sg_arrays(p::_SlabGyroStage) = (p.X..., p.vpar)

"Scalar field receipt; saturated displayed values do not replace exact acceptance comparisons."
struct SlabGyrokineticFieldReport{T}
    residual::T
    tolerance::T
    source_norm::T
    gauge::T
    gauge_tolerance::T
    ion_dc_mean::T
    ion_dc_imaginary::T
    ion_dc_inventory_defect::T
    charge_defect::T
    electron_response_amplitude::T
    range_resolved::Bool
end

struct _SlabGyroWork{T}
    old::Vector{_SlabGyroStage{T}}
    trial::Vector{_SlabGyroStage{T}}
    endpoint::Vector{_SlabGyroStage{T}}
    rates::NTuple{4,Vector{_SlabGyroStage{T}}}
    rho::Array{T,3}
    phi::Array{T,3}
    source::Array{Complex{T},3}
    phihat::Array{Complex{T},3}
end

"""
    SlabGyrokineticPIC(g, populations; B0, Te, e=1, nring=16,
                      field_rtol=sqrt(eps(T)), field_atol=0)

Collective 3D periodic CPU full-f positive-ion gyrokinetics with fixed physical
mu/number weights, equilibrium long-wavelength polarization and modified
adiabatic electrons. Call `init!`, then `step!(state, populations, dt)`.
The four RK4 stages and completed endpoint each solve the common potential.
See `docs/SLAB_GYROKINETIC.md` for units, nullspace and acceptance limitations.
"""
mutable struct SlabGyrokineticPIC{T,G}
    const g::G
    const owners::Vector{GyrocentreSet{T}}
    const B0::T
    const Te::T
    const e::T
    const nring::Int
    const field_rtol::T
    const field_atol::T
    const chi::T
    const alpha::T
    const ne0_display::T
    const rho_gc::Array{T,3}
    const phi::Array{T,3}
    const _counts::Vector{Int}
    const _qm::Vector{NTuple{2,T}}
    const _mu::Vector{Vector{T}}
    const _weight::Vector{Vector{T}}
    const _radius::Vector{Vector{T}}
    const _components::Vector{Tuple}
    const _phases::Vector{NTuple{2,T}}
    const _lambda::Array{T,3}
    const _gain::T
    const work::_SlabGyroWork{T}
    initialized::Bool
    time::T
    step::Int
    last_dt::T
    last_report::Union{Nothing,SlabGyrokineticFieldReport{T}}
    busy::Bool
end

_sg_volume(g) = prod(_exact_rational, g.L)
_sg_cell_volume(g) = _sg_volume(g) / length(g.cbuf)
function _sg_inventory(populations, g, e)
    Q = _exact_rational
    volume = _sg_volume(g)
    charge = sum(
        (Q(p.q) * _exact_float_sum(p.weight) for p in populations);
        init = zero(Rational{BigInt}),
    )
    mass = sum(
        (Q(p.m) * _exact_float_sum(p.weight) for p in populations);
        init = zero(Rational{BigInt}),
    )
    return (charge = charge, mass = mass, ne0 = charge / (Q(e) * volume))
end
function _sg_population(p::GyrocentreSet{T}; values = true) where {T}
    N = nparticles(p)
    all(a -> length(a) == N, _sg_arrays(p)) ||
        throw(DimensionMismatch("gyrocentre component counts differ"))
    isfinite(p.q) && p.q > 0 && isfinite(p.m) && p.m > 0 ||
        throw(ArgumentError("gyrocentre charge and mass must be finite positive stored values"))
    if values
        all(a -> all(isfinite, a), (p.X..., p.vpar, p.mu, p.weight)) ||
            throw(ArgumentError("gyrocentre physical values must be finite"))
        all(>=(zero(T)), p.mu) && all(>=(zero(T)), p.weight) ||
            throw(ArgumentError("gyrocentre mu and full-f weights must be nonnegative"))
    end
    return nothing
end
function _sg_grid(g::FourierGrid{3,T}) where {T}
    all(>=(3), g.n) ||
        throw(ArgumentError("slab gyrokinetics requires at least three cells per axis"))
    metadata = SpectralOperators._fourier_metadata(g.n, g.L)
    isequal(g.dx, metadata.dx) &&
        isequal(g.midx, metadata.midx) &&
        isequal(g.kfull, metadata.kfull) &&
        isequal(g.kvec, metadata.kvec) &&
        isequal(g.ik, metadata.ik) || throw(ArgumentError("gyrokinetic Fourier metadata changed"))
    all(a -> a isa Array{Complex{T},3} && size(a) == g.n, (g.cbuf, g.tbuf, g.abuf, g.sbuf)) ||
        throw(ArgumentError("gyrokinetic Fourier buffers must be matching ordinary Arrays"))
    return nothing
end
function _sg_work_arrays(w::_SlabGyroWork)
    arrays = Any[w.rho, w.phi, w.source, w.phihat]
    for collection in (w.old, w.trial, w.endpoint, w.rates...), p in collection
        append!(arrays, _sg_arrays(p))
    end
    return arrays
end
function _sg_storage_arrays(st)
    arrays = Any[st.rho_gc, st.phi, st._counts, st._qm, st._phases, st._lambda]
    for p in st.owners
        append!(arrays, _sg_arrays(p))
    end
    append!(arrays, st._mu)
    append!(arrays, st._weight)
    append!(arrays, st._radius)
    append!(arrays, _sg_work_arrays(st.work))
    append!(
        arrays,
        (
            st.g.cbuf,
            st.g.tbuf,
            st.g.abuf,
            st.g.sbuf,
            st.g.midx...,
            st.g.kfull...,
            st.g.kvec...,
            st.g.ik...,
        ),
    )
    return arrays
end
function _sg_storage(st)
    arrays = _sg_storage_arrays(st)
    for i in eachindex(arrays), j = 1:(i-1)
        _particle_arrays_mightalias(arrays[i], arrays[j]) && throw(
            ArgumentError("gyrokinetic public/private/metadata storage roles must be independent"),
        )
    end
    # FFT plans may retain arrays in addition to the four explicit grid buffers.
    _visit_checkpoint_objects((st.g.plan, st.g.iplan)) do value
        if value isa AbstractArray
            for a in arrays
                any(x -> value === x, (st.g.cbuf, st.g.tbuf, st.g.abuf, st.g.sbuf)) && continue
                _particle_arrays_mightalias(value, a) && throw(
                    ArgumentError(
                        "gyrokinetic physical/configuration storage overlaps FFT plan storage",
                    ),
                )
            end
        end
    end
    return nothing
end
function _sg_geometry(st)
    N = length(st.owners)
    all(
        a -> length(a) == N,
        (st._counts, st._qm, st._mu, st._weight, st._radius, st._components),
    ) || throw(DimensionMismatch("gyrokinetic fixed inventory geometry changed"))
    length(st._phases) == st.nring && all(phase -> all(isfinite, phase), st._phases) ||
        throw(ArgumentError("gyrokinetic fixed quadrature geometry changed"))
    all(
        a -> size(a) == st.g.n,
        (st.rho_gc, st.phi, st._lambda, st.work.rho, st.work.phi, st.work.source, st.work.phihat),
    ) || throw(DimensionMismatch("gyrokinetic field/work geometry changed"))
    for k in eachindex(st.owners)
        all(a -> length(a) == st._counts[k], (st._mu[k], st._weight[k], st._radius[k])) ||
            throw(DimensionMismatch("gyrokinetic invariant/radius geometry changed"))
    end
    for collection in (st.work.old, st.work.trial, st.work.endpoint, st.work.rates...)
        length(collection) == N ||
            throw(DimensionMismatch("gyrokinetic stage population count changed"))
        for k in eachindex(collection)
            all(a -> length(a) == st._counts[k], _sg_arrays(collection[k])) ||
                throw(DimensionMismatch("gyrokinetic stage component geometry changed"))
        end
    end
    return nothing
end
function _sg_bind(st::SlabGyrokineticPIC{T}, populations; values = true) where {T}
    _boris_accuracy_mode() ||
        throw(ArgumentError("slab gyrokinetics requires nearest rounding and gradual underflow"))
    populations isa AbstractVector ||
        throw(ArgumentError("gyrocentre populations must be an ordered vector"))
    Base.require_one_based_indexing(populations)
    length(populations) == length(st.owners) ||
        throw(DimensionMismatch("gyrocentre population count changed"))
    _sg_grid(st.g)
    _sg_geometry(st)
    for k in eachindex(st.owners)
        p = populations[k]
        p === st.owners[k] || throw(ArgumentError("gyrocentre population identity/order changed"))
        _sg_population(p; values)
        nparticles(p) == st._counts[k] ||
            throw(DimensionMismatch("gyrocentre particle count changed"))
        all(pair -> pair[1] === pair[2], zip(_sg_arrays(p), st._components[k])) ||
            throw(ArgumentError("gyrocentre component identities changed"))
        isequal((p.q, p.m), st._qm[k]) &&
            isequal(p.mu, st._mu[k]) &&
            isequal(p.weight, st._weight[k]) ||
            throw(ArgumentError("gyrokinetic fixed charges/masses/mu/weights changed"))
    end
    _sg_storage(st)
    return nothing
end
function _sg_acquire(st)
    st.busy && throw(ArgumentError("gyrokinetic state/grid workspace is already in use"))
    st.busy = true
    return nothing
end
function _sg_phase(initialized, time::T, step::Int, h::T, report) where {T}
    initialized && isfinite(time) && time >= 0 && step >= 0 && isfinite(h) && h >= 0 ||
        throw(ArgumentError("invalid synchronous gyrokinetic phase"))
    (step == 0 ? iszero(time) && iszero(h) : time > 0 && 0 < h <= time) ||
        throw(ArgumentError("gyrokinetic clock/completed interval disagree"))
    report isa SlabGyrokineticFieldReport{T} ||
        throw(ArgumentError("gyrokinetic boundary requires its field receipt"))
    return nothing
end
function _sg_boundary(st)
    _sg_phase(st.initialized, st.time, st.step, st.last_dt, st.last_report)
    all(isfinite, st.rho_gc) && all(isfinite, st.phi) ||
        throw(ArgumentError("gyrokinetic boundary fields must be finite"))
    for p in st.owners, d = 1:3, x in p.X[d]
        0 <= x < st.g.L[d] ||
            throw(ArgumentError("gyrocentre boundary positions must be periodic representatives"))
    end
    return nothing
end

function SlabGyrokineticPIC(
    g::FourierGrid{3,T},
    populations::AbstractVector;
    B0,
    Te,
    e = one(T),
    nring = 16,
    field_rtol = sqrt(eps(T)),
    field_atol = zero(T),
) where {T<:Union{Float32,Float64}}
    _boris_accuracy_mode() ||
        throw(ArgumentError("slab gyrokinetics requires nearest rounding and gradual underflow"))
    _sg_grid(g)
    Base.require_one_based_indexing(populations)
    isempty(populations) && throw(ArgumentError("at least one gyrocentre population is required"))
    all(p -> p isa GyrocentreSet{T}, populations) ||
        throw(ArgumentError("gyrocentre populations must match grid precision"))
    owners = GyrocentreSet{T}[p for p in populations]
    foreach(_sg_population, owners)
    B, temp, charge = (
        _require_finite_positive_real(name, value, T) for
        (name, value) in (("B0", B0), ("Te", Te), ("e", e))
    )
    count = _require_positive_intlike("nring", nring)
    count >= 3 && Float64(count) == count || throw(
        ArgumentError("nring must be at least three and exactly representable in phase precision"),
    )
    rtol = _require_finite_nonnegative_real("field_rtol", field_rtol, T)
    atol = _require_finite_nonnegative_real("field_atol", field_atol, T)
    rtol < 1 && (!iszero(rtol) || !iszero(atol)) ||
        throw(ArgumentError("field tolerances require rtol<1 and a positive tolerance"))
    Q = _exact_rational
    inventory = _sg_inventory(owners, g, charge)
    inventory.charge > 0 ||
        throw(ArgumentError("gyrocentre collection must carry positive full-f number"))
    chi = _from_exact(T, inventory.mass / (_sg_volume(g) * Q(B)^2))
    alpha = _from_exact(T, Q(charge)^2 * inventory.ne0 / Q(temp))
    display = _from_exact(T, inventory.ne0)
    all(x -> isfinite(x) && x >= floatmin(T), (chi, alpha, display)) || throw(
        ArgumentError(
            "gyrokinetic susceptibilities and displayed equilibrium density must be positive normal values",
        ),
    )
    counts = map(nparticles, owners)
    radii = map(owners) do p
        map(p.mu) do mu
            radius = _sqrt_from_exact(T, 2Q(p.m) * Q(mu) / (Q(p.q)^2 * Q(B)))
            isfinite(radius) && (iszero(mu) || !iszero(radius)) || throw(
                ArgumentError("gyrocentre radius must be finite and positive when mu is positive"),
            )
            radius
        end
    end
    phases = NTuple{2,T}[]
    for k = 0:(count-1)
        sine, cosine = sincospi(2 * (Float64(k) / Float64(count)))
        push!(phases, (T(cosine), T(sine)))
    end
    lambda = zeros(T, g.n)
    gain = one(T)
    for I in CartesianIndices(lambda)
        all(==(1), Tuple(I)) && continue
        symbol = Q(chi) * (Q(g.kfull[1][I[1]])^2 + Q(g.kfull[2][I[2]])^2)
        (g.midx[2][I[2]] != 0 || g.midx[3][I[3]] != 0) && (symbol += Q(alpha))
        value = _from_exact(T, symbol)
        isfinite(value) && value >= floatmin(T) || throw(
            ArgumentError("gyrokinetic non-DC field symbols must be positive finite normal values"),
        )
        lambda[I] = value
        gain = max(gain, inv(value))
    end
    isfinite(gain) || throw(ArgumentError("gyrokinetic inverse field gain must be finite"))
    stages() = [_SlabGyroStage(T, N) for N in counts]
    work = _SlabGyroWork{T}(
        stages(),
        stages(),
        stages(),
        ntuple(_ -> stages(), 4),
        zeros(T, g.n),
        zeros(T, g.n),
        zeros(Complex{T}, g.n),
        zeros(Complex{T}, g.n),
    )
    st = SlabGyrokineticPIC{T,typeof(g)}(
        g,
        owners,
        B,
        temp,
        charge,
        count,
        rtol,
        atol,
        chi,
        alpha,
        display,
        zeros(T, g.n),
        zeros(T, g.n),
        counts,
        [(p.q, p.m) for p in owners],
        [copy(p.mu) for p in owners],
        [copy(p.weight) for p in owners],
        radii,
        Tuple[_sg_arrays(p) for p in owners],
        phases,
        lambda,
        gain,
        work,
        false,
        zero(T),
        0,
        zero(T),
        nothing,
        false,
    )
    _sg_bind(st, populations)
    return st
end

function _sg_ring(st, stage, k, p)
    X = ntuple(d -> stage[k].X[d][p], 3)
    radius = st._radius[k][p]
    B = (zero(st.B0), zero(st.B0), st.B0)
    return (X, radius, B, _gyro_ring_geometry(X, radius, B))
end
function _sg_ring_point(ring, phase)
    point = _gyro_ring_point(ring..., phase...)
    all(isfinite, point) || throw(ArgumentError("gyrocentre ring point is not representable"))
    return point
end
_sg_stencil(point, g) =
    ntuple(d -> _stencil1d(TSC(), _periodic_cell_position(point[d], g.L[d], g.dx[d])), 3)

function _sg_deposit!(st::SlabGyrokineticPIC{T}, stage) where {T}
    out, g = st.work.rho, st.g
    fill!(out, zero(T))
    lost = false
    for (k, pset) in enumerate(st.owners), p in eachindex(pset.weight)
        ring = _sg_ring(st, stage, k, p)
        charge = pset.q * pset.weight[p]
        lost |= _deposition_product_lost(pset.q, pset.weight[p], charge)
        for phase in st._phases
            stencils = _sg_stencil(_sg_ring_point(ring, phase), g)
            for O in CartesianIndices((3, 3, 3))
                index = ntuple(d -> _wrapped_cell(stencils[d][1] + O[d] - 1, g.n[d]), 3)
                term = charge
                for d = 1:3
                    factor = stencils[d][2][O[d]]
                    next = term * factor
                    lost |= _deposition_product_lost(term, factor, next)
                    term = next
                end
                out[index...] += term
                lost |= !isfinite(out[index...])
            end
        end
    end
    volume = _from_exact(T, _sg_cell_volume(g))
    denominator = T(st.nring) * volume
    lost |= T(st.nring) != st.nring || !isfinite(denominator) || denominator < floatmin(T)
    if !lost
        for I in eachindex(out)
            original = out[I]
            out[I] /= denominator
            lost |= _deposition_range_lost(out[I]) || (!iszero(original) && iszero(out[I]))
        end
    end
    if lost
        Q = _exact_rational
        sums = fill(zero(Rational{BigInt}), g.n)
        for (k, pset) in enumerate(st.owners), p in eachindex(pset.weight)
            ring = _sg_ring(st, stage, k, p)
            charge = Q(pset.q) * Q(pset.weight[p])
            for phase in st._phases
                stencils = _sg_stencil(_sg_ring_point(ring, phase), g)
                for O in CartesianIndices((3, 3, 3))
                    index = ntuple(d -> _wrapped_cell(stencils[d][1] + O[d] - 1, g.n[d]), 3)
                    term = charge * prod(d -> Q(stencils[d][2][O[d]]), 1:3)
                    sums[index...] += term
                end
            end
        end
        divisor = st.nring * _sg_cell_volume(g)
        for I in eachindex(out)
            out[I] = _from_exact(T, sums[I] / divisor)
        end
    end
    all(isfinite, out) || throw(ArgumentError("gyroaveraged charge density is not representable"))
    return out
end

# Range-checked local TSC potential and derivative. Centered differences exactly
# annihilate constant represented potentials; cold arithmetic retains delta.
function _sg_sample(field, point, g::FourierGrid{3,T}) where {T}
    stencils = _sg_stencil(point, g)
    pairs = ntuple(
        d -> _tsc_derivative_pair(
            _periodic_cell_position(point[d], g.L[d], g.dx[d]),
            stencils[d][1],
        ),
        3,
    )
    value = zero(T)
    gradient = ntuple(_ -> zero(T), 3)
    lost = false
    for O in CartesianIndices((3, 3, 3))
        index = ntuple(d -> _wrapped_cell(stencils[d][1] + O[d] - 1, g.n[d]), 3)
        term = field[index...]
        for d = 1:3
            factor = stencils[d][2][O[d]]
            next = term * factor
            lost |= _deposition_product_lost(term, factor, next)
            term = next
        end
        next = value + term
        lost |=
            !isfinite(next) ||
            issubnormal(next) ||
            (iszero(next) && !iszero(value) && !iszero(term))
        value = next
    end
    for d = 1:3
        total = zero(T)
        others = filter(!=(d), (1, 2, 3))
        for a = 1:3, b = 1:3
            index = ntuple(
                j ->
                    j == d ? _wrapped_cell(stencils[j][1] + 1, g.n[j]) :
                    _wrapped_cell(stencils[j][1] + (j == others[1] ? a : b) - 1, g.n[j]),
                3,
            )
            left = Base.setindex(index, _wrapped_cell(stencils[d][1], g.n[d]), d)
            right = Base.setindex(index, _wrapped_cell(stencils[d][1] + 2, g.n[d]), d)
            dl, dr = field[left...] - field[index...], field[right...] - field[index...]
            u, v = pairs[d].left * dl, pairs[d].right * dr
            line = u + v
            lost |= _deposition_range_lost(line) || (iszero(line) && !iszero(u) && !iszero(v))
            term = line
            for factor in (stencils[others[1]][2][a], stencils[others[2]][2][b])
                next = term * factor
                lost |= _deposition_product_lost(term, factor, next)
                term = next
            end
            next = total + term
            lost |=
                any(x -> !isfinite(x) || issubnormal(x), (dl, dr, u, v, term, next)) ||
                (iszero(u) && !iszero(pairs[d].left) && !iszero(dl)) ||
                (iszero(v) && !iszero(pairs[d].right) && !iszero(dr)) ||
                (iszero(next) && !iszero(total) && !iszero(term))
            total = next
        end
        result = total / g.dx[d]
        lost |= !isfinite(result) || issubnormal(result) || (!iszero(total) && iszero(result))
        gradient = Base.setindex(gradient, result, d)
    end
    if lost
        Q = _exact_rational
        exact_value = zero(Rational{BigInt})
        for O in CartesianIndices((3, 3, 3))
            index = ntuple(d -> _wrapped_cell(stencils[d][1] + O[d] - 1, g.n[d]), 3)
            exact_value += Q(field[index...]) * prod(d -> Q(stencils[d][2][O[d]]), 1:3)
        end
        value = _from_exact(T, exact_value)
        gradient = ntuple(3) do d
            total = zero(Rational{BigInt})
            others = filter(!=(d), (1, 2, 3))
            delta = Q(pairs[d].delta)
            for a = 1:3, b = 1:3
                index = ntuple(
                    j ->
                        j == d ? _wrapped_cell(stencils[j][1] + 1, g.n[j]) :
                        _wrapped_cell(stencils[j][1] + (j == others[1] ? a : b) - 1, g.n[j]),
                    3,
                )
                left = Base.setindex(index, _wrapped_cell(stencils[d][1], g.n[d]), d)
                right = Base.setindex(index, _wrapped_cell(stencils[d][1] + 2, g.n[d]), d)
                line =
                    (delta - 1 // 2) * (Q(field[left...]) - Q(field[index...])) +
                    (delta + 1 // 2) * (Q(field[right...]) - Q(field[index...]))
                total += line * Q(stencils[others[1]][2][a]) * Q(stencils[others[2]][2][b])
            end
            _from_exact(T, total / Q(g.dx[d]))
        end
    end
    all(isfinite, (value, gradient...)) ||
        throw(ArgumentError("gyroaveraged potential sample/gradient is not representable"))
    return (value, gradient...)
end
function _sg_stream_add(state, sample)
    parts, exact = state
    if exact === nothing
        updated = _gyro_ring_add(parts, sample)
        return updated[4] ? (updated, nothing) :
               (parts, _gyro_exact_parts(parts) + _exact_rational(sample))
    end
    return parts, exact + _exact_rational(sample)
end
function _sg_gather(st::SlabGyrokineticPIC{T}, stage, k, p) where {T}
    sums = ntuple(_ -> ((zero(T), zero(T), 0, true), nothing), 4)
    ring = _sg_ring(st, stage, k, p)
    for phase in st._phases
        sample = _sg_sample(st.work.phi, _sg_ring_point(ring, phase), st.g)
        sums = (
            _sg_stream_add(sums[1], sample[1]),
            _sg_stream_add(sums[2], sample[2]),
            _sg_stream_add(sums[3], sample[3]),
            _sg_stream_add(sums[4], sample[4]),
        )
    end
    values = (
        _gyro_ring_mean(sums[1]..., st.nring, T),
        _gyro_ring_mean(sums[2]..., st.nring, T),
        _gyro_ring_mean(sums[3]..., st.nring, T),
        _gyro_ring_mean(sums[4]..., st.nring, T),
    )
    all(isfinite, values) ||
        throw(ArgumentError("gyroaveraged potential/gradient is not representable"))
    return values
end

function _sg_fft!(out, field, g::FourierGrid{3,T}, gain::T) where {T}
    magnitude = maximum(abs, field)
    isfinite(magnitude) || throw(ArgumentError("gyrokinetic FFT input must be finite"))
    safe = SpectralOperators._fft_safe_magnitude(T, length(field), gain)
    scale = SpectralOperators._fft_input_scale(magnitude, safe)
    SpectralOperators._copy_fft_input!(out, field, scale)
    g.plan * out
    all(isfinite, out) || throw(ArgumentError("gyrokinetic FFT data are not finite"))
    return scale
end
function _sg_symbol(st, I)
    Q = _exact_rational
    value = Q(st.chi) * (Q(st.g.kfull[1][I[1]])^2 + Q(st.g.kfull[2][I[2]])^2)
    (st.g.midx[2][I[2]] != 0 || st.g.midx[3][I[3]] != 0) && (value += Q(st.alpha))
    return value
end
function _sg_field_receipt(st::SlabGyrokineticPIC{T}, source_scale, phi_scale) where {T}
    Q = _exact_rational
    g, w = st.g, st.work
    N = length(w.rho)
    residual = source_norm = potential_norm = zero(Rational{BigInt})
    minimum_symbol = nothing
    for I in CartesianIndices(w.source)
        for part in (real, imag)
            phi = Q(part(w.phihat[I])) * Q(phi_scale)
            potential_norm = max(potential_norm, abs(phi) / N)
            all(==(1), Tuple(I)) && continue
            source = Q(part(w.source[I])) * Q(source_scale)
            symbol = _sg_symbol(st, I)
            minimum_symbol = minimum_symbol === nothing ? symbol : min(minimum_symbol, symbol)
            residual = max(residual, abs(symbol * phi - source) / N)
            source_norm = max(source_norm, abs(source) / N)
        end
    end
    tolerance = Q(st.field_atol) + Q(st.field_rtol) * source_norm
    residual <= tolerance ||
        throw(ArgumentError("gyrokinetic potential failed its fresh modal field residual"))
    gauge = max(abs(Q(real(w.phihat[1]))), abs(Q(imag(w.phihat[1])))) * Q(phi_scale) / N
    gauge_tolerance = Q(st.field_atol) / minimum_symbol + Q(st.field_rtol) * potential_norm
    gauge <= gauge_tolerance ||
        throw(ArgumentError("gyrokinetic potential failed its zero-mean gauge gate"))
    inventory = _sg_inventory(st.owners, g, st.e)
    amplitude = zero(Rational{BigInt})
    for i = 1:g.n[1]
        mean_phi =
            sum((Q(w.phi[i, j, k]) for j = 1:g.n[2], k = 1:g.n[3]); init = zero(Rational{BigInt})) /
            (g.n[2] * g.n[3])
        for j = 1:g.n[2], k = 1:g.n[3]
            response = Q(w.phi[i, j, k]) - mean_phi
            density = inventory.ne0 + Q(st.alpha) * response / Q(st.e)
            materialized = _from_exact(T, density)
            density > 0 && isfinite(materialized) && materialized > 0 || throw(
                ArgumentError(
                    "linear adiabatic electron density must remain finite and strictly positive",
                ),
            )
            amplitude = max(amplitude, abs(Q(st.e) * response / Q(st.Te)))
        end
    end
    mean_real = Q(real(w.source[1])) * Q(source_scale) / N
    mean_imag = Q(imag(w.source[1])) * Q(source_scale) / N
    dc_defect = mean_real - inventory.charge / _sg_volume(g)
    charge_defect = _sg_cell_volume(g) * _exact_float_sum(w.rho) - inventory.charge
    values = map(
        x -> _from_exact(T, x),
        (
            residual,
            tolerance,
            source_norm,
            gauge,
            gauge_tolerance,
            mean_real,
            mean_imag,
            dc_defect,
            charge_defect,
            amplitude,
        ),
    )
    any(isnan, values) && throw(ArgumentError("invalid gyrokinetic field diagnostic"))
    return SlabGyrokineticFieldReport{T}(values..., all(isfinite, values))
end
function _sg_field!(st::SlabGyrokineticPIC{T}) where {T}
    g, w = st.g, st.work
    source_scale = _sg_fft!(w.source, w.rho, g, st._gain)
    for I in CartesianIndices(g.cbuf)
        if all(==(1), Tuple(I))
            g.cbuf[I] = zero(Complex{T})
        else
            source, divisor = w.source[I], st._lambda[I]
            value = source / divisor
            # An active component lost before inverse rescaling cannot be
            # recovered by that scale. Reject this admitted-range boundary.
            all(part -> iszero(part(source)) || !iszero(part(value)), (real, imag)) ||
                throw(ArgumentError("gyrokinetic modal division lost a nonzero component"))
            g.cbuf[I] = value
        end
    end
    all(isfinite, g.cbuf) || throw(ArgumentError("gyrokinetic modal potential is not finite"))
    g.iplan * g.cbuf
    SpectralOperators._store_fft_output!(w.phi, g.cbuf, source_scale)
    all(isfinite, w.phi) || throw(ArgumentError("gyrokinetic potential is not representable"))
    phi_scale = _sg_fft!(w.phihat, w.phi, g, one(T))
    return _sg_field_receipt(st, source_scale, phi_scale)
end
function _sg_rhs!(out, st::SlabGyrokineticPIC{T}, stage) where {T}
    _sg_deposit!(st, stage)
    _sg_field!(st)
    for (k, pset) in enumerate(st.owners), p in eachindex(pset.weight)
        _, gx, gy, gz = _sg_gather(st, stage, k, p)
        rate = (
            _gyro_scaled_ratio((-gy,), (st.B0,)),
            _gyro_scaled_ratio((gx,), (st.B0,)),
            stage[k].vpar[p],
        )
        acceleration = _gyro_scaled_ratio((-pset.q, gz), (pset.m,))
        all(isfinite, (rate..., acceleration)) ||
            throw(ArgumentError("gyrocentre characteristic rate is not representable"))
        for d = 1:3
            out[k].X[d][p] = rate[d]
        end
        out[k].vpar[p] = acceleration
    end
    return nothing
end
function _sg_copy!(dest, source)
    for k in eachindex(dest)
        for d = 1:3
            copyto!(dest[k].X[d], source[k].X[d])
        end
        copyto!(dest[k].vpar, source[k].vpar)
    end
    return nothing
end
function _sg_affine!(dest, old, rates, op)
    for k in eachindex(dest), c = 1:4
        target = c == 4 ? dest[k].vpar : dest[k].X[c]
        initial = c == 4 ? old[k].vpar : old[k].X[c]
        vectors = map(r -> c == 4 ? r[k].vpar : r[k].X[c], rates)
        for p in eachindex(target)
            target[p] = op(initial[p], map(a -> a[p], vectors)...)
            isfinite(target[p]) || throw(ArgumentError("gyrocentre RK stage is not representable"))
        end
    end
    return nothing
end
function _sg_wrap!(stage, g)
    for pop in stage, d = 1:3, p in eachindex(pop.vpar)
        x = pop.X[d][p]
        isfinite(x) || throw(ArgumentError("gyrocentre endpoint must be finite before wrapping"))
        pop.X[d][p] = _periodic_coordinate(x, zero(x), g.L[d], g.L[d])
    end
    return nothing
end
function _sg_publish!(st, stage, rho, phi, time, step, h, report)
    for k in eachindex(st.owners)
        for d = 1:3, p in eachindex(st.owners[k].vpar)
            @inbounds st.owners[k].X[d][p] = stage[k].X[d][p]
        end
        for p in eachindex(st.owners[k].vpar)
            @inbounds st.owners[k].vpar[p] = stage[k].vpar[p]
        end
    end
    for I in eachindex(st.rho_gc)
        @inbounds st.rho_gc[I] = rho[I]
        @inbounds st.phi[I] = phi[I]
    end
    st.initialized = true
    st.time, st.step, st.last_dt, st.last_report = time, step, h, report
    return st
end
function init!(st::SlabGyrokineticPIC{T}, populations::AbstractVector = st.owners) where {T}
    _sg_acquire(st)
    try
        _sg_bind(st, populations)
        _sg_copy!(st.work.endpoint, populations)
        _sg_wrap!(st.work.endpoint, st.g)
        _sg_deposit!(st, st.work.endpoint)
        report = _sg_field!(st)
        return _sg_publish!(
            st,
            st.work.endpoint,
            st.work.rho,
            st.work.phi,
            zero(T),
            0,
            zero(T),
            report,
        )
    finally
        st.busy = false
    end
end
function step!(st::SlabGyrokineticPIC{T}, populations::AbstractVector, dt::Real) where {T}
    _sg_acquire(st)
    try
        _sg_bind(st, populations)
        _sg_boundary(st)
        h = _require_finite_real("gyrokinetic dt", dt, T)
        h >= 0 && (iszero(dt) || !iszero(h)) || throw(
            ArgumentError(
                "gyrokinetic dt must be nonnegative and nonzero inputs must survive conversion",
            ),
        )
        iszero(h) && return st
        time, count = _validated_step_clock(st.time, st.step, h)
        time > st.time || throw(ArgumentError("gyrokinetic clock must advance at stored precision"))
        w = st.work
        _sg_copy!(w.old, populations)
        _sg_rhs!(w.rates[1], st, w.old)
        _sg_affine!(w.trial, w.old, (w.rates[1],), _RK4Affine(h, 2))
        _sg_rhs!(w.rates[2], st, w.trial)
        _sg_affine!(w.trial, w.old, (w.rates[2],), _RK4Affine(h, 2))
        _sg_rhs!(w.rates[3], st, w.trial)
        _sg_affine!(w.trial, w.old, (w.rates[3],), _RK4Affine(h, 1))
        _sg_rhs!(w.rates[4], st, w.trial)
        _sg_affine!(w.endpoint, w.old, w.rates, _RK4Affine(h, 6))
        _sg_wrap!(w.endpoint, st.g)
        _sg_deposit!(st, w.endpoint)
        report = _sg_field!(st)
        return _sg_publish!(st, w.endpoint, w.rho, w.phi, time, count, h, report)
    finally
        st.busy = false
    end
end

"""Synchronous physical reservoirs and field receipt; no finite-step conservation theorem."""
function slab_gyrokinetic_diagnostics(st::SlabGyrokineticPIC{T}, populations = st.owners) where {T}
    _sg_acquire(st)
    try
        _sg_bind(st, populations)
        _sg_boundary(st)
        Q = _exact_rational
        Kpar = Kmu = momentum = zero(Rational{BigInt})
        number = Rational{BigInt}[]
        charge = Rational{BigInt}[]
        for p in st.owners
            n = _exact_float_sum(p.weight)
            push!(number, n)
            push!(charge, Q(p.q) * n)
            for i in eachindex(p.weight)
                Kpar += Q(p.weight[i]) * Q(p.m) * Q(p.vpar[i])^2 / 2
                Kmu += Q(p.weight[i]) * Q(p.mu[i]) * Q(st.B0)
                momentum += Q(p.weight[i]) * Q(p.m) * Q(p.vpar[i])
            end
        end
        scale = _sg_fft!(st.work.phihat, st.phi, st.g, one(T))
        Upol = Ue = zero(Rational{BigInt})
        for I in CartesianIndices(st.work.phihat)
            z = st.work.phihat[I]
            square = (Q(real(z))^2 + Q(imag(z))^2) * Q(scale)^2
            Upol += Q(st.chi) * (Q(st.g.kfull[1][I[1]])^2 + Q(st.g.kfull[2][I[2]])^2) * square
            (st.g.midx[2][I[2]] != 0 || st.g.midx[3][I[3]] != 0) && (Ue += Q(st.alpha) * square)
        end
        factor = _sg_cell_volume(st.g) / (2length(st.phi))
        Upol *= factor
        Ue *= factor
        pairing =
            _sg_cell_volume(st.g) *
            sum((Q(a) * Q(b) for (a, b) in zip(st.rho_gc, st.phi)); init = zero(Rational{BigInt})) /
            2
        total = Kpar + Kmu + Upol + Ue
        radius = maximum((r for rs in st._radius for r in rs); init = zero(T))
        k2 =
            maximum(_exact_rational(x)^2 for x in st.g.kfull[1]) +
            maximum(_exact_rational(y)^2 for y in st.g.kfull[2])
        coverage = _sqrt_from_exact(T, Q(radius)^2 * k2)
        reservoirs = map(x -> _from_exact(T, x), (Kpar, Kmu, Upol, Ue, total, Upol + Ue - pairing))
        displayed_momentum = _from_exact(T, momentum)
        resolved =
            all(isfinite, (reservoirs..., displayed_momentum, radius, coverage)) &&
            st.last_report.range_resolved
        return (
            time = st.time,
            step = st.step,
            last_dt = st.last_dt,
            phase = :synchronous_integer,
            velocity_time_offset = zero(T),
            number = number,
            charge = charge,
            electron_background_charge = -sum(charge),
            parallel_momentum = displayed_momentum,
            energy = (
                parallel = reservoirs[1],
                magnetic_moment = reservoirs[2],
                polarization = reservoirs[3],
                electron_response_free_energy = reservoirs[4],
                total = reservoirs[5],
                field_pairing_defect = reservoirs[6],
                range_resolved = all(isfinite, reservoirs),
            ),
            max_radius = radius,
            max_resolved_kperp_radius = coverage,
            last_field = st.last_report,
            range_resolved = resolved,
        )
    finally
        st.busy = false
    end
end
