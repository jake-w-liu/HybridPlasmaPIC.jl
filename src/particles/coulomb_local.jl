# Cell-local scheduling for the canonical Coulomb angle/rotation kernel.
# Species and RNG are caller-owned; every velocity update is staged here.

struct _CoulombCellWork{T}
    velocity::NTuple{3,Vector{T}}
    cell::Vector{Int}
    index::Vector{Int}
    offset::Vector{Int}
    cursor::Vector{Int}
end

"""
    CoulombCollisionWorkspace(species, grid)

Reusable CPU storage for cell-local [`collide_coulomb!`](@ref). The workspace
contains derived cell bins, pairing indices and staged velocities, and holds no
RNG or simulation state. Use one workspace per concurrent collision call.
Particle-count changes resize its buffers; species count and grid cell count
must match the constructor. A restarted calculation can recreate this scratch.
"""
struct CoulombCollisionWorkspace{T}
    populations::Vector{_CoulombCellWork{T}}
    ncells::Int
end

@inline _coulomb_index_arrays(w::_CoulombCellWork) = (w.cell, w.index, w.offset, w.cursor)

@inline function _require_coulomb_cpu_array(array)
    array isa StridedVector ||
        throw(ArgumentError("cell-local collide_coulomb! requires CPU strided vectors"))
    return nothing
end

function _validate_local_coulomb_species(species, g::FourierGrid{D,T}) where {D,T}
    if species isa AbstractVector
        axes(species, 1) == Base.OneTo(length(species)) ||
            throw(ArgumentError("collide_coulomb!: species must use one-based axes"))
    end
    T <: Union{Float32,Float64} ||
        throw(ArgumentError("cell-local collide_coulomb! requires Float32 or Float64"))
    for (a, ps) in enumerate(species)
        ps isa ParticleSet{D,T} || throw(
            ArgumentError("collide_coulomb!: species $a must match grid dimension and precision"),
        )
        _check_particle_storage_axes(ps)
        for array in ps.x
            _require_coulomb_cpu_array(array)
        end
        for array in ps.v
            _require_coulomb_cpu_array(array)
        end
        _require_coulomb_cpu_array(ps.weight)
        _require_coulomb_cpu_array(ps.id)
        _require_coulomb_cpu_array(ps.tag)
        _particle_components_mightalias(ps) &&
            throw(ArgumentError("collide_coulomb!: particle components must not overlap"))
        for b = 1:(a-1)
            _particle_sets_mightalias(ps, species[b]) &&
                throw(ArgumentError("collide_coulomb!: species must have independent storage"))
        end
        _require_finite_real("species charge", ps.q, T)
        _require_finite_positive_real("species mass", ps.m, T)
        _require_finite_particle_positions(ps)
        _require_finite_particle_velocities(ps, "collide_coulomb!")
        for wp in ps.weight
            (isfinite(wp) && wp >= zero(T)) ||
                throw(ArgumentError("collide_coulomb!: weights must be finite and non-negative"))
        end
    end
    return nothing
end

function CoulombCollisionWorkspace(
    species::Union{Tuple,AbstractVector},
    g::FourierGrid{D,T},
) where {D,T}
    _validate_local_coulomb_species(species, g)
    ncells = prod(g.n)
    populations = Vector{_CoulombCellWork{T}}(undef, length(species))
    for (a, ps) in enumerate(species)
        n = nparticles(ps)
        populations[a] = _CoulombCellWork(
            ntuple(_ -> Vector{T}(undef, n), 3),
            Vector{Int}(undef, n),
            Vector{Int}(undef, n),
            Vector{Int}(undef, ncells + 1),
            Vector{Int}(undef, ncells),
        )
    end
    return CoulombCollisionWorkspace(populations, ncells)
end

function _validate_coulomb_work_array(array, species, coulomb_log)
    for ps in species
        _mightalias_particle_storage(array, ps) &&
            throw(ArgumentError("collide_coulomb!: workspace must not alias particles"))
    end
    if coulomb_log isa AbstractArray
        Base.mightalias(array, coulomb_log) &&
            throw(ArgumentError("collide_coulomb!: workspace must not alias Coulomb logarithms"))
    end
    return nothing
end

@inline function _require_coulomb_independent(array, other)
    _particle_arrays_mightalias(array, other) &&
        throw(ArgumentError("collide_coulomb!: workspace arrays must not overlap"))
    return nothing
end

function _require_coulomb_work_independent(array, other::_CoulombCellWork)
    for buffer in other.velocity
        _require_coulomb_independent(array, buffer)
    end
    for buffer in _coulomb_index_arrays(other)
        _require_coulomb_independent(array, buffer)
    end
    return nothing
end

function _validate_coulomb_workspace(work, species, g::FourierGrid{D,T}, coulomb_log) where {D,T}
    work isa CoulombCollisionWorkspace{T} ||
        throw(ArgumentError("collide_coulomb!: workspace precision must match the grid"))
    length(work.populations) == length(species) ||
        throw(DimensionMismatch("collide_coulomb!: workspace species count does not match"))
    work.ncells == prod(g.n) ||
        throw(DimensionMismatch("collide_coulomb!: workspace grid cell count does not match"))
    for (a, buffers) in enumerate(work.populations)
        length(buffers.offset) == work.ncells + 1 && length(buffers.cursor) == work.ncells ||
            throw(DimensionMismatch("collide_coulomb!: workspace cell buffers have wrong lengths"))
        # Homogeneous groups avoid repeatedly boxing mixed-element tuples while
        # retaining every within-buffer, cross-population and live-state check.
        for c = 1:3
            array = buffers.velocity[c]
            _validate_coulomb_work_array(array, species, coulomb_log)
            for d = 1:(c-1)
                _require_coulomb_independent(array, buffers.velocity[d])
            end
            for b = 1:(a-1)
                _require_coulomb_work_independent(array, work.populations[b])
            end
        end
        indices = _coulomb_index_arrays(buffers)
        for c = 1:4
            array = indices[c]
            _validate_coulomb_work_array(array, species, coulomb_log)
            for buffer in buffers.velocity
                _require_coulomb_independent(array, buffer)
            end
            for d = 1:(c-1)
                _require_coulomb_independent(array, indices[d])
            end
            for b = 1:(a-1)
                _require_coulomb_work_independent(array, work.populations[b])
            end
        end
    end
    return nothing
end

@inline _coulomb_log_value(value::Real, a, b, ::Type{T}) where {T} = T(value)
@inline _coulomb_log_value(value::AbstractMatrix, a, b, ::Type{T}) where {T} = T(value[a, b])

function _validate_coulomb_logs(value, count::Int, ::Type{T}) where {T}
    if value isa Real
        _require_finite_nonnegative_real("coulomb_log", value, T)
    elseif value isa AbstractMatrix{<:Real}
        axes(value) == (Base.OneTo(count), Base.OneTo(count)) ||
            throw(DimensionMismatch("coulomb_log must be a species-by-species matrix"))
        for a = 1:count, b = 1:count
            v = _require_finite_nonnegative_real("coulomb_log", value[a, b], T)
            v == T(value[b, a]) || throw(ArgumentError("coulomb_log must be symmetric"))
        end
    else
        throw(ArgumentError("coulomb_log must be a real scalar or symmetric real matrix"))
    end
    return nothing
end

function _bin_coulomb_particles!(w::_CoulombCellWork, ps, g, strides)
    n = nparticles(ps)
    for array in w.velocity
        length(array) == n || resize!(array, n)
    end
    length(w.cell) == n || resize!(w.cell, n)
    length(w.index) == n || resize!(w.index, n)
    for c = 1:3
        copyto!(w.velocity[c], ps.v[c])
    end
    fill!(w.offset, 0)
    @inbounds for p = 1:n
        cell = iszero(ps.weight[p]) ? 0 : _particle_linear_cell(ps, g, p, strides)
        w.cell[p] = cell
        iszero(cell) || (w.offset[cell+1] += 1)
    end
    w.offset[1] = 1
    @inbounds for cell in eachindex(w.cursor)
        w.offset[cell+1] += w.offset[cell]
        w.cursor[cell] = w.offset[cell]
    end
    @inbounds for p = 1:n
        cell = w.cell[p]
        iszero(cell) && continue
        w.index[w.cursor[cell]] = p
        w.cursor[cell] += 1
    end
    return nothing
end

function _shuffle_coulomb_cell!(index, first::Int, last::Int, rng)
    @inbounds for k = last:-1:(first+1)
        j = rand(rng, first:k)
        index[k], index[j] = index[j], index[k]
    end
    return nothing
end

# Scale the masses before division; ma+mb can overflow despite finite fractions.
@inline function _coulomb_mass_fractions(ma::T, mb::T) where {T}
    if ma >= mb
        ratio = mb / ma
        fa = inv(one(T) + ratio)
        return fa, ratio * fa
    end
    ratio = ma / mb
    fb = inv(one(T) + ratio)
    return ratio * fb, fb
end

# Original physical factors accompany the scaled approximation through the
# variance-one decision. The masses and count/reuse factors must not be replaced
# by an already rounded reduced mass or quotient in that comparison.
struct _CoulombLocalRate{D,T}
    dt::T
    number_scale::T
    logarithm::T
    qa::T
    qb::T
    ma::T
    mb::T
    epsilon::T
    widths::NTuple{D,T}
    fraction::T
    power::Int
end

struct _CoulombLocalPairRate{D,T}
    channel::_CoulombLocalRate{D,T}
    count::Int
    reuse::Int
    triplet::Bool
end

@inline function _coulomb_rate_product(fraction, power, value, divide::Bool)
    part, shift = frexp(value)
    fraction = divide ? fraction / part : fraction * part
    fraction, extra = frexp(fraction)
    return fraction, power + (divide ? -shift : shift) + extra
end

function _coulomb_scaled_channel(dt, number_scale, logarithm, qa, qb, ma, mb, epsilon, widths)
    # Mantissas stay bounded while integer exponents carry every scale.
    lesser, greater = Float64(min(ma, mb)), Float64(max(ma, mb))
    correction = 1 + lesser / greater
    fraction, power = 0.5, 1
    for value in (
        Float64(dt),
        Float64(number_scale),
        Float64(logarithm),
        abs(Float64(qa)),
        abs(Float64(qa)),
        abs(Float64(qb)),
        abs(Float64(qb)),
        correction,
        correction,
    )
        fraction, power = _coulomb_rate_product(fraction, power, value, false)
    end
    for value in (8.0, Float64(pi), Float64(epsilon), Float64(epsilon), lesser, lesser)
        fraction, power = _coulomb_rate_product(fraction, power, value, true)
    end
    for width in widths
        fraction, power = _coulomb_rate_product(fraction, power, Float64(width), true)
    end
    return fraction, power
end

function _coulomb_channel_rate(
    pa,
    pb,
    g::FourierGrid{D,T},
    dt,
    epsilon,
    number_scale,
    logarithm,
) where {D,T}
    fraction, power =
        _coulomb_scaled_channel(dt, number_scale, logarithm, pa.q, pb.q, pa.m, pb.m, epsilon, g.dx)
    # Preserve the existing native coefficient and its rounding point.
    native, extra = frexp(T(fraction))
    return _CoulombLocalRate{D,T}(
        dt,
        number_scale,
        logarithm,
        pa.q,
        pb.q,
        pa.m,
        pb.m,
        epsilon,
        g.dx,
        native,
        power + extra,
    )
end

@inline function _coulomb_scaled_density(fraction, power, count, reuse, triplet)
    fraction, power = _coulomb_rate_product(fraction, power, Float64(count), false)
    fraction, power = _coulomb_rate_product(fraction, power, Float64(reuse), true)
    triplet && (power -= 1)
    return fraction, power
end

@inline function _coulomb_local_scaled_coefficient(rate::_CoulombLocalPairRate{D,T}) where {D,T}
    fraction, power = _coulomb_scaled_density(
        Float64(rate.channel.fraction),
        rate.channel.power,
        rate.count,
        rate.reuse,
        rate.triplet,
    )
    native, extra = frexp(T(fraction))
    return native, power + extra
end

@inline function _coulomb_scaled_weight_speed(fraction, power, wmax, wbar, speed)
    fraction, power = _coulomb_rate_product(fraction, power, Float64(wmax), false)
    fraction, power = _coulomb_rate_product(fraction, power, Float64(wbar), true)
    for _ = 1:3
        fraction, power = _coulomb_rate_product(fraction, power, Float64(speed), true)
    end
    return fraction, power
end

function _coulomb_scaled_variance(
    rate::_CoulombLocalPairRate{D,T},
    wmax::T,
    wbar,
    speed::T,
) where {D,T}
    native, power = _coulomb_local_scaled_coefficient(rate)
    fraction, power = _coulomb_scaled_weight_speed(Float64(native), power, wmax, wbar, speed)
    return T(fraction), power
end

function _coulomb_local_variance64(rate::_CoulombLocalPairRate, wmax, wbar, speed)
    channel = rate.channel
    fraction, power = _coulomb_scaled_channel(
        channel.dt,
        channel.number_scale,
        channel.logarithm,
        channel.qa,
        channel.qb,
        channel.ma,
        channel.mb,
        channel.epsilon,
        channel.widths,
    )
    fraction, power = _coulomb_scaled_density(fraction, power, rate.count, rate.reuse, rate.triplet)
    fraction, power = _coulomb_scaled_weight_speed(fraction, power, wmax, wbar, speed)
    return ldexp(fraction, power)
end

# Includes channel products/divisions, reduced-mass correction, native casts,
# count/reuse and the shared sampler's eight operations. The doubled margin
# converts the relative forward-error bound to an absolute guard near one.
@inline _coulomb_variance_roundoff(::_CoulombLocalPairRate{D}, ::Type{T}) where {D,T} =
    (64 + 4D) * eps(T)

# Machin's identity and alternating remainders give rational lower/upper bounds
# without changing or depending on process-wide MPFR precision/exponent limits.
# A rational cannot equal pi; refine until the comparison has a certificate.
@noinline function _coulomb_above_pi(value::Rational{BigInt})
    a, b = BigInt(0) // BigInt(1), BigInt(0) // BigInt(1)
    ta, tb = BigInt(1) // BigInt(5), BigInt(1) // BigInt(239)
    n = 0
    while true
        a += ta
        b += tb
        ta *= -(2n + 1) // (25 * (2n + 3))
        tb *= -(2n + 1) // (239^2 * (2n + 3))
        low = 16min(a, a + ta) - 4max(b, b + tb)
        high = 16max(a, a + ta) - 4min(b, b + tb)
        value <= low && return false
        value >= high && return true
        n += 1
    end
end

@noinline function _coulomb_large_angle_exact(
    rate::_CoulombLocalPairRate{D,T},
    wmax,
    wbar,
    speed,
) where {D,T}
    if T === Float32 && wmax isa Float32 && speed isa Float32 && wbar isa Union{Float32,Float64}
        # Reevaluate original factors, not the narrowed channel coefficient.
        # The same 32+D rounding-factor bound applies with binary64 precision;
        # only a certified side of one bypasses the rational/pi comparison.
        variance = _coulomb_local_variance64(rate, wmax, wbar, speed)
        abs(variance - 1.0) > _coulomb_variance_roundoff(rate, Float64) && return variance > 1.0
    end
    channel = rate.channel
    ma, mb = _exact_rational(channel.ma), _exact_rational(channel.mb)
    numerator =
        _exact_rational(channel.dt) *
        _exact_rational(channel.number_scale) *
        _exact_rational(channel.logarithm) *
        _exact_rational(channel.qa)^2 *
        _exact_rational(channel.qb)^2 *
        (ma + mb)^2 *
        _exact_rational(wmax) *
        rate.count
    denominator =
        8 *
        _exact_rational(channel.epsilon)^2 *
        ma^2 *
        mb^2 *
        _exact_rational(wbar) *
        _exact_rational(speed)^3 *
        rate.reuse
    rate.triplet && (denominator *= 2)
    for width in channel.widths
        denominator *= _exact_rational(width)
    end
    return _coulomb_above_pi(numerator / denominator)
end

@inline function _coulomb_local_pair!(
    wa::_CoulombCellWork{T},
    wb::_CoulombCellWork{T},
    pa,
    pb,
    i::Int,
    j::Int,
    rate::_CoulombLocalPairRate,
    uf::T,
    fa::T,
    fb::T,
    rng,
) where {T}
    @inbounds begin
        va, vb = wa.velocity, wb.velocity
        gx, gy, gz = va[1][i] - vb[1][j], va[2][i] - vb[2][j], va[3][i] - vb[3][j]
        (isfinite(gx) && isfinite(gy) && isfinite(gz)) ||
            throw(ArgumentError("collide_coulomb!: relative velocity exceeds numeric range"))
        speed = hypot(gx, gy, gz)
        isfinite(speed) ||
            throw(ArgumentError("collide_coulomb!: relative speed exceeds numeric range"))
        iszero(speed) && return nothing
        wi, wj = pa.weight[i], pb.weight[j]
        maximum_weight = max(wi, wj)
        fraction, power = _coulomb_local_scaled_coefficient(rate)
        sample = _coulomb_sample_increment(
            gx,
            gy,
            gz,
            speed,
            ldexp(fraction, power),
            maximum_weight,
            one(T),
            uf,
            2 * T(π),
            rng,
            rate,
        )
        accept_a = wi == wj || rand(rng, T) < wj / maximum_weight
        accept_b = wi == wj || rand(rng, T) < wi / maximum_weight
        sample.identity && return nothing
        original_a = (va[1][i], va[2][i], va[3][i])
        original_b = (vb[1][j], vb[2][j], vb[3][j])
        increment = sample.increment
        complete =
            sample.complete || any(
                c -> _coulomb_complete_needed(original_a[c], original_b[c], increment[c], fa, fb),
                1:3,
            )
        if complete
            newa, newb = _coulomb_complete_update(original_a, original_b, pa.m, pb.m, speed, sample)
            ((!accept_a || all(isfinite, newa)) && (!accept_b || all(isfinite, newb))) ||
                throw(ArgumentError("collide_coulomb!: scattered velocity must be finite"))
            for c = 1:3
                accept_a && (va[c][i] = newa[c])
                accept_b && (vb[c][j] = newb[c])
            end
            return nothing
        end
        for c = 1:3
            na, nb = _coulomb_mass_update(va[c][i], vb[c][j], pa.m, pb.m, fa, fb, increment[c])
            ((!accept_a || isfinite(na)) && (!accept_b || isfinite(nb))) ||
                throw(ArgumentError("collide_coulomb!: scattered velocity must be finite"))
            accept_a && (va[c][i] = na)
            accept_b && (vb[c][j] = nb)
        end
    end
    return nothing
end

function _coulomb_self_cell!(
    work,
    ps,
    first::Int,
    last::Int,
    channel::_CoulombLocalRate{D,T},
    uf::T,
    rng,
) where {D,T}
    n = last - first + 1
    n < 2 && return nothing
    _shuffle_coulomb_cell!(work.index, first, last, rng)
    rate = _CoulombLocalPairRate(channel, n - 1, 1, false)
    start = first
    if isodd(n)
        i, j, k = work.index[first], work.index[first+1], work.index[first+2]
        for (a, b) in ((i, j), (j, k), (k, i))
            _coulomb_local_pair!(
                work,
                work,
                ps,
                ps,
                a,
                b,
                _CoulombLocalPairRate(channel, n - 1, 1, true),
                uf,
                T(0.5),
                T(0.5),
                rng,
            )
        end
        start += 3
    end
    @inbounds for k = start:2:last
        _coulomb_local_pair!(
            work,
            work,
            ps,
            ps,
            work.index[k],
            work.index[k+1],
            rate,
            uf,
            T(0.5),
            T(0.5),
            rng,
        )
    end
    return nothing
end

function _coulomb_cross_cell!(
    wa,
    wb,
    pa,
    pb,
    cell::Int,
    channel::_CoulombLocalRate{D,T},
    uf::T,
    fa::T,
    fb::T,
    rng,
) where {D,T}
    firsta, lasta = wa.offset[cell], wa.offset[cell+1] - 1
    firstb, lastb = wb.offset[cell], wb.offset[cell+1] - 1
    na, nb = lasta - firsta + 1, lastb - firstb + 1
    (iszero(na) || iszero(nb)) && return nothing
    _shuffle_coulomb_cell!(wa.index, firsta, lasta, rng)
    _shuffle_coulomb_cell!(wb.index, firstb, lastb, rng)
    larger, smaller = max(na, nb), min(na, nb)
    repeats, extra = divrem(larger, smaller)
    @inbounds for k = 0:(larger-1)
        small = mod(k, smaller)
        reuse = repeats + (small < extra)
        i = wa.index[firsta+(na >= nb ? k : small)]
        j = wb.index[firstb+(nb >= na ? k : small)]
        _coulomb_local_pair!(
            wa,
            wb,
            pa,
            pb,
            i,
            j,
            _CoulombLocalPairRate(channel, larger, reuse, false),
            uf,
            fa,
            fb,
            rng,
        )
    end
    return nothing
end

"""
    collide_coulomb!(species, grid, dt; coulomb_log, epsilon0, weight_to_number,
                     u_floor=1e-3, rng=Random.default_rng(), workspace=nothing)

Nonrelativistic, cell-local self and interspecies Coulomb scattering on a
periodic `FourierGrid`. All populations must use its dimension and Float32/64
precision with independent CPU strided-vector storage and one-based species
collection axes. Only positive-weight
particles pair. A required nonnegative scalar or symmetric matrix `coulomb_log`
sets each unordered species channel; zero entries and neutral charges disable
channels. Masses, `epsilon0`, `weight_to_number` and `u_floor` must be finite and
positive. Charges, positions and velocities must be finite; weights and `dt`
must be finite and nonnegative.
`u_floor` has the same velocity units as the particles and must be rescaled
when converting a calculation between unit systems.

`weight_to_number` converts stored macro-weight into real-particle counts in
the chosen units. Permittivity and this conversion are explicit: normalized
density weights are not physical counts. Use [`coulomb_normalization`](@ref)
with `PlasmaUnits` for the package's normalized hybrid convention. In D<3,
cell volume has unit transverse extent in the chosen coordinates.

The physical reduced mass controls each pair's center-of-mass rotation. Pair
variance is `(qa*qb)^2*logLambda*dt/(8π*epsilon0^2*mu^2*max(abs(g),u_floor)^3)`
times the weighted partner density and pairing correction. Self collisions
exclude self density (`N-1` partners), with half-strength odd-count triplets.
Unequal cross counts correct every smaller-list particle's actual reuse count.
The prefactor uses mathematical π. Binary scaling preserves physical rate
factors across intermediate range limits; near variance one, original masses,
integer counts and certified π bounds determine the scattering regime.
The Higginson–Holod–Link weight rejection preserves momentum and energy in
expectation; equal weights conserve them pairwise to roundoff. The existing
Takizuka–Abe angle law and its isotropic strong-scattering branch are shared
with the 0-D method. Resolve collision timesteps; finite-step accuracy is not
implied by conservation.

All velocities stage until the complete operation succeeds. Parameter or
sampling errors leave populations unchanged; RNG and workspace may advance.
Positions and particle metadata remain unchanged. Identical species order,
initial state, RNG state and environment reproduce the result. Existing
simulation checkpoints do not save an external RNG. Returns `species`.
"""
function collide_coulomb!(
    species::Union{Tuple,AbstractVector},
    g::FourierGrid{D,T},
    dt::Real;
    coulomb_log,
    epsilon0::Real,
    weight_to_number::Real,
    u_floor::Real = 1e-3,
    rng = Random.default_rng(),
    workspace = nothing,
) where {D,T}
    dtT = _require_finite_nonnegative_real("dt", dt, T)
    epsilon = _require_finite_positive_real("epsilon0", epsilon0, T)
    number_scale = _require_finite_positive_real("weight_to_number", weight_to_number, T)
    uf = _require_finite_positive_real("u_floor", u_floor, T)
    _validate_local_coulomb_species(species, g)
    _validate_coulomb_logs(coulomb_log, length(species), T)
    for dx in g.dx
        _require_finite_positive_real("grid cell width", dx, T)
    end
    workspace === nothing || _validate_coulomb_workspace(workspace, species, g, coulomb_log)
    (isempty(species) || iszero(dtT)) && return species
    active = false
    for a in eachindex(species), b = a:length(species)
        active |=
            !iszero(species[a].q) &&
            !iszero(species[b].q) &&
            !iszero(_coulomb_log_value(coulomb_log, a, b, T))
    end
    active || return species
    work = workspace === nothing ? CoulombCollisionWorkspace(species, g) : workspace
    strides = _particle_cell_strides(g)
    for a in eachindex(species)
        _bin_coulomb_particles!(work.populations[a], species[a], g, strides)
    end
    for a in eachindex(species), b = a:length(species)
        pa, pb = species[a], species[b]
        logarithm = _coulomb_log_value(coulomb_log, a, b, T)
        (iszero(pa.q) || iszero(pb.q) || iszero(logarithm)) && continue
        channel = _coulomb_channel_rate(pa, pb, g, dtT, epsilon, number_scale, logarithm)
        fa, fb = _coulomb_mass_fractions(pa.m, pb.m)
        wa, wb = work.populations[a], work.populations[b]
        for cell = 1:work.ncells
            if a == b
                _coulomb_self_cell!(
                    wa,
                    pa,
                    wa.offset[cell],
                    wa.offset[cell+1] - 1,
                    channel,
                    uf,
                    rng,
                )
            else
                _coulomb_cross_cell!(wa, wb, pa, pb, cell, channel, uf, fa, fb, rng)
            end
        end
    end
    for a in eachindex(species), c = 1:3
        copyto!(species[a].v[c], work.populations[a].velocity[c])
    end
    return species
end
