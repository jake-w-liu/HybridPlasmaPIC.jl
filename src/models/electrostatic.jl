# Periodic electrostatic PIC with a shared spectral Poisson field. Scalar states
# retain the normalized electron/background API; dimension-parametric collection
# states support arbitrary kinetic charges and masses with a common leapfrog phase.

"""
    Electrostatic1D(g, Nparticles; n0=1.0)

State for a 1D electrostatic PIC on periodic grid `g`: electric field, electron
number density, a per-particle gather buffer, and the immobile-ion background
density `n0` (charge `ρ = n0 − n_e`).
"""
struct Electrostatic1D{T,G}
    g::G
    n0::T
    E::Vector{T}
    ne::Vector{T}
    Ep::Vector{T}
    primed::Base.RefValue{Bool}        # false until the one-time leapfrog velocity priming
    last_dt::Base.RefValue{T}          # preceding drift interval sets the velocity half level
end

function Electrostatic1D(g::FourierGrid{1,T}, Nparticles::Integer; n0 = 1.0) where {T}
    n = g.n[1]
    Np = _particle_length(Nparticles)
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    Electrostatic1D{T,typeof(g)}(
        g,
        n0T,
        zeros(T, n),
        zeros(T, n),
        zeros(T, Np),
        Ref(false),
        Ref(zero(T)),
    )
end

# Collection work owns proposed particles and fields. Public particles/fields are
# published only after the complete collection step succeeds; scratch may change.
const _ElectrostaticHostParticle{D,T} =
    ParticleSet{D,T,Vector{T},Vector{T},Vector{T},Vector{UInt64},Vector{UInt32}}
struct _ElectrostaticCollectionWork{D,T,P<:ParticleSet{D,T}}
    neutrality::Symbol
    charge_rtol::T
    owners::Vector{P}
    staged::Vector{_ElectrostaticHostParticle{D,T}}
    E::NTuple{3,Array{T,D}}
    ne::Array{T,D}
    density::Array{T,D}
    initialized::Base.RefValue{Bool}
end

"""
    ElectrostaticPIC(g, Nparticles; n0=1.0)

Dimension-parametric electrostatic PIC on a periodic `FourierGrid{D}`. The state
stores a 3-component electric field tuple over the `D`-dimensional mesh, electron
number density, per-particle gather buffers, and an immobile ion background
density `n0` (charge `ρ = n0 - n_e`). Components `D+1:3` of `E` are zero.
"""
struct ElectrostaticPIC{D,T,G,W<:Union{Nothing,_ElectrostaticCollectionWork{D,T}}}
    g::G
    n0::T
    E::NTuple{3,Array{T,D}}
    ne::Array{T,D}
    Ep::NTuple{3,Vector{T}}
    primed::Base.RefValue{Bool}        # false until the one-time leapfrog velocity priming
    last_dt::Base.RefValue{T}          # preceding drift interval sets the velocity half level
    collection::W
end

function ElectrostaticPIC(g::FourierGrid{D,T}, Nparticles::Integer; n0 = 1.0) where {D,T}
    1 <= D <= 3 || throw(ArgumentError("ElectrostaticPIC supports D = 1, 2, or 3"))
    Np = _particle_length(Nparticles)
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    return ElectrostaticPIC{D,T,typeof(g),Nothing}(
        g,
        n0T,
        ntuple(_ -> zeros(T, g.n), 3),
        zeros(T, g.n),
        ntuple(_ -> zeros(T, Np), 3),
        Ref(false),
        Ref(zero(T)),
        nothing,
    )
end

function _resize_espic_particle_workspaces!(es::Electrostatic1D, N::Integer)
    length(es.Ep) == N || resize!(es.Ep, N)
    return es
end

function _resize_espic_particle_workspaces!(es::ElectrostaticPIC, N::Integer)
    for c = 1:3
        length(es.Ep[c]) == N || resize!(es.Ep[c], N)
    end
    return es
end

"Solve −∂²φ = ρ/ε0, E = −∂φ, spectrally: Ê_m = −i ρ̂_m/(ε0 k_m), Ê_0 = 0."
function poisson_E!(es::Electrostatic1D{T}) where {T}
    g = es.g
    n = g.n[1]
    cb = g.cbuf
    @inbounds for i = 1:n
        cb[i] = Complex{T}(es.n0 - es.ne[i])
    end
    g.plan * cb
    k = g.kvec[1]
    @inbounds for m = 1:n
        km = k[m]
        cb[m] = km == 0 ? zero(Complex{T}) : (-im * cb[m]) / km
    end
    g.iplan * cb
    @inbounds for i = 1:n
        es.E[i] = real(cb[i])
    end
    return es.E
end

# Longitudinal Fourier modes share exact-zero classification and one cold
# complete-expression projection. FFT input/output range remains owned by the
# caller: this path starts with the represented finite Fourier coefficients.
@inline function _spectral_wavevector(g::FourierGrid{D,T}, I::CartesianIndex{D}) where {D,T}
    return ntuple(d -> @inbounds(g.kvec[d][I[d]]), Val(D))
end

# A subnormal intermediate can lose bits that a later large factor restores.
# Keep the original operation order when every intermediate retains its range.
@inline _spectral_product_lost_range(a::T, b::T, p::T) where {T} =
    !isfinite(p) || (!iszero(a) && !iszero(b) && abs(p) < floatmin(T))
@inline _spectral_quotient_lost_range(a::T, b::T, p::T) where {T} =
    !isfinite(p) || (!iszero(a) && abs(p) < floatmin(T))
@inline _spectral_mode_is_null(k) = all(iszero, k)
function _spectral_norm2(k::NTuple{D,T}) where {D,T}
    total = zero(T)
    lost = false
    for kd in k
        square = kd^2
        lost |= _spectral_product_lost_range(kd, kd, square)
        total += square
    end
    return total, lost || !isfinite(total)
end
_spectral_exact_complex(z::Complex) = complex(_exact_rational(real(z)), _exact_rational(imag(z)))
function _spectral_exact_projection(
    a::NTuple{D,Complex{T}},
    k::NTuple{D,T},
    source::Complex{T},
    divisor::T,
) where {D,T}
    # The homogeneous zero mode needs no exact allocations, even at extreme k.
    iszero(source) && all(iszero, a) && return a
    kr = map(_exact_rational, k)
    ar = map(_spectral_exact_complex, a)
    k2 = sum(x -> x * x, kr)
    # Source/divisor and the original components remain unrounded until the
    # complete longitudinal correction is formed. In particular, dt is not
    # consumed by an overflowing/underflowing target before reaching here.
    defect =
        _spectral_exact_complex(source) / _exact_rational(divisor) -
        im * sum(kr[d] * ar[d] for d = 1:D)
    return map(ar, kr) do ad, kd
        result = ad - im * kd * defect / k2
        complex(_from_exact(T, real(result)), _from_exact(T, imag(result)))
    end
end
# Preserve callers' existing nonfinite-input behavior; the exact path accepts
# only finite represented operands and a nonzero divisor.
@inline _spectral_projection_eligible(a, k, source, divisor) =
    all(isfinite, a) &&
    all(isfinite, k) &&
    isfinite(source) &&
    isfinite(divisor) &&
    !iszero(divisor)

@inline function _spectral_longitudinal_current(source::Complex{T}, k::T, dt::T) where {T}
    # Legacy continuity uses -source/(dt*ik) on a nonzero derivative mode.
    # Keep its ordinary two-quotient order, retaining original operands when
    # either intermediate loses range before the completed modal coefficient.
    scaled = source / dt
    result = -scaled / Complex{T}(0, k)
    # Each source component follows its scalar quotient through the pure
    # imaginary divisor. A nonzero component must remain normal at both stages.
    real_normal = (abs(real(scaled)) >= floatmin(T)) & (abs(imag(result)) >= floatmin(T))
    imaginary_normal = (abs(imag(scaled)) >= floatmin(T)) & (abs(real(result)) >= floatmin(T))
    safe =
        isfinite(scaled) &
        isfinite(result) &
        (iszero(real(source)) | real_normal) &
        (iszero(imag(source)) | imaginary_normal)
    zero_mode = (zero(Complex{T}),)
    cold = !safe && _spectral_projection_eligible(zero_mode, (k,), -source, dt)
    return cold ? only(_spectral_exact_projection(zero_mode, (k,), -source, dt)) : result
end

function _spectral_project_mode(
    a::NTuple{D,Complex{T}},
    k::NTuple{D,T},
    source::Complex{T},
    divisor::T,
) where {D,T}
    _spectral_mode_is_null(k) && return a
    k2, lost = _spectral_norm2(k)
    divergence = zero(Complex{T})
    for d = 1:D
        term = Complex{T}(0, k[d]) * a[d]
        lost |=
            _spectral_product_lost_range(k[d], imag(a[d]), real(term)) ||
            _spectral_product_lost_range(k[d], real(a[d]), imag(term))
        divergence += term
        lost |= !isfinite(divergence)
    end
    target = source / divisor
    lost |=
        _spectral_quotient_lost_range(real(source), divisor, real(target)) ||
        _spectral_quotient_lost_range(imag(source), divisor, imag(target))
    defect = target - divergence
    lost |= !isfinite(defect)
    coefficients = map(kd -> (-im * kd / k2), k)
    for d = 1:D
        lost |= _spectral_quotient_lost_range(k[d], k2, imag(coefficients[d]))
    end
    corrections = map(coefficient -> coefficient * defect, coefficients)
    for d = 1:D
        lost |=
            _spectral_product_lost_range(
                imag(coefficients[d]),
                imag(defect),
                real(corrections[d]),
            ) ||
            _spectral_product_lost_range(imag(coefficients[d]), real(defect), imag(corrections[d]))
    end
    result = map(+, a, corrections)
    lost |= !all(isfinite, result)
    cold = lost && _spectral_projection_eligible(a, k, source, divisor)
    return cold ? _spectral_exact_projection(a, k, source, divisor) : result
end
function _spectral_poisson_mode(k::NTuple{D,T}, source::Complex{T}) where {D,T}
    zeros_mode = map(_ -> zero(Complex{T}), k)
    _spectral_mode_is_null(k) && return zeros_mode
    k2, lost = _spectral_norm2(k)
    numerator = map(kd -> (-im * kd * source), k)
    result = map(value -> value / k2, numerator)
    for d = 1:D
        lost |=
            _spectral_product_lost_range(k[d], imag(source), real(numerator[d])) ||
            _spectral_product_lost_range(k[d], real(source), imag(numerator[d]))
        lost |=
            _spectral_quotient_lost_range(real(numerator[d]), k2, real(result[d])) ||
            _spectral_quotient_lost_range(imag(numerator[d]), k2, imag(result[d]))
    end
    cold = lost && _spectral_projection_eligible(zeros_mode, k, source, one(T))
    return cold ? _spectral_exact_projection(zeros_mode, k, source, one(T)) : result
end

"""
    poisson_E!(es::ElectrostaticPIC)

Solve `-∇²φ = ρ`, `E = -∇φ` spectrally in 1D, 2D, or 3D.
For a collection state, refresh the charge from its initialized populations,
validate the configured physical neutrality policy, and publish the staged field
without changing the leapfrog phase. `ne` is derived from those populations.
Scalar electron states instead solve from the caller-supplied `n0-ne` array.
"""
function poisson_E!(es::ElectrostaticPIC)
    if es.collection === nothing
        return _espic_poisson_E!(es.E, es.ne, es.n0, es.g)
    end
    work = _espic_collection(es)
    work.initialized[] || throw(
        ArgumentError(
            "initialize the electrostatic species collection before refreshing its field",
        ),
    )
    _espic_validate_species(es, work.owners)
    _espic_collection_fields!(es, work.owners)
    _espic_publish_fields!(es, work)
    return es.E
end

function _espic_poisson_E!(E, ne, n0::T, g::FourierGrid{D,T}) where {D,T}
    cb = g.cbuf
    @inbounds for I in CartesianIndices(cb)
        cb[I] = Complex{T}(n0 - ne[I])
    end
    g.plan * cb
    for c = 1:3
        fill!(E[c], zero(T))
    end
    buffers = (g.tbuf, g.abuf, g.sbuf)
    @inbounds for I in CartesianIndices(cb)
        values = _spectral_poisson_mode(_spectral_wavevector(g, I), cb[I])
        for d = 1:D
            buffers[d][I] = values[d]
        end
    end
    for d = 1:D
        g.iplan * buffers[d]
        for I in CartesianIndices(cb)
            E[d][I] = real(buffers[d][I])
        end
    end
    return E
end

@inline function _require_espic_electrons(e::ParticleSet{D,T}) where {D,T}
    q = _require_finite_real("electron charge q", e.q, T)
    m = _require_finite_real("electron mass m", e.m, T)
    q == -one(T) ||
        throw(ArgumentError("electrostatic PIC requires electron ParticleSet with q = -1"))
    m == one(T) ||
        throw(ArgumentError("electrostatic PIC requires electron ParticleSet with m = 1"))
    return nothing
end

"Deposit electron density and solve the field (call once after loading)."
function init_espic!(es::Electrostatic1D{T}, e::ParticleSet{1,T}) where {T}
    _require_espic_electrons(e)
    _resize_espic_particle_workspaces!(es, nparticles(e))
    density!(es.ne, e, es.g, CIC())
    poisson_E!(es)
    es.primed[] = false                 # re-arm the one-time priming for the next step
    es.last_dt[] = zero(T)
    return es
end

function init_espic!(es::ElectrostaticPIC{D,T}, e::ParticleSet{D,T}) where {D,T}
    _espic_require_scalar(es)
    _require_espic_electrons(e)
    _resize_espic_particle_workspaces!(es, nparticles(e))
    density!(es.ne, e, es.g, CIC())
    poisson_E!(es)
    es.primed[] = false                 # re-arm the one-time priming for the next step
    es.last_dt[] = zero(T)
    return es
end

"""
    step_espic!(es, electrons, dt)

Advance one leapfrog step: gather E → Boris kick (B=0) → drift → periodic wrap →
deposit density → solve Poisson. When `dt` changes, the velocity kick spans
`(last_dt + dt)/2` between the preceding and next half levels; positions advance
by `dt`. A zero timestep leaves both the state and its time staggering unchanged.
A nonzero `dt` below `2*floatmin(T)` throws an `ArgumentError` before any state
changes. Priming, kicks and drifts follow the normwise range contract of
[`boris_kick`](@ref); a primed velocity that exceeds `T` throws an `ArgumentError`.
"""
function step_espic!(es::Electrostatic1D{T}, e::ParticleSet{1,T}, dt::Real) where {T}
    _require_espic_electrons(e)
    dtT =
        _require_resolved_dt("step_espic!", _validated_nonnegative_dt(T, dt; name = "step_espic!"))
    iszero(dtT) && return es            # dt=0 no-op: do not consume the one-time priming
    _resize_espic_particle_workspaces!(es, nparticles(e))
    g = es.g
    gather_scalar!(es.Ep, es.E, e, g, CIC())
    _espic_push!(e, (es.Ep,), g, dtT, es.primed[], es.last_dt[])
    es.primed[] = true
    density!(es.ne, e, g, CIC())
    poisson_E!(es)
    es.last_dt[] = dtT
    return es
end

function step_espic!(es::ElectrostaticPIC{D,T}, e::ParticleSet{D,T}, dt::Real) where {D,T}
    _espic_require_scalar(es)
    _require_espic_electrons(e)
    dtT =
        _require_resolved_dt("step_espic!", _validated_nonnegative_dt(T, dt; name = "step_espic!"))
    iszero(dtT) && return es            # dt=0 no-op: do not consume the one-time priming
    _resize_espic_particle_workspaces!(es, nparticles(e))
    g = es.g
    gather_vector!(es.Ep, es.E, e, g, CIC())
    _espic_push!(e, es.Ep, g, dtT, es.primed[], es.last_dt[])
    es.primed[] = true
    density!(es.ne, e, g, CIC())
    poisson_E!(es)
    es.last_dt[] = dtT
    return es
end

"Electrostatic field energy ∫ ½ E² dx, retaining finite complete-expression range."
field_energy(es::Electrostatic1D) = _espic_field_energy((es.E,), es.g.dx)
field_energy(es::ElectrostaticPIC) = _espic_field_energy(es.E, es.g.dx)

function _espic_field_energy(E, dx::NTuple{D,T}) where {D,T}
    if T === Float32
        value = _field_energy_float32(E, (), 1.0, dx)
        value === nothing || return value
    end
    # Retain the scalar/dimension-parametric historical summation order. The
    # shared field-energy finisher repairs only lost range for Float64.
    s = sum(sum(abs2, component) for component in E)
    minimum_e2 = T(Inf)
    for component in E, value in component
        iszero(value) || (minimum_e2 = min(minimum_e2, value^2))
    end
    return _em_field_energy_finish(E, (), one(T), dx, s, minimum_e2, T(Inf))
end

# Embed the longitudinal compatibility buffer in the physical three-vector.
# Separate tuple arities make every accessed field component explicit.
@inline _espic_particle_field(Ep::Tuple{V}, p, ::Type{T}) where {V,T} = (Ep[1][p], zero(T), zero(T))
@inline _espic_particle_field(Ep::NTuple{3}, p, ::Type{T}) where {T} =
    (Ep[1][p], Ep[2][p], Ep[3][p])

# The scalar electron and arbitrary-population paths share exactly this kick,
# priming and drift. The 1D compatibility state carries only longitudinal E.
# Priming, kick and drift keep their native order for ordinary inputs and follow
# the normwise range contract of kinematics_range.jl; a primed velocity that
# exceeds T throws an ArgumentError.
function _espic_push!(ps::ParticleSet{D,T}, Ep::NTuple{C}, g, dt::T, primed, last_dt) where {D,T,C}
    qm = ps.q / ps.m
    kick_dt = primed ? last_dt / 2 + dt / 2 : dt
    indices = eachindex(ps.weight)
    if !primed
        h = dt / 2
        hq = h * qm
        if _es_prime_batch_certificate(ps.v, Ep, h, qm, hq, C, indices)
            _espic_prime!(ps, Ep, h, qm, hq, Val(false))
        else
            _espic_prime!(ps, Ep, h, qm, hq, Val(true))
        end
    end
    z = zero(T)
    hk = qm * kick_dt / 2
    velocities = ntuple(c -> c <= C ? ps.v[c] : z, Val(3))
    scales = _boris_batch_scales(velocities, _espic_kick_fields(Ep, z), (z, z, z), hk, indices)
    kick, vbound = _boris_batch_certificate(hk, qm, kick_dt, scales...)
    if kick && _coordinate_batch_certified(ps, _coordinate_batch_max(ps.x, T), dt, vbound)
        _espic_kick_drift!(ps, Ep, qm, dt, kick_dt, Val(false))
    else
        _espic_kick_drift!(ps, Ep, qm, dt, kick_dt, Val(true))
    end
    _require_finite_particle_velocities(ps, "electrostatic stepping")
    # The boundary owner validates all proposed positions before wrapping any.
    apply_periodic!(ps, ntuple(_ -> zero(T), D), g.L)
    return ps
end

# Field components of the kick as a three-vector of arrays (absent longitudinal-
# compatibility components are zero, as in the kick loop).
@inline _espic_kick_fields(Ep::Tuple{V}, z) where {V} = (Ep[1], z, z)
@inline _espic_kick_fields(Ep::NTuple{3}, z) = Ep

function _espic_prime!(
    ps::ParticleSet{D,T},
    Ep::NTuple{C},
    h,
    qm,
    hq,
    ::Val{guard},
) where {D,T,C,guard}
    @inbounds for c = 1:C, p in eachindex(ps.weight)
        ps.v[c][p] =
            guard ? _es_prime_component(ps.v[c][p], Ep[c][p], h, qm, hq) :
            ps.v[c][p] - hq * Ep[c][p]
    end
    return ps
end

function _espic_kick_drift!(
    ps::ParticleSet{D,T},
    Ep::NTuple{C},
    qm,
    dt,
    kick_dt,
    guard::Val,
) where {D,T,C}
    @inbounds for p in eachindex(ps.weight)
        # Boris with B=0 retains the established two-half-kick rounding order.
        v = ntuple(c -> c <= C ? ps.v[c][p] : zero(T), 3)
        E = _espic_particle_field(Ep, p, T)
        moved = _kick_eval(guard, v..., E..., zero(T), zero(T), zero(T), qm, kick_dt)
        for c = 1:C
            ps.v[c][p] = moved[c]
        end
        for d = 1:D
            ps.x[d][p] = _coordinate_eval(guard, ps.x[d][p], ps.v[d][p], dt)
        end
    end
    return ps
end

_espic_require_scalar(es::ElectrostaticPIC) =
    es.collection === nothing ? nothing :
    throw(ArgumentError("collection ElectrostaticPIC requires the ordered species vector"))

"""
    ElectrostaticPIC(g, species::AbstractVector; background_charge=0,
                     neutrality=:require, charge_rtol=64eps(T))

Periodic CPU electrostatic PIC for ordered `ParticleSet{D,T}` populations, with
finite arbitrary charges, positive masses, representable nonzero q/m for charged
populations, and independent strided storage.
`T` is Float32 or Float64. All populations share one electric field and one
leapfrog time level; call `init_espic!(es, species)` before stepping.

`background_charge` is a fixed uniform charge density (`ε0=1`). `:require`
rejects total physical charge larger than `charge_rtol` times the sum of absolute
particle/background charges. The default tolerance covers floating-point
summation/deposition roundoff. `:neutralize` explicitly adds a uniform background
that cancels the mean source; `electrostatic_diagnostics` reports that charge.
Neither policy alters particle weights. The integer-count electron constructor
retains its historical neutralizing-background behavior.

`ne` stores electron-equivalent charge density `-sum(q_s*n_s)`, so `n0-ne` is
the total physical source. The collection and order are bound by initialization;
call `init_espic!` to replace/reorder populations. Initialization resets the phase
without converting velocities: supply physical full-time velocities when calling
it. Replacement populations must fit the constructor collection element type;
construct a new state when changing its storage types.
Particle growth requires new velocities at the
existing half level. There is no subcycling or independent species timestep.

Deposited charge and fields must fit T. Physical charge totals and neutrality
checks retain the exact represented weights, charges and box lengths.
Input failures leave physical fields, particles and phase unchanged. Collection
updates are staged, including arithmetic failures; scratch/FFT work may change.
Concurrent calls require independent solver states, populations and grids.
"""
function ElectrostaticPIC(
    g::FourierGrid{D,T},
    species::AbstractVector;
    background_charge = zero(T),
    neutrality::Symbol = :require,
    charge_rtol = 64eps(T),
) where {D,T}
    1 <= D <= 3 || throw(ArgumentError("ElectrostaticPIC supports D = 1, 2, or 3"))
    T in (Float32, Float64) ||
        throw(ArgumentError("electrostatic collections require Float32 or Float64"))
    neutrality in (:require, :neutralize) ||
        throw(ArgumentError("neutrality must be :require or :neutralize"))
    background = _require_finite_real("background_charge", background_charge, T)
    rtol = _require_finite_nonnegative_real("charge_rtol", charge_rtol, T)
    _espic_validate_collection_type(species, g)
    P = eltype(species)
    work = _ElectrostaticCollectionWork{D,T,P}(
        neutrality,
        rtol,
        P[],
        _ElectrostaticHostParticle{D,T}[],
        ntuple(_ -> zeros(T, g.n), 3),
        zeros(T, g.n),
        zeros(T, g.n),
        Ref(false),
    )
    es = ElectrostaticPIC{D,T,typeof(g),typeof(work)}(
        g,
        background,
        ntuple(_ -> zeros(T, g.n), 3),
        zeros(T, g.n),
        ntuple(_ -> T[], 3),
        Ref(false),
        Ref(zero(T)),
        work,
    )
    _espic_validate_species(es, species)
    return es
end

function _espic_collection(es::ElectrostaticPIC)
    work = es.collection
    work === nothing &&
        throw(ArgumentError("species vectors require ElectrostaticPIC(g, species; ...)"))
    return work
end

# Exact stored charge is needed to distinguish genuinely neutral states from
# products that merely round to the same value. The shared dyadic sum uses a
# bounded integer accumulator when the exponent spread permits it.
function _espic_charges(es::ElectrostaticPIC, species)
    particle = zero(Rational{BigInt})
    absolute = zero(Rational{BigInt})
    for ps in species
        term = _exact_rational(ps.q) * _exact_float_sum(ps.weight)
        particle += term
        absolute += abs(term) # Validated particle weights are nonnegative.
    end
    background = _exact_rational(es.n0) * prod(_exact_rational, es.g.L)
    return (; particle, background, net = particle + background, scale = absolute + abs(background))
end

function _espic_check_neutrality(work, charges)
    work.neutrality === :require &&
        abs(charges.net) > _exact_rational(work.charge_rtol) * charges.scale &&
        throw(
            ArgumentError(
                "periodic Poisson requires neutral total charge; supply background_charge or neutrality=:neutralize",
            ),
        )
    return nothing
end

function _espic_validate_collection_type(species, g::FourierGrid{D,T}) where {D,T}
    species isa AbstractVector{<:ParticleSet{D,T}} || throw(
        ArgumentError("species must be a typed vector matching the grid dimension and precision"),
    )
    axes(species, 1) == Base.OneTo(length(species)) ||
        throw(ArgumentError("species vector must use one-based axes"))
    return nothing
end

function _espic_validate_storage(es::ElectrostaticPIC{D,T}, species) where {D,T}
    work = _espic_collection(es)
    _espic_validate_collection_type(species, es.g)
    _require_grid_tuple(:E, es.E, es.g)
    _require_grid_array(:ne, es.ne, es.g)
    _require_grid_tuple(:staged_E, work.E, es.g)
    _require_grid_array(:staged_ne, work.ne, es.g)
    _require_grid_array(:density, work.density, es.g)
    for (i, ps) in enumerate(species)
        ps isa ParticleSet{D,T} ||
            throw(ArgumentError("electrostatic species must match grid dimension and precision"))
        _check_particle_storage_axes(ps)
        all(a -> a isa StridedVector, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)) ||
            throw(ArgumentError("electrostatic collections require CPU strided particle storage"))
        _particle_components_mightalias(ps) &&
            throw(ArgumentError("electrostatic particle components must not alias"))
        for j = 1:(i-1)
            _particle_sets_mightalias(ps, species[j]) &&
                throw(ArgumentError("electrostatic populations must have independent storage"))
        end
        for a in (
            es.E...,
            es.ne,
            es.Ep...,
            work.E...,
            work.ne,
            work.density,
            es.g.kvec...,
            es.g.cbuf,
            es.g.tbuf,
            es.g.abuf,
            es.g.sbuf,
        )
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("electrostatic particles must not alias solver storage"))
        end
        for staged in work.staged
            _particle_sets_mightalias(ps, staged) &&
                throw(ArgumentError("electrostatic particles must not alias staged storage"))
        end
    end
    return nothing
end

function _espic_validate_species(es::ElectrostaticPIC{D,T}, species) where {D,T}
    _espic_validate_storage(es, species)
    for ps in species
        qm = _validated_qm(ps)
        !iszero(ps.q) &&
            iszero(qm) &&
            throw(ArgumentError("electrostatic charge-to-mass ratio underflows; rescale the units"))
        _validate_particle_positions(ps)
        _require_finite_particle_velocities(ps, "electrostatic stepping")
        all(w -> isfinite(w) && w >= zero(T), ps.weight) ||
            throw(ArgumentError("electrostatic weights must be finite and nonnegative"))
    end
    charges = _espic_charges(es, species)
    _espic_check_neutrality(_espic_collection(es), charges)
    return charges
end

function _espic_prepare_stages!(es::ElectrostaticPIC{D,T}, species) where {D,T}
    work = _espic_collection(es)
    while length(work.staged) < length(species)
        push!(work.staged, ParticleSet{D,T}(0))
    end
    resize!(work.staged, length(species))
    for (stage, ps) in zip(work.staged, species)
        for (out, input) in zip(
            (stage.x..., stage.v..., stage.weight, stage.id, stage.tag),
            (ps.x..., ps.v..., ps.weight, ps.id, ps.tag),
        )
            resize!(out, length(input))
            copyto!(out, input)
        end
        stage.q, stage.m = ps.q, ps.m
    end
    return work
end

function _espic_collection_fields!(es, species)
    work = _espic_collection(es)
    _charge_density_sum!(work.ne, work.density, species, es.g, CIC(); sign = -one(eltype(es.ne)))
    all(isfinite, work.ne) ||
        throw(ArgumentError("electrostatic deposited charge exceeds numeric range"))
    _espic_poisson_E!(work.E, work.ne, es.n0, es.g)
    all(a -> all(isfinite, a), work.E) ||
        throw(ArgumentError("electrostatic field exceeds numeric range"))
    return work
end

function _espic_publish_fields!(es, work)
    copyto!(es.ne, work.ne)
    foreach(copyto!, es.E, work.E)
    return es
end

function _espic_validate_owner_types(work, species)
    all(ps -> ps isa eltype(work.owners), species) || throw(
        ArgumentError(
            "electrostatic population storage type changed; construct a state for the new storage",
        ),
    )
    return nothing
end

function init_espic!(es::ElectrostaticPIC{D,T}, species::AbstractVector) where {D,T}
    _espic_validate_owner_types(_espic_collection(es), species)
    _espic_validate_species(es, species)
    work = _espic_collection_fields!(es, species)
    _espic_publish_fields!(es, work)
    empty!(work.owners)
    append!(work.owners, species)
    work.initialized[] = true
    es.primed[] = false
    es.last_dt[] = zero(T)
    return es
end

function step_espic!(es::ElectrostaticPIC{D,T}, species::AbstractVector, dt::Real) where {D,T}
    dtT =
        _require_resolved_dt("step_espic!", _validated_nonnegative_dt(T, dt; name = "step_espic!"))
    _espic_validate_species(es, species)
    work = _espic_collection(es)
    work.initialized[] ||
        throw(ArgumentError("initialize the electrostatic species collection before stepping"))
    _espic_require_bound_species(work, species)
    isfinite(es.last_dt[]) && (es.primed[] ? es.last_dt[] > zero(T) : iszero(es.last_dt[])) ||
        throw(ArgumentError("electrostatic last_dt must be positive exactly when primed"))
    all(a -> all(isfinite, a), es.E) || throw(ArgumentError("electrostatic field must be finite"))
    iszero(dtT) && return es
    _espic_prepare_stages!(es, species)
    for ps in work.staged
        _resize_espic_particle_workspaces!(es, nparticles(ps))
        gather_vector!(es.Ep, es.E, ps, es.g, CIC())
        _espic_push!(ps, es.Ep, es.g, dtT, es.primed[], es.last_dt[])
    end
    _espic_collection_fields!(es, work.staged)
    for (ps, stage) in zip(species, work.staged)
        foreach(copyto!, ps.x, stage.x)
        foreach(copyto!, ps.v, stage.v)
    end
    _espic_publish_fields!(es, work)
    es.primed[] = true
    es.last_dt[] = dtT
    return es
end

"""
    electrostatic_diagnostics(es::ElectrostaticPIC, species)

Particle/background charge, explicitly neutralizing charge, deposited charge,
particle momentum, kinetic energy, shared field energy (once), and total energy.
Particle/background charge totals are exact Rational{BigInt} sums of the stored
inputs; deposited charge also retains the complete represented grid integral.
`velocity_time_offset` is `-last_dt/2` after stepping: energies use stored
half-time velocities and full-time E, without silently centering velocities.
The integrator converges in energy; this is not an exact discrete energy law.
"""
function electrostatic_diagnostics(es::ElectrostaticPIC{D,T}, species::AbstractVector) where {D,T}
    charges = _espic_validate_species(es, species)
    work = _espic_collection(es)
    work.initialized[] ||
        throw(ArgumentError("initialize the electrostatic species collection before diagnostics"))
    _espic_require_bound_species(work, species)
    kinetic = kinetic_energy(species)
    momentum = isempty(species) ? ntuple(_ -> zero(T), 3) : total_momentum(species)
    field = field_energy(es)
    neutralizing = work.neutrality === :neutralize ? -charges.net : zero(charges.net)
    return (
        particle_charge = charges.particle,
        background_charge = charges.background,
        neutralizing_charge = neutralizing,
        net_charge = charges.net + neutralizing,
        deposited_charge = (_exact_rational(es.n0) * length(es.ne) - _exact_float_sum(es.ne)) *
                           prod(_exact_rational, es.g.dx),
        momentum = momentum,
        kinetic = kinetic,
        field = field,
        total = kinetic + field,
        velocity_time_offset = es.primed[] ? -es.last_dt[] / 2 : zero(T),
    )
end

function _espic_require_bound_species(work, species)
    length(work.owners) == length(species) &&
        all(a === b for (a, b) in zip(work.owners, species)) ||
        throw(ArgumentError("electrostatic population membership/order changed; call init_espic!"))
    return nothing
end
