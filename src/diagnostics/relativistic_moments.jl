# Complete relativistic particle moments reuse the existing residual-ball and
# exact rational-radical owners. No separate Lorentz-factor approximation lives
# here; ordinary Float64 error bounds include every original source operand.

function _relativistic_particle_moments_exact(species, c::T) where {T}
    energy = _relativistic_moment_exact(species, c, Val(0))
    momentum = ntuple(d -> _relativistic_moment_exact(species, c, Val(d)), Val(3))
    return energy, momentum
end

function _relativistic_particle_moments(species, c::T) where {T}
    for ps in species
        _check_particle_storage_axes(ps)
        all(a -> a isa StridedVector, (ps.v..., ps.weight)) ||
            throw(ArgumentError("relativistic moments require CPU strided velocity/weight storage"))
        isfinite(ps.m) && all(isfinite, ps.weight) ||
            throw(ArgumentError("relativistic moments require finite masses and weights"))
        _require_finite_particle_velocities(ps, "relativistic moments")
    end
    _boris_accuracy_mode() || return _relativistic_particle_moments_exact(species, c)
    _rel_prime_dd_operand(c) || return _relativistic_particle_moments_exact(species, c)
    sums = ntuple(_ -> _rel_prime_ball(0.0), Val(4))
    active = ntuple(_ -> false, Val(4))
    for ps in species
        _rel_prime_dd_operand(ps.m) || return _relativistic_particle_moments_exact(species, c)
        @inbounds for p in eachindex(ps.weight)
            velocity = (ps.v[1][p], ps.v[2][p], ps.v[3][p])
            weight = ps.weight[p]
            all(_rel_prime_dd_operand, (velocity..., weight)) ||
                return _relativistic_particle_moments_exact(species, c)
            # Zero velocity is an exactly valid physical state at positive c.
            all(iszero, velocity) && continue
            gamma = _rel_prime_gamma_ball(velocity, c)
            _rel_prime_ball_valid(gamma) || return _relativistic_particle_moments_exact(species, c)
            (iszero(ps.m) || iszero(weight)) && continue
            weighted_mass = _rel_prime_ball_mul(_rel_prime_ball(ps.m), _rel_prime_ball(weight))
            speed2 = _rel_prime_ball(0.0)
            for x in velocity
                iszero(x) && continue
                term = _rel_prime_ball_mul(_rel_prime_ball(x), _rel_prime_ball(x))
                speed2 = _rel_prime_ball_add(speed2, term)
            end
            # (gamma-1)c² = |v|²*gamma*(gamma/(gamma+1)). This order
            # avoids a small-speed subtraction and never squares a large gamma.
            ratio = _rel_prime_ball_div(gamma, _rel_prime_ball_add(gamma, _rel_prime_ball(1.0)))
            energy = _rel_prime_ball_mul(
                _rel_prime_ball_mul(weighted_mass, speed2),
                _rel_prime_ball_mul(gamma, ratio),
            )
            sums = Base.setindex(sums, _rel_prime_ball_add(sums[1], energy), 1)
            active = Base.setindex(active, true, 1)
            for d = 1:3
                iszero(velocity[d]) && continue
                term = _rel_prime_ball_mul(
                    weighted_mass,
                    _rel_prime_ball_mul(gamma, _rel_prime_ball(velocity[d])),
                )
                sums = Base.setindex(sums, _rel_prime_ball_add(sums[d+1], term), d + 1)
                active = Base.setindex(active, true, d + 1)
            end
        end
    end
    values = (
        _rel_moment_entry(T, species, c, active[1], sums[1], 1),
        _rel_moment_entry(T, species, c, active[2], sums[2], 2),
        _rel_moment_entry(T, species, c, active[3], sums[3], 3),
        _rel_moment_entry(T, species, c, active[4], sums[4], 4),
    )
    return values[1], (values[2], values[3], values[4])
end

function _rel_moment_entry(::Type{T}, species, c, isactive, sum, i) where {T}
    isactive || return zero(T)
    value, resolved = _rel_prime_round_cell(T, sum[1], sum[2])
    return resolved ? value : _relativistic_moment_exact(species, c, Val(i - 1))
end

"""
    total_momentum_relativistic(particles, c) -> NTuple{3}
    total_momentum_relativistic(species::AbstractVector, c) -> NTuple{3}

Physical particle momentum `Σ w*m*γ*v`, with `γ = 1/sqrt(1-|v|²/c²)`.
Supports CPU Float32/Float64 populations with finite masses, weights and stored
subluminal velocities. Signed finite mass/weight factors retain their algebraic
meaning. A typed empty collection returns zero. This diagnostic's physical
`|v|<c` domain is wider than the mover's supported rounded-output domain.

The existing bounded residual-ball owner encloses original sources and whole
collection sums; unresolved operands or rounding cells use the same exact
radical owner as kinetic energy. Overflowing population momenta can cancel
before final rounding. Velocities use their stored time level; this function
does not synchronize particles or fields.
"""
function total_momentum_relativistic(
    ps::ParticleSet{D,T},
    c::Real,
) where {D,T<:Union{Float32,Float64}}
    cT = _require_finite_positive_real("c", c, T)
    return last(_relativistic_particle_moments((ps,), cT))
end

function total_momentum_relativistic(
    species::AbstractVector{<:ParticleSet{D,T}},
    c::Real,
) where {D,T<:Union{Float32,Float64}}
    cT = _require_finite_positive_real("c", c, T)
    return last(_relativistic_particle_moments(species, cT))
end
