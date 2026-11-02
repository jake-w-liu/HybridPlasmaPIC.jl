# Two-face full-time transport/source orchestration; the flux and Boris kernels
# remain canonical. Every object returned here is still a private candidate.
function _open_hall_birth_plan(st, h)
    ns = length(st.owners)
    counts = zeros(Int, 2ns)
    alpha = copy(st.alpha)
    next = copy(st.next_index)
    total = BigInt(0)
    area = prod(_exact_rational, Base.tail(st.op.L))
    for s = 1:ns
        produced = BigInt(0)
        for face = 1:2
            i = 2(s - 1) + face
            record = st.sources[s][face]
            u = alpha[i] + _exact_rational(record.density) * _exact_rational(st.fluxes[s][face]) *
                area * _exact_rational(h) / _exact_rational(record.weight)
            number = fld(numerator(u), denominator(u))
            number >= 0 && number <= st.controls.max_births - total ||
                throw(ArgumentError("open Hall predicted births exceed max_births; no source count is capped"))
            counts[i] = Int(number)
            alpha[i] = u - number
            total += number
            produced += number
        end
        BigInt(st.next_index[s]) + produced <= BigInt(_PARTICLE_ID_MAX_INDEX) + 1 ||
            throw(ArgumentError("open Hall particle provenance index is exhausted"))
        next[s] = UInt64(BigInt(st.next_index[s]) + produced)
    end
    old = sum(BigInt(nparticles(ps)) for ps in st.owners)
    old + total <= st.controls.max_particles ||
        throw(ArgumentError("open Hall old plus all births exceeds max_particles before transport"))
    return (counts = counts, alpha = alpha, next_index = next, total = Int(total))
end

function _open_hall_event(ps, weight, velocity)
    Q = _exact_rational
    mass, charge = Q(ps.m) * Q(weight), Q(ps.q) * Q(weight)
    v = map(Q, velocity)
    return _OpenHallExchange(BigInt(1), mass, charge, ntuple(c -> mass * v[c], 3), mass * sum(abs2, v) / 2)
end

# Return the first normal exit or a once-rounded retained endpoint. Periodic
# wrapping composes the complete displacement before one exact mod and round.
function _open_hall_ballistic(position, velocity, age, op::MixedShockOperators{D,T}) where {D,T}
    Q = _exact_rational
    x, vx = Q(position[1]), Q(velocity[1])
    if vx != 0
        face = vx > 0 ? 2 : 1
        boundary = face == 1 ? zero(_OpenHallQ) : Q(op.L[1])
        exit_time = (boundary - x) / vx
        0 <= exit_time <= age && return (face = face, position = nothing)
    end
    endpoint = ntuple(Val(D)) do d
        value = _particle_coordinate_native(Q(position[d]), Q(velocity[d]), age)
        d > 1 && (value = mod(value, Q(op.L[d])))
        stored = _from_exact(T, value)
        isfinite(stored) || throw(ArgumentError("open Hall retained ballistic coordinate exceeds its finite range"))
        d > 1 && stored == op.L[d] && (stored = zero(T))
        (d == 1 ? 0 <= stored <= op.L[d] : 0 <= stored < op.L[d]) ||
            throw(ArgumentError("open Hall retained coordinate crosses its physical domain"))
        return stored
    end
    return (face = 0, position = endpoint)
end

function _open_hall_put_particle!(ps, p, position, velocity, weight, id, tag)
    for d in eachindex(ps.x)
        ps.x[d][p] = position[d]
    end
    for c = 1:3
        ps.v[c][p] = velocity[c]
    end
    ps.weight[p], ps.id[p], ps.tag[p] = weight, id, tag
    return nothing
end

function _open_hall_transport!(st::OpenHallHybrid{D,T}, h, plan) where {D,T}
    # Copy only after every source count/capacity/ID range has been admitted.
    copy!(st.rng_candidate, st.rng)
    _validate_collisional_rng_pair(st.rng, st.rng_candidate, "open Hall source")
    draws = _ReservoirDrawView(st.rng_candidate, st.controls.max_draws, 0)
    ns = length(st.owners)
    prepared = _ElectrostaticHostParticle{D,T}[]
    retained_counts = Int[]
    born = [_OpenHallExchange() for _ = 1:2ns]
    exited = [_OpenHallExchange() for _ = 1:2ns]
    kick_energy = kick_work = zero(_OpenHallQ)
    kick_momentum = ntuple(_ -> zero(_OpenHallQ), 3)
    Q = _exact_rational
    for s = 1:ns
        ps = st.owners[s]
        n = nparticles(ps)
        capacity = Base.checked_add(n, Base.checked_add(plan.counts[2s - 1], plan.counts[2s]))
        candidate = ParticleSet{D,T}(capacity; q = ps.q, m = ps.m)
        _open_hall_gather!(st, ps, st.work.force)
        factor = _validated_boris_factor(ps)
        hk = _boris_native_factor(factor) * h / T(2)
        scales = _boris_batch_scales(ps.v, st.work.Ep, st.work.Bp, hk, 1:n)
        native, _ = _boris_batch_admission(ps, hk, factor, h, scales...)
        guard = Val(!native)
        # Kick every old record first. No preceding kick is assigned to newborns.
        for p = 1:n
            old = ntuple(c -> ps.v[c][p], 3)
            E = ntuple(c -> st.work.Ep[c][p], 3)
            B = ntuple(c -> st.work.Bp[c][p], 3)
            velocity = _kick_eval(guard, old..., E..., B..., factor, h)
            all(isfinite, velocity) || throw(ArgumentError("open Hall Boris velocity exceeds its finite materialized range"))
            _open_hall_put_particle!(candidate, p, ntuple(d -> ps.x[d][p], D), velocity,
                ps.weight[p], ps.id[p], ps.tag[p])
            mass, charge = Q(ps.m) * Q(ps.weight[p]), Q(ps.q) * Q(ps.weight[p])
            kick_energy += mass * sum(Q(velocity[c])^2 - Q(old[c])^2 for c = 1:3) / 2
            kick_work += Q(h) * charge * sum((Q(old[c]) + Q(velocity[c])) * Q(E[c]) for c = 1:3) / 2
            kick_momentum = ntuple(c -> kick_momentum[c] + mass * (Q(velocity[c]) - Q(old[c])), 3)
        end
        retained = 0
        # Stable forward compaction is safe: every written index is <= p.
        for p = 1:n
            velocity = ntuple(c -> candidate.v[c][p], 3)
            result = _open_hall_ballistic(ntuple(d -> candidate.x[d][p], D), velocity, Q(h), st.op)
            if result.face == 0
                retained += 1
                _open_hall_put_particle!(candidate, retained, result.position, velocity,
                    candidate.weight[p], candidate.id[p], candidate.tag[p])
            else
                i = 2(s - 1) + result.face
                exited[i] = _open_hall_exchange_add(exited[i], _open_hall_event(ps, candidate.weight[p], velocity))
            end
        end
        next = st.next_index[s]
        for face = 1:2
            i = 2(s - 1) + face
            source = st.sources[s][face]
            for birth = 1:plan.counts[i]
                speed = flux_speed(draws, source.drift, source.spread)
                age_draw = rand(draws, T)
                isfinite(speed) && speed >= 0 && isfinite(age_draw) && 0 <= age_draw < 1 ||
                    throw(ArgumentError("open Hall reservoir produced an invalid normal speed or age draw"))
                position = ntuple(Val(D)) do d
                    d == 1 && return face == 1 ? zero(T) : st.op.L[1]
                    uniform = rand(draws, T)
                    isfinite(uniform) && 0 <= uniform < 1 || throw(ArgumentError("open Hall reservoir position draw is outside[0,1)"))
                    value = _from_exact(T, Q(st.op.L[d]) * Q(uniform))
                    return value == st.op.L[d] ? zero(T) : value
                end
                normals = (randn(draws, T), randn(draws, T))
                velocity = (face == 1 ? speed : -speed,
                    source.tangent_drift[1] + source.tangent_spread[1] * normals[1],
                    source.tangent_drift[2] + source.tangent_spread[2] * normals[2])
                all(isfinite, velocity) || throw(ArgumentError("open Hall sampled tangent velocity exceeds its finite range"))
                id = global_particle_id(next; species = st.namespaces[s])
                next += UInt64(1) # The complete two-face reservation already fits.
                event = _open_hall_event(ps, source.weight, velocity)
                born[i] = _open_hall_exchange_add(born[i], event)
                result = _open_hall_ballistic(position, velocity, Q(h) * Q(age_draw), st.op)
                if result.face == 0
                    retained += 1
                    _open_hall_put_particle!(candidate, retained, result.position, velocity,
                        source.weight, id, source.tag)
                else
                    out = 2(s - 1) + result.face
                    exited[out] = _open_hall_exchange_add(exited[out], event)
                end
            end
        end
        next == plan.next_index[s] || error("open Hall source reservation and draw traversal disagree")
        for a in (candidate.x..., candidate.v..., candidate.weight, candidate.id, candidate.tag)
            resize!(a, retained)
        end
        push!(prepared, candidate)
        push!(retained_counts, retained)
    end
    _open_hall_particles(st, prepared)
    before, after = _open_hall_particle_totals(st.owners), _open_hall_particle_totals(prepared)
    births = foldl(_open_hall_exchange_add, born; init = _OpenHallExchange())
    exits = foldl(_open_hall_exchange_add, exited; init = _OpenHallExchange())
    expected = _open_hall_exchange_add(_open_hall_exchange_add(before, births), exits, -1)
    after.count == expected.count && after.mass == expected.mass && after.charge == expected.charge &&
        after.momentum == ntuple(c -> expected.momentum[c] + kick_momentum[c], 3) &&
        after.kinetic == expected.kinetic + kick_energy || error("open Hall retained-record event/kick identities disagree")
    return (particles = prepared, counts = retained_counts, born = born, exited = exited,
        kick_energy = kick_energy, kick_work = kick_work, kick_momentum = kick_momentum,
        event_energy = births.kinetic - exits.kinetic, draws = draws.used,
        births = Int(births.count), exits = Int(exits.count), old_totals = before, new_totals = after)
end
