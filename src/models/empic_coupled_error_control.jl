# Full-time coupled projected Euler/Heun control. Legacy half phases are separate.
const _EMPICCoupledQ = Rational{BigInt}
const _EC_CONTROL_KEYS = (
    :pair,
    :metric,
    :domain,
    :atol,
    :rtol,
    :h_min,
    :h_max,
    :max_trials,
    :max_force_evals,
    :max_deposits,
    :max_fft_calls,
    :max_report_cells,
)
const _EC_COST_KEYS = (:trials, :force_evals, :deposits, :fft_calls, :report_cells)
const _EC_PHASE_KEYS = (
    :time,
    :elapsed,
    :interval,
    :last,
    :next,
    :next_h,
    :initial_energy,
    :field_defect,
    :particle_defect,
    :coupling_defect,
    :admission_control,
    :costs,
    :report,
)
struct _EMPICCoupledPhase{T}
    time::_EMPICCoupledQ
    elapsed::_EMPICCoupledQ
    interval::Vector{T}
    last::Vector{_EMPICCoupledQ}
    next::Vector{_EMPICCoupledQ}
    next_h::T
    initial_energy::_EMPICCoupledQ
    field_defect::_EMPICCoupledQ
    particle_defect::_EMPICCoupledQ
    coupling_defect::_EMPICCoupledQ
    admission_control::NamedTuple
    costs::NamedTuple
    report::Union{Nothing,NamedTuple}
end
mutable struct _EMPICCoupledSchedule{T}
    proposal::NamedTuple
    control::NamedTuple
    phase::_EMPICCoupledPhase{T}
    candidate::_EMPICCoupledPhase{T}
    shadow::Any
    candidate_shadow::Any
    signature::Union{Nothing,String}
    busy::Threads.Atomic{Bool}
end
mutable struct _EMPICCoupledMeter
    control::NamedTuple
    trials::Int
    force_evals::Int
    deposits::Int
    fft_calls::Int
    report_cells::Int
    continuity::Union{Nothing,NamedTuple}
end
struct _EMPICCoupledReject <: Exception
    reason::Symbol
    stage::Symbol
    species::Int
    particle::Int
    defect::Union{Nothing,_EMPICCoupledQ}
end
Base.showerror(io::IO, e::_EMPICCoupledReject) = print(
    io,
    "coupled trial rejected: ",
    e.reason,
    " at ",
    e.stage,
    " species=",
    e.species,
    " particle=",
    e.particle,
    " defect=",
    e.defect,
)
_ec_q(x) = x isa _EMPICCoupledQ && denominator(x) > 0
function _ec_identical(left, right)
    typeof(left) === typeof(right) || return false
    if left isa NamedTuple || left isa Tuple
        return length(left) == length(right) &&
               all(pair -> _ec_identical(pair[1], pair[2]), zip(left, right))
    elseif left isa AbstractArray
        return axes(left) == axes(right) &&
               all(pair -> _ec_identical(pair[1], pair[2]), zip(left, right))
    end
    return isequal(left, right)
end
_ec_zero_costs() = (trials = 0, force_evals = 0, deposits = 0, fft_calls = 0, report_cells = 0)
_ec_meter(control) = _EMPICCoupledMeter(control, 0, 0, 0, 0, 0, nothing)
_ec_costs(m::_EMPICCoupledMeter) =
    NamedTuple{_EC_COST_KEYS}(Tuple(getfield(m, k) for k in _EC_COST_KEYS))
function _ec_costs_check(costs)
    costs isa NamedTuple &&
        keys(costs) == _EC_COST_KEYS &&
        all(x -> x isa Int && x >= 0, values(costs)) ||
        throw(ArgumentError("coupled cost layout differs"))
    return nothing
end
_ec_cost_add(a, b) = NamedTuple{_EC_COST_KEYS}(
    Tuple(Base.checked_add(getproperty(a, k), getproperty(b, k)) for k in _EC_COST_KEYS),
)
_ec_cost_sub(a, b) =
    NamedTuple{_EC_COST_KEYS}(Tuple(getproperty(a, k) - getproperty(b, k) for k in _EC_COST_KEYS))
function _ec_charge!(meter::_EMPICCoupledMeter, role::Symbol, count::Int)
    count >= 0 || throw(ArgumentError("coupled work charge must be nonnegative"))
    role in _EC_COST_KEYS || throw(ArgumentError("unknown coupled work role"))
    limit = getproperty(meter.control, Symbol("max_", role))::Int
    next = Base.checked_add(getfield(meter, role)::Int, count)
    next <= limit ||
        throw(ArgumentError("coupled $(role) work budget exhausted: required=$next, cap=$limit"))
    setfield!(meter, role, next)
    return nothing
end
_empic_fft_charge!(meter::_EMPICCoupledMeter, count::Int) = _ec_charge!(meter, :fft_calls, count)
function _empic_continuity_observe!(meter::_EMPICCoupledMeter, absolute, scale)
    meter.continuity =
        isfinite(absolute) && isfinite(scale) ?
        (absolute = _exact_rational(absolute), scale = _exact_rational(scale)) : nothing
    return nothing
end
_ec_reject(reason, stage, s = 0, p = 0, defect = nothing) =
    throw(_EMPICCoupledReject(reason, stage, s, p, defect))
function _ec_mode(::Type{T}) where {T}
    rounding(T) === RoundNearest && rounding(Float64) === RoundNearest && !get_zero_subnormals() ||
        throw(ArgumentError("coupled control requires nearest rounding and gradual underflow"))
    return nothing
end
function _ec_control(proposal, ::Type{T}) where {T}
    proposal isa NamedTuple && keys(proposal) == _EC_CONTROL_KEYS ||
        throw(ArgumentError("coupled error_control requires the closed projected-pair proposal"))
    proposal.pair === :projected_euler_heun_v1 &&
        proposal.metric === :component_max_v1 &&
        proposal.domain === :cic_tsc_continuous_v1 ||
        throw(ArgumentError("unsupported coupled pair/metric/domain"))
    roles = (:x, :v, :E, :B, :clock)
    for name in (:atol, :rtol)
        record = getproperty(proposal, name)
        record isa NamedTuple && keys(record) == roles ||
            throw(ArgumentError("coupled tolerances require x/v/E/B/clock roles"))
    end
    atols = NamedTuple{roles}(
        map(values(proposal.atol)) do v
            result = _empic_event_scalar("coupled absolute tolerance", v, T)
            _stored_positive(result) ||
                throw(ArgumentError("coupled absolute tolerances must be positive"))
            result
        end,
    )
    rtols = NamedTuple{roles}(
        map(values(proposal.rtol)) do v
            result = _empic_event_scalar("coupled relative tolerance", v, T)
            _stored_nonnegative(result) && result < one(T) ||
                throw(ArgumentError("coupled relative tolerances must be in[0,1)"))
            result
        end,
    )
    minimum = _empic_event_scalar("coupled h_min", proposal.h_min, T)
    maximum = _empic_event_scalar("coupled h_max", proposal.h_max, T)
    0 < minimum <= maximum ||
        throw(ArgumentError("coupled interval bounds must be positive and ordered"))
    _require_resolved_dt("coupled h_min", minimum)
    _require_resolved_dt("coupled h_max", maximum)
    names = (:max_trials, :max_force_evals, :max_deposits, :max_fft_calls, :max_report_cells)
    limits = NamedTuple{names}(
        map(names) do name
            v = getproperty(proposal, name)
            v isa Integer && !(v isa Bool) && 1 <= v <= typemax(Int) ||
                throw(ArgumentError("coupled work caps must be positive fitting Int values"))
            Int(v)
        end,
    )
    return merge(
        (
            pair = proposal.pair,
            metric = proposal.metric,
            domain = proposal.domain,
            atol = atols,
            rtol = rtols,
            h_min = minimum,
            h_max = maximum,
        ),
        limits,
    )
end
function _ec_zero_phase(
    ::Type{T},
    count,
    control;
    energy = zero(_EMPICCoupledQ),
    costs = _ec_zero_costs(),
) where {T}
    z = zero(_EMPICCoupledQ)
    return _EMPICCoupledPhase{T}(
        z,
        z,
        zeros(T, count),
        fill(z, count),
        fill(z, count),
        control.h_max,
        energy,
        z,
        z,
        z,
        control,
        costs,
        nothing,
    )
end
function _empic_coupled_schedule(
    mobile,
    ns,
    stride,
    cap,
    adaptive,
    minimum,
    limits,
    budget,
    ::Type{T},
    control,
) where {T}
    cap === nothing && throw(
        ArgumentError("coupled error control requires species_dt and event_budget at construction"),
    )
    stride === nothing || throw(ArgumentError("coupled control cannot use species_stride"))
    proposal = _empic_event_proposal(mobile, cap, adaptive, minimum, limits, budget, T)
    ns <= proposal.budget ||
        throw(ArgumentError("base field partition exceeds coupled event budget"))
    closed = _ec_control(control, T)
    phase = _ec_zero_phase(T, length(mobile), closed)
    return _EMPICCoupledSchedule{T}(
        proposal,
        closed,
        phase,
        _ec_zero_phase(T, length(mobile), closed),
        nothing,
        nothing,
        nothing,
        Threads.Atomic{Bool}(false),
    )
end
_empic_coupled_bind!(es, schedule) = nothing
function _ec_spectral_signature(g)
    expected = SpectralOperators._fourier_metadata(g.n, g.L)
    for role in (:n, :L, :dx, :midx, :kfull, :ik, :kvec)
        isequal(getproperty(g, role), getproperty(expected, role)) ||
            throw(ArgumentError("coupled Fourier metadata is not canonical"))
    end
    bytes = UInt8[]
    for vector in
        (collect(g.n), collect(g.L), collect(g.dx), g.midx..., g.kfull..., g.ik..., g.kvec...)
        append!(bytes, reinterpret(UInt8, vector))
    end
    return bytes2hex(SHA.sha256(bytes))
end
function _empic_coupled_bind!(es, schedule::_EMPICCoupledSchedule)
    schedule.signature = _ec_spectral_signature(es.g)
    return nothing
end
function _ec_acquire(schedule)
    Threads.atomic_cas!(schedule.busy, false, true) == false ||
        throw(ArgumentError("coupled EMPIC already has an active operation"))
    return nothing
end
_ec_release(schedule) = (Threads.atomic_xchg!(schedule.busy, false); nothing)
function _ec_copy_fields(fields)
    return NamedTuple{keys(fields)}(map(v -> v isa Tuple ? map(copy, v) : copy(v), values(fields)))
end
function _ec_copy_species(species, ::Val{D}, ::Type{T}) where {D,T}
    result = _ElectrostaticHostParticle{D,T}[]
    _pic_prepare_species_stages!(result, species)
    return result
end
function _ec_new_stage(
    es::EMPIC{D,T};
    c = es.c,
    shape = es.shape,
    relativistic = es.relativistic,
    n_sub = es.n_sub,
) where {D,T}
    return _new_empic_state(es.g, 0, es.n0, c, shape, relativistic, false, one(T), n_sub, nothing)
end
function _ec_new_stage(
    es::EMPIC1D{T};
    c = es.c,
    shape = es.shape,
    relativistic = es.relativistic,
    n_sub = es.n_sub,
) where {T}
    return _new_empic1d_state(es.g, 0, es.n0, c, shape, relativistic, false, one(T), n_sub, nothing)
end
function _ec_restore_fields!(stage::Union{EMPIC,EMPIC1D}, fields::NamedTuple)
    _restore_checkpoint_fields!(_empic_event_fields(stage), fields)
    return stage
end
function _ec_snapshot(stage, species)
    return (
        species = NamedTuple[deepcopy(_checkpoint_particle_state(ps)) for ps in species],
        fields = _ec_copy_fields(_empic_event_fields(stage)),
    )
end
_ec_snapshot_populations(record, ::Val{D}, ::Type{T}) where {D,T} =
    _ElectrostaticHostParticle{D,T}[_checkpoint_species(ps, Val(D), T) for ps in record.species]
function _ec_energy(stage, species)
    U = _empic_event_energy(stage)
    K = if stage.relativistic
        value = first(_relativistic_particle_moments(species, stage.c))
        isfinite(value) ||
            throw(ArgumentError("coupled reported relativistic kinetic energy is nonfinite"))
        _exact_rational(value)
    else
        _kinetic_energy_rational(species)
    end
    return K, U
end
function _ec_count_particles(species, mobile = nothing)
    count = 0
    for (i, ps) in enumerate(species)
        mobile === nothing || mobile[i] || continue
        count = Base.checked_add(count, nparticles(ps))
    end
    return count
end
function _ec_source!(
    rho,
    stage::Union{EMPIC,EMPIC1D},
    species::AbstractVector,
    target,
    meter::_EMPICCoupledMeter,
)
    _ec_charge!(meter, :deposits, _ec_count_particles(species))
    return _empic_charge_sum!(rho, stage.ne, species, stage.g, stage.shape, stage.n0, target)
end
function _ec_current!(stage::EMPIC, species, mobile, meter)
    _ec_charge!(meter, :deposits, Base.checked_mul(3, _ec_count_particles(species, mobile)))
    midpoints = [map(_ -> eltype(stage.rho_n)[], ps.x) for ps in species]
    _empic_current_sum!(
        stage.Jraw,
        (work = stage.worke, density = stage.ne),
        species,
        midpoints,
        mobile,
        stage.g,
        stage.shape,
        Val(:stored_positions),
    )
    return (current = map(copy, stage.Jraw), dc = nothing)
end
function _ec_current!(stage::EMPIC1D{T}, species, mobile, meter) where {T}
    _ec_charge!(meter, :deposits, Base.checked_mul(2, _ec_count_particles(species, mobile)))
    Jy = similar(stage.Jy)
    _empic_exact_current_sum!(
        Jy,
        species,
        fill(nothing, length(species)),
        mobile,
        stage.g,
        stage.shape,
        2,
        Val(:stored_positions),
    )
    all(isfinite, Jy) || throw(ArgumentError("coupled raw transverse current is nonfinite"))
    dc = _empic1d_collection_Jx0(stage, (staged = species, mobile = mobile))
    return (current = (Jy,), dc = dc)
end
function _ec_average_current(a, b, ::Type{T}) where {T}
    current = map(a.current, b.current) do left, right
        map((x, y) -> _from_exact(T, (_exact_rational(x) + _exact_rational(y)) / 2), left, right)
    end
    if a.dc === nothing
        b.dc === nothing ||
            throw(ArgumentError("coupled branch current layouts differ"))
        all(x -> all(isfinite, x), current) ||
            throw(ArgumentError("coupled averaged current is nonfinite"))
        return (current = current, dc = nothing)
    end
    b.dc === nothing && throw(ArgumentError("coupled branch current layouts differ"))
    dc = _from_exact(T, (_exact_rational(a.dc) + _exact_rational(b.dc)) / 2)
    all(x -> all(isfinite, x), current) && isfinite(dc) ||
        throw(ArgumentError("coupled averaged current is nonfinite"))
    return (current = current, dc = dc)
end
function _ec_correct!(stage::EMPIC, raw, h, meter, which)
    foreach(copyto!, stage.Jraw, raw.current)
    _empic_correct_current!(stage.J, stage.Jraw, stage.rho_n, stage.rho_np1, stage.g, h, meter)
    all(a -> all(isfinite, a), stage.J) || _ec_reject(:nonfinite_current, which)
    return nothing
end
function _ec_correct!(stage::EMPIC1D, raw, h, meter, which)
    copyto!(stage.Jy, only(raw.current))
    _zero_nyquist_mode!(stage.Jy)
    _esirkepov_Jx!(stage, raw.dc, h, meter)
    all(a -> all(isfinite, a), (stage.Jx, stage.Jy)) || _ec_reject(:nonfinite_current, which)
    return nothing
end
function _ec_force!(stage, species, mobile, H, meter)
    rows = Union{Nothing,NamedTuple}[]
    for (i, ps) in enumerate(species)
        if !mobile[i]
            push!(rows, nothing)
            continue
        end
        _ec_charge!(meter, :force_evals, nparticles(ps))
        _empic_collection_gather!(stage, ps; initial = true)
        T = eltype(ps)
        E, B = if stage isa EMPIC
            map(copy, stage.Ep), map(copy, stage.Bp)
        else
            (
                (copy(stage.Exp), copy(stage.Eyp), zeros(T, nparticles(ps))),
                (zeros(T, nparticles(ps)), zeros(T, nparticles(ps)), copy(stage.Bzp)),
            )
        end
        all(a -> all(isfinite, a), (E..., B...)) || _ec_reject(:nonfinite_gather, :force, i)
        factor = _validated_boris_factor(ps)
        increments = NamedTuple[
            _em_lorentz_scaled_parts(
                ntuple(c -> ps.v[c][p], 3),
                ntuple(c -> E[c][p], 3),
                ntuple(c -> B[c][p], 3),
                factor,
                H,
                stage.c,
                Val(stage.relativistic),
            ) for p in eachindex(ps.weight)
        ]
        push!(rows, (E = E, B = B, increments = increments))
    end
    return rows
end
function _ec_parts(derivative, parts)
    if parts === nothing
        all(isfinite, derivative) || _ec_reject(:nonfinite_derivative, :field)
        return [(copy(derivative), one(_EMPICCoupledQ))]
    end
    all(pair -> all(isfinite, first(pair)), parts) || _ec_reject(:nonfinite_derivative, :field)
    return [(copy(array), scale) for (array, scale) in parts]
end
function _ec_derivatives!(stage::EMPIC, meter)
    ep = _em_curl_parts!(stage.curlE, stage.E, stage.g, meter)
    bp = _em_curl_parts!(stage.curlB, stage.B, stage.g, meter)
    return (
        electric = ntuple(c -> _ec_parts(stage.curlE[c], ep === nothing ? nothing : ep[c]), 3),
        magnetic = ntuple(c -> _ec_parts(stage.curlB[c], bp === nothing ? nothing : bp[c]), 3),
    )
end
function _ec_derivatives!(stage::EMPIC1D{T}, meter) where {T}
    ep = _em_derivative_parts!(stage.dEy, stage.Ey, stage.g, 1, meter)
    bp = _em_derivative_parts!(stage.dBz, stage.Bz, stage.g, 1, meter)
    electric = (_ec_parts(stage.dEy, ep),)
    transverse = [(array, -scale) for (array, scale) in _ec_parts(stage.dBz, bp)]
    longitudinal = Tuple{Array{T,1},_EMPICCoupledQ}[]
    return (electric = electric, magnetic = (longitudinal, transverse))
end
function _ec_cut_receipt(x, y, L, dx, shape)
    Q = _exact_rational
    period = Q(L)
    spacing = Q(dx)
    iszero(y - x) && return (direction = 0, seam = false, segments = ())
    direction = y > x ? 1 : -1
    parameters = _EMPICCoupledQ[0, 1]
    first_plane = ceil(BigInt, min(x, y) / period)
    last_plane = floor(BigInt, max(x, y) / period)
    last_plane - first_plane <= 1 ||
        throw(ArgumentError("coupled branch crossed its admitted travel range"))
    for plane = first_plane:last_plane
        boundary = plane * period
        theta = (boundary - x) / (y - x)
        0 < theta < 1 && push!(parameters, theta)
    end
    sort!(unique!(parameters))
    offset = shape isa CIC ? zero(_EMPICCoupledQ) : BigInt(1) // BigInt(2)
    pieces = NamedTuple[]
    for j = 1:length(parameters)-1
        a, b = parameters[j], parameters[j+1]
        tile = floor(BigInt, (x + (a + b) * (y - x) / 2) / period)
        left = x + a * (y - x) - tile * period
        right = x + b * (y - x) - tile * period
        first = ceil(BigInt, min(left, right) / spacing - offset)
        last = floor(BigInt, max(left, right) / spacing - offset)
        push!(
            pieces,
            (
                tile = tile,
                start = left,
                finish = right,
                first = first,
                last = last,
                offset = offset,
            ),
        )
    end
    return (
        direction = direction,
        seam = floor(BigInt, x / period) != floor(BigInt, y / period),
        segments = Tuple(pieces),
    )
end
function _ec_stencil_receipt(ps, g, shape)
    D = length(g.n)
    return [
        ntuple(d -> begin
            coordinate = _particle_cell_position(ps, g, d, p)
            base, weights = _stencil1d(shape, coordinate)
            (coordinate = coordinate, base = base, weights = weights)
        end, D) for p in eachindex(ps.weight)
    ]
end
function _ec_particle_branch!(
    destination,
    origin,
    euler,
    forces0,
    forces1,
    stage::Union{EMPIC,EMPIC1D},
    mobile,
    h,
    which,
)
    T = eltype(stage.rho_n)
    D = length(stage.g.n)
    H = _exact_rational(h)
    rows = Union{Nothing,NamedTuple}[]
    for (s, (before, after)) in enumerate(zip(origin, destination))
        if !mobile[s]
            push!(rows, nothing)
            continue
        end
        n = nparticles(before)
        ideal = ntuple(_ -> Vector{_EMPICCoupledQ}(undef, n), D)
        lifted = ntuple(_ -> Vector{_EMPICCoupledQ}(undef, n), D)
        winding = ntuple(_ -> Vector{BigInt}(undef, n), D)
        defect = ntuple(_ -> Vector{_EMPICCoupledQ}(undef, n), D)
        wrapped = ntuple(_ -> Vector{T}(undef, n), D)
        cuts = Vector{NTuple{D,NamedTuple}}(undef, n)
        for p in eachindex(before.weight)
            v0 = ntuple(c -> before.v[c][p], 3)
            p1 = forces1 === nothing ? nothing : forces1[s].increments[p]
            output = _em_lorentz_pair_output(
                v0,
                forces0[s].increments[p],
                p1,
                stage.c,
                Val(stage.relativistic),
            )
            all(isfinite, output) || _ec_reject(:nonfinite_velocity, which, s, p)
            stage.relativistic &&
                !_supported_relativistic_output(output, stage.c) &&
                _ec_reject(:unsupported_velocity, which, s, p)
            for c = 1:3
                after.v[c][p] = output[c]
            end
            for d = 1:D
                x = _exact_rational(before.x[d][p])
                velocity =
                    forces1 === nothing ? _exact_rational(v0[d]) :
                    (_exact_rational(v0[d]) + _exact_rational(euler[s].v[d][p])) / 2
                target = _particle_coordinate(x, velocity, H)
                represented = _periodic_represented_lift(target, stage.g.L[d])
                bound = _exact_rational(stage.g.L[d]) / 4
                max(abs(target - x), abs(represented.lifted - x)) < bound || _ec_reject(
                    :travel,
                    which,
                    s,
                    p,
                    max(abs(target - x), abs(represented.lifted - x)) / bound,
                )
                after.x[d][p] = represented.wrapped
                ideal[d][p], lifted[d][p] = target, represented.lifted
                winding[d][p], defect[d][p], wrapped[d][p] =
                    represented.winding, represented.defect, represented.wrapped
            end
            cuts[p] = ntuple(
                d -> _ec_cut_receipt(
                    _exact_rational(before.x[d][p]),
                    lifted[d][p],
                    stage.g.L[d],
                    stage.g.dx[d],
                    stage.shape,
                ),
                D,
            )
        end
        push!(
            rows,
            (
                ideal = ideal,
                lifted = lifted,
                winding = winding,
                defect = defect,
                wrapped = wrapped,
                cuts = cuts,
                stencil = _ec_stencil_receipt(after, stage.g, stage.shape),
            ),
        )
    end
    return rows
end
function _ec_constraint(stage::Union{EMPIC,EMPIC1D}, meter::_EMPICCoupledMeter)
    T = eltype(stage.rho_n)
    E = stage isa EMPIC ? stage.E : (stage.Ex, stage.Ey)
    B = stage isa EMPIC ? stage.B : (stage.Bz,)
    projectedE =
        stage isa EMPIC ? map(copy, E) : (copy(stage.Ex), copy(stage.Ey), zeros(T, stage.g.n))
    projectedB =
        stage isa EMPIC ? map(copy, B) : (zeros(T, stage.g.n), zeros(T, stage.g.n), copy(stage.Bz))
    _empic_fft_charge!(meter, 7)
    _spectral_divergence_project!(projectedE, stage.rho_np1, stage.g)
    _empic_fft_charge!(meter, 6)
    _spectral_divergence_project!(projectedB, nothing, stage.g)
    e = ntuple(
        c -> map((a, b) -> _exact_rational(a) - _exact_rational(b), projectedE[c], E[c]),
        length(E),
    )
    b = ntuple(
        c -> begin
            projected = stage isa EMPIC ? projectedB[c] : projectedB[3]
            map((a, x) -> _exact_rational(a) - _exact_rational(x), projected, B[c])
        end,
        length(B),
    )
    return (electric = e, magnetic = b)
end
function _ec_component_norm(difference, start, low, high, absolute, relative)
    result = zero(_EMPICCoupledQ)
    for i in eachindex(difference)
        scale =
            _exact_rational(absolute) +
            _exact_rational(relative) * max(
                abs(_exact_rational(start[i])),
                abs(_exact_rational(low[i])),
                abs(_exact_rational(high[i])),
            )
        result = max(result, abs(difference[i]) / scale)
    end
    return result
end
function _ec_field_norm(difference, start, low, high, control, role)
    result = zero(_EMPICCoupledQ)
    for c in eachindex(difference)
        result = max(
            result,
            _ec_component_norm(
                difference[c],
                start[c],
                low[c],
                high[c],
                getproperty(control.atol, role),
                getproperty(control.rtol, role),
            ),
        )
    end
    return result
end
function _ec_observed_constraint_norm(constraint, origin, low, high, control)
    E0 = _empic_event_field_E(low, origin)
    B0 = _empic_event_field_B(low, origin)
    return (
        electric = _ec_field_norm(
            constraint.electric,
            E0,
            _empic_event_E(low),
            _empic_event_E(high),
            control,
            :E,
        ),
        magnetic = _ec_field_norm(
            constraint.magnetic,
            B0,
            _empic_event_B(low),
            _empic_event_B(high),
            control,
            :B,
        ),
    )
end
function _ec_particle_norm(origin, low, high, lifts, mobile, g, control)
    position, velocity = zero(_EMPICCoupledQ), zero(_EMPICCoupledQ)
    for (s, (before, left, right)) in enumerate(zip(origin, low, high))
        mobile[s] || continue
        for p in eachindex(before.weight)
            for d in eachindex(before.x)
                delta = lifts.high[s].lifted[d][p] - lifts.low[s].lifted[d][p]
                abs(delta) < _exact_rational(g.L[d]) / 2 ||
                    _ec_reject(:ambiguous_position, :pair, s, p)
                position = max(
                    position,
                    abs(delta) / (
                        _exact_rational(control.atol.x) +
                        _exact_rational(control.rtol.x) * _exact_rational(g.L[d])
                    ),
                )
            end
            for c = 1:3
                scale =
                    _exact_rational(control.atol.v) +
                    _exact_rational(control.rtol.v) * max(
                        abs(_exact_rational(before.v[c][p])),
                        abs(_exact_rational(left.v[c][p])),
                        abs(_exact_rational(right.v[c][p])),
                    )
                velocity = max(
                    velocity,
                    abs(_exact_rational(right.v[c][p]) - _exact_rational(left.v[c][p])) / scale,
                )
            end
        end
    end
    return position, velocity
end
function _ec_row_receipts(origin, out, stages, weights, h)
    H = _exact_rational(h)
    c2 = _exact_rational(out.c)^2
    E0, B0 = _empic_event_field_E(out, origin), _empic_event_field_B(out, origin)
    E, B, J = _empic_event_E(out), _empic_event_B(out), _empic_event_J(out)
    etaE = ntuple(length(E)) do c
        parts = _empic_weighted_derivative(stages, weights, :magnetic, c)
        result = similar(E[c], _EMPICCoupledQ)
        for i in eachindex(result)
            result[i] =
                _exact_rational(E[c][i]) - _exact_rational(E0[c][i]) -
                H * c2 * _em_exact_derivative(parts, i) + H * _exact_rational(J[c][i])
        end
        result
    end
    etaB = ntuple(length(B)) do c
        parts = _empic_weighted_derivative(stages, weights, :electric, c)
        result = similar(B[c], _EMPICCoupledQ)
        for i in eachindex(result)
            result[i] =
                _exact_rational(B[c][i]) - _exact_rational(B0[c][i]) +
                H * _em_exact_derivative(parts, i)
        end
        result
    end
    return (electric = etaE, magnetic = etaB)
end
function _ec_continuity(stage, h, meter)
    meter.continuity = nothing
    relative = charge_conservation_residual(stage, h, meter)
    isfinite(relative) && relative >= 0 && meter.continuity !== nothing ||
        _ec_reject(:unresolved_continuity, :field)
    return merge((relative = relative,), meter.continuity)
end
function _ec_stage_power(species, forces, mobile)
    result = zero(_EMPICCoupledQ)
    for (s, ps) in enumerate(species)
        mobile[s] || continue
        for p in eachindex(ps.weight)
            result +=
                _exact_rational(ps.q) *
                _exact_rational(ps.weight[p]) *
                sum(c -> _exact_rational(ps.v[c][p]) * _exact_rational(forces[s].E[c][p]), 1:3)
        end
    end
    return result
end
function _ec_interval_work(
    stage::Union{EMPIC,EMPIC1D},
    origin_species,
    low_species,
    high_species,
    origin_fields,
    forces0,
    forces1,
    raw,
    h,
    mobile,
)
    H = _exact_rational(h)
    volume = prod(_exact_rational, stage.g.dx)
    E0 = _empic_event_field_E(stage, origin_fields)
    E1 = _empic_event_E(stage)
    J = _empic_event_J(stage)
    used =
        volume *
        H *
        sum(
            c -> sum(
                i ->
                    _exact_rational(J[c][i]) *
                    (_exact_rational(E0[c][i]) + _exact_rational(E1[c][i])) / 2,
                eachindex(J[c]);
                init = zero(_EMPICCoupledQ),
            ),
            eachindex(J);
            init = zero(_EMPICCoupledQ),
        )
    rawwork = if stage isa EMPIC
        volume *
        H *
        sum(
            c -> sum(
                i ->
                    _exact_rational(raw.current[c][i]) *
                    (_exact_rational(E0[c][i]) + _exact_rational(E1[c][i])) / 2,
                eachindex(raw.current[c]);
                init = zero(_EMPICCoupledQ),
            ),
            eachindex(raw.current);
            init = zero(_EMPICCoupledQ),
        )
    else
        volume *
        H *
        sum(
            i ->
                _exact_rational(only(raw.current)[i]) *
                (_exact_rational(E0[2][i]) + _exact_rational(E1[2][i])) / 2,
            eachindex(stage.Jy);
            init = zero(_EMPICCoupledQ),
        )
    end
    particle =
        H * (
            _ec_stage_power(origin_species, forces0, mobile) +
            _ec_stage_power(low_species, forces1, mobile)
        ) / 2
    K1, U1 = _ec_energy(stage, high_species)
    K0 =
        stage.relativistic ?
        begin
            value = first(_relativistic_particle_moments(origin_species, stage.c))
            isfinite(value) ||
                throw(ArgumentError("coupled initial reported kinetic energy is nonfinite"))
            _exact_rational(value)
        end : _kinetic_energy_rational(origin_species)
    U0 = _field_energy_rational(
        E0,
        _empic_event_field_B(stage, origin_fields),
        _exact_rational(stage.c)^2,
        stage.g.dx,
    )
    df, dp, dc = U1 - U0 + used, K1 - K0 - particle, particle - used
    return (
        work_used = used,
        work_raw = stage isa EMPIC ? rawwork : nothing,
        work_raw_transverse = stage isa EMPIC1D ? rawwork : nothing,
        particle_work = particle,
        field_before = U0,
        field_after = U1,
        kinetic_before = K0,
        kinetic_after = K1,
        field_defect = df,
        particle_defect = dp,
        coupling_defect = dc,
        total_defect = U1 - U0 + K1 - K0,
    )
end
function _ec_pair!(
    stage::Union{EMPIC,EMPIC1D},
    origin_species::AbstractVector,
    origin_fields::NamedTuple,
    mobile::BitVector,
    target::T,
    h::T,
    delta::_EMPICCoupledQ,
    time::T,
    next_time::T,
    source::_EMPICCoupledQ,
    elapsed::_EMPICCoupledQ,
    control::NamedTuple,
    meter::_EMPICCoupledMeter,
) where {T<:AbstractFloat}
    T = eltype(stage.rho_n)
    D = length(stage.g.n)
    H = _exact_rational(h)
    before_cost = _ec_costs(meter)
    _ec_charge!(meter, :trials, 1)
    _ec_restore_fields!(stage, origin_fields)
    rho0 = similar(stage.rho_n)
    _ec_source!(rho0, stage, origin_species, target, meter)
    isequal(rho0, origin_fields.rho_np1) ||
        throw(ArgumentError("coupled origin charge differs from its committed source"))
    force0 = _ec_force!(stage, origin_species, mobile, H, meter)
    raw0 = _ec_current!(stage, origin_species, mobile, meter)
    derivatives0 = _ec_derivatives!(stage, meter)
    low = _ec_new_stage(stage)
    _ec_restore_fields!(low, origin_fields)
    low_species = _ec_copy_species(origin_species, Val(D), T)
    low_lift = _ec_particle_branch!(
        low_species,
        origin_species,
        nothing,
        force0,
        nothing,
        stage,
        mobile,
        h,
        :low,
    )
    copyto!(low.rho_n, rho0)
    _ec_source!(low.rho_np1, low, low_species, target, meter)
    _ec_correct!(low, raw0, h, meter, :low)
    _empic_full_time_field_update!(low, origin_fields, (derivatives0,), (one(_EMPICCoupledQ),), h)
    all(
        a -> all(isfinite, a),
        (_empic_event_E(low)..., _empic_event_B(low)..., _empic_event_J(low)...),
    ) || _ec_reject(:nonfinite_field, :low)
    constraint_low = _ec_constraint(low, meter)
    cl = _ec_observed_constraint_norm(constraint_low, origin_fields, low, low, control)
    max(cl.electric, cl.magnetic) <= 1 ||
        _ec_reject(:field_constraint, :low, 0, 0, max(cl.electric, cl.magnetic))
    continuity_low = _ec_continuity(low, h, meter)
    force1 = _ec_force!(low, low_species, mobile, H, meter)
    raw1 = _ec_current!(low, low_species, mobile, meter)
    derivatives1 = _ec_derivatives!(low, meter)
    # Gathering/current scratch cannot replace the low branch's consumed current.
    low_snapshot = _ec_snapshot(low, low_species)
    if low isa EMPIC
        foreach(copyto!, low.Jraw, raw0.current)
        low_snapshot = merge(
            low_snapshot,
            (fields = merge(low_snapshot.fields, (Jraw = map(copy, raw0.current),)),),
        )
    end
    high_species = _ec_copy_species(origin_species, Val(D), T)
    high_lift = _ec_particle_branch!(
        high_species,
        origin_species,
        low_species,
        force0,
        force1,
        stage,
        mobile,
        h,
        :high,
    )
    rawhigh = _ec_average_current(raw0, raw1, T)
    copyto!(stage.rho_n, rho0)
    _ec_source!(stage.rho_np1, stage, high_species, target, meter)
    _ec_correct!(stage, rawhigh, h, meter, :high)
    weights = (BigInt(1) // BigInt(2), BigInt(1) // BigInt(2))
    stages = (derivatives0, derivatives1)
    _empic_full_time_field_update!(stage, origin_fields, stages, weights, h)
    all(
        a -> all(isfinite, a),
        (_empic_event_E(stage)..., _empic_event_B(stage)..., _empic_event_J(stage)...),
    ) || _ec_reject(:nonfinite_field, :high)
    constraint_high = _ec_constraint(stage, meter)
    continuity_high = _ec_continuity(stage, h, meter)
    lifts = (low = low_lift, high = high_lift)
    nx, nv = _ec_particle_norm(
        origin_species,
        low_species,
        high_species,
        lifts,
        mobile,
        stage.g,
        control,
    )
    E0 = _empic_event_field_E(stage, origin_fields)
    B0 = _empic_event_field_B(stage, origin_fields)
    de = map(
        (a, b) -> map((x, y) -> _exact_rational(x) - _exact_rational(y), a, b),
        _empic_event_E(stage),
        _empic_event_E(low),
    )
    db = map(
        (a, b) -> map((x, y) -> _exact_rational(x) - _exact_rational(y), a, b),
        _empic_event_B(stage),
        _empic_event_B(low),
    )
    ne = _ec_field_norm(de, E0, _empic_event_E(low), _empic_event_E(stage), control, :E)
    nb = _ec_field_norm(db, B0, _empic_event_B(low), _empic_event_B(stage), control, :B)
    scale =
        _exact_rational(control.atol.clock) +
        _exact_rational(control.rtol.clock) * max(
            abs(_exact_rational(time)),
            abs(_exact_rational(next_time)),
            abs(source),
            abs(source + delta),
        )
    duration_defect = H - delta
    addition_defect = _exact_rational(next_time) - _exact_rational(time) - H
    nc =
        max(abs(duration_defect), abs(addition_defect), abs(duration_defect + addition_defect)) /
        scale
    ch = _ec_observed_constraint_norm(constraint_high, origin_fields, low, stage, control)
    norms = (
        x = nx,
        v = nv,
        E = ne,
        B = nb,
        clock = nc,
        constraint_E_low = cl.electric,
        constraint_B_low = cl.magnetic,
        constraint_E_high = ch.electric,
        constraint_B_high = ch.magnetic,
    )
    norm = maximum(values(norms))
    norm <= 1 || _ec_reject(:embedded_defect, :pair, 0, 0, norm)
    work = _ec_interval_work(
        stage,
        origin_species,
        low_species,
        high_species,
        origin_fields,
        force0,
        force1,
        rawhigh,
        h,
        mobile,
    )
    interval = (
        dt = h,
        source_before = source,
        source_after = source + delta,
        elapsed_before = elapsed,
        elapsed_after = elapsed + H,
        public_before = time,
        public_after = next_time,
        duration_defect = duration_defect,
        addition_defect = addition_defect,
        origin = (
            species = NamedTuple[deepcopy(_checkpoint_particle_state(ps)) for ps in origin_species],
            fields = _ec_copy_fields(origin_fields),
        ),
        low = low_snapshot,
        high = _ec_snapshot(stage, high_species),
        stage0 = (
            forces = force0,
            raw = raw0,
            derivatives = derivatives0,
            stencil = [_ec_stencil_receipt(ps, stage.g, stage.shape) for ps in origin_species],
        ),
        stage1 = (
            forces = force1,
            raw = raw1,
            derivatives = derivatives1,
            stencil = [_ec_stencil_receipt(ps, stage.g, stage.shape) for ps in low_species],
        ),
        raw_effective = rawhigh,
        lifts = lifts,
        norms = norms,
        norm = norm,
        constraint_low = constraint_low,
        constraint_high = constraint_high,
        continuity = (low = continuity_low, high = continuity_high),
        rows = (
            low = _ec_row_receipts(origin_fields, low, (derivatives0,), (one(_EMPICCoupledQ),), h),
            high = _ec_row_receipts(origin_fields, stage, stages, weights, h),
        ),
        work = work,
        costs = _ec_cost_sub(_ec_costs(meter), before_cost),
    )
    return high_species, interval
end
function _ec_role_arrays(value, result = AbstractArray[], descended = IdDict{Any,Bool}(); depth = 0)
    depth <= 32 || throw(ArgumentError("coupled closed receipt exceeds its schema nesting bound"))
    if value isa AbstractArray
        push!(result, value)
        # Repeated roles remain visible to independence checks, but a cycle
        # must not recursively enumerate the same container forever.
        haskey(descended, value) && return result
        descended[value] = true
        if !(eltype(value) <: Number)
            for entry in value
                _ec_role_arrays(entry, result, descended; depth = depth + 1)
            end
        end
    elseif value isa NamedTuple || value isa Tuple
        for entry in value
            _ec_role_arrays(entry, result, descended; depth = depth + 1)
        end
    elseif value isa ParticleSet
        append!(result, _empic_event_arrays(value))
    end
    return result
end
function _ec_independent(arrays)
    for (i, a) in enumerate(arrays), j = 1:i-1
        b = arrays[j]
        left, right = _checkpoint_packed_chunks(a), _checkpoint_packed_chunks(b)
        (
            _particle_arrays_mightalias(a, b) ||
            (left !== nothing && _particle_arrays_mightalias(left, b)) ||
            (right !== nothing && _particle_arrays_mightalias(a, right)) ||
            (left !== nothing && right !== nothing && _particle_arrays_mightalias(left, right))
        ) && throw(ArgumentError("coupled particle/field/controller/report storage roles overlap"))
    end
    return nothing
end
function _ec_report_cells(value)
    count = 0
    for array in _ec_role_arrays(value)
        count = Base.checked_add(count, length(array))
    end
    return count
end
function _empic_collection_schedule_arrays(schedule::_EMPICCoupledSchedule)
    arrays =
        AbstractArray[schedule.proposal.cap, schedule.proposal.adaptive, schedule.proposal.minimum]
    for phase in (schedule.phase, schedule.candidate)
        append!(arrays, (phase.interval, phase.last, phase.next))
        append!(arrays, _ec_role_arrays(phase.report))
    end
    append!(arrays, _empic_event_shadow_arrays(schedule.shadow))
    append!(arrays, _empic_event_shadow_arrays(schedule.candidate_shadow))
    return arrays
end
_empic_collection_schedule_axes(::_EMPICCoupledSchedule, array, count) = nothing
function _empic_collection_validate_schedule(es, schedule::_EMPICCoupledSchedule{T}) where {T}
    _ec_mode(T)
    es.shape isa Union{CIC,TSC} || throw(ArgumentError("coupled error control admits CIC/TSC only"))
    _ec_identical(schedule.control, _ec_control(schedule.control, T)) ||
        throw(ArgumentError("coupled control representation changed"))
    p = schedule.proposal
    _ec_identical(
        p,
        _empic_event_proposal(
            es.collection.mobile,
            p.cap,
            p.adaptive,
            p.minimum,
            p.limits,
            p.budget,
            T,
        ),
    ) || throw(ArgumentError("coupled nominal proposal representation changed"))
    es.n_sub isa Int && 1 <= es.n_sub <= p.budget ||
        throw(ArgumentError("coupled base/event budget differs"))
    schedule.signature isa String && schedule.signature == _ec_spectral_signature(es.g) ||
        throw(ArgumentError("coupled spectral symbol changed"))
    for phase in (schedule.phase, schedule.candidate)
        all(
            a -> axes(a, 1) == Base.OneTo(length(es.collection.mobile)),
            (phase.interval, phase.last, phase.next),
        ) || throw(ArgumentError("coupled calendar row layout changed"))
        phase.report === nothing ||
            phase.report isa NamedTuple ||
            throw(ArgumentError("coupled report container differs"))
    end
    return nothing
end
function _ec_configuration(
    es;
    proposal = es.collection.schedule.proposal,
    control = es.collection.schedule.control,
    n_sub = es.n_sub,
    c = es.c,
    relativistic = es.relativistic,
    signature = es.collection.schedule.signature,
)
    return merge(
        _empic_event_configuration(es; c, relativistic),
        (
            n_sub = n_sub,
            proposal = (
                cap = Tuple(proposal.cap),
                adaptive = Tuple(proposal.adaptive),
                minimum = Tuple(proposal.minimum),
                limits = proposal.limits,
                budget = proposal.budget,
            ),
            control = control,
            signature = signature,
        ),
    )
end
function _ec_shadow(es, populations, fields, bindings, clock; kwargs...)
    return _empic_event_shadow_values(
        es,
        populations,
        fields,
        bindings,
        _ec_configuration(es; kwargs...),
        clock,
    )
end
function _ec_shadow_matches(es, schedule)
    shadow = schedule.shadow
    work = es.collection
    shadow isa NamedTuple &&
        keys(shadow) == (:species, :fields, :bindings, :configuration, :clock) ||
        throw(ArgumentError("coupled committed shadow is absent"))
    _ec_identical(shadow.configuration, _ec_configuration(es)) &&
        isequal(shadow.clock, (es.time[], es.step[], es.last_dt[], es.last_sub_dt[])) &&
        isequal(shadow.fields, _empic_event_fields(es)) &&
        length(shadow.species) == length(shadow.bindings) == length(work.owners) || throw(
        ArgumentError(
            "coupled source/configuration/physical field/clock changed; initialize full-time data",
        ),
    )
    for (ps, saved, binding) in zip(work.owners, shadow.species, shadow.bindings)
        length(binding) == length(_empic_event_arrays(ps)) &&
            all(p -> p[1] === p[2], zip(_empic_event_arrays(ps), binding)) &&
            isequal((ps.q, ps.m), (saved.q, saved.m)) &&
            all(p -> isequal(p[1], p[2]), zip(_empic_event_arrays(ps), _empic_event_arrays(saved))) ||
            throw(
                ArgumentError(
                    "coupled source/provenance/component identity changed; initialize full-time data",
                ),
            )
    end
    return nothing
end
function _ec_phase_scalar_check(es, phase, population_count::Int, initialized, clock)
    T = eltype(es.rho_n)
    time, step, full, fast = clock
    phase isa _EMPICCoupledPhase{T} && initialized isa Bool ||
        throw(ArgumentError("coupled phase type differs"))
    _validate_empic_clock_phase(time, step, full, fast)
    all(
            _ec_q,
            (
                phase.time,
                phase.elapsed,
                phase.initial_energy,
                phase.field_defect,
                phase.particle_defect,
                phase.coupling_defect,
                phase.last...,
                phase.next...,
            ),
        ) &&
        phase.time >= 0 &&
        phase.elapsed >= 0 &&
        phase.initial_energy >= 0 &&
        isfinite(phase.next_h) &&
        _stored_positive(phase.next_h) || throw(ArgumentError("coupled phase exact scalars differ"))
    _require_resolved_dt("coupled next interval", phase.next_h)
    _ec_costs_check(phase.costs)
    _ec_identical(phase.admission_control, _ec_control(phase.admission_control, T)) ||
        throw(ArgumentError("coupled historical admission control differs"))
    count = length(es.collection.mobile)
    population_count == count &&
        all(a -> axes(a, 1) == Base.OneTo(count), (phase.interval, phase.last, phase.next)) ||
        throw(ArgumentError("coupled phase population rows differ"))
    all(isfinite, phase.interval) || throw(ArgumentError("coupled nominal intervals are nonfinite"))
    if step == 0
        _particle_stored_zero(time) &&
            iszero(phase.time) &&
            iszero(phase.elapsed) &&
            phase.report === nothing &&
            all(_particle_stored_zero, phase.interval) &&
            all(iszero, (phase.last..., phase.next...)) &&
            all(iszero, (phase.field_defect, phase.particle_defect, phase.coupling_defect)) ||
            throw(ArgumentError("coupled step-zero phase is advanced"))
    else
        initialized && phase.report !== nothing && phase.time > 0 && phase.elapsed > 0 ||
            throw(ArgumentError("coupled advanced phase lacks a complete report"))
        for i in eachindex(es.collection.mobile)
            if es.collection.mobile[i]
                phase.interval[i] >= 2floatmin(T) &&
                    phase.last[i] <= phase.time <= phase.next[i] &&
                    phase.next[i] - phase.last[i] == _exact_rational(phase.interval[i]) ||
                    throw(ArgumentError("coupled nominal deadline history differs"))
            else
                _particle_stored_zero(phase.interval[i]) &&
                    iszero(phase.last[i]) &&
                    iszero(phase.next[i]) ||
                    throw(ArgumentError("fixed coupled population has an active calendar"))
            end
        end
    end
    return nothing
end
function _ec_validate_phase(
    es::Union{EMPIC,EMPIC1D},
    phase::_EMPICCoupledPhase,
    species::AbstractVector,
    fields::NamedTuple,
    initialized::Bool,
    clock::Tuple,
    meter::_EMPICCoupledMeter;
    replay::Bool = true,
)
    _ec_phase_scalar_check(es, phase, length(species), initialized, clock)
    if initialized
        if iszero(clock[2])
            currents = es isa EMPIC ? (fields.J..., fields.Jraw...) : (fields.Jx, fields.Jy)
            _ec_identical(fields.rho_n, fields.rho_np1) &&
                all(array -> all(_particle_stored_zero, array), currents) || throw(
                ArgumentError(
                    "initialized coupled step-zero source/current inventory differs from full-time initialization",
                ),
            )
        end
        probe = _ec_new_stage(es)
        _ec_restore_fields!(probe, fields)
        K, U = _ec_energy(probe, species)
        K + U - phase.initial_energy ==
        phase.field_defect + phase.particle_defect + phase.coupling_defect ||
            throw(ArgumentError("coupled cumulative measured-energy ledger differs"))
        charges = _pic_collection_charges(es.g, es.n0, species)
        target = _empic_collection_charge_target(es, es.collection, charges)
        rho = similar(probe.rho_np1)
        _ec_source!(rho, probe, species, target, meter)
        isequal(rho, fields.rho_np1) || throw(ArgumentError("coupled ending charge/source differs"))
        constraint = _ec_constraint(probe, meter)
        control = phase.admission_control
        observed = _ec_observed_constraint_norm(constraint, fields, probe, probe, control)
        max(observed.electric, observed.magnetic) <= 1 ||
            throw(ArgumentError("coupled ending field constraint differs"))
        phase.report === nothing ||
            _ec_report_check(es, phase, species, fields, clock, meter; replay)
    end
    return nothing
end
function _empic_collection_validate_species_history(es, schedule::_EMPICCoupledSchedule)
    # Public controlled operations own the full observed/replayed admission.
    count =
        es.collection.initialized[] ? length(es.collection.owners) : length(es.collection.mobile)
    _ec_phase_scalar_check(
        es,
        schedule.phase,
        count,
        es.collection.initialized[],
        (es.time[], es.step[], es.last_dt[], es.last_sub_dt[]),
    )
    es.collection.initialized[] && _ec_shadow_matches(es, schedule)
    return nothing
end
function _ec_ready(
    es::Union{EMPIC,EMPIC1D},
    species::AbstractVector,
    meter::_EMPICCoupledMeter,
)
    charges = _empic_collection_validate_species(es, species)
    work = es.collection
    work.initialized[] ||
        throw(ArgumentError("initialize coupled full-time data before evolution or observation"))
    _empic_collection_bound(work, species)
    _ec_shadow_matches(es, work.schedule)
    _empic_collection_require_finite_fields(es)
    _ec_validate_phase(
        es,
        work.schedule.phase,
        species,
        _empic_event_fields(es),
        true,
        (es.time[], es.step[], es.last_dt[], es.last_sub_dt[]),
        meter,
    )
    return charges
end
function _ec_sync_magnetic!(stage::EMPIC)
    foreach(copyto!, stage.Bc, stage.B)
    return nothing
end
function _ec_sync_magnetic!(stage::EMPIC1D)
    copyto!(stage.Bzc, stage.Bz)
    return nothing
end
function _ec_publish!(es, work, phase, shadow; particles)
    schedule = work.schedule
    reset = _ec_zero_phase(eltype(es.rho_n), length(work.mobile), schedule.control)
    schedule.candidate, schedule.candidate_shadow = phase, shadow
    if particles
        for (ps, private) in zip(work.owners, work.staged)
            foreach(copyto!, ps.x, private.x)
            foreach(copyto!, ps.v, private.v)
        end
    end
    foreach(copyto!, _empic_transaction_fields(es), _empic_transaction_fields(work.stage))
    es.time[], es.step[] = work.stage.time[], work.stage.step[]
    es.last_dt[], es.last_sub_dt[] = work.stage.last_dt[], work.stage.last_sub_dt[]
    schedule.phase, schedule.shadow = phase, shadow
    schedule.candidate, schedule.candidate_shadow = reset, nothing
    return es
end
function _empic_coupled_init!(es::Union{EMPIC,EMPIC1D}, species::AbstractVector)
    work = es.collection
    schedule = work.schedule
    _ec_acquire(schedule)
    try
        T = eltype(es.rho_n)
        all(ps -> ps isa eltype(work.owners), species) ||
            throw(ArgumentError("coupled population storage type differs"))
        charges = _empic_collection_validate_species(es, species)
        all(a -> all(isfinite, a), _empic_collection_seed_fields(es)) ||
            throw(ArgumentError("coupled full-time seed fields must be finite"))
        _empic_collection_prepare!(es, species)
        meter = _ec_meter(schedule.control)
        stage = work.stage
        # EMPIC1D derives Ex from charge, which would erase a seeded uniform
        # (zero-k) longitudinal mode. The controller contract retains excluded
        # homogeneous modes, so preserve the original mean and restore it after
        # the charge-gauge projection. Nonzero-k Ex remains charge-derived.
        seeded_dc =
            stage isa EMPIC1D ? _from_exact(T, _exact_float_sum(stage.Ex) / length(stage.Ex)) :
            zero(T)
        target = _empic_collection_charge_target(es, work, charges)
        _ec_source!(stage.rho_n, stage, work.staged, target, meter)
        _empic_fft_charge!(meter, stage isa EMPIC ? 7 : 2)
        _empic_collection_initialize_fields!(stage)
        if stage isa EMPIC1D && !iszero(seeded_dc)
            @. stage.Ex += seeded_dc
        end
        stage.time[], stage.step[] = zero(T), 0
        stage.last_dt[], stage.last_sub_dt[] = zero(T), zero(T)
        K, U = _ec_energy(stage, work.staged)
        phase = _ec_zero_phase(
            T,
            length(species),
            schedule.control;
            energy = K + U,
            costs = _ec_costs(meter),
        )
        _ec_validate_phase(
            es,
            phase,
            work.staged,
            _empic_event_fields(stage),
            true,
            (zero(T), 0, zero(T), zero(T)),
            meter,
        )
        phase = _EMPICCoupledPhase{T}(
            phase.time,
            phase.elapsed,
            phase.interval,
            phase.last,
            phase.next,
            phase.next_h,
            phase.initial_energy,
            phase.field_defect,
            phase.particle_defect,
            phase.coupling_defect,
            phase.admission_control,
            _ec_costs(meter),
            nothing,
        )
        shadow = _ec_shadow(
            es,
            work.staged,
            _empic_event_fields(stage),
            species,
            (zero(T), 0, zero(T), zero(T)),
        )
        _ec_independent(
            vcat(
                _ec_role_arrays(phase.report),
                _empic_event_shadow_arrays(shadow),
                _ec_role_arrays((phase.interval, phase.last, phase.next)),
            ),
        )
        _empic_collection_validate_storage(es, species)
        sizehint!(work.owners, length(species))
        _ec_publish!(es, work, phase, shadow; particles = false)
        empty!(work.owners)
        append!(work.owners, species)
        work.initialized[] = true
        return es
    finally
        _ec_release(schedule)
    end
end
function _empic_coupled_set!(es, n_sub, stride, cap, adaptive, minimum, limits, budget, control)
    schedule = es.collection.schedule
    _ec_acquire(schedule)
    try
        stride === nothing || throw(ArgumentError("coupled control cannot use species_stride"))
        n_sub isa Integer && !(n_sub isa Bool) && 1 <= n_sub <= typemax(Int) ||
            throw(ArgumentError("coupled n_sub must fit a positive Int"))
        T = eltype(es.rho_n)
        old = schedule.proposal
        proposed = _empic_event_proposal(
            es.collection.mobile,
            cap === nothing ? old.cap : cap,
            adaptive === nothing ? old.adaptive : adaptive,
            minimum === nothing ? old.minimum : minimum,
            limits === nothing ? old.limits : limits,
            budget === nothing ? old.budget : budget,
            T,
        )
        n_sub <= proposed.budget ||
            throw(ArgumentError("coupled base field partition exceeds event budget"))
        closed = control === nothing ? schedule.control : _ec_control(control, T)
        _empic_collection_validate_schedule(es, schedule)
        es.collection.initialized[] &&
            _ec_ready(es, es.collection.owners, _ec_meter(schedule.phase.admission_control))
        phase = schedule.phase
        next_h = _ec_future_hint(phase, closed)
        candidate = _EMPICCoupledPhase{T}(
            phase.time,
            phase.elapsed,
            copy(phase.interval),
            deepcopy(phase.last),
            deepcopy(phase.next),
            next_h,
            phase.initial_energy,
            phase.field_defect,
            phase.particle_defect,
            phase.coupling_defect,
            phase.admission_control,
            phase.costs,
            deepcopy(phase.report),
        )
        shadow =
            es.collection.initialized[] ?
            _ec_shadow(
                es,
                es.collection.owners,
                _empic_event_fields(es),
                es.collection.owners,
                (es.time[], es.step[], es.last_dt[], es.last_sub_dt[]);
                proposal = proposed,
                control = closed,
                n_sub = Int(n_sub),
            ) : nothing
        _ec_independent(
            vcat(
                _ec_role_arrays((
                    proposed,
                    candidate.interval,
                    candidate.last,
                    candidate.next,
                    candidate.report,
                )),
                _empic_event_shadow_arrays(shadow),
            ),
        )
        es.n_sub = Int(n_sub)
        schedule.proposal, schedule.control = proposed, closed
        schedule.phase, schedule.shadow = candidate, shadow
        return es
    finally
        _ec_release(schedule)
    end
end
function _ec_recommend(h, norm, control)
    exact = norm <= BigInt(1) // BigInt(4) ? 2 * _exact_rational(h) : _exact_rational(h)
    return _empic_event_floor(
        typeof(h),
        min(_exact_rational(control.h_max), max(_exact_rational(control.h_min), exact)),
    )
end
function _ec_future_hint(phase, control)
    hint =
        phase.report === nothing ? phase.admission_control.h_max :
        _ec_recommend(
            last(phase.report.intervals).dt,
            last(phase.report.intervals).norm,
            phase.report.control,
        )
    return min(control.h_max, max(control.h_min, hint))
end
const _EC_REPORT_KEYS = (
    :kind,
    :control,
    :event_proposal,
    :configuration,
    :n_sub,
    :start,
    :finish,
    :requested,
    :elapsed_before,
    :elapsed_after,
    :public_before,
    :public_after,
    :step_before,
    :step_after,
    :initial_calendar,
    :initial_next_h,
    :cost_before,
    :costs,
    :cost_after,
    :admission_costs,
    :verification_costs,
    :initial_energy,
    :defects_before,
    :intervals,
    :rejections,
    :decisions,
    :work_used,
    :work_raw,
    :work_raw_transverse,
    :particle_work,
    :field_before,
    :field_after,
    :kinetic_before,
    :kinetic_after,
    :field_defect,
    :particle_defect,
    :coupling_defect,
    :total_defect,
    :clock_defect,
)
const _EC_INTERVAL_KEYS = (
    :dt,
    :source_before,
    :source_after,
    :elapsed_before,
    :elapsed_after,
    :public_before,
    :public_after,
    :duration_defect,
    :addition_defect,
    :origin,
    :low,
    :high,
    :stage0,
    :stage1,
    :raw_effective,
    :lifts,
    :norms,
    :norm,
    :constraint_low,
    :constraint_high,
    :continuity,
    :rows,
    :work,
    :costs,
)
function _ec_nominal_decision!(
    stage,
    species,
    mobile,
    proposal,
    i,
    time,
    previous,
    bootstrap,
    meter,
)
    ps = species[i]
    _ec_charge!(meter, :force_evals, nparticles(ps))
    _empic_collection_gather!(stage, ps; initial = true)
    choice = _empic_event_interval_choice(stage, proposal, ps, i, previous, bootstrap)
    return merge((species = i, time = time, previous = previous, bootstrap = bootstrap), choice)
end
function _ec_sum_work(intervals, kind)
    first_work = first(intervals).work
    last_work = last(intervals).work
    total(name) = sum(i -> getproperty(i.work, name), intervals; init = zero(_EMPICCoupledQ))
    return (
        work_used = total(:work_used),
        work_raw = kind === :vector ? total(:work_raw) : nothing,
        work_raw_transverse = kind === :planar ? total(:work_raw_transverse) : nothing,
        particle_work = total(:particle_work),
        field_before = first_work.field_before,
        field_after = last_work.field_after,
        kinetic_before = first_work.kinetic_before,
        kinetic_after = last_work.kinetic_after,
        field_defect = total(:field_defect),
        particle_defect = total(:particle_defect),
        coupling_defect = total(:coupling_defect),
        total_defect = total(:total_defect),
    )
end
function _empic_coupled_step!(
    es::Union{EMPIC,EMPIC1D},
    species::AbstractVector,
    dt::Real,
)
    work = es.collection
    schedule = work.schedule
    _ec_acquire(schedule)
    try
        T = eltype(es.rho_n)
        D = length(es.g.n)
        duration = _require_resolved_dt(
            "coupled step",
            _validated_nonnegative_dt(T, dt; name = "coupled EMPIC"),
        )
        isfinite(dt) && dt >= 0 && (iszero(dt) || !_particle_stored_zero(duration)) ||
            throw(ArgumentError("coupled nonzero duration must survive conversion"))
        meter = _ec_meter(iszero(duration) ? schedule.phase.admission_control : schedule.control)
        charges = _ec_ready(es, species, meter)
        _particle_stored_zero(duration) && return es
        original = schedule.phase
        control = schedule.control
        proposal = schedule.proposal
        admission_costs = _ec_costs(meter)
        _empic_collection_prepare!(es, species)
        stage = work.stage
        start, finish = original.time, original.time + _exact_rational(duration)
        source, elapsed, clock_value = start, original.elapsed, stage.time[]
        H = copy(original.interval)
        a = deepcopy(original.last)
        b = deepcopy(original.next)
        _ec_charge!(meter, :report_cells, Base.checked_mul(6, length(work.mobile)))
        historical_proposal = deepcopy(proposal)
        initial_calendar = (interval = copy(H), last = deepcopy(a), next = deepcopy(b))
        intervals = NamedTuple[]
        rejections = NamedTuple[]
        decisions = NamedTuple[]
        base = 1
        next_h = original.next_h
        target = _empic_collection_charge_target(es, work, charges)
        while source < finish
            origin_species = _ec_copy_species(work.staged, Val(D), T)
            origin_fields = _ec_copy_fields(_empic_event_fields(stage))
            for i in eachindex(work.mobile)
                work.mobile[i] || continue
                bootstrap = es.step[] == 0 && isempty(intervals) && iszero(H[i])
                if bootstrap || b[i] == source
                    _ec_restore_fields!(stage, origin_fields)
                    decision = _ec_nominal_decision!(
                        stage,
                        origin_species,
                        work.mobile,
                        proposal,
                        i,
                        source,
                        H[i],
                        bootstrap,
                        meter,
                    )
                    length(decisions) < proposal.budget - length(intervals) ||
                        throw(ArgumentError("coupled nominal decision/event budget exhausted"))
                    _ec_charge!(meter, :report_cells, 1)
                    push!(decisions, decision)
                    H[i], a[i], b[i] =
                        decision.selected, source, source + _exact_rational(decision.selected)
                end
            end
            calendar = (time = start, next = b)
            edge, base = _empic_event_next(calendar, work.mobile, source, finish, base, es.n_sub)
            delta = min(_exact_rational(next_h), edge - source)
            while true
                delta >= _exact_rational(control.h_min) ||
                    throw(ArgumentError("coupled source gap fell below h_min"))
                h = _empic_event_interval(T, delta)
                control.h_min <= h <= control.h_max ||
                    throw(ArgumentError("coupled actual interval is outside h_min/h_max"))
                next_public, _ = _validated_step_clock(clock_value, es.step[], h)
                next_public > clock_value ||
                    throw(ArgumentError("coupled public clock cannot advance"))
                before = _ec_costs(meter)
                outcome = try
                    _ec_pair!(
                        stage,
                        origin_species,
                        origin_fields,
                        work.mobile,
                        target,
                        h,
                        delta,
                        clock_value,
                        next_public,
                        source,
                        elapsed,
                        control,
                        meter,
                    )
                catch failure
                    failure isa _EMPICCoupledReject || rethrow()
                    rejected = (
                        dt = h,
                        source = source,
                        delta = delta,
                        clock_value = clock_value,
                        reason = failure.reason,
                        stage = failure.stage,
                        species = failure.species,
                        particle = failure.particle,
                        defect = failure.defect,
                        costs = _ec_cost_sub(_ec_costs(meter), before),
                    )
                    _ec_charge!(meter, :report_cells, 1)
                    push!(rejections, rejected)
                    nothing
                end
                if outcome === nothing
                    delta /= 2
                    continue
                end
                high, interval = outcome
                length(intervals) < proposal.budget - length(decisions) ||
                    throw(ArgumentError("coupled accepted interval/event budget exhausted"))
                cells = Base.checked_add(1, _ec_report_cells(interval))
                _ec_charge!(meter, :report_cells, cells)
                push!(intervals, interval)
                empty!(work.staged)
                append!(work.staged, high)
                source, elapsed, clock_value =
                    source + delta, elapsed + _exact_rational(h), next_public
                stage.time[] = clock_value
                next_h = _ec_recommend(h, interval.norm, control)
                break
            end
        end
        _, count = _validated_step_clock(es.time[], es.step[], zero(T))
        stage.step[], stage.last_dt[], stage.last_sub_dt[] = count, duration, last(intervals).dt
        _ec_sync_magnetic!(stage)
        kind = es isa EMPIC ? :vector : :planar
        summary = _ec_sum_work(intervals, kind)
        report = (
            kind = kind,
            control = control,
            event_proposal = historical_proposal,
            configuration = _empic_event_configuration(es),
            n_sub = es.n_sub,
            start = start,
            finish = finish,
            requested = duration,
            elapsed_before = original.elapsed,
            elapsed_after = elapsed,
            public_before = es.time[],
            public_after = clock_value,
            step_before = es.step[],
            step_after = count,
            initial_calendar = initial_calendar,
            initial_next_h = original.next_h,
            cost_before = original.costs,
            costs = _ec_costs(meter),
            cost_after = _ec_cost_add(original.costs, _ec_costs(meter)),
            admission_costs = admission_costs,
            verification_costs = _ec_zero_costs(),
            initial_energy = original.initial_energy,
            defects_before = (
                field = original.field_defect,
                particle = original.particle_defect,
                coupling = original.coupling_defect,
            ),
            intervals = intervals,
            rejections = rejections,
            decisions = decisions,
            work_used = summary.work_used,
            work_raw = summary.work_raw,
            work_raw_transverse = summary.work_raw_transverse,
            particle_work = summary.particle_work,
            field_before = summary.field_before,
            field_after = summary.field_after,
            kinetic_before = summary.kinetic_before,
            kinetic_after = summary.kinetic_after,
            field_defect = summary.field_defect,
            particle_defect = summary.particle_defect,
            coupling_defect = summary.coupling_defect,
            total_defect = summary.total_defect,
            clock_defect = _exact_rational(clock_value) - _exact_rational(es.time[]) -
                           _exact_rational(duration),
        )
        _ec_report_cells(report) == meter.report_cells ||
            throw(ArgumentError("coupled report-cell reservation differs"))
        phase = _EMPICCoupledPhase{T}(
            finish,
            elapsed,
            H,
            a,
            b,
            next_h,
            original.initial_energy,
            original.field_defect + summary.field_defect,
            original.particle_defect + summary.particle_defect,
            original.coupling_defect + summary.coupling_defect,
            control,
            report.cost_after,
            report,
        )
        before_verification = _ec_costs(meter)
        _ec_validate_phase(
            es,
            phase,
            work.staged,
            _empic_event_fields(stage),
            true,
            (clock_value, count, duration, last(intervals).dt),
            meter,
        )
        verification = _ec_cost_sub(_ec_costs(meter), before_verification)
        report = merge(
            report,
            (
                costs = _ec_costs(meter),
                cost_after = _ec_cost_add(original.costs, _ec_costs(meter)),
                verification_costs = verification,
            ),
        )
        phase = _EMPICCoupledPhase{T}(
            finish,
            elapsed,
            H,
            a,
            b,
            next_h,
            original.initial_energy,
            phase.field_defect,
            phase.particle_defect,
            phase.coupling_defect,
            control,
            report.cost_after,
            report,
        )
        _ec_report_check(
            es,
            phase,
            work.staged,
            _empic_event_fields(stage),
            (clock_value, count, duration, last(intervals).dt),
            nothing;
            replay = false,
        )
        shadow = _ec_shadow(
            es,
            work.staged,
            _empic_event_fields(stage),
            species,
            (clock_value, count, duration, last(intervals).dt),
        )
        arrays = vcat(
            _ec_role_arrays((phase.interval, phase.last, phase.next, phase.report)),
            _empic_event_shadow_arrays(shadow),
        )
        _ec_independent(arrays)
        _empic_collection_validate_storage(es, species)
        for array in arrays, ps in species
            _mightalias_particle_storage(array, ps) &&
                throw(ArgumentError("coupled candidate/report aliases public particles"))
        end
        return _ec_publish!(es, work, phase, shadow; particles = true)
    finally
        _ec_release(schedule)
    end
end
function _ec_closed(value, arrays = IdDict{Any,Bool}(); depth = 0)
    # Schema4 has a fixed nonrecursive packet hierarchy;32 exceeds its actual
    # stage/particle/cut nesting while bounding malformed input traversal.
    depth <= 32 || throw(ArgumentError("coupled receipt exceeds its closed schema nesting bound"))
    if value === nothing || value isa Union{Bool,Symbol,String,Int,UInt64,UInt32,BigInt,CIC,TSC}
        return nothing
    elseif value isa _EMPICCoupledQ
        _ec_q(value) || throw(ArgumentError("coupled exact numeric receipt is invalid"))
    elseif value isa Union{Float32,Float64}
        isfinite(value) || throw(ArgumentError("coupled numeric receipt is nonfinite"))
    elseif value isa NamedTuple || value isa Tuple
        foreach(entry -> _ec_closed(entry, arrays; depth = depth + 1), value)
    elseif value isa Array || value isa BitVector
        haskey(arrays, value) &&
            throw(ArgumentError("coupled receipt repeats an array role or contains a cycle"))
        arrays[value] = true
        all(i -> isassigned(value, i), eachindex(value)) ||
            throw(ArgumentError("coupled receipt has unassigned entries"))
        foreach(entry -> _ec_closed(entry, arrays; depth = depth + 1), value)
    else
        throw(ArgumentError("coupled receipt contains an unsupported value type $(typeof(value))"))
    end
    return nothing
end
function _ec_snapshot_check(es, record)
    T = eltype(es.rho_n)
    D = length(es.g.n)
    record isa NamedTuple &&
        keys(record) == (:species, :fields) &&
        record.species isa Vector &&
        length(record.species) == length(es.collection.mobile) ||
        throw(ArgumentError("coupled stage snapshot layout differs"))
    for ps in record.species
        ps isa NamedTuple &&
            keys(ps) == _CHECKPOINT_PARTICLE_FIELDS &&
            ps.x isa NTuple{D,Vector{T}} &&
            ps.v isa NTuple{3,Vector{T}} &&
            ps.weight isa Vector{T} &&
            ps.id isa Vector{UInt64} &&
            ps.tag isa Vector{UInt32} &&
            ps.q isa T &&
            ps.m isa T || throw(
            ArgumentError("coupled stage particle record requires canonical ordinary vectors"),
        )
        count = length(ps.weight)
        all(a -> axes(a, 1) == Base.OneTo(count), (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)) ||
            throw(DimensionMismatch("coupled stage particle sizes differ"))
    end
    _validate_checkpoint_fields(record.fields, _empic_event_fields(es))
    keys(record.fields) == keys(_empic_event_fields(es)) ||
        throw(ArgumentError("coupled stage field inventory has extra fields"))
    for value in values(record.fields), array in (value isa Tuple ? value : (value,))
        array isa Array{T,D} && size(array) == es.g.n && all(isfinite, array) ||
            throw(ArgumentError("coupled stage fields are invalid"))
    end
    populations = _ec_snapshot_populations(record, Val(D), T)
    _empic_collection_validate_population_values(es, populations, es.c, es.relativistic)
    return populations
end
function _ec_install_samples!(stage::EMPIC, ps, sample)
    _ensure_empic_electron_buffers!(stage, nparticles(ps))
    foreach(copyto!, stage.Ep, sample.E)
    foreach(copyto!, stage.Bp, sample.B)
    return nothing
end
function _ec_install_samples!(stage::EMPIC1D, ps, sample)
    _ensure_electron_buffers!(stage, nparticles(ps))
    copyto!(stage.Exp, sample.E[1])
    copyto!(stage.Eyp, sample.E[2])
    copyto!(stage.Bzp, sample.B[3])
    return nothing
end
function _ec_report_check(
    es::Union{EMPIC,EMPIC1D},
    phase::_EMPICCoupledPhase,
    species::AbstractVector,
    fields::NamedTuple,
    clock::Tuple,
    meter::Union{Nothing,_EMPICCoupledMeter};
    replay::Bool = true,
)
    r = phase.report
    T = eltype(es.rho_n)
    D = length(es.g.n)
    r isa NamedTuple && keys(r) == _EC_REPORT_KEYS ||
        throw(ArgumentError("coupled report has noncanonical fields"))
    _ec_closed(r)
    r.kind === (es isa EMPIC ? :vector : :planar) &&
        _ec_identical(r.control, _ec_control(r.control, T)) &&
        _ec_identical(r.control, phase.admission_control) &&
        _ec_identical(r.configuration, _empic_event_configuration(es)) ||
        throw(ArgumentError("coupled historical report policy/source differs"))
    p = r.event_proposal
    p isa NamedTuple &&
        keys(p) == (:cap, :adaptive, :minimum, :limits, :budget) &&
        p.cap isa Vector{T} &&
        p.minimum isa Vector{T} &&
        p.adaptive isa BitVector ||
        throw(ArgumentError("coupled historical calendar proposal differs"))
    _empic_event_proposal(es.collection.mobile, p.cap, p.adaptive, p.minimum, p.limits, p.budget, T)
    r.n_sub isa Int &&
        1 <= r.n_sub <= p.budget &&
        r.requested isa T &&
        r.requested > 0 &&
        r.public_before isa T &&
        r.public_after isa T &&
        r.public_before >= 0 &&
        r.public_after > r.public_before &&
        r.step_before isa Int &&
        r.step_after isa Int &&
        0 <= r.step_before < r.step_after &&
        r.step_after == Base.checked_add(r.step_before, 1) &&
        all(
            _ec_q,
            (
                r.start,
                r.finish,
                r.elapsed_before,
                r.elapsed_after,
                r.initial_energy,
                r.clock_defect,
            ),
        ) &&
        0 <= r.start < r.finish &&
        r.finish - r.start == _exact_rational(r.requested) &&
        0 <= r.elapsed_before < r.elapsed_after ||
        throw(ArgumentError("coupled report clocks differ"))
    time, step, full, fast = clock
    isequal((r.public_after, r.step_after, r.requested), (time, step, full)) &&
        r.finish == phase.time &&
        r.elapsed_after == phase.elapsed &&
        r.initial_energy == phase.initial_energy ||
        throw(ArgumentError("coupled report differs from complete phase"))
    r.initial_calendar isa NamedTuple &&
        keys(r.initial_calendar) == (:interval, :last, :next) &&
        r.initial_calendar.interval isa Vector{T} &&
        r.initial_calendar.last isa Vector{_EMPICCoupledQ} &&
        r.initial_calendar.next isa Vector{_EMPICCoupledQ} &&
        all(a -> axes(a, 1) == Base.OneTo(length(species)), values(r.initial_calendar)) &&
        r.initial_next_h isa T &&
        r.control.h_min <= r.initial_next_h <= r.control.h_max ||
        throw(ArgumentError("coupled report initial calendar/controller differs"))
    for costs in (r.cost_before, r.costs, r.cost_after, r.admission_costs, r.verification_costs)
        _ec_costs_check(costs)
    end
    r.cost_after == phase.costs == _ec_cost_add(r.cost_before, r.costs) ||
        throw(ArgumentError("coupled cumulative work counters differ"))
    for name in _EC_COST_KEYS
        getproperty(r.costs, name) <= getproperty(r.control, Symbol("max_", name)) ||
            throw(ArgumentError("coupled operation exceeded its historical work cap"))
    end
    r.intervals isa Vector{NamedTuple} &&
        !isempty(r.intervals) &&
        r.rejections isa Vector{NamedTuple} &&
        r.decisions isa Vector{NamedTuple} ||
        throw(ArgumentError("coupled trial history containers differ"))
    isequal(phase.next_h, _ec_future_hint(phase, es.collection.schedule.control)) ||
        throw(ArgumentError("coupled future interval hint differs from historical admission"))
    length(r.intervals) + length(r.decisions) <= p.budget &&
        _ec_report_cells(r) == r.costs.report_cells ||
        throw(ArgumentError("coupled report event/cell bound differs"))
    _ec_independent(_ec_role_arrays(r))
    H = copy(r.initial_calendar.interval)
    a = deepcopy(r.initial_calendar.last)
    b = deepcopy(r.initial_calendar.next)
    for i in eachindex(es.collection.mobile)
        if r.step_before == 0 || !es.collection.mobile[i]
            _particle_stored_zero(H[i]) && iszero(a[i]) && iszero(b[i]) ||
                throw(ArgumentError("coupled initial inactive calendar differs"))
        else
            H[i] >= 2floatmin(T) &&
                a[i] <= r.start <= b[i] &&
                b[i] - a[i] == _exact_rational(H[i]) ||
                throw(ArgumentError("coupled initial pending calendar differs"))
        end
    end
    source, elapsed, clock_value = r.start, r.elapsed_before, r.public_before
    hint = r.initial_next_h
    base = 1
    dcursor = 1
    rcursor = 1
    trial_cost = _ec_zero_costs()
    decision_force = 0
    previous_high = nothing
    probe = _ec_new_stage(es)
    reasons = (
        :nonfinite_gather,
        :nonfinite_derivative,
        :nonfinite_velocity,
        :unsupported_velocity,
        :travel,
        :ambiguous_position,
        :nonfinite_current,
        :nonfinite_field,
        :field_constraint,
        :unresolved_continuity,
        :embedded_defect,
    )
    for interval in r.intervals
        interval isa NamedTuple && keys(interval) == _EC_INTERVAL_KEYS ||
            throw(ArgumentError("coupled interval fields differ"))
        before = _ec_snapshot_check(es, interval.origin)
        _ec_snapshot_check(es, interval.low)
        _ec_snapshot_check(es, interval.high)
        previous_high === nothing ||
            isequal(previous_high, interval.origin) ||
            throw(ArgumentError("coupled accepted physical history is discontinuous"))
        interval.source_before == source &&
            interval.elapsed_before == elapsed &&
            isequal(interval.public_before, clock_value) ||
            throw(ArgumentError("coupled interval origins differ"))
        for i in eachindex(es.collection.mobile)
            es.collection.mobile[i] || continue
            bootstrap = r.step_before == 0 && source == r.start && iszero(H[i])
            if bootstrap || b[i] == source
                dcursor <= length(r.decisions) ||
                    throw(ArgumentError("coupled report omitted a due nominal decision"))
                decision = r.decisions[dcursor]
                keys(decision) == (
                        :species,
                        :time,
                        :previous,
                        :bootstrap,
                        :selected,
                        :cap,
                        :upper,
                        :minimum,
                        :adaptive,
                        :limits,
                        :q,
                        :m,
                        :E1,
                        :B1,
                        :V1,
                        :halvings,
                    ) &&
                    decision.species == i &&
                    decision.time == source &&
                    isequal(decision.previous, H[i]) &&
                    decision.bootstrap === bootstrap ||
                    throw(ArgumentError("coupled nominal decision ordering differs"))
                sample = interval.stage0.forces[i]
                sample === nothing && throw(
                    ArgumentError("coupled nominal decision has no actual force-stage sample"),
                )
                _ec_install_samples!(probe, before[i], sample)
                choice = _empic_event_interval_choice(probe, p, before[i], i, H[i], bootstrap)
                expected = merge(
                    (species = i, time = source, previous = H[i], bootstrap = bootstrap),
                    choice,
                )
                _ec_identical(decision, expected) || throw(
                    ArgumentError("coupled nominal decision differs from canonical sampled policy"),
                )
                decision_force = Base.checked_add(decision_force, nparticles(before[i]))
                H[i], a[i], b[i] =
                    decision.selected, source, source + _exact_rational(decision.selected)
                dcursor += 1
            end
        end
        edge, base = _empic_event_next(
            (time = r.start, next = b),
            es.collection.mobile,
            source,
            r.finish,
            base,
            r.n_sub,
        )
        delta = min(_exact_rational(hint), edge - source)
        while rcursor <= length(r.rejections) && r.rejections[rcursor].source == source
            rejected = r.rejections[rcursor]
            keys(rejected) == (
                    :dt,
                    :source,
                    :delta,
                    :clock_value,
                    :reason,
                    :stage,
                    :species,
                    :particle,
                    :defect,
                    :costs,
                ) &&
                rejected.dt isa T &&
                r.control.h_min <= rejected.dt <= r.control.h_max &&
                rejected.delta == delta &&
                isequal(rejected.dt, _empic_event_interval(T, delta)) &&
                isequal(rejected.clock_value, clock_value) &&
                rejected.reason in reasons &&
                rejected.stage isa Symbol &&
                rejected.species isa Int &&
                rejected.particle isa Int &&
                rejected.species >= 0 &&
                rejected.particle >= 0 &&
                (rejected.defect === nothing || _ec_q(rejected.defect)) ||
                throw(ArgumentError("coupled rejected-trial receipt differs"))
            _ec_costs_check(rejected.costs)
            rejected.costs.trials == 1 && rejected.costs.report_cells == 0 ||
                throw(ArgumentError("coupled rejected-trial counters differ"))
            rejected.reason in (:embedded_defect, :field_constraint) &&
                !(rejected.defect !== nothing && rejected.defect > 1) &&
                throw(ArgumentError("coupled rejection has no failed measured defect"))
            rejected.reason === :travel &&
                !(rejected.defect !== nothing && rejected.defect >= 1) &&
                throw(ArgumentError("coupled travel rejection has no failed travel ratio"))
            trial_cost = _ec_cost_add(trial_cost, rejected.costs)
            delta /= 2
            rcursor += 1
        end
        delta >= _exact_rational(r.control.h_min) &&
            interval.source_after == source + delta &&
            interval.dt isa T &&
            isequal(interval.dt, _empic_event_interval(T, delta)) &&
            r.control.h_min <= interval.dt <= r.control.h_max &&
            interval.norm isa _EMPICCoupledQ &&
            0 <= interval.norm <= 1 &&
            interval.norm == maximum(values(interval.norms)) ||
            throw(ArgumentError("coupled accepted trial/halving/defect differs"))
        next_clock, _ = _validated_step_clock(clock_value, r.step_before, interval.dt)
        isequal(next_clock, interval.public_after) &&
            next_clock > clock_value &&
            interval.duration_defect == _exact_rational(interval.dt) - delta &&
            interval.addition_defect ==
            _exact_rational(next_clock) - _exact_rational(clock_value) -
            _exact_rational(interval.dt) &&
            interval.elapsed_after == elapsed + _exact_rational(interval.dt) ||
            throw(ArgumentError("coupled actual interval clock/materialization differs"))
        _ec_costs_check(interval.costs)
        interval.costs.trials == 1 && interval.costs.report_cells == 0 ||
            throw(ArgumentError("coupled accepted-trial work counters differ"))
        if replay
            charges = _pic_collection_charges(es.g, es.n0, before)
            target = _empic_collection_charge_target(es, es.collection, charges)
            fresh = try
                _, packet = _ec_pair!(
                    probe,
                    before,
                    interval.origin.fields,
                    es.collection.mobile,
                    target,
                    interval.dt,
                    delta,
                    clock_value,
                    next_clock,
                    source,
                    elapsed,
                    r.control,
                    meter,
                )
                packet
            catch failure
                failure isa _EMPICCoupledReject || rethrow()
                throw(
                    ArgumentError(
                        "saved coupled accepted trial now fails its own historical admission: $(failure.reason)",
                    ),
                )
            end
            _ec_identical(interval, fresh) || throw(
                ArgumentError(
                    "coupled accepted trial differs from canonical source/force/field/rounding replay",
                ),
            )
        end
        trial_cost = _ec_cost_add(trial_cost, interval.costs)
        source, elapsed, clock_value = source + delta, interval.elapsed_after, next_clock
        hint = _ec_recommend(interval.dt, interval.norm, r.control)
        previous_high = interval.high
    end
    dcursor == length(r.decisions) + 1 &&
        rcursor == length(r.rejections) + 1 &&
        source == r.finish &&
        elapsed == r.elapsed_after &&
        isequal(clock_value, r.public_after) &&
        isequal(last(r.intervals).dt, fast) &&
        H == phase.interval &&
        a == phase.last &&
        b == phase.next ||
        throw(ArgumentError("coupled final clock/calendar/retry history differs"))
    previous_high !== nothing && isequal(previous_high.fields, fields) ||
        throw(ArgumentError("coupled report ending physical fields differ"))
    _ec_species_records_equal(previous_high.species, species) ||
        throw(ArgumentError("coupled report ending particles differ"))
    expected_cost = _ec_cost_add(_ec_cost_add(r.admission_costs, r.verification_costs), trial_cost)
    expected_cost = merge(
        expected_cost,
        (
            force_evals = Base.checked_add(expected_cost.force_evals, decision_force),
            report_cells = r.costs.report_cells,
        ),
    )
    expected_cost == r.costs || throw(
        ArgumentError(
            "coupled operation counters omit actual admission/verification/trial/decision work",
        ),
    )
    summary = _ec_sum_work(r.intervals, r.kind)
    for name in keys(summary)
        isequal(getproperty(r, name), getproperty(summary, name)) ||
            throw(ArgumentError("coupled accepted work aggregate differs"))
    end
    r.defects_before isa NamedTuple &&
        keys(r.defects_before) == (:field, :particle, :coupling) &&
        all(_ec_q, values(r.defects_before)) &&
        r.field_before + r.kinetic_before - r.initial_energy == sum(values(r.defects_before)) &&
        r.total_defect == r.field_defect + r.particle_defect + r.coupling_defect &&
        phase.field_defect == r.defects_before.field + r.field_defect &&
        phase.particle_defect == r.defects_before.particle + r.particle_defect &&
        phase.coupling_defect == r.defects_before.coupling + r.coupling_defect &&
        r.clock_defect ==
        _exact_rational(r.public_after) - _exact_rational(r.public_before) -
        _exact_rational(r.requested) ||
        throw(ArgumentError("coupled cumulative work/clock receipt differs"))
    return nothing
end
function _ec_species_records_equal(records, species)
    length(records) == length(species) || return false
    return all(zip(records, species)) do pair
        record, ps = pair
        isequal((record.q, record.m), (ps.q, ps.m)) && all(
            p -> isequal(p[1], p[2]),
            zip(
                (record.x..., record.v..., record.weight, record.id, record.tag),
                _empic_event_arrays(ps),
            ),
        )
    end
end
