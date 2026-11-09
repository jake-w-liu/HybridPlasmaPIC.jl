struct _TAHStage{T}
    species::Vector{_ElectrostaticHostParticle{3,T}}
    beta::NTuple{3,Array{T,3}}
    pe::Array{T,3}
end
struct _TAHParticleRate{T}
    x::NTuple{3,Vector{T}}
    v::NTuple{3,Vector{T}}
end
struct _TAHRate{T}
    particles::Vector{_TAHParticleRate{T}}
    beta::NTuple{3,Array{T,3}}
    pe::Array{T,3}
end
mutable struct _TAHWork{T,F}
    old::_TAHStage{T}
    trial::_TAHStage{T}
    endpoint::_TAHStage{T}
    rates::NTuple{4,_TAHRate{T}}
    endpoint_rate::_TAHRate{T}
    fields::F
    moments::NTuple{4,Array{_CompositeQ,3}}
    Bcov::NTuple{3,Array{T,3}}
    curlB::NTuple{3,Array{T,3}}
    alpha::NTuple{3,Array{T,3}}
    gradient::NTuple{3,Array{T,3}}
    U::NTuple{3,Array{T,3}}
    curlE::NTuple{3,Array{T,3}}
    zero_current::NTuple{3,Array{T,3}}
    divue::Array{T,3}
    divbeta::Array{T,3}
end

struct ToroidalAnnulusStepReport{T}
    time::T
    step::Int
    dt::T
    before::NTuple{3,_CompositeQ}
    after::NTuple{3,_CompositeQ}
    stages::NTuple{4,ToroidalAnnulusStageReport{T}}
    endpoint::ToroidalAnnulusStageReport{T}
    outward_work::_CompositeQ
    controller_work::_CompositeQ
    stage_defect::_CompositeQ
    rk_defect::_CompositeQ
end

Base.isequal(a::ToroidalAnnulusStepReport, b::ToroidalAnnulusStepReport) =
    typeof(a) === typeof(b) && all(name -> isequal(getfield(a, name), getfield(b, name)), fieldnames(typeof(a)))
Base.hash(report::ToroidalAnnulusStepReport, seed::UInt) =
    hash((typeof(report), ntuple(i -> getfield(report, i), fieldcount(typeof(report)))), seed)

"""
    ToroidalAnnulusHybridPIC(mesh, species; closure, beta0, pe0=0,
        pressure_reservoir=(0,0), div_rtol=sqrt(eps(T)), div_atol=0,
        refresh_rtol=sqrt(eps(T)))

Coupled full-f ions, magnetic flux and adiabatic electron pressure on a
positive-J annulus. Fixed electrical contacts and incoming pressure reservoirs
have explicit measured work. This is not solid-axis or material-wall evolution.
"""
mutable struct ToroidalAnnulusHybridPIC{T,M,C,F}
    const mesh::M
    const closure::C
    const owners::Vector{_ElectrostaticHostParticle{3,T}}
    const beta::NTuple{3,Array{T,3}}
    const fields::F
    const pressure_reservoir::NTuple{2,T}
    const div_rtol::T
    const div_atol::T
    const refresh_rtol::T
    const _counts::Vector{Int}
    const _qm::Vector{NTuple{2,T}}
    const _weights::Vector{Vector{T}}
    const _components::Vector{Tuple}
    const _normal_flux::Tuple{Tuple,Tuple}
    const work::_TAHWork{T,F}
    _committed::_TAHStage{T}
    _configuration::NamedTuple
    initialized::Bool
    time::T
    step::Int
    last_dt::T
    initial_energy::_CompositeQ
    outward_work::_CompositeQ
    controller_work::_CompositeQ
    stage_defect::_CompositeQ
    rk_defect::_CompositeQ
    last_report::Union{Nothing,ToroidalAnnulusStepReport{T}}
    last_endpoint::Union{Nothing,ToroidalAnnulusStageReport{T}}
    _busy::Threads.Atomic{Bool}
end

function _tah_configuration(st)
    m = st.mesh
    return (precision = eltype(st.beta[1]), n = m.n, R0 = m.R0, r_inner = m.r_inner,
        r_outer = m.r_outer, tau = m.tau, L = m.L, dx = m.dx,
        geometry = :circular_annulus_normalized_turns_sigma_minus_one,
        radii = Tuple(m.r), major_coordinates = Tuple(vec(m.R)),
        trigonometry = (Tuple(m.su), Tuple(m.cu), Tuple(m.sv), Tuple(m.cv)),
        sbp_norm = Tuple(m.op.sbp.H),
        transverse_symbols = Tuple(Tuple(work.ky) for work in m.op.transverse),
        metric = :original_tau_radius_major_coordinate_diagonal_v1,
        quadrature = :positive_J_times_stored_SBP_norm,
        coordinates = :logical_positions_cartesian_velocities,
        normalization = :Omega_ci_e1_mu0_1,
        closure = st.closure, shape = :complete_CIC_fold_x_periodic_turns,
        moments = :whole_collection_charge_cartesian_current,
        source_policy = :reject_ohm_source_risk,
        boundary = :zero_tangential_covector_incoming_pressure_SAT_v1,
        pressure_reservoir = st.pressure_reservoir, normal_flux = st._normal_flux,
        sat_coefficient = 1, sat_bound = :h_nu_at_most_one_only,
        time_policy = :collective_RK4_all_sources_endpoint_refresh,
        count = Tuple(st._counts), charges_masses = Tuple(st._qm),
        weights = Tuple(Tuple(w) for w in st._weights), div_rtol = st.div_rtol,
        div_atol = st.div_atol, refresh_rtol = st.refresh_rtol)
end

function _tah_normal_flux(st, beta)
    mesh = st.mesh
    for (face, side) in enumerate((1, mesh.n[1])), k in 1:mesh.n[3], j in 1:mesh.n[2]
        index = j + (k - 1) * mesh.n[2]
        _exact_rational(beta[1][side, j, k]) == _exact_rational(st._normal_flux[face][index]) ||
            throw(ArgumentError("annular normal magnetic boundary data changed"))
    end
    return nothing
end

function _tah_stage_arrays(stage)
    arrays = Any[stage.species, stage.beta..., stage.pe]
    for p in stage.species
        append!(arrays, _tah_particle_arrays(p))
    end
    return arrays
end
function _tah_rate_arrays(rate)
    arrays = Any[rate.particles, rate.beta..., rate.pe]
    for p in rate.particles
        append!(arrays, (p.x..., p.v...))
    end
    return arrays
end
function _tah_owned_arrays(st)
    arrays = Any[st.owners, st.beta..., _tah_field_arrays(st.fields)...,
        st._counts, st._qm, st._weights, st._components]
    append!(arrays, st._weights)
    for p in st.owners
        append!(arrays, _tah_particle_arrays(p))
    end
    for stage in (st._committed, st.work.old, st.work.trial, st.work.endpoint)
        append!(arrays, _tah_stage_arrays(stage))
    end
    for rate in (st.work.rates..., st.work.endpoint_rate)
        append!(arrays, _tah_rate_arrays(rate))
    end
    w = st.work
    append!(arrays, (_tah_field_arrays(w.fields)..., w.moments..., w.Bcov..., w.curlB...,
        w.alpha..., w.gradient..., w.U..., w.curlE..., w.zero_current..., w.divue, w.divbeta))
    append!(arrays, _tah_mesh_arrays(st.mesh))
    return arrays
end

function _tah_storage(st::ToroidalAnnulusHybridPIC{T}, species; values = true) where {T}
    _tah_geometry(st.mesh)
    species isa AbstractVector || throw(ArgumentError("annular species must be ordered"))
    Base.require_one_based_indexing(species)
    length(species) == length(st.owners) == length(st._counts) == length(st._qm) ==
        length(st._weights) == length(st._components) ||
        throw(DimensionMismatch("annular fixed population layout changed"))
    isequal(st._configuration, _tah_configuration(st)) ||
        throw(ArgumentError("annular fixed configuration changed"))
    for (k, p) in enumerate(species)
        _tah_particle_storage(p, T)
        p === st.owners[k] && length(p.weight) == st._counts[k] &&
            length(st._components[k]) == 9 &&
            all(pair -> pair[1] === pair[2], zip(_tah_particle_arrays(p), st._components[k])) ||
            throw(ArgumentError("annular owner/component/count binding changed"))
        if values
            _tah_particle_values(p, st.mesh)
            isequal((p.q, p.m), st._qm[k]) && isequal(p.weight, st._weights[k]) ||
                throw(ArgumentError("annular original particle factors/weights changed"))
        end
    end
    n = st.mesh.n
    for stage in (st._committed, st.work.old, st.work.trial, st.work.endpoint)
        length(stage.species) == length(st.owners) || throw(DimensionMismatch("annular stage population layout changed"))
        all(a -> a isa Array{T,3} && size(a) == n, (stage.beta..., stage.pe)) ||
            throw(DimensionMismatch("annular stage field geometry changed"))
        for (k, p) in enumerate(stage.species)
            _tah_particle_storage(p, T)
            length(p.weight) == st._counts[k] || throw(DimensionMismatch("annular stage count changed"))
        end
    end
    for rate in (st.work.rates..., st.work.endpoint_rate)
        length(rate.particles) == length(st.owners) || throw(DimensionMismatch("annular rate population layout changed"))
        all(a -> a isa Array{T,3} && size(a) == n, (rate.beta..., rate.pe)) ||
            throw(DimensionMismatch("annular rate field geometry changed"))
        for (k, p) in enumerate(rate.particles)
            all(a -> a isa Vector{T} && length(a) == st._counts[k], (p.x..., p.v...)) ||
                throw(DimensionMismatch("annular particle-rate geometry changed"))
        end
    end
    w = st.work
    all(a -> a isa Array{T,3} && size(a) == n,
        (st.beta..., _tah_field_arrays(st.fields)..., _tah_field_arrays(w.fields)...,
            w.Bcov..., w.curlB..., w.alpha..., w.gradient..., w.U..., w.curlE...,
            w.zero_current..., w.divue, w.divbeta)) &&
        all(a -> a isa Array{_CompositeQ,3} && size(a) == n, w.moments) ||
        throw(DimensionMismatch("annular work/field geometry changed"))
    all(a -> all(_particle_stored_zero, a), w.zero_current) ||
        throw(ArgumentError("annular inactive transport scratch changed"))
    _tah_independent(_tah_owned_arrays(st))
    return nothing
end

function _tah_make_stage(species, beta, pe, ::Type{T}) where {T}
    return _TAHStage{T}(_tah_clone_species(species, T), map(copy, beta), copy(pe))
end
function _tah_make_rate(counts, n, ::Type{T}) where {T}
    particles = [_TAHParticleRate{T}(ntuple(_ -> zeros(T, count), 3),
        ntuple(_ -> zeros(T, count), 3)) for count in counts]
    return _TAHRate{T}(particles, ntuple(_ -> zeros(T, n), 3), zeros(T, n))
end

function ToroidalAnnulusHybridPIC(mesh::ToroidalAnnulusGrid{T}, species::AbstractVector;
    closure = EvolvingElectrons{T}(T(5) / T(3)), beta0 = nothing, pe0 = zero(T),
    pressure_reservoir = (zero(T), zero(T)), div_rtol = sqrt(eps(T)), div_atol = zero(T),
    refresh_rtol = sqrt(eps(T))) where {T}
    _tah_geometry(mesh)
    closure isa EvolvingElectrons{T} || throw(ArgumentError("annular hybrid requires grid-precision EvolvingElectrons"))
    closure.γ > one(T) && isfinite(closure.γ) || throw(ArgumentError("annular gamma must exceed one"))
    Base.require_one_based_indexing(species)
    isempty(species) && throw(ArgumentError("annular hybrid requires at least one population"))
    for p in species
        _tah_particle_storage(p, T)
        _tah_particle_values(p, mesh)
    end
    initial_beta = beta0 === nothing ? toroidal_annulus_vacuum_flux(mesh; strength = one(T)) : beta0
    initial_beta isa NTuple{3,Array{T,3}} && all(a -> size(a) == mesh.n && all(isfinite, a), initial_beta) ||
        throw(ArgumentError("annular beta0 requires matching finite ordinary arrays"))
    pressure = if pe0 isa Real
        p = _require_finite_nonnegative_real("annular initial pressure", pe0, T)
        (iszero(pe0) || !iszero(p)) || throw(ArgumentError("nonzero annular pressure disappears on conversion"))
        fill(p, mesh.n)
    else
        pe0 isa Array{T,3} && size(pe0) == mesh.n && all(p -> isfinite(p) && p >= 0, pe0) ||
            throw(ArgumentError("annular pe0 requires matching nonnegative finite pressure"))
        pe0
    end
    pressure_reservoir isa Union{Tuple,AbstractVector} && length(pressure_reservoir) == 2 ||
        throw(DimensionMismatch("annular pressure reservoirs require two values"))
    pressure_reservoir isa AbstractVector && Base.require_one_based_indexing(pressure_reservoir)
    all(x -> x isa Real, pressure_reservoir) || throw(ArgumentError("annular pressure reservoirs must be real scalars"))
    reservoirs = ntuple(2) do side
        value = _require_finite_nonnegative_real("annular pressure reservoir", pressure_reservoir[side], T)
        (iszero(pressure_reservoir[side]) || !iszero(value)) ||
            throw(ArgumentError("nonzero pressure reservoir disappears on conversion"))
        value
    end
    r = _require_finite_nonnegative_real("annular div_rtol", div_rtol, T)
    a = _require_finite_nonnegative_real("annular div_atol", div_atol, T)
    refresh = _require_finite_positive_real("annular refresh_rtol", refresh_rtol, T)
    r < 1 && (r > 0 || a > 0) && refresh < 1 || throw(ArgumentError("annular relative tolerances must be below one with positive allowance"))
    original = Any[initial_beta..., pressure, _tah_mesh_arrays(mesh)...]
    for p in species
        append!(original, _tah_particle_arrays(p))
    end
    _tah_independent(original)
    owned = ToroidalAnnulusGrid(mesh.n; R0 = mesh.R0, r_inner = mesh.r_inner, r_outer = mesh.r_outer, T = T)
    beta = map(copy, initial_beta)
    fields = _tah_fields(T, owned.n)
    copyto!(fields.pe, pressure)
    counts = Int[length(p.weight) for p in species]
    stage() = _tah_make_stage(species, beta, fields.pe, T)
    vector() = ntuple(_ -> zeros(T, owned.n), 3)
    work = _TAHWork{T,typeof(fields)}(stage(), stage(), stage(),
        ntuple(_ -> _tah_make_rate(counts, owned.n, T), 4), _tah_make_rate(counts, owned.n, T),
        _tah_fields(T, owned.n), ntuple(_ -> fill(zero(_CompositeQ), owned.n), 4),
        vector(), vector(), vector(), vector(), vector(), vector(), vector(),
        zeros(T, owned.n), zeros(T, owned.n))
    normal = ntuple(side -> Tuple(vec(copy(view(beta[1], side == 1 ? 1 : owned.n[1], :, :)))), 2)
    z = zero(_CompositeQ)
    st = ToroidalAnnulusHybridPIC{T,typeof(owned),typeof(closure),typeof(fields)}(
        owned, closure, _ElectrostaticHostParticle{3,T}[p for p in species], beta, fields,
        reservoirs, r, a, refresh, counts, [(p.q, p.m) for p in species],
        [copy(p.weight) for p in species], Tuple[_tah_particle_arrays(p) for p in species],
        normal, work, stage(), NamedTuple(), false, zero(T), 0, zero(T), z, z, z, z, z,
        nothing, nothing, Threads.Atomic{Bool}(false))
    st._configuration = _tah_configuration(st)
    _tah_storage(st, species)
    return st
end

function _tah_acquire(st)
    Threads.atomic_cas!(st._busy, false, true) == false || throw(ArgumentError("annular hybrid already has an active operation"))
    return nothing
end
_tah_release(st) = (Threads.atomic_xchg!(st._busy, false); nothing)

function _tah_copy!(destination, species, beta, pe)
    _pic_prepare_species_stages!(destination.species, species)
    foreach(copyto!, destination.beta, beta)
    copyto!(destination.pe, pe)
    return nothing
end

function _tah_rhs!(rate, st::ToroidalAnnulusHybridPIC{T}, stage, h::T) where {T}
    data = _tah_prepare_fields!(st, stage, rate, h)
    kinetic = particle = zero(_CompositeQ)
    Q = _exact_rational
    for (p, output) in zip(stage.species, rate.particles), marker in eachindex(p.weight)
        E, B = _tah_gather(st, p, marker)
        velocity = ntuple(c -> p.v[c][marker], 3)
        acceleration = _sheared_lorentz_rate(p.q, p.m, E, velocity, B)
        geo = _tah_point_geometry(st.mesh, ntuple(c -> p.x[c][marker], 3))
        for c in 1:3
            output.x[c][marker] = _tah_round(T, _tah_dot(geo.basis[c], map(Q, velocity)) / geo.h[c],
                "annular logical position rate"; preserve_nonzero = true)
            output.v[c][marker] = acceleration[c]
        end
        kinetic += Q(p.m) * Q(p.weight[marker]) * _tah_dot(map(Q, velocity), map(Q, acceleration))
        particle += Q(p.q) * Q(p.weight[marker]) * _tah_dot(map(Q, velocity), map(Q, E))
    end
    return _tah_stage_report(T, h, data, kinetic, particle)
end

function _tah_affine!(destination, original, rates, affine)
    _tah_copy!(destination, original.species, original.beta, original.pe)
    for k in eachindex(destination.species), c in 1:3, name in (:x, :v)
        output = getproperty(destination.species[k], name)[c]
        initial = getproperty(original.species[k], name)[c]
        for marker in eachindex(output)
            output[marker] = affine(initial[marker], map(r -> getproperty(r.particles[k], name)[c][marker], rates)...)
            isfinite(output[marker]) || throw(ArgumentError("annular RK particle stage is not representable"))
        end
    end
    for c in 1:3, I in eachindex(destination.beta[c])
        destination.beta[c][I] = affine(original.beta[c][I], map(r -> r.beta[c][I], rates)...)
        isfinite(destination.beta[c][I]) || throw(ArgumentError("annular RK magnetic stage is not representable"))
    end
    for I in eachindex(destination.pe)
        destination.pe[I] = affine(original.pe[I], map(r -> r.pe[I], rates)...)
        isfinite(destination.pe[I]) && destination.pe[I] >= 0 || throw(ArgumentError("annular RK pressure must remain nonnegative"))
    end
    return nothing
end

function _tah_wrap!(stage)
    for p in stage.species, d in 2:3, marker in eachindex(p.weight)
        p.x[d][marker] = _periodic_coordinate(p.x[d][marker], zero(eltype(p.weight)), one(eltype(p.weight)), one(eltype(p.weight)))
    end
    return nothing
end

function _tah_reservoirs(st, species, beta, pe)
    Q = _exact_rational
    magnetic = thermal = zero(_CompositeQ)
    gamma = Q(st.closure.γ)
    for I in CartesianIndices(pe)
        metric, J = _tah_metric(st.mesh, I)
        volume = _tah_volume(st.mesh, I)
        magnetic += volume * sum(metric[c]^2 * Q(beta[c][I])^2 / J^2 for c in 1:3) / 2
        thermal += volume * Q(pe[I]) / (gamma - 1)
    end
    return (_kinetic_energy_rational(species), magnetic, thermal)
end

function _tah_report_valid(report::ToroidalAnnulusStageReport{T}) where {T}
    all(_sce_q_finite, (report.divergence, report.divergence_tolerance, report.nu...,
        report.kinetic_power, report.particle_work, report.grid_work, report.current_work,
        report.magnetic_power, report.thermal_power, report.pressure_work, report.poynting...,
        report.enthalpy..., report.electric_controller..., report.thermal_controller...,
        report.lorentz_defect, report.shape_defect, report.ohm_defect, report.magnetic_defect,
        report.thermal_defect, report.balance_defect, report.logical_curl_defect,
        report.pressure_assembly_defect)) && isfinite(report.dt) && report.dt >= 0 &&
        0 <= report.divergence <= report.divergence_tolerance && all(>=(0), report.nu) &&
        all(f -> _exact_rational(report.dt) * f <= 1, report.nu) &&
        report.lorentz_defect == report.kinetic_power - report.particle_work &&
        report.shape_defect == report.particle_work - report.grid_work &&
        report.ohm_defect == report.grid_work - report.current_work + report.pressure_work - sum(report.electric_controller) &&
        report.magnetic_defect == report.magnetic_power + report.current_work + sum(report.poynting) &&
        report.thermal_defect == report.thermal_power - sum(report.thermal_controller) - report.pressure_work + sum(report.enthalpy) &&
        report.balance_defect == report.lorentz_defect + report.shape_defect + report.ohm_defect + report.magnetic_defect + report.thermal_defect ||
        throw(ArgumentError("annular stage receipt has inconsistent source/work/admission values"))
    return nothing
end

function _tah_integral(stages, h, value)
    Q = _exact_rational
    return Q(h) * sum((weight * value(stage) for (weight, stage) in zip((1, 2, 2, 1), stages)); init = zero(_CompositeQ)) / 6
end

function _tah_step_report(time, count, h, before, after, stages, endpoint)
    outward = _tah_integral(stages, h, r -> sum(r.poynting) + sum(r.enthalpy))
    controller = _tah_integral(stages, h, r -> sum(r.electric_controller) + sum(r.thermal_controller))
    defect = _tah_integral(stages, h, r -> r.balance_defect)
    rk = sum(after) - sum(before) - _tah_integral(stages, h, r -> r.kinetic_power + r.magnetic_power + r.thermal_power)
    return ToroidalAnnulusStepReport(time, count, h, before, after, stages, endpoint, outward, controller, defect, rk)
end

function _tah_step_valid(report::ToroidalAnnulusStepReport{T}, time, count, h) where {T}
    isequal(report.time, time) && report.step == count && isequal(report.dt, h) &&
        isfinite(h) && h > 0 && all(x -> _sce_q_finite(x) && x >= 0, (report.before..., report.after...)) ||
        throw(ArgumentError("annular physical step report has invalid phase/reservoirs"))
    foreach(_tah_report_valid, (report.stages..., report.endpoint))
    all(r -> isequal(r.dt, h), (report.stages..., report.endpoint)) || throw(ArgumentError("annular stage durations differ"))
    expected = _tah_step_report(time, count, h, report.before, report.after, report.stages, report.endpoint)
    all(name -> getfield(report, name) == getfield(expected, name),
        (:outward_work, :controller_work, :stage_defect, :rk_defect)) &&
        sum(report.after) - sum(report.before) + report.outward_work - report.controller_work == report.stage_defect + report.rk_defect ||
        throw(ArgumentError("annular complete step energy ledger is inconsistent"))
    return nothing
end

function _tah_phase(st, reservoirs, initialized, time, count, h, initial, outward, controller, defect, rk, report, endpoint)
    T = eltype(st.beta[1])
    reservoirs isa NTuple{3,_CompositeQ} && all(x -> _sce_q_finite(x) && x >= 0, reservoirs) ||
        throw(ArgumentError("annular current physical reservoirs are invalid"))
    energy = sum(reservoirs)
    initialized === true && time isa T && count isa Int && h isa T && count >= 0 &&
        isfinite(time) && time >= 0 && isfinite(h) && h >= 0 &&
        all(_sce_q_finite, (initial, outward, controller, defect, rk)) && initial >= 0 &&
        energy - initial + outward - controller == defect + rk &&
        endpoint isa ToroidalAnnulusStageReport{T} || throw(ArgumentError("annular synchronous phase/accounting is invalid"))
    _tah_report_valid(endpoint)
    isequal(endpoint.dt, h) || throw(ArgumentError("annular endpoint duration differs from completed interval"))
    if count == 0
        _particle_stored_zero(time) && _particle_stored_zero(h) && report === nothing &&
            all(iszero, (outward, controller, defect, rk)) || throw(ArgumentError("unadvanced annular phase carries completed-step history"))
    else
        time > 0 && 0 < h <= time && report isa ToroidalAnnulusStepReport{T} ||
            throw(ArgumentError("advanced annular phase lacks a real positive completed interval"))
        _tah_step_valid(report, time, count, h)
        report.after == reservoirs && isequal(report.endpoint, endpoint) ||
            throw(ArgumentError("annular preceding report differs from its current primary energy/endpoint"))
        if count == 1
            isequal(time, h) && sum(report.before) == initial &&
                outward == report.outward_work && controller == report.controller_work &&
                defect == report.stage_defect && rk == report.rk_defect ||
                throw(ArgumentError("annular first completed interval/history is inconsistent"))
        end
    end
    return nothing
end

function _tah_same_primary(st, species)
    current = st._committed
    all(pair -> isequal(pair[1], pair[2]), zip(st.beta, current.beta)) &&
        isequal(st.fields.pe, current.pe) || throw(ArgumentError("annular primary fields changed; call init!"))
    for (p, saved) in zip(species, current.species)
        isequal((p.q, p.m), (saved.q, saved.m)) &&
            all(pair -> isequal(pair[1], pair[2]), zip(_tah_particle_arrays(p), _tah_particle_arrays(saved))) ||
            throw(ArgumentError("annular physical particle values/provenance changed; call init!"))
        all(d -> all(x -> 0 <= x < 1, p.x[d]), 2:3) || throw(ArgumentError("annular completed turns must be half-open representatives"))
    end
    return nothing
end

function _tah_refresh(saved, prepared, tolerance)
    for name in keys(saved)
        left, right = getproperty(saved, name), getproperty(prepared, name)
        if name in (:n, :Ji, :pe)
            isequal(left, right) || throw(ArgumentError("annular source/pressure field differs from actual complete populations"))
        else
            _sh_relative_block(left, right, tolerance)
        end
    end
    return nothing
end

function _tah_ready(st, species)
    _tah_storage(st, species)
    st.initialized || throw(ArgumentError("initialize the annular hybrid before evolution or save"))
    _tah_same_primary(st, species)
    _tah_normal_flux(st, st.beta)
    reservoirs = _tah_reservoirs(st, species, st.beta, st.fields.pe)
    _tah_phase(st, reservoirs, true, st.time, st.step, st.last_dt, st.initial_energy,
        st.outward_work, st.controller_work, st.stage_defect, st.rk_defect, st.last_report, st.last_endpoint)
    observed = _tah_rhs!(st.work.endpoint_rate, st, st._committed, st.last_dt)
    _tah_refresh(st.fields, st.work.fields, st.refresh_rtol)
    _tah_contact_fields(st, st.fields)
    _tah_endpoint_receipt(st.last_endpoint, observed, st.fields, st)
    return observed
end

function _tah_publish!(st, stage, fields, time, count, h, initial, outward, controller, defect, rk, report, endpoint;
    provenance = false)
    # All graph/range/source/boundary/ledger checks and allocations precede these stores.
    for (destination, prepared) in zip(st.owners, stage.species)
        foreach(copyto!, destination.x, prepared.x)
        foreach(copyto!, destination.v, prepared.v)
        if provenance
            copyto!(destination.weight, prepared.weight)
            copyto!(destination.id, prepared.id)
            copyto!(destination.tag, prepared.tag)
            destination.q, destination.m = prepared.q, prepared.m
        end
    end
    foreach(copyto!, st.beta, stage.beta)
    foreach(copyto!, _tah_field_arrays(st.fields), _tah_field_arrays(fields))
    old = st._committed
    st._committed = stage
    st.work.endpoint = old
    st.time, st.step, st.last_dt = time, count, h
    st.initial_energy, st.outward_work, st.controller_work = initial, outward, controller
    st.stage_defect, st.rk_defect = defect, rk
    st.last_report, st.last_endpoint = report, endpoint
    st.initialized = true
    return st
end

function init!(st::ToroidalAnnulusHybridPIC{T}, species::AbstractVector = st.owners) where {T}
    _tah_acquire(st)
    try
        _tah_storage(st, species)
        _tah_copy!(st.work.endpoint, species, st.beta, st.fields.pe)
        _tah_wrap!(st.work.endpoint)
        endpoint = _tah_rhs!(st.work.endpoint_rate, st, st.work.endpoint, zero(T))
        reservoirs = _tah_reservoirs(st, st.work.endpoint.species, st.work.endpoint.beta, st.work.endpoint.pe)
        initial, z = sum(reservoirs), zero(_CompositeQ)
        _tah_phase(st, reservoirs, true, zero(T), 0, zero(T), initial, z, z, z, z, nothing, endpoint)
        _tah_storage(st, species)
        return _tah_publish!(st, st.work.endpoint, st.work.fields, zero(T), 0, zero(T), initial, z, z, z, z, nothing, endpoint)
    finally
        _tah_release(st)
    end
end

function step!(st::ToroidalAnnulusHybridPIC{T}, species::AbstractVector, dt::Real) where {T}
    _tah_acquire(st)
    try
        _tah_ready(st, species)
        h = _require_finite_real("annular dt", dt, T)
        h >= 0 && (iszero(dt) || !iszero(h)) || throw(ArgumentError("annular dt must be nonnegative and nonzero inputs must survive conversion"))
        iszero(h) && return st
        time, count = _validated_step_clock(st.time, st.step, h)
        time > st.time || throw(ArgumentError("annular global clock must advance"))
        w = st.work
        before = _tah_reservoirs(st, species, st.beta, st.fields.pe)
        _tah_copy!(w.old, species, st.beta, st.fields.pe)
        p1 = _tah_rhs!(w.rates[1], st, w.old, h)
        _tah_affine!(w.trial, w.old, (w.rates[1],), _RK4Affine(h, 2))
        p2 = _tah_rhs!(w.rates[2], st, w.trial, h)
        _tah_affine!(w.trial, w.old, (w.rates[2],), _RK4Affine(h, 2))
        p3 = _tah_rhs!(w.rates[3], st, w.trial, h)
        _tah_affine!(w.trial, w.old, (w.rates[3],), _RK4Affine(h, 1))
        p4 = _tah_rhs!(w.rates[4], st, w.trial, h)
        _tah_affine!(w.endpoint, w.old, w.rates, _RK4Affine(h, 6))
        _tah_wrap!(w.endpoint)
        endpoint = _tah_rhs!(w.endpoint_rate, st, w.endpoint, h)
        after = _tah_reservoirs(st, w.endpoint.species, w.endpoint.beta, w.endpoint.pe)
        report = _tah_step_report(time, count, h, before, after, (p1, p2, p3, p4), endpoint)
        outward, controller = st.outward_work + report.outward_work, st.controller_work + report.controller_work
        defect, rk = st.stage_defect + report.stage_defect, st.rk_defect + report.rk_defect
        _tah_phase(st, after, true, time, count, h, st.initial_energy, outward, controller, defect, rk, report, endpoint)
        _tah_storage(st, species)
        return _tah_publish!(st, w.endpoint, w.fields, time, count, h, st.initial_energy, outward, controller, defect, rk, report, endpoint)
    finally
        _tah_release(st)
    end
end
