# Fully coupled synchronous ideal kinetic-ion / isothermal-electron hybrid.
# Logical positions, Cartesian velocities, contravariant magnetic flux density.
function _sheared_fields(::Type{T}, n) where {T}
    scalar() = zeros(T, n)
    vector() = ntuple(_ -> scalar(), 3)
    return (n = scalar(), Ji = vector(), ui = vector(), B = vector(), E = vector(),
        J = vector(), pe = scalar(), gradp = vector(), ninv = scalar())
end
_sh_field_arrays(f) = (f.n, f.Ji..., f.ui..., f.B..., f.E..., f.J..., f.pe, f.gradp..., f.ninv)
_sh_particle_arrays(p) = (p.x..., p.v..., p.weight, p.id, p.tag)

struct _ShearedStage{T}
    species::Vector{_ElectrostaticHostParticle{3,T}}
    beta::NTuple{3,Array{T,3}}
end
struct _ShearedParticleRate{T}
    x::NTuple{3,Vector{T}}
    v::NTuple{3,Vector{T}}
end
struct _ShearedRate{T}
    particles::Vector{_ShearedParticleRate{T}}
    beta::NTuple{3,Array{T,3}}
end
struct _ShearedWork{T,F}
    old::_ShearedStage{T}
    trial::_ShearedStage{T}
    endpoint::_ShearedStage{T}
    rates::NTuple{4,_ShearedRate{T}}
    fields::F
    moments::NTuple{4,Array{Rational{BigInt},3}}
    Bcov::NTuple{3,Array{T,3}}
    Ecov::NTuple{3,Array{T,3}}
    curl::NTuple{3,Array{T,3}}
    gradient::NTuple{3,Array{T,3}}
    lapJ::NTuple{3,Array{T,3}}
    hats::NTuple{3,Array{Complex{T},3}}
end

"""
    ShearedPeriodicHybridPIC(mesh, species; Te, B0=(0,0,1),
        div_rtol=sqrt(eps(T)), div_atol=0, refresh_rtol=sqrt(eps(T)))

Real full-f ideal massless-electron hybrid on one periodic nonorthogonal shear.
ParticleSet x is logical and v is Cartesian in this binding. Primary beta is
contravariant flux; fields are Cartesian. init! then coupled RK4 step! own the
complete synchronous phase. This family does not implement a toroidal axis/wall.
"""
mutable struct ShearedPeriodicHybridPIC{T,G,M,F}
    const g::G
    const mesh::M
    const closure::IsothermalElectrons{T}
    const owners::Vector{_ElectrostaticHostParticle{3,T}}
    const beta::NTuple{3,Array{T,3}}
    const fields::F
    const div_rtol::T
    const div_atol::T
    const refresh_rtol::T
    const _counts::Vector{Int}
    const _qm::Vector{NTuple{2,T}}
    const _weights::Vector{Vector{T}}
    const _components::Vector{Tuple}
    const work::_ShearedWork{T,F}
    initialized::Bool
    time::T
    step::Int
    last_dt::T
    last_report::Union{Nothing,ShearedHybridFieldReport{T}}
    busy::Bool
end

function _sh_clone(p::ParticleSet{3,T}) where {T}
    return ParticleSet{3,T}(map(copy, p.x), map(copy, p.v), copy(p.weight),
        copy(p.id), copy(p.tag), p.q, p.m)
end
function _sh_particle(p::ParticleSet{3,T}; values = true) where {T}
    p isa _ElectrostaticHostParticle{3,T} || throw(ArgumentError("sheared ions require ordinary Vector leaves"))
    N = nparticles(p)
    all(a -> length(a) == N, _sh_particle_arrays(p)) || throw(DimensionMismatch("sheared particle component counts differ"))
    isfinite(p.q) && p.q > 0 && isfinite(p.m) && p.m > 0 ||
        throw(ArgumentError("sheared ion q and m must be finite positive stored values"))
    all(x -> isfinite(x) && x >= 0, p.weight) || throw(ArgumentError("sheared full-f weights must be finite nonnegative values"))
    if values
        all(a -> all(isfinite, a), (p.x..., p.v...)) || throw(ArgumentError("sheared phase-space values must be finite"))
    end
    return nothing
end

function _sh_work_arrays(w)
    out = Any[_sh_field_arrays(w.fields)..., w.moments..., w.Bcov..., w.Ecov...,
        w.curl..., w.gradient..., w.lapJ..., w.hats...]
    for stage in (w.old, w.trial, w.endpoint)
        push!(out, stage.species)
        append!(out, stage.beta)
        for p in stage.species
            append!(out, _sh_particle_arrays(p))
        end
    end
    for rate in w.rates
        push!(out, rate.particles)
        append!(out, rate.beta)
        for p in rate.particles
            append!(out, (p.x..., p.v...))
        end
    end
    return out
end
function _sh_storage_arrays(st)
    out = Any[st.owners, st._counts, st._qm, st._weights, st._components,
        st.mesh._s, st.beta..., _sh_field_arrays(st.fields)...]
    for p in st.owners
        append!(out, _sh_particle_arrays(p))
    end
    append!(out, st._weights)
    append!(out, _sh_work_arrays(st.work))
    append!(out, (st.g.cbuf, st.g.tbuf, st.g.abuf, st.g.sbuf,
        st.g.midx..., st.g.kfull..., st.g.kvec..., st.g.ik...))
    return out
end
function _sh_storage(st)
    arrays = _sh_storage_arrays(st)
    for i in eachindex(arrays), j in 1:(i - 1)
        _particle_arrays_mightalias(arrays[i], arrays[j]) &&
            throw(ArgumentError("sheared public/private/readonly/FFT roles must have independent storage"))
    end
    _visit_checkpoint_objects((st.g.plan, st.g.iplan)) do value
        if value isa AbstractArray && !any(a -> a === value, (st.g.cbuf, st.g.tbuf, st.g.abuf, st.g.sbuf))
            any(a -> _particle_arrays_mightalias(value, a), arrays) &&
                throw(ArgumentError("sheared storage overlaps retained FFT plan storage"))
        end
    end
    return nothing
end
function _sh_work_geometry(st)
    N = length(st.owners)
    all(a -> length(a) == N, (st._counts, st._qm, st._weights, st._components)) ||
        throw(DimensionMismatch("sheared fixed inventory geometry changed"))
    all(a -> size(a) == st.g.n, (st.beta..., _sh_field_arrays(st.fields)...,
        _sh_field_arrays(st.work.fields)..., st.work.moments..., st.work.Bcov...,
        st.work.Ecov..., st.work.curl..., st.work.gradient..., st.work.lapJ..., st.work.hats...)) ||
        throw(DimensionMismatch("sheared field/work geometry changed"))
    for stage in (st.work.old, st.work.trial, st.work.endpoint)
        length(stage.species) == N && all(a -> size(a) == st.g.n, stage.beta) ||
            throw(DimensionMismatch("sheared private stage geometry changed"))
        for k in eachindex(stage.species)
            p = stage.species[k]
            _sh_particle(p; values = false)
            nparticles(p) == st._counts[k] && isequal((p.q, p.m), st._qm[k]) &&
                isequal(p.weight, st._weights[k]) || throw(ArgumentError("sheared private stage inventory changed"))
        end
    end
    for rate in st.work.rates
        length(rate.particles) == N && all(a -> size(a) == st.g.n, rate.beta) ||
            throw(DimensionMismatch("sheared private rate geometry changed"))
        for k in eachindex(rate.particles)
            p = rate.particles[k]
            all(a -> length(a) == st._counts[k], (p.x..., p.v...)) ||
                throw(DimensionMismatch("sheared private rate component count changed"))
        end
    end
    return nothing
end
function _sh_bind(st::ShearedPeriodicHybridPIC, species; values = true)
    _boris_accuracy_mode() || throw(ArgumentError("sheared hybrid requires nearest rounding and gradual underflow"))
    species isa AbstractVector || throw(ArgumentError("sheared populations must be an ordered vector"))
    Base.require_one_based_indexing(species)
    st.g === st.mesh.logical || throw(ArgumentError("sheared logical grid binding changed"))
    _shear_geometry(st.mesh)
    length(species) == length(st.owners) || throw(DimensionMismatch("sheared population count changed"))
    _sh_work_geometry(st)
    for k in eachindex(st.owners)
        p = species[k]
        p === st.owners[k] || throw(ArgumentError("sheared population identity/order changed"))
        _sh_particle(p; values)
        nparticles(p) == st._counts[k] && isequal((p.q, p.m), st._qm[k]) &&
            isequal(p.weight, st._weights[k]) || throw(ArgumentError("sheared fixed population inventory changed"))
        all(pair -> pair[1] === pair[2], zip(_sh_particle_arrays(p), st._components[k])) ||
            throw(ArgumentError("sheared public particle leaf identities changed"))
    end
    _sh_storage(st)
    return nothing
end
function _sh_acquire(st)
    st.busy && throw(ArgumentError("sheared state/grid workspace is already in use"))
    st.busy = true
    return nothing
end
function _sh_report(report::ShearedHybridFieldReport{T}) where {T}
    values = (report.divergence, report.tolerance, report.source_norm)
    all(x -> !isnan(x) && x >= 0, values) &&
        report.divergence <= report.tolerance && !report.source_risk &&
        report.range_resolved === all(isfinite, values) ||
        throw(ArgumentError("sheared accepted field receipt is malformed or not Faraday ready"))
    return nothing
end
function _sh_phase(time::T, step::Int, h::T, report) where {T}
    isfinite(time) && time >= 0 && step >= 0 && isfinite(h) && h >= 0 &&
        (step == 0 ? iszero(time) && iszero(h) : time > 0 && 0 < h <= time) ||
        throw(ArgumentError("sheared synchronous clock/completed interval disagree"))
    report isa ShearedHybridFieldReport{T} || throw(ArgumentError("sheared synchronous phase requires its field receipt"))
    _sh_report(report)
    return nothing
end
function _sh_boundary(st)
    st.initialized || throw(ArgumentError("initialize the sheared hybrid before stepping or saving"))
    _sh_phase(st.time, st.step, st.last_dt, st.last_report)
    all(a -> all(isfinite, a), (st.beta..., _sh_field_arrays(st.fields)...)) &&
        all(>(0), st.fields.n) && all(>(0), st.fields.ninv) && all(>=(0), st.fields.pe) ||
        throw(ArgumentError("sheared completed physical fields must be finite with positive actual density/reciprocal"))
    for p in st.owners, d in 1:3, x in p.x[d]
        0 <= x < st.g.L[d] || throw(ArgumentError("sheared completed logical positions must be half-open periodic representatives"))
    end
    return nothing
end

function ShearedPeriodicHybridPIC(mesh::PeriodicShearGrid{T}, species::AbstractVector;
    Te, B0 = (0, 0, 1), div_rtol = sqrt(eps(T)), div_atol = zero(T),
    refresh_rtol = sqrt(eps(T))) where {T<:Union{Float32,Float64}}
    _boris_accuracy_mode() || throw(ArgumentError("sheared hybrid requires nearest rounding and gradual underflow"))
    _shear_geometry(mesh)
    Base.require_one_based_indexing(species)
    isempty(species) && throw(ArgumentError("sheared hybrid requires at least one population"))
    all(p -> p isa _ElectrostaticHostParticle{3,T}, species) || throw(ArgumentError("sheared ions must be ordinary host populations at grid precision"))
    owners = _ElectrostaticHostParticle{3,T}[p for p in species]
    foreach(_sh_particle, owners)
    any(p -> any(>(0), p.weight), owners) || throw(ArgumentError("sheared collection must carry positive full-f charge"))
    temp = _require_finite_nonnegative_real("sheared Te", Te, T)
    (iszero(Te) || !iszero(temp)) || throw(ArgumentError("nonzero sheared Te must survive conversion"))
    r = _require_finite_nonnegative_real("div_rtol", div_rtol, T)
    a = _require_finite_nonnegative_real("div_atol", div_atol, T)
    rr = _require_finite_positive_real("refresh_rtol", refresh_rtol, T)
    r < 1 && (!iszero(r) || !iszero(a)) && rr < 1 || throw(ArgumentError("sheared tolerances require positive allowance and relative values below one"))
    length(B0) == 3 || throw(DimensionMismatch("uniform physical B0 requires three components"))
    initial = ntuple(c -> _require_finite_real("B0", B0[c], T), 3)
    g, counts = mesh.logical, map(nparticles, owners)
    vector() = ntuple(_ -> zeros(T, g.n), 3)
    stage() = _ShearedStage{T}([_sh_clone(p) for p in owners], vector())
    function rate()
        particles = [_ShearedParticleRate{T}(ntuple(_ -> zeros(T, N), 3), ntuple(_ -> zeros(T, N), 3)) for N in counts]
        return _ShearedRate{T}(particles, vector())
    end
    fields = _sheared_fields(T, g.n)
    work = _ShearedWork{T,typeof(fields)}(stage(), stage(), stage(), ntuple(_ -> rate(), 4),
        _sheared_fields(T, g.n), ntuple(_ -> fill(zero(Rational{BigInt}), g.n), 4),
        vector(), vector(), vector(), vector(), vector(), ntuple(_ -> zeros(Complex{T}, g.n), 3))
    beta = vector()
    for I in CartesianIndices(beta[1])
        beta[1][I] = _shear_linear(initial[1], -mesh._s[I[2]], initial[2])
        beta[2][I], beta[3][I] = initial[2], initial[3]
    end
    st = ShearedPeriodicHybridPIC{T,typeof(g),typeof(mesh),typeof(fields)}(g, mesh,
        IsothermalElectrons{T}(temp), owners, beta, fields, r, a, rr, counts,
        [(p.q, p.m) for p in owners], [copy(p.weight) for p in owners],
        Tuple[_sh_particle_arrays(p) for p in owners], work, false, zero(T), 0, zero(T), nothing, false)
    _sh_bind(st, species)
    return st
end

function _sh_copy!(stage, species, beta)
    for k in eachindex(stage.species), c in 1:3
        copyto!(stage.species[k].x[c], species[k].x[c])
        copyto!(stage.species[k].v[c], species[k].v[c])
    end
    for c in 1:3
        copyto!(stage.beta[c], beta[c])
    end
    return nothing
end
function _sh_wrap!(stage, g)
    for p in stage.species, d in 1:3, i in eachindex(p.weight)
        x = p.x[d][i]
        isfinite(x) || throw(ArgumentError("sheared endpoint coordinate must be finite before wrapping"))
        p.x[d][i] = _periodic_coordinate(x, zero(x), g.L[d], g.L[d])
    end
    return nothing
end
function _sh_prepare!(st, stage)
    all(a -> all(isfinite, a), stage.beta) || throw(ArgumentError("sheared primary magnetic flux must be finite"))
    _sheared_moments!(st, stage.species)
    return _sheared_complete_fields!(st, stage.beta)
end

# True continuous Lorentz derivative: no duration, finite-step map, rounded
# q/m or prematurely rounded/cancelled cross product enters this owner.
function _sheared_lorentz_rate(q::T, m::T, E, v, B) where {T}
    all(isfinite, (q, m, E..., v..., B...)) && m > 0 ||
        throw(ArgumentError("sheared Lorentz derivative requires finite original operands and positive mass"))
    Q = _exact_rational
    return ntuple(3) do c
        j, k = mod1(c + 1, 3), mod1(c + 2, 3)
        value = Q(q) * (Q(E[c]) + Q(v[j]) * Q(B[k]) - Q(v[k]) * Q(B[j])) / Q(m)
        _shear_round(T, value, "continuous Lorentz rate"; preserve_nonzero = true)
    end
end
function _sh_rhs!(out, st::ShearedPeriodicHybridPIC{T}, stage) where {T}
    _sh_prepare!(st, stage)
    for c in 1:3
        copyto!(out.beta[c], st.work.curl[c])
    end
    Q = _exact_rational
    for k in eachindex(stage.species)
        p, rate = stage.species[k], out.particles[k]
        for i in eachindex(p.weight)
            E, B = _sheared_gather(st, p, i)
            v = ntuple(c -> p.v[c][i], 3)
            acceleration = _sheared_lorentz_rate(p.q, p.m, E, v, B)
            slope = _shear_slope(st.mesh, p.x[2][i])
            rate.x[1][i] = _shear_round(T, Q(v[1]) - Q(slope) * Q(v[2]), "logical position rate"; preserve_nonzero = true)
            rate.x[2][i], rate.x[3][i] = v[2], v[3]
            for c in 1:3
                rate.v[c][i] = acceleration[c]
            end
        end
    end
    return nothing
end
function _sh_affine!(destination, original, rates, affine)
    for k in eachindex(destination.species), c in 1:3
        for component in (:x, :v)
            dst = getproperty(destination.species[k], component)[c]
            old = getproperty(original.species[k], component)[c]
            sources = map(r -> getproperty(r.particles[k], component)[c], rates)
            for i in eachindex(dst)
                dst[i] = affine(old[i], map(a -> a[i], sources)...)
                isfinite(dst[i]) || throw(ArgumentError("sheared RK particle stage is not representable"))
            end
        end
    end
    for c in 1:3, I in eachindex(destination.beta[c])
        destination.beta[c][I] = affine(original.beta[c][I], map(r -> r.beta[c][I], rates)...)
        isfinite(destination.beta[c][I]) || throw(ArgumentError("sheared RK magnetic stage is not representable"))
    end
    return nothing
end
function _sh_publish!(st, stage, fields, time, count, h, report, provenance = nothing)
    destinations, sources = _sh_field_arrays(st.fields), _sh_field_arrays(fields)
    # Every source/destination role and scalar type is prepared before this point.
    if provenance !== nothing
        for k in eachindex(st.owners), i in eachindex(st.owners[k].weight)
            @inbounds st.owners[k].weight[i] = provenance[k].weight[i]
            @inbounds st.owners[k].id[i] = provenance[k].id[i]
            @inbounds st.owners[k].tag[i] = provenance[k].tag[i]
        end
    end
    for k in eachindex(st.owners), c in 1:3, i in eachindex(st.owners[k].weight)
        @inbounds st.owners[k].x[c][i] = stage.species[k].x[c][i]
        @inbounds st.owners[k].v[c][i] = stage.species[k].v[c][i]
    end
    for c in 1:3, I in eachindex(st.beta[c])
        @inbounds st.beta[c][I] = stage.beta[c][I]
    end
    for k in eachindex(destinations), I in eachindex(destinations[k])
        @inbounds destinations[k][I] = sources[k][I]
    end
    st.initialized = true
    st.time = time
    st.step = count
    st.last_dt = h
    st.last_report = report
    return st
end

function init!(st::ShearedPeriodicHybridPIC{T}, species::AbstractVector = st.owners) where {T}
    _sh_acquire(st)
    try
        _sh_bind(st, species)
        _sh_copy!(st.work.endpoint, species, st.beta)
        _sh_wrap!(st.work.endpoint, st.g)
        report = _sh_prepare!(st, st.work.endpoint)
        return _sh_publish!(st, st.work.endpoint, st.work.fields, zero(T), 0, zero(T), report)
    finally
        st.busy = false
    end
end
function step!(st::ShearedPeriodicHybridPIC{T}, species::AbstractVector, dt::Real) where {T}
    _sh_acquire(st)
    try
        _sh_bind(st, species)
        _sh_boundary(st)
        h = _require_finite_real("sheared dt", dt, T)
        h >= 0 && (iszero(dt) || !iszero(h)) || throw(ArgumentError("sheared dt must be nonnegative and nonzero inputs must survive conversion"))
        iszero(h) && return st
        time, count = _validated_step_clock(st.time, st.step, h)
        time > st.time || throw(ArgumentError("sheared clock must advance at stored precision"))
        w = st.work
        _sh_copy!(w.old, species, st.beta)
        _sh_rhs!(w.rates[1], st, w.old)
        _sh_affine!(w.trial, w.old, (w.rates[1],), _RK4Affine(h, 2))
        _sh_rhs!(w.rates[2], st, w.trial)
        _sh_affine!(w.trial, w.old, (w.rates[2],), _RK4Affine(h, 2))
        _sh_rhs!(w.rates[3], st, w.trial)
        _sh_affine!(w.trial, w.old, (w.rates[3],), _RK4Affine(h, 1))
        _sh_rhs!(w.rates[4], st, w.trial)
        _sh_affine!(w.endpoint, w.old, w.rates, _RK4Affine(h, 6))
        _sh_wrap!(w.endpoint, st.g)
        report = _sh_prepare!(st, w.endpoint)
        return _sh_publish!(st, w.endpoint, w.fields, time, count, h, report)
    finally
        st.busy = false
    end
end

"""Synchronous physical reservoirs and actual work defects; no finite-step conservation theorem."""
function sheared_hybrid_diagnostics(st::ShearedPeriodicHybridPIC{T}, species = st.owners) where {T}
    _sh_acquire(st)
    try
        _sh_bind(st, species)
        _sh_boundary(st)
        Q = _exact_rational
        volume = prod(Q, st.g.dx)
        kinetic = _kinetic_energy_rational(st.owners)
        magnetic = _field_square_sum_exact(st.fields.B) * volume / 2
        electron = _isothermal_exact_energy(st.fields.n, st.closure.Te, st.g.dx)
        total = _isothermal_exact_energy(st.fields.n, st.closure.Te, st.g.dx, kinetic + magnetic)
        number = [_exact_float_sum(p.weight) for p in st.owners]
        charge = [Q(p.q) * n for (p, n) in zip(st.owners, number)]
        momentum = ntuple(c -> sum((Q(p.m) * Q(p.weight[i]) * Q(p.v[c][i])
            for p in st.owners for i in eachindex(p.weight)); init = zero(Rational{BigInt})), 3)
        # Reprepare the same private field to observe work at this physical phase;
        # saved/public fields remain unchanged and are never replaced by this call.
        _sh_copy!(st.work.endpoint, species, st.beta)
        report = _sh_prepare!(st, st.work.endpoint)
        ion_work = zero(Rational{BigInt})
        for p in st.work.endpoint.species, i in eachindex(p.weight)
            E, _ = _sheared_gather(st, p, i)
            ion_work += Q(p.q) * Q(p.weight[i]) * sum(c -> Q(p.v[c][i]) * Q(E[c]), 1:3)
        end
        f = st.work.fields
        grid_work = volume * sum((Q(f.Ji[c][I]) * Q(f.E[c][I]) for c in 1:3 for I in eachindex(f.n)); init = zero(Rational{BigInt}))
        for c in 1:3
            copyto!(st.work.rates[1].beta[c], st.work.curl[c])
        end
        magnetic_power = volume * sum((Q(st.work.Bcov[c][I]) * Q(st.work.rates[1].beta[c][I])
            for c in 1:3 for I in eachindex(f.n)); init = zero(Rational{BigInt}))
        current_work = volume * sum((Q(f.J[c][I]) * Q(f.E[c][I]) for c in 1:3 for I in eachindex(f.n)); init = zero(Rational{BigInt}))
        beta_mean = ntuple(c -> _exact_float_sum(st.beta[c]) / length(st.beta[c]), 3)
        flux = ntuple(c -> beta_mean[c] * prod(d -> d == c ? one(Rational{BigInt}) : Q(st.g.dx[d]) * st.g.n[d], 1:3), 3)
        scalars = map(x -> _from_exact(T, x), (kinetic, magnetic, ion_work - grid_work, magnetic_power + current_work))
        displayed_momentum = map(x -> _from_exact(T, x), momentum)
        return (time = st.time, step = st.step, last_dt = st.last_dt,
            phase = :synchronous_integer, velocity_time_offset = zero(T), number = number,
            charge = charge, momentum = displayed_momentum,
            energy = (kinetic = scalars[1], magnetic = scalars[2], electron_free_energy = electron, total = total),
            ion_work_pairing_defect = scalars[3], magnetic_curl_work_defect = scalars[4],
            minimum_density = minimum(st.fields.n), maximum_shear = maximum(abs, st.mesh._s),
            quadrature_volume_defect = length(st.fields.n) * volume - prod(Q, st.g.L),
            deposited_charge_defect = volume * _exact_float_sum(st.fields.n) - sum(charge),
            mean_beta = beta_mean, mean_cross_section_flux = flux, last_field = st.last_report,
            observed_field = report,
            range_resolved = all(isfinite, (scalars..., displayed_momentum..., electron, total)) && report.range_resolved && st.last_report.range_resolved)
    finally
        st.busy = false
    end
end
