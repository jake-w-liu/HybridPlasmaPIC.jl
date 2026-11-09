# One composite Neumann solve and synchronous original-force Verlet map. All
# solve/report candidates are private to a call; only particle stage leaves persist.
struct StaticCompositeESReport{T<:Union{Float32,Float64}}
    dt::T
    source_work::_CompositeQ
    field_change::_CompositeQ
    kinetic_change::_CompositeQ
    field_solve_defect::_CompositeQ
    particle_work_defect::_CompositeQ
    total_energy_defect::_CompositeQ
    max_gauss_residual::_CompositeQ
    max_constitutive_residual::_CompositeQ
    interface_crossings::Vector{Int}
end

struct _CompositeESBinding{T}
    owners::Vector{_ElectrostaticHostParticle{1,T}}
    components::Vector{Tuple}
    candidates::Vector{_ElectrostaticHostParticle{1,T}}
    counts::Vector{Int}
    bins::_CompositeParticleBins

    function _CompositeESBinding{T}(
        owners::Vector{_ElectrostaticHostParticle{1,T}},
        components::Vector{Tuple},
        candidates::Vector{_ElectrostaticHostParticle{1,T}},
        counts::Vector{Int},
        bins::_CompositeParticleBins,
    ) where {T}
        return new{T}(owners, components, candidates, counts, bins)
    end
end

struct StaticCompositeElectrostaticPIC{
    T<:Union{Float32,Float64},
    M<:StaticAMRHierarchy1D{T},
    S<:YeeTriangularCloud,
    B,
    P,
    F<:Tuple,
}
    mesh::M
    shape::S
    background_charge::B
    left_plate_charge::P
    E::Vector{T}
    phi::Vector{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    initialized::Base.RefValue{Bool}
    last_dt::Base.RefValue{T}
    right_plate_charge::Base.RefValue{_CompositeQ}
    last_impulse::Vector{_CompositeQ}
    last_report::Base.RefValue{Union{Nothing,StaticCompositeESReport{T}}}
    _factors::F
    _binding::Base.RefValue{_CompositeESBinding{T}}
end

function _sce_mode()
    rounding(Float32) === RoundNearest &&
        rounding(Float64) === RoundNearest &&
        !get_zero_subnormals() ||
        throw(ArgumentError("static composite ES requires nearest rounding and gradual underflow"))
    return nothing
end
_sce_particle_arrays(ps) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
_sce_q_finite(x) = x isa _CompositeQ && denominator(x) > 0
_sce_ieee_word(x::Float32) = reinterpret(UInt32, x)
_sce_ieee_word(x::Float64) = reinterpret(UInt64, x)
_sce_ieee_config(x) = (typeof(x), _sce_ieee_word(x))

function _sce_population_storage(species, ::Type{T}) where {T}
    species isa AbstractVector ||
        throw(ArgumentError("composite species must be an ordered vector"))
    Base.require_one_based_indexing(species)
    for ps in species
        ps isa _ElectrostaticHostParticle{1,T} || throw(
            ArgumentError("composite ES requires 1D grid-precision ordinary Vector particles"),
        )
        _check_particle_storage_axes(ps)
        _require_particle_ordinary_stores(ps, "static composite ES")
    end
    return nothing
end

function _sce_values(mesh, species, factors)
    length(species) == length(factors) ||
        throw(ArgumentError("composite ES species count is constructor-bound"))
    for (k, ps) in enumerate(species)
        q, m = _validated_charge_mass(ps)
        isequal((q, m), factors[k]) || throw(
            ArgumentError("composite ES original q/m are constructor-bound; construct a new model"),
        )
        for p in eachindex(ps.weight)
            w = ps.weight[p]
            isfinite(w) && _stored_nonnegative(w) ||
                throw(ArgumentError("composite particle weights must be finite and nonnegative"))
            all(c -> isfinite(ps.v[c][p]), 1:3) && isfinite(ps.x[1][p]) ||
                throw(ArgumentError("composite particle positions and velocities must be finite"))
            _particle_stored_zero(ps.v[2][p]) && _particle_stored_zero(ps.v[3][p]) ||
                throw(ArgumentError("composite ES requires stored-zero transverse velocities"))
            _static_amr_leaf(mesh, _exact_rational(ps.x[1][p]))
        end
    end
    return nothing
end

function _sce_make_binding(species, sources, mesh, ::Type{T}; components = species) where {T}
    owners = _ElectrostaticHostParticle{1,T}[]
    for ps in species
        ps isa _ElectrostaticHostParticle{1,T} ||
            throw(ArgumentError("composite ES requires ordinary 1D grid-precision Vector particles"))
        push!(owners, ps)
    end
    arrays = Tuple[_sce_particle_arrays(ps) for ps in components]
    candidates = _ElectrostaticHostParticle{1,T}[]
    _pic_prepare_species_stages!(candidates, sources)
    return _CompositeESBinding{T}(
        owners,
        arrays,
        candidates,
        Int[length(ps.weight) for ps in sources],
        _sce_bins(sources, mesh),
    )
end

function _sce_arrays(st, species)
    writable = Any[st.E, st.phi, st.last_impulse]
    for ps in species
        append!(writable, _sce_particle_arrays(ps))
    end
    for ps in st._binding[].candidates
        _check_particle_storage_axes(ps)
        append!(writable, _sce_particle_arrays(ps))
    end
    readonly = _static_amr_arrays(st.mesh)
    binding = st._binding[]
    append!(
        readonly,
        (binding.counts, binding.bins.offsets, binding.bins.species, binding.bins.particle),
    )
    st.last_report[] isa StaticCompositeESReport &&
        push!(readonly, st.last_report[].interface_crossings)
    return writable, readonly
end

function _sce_storage(st::StaticCompositeElectrostaticPIC{T}, species) where {T}
    _static_amr_validate(st.mesh)
    _sce_population_storage(species, T)
    length(species) == length(st._factors) ||
        throw(ArgumentError("composite ES species count is constructor-bound"))
    n = length(st.mesh.widths)
    length(st.E) == n && length(st.phi) == n + 1 && length(st.last_impulse) == n ||
        throw(DimensionMismatch("composite physical field/history sizes changed"))
    binding = st._binding[]
    length(binding.owners) ==
    length(binding.components) ==
    length(binding.candidates) ==
    length(binding.counts) ==
    length(st._factors) || throw(DimensionMismatch("composite binding storage changed"))
    length(binding.bins.offsets) == n + 1 &&
        length(binding.bins.species) == length(binding.bins.particle) ||
        throw(DimensionMismatch("composite ownership-bin storage changed"))
    _sce_population_storage(binding.candidates, T)
    all(a -> length(a) == 7, binding.components) ||
        throw(DimensionMismatch("composite component binding storage changed"))
    refs = (
        st.time,
        st.step,
        st.initialized,
        st.last_dt,
        st.right_plate_charge,
        st.last_report,
        st._binding,
    )
    for i in eachindex(refs), j = 1:(i-1)
        refs[i] === refs[j] &&
            throw(ArgumentError("composite scalar state roles must be independent"))
    end
    writable, readonly = _sce_arrays(st, species)
    for i in eachindex(writable)
        for j = 1:(i-1)
            _particle_arrays_mightalias(writable[i], writable[j]) && throw(
                ArgumentError("composite field/particle/history/stage storage must be independent"),
            )
        end
        for a in readonly
            _particle_arrays_mightalias(writable[i], a) && throw(
                ArgumentError(
                    "composite writable storage must not alias geometry or ownership metadata",
                ),
            )
        end
    end
    return nothing
end

function _sce_binding(st, species; values = true, counts = true)
    binding = st._binding[]
    length(species) == length(binding.owners) ||
        throw(ArgumentError("composite population membership changed; call init!"))
    for (k, ps) in enumerate(species)
        ps === binding.owners[k] ||
            throw(ArgumentError("composite population order/identity changed; call init!"))
        all(pair -> pair[1] === pair[2], zip(_sce_particle_arrays(ps), binding.components[k])) ||
            throw(ArgumentError("composite particle component identity changed; call init!"))
        counts &&
            length(ps.weight) != binding.counts[k] &&
            throw(ArgumentError("composite particle count changed; call init!"))
        if values
            frozen = binding.candidates[k]
            isequal((ps.q, ps.m), (frozen.q, frozen.m)) &&
                isequal(ps.weight, frozen.weight) &&
                isequal(ps.id, frozen.id) &&
                isequal(ps.tag, frozen.tag) ||
                throw(ArgumentError("composite particle sources/provenance changed; call init!"))
        end
    end
    return nothing
end

function _sce_configuration(st::StaticCompositeElectrostaticPIC{T}) where {T}
    all(pair -> all(isfinite, pair) && _stored_positive(pair[2]), st._factors) ||
        throw(ArgumentError("invalid composite constructor factor configuration"))
    isfinite(st.background_charge) &&
        isfinite(st.left_plate_charge) &&
        isfinite(st.shape.half_width) &&
        _stored_positive(st.shape.half_width) ||
        throw(ArgumentError("invalid composite scalar source configuration"))
    return (
        T = T,
        nbase = st.mesh.nbase,
        length = _sce_ieee_config(st.mesh.L),
        refinements = st.mesh.refinements,
        ratio = 2,
        ghosts = 0,
        shape = (:triangular_physical_cloud, _sce_ieee_config(st.shape.half_width)),
        background = _sce_ieee_config(st.background_charge),
        left_plate = _sce_ieee_config(st.left_plate_charge),
        factors = map(pair -> map(_sce_ieee_config, pair), st._factors),
        field_policy = :original_source_composite_neumann_rounded_v1,
        force_policy = :original_geometry_boris_half_v1,
        phase = :synchronous_velocity_verlet,
        timestep = :fixed_after_first_success,
        event_policy = :strict_event_free,
        arithmetic = :nearest_gradual,
    )
end

function _sce_solution(mesh, shape, background, left, right, species, ::Type{T}) where {T}
    n = length(mesh.widths)
    rhs = fill(zero(_CompositeQ), n + 1)
    _sce_source_rhs!(rhs, species, mesh, shape, background, left, right)
    E, phi = Vector{T}(undef, n), Vector{T}(undef, n + 1)
    phi[1] = zero(T)
    field, potential = zero(_CompositeQ), zero(_CompositeQ)
    for cell = 1:n
        field += rhs[cell]
        potential -= mesh.widths[cell] * field
        E[cell], phi[cell+1] = _from_exact(T, field), _from_exact(T, potential)
        isfinite(E[cell]) || throw(
            ArgumentError("composite electric field materialization is nonfinite at leaf $cell"),
        )
        isfinite(phi[cell+1]) || throw(
            ArgumentError("composite potential materialization is nonfinite at node $(cell+1)"),
        )
    end
    return (; E, phi, rhs)
end

function _sce_residuals(mesh, E, phi, rhs)
    n = length(E)
    gauss = Vector{_CompositeQ}(undef, n + 1)
    constitutive = Vector{_CompositeQ}(undef, n)
    gauss[1] = _exact_rational(E[1]) - rhs[1]
    for i = 2:n
        gauss[i] = _exact_rational(E[i]) - _exact_rational(E[i-1]) - rhs[i]
    end
    gauss[end] = -_exact_rational(E[end]) - rhs[end]
    for cell = 1:n
        constitutive[cell] =
            _exact_rational(E[cell]) +
            (_exact_rational(phi[cell+1]) - _exact_rational(phi[cell])) / mesh.widths[cell]
    end
    return (;
        gauss,
        constitutive,
        max_gauss = maximum(abs, gauss),
        max_constitutive = maximum(abs, constitutive),
    )
end

function _sce_report_values(st::StaticCompositeElectrostaticPIC{T}, report, dt, total) where {T}
    report isa StaticCompositeESReport{T} && isequal(report.dt, dt) ||
        throw(ArgumentError("invalid composite step report or duration"))
    for name in (
        :source_work,
        :field_change,
        :kinetic_change,
        :field_solve_defect,
        :particle_work_defect,
        :total_energy_defect,
        :max_gauss_residual,
        :max_constitutive_residual,
    )
        _sce_q_finite(getfield(report, name)) ||
            throw(ArgumentError("composite report $name must be a finite exact scalar"))
    end
    report.field_change + report.source_work == report.field_solve_defect &&
        report.kinetic_change - report.source_work == report.particle_work_defect &&
        report.field_change + report.kinetic_change == report.total_energy_defect &&
        report.field_solve_defect + report.particle_work_defect == report.total_energy_defect ||
        throw(ArgumentError("composite signed energy report identities are inconsistent"))
    report.max_gauss_residual >= 0 && report.max_constitutive_residual >= 0 ||
        throw(ArgumentError("composite residual magnitudes must be nonnegative"))
    length(report.interface_crossings) == length(st.mesh.interfaces) &&
        all(x -> 0 <= x <= total, report.interface_crossings) ||
        throw(ArgumentError("invalid composite interface crossing counts"))
    return nothing
end

function _sce_phase_values(
    st::StaticCompositeElectrostaticPIC{T},
    initialized,
    time,
    step,
    dt,
    impulse,
    report,
    total,
) where {T}
    initialized isa Bool &&
        time isa T &&
        step isa Int &&
        step >= 0 &&
        dt isa T &&
        isfinite(time) &&
        _stored_nonnegative(time) &&
        isfinite(dt) &&
        _stored_nonnegative(dt) || throw(ArgumentError("invalid composite phase or clock"))
    impulse isa Vector{_CompositeQ} &&
        length(impulse) == length(st.mesh.widths) &&
        all(_sce_q_finite, impulse) || throw(
        ArgumentError("composite integrated-current history must contain finite exact scalars"),
    )
    if step == 0
        _particle_stored_zero(time) &&
            _particle_stored_zero(dt) &&
            all(iszero, impulse) &&
            report === nothing ||
            throw(ArgumentError("unadvanced composite clock/current/report must be zero"))
    else
        initialized && _stored_positive(time) && _stored_positive(dt) && dt <= time || throw(
            ArgumentError(
                "advanced composite state requires positive finite duration not exceeding time",
            ),
        )
        _sce_report_values(st, report, dt, total)
    end
    return nothing
end

function _sce_state_values(
    st::StaticCompositeElectrostaticPIC{T},
    species,
    E,
    phi,
    initialized,
    time,
    step,
    dt,
    right,
    impulse,
    report;
    map = true,
) where {T}
    E isa Vector{T} &&
        phi isa Vector{T} &&
        length(E) == length(st.mesh.widths) &&
        length(phi) == length(st.mesh.vertices) &&
        all(isfinite, E) &&
        all(isfinite, phi) ||
        throw(ArgumentError("composite fields must be finite ordinary matching Vectors"))
    _sce_values(st.mesh, species, st._factors)
    total = 0
    for ps in species
        total = Base.checked_add(total, length(ps.weight))
    end
    _sce_phase_values(st, initialized, time, step, dt, impulse, report, total)
    _sce_q_finite(right) &&
        right ==
        _yee_closed_right_plate(st.left_plate_charge, st.background_charge, st.mesh.L, species) ||
        throw(ArgumentError("composite derived right plate does not match original sources"))
    if map
        solution = _sce_solution(
            st.mesh,
            st.shape,
            st.background_charge,
            st.left_plate_charge,
            right,
            species,
            T,
        )
        isequal(E, solution.E) && isequal(phi, solution.phi) || throw(
            ArgumentError(
                "composite fields do not match the declared rounded original-source solve",
            ),
        )
        residuals = _sce_residuals(st.mesh, E, phi, solution.rhs)
        if step > 0
            report.max_gauss_residual == residuals.max_gauss &&
                report.max_constitutive_residual == residuals.max_constitutive ||
                throw(ArgumentError("composite report residuals do not match current fields"))
        end
        return solution
    end
    return nothing
end

function _sce_ready(st, species; map = true)
    _sce_mode()
    _sce_storage(st, species)
    _sce_configuration(st)
    st.initialized[] || throw(ArgumentError("call init! before composite stepping or diagnostics"))
    _sce_binding(st, species)
    _sce_values(st.mesh, species, st._factors)
    _sce_validate_bins(st._binding[].bins, species, st.mesh)
    return _sce_state_values(
        st,
        species,
        st.E,
        st.phi,
        st.initialized[],
        st.time[],
        st.step[],
        st.last_dt[],
        st.right_plate_charge[],
        st.last_impulse,
        st.last_report[];
        map,
    )
end

"""
    StaticCompositeElectrostaticPIC(mesh, species; shape, background_charge=0,
                                    left_plate_charge=0)

Static CPU 1D composite electrostatic PIC with physical triangular clouds and
one global Neumann solve on unique leaves/shared nodes. Positive-ion or signed
species use original finite q and positive finite m; factors are constructor-
bound. `init!` can atomically rebind matching-factor populations/counts/sources.
Fields and velocities are synchronous; positive dt is fixed after first success.
Physical walls are event-free, while coarse/fine interfaces may be crossed.
Only ordinary Vector particles are supported; one call at a time. Alias checks
require storage whose overlap Base.mightalias faithfully reports. Source,
compiler/runtime, self-force convergence and resource acceptance are distinct.
See `docs/STATIC_COMPOSITE_ELECTROSTATIC.md` for scope and exact defect ledgers.
"""
function StaticCompositeElectrostaticPIC(
    mesh::StaticAMRHierarchy1D{T},
    species::AbstractVector;
    shape::YeeTriangularCloud,
    background_charge = zero(T),
    left_plate_charge = zero(T),
) where {T}
    _sce_mode()
    _static_amr_validate(mesh)
    _sce_population_storage(species, T)
    factors = Tuple(_validated_charge_mass(ps) for ps in species)
    _sce_values(mesh, species, factors)
    b = _yl_source_scalar("background_charge", background_charge, T)
    left = _yl_source_scalar("left_plate_charge", left_plate_charge, T)
    right = _yee_closed_right_plate(left, b, mesh.L, species)
    binding = _sce_make_binding(species, species, mesh, T)
    report = Ref{Union{Nothing,StaticCompositeESReport{T}}}(nothing)
    st = StaticCompositeElectrostaticPIC{
        T,
        typeof(mesh),
        typeof(shape),
        typeof(b),
        typeof(left),
        typeof(factors),
    }(
        mesh,
        shape,
        b,
        left,
        zeros(T, length(mesh.widths)),
        zeros(T, length(mesh.vertices)),
        Ref(zero(T)),
        Ref(0),
        Ref(false),
        Ref(zero(T)),
        Ref(right),
        fill(zero(_CompositeQ), length(mesh.widths)),
        report,
        factors,
        Ref(binding),
    )
    _sce_storage(st, species)
    _sce_configuration(st)
    return st
end

function init!(st::StaticCompositeElectrostaticPIC{T}, species::AbstractVector) where {T}
    _sce_mode()
    _sce_storage(st, species)
    _sce_configuration(st)
    _sce_values(st.mesh, species, st._factors)
    binding = _sce_make_binding(species, species, st.mesh, T)
    right = _yee_closed_right_plate(st.left_plate_charge, st.background_charge, st.mesh.L, species)
    solution = _sce_solution(
        st.mesh,
        st.shape,
        st.background_charge,
        st.left_plate_charge,
        right,
        species,
        T,
    )
    zero_impulse = zero(_CompositeQ)
    copyto!(st.E, solution.E)
    copyto!(st.phi, solution.phi)
    fill!(st.last_impulse, zero_impulse)
    st._binding[] = binding
    st.right_plate_charge[] = right
    st.time[], st.step[], st.last_dt[] = zero(T), 0, zero(T)
    st.last_report[], st.initialized[] = nothing, true
    return st
end

function _sce_admit(st, species, tau)
    duration = _exact_rational(tau)
    a, hmin = _exact_rational(st.shape.half_width), minimum(st.mesh.widths)
    stiffness, charge = zero(_CompositeQ), zero(_CompositeQ)
    for ps in species
        q, m = _exact_rational(ps.q), _exact_rational(ps.m)
        for w in ps.weight
            weight = _exact_rational(w)
            stiffness += q^2 * weight / m
            charge += abs(q * weight)
        end
    end
    duration^2 * min(inv(hmin), inv(a)) * stiffness < 4 ||
        throw(ArgumentError("dt violates the strict composite zero-field rest-mode bound"))
    bound = max(
        maximum(x -> abs(_exact_rational(x)), st.E),
        abs(_exact_rational(st.left_plate_charge)) +
        abs(_exact_rational(st.background_charge)) * _exact_rational(st.mesh.L) +
        charge,
    )
    L = _exact_rational(st.mesh.L)
    for ps in species
        acceleration = abs(_boris_exact_factor(_validated_boris_factor(ps))) * bound
        for p in eachindex(ps.weight)
            x, v = _exact_rational(ps.x[1][p]), _exact_rational(ps.v[1][p])
            abs(v) * duration + acceleration * duration^2 / 2 < min(x, L - x) ||
                throw(ArgumentError("dt fails the strict composite physical-wall exclusion tube"))
        end
    end
    return nothing
end

function _sce_half_kick(ps, p, field, factor, half)
    v = ntuple(c -> ps.v[c][p], 3)
    z = zero(_CompositeQ)
    out = _yee_particle_half_kick(v, (field, z, z), (z, z, z), factor, half)
    all(isfinite, out) && _particle_stored_zero(out[2]) && _particle_stored_zero(out[3]) ||
        throw(ArgumentError("composite original-force kick must round finite and longitudinal"))
    return out
end

function _sce_field_energy(mesh, E)
    return sum(
        (mesh.widths[i] * _exact_rational(E[i])^2 / 2 for i in eachindex(E));
        init = zero(_CompositeQ),
    )
end

function step!(st::StaticCompositeElectrostaticPIC{T}, species::AbstractVector, dt::Real) where {T}
    _sce_mode()
    isfinite(dt) && dt >= zero(dt) ||
        throw(ArgumentError("composite dt must be finite and nonnegative"))
    tau = _validated_nonnegative_dt(T, dt; name = "static composite ES")
    _particle_stored_zero(tau) &&
        !iszero(dt) &&
        throw(ArgumentError("nonzero composite duration must not round to zero"))
    # The true-zero route checks read-only storage/phase but does not solve fields
    # or mutate any private numerical candidate.
    if _particle_stored_zero(tau)
        _sce_ready(st, species; map = false)
        return st
    end
    st.step[] == 0 ||
        isequal(tau, st.last_dt[]) ||
        throw(ArgumentError("composite dt is fixed after first success"))
    next_time, next_step = _validated_step_clock(st.time[], st.step[], tau)
    next_time > st.time[] || throw(ArgumentError("composite time must advance finitely"))
    _sce_ready(st, species)
    _sce_admit(st, species, tau)
    binding = st._binding[]
    candidates = binding.candidates
    half = _exact_rational(tau) / 2
    factors = [_boris_exact_factor(_validated_boris_factor(ps)) for ps in species]
    oldkinetic = _kinetic_energy_rational(species)
    for (ps, candidate, factor) in zip(species, candidates, factors), p in eachindex(ps.weight)
        x = ps.x[1][p]
        field = _sce_cloud_gather(st.E, _exact_rational(x), st.mesh, st.shape)
        velocity = _sce_half_kick(ps, p, field, factor, half)
        y = _particle_coordinate(x, velocity[1], tau)
        isfinite(y) || throw(ArgumentError("composite drift must round finite"))
        _static_amr_leaf(st.mesh, _exact_rational(y))
        candidate.x[1][p] = y
        for c = 1:3
            candidate.v[c][p] = velocity[c]
        end
    end
    _sce_values(st.mesh, candidates, st._factors)
    right = st.right_plate_charge[]
    solution = _sce_solution(
        st.mesh,
        st.shape,
        st.background_charge,
        st.left_plate_charge,
        right,
        candidates,
        T,
    )
    for (candidate, factor) in zip(candidates, factors), p in eachindex(candidate.weight)
        field = _sce_cloud_gather(solution.E, _exact_rational(candidate.x[1][p]), st.mesh, st.shape)
        velocity = _sce_half_kick(candidate, p, field, factor, half)
        for c = 1:3
            candidate.v[c][p] = velocity[c]
        end
    end
    _sce_values(st.mesh, candidates, st._factors)
    impulse = fill(zero(_CompositeQ), length(st.E))
    _sce_cloud_impulse!(impulse, species, candidates, st.mesh, st.shape)
    work = sum(
        (
            st.mesh.widths[i] *
            (_exact_rational(st.E[i]) + _exact_rational(solution.E[i])) *
            impulse[i] / 2 for i in eachindex(impulse)
        );
        init = zero(_CompositeQ),
    )
    fieldchange = _sce_field_energy(st.mesh, solution.E) - _sce_field_energy(st.mesh, st.E)
    kineticchange = _kinetic_energy_rational(candidates) - oldkinetic
    residuals = _sce_residuals(st.mesh, solution.E, solution.phi, solution.rhs)
    report = StaticCompositeESReport{T}(
        tau,
        work,
        fieldchange,
        kineticchange,
        fieldchange + work,
        kineticchange - work,
        fieldchange + kineticchange,
        residuals.max_gauss,
        residuals.max_constitutive,
        _sce_crossings(species, candidates, st.mesh),
    )
    bins = _sce_bins(candidates, st.mesh)
    prepared =
        _CompositeESBinding{T}(binding.owners, binding.components, candidates, binding.counts, bins)
    _sce_state_values(
        st,
        candidates,
        solution.E,
        solution.phi,
        true,
        next_time,
        next_step,
        tau,
        right,
        impulse,
        report;
        map = false,
    )
    # Every allocation, source solve, range/ownership/history check precedes this
    # ordinary-store publication. No public lengths or particle identities change.
    for (ps, candidate) in zip(species, candidates)
        copyto!(ps.x[1], candidate.x[1])
        foreach(copyto!, ps.v, candidate.v)
    end
    copyto!(st.E, solution.E)
    copyto!(st.phi, solution.phi)
    copyto!(st.last_impulse, impulse)
    st._binding[], st.last_report[] = prepared, report
    st.time[], st.step[], st.last_dt[] = next_time, next_step, tau
    return st
end

"""
    composite_electrostatic_diagnostics(state, species)

Exact physical charge/plate, unique-leaf kinetic/field energies, all-node Gauss
and constitutive residuals, source current and signed step defects. The report
separates field-solve rounding from particle work error; neither is forced zero.
No interface self-force, nonlinear stability or cost acceptance is implied.
"""
function composite_electrostatic_diagnostics(
    st::StaticCompositeElectrostaticPIC,
    species::AbstractVector,
)
    solution = _sce_ready(st, species)
    residuals = _sce_residuals(st.mesh, st.E, st.phi, solution.rhs)
    charge = sum(
        (_exact_rational(ps.q) * _exact_float_sum(ps.weight) for ps in species);
        init = zero(_CompositeQ),
    )
    kinetic, field = _kinetic_energy_rational(species), _sce_field_energy(st.mesh, st.E)
    volume_charge = copy(solution.rhs)
    volume_charge[1] -= _exact_rational(st.left_plate_charge)
    volume_charge[end] -= st.right_plate_charge[]
    density =
        _CompositeQ[volume_charge[i] / st.mesh.dual_volumes[i] for i in eachindex(volume_charge)]
    return (
        particle_charge = charge,
        background_charge = _exact_rational(st.background_charge) * _exact_rational(st.mesh.L),
        left_plate_charge = _exact_rational(st.left_plate_charge),
        right_plate_charge = st.right_plate_charge[],
        nodal_source = copy(solution.rhs),
        volume_nodal_charge = volume_charge,
        charge_density = density,
        dual_volumes = copy(st.mesh.dual_volumes),
        leaf_volumes = copy(st.mesh.widths),
        gauss_residual = residuals.gauss,
        constitutive_residual = residuals.constitutive,
        kinetic_energy = kinetic,
        field_energy = field,
        total_energy = kinetic + field,
        last_impulse = copy(st.last_impulse),
        last_report = deepcopy(st.last_report[]),
        velocity_time_offset = zero(st.time[]),
    )
end

function patch_electric_field(st::StaticCompositeElectrostaticPIC{T}, patch_id::Integer) where {T}
    _sce_ready(st, st._binding[].owners)
    1 <= patch_id <= length(st.mesh.patches) || throw(BoundsError(st.mesh.patches, patch_id))
    patch = st.mesh.patches[Int(patch_id)]
    out = Vector{T}(undef, patch.ncell)
    for cell = 1:patch.ncell
        lo, hi = patch.leaf_first[cell], patch.leaf_last[cell]
        value =
            sum(
                (st.mesh.widths[e] * _exact_rational(st.E[e]) for e = lo:hi);
                init = zero(_CompositeQ),
            ) / (st.mesh.vertices[hi+1] - st.mesh.vertices[lo])
        out[cell] = _from_exact(T, value)
        isfinite(out[cell]) || throw(ArgumentError("derived composite patch field is nonfinite"))
    end
    return out
end
