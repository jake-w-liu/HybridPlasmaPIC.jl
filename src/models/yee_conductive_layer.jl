const _YeeLayerQ = Rational{BigInt}
const _YEE_LAYER_LEDGER_FIELDS = (
    :initial_energy, :last_loss, :total_loss, :last_work, :total_work,
    :last_defect, :total_defect,
)

"""
    YeeConductiveLayer(vacuum::MaxwellVacuum; dt, cells, peak_rate,
                       profile=:quadratic)

Fixed-step CPU Float32/64 scalar electric/magnetic conductivity absorber in
D=1:3. `cells` and `peak_rate` are tuples of (lower,upper) pairs, one per axis.
Profiles are `:quadratic` or `:cubic`, sampled at each component's Yee staggering;
rates add at edges/corners. Periodic-axis layers are disabled. This artificial
absorber is not PML/CPML or an exact outgoing boundary. Exterior PEC/periodic
conditions retain their meanings, and damping can change divergence in the layer.

Seed `vacuum.fields` and call `init!(layer)`. The supplied unadvanced vacuum is
borrowed exclusively; do not advance or mutate it outside this wrapper. Each
positive fixed step applies conductivity midpoint half-steps around the unchanged
synchronized Yee KDK candidate. Original rational factors round only final
fields. Strict `dt*sum(maximum face rate)/4 < 1` excludes midpoint sign reversal;
the vacuum CFL is also required. Neither condition is a global stability or
finite-grid reflection theorem. Nearest rounding and gradual underflow are required.

Actual removed field energy and prescribed-source work are retained separately
from the signed physical Yee balance defect. All fields, clock and exact budget
ledgers stage before publication. An expected input/arithmetic error changes only
private scratch. All-zero profiles retain the canonical vacuum field map.
Prescribed currents retain D1-only support; no particle/source-charge model is
supplied. `init!` resets the budget baseline from current finite PEC fields.

One instance requires exclusive access. Ordinary publication excludes resource
exhaustion, protected memory, custom setters and concurrent external mutation.
Checkpoint/archive state includes the complete synchronized fields and exact
ledgers; derived caches/staging are omitted. Runtime/resource acceptance is deferred.
"""
mutable struct YeeConductiveLayer{D,T}
    const inner::MaxwellVacuum{D,T}
    const dt::T
    const cells::NTuple{D,Tuple{Int,Int}}
    const peak_rate::NTuple{D,Tuple{T,T}}
    const profile::Symbol
    const _predamped::YeeFields{D,T}
    _axis_r::NTuple{D,Tuple{Vector{_YeeLayerQ},Vector{_YeeLayerQ}}}
    const busy::Threads.Atomic{Bool}
    phase::Symbol
    initialized::Bool
    ledger_step::Int
    initial_energy::_YeeLayerQ
    last_loss::_YeeLayerQ
    total_loss::_YeeLayerQ
    last_work::_YeeLayerQ
    total_work::_YeeLayerQ
    last_defect::_YeeLayerQ
    total_defect::_YeeLayerQ
end

function _yee_layer_mode()
    _boris_accuracy_mode() ||
        throw(ArgumentError("YeeConductiveLayer requires nearest rounding and gradual underflow"))
    return nothing
end

function _yee_layer_config(g::YeeGrid{D,T}, dt, cells, rates, profile) where {D,T}
    h = _require_finite_positive_real("dt", dt, T)
    widths, peaks = _yee_profile_config(g, cells, rates, profile)
    bound = _exact_rational(h) * sum(d -> max(_exact_rational(peaks[d][1]), _exact_rational(peaks[d][2])), 1:D) / 4
    bound < 1 ||
        throw(ArgumentError("layer dt*sum(maximum face rate)/4 must be strictly below one"))
    return h, widths, peaks
end

function _yee_profile_config(g::YeeGrid{D,T}, cells, rates, profile) where {D,T}
    profile in (:quadratic, :cubic) ||
        throw(ArgumentError("layer profile must be :quadratic or :cubic"))
    cells isa Tuple && length(cells) == D && rates isa Tuple && length(rates) == D ||
        throw(ArgumentError("cells and peak_rate need one (lower,upper) pair per axis"))
    widths = ntuple(D) do d
        pair = cells[d]
        pair isa Tuple && length(pair) == 2 ||
            throw(ArgumentError("each layer thickness must be a (lower,upper) pair"))
        all(n -> n isa Integer && !(n isa Bool) && 0 <= n < g.n[d], pair) ||
            throw(ArgumentError("layer cell counts must be nonnegative integers smaller than the axis"))
        low, high = Int(pair[1]), Int(pair[2])
        low < g.n[d] - high ||
            throw(ArgumentError("layer thicknesses must retain at least one interior cell"))
        (low, high)
    end
    peaks = ntuple(D) do d
        pair = rates[d]
        pair isa Tuple && length(pair) == 2 ||
            throw(ArgumentError("each peak_rate must be a (lower,upper) pair"))
        all(rate -> rate isa Real, pair) ||
            throw(ArgumentError("layer peak rates must be real values"))
        ntuple(2) do side
            rate = _require_finite_nonnegative_real("peak_rate", pair[side], T)
            active = widths[d][side] > 0
            (active ? rate > zero(T) : _particle_stored_zero(rate)) ||
                throw(ArgumentError("positive layer thickness requires positive rate; disabled faces require zero rate"))
            !g.periodic[d] || !active ||
                throw(ArgumentError("periodic axes cannot have conductive layers"))
            rate
        end
    end
    return widths, peaks
end

function _yee_layer_axis_value(g, cells, peaks, profile, hquarter, d, node, i)
    # BigInt formation precedes doubling, so valid Int lattice indices cannot
    # overflow while forming the exact half-index.
    xi = node ? BigInt(i - 1) // BigInt(1) : (2 * BigInt(i) - 1) // BigInt(2)
    power = profile === :quadratic ? 2 : 3
    value = zero(_YeeLayerQ)
    for side = 1:2
        width = cells[d][side]
        width == 0 && continue
        distance = side == 1 ? width - xi : xi - (g.n[d] - width)
        distance <= 0 && continue
        eta = distance / width
        value += _exact_rational(peaks[d][side]) * eta^power
    end
    return hquarter * value
end

function _yee_layer_axis_cache(g::YeeGrid{D}, h, cells, peaks, profile) where {D}
    hquarter = _exact_rational(h) / 4
    return _yee_profile_axis_cache(g, cells, peaks, profile, hquarter)
end

function _yee_profile_axis_cache(g::YeeGrid{D}, cells, peaks, profile, hquarter) where {D}
    return ntuple(D) do d
        nodes = g.n[d] + Int(!g.periodic[d])
        nodal = _YeeLayerQ[_yee_layer_axis_value(g, cells, peaks, profile, hquarter, d, true, i) for i = 1:nodes]
        centered = _YeeLayerQ[_yee_layer_axis_value(g, cells, peaks, profile, hquarter, d, false, i) for i = 1:g.n[d]]
        (nodal, centered)
    end
end

function YeeConductiveLayer(
    vacuum::MaxwellVacuum{D,T}; dt, cells, peak_rate, profile::Symbol = :quadratic,
) where {D,T}
    _yee_layer_mode()
    h, widths, peaks = _yee_layer_config(vacuum.g, dt, cells, peak_rate, profile)
    _require_finite_positive_real("c", vacuum.c, T)
    isequal(maxwell_cfl(vacuum), _yee_cfl(vacuum.g, vacuum.c)) ||
        throw(ArgumentError("layer vacuum CFL cache must match its grid and speed"))
    h < maxwell_cfl(vacuum) || throw(ArgumentError("layer dt must be strictly sub-CFL"))
    _yee_validate_storage(vacuum)
    _yee_validate_fields(vacuum.fields, vacuum.g)
    _particle_stored_zero(vacuum.time[]) && vacuum.step[] == 0 ||
        throw(ArgumentError("construct a conductive layer from an unadvanced vacuum"))
    z = zero(_YeeLayerQ)
    model = YeeConductiveLayer(
        vacuum, h, widths, peaks, profile, YeeFields(vacuum.g),
        _yee_layer_axis_cache(vacuum.g, h, widths, peaks, profile),
        Threads.Atomic{Bool}(false), :uninitialized, false, 0,
        z, z, z, z, z, z, z,
    )
    _yee_layer_validate_storage(model)
    return model
end

maxwell_cfl(model::YeeConductiveLayer) = maxwell_cfl(model.inner)

function _yee_layer_validate_config(model::YeeConductiveLayer; cache::Bool = true)
    _yee_layer_mode()
    inner = model.inner
    _require_finite_positive_real("c", inner.c, typeof(model.dt))
    isequal(maxwell_cfl(inner), _yee_cfl(inner.g, inner.c)) ||
        throw(ArgumentError("layer vacuum CFL cache must match its grid and speed"))
    _yee_layer_config(model.inner.g, model.dt, model.cells, model.peak_rate, model.profile)
    model.dt < maxwell_cfl(model) || throw(ArgumentError("layer dt must be strictly sub-CFL"))
    if cache
        g = model.inner.g
        hquarter = _exact_rational(model.dt) / 4
        for d = 1:length(g.n), slot = 1:2
            table = model._axis_r[d][slot]
            length(table) == g.n[d] + Int(slot == 1 && !g.periodic[d]) ||
                throw(DimensionMismatch("conductive layer axis cache has the wrong length"))
            for i in eachindex(table)
                expected = _yee_layer_axis_value(g, model.cells, model.peak_rate, model.profile, hquarter, d, slot == 1, i)
                isequal(table[i], expected) ||
                    throw(ArgumentError("conductive layer cache does not match its fixed profile"))
            end
        end
    end
    return nothing
end

function _yee_layer_validate_storage(model::YeeConductiveLayer; check_pec::Bool = true, cache::Bool = true)
    inner, g = model.inner, model.inner.g
    _yee_validate_storage(inner; check_pec)
    _yee_check_shape(model._predamped.E, g, true)
    _yee_check_shape(model._predamped.B, g, false)
    caches = cache ? Tuple(a for pair in model._axis_r for a in pair) : ()
    _yee_check_independent_arrays((
        inner.fields.E..., inner.fields.B..., inner._candidate.E..., inner._candidate.B...,
        model._predamped.E..., model._predamped.B..., caches...,
    ))
    return nothing
end

_yee_layer_energy(fields, inner) =
    _field_energy_rational(fields.E, fields.B, _exact_rational(inner.c)^2, inner.g.dx)

function _yee_layer_validate_budget(s, energy, time, step, dt::T, ::Type{T}) where {T}
    s.phase === :ready && s.initialized === true ||
        throw(ArgumentError("conductive layer must be initialized at a completed boundary"))
    _validate_checkpoint_clock_state((; time, step), T)
    s.ledger_step isa Int && s.ledger_step == step ||
        throw(ArgumentError("conductive layer ledger count must match its field clock"))
    for name in _YEE_LAYER_LEDGER_FIELDS
        value = getproperty(s, name)
        value isa _YeeLayerQ && denominator(value) > 0 ||
            throw(ArgumentError("conductive layer ledgers must be finite exact rational values"))
    end
    0 <= s.last_loss <= s.total_loss && s.initial_energy >= 0 ||
        throw(ArgumentError("conductive layer losses and initial energy must be nonnegative"))
    if step == 0
        _particle_stored_zero(time) && all(name -> iszero(getproperty(s, name)), _YEE_LAYER_LEDGER_FIELDS[2:end]) ||
            throw(ArgumentError("unadvanced layer clock and loss/work/defect ledgers must be zero"))
    else
        time > zero(T) && dt <= time || throw(ArgumentError("advanced layer time must be positive"))
        if step == 1
            s.last_loss == s.total_loss && s.last_work == s.total_work && s.last_defect == s.total_defect ||
                throw(ArgumentError("first-step layer ledgers must match their cumulative values"))
        end
    end
    energy + s.total_loss + s.total_work - s.initial_energy == s.total_defect ||
        throw(ArgumentError("conductive layer field and budget ledgers are inconsistent"))
    return nothing
end

function _yee_layer_boundary(model::YeeConductiveLayer{D,T}) where {D,T}
    _yee_validate_fields(model.inner.fields, model.inner.g)
    energy = _yee_layer_energy(model.inner.fields, model.inner)
    _yee_layer_validate_budget(model, energy, model.inner.time[], model.inner.step[], model.dt, T)
    return energy
end

function _yee_layer_acquire!(model::YeeConductiveLayer)
    Threads.atomic_cas!(model.busy, false, true) &&
        throw(ArgumentError("concurrent or reentrant conductive layer operation"))
    return nothing
end

function init!(model::YeeConductiveLayer{D,T}) where {D,T}
    _yee_layer_acquire!(model)
    try
        _yee_layer_validate_config(model; cache = false)
        _yee_layer_validate_storage(model; cache = false)
        _yee_validate_fields(model.inner.fields, model.inner.g)
        energy = _yee_layer_energy(model.inner.fields, model.inner)
        cache = _yee_layer_axis_cache(model.inner.g, model.dt, model.cells, model.peak_rate, model.profile)
        z = zero(_YeeLayerQ)
        model.inner.time[] = zero(T)
        model.inner.step[] = 0
        model._axis_r = cache
        model.initial_energy = energy
        model.last_loss = z
        model.total_loss = z
        model.last_work = z
        model.total_work = z
        model.last_defect = z
        model.total_defect = z
        model.ledger_step = 0
        model.initialized = true
        model.phase = :ready
        return model
    finally
        model.busy[] = false
    end
end

function _yee_layer_damp!(destination, original, model::YeeConductiveLayer{D,T}) where {D,T}
    for (electric, output, input) in ((true, destination.E, original.E), (false, destination.B, original.B))
        for c = 1:3, I in CartesianIndices(input[c])
            r = zero(_YeeLayerQ)
            for d = 1:D
                slot = _yee_node(electric, c, d) ? 1 : 2
                r += model._axis_r[d][slot][I[d]]
            end
            old = input[c][I]
            # Preserve original signed zeros and exact zero-profile identity.
            value = iszero(r) || _particle_stored_zero(old) ? old :
                    _from_exact(T, (1 - r) * _exact_rational(old) / (1 + r))
            isfinite(value) || throw(ArgumentError("conductive layer produces a nonfinite field"))
            output[c][I] = value
        end
    end
    return nothing
end

step_maxwell!(model::YeeConductiveLayer, dt::Real) = _yee_layer_step!(model, dt, nothing)
function step_maxwell!(model::YeeConductiveLayer{1}, dt::Real, impulse)
    impulse === nothing && throw(ArgumentError("impulse must be a tuple of three current arrays"))
    return _yee_layer_step!(model, dt, impulse)
end
step_maxwell!(model::YeeConductiveLayer, dt::Real, impulse) =
    throw(ArgumentError("prescribed conductive-layer currents support one spatial dimension"))

function _yee_layer_step!(model::YeeConductiveLayer{D,T}, dt::Real, impulse) where {D,T}
    h = _require_finite_positive_real("dt", dt, T)
    _yee_layer_acquire!(model)
    try
        _yee_layer_validate_config(model)
        _yee_layer_validate_storage(model)
        h == model.dt || throw(ArgumentError("conductive layer requires its fixed converted dt"))
        original_energy = _yee_layer_boundary(model)
        inner = model.inner
        next_time, next_step = _validated_step_clock(inner.time[], inner.step[], h)
        next_time > inner.time[] || throw(ArgumentError("conductive layer time must advance finitely"))
        if impulse !== nothing
            _yee_validate_impulse(impulse, inner.g)
            caches = Tuple(a for pair in model._axis_r for a in pair)
            _yee_check_independent_arrays((
                inner.fields.E..., inner.fields.B..., inner._candidate.E..., inner._candidate.B...,
                model._predamped.E..., model._predamped.B..., caches..., impulse...,
            ))
        end
        previous_phase = model.phase
        model.phase = :inflight
        completed = false
        try
            active = any(pair -> any(rate -> !_particle_stored_zero(rate), pair), model.peak_rate)
            before = active ? model._predamped : inner.fields
            if active
                _yee_layer_damp!(before, inner.fields, model)
            end
            before_energy = active ? _yee_layer_energy(before, inner) : original_energy
            candidate = inner._candidate
            _yee_candidate_step!(candidate, before, inner.g, h, inner.c, impulse)
            middle_energy = _yee_layer_energy(candidate, inner)
            work = impulse === nothing ? zero(_YeeLayerQ) :
                   _yee_source_work_rational(before.E, candidate.E, impulse, inner.g.dx[1])
            if active
                _yee_layer_damp!(candidate, candidate, model)
            end
            final_energy = active ? _yee_layer_energy(candidate, inner) : middle_energy
            loss = (original_energy - before_energy) + (middle_energy - final_energy)
            loss >= 0 || throw(ArgumentError("conductive layer damping must not add field energy"))
            defect = middle_energy - before_energy + work
            total_loss = model.total_loss + loss
            total_work = model.total_work + work
            total_defect = model.total_defect + defect
            budget = (
                phase = :ready, initialized = true, ledger_step = next_step,
                initial_energy = model.initial_energy,
                last_loss = loss, total_loss = total_loss,
                last_work = work, total_work = total_work,
                last_defect = defect, total_defect = total_defect,
            )
            _yee_layer_validate_budget(budget, final_energy, next_time, next_step, model.dt, T)
            # All field arithmetic, exact budgets and clock calculations passed.
            # Shared fixed-size copies and typed assignments only follow.
            _yee_publish_candidate!(inner, candidate, inner.fields, next_time, next_step)
            model.last_loss = loss
            model.total_loss = total_loss
            model.last_work = work
            model.total_work = total_work
            model.last_defect = defect
            model.total_defect = total_defect
            model.ledger_step = next_step
            completed = true
            return model
        finally
            model.phase = completed ? :ready : previous_phase
        end
    finally
        model.busy[] = false
    end
end

function _yee_layer_diagnostic(f, model::YeeConductiveLayer)
    _yee_layer_acquire!(model)
    try
        _yee_layer_validate_config(model)
        _yee_layer_validate_storage(model)
        _yee_layer_boundary(model)
        return f(model.inner)
    finally
        model.busy[] = false
    end
end

em_field_energy(model::YeeConductiveLayer) = _yee_layer_diagnostic(em_field_energy, model)
yee_poynting(model::YeeConductiveLayer) = _yee_layer_diagnostic(yee_poynting, model)
yee_poynting_flux(model::YeeConductiveLayer) = _yee_layer_diagnostic(yee_poynting_flux, model)

"""
    absorber_diagnostics(model::YeeConductiveLayer)

Return canonical rounded physical field energy, exact field energy and actual
last/cumulative damping loss, external-source work (positive out), and signed
physical Yee-map balance defect. Exact bookkeeping residual is distinct from
that generally nonzero defect. No boundary radiation or PML result is inferred.
"""
function absorber_diagnostics(model::YeeConductiveLayer)
    _yee_layer_acquire!(model)
    try
        _yee_layer_validate_config(model)
        _yee_layer_validate_storage(model)
        energy = _yee_layer_boundary(model)
        return (
            time = model.inner.time[], step = model.inner.step[],
            field_energy = em_field_energy(model.inner), field_energy_exact = energy,
            initial_energy = model.initial_energy,
            last_loss = model.last_loss, total_loss = model.total_loss,
            last_source_work = model.last_work, total_source_work = model.total_work,
            last_field_defect = model.last_defect, total_field_defect = model.total_defect,
            bookkeeping_residual = energy + model.total_loss + model.total_work - model.initial_energy - model.total_defect,
        )
    finally
        model.busy[] = false
    end
end
