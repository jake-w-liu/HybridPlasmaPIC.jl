const _YEE_PML_LEDGER_FIELDS = (
    :initial_energy, :last_exchange, :total_exchange,
    :last_field_defect, :total_field_defect,
)

struct _YeePMLMemory{D,T}
    E::NTuple{3,NTuple{2,Union{Nothing,Array{T,D}}}}
    B::NTuple{3,NTuple{2,Union{Nothing,Array{T,D}}}}
end

"""
    YeePML(vacuum::MaxwellVacuum; dt, cells, peak_rate, profile=:quadratic)

Source-free serial Float32/64 directional ADE PML in D1–3, with kappa=1 and
alpha=0. Each enabled curl derivative retains its own causal memory at its
E/B target staggering. This is not the scalar-conductivity layer. Profiles
use (lower,upper) cell/rate pairs per axis; periodic faces must be disabled.

Seed the unadvanced, exclusively borrowed vacuum and call `init!(pml)` to reset
all memories to zero, clock and budget. This declares zero prehistory, without
inferring prior radiation from seeded fields. Fixed positive converted steps
use original-duration midpoint B-half/E-full/B-half maps. Strict sub-CFL and
h*maximum_axis_rate/2<1 are admission policies, not a coupled PML stability or
finite-grid reflection theorem. Nearest rounding/gradual underflow are required.
D1 Ex/Bx have no derivative/memory and retain their original values.

Physical fields, memories, clock and signed exchange/defect budgets stage before
ordinary publication. Auxiliary exchange is not nonnegative heat. Conductivity
can change ordinary divergence in the artificial layer; no particle/source,
projection, general CPML or material model is supplied. Checkpoints retain every
advanced memory verbatim. One instance requires exclusive access; protected
memory, custom setters, resource exhaustion and concurrent external mutation are
outside ordinary-copy guarantees. Runtime/compiler/resource acceptance is deferred.
"""
mutable struct YeePML{D,T}
    const inner::MaxwellVacuum{D,T}
    const dt::T
    const cells::NTuple{D,Tuple{Int,Int}}
    const peak_rate::NTuple{D,Tuple{T,T}}
    const profile::Symbol
    const _memory::_YeePMLMemory{D,T}
    const _candidate_memory::_YeePMLMemory{D,T}
    _axis_rate::NTuple{D,Tuple{Vector{_YeeLayerQ},Vector{_YeeLayerQ}}}
    const busy::Threads.Atomic{Bool}
    phase::Symbol
    initialized::Bool
    ledger_step::Int
    initial_energy::_YeeLayerQ
    last_exchange::_YeeLayerQ
    total_exchange::_YeeLayerQ
    last_field_defect::_YeeLayerQ
    total_field_defect::_YeeLayerQ
end

function _yp_mode()
    _boris_accuracy_mode() ||
        throw(ArgumentError("YeePML requires nearest rounding and gradual underflow"))
    return nothing
end

function _yp_config(g::YeeGrid{D,T}, dt, cells, rates, profile) where {D,T}
    h = _require_finite_positive_real("dt", dt, T)
    widths, peaks = _yee_profile_config(g, cells, rates, profile)
    bound = _exact_rational(h) * maximum(d -> max(_exact_rational(peaks[d][1]), _exact_rational(peaks[d][2])), 1:D) / 2
    bound < 1 || throw(ArgumentError("PML dt*maximum axis rate/2 must be strictly below one"))
    return h, widths, peaks
end

function _yp_make_memory(g::YeeGrid{D,T}, cells) where {D,T}
    families = map((true, false)) do electric
        ntuple(3) do c
            ntuple(2) do slot
                axis = mod1(c + slot, 3)
                axis <= D && any(>(0), cells[axis]) ? zeros(T, _yee_shape(g, electric, c)) : nothing
            end
        end
    end
    return _YeePMLMemory{D,T}(families[1], families[2])
end

_yp_arrays(memory::_YeePMLMemory) =
    Tuple(a for family in (memory.E, memory.B) for pair in family for a in pair if a !== nothing)

function _yp_memory_storage(memory::_YeePMLMemory{D,T}, g::YeeGrid{D,T}, cells) where {D,T}
    for (electric, family) in ((true, memory.E), (false, memory.B)), c = 1:3, slot = 1:2
        axis = mod1(c + slot, 3)
        a = family[c][slot]
        present = axis <= D && any(>(0), cells[axis])
        if present
            a isa Array{T,D} && size(a) == _yee_shape(g, electric, c) ||
                throw(DimensionMismatch("PML memory must match its target Yee component"))
        else
            a === nothing || throw(ArgumentError("suppressed or disabled PML directions have no memory"))
        end
    end
    return nothing
end

function _yp_memory_values(memory, g, rates, step)
    for (electric, family) in ((true, memory.E), (false, memory.B)), c = 1:3, slot = 1:2
        a = family[c][slot]
        a === nothing && continue
        _yee_all_finite(a) || throw(ArgumentError("PML memories must be finite"))
        _yee_walls_zero(a, g, electric, c) ||
            throw(ArgumentError("PML memories must vanish on constrained PEC planes"))
        axis = mod1(c + slot, 3)
        table = rates[axis][_yee_node(electric, c, axis) ? 1 : 2]
        for I in CartesianIndices(a)
            if step == 0 || iszero(table[I[axis]])
                _particle_stored_zero(a[I]) ||
                    throw(ArgumentError("PML memories must vanish at step zero and at zero-rate samples"))
            end
        end
    end
    return nothing
end

function _yp_validate_config(model::YeePML; cache::Bool = true)
    _yp_mode()
    inner = model.inner
    _require_finite_positive_real("c", inner.c, typeof(model.dt))
    isequal(maxwell_cfl(inner), _yee_cfl(inner.g, inner.c)) ||
        throw(ArgumentError("PML vacuum CFL cache must match its grid and speed"))
    _yp_config(inner.g, model.dt, model.cells, model.peak_rate, model.profile)
    model.dt < maxwell_cfl(inner) || throw(ArgumentError("PML dt must be strictly sub-CFL"))
    if cache
        for d = 1:length(inner.g.n), slot = 1:2
            table = model._axis_rate[d][slot]
            length(table) == inner.g.n[d] + Int(slot == 1 && !inner.g.periodic[d]) ||
                throw(DimensionMismatch("PML axis cache has the wrong length"))
            for i in eachindex(table)
                expected = _yee_layer_axis_value(inner.g, model.cells, model.peak_rate, model.profile, one(_YeeLayerQ), d, slot == 1, i)
                isequal(table[i], expected) || throw(ArgumentError("PML cache does not match its fixed profile"))
            end
        end
    end
    return nothing
end

function _yp_validate_storage(model::YeePML; check_pec::Bool = true, cache::Bool = true)
    inner = model.inner
    _yee_validate_storage(inner; check_pec)
    _yp_memory_storage(model._memory, inner.g, model.cells)
    _yp_memory_storage(model._candidate_memory, inner.g, model.cells)
    caches = cache ? Tuple(a for pair in model._axis_rate for a in pair) : ()
    _yee_check_independent_arrays((
        inner.fields.E..., inner.fields.B..., inner._candidate.E..., inner._candidate.B...,
        _yp_arrays(model._memory)..., _yp_arrays(model._candidate_memory)..., caches...,
    ))
    return nothing
end

function YeePML(vacuum::MaxwellVacuum{D,T}; dt, cells, peak_rate, profile::Symbol = :quadratic) where {D,T}
    _yp_mode()
    h, widths, peaks = _yp_config(vacuum.g, dt, cells, peak_rate, profile)
    z = zero(_YeeLayerQ)
    model = YeePML(
        vacuum, h, widths, peaks, profile,
        _yp_make_memory(vacuum.g, widths), _yp_make_memory(vacuum.g, widths),
        _yee_profile_axis_cache(vacuum.g, widths, peaks, profile, one(_YeeLayerQ)),
        Threads.Atomic{Bool}(false), :uninitialized, false, 0, z, z, z, z, z,
    )
    _yp_validate_config(model)
    _yp_validate_storage(model)
    _yee_validate_fields(vacuum.fields, vacuum.g)
    _particle_stored_zero(vacuum.time[]) && vacuum.step[] == 0 ||
        throw(ArgumentError("construct PML from an unadvanced vacuum"))
    return model
end

maxwell_cfl(model::YeePML) = maxwell_cfl(model.inner)

function _yp_validate_budget(s, energy, time, step, ::Type{T}) where {T}
    s.phase === :ready && s.initialized === true ||
        throw(ArgumentError("PML requires an initialized completed boundary"))
    _validate_checkpoint_clock_state((; time, step), T)
    s.ledger_step isa Int && s.ledger_step == step ||
        throw(ArgumentError("PML ledger count must match its field clock"))
    for name in _YEE_PML_LEDGER_FIELDS
        value = getproperty(s, name)
        value isa _YeeLayerQ && denominator(value) > 0 ||
            throw(ArgumentError("PML ledgers must be finite exact rational values"))
    end
    s.initial_energy >= 0 || throw(ArgumentError("PML initial field norm must be nonnegative"))
    if step == 0
        _particle_stored_zero(time) && all(name -> iszero(getproperty(s, name)), _YEE_PML_LEDGER_FIELDS[2:end]) ||
            throw(ArgumentError("unadvanced PML clock and exchange/defect ledgers must be zero"))
    else
        time > zero(T) || throw(ArgumentError("advanced PML time must be positive"))
        if step == 1
            s.last_exchange == s.total_exchange && s.last_field_defect == s.total_field_defect ||
                throw(ArgumentError("first-step PML ledgers must match cumulative values"))
        end
    end
    energy + s.total_exchange - s.initial_energy == s.total_field_defect ||
        throw(ArgumentError("PML fields and exact budget ledgers are inconsistent"))
    return nothing
end

function _yp_boundary(model::YeePML{D,T}) where {D,T}
    _yee_validate_fields(model.inner.fields, model.inner.g)
    _yp_memory_values(model._memory, model.inner.g, model._axis_rate, model.inner.step[])
    energy = _yee_layer_energy(model.inner.fields, model.inner)
    _yp_validate_budget(model, energy, model.inner.time[], model.inner.step[], T)
    return energy
end

function _yp_acquire!(model::YeePML)
    Threads.atomic_cas!(model.busy, false, true) &&
        throw(ArgumentError("concurrent or reentrant PML operation"))
    return nothing
end

function _yp_copy_memory!(destination, source)
    for (out, input) in ((destination.E, source.E), (destination.B, source.B)), c = 1:3, slot = 1:2
        a = input[c][slot]
        a === nothing || copyto!(out[c][slot], a)
    end
    return nothing
end

function init!(model::YeePML{D,T}) where {D,T}
    _yp_acquire!(model)
    try
        _yp_validate_config(model; cache = false)
        _yp_validate_storage(model; cache = false)
        _yee_validate_fields(model.inner.fields, model.inner.g)
        energy = _yee_layer_energy(model.inner.fields, model.inner)
        rates = _yee_profile_axis_cache(model.inner.g, model.cells, model.peak_rate, model.profile, one(_YeeLayerQ))
        z, time = zero(_YeeLayerQ), zero(T)
        for a in _yp_arrays(model._candidate_memory)
            fill!(a, time)
        end
        _yp_memory_values(model._candidate_memory, model.inner.g, rates, 0)
        budget = (
            phase = :ready, initialized = true, ledger_step = 0,
            initial_energy = energy, last_exchange = z, total_exchange = z,
            last_field_defect = z, total_field_defect = z,
        )
        _yp_validate_budget(budget, energy, time, 0, T)
        _yp_copy_memory!(model._memory, model._candidate_memory)
        model.inner.time[] = time
        model.inner.step[] = 0
        model._axis_rate = rates
        model.initial_energy = energy
        model.last_exchange = z
        model.total_exchange = z
        model.last_field_defect = z
        model.total_field_defect = z
        model.ledger_step = 0
        model.initialized = true
        model.phase = :ready
        return model
    finally
        model.busy[] = false
    end
end

# A frozen-source midpoint subflow. Canonical exact differences retain original
# entries/dx, and delta is Q(h) or Q(h)/2, never an independently rounded duration.
function _yp_stage!(out, old, source, out_memory, old_memory, model::YeePML{D,T}, delta, c2, volume, electric) where {D,T}
    g = model.inner.g
    exchange = zero(_YeeLayerQ)
    beta = electric ? c2 : -one(_YeeLayerQ)
    for c = 1:3
        ja, jb = mod1(c + 1, 3), mod1(c + 2, 3)
        if ja > D && jb > D
            # D1 longitudinal components have no derivative or memory. Preserve
            # original values directly, including their stored signed zeros.
            copyto!(out[c], old[c])
            continue
        end
        for I in CartesianIndices(out[c])
            original = _exact_rational(old[c][I])
            terms = ntuple(2) do slot
                axis = slot == 1 ? ja : jb
                if axis > D
                    return (flux = zero(_YeeLayerQ), before = zero(_YeeLayerQ), after = nothing)
                end
                source_component = slot == 1 ? jb : ja
                derivative = _yee_difference_exact(source[source_component], I, axis, g, !electric)
                a = old_memory[c][slot]
                psi = a === nothing ? zero(_YeeLayerQ) : _exact_rational(a[I])
                lambda = model._axis_rate[axis][_yee_node(electric, c, axis) ? 1 : 2][I[axis]]
                r = delta * lambda / 2
                flux = (derivative + psi) / (1 + r)
                value = if a === nothing
                    nothing
                else
                    next = _from_exact(T, ((1 - r) * psi - 2r * derivative) / (1 + r))
                    isfinite(next) || throw(ArgumentError("PML stage produces a nonfinite memory"))
                    next
                end
                (flux = flux, before = psi, after = value)
            end
            value = _yee_source_row_exact(original + beta * delta * (terms[1].flux - terms[2].flux), nothing, I, T)
            average_field = (original + _exact_rational(value)) / 2
            average_memory = map(terms) do term
                term.after === nothing ? zero(_YeeLayerQ) : (term.before + _exact_rational(term.after)) / 2
            end
            exchange += (electric ? -c2 : c2) * delta * average_field * (average_memory[1] - average_memory[2])
            out[c][I] = value
            for slot = 1:2
                next = terms[slot].after
                next === nothing || (out_memory[c][slot][I] = next)
            end
        end
    end
    return volume * exchange
end

step_maxwell!(model::YeePML, dt::Real) = _yp_step!(model, dt)
step_maxwell!(model::YeePML, dt::Real, impulse) =
    throw(ArgumentError("this Yee PML supports source-free evolution only"))

function _yp_step!(model::YeePML{D,T}, dt::Real) where {D,T}
    h = _require_finite_positive_real("dt", dt, T)
    _yp_acquire!(model)
    try
        _yp_validate_config(model)
        _yp_validate_storage(model)
        h == model.dt || throw(ArgumentError("PML requires its fixed converted dt"))
        initial = _yp_boundary(model)
        inner = model.inner
        next_time, next_step = _validated_step_clock(inner.time[], inner.step[], h)
        next_time > inner.time[] || throw(ArgumentError("PML time must advance finitely"))
        duration = _exact_rational(h)
        half = duration / 2
        c2 = _exact_rational(inner.c)^2
        volume = prod(_exact_rational, inner.g.dx)
        previous_phase = model.phase
        completed = false
        model.phase = :inflight
        try
            candidate, old = inner._candidate, inner.fields
            memory, staged = model._memory, model._candidate_memory
            if any(pair -> any(>(0), pair), model.cells)
                exchange = _yp_stage!(candidate.B, old.B, old.E, staged.B, memory.B, model, half, c2, volume, false)
                exchange += _yp_stage!(candidate.E, old.E, candidate.B, staged.E, memory.E, model, duration, c2, volume, true)
                exchange += _yp_stage!(candidate.B, candidate.B, candidate.E, staged.B, staged.B, model, half, c2, volume, false)
            else
                _yee_candidate_step!(candidate, old, inner.g, h, inner.c, nothing)
                exchange = zero(_YeeLayerQ)
            end
            _yp_memory_values(staged, inner.g, model._axis_rate, next_step)
            energy = _yee_layer_energy(candidate, inner)
            defect = energy - initial + exchange
            total_exchange = model.total_exchange + exchange
            total_defect = model.total_field_defect + defect
            budget = (
                phase = :ready, initialized = true, ledger_step = next_step,
                initial_energy = model.initial_energy,
                last_exchange = exchange, total_exchange = total_exchange,
                last_field_defect = defect, total_field_defect = total_defect,
            )
            _yp_validate_budget(budget, energy, next_time, next_step, T)
            # Complete candidate fields/memories, clock and exact ledgers passed.
            # Only prevalidated fixed-size ordinary copies and typed stores follow.
            _yp_copy_memory!(memory, staged)
            _yee_publish_candidate!(inner, candidate, old, next_time, next_step)
            model.last_exchange = exchange
            model.total_exchange = total_exchange
            model.last_field_defect = defect
            model.total_field_defect = total_defect
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

function _yp_diagnostic(f, model::YeePML)
    _yp_acquire!(model)
    try
        _yp_validate_config(model)
        _yp_validate_storage(model)
        _yp_boundary(model)
        return f(model.inner)
    finally
        model.busy[] = false
    end
end

em_field_energy(model::YeePML) = _yp_diagnostic(em_field_energy, model)
yee_poynting(model::YeePML) = _yp_diagnostic(yee_poynting, model)
yee_poynting_flux(model::YeePML) = _yp_diagnostic(yee_poynting_flux, model)

"""
    pml_diagnostics(model::YeePML)

Return canonical E/B field norm, signed last/cumulative memory exchange and
physical field-balance defect, plus exact memory square sums. The latter have
derivative units and are not positive augmented material energy or heat. The
exact bookkeeping residual does not imply zero field defect or finite-grid
reflectionlessness. Outer Poynting flux keeps actual PEC/periodic semantics.
"""
function pml_diagnostics(model::YeePML)
    _yp_acquire!(model)
    try
        _yp_validate_config(model)
        _yp_validate_storage(model)
        energy = _yp_boundary(model)
        memory_norm(family) = _field_square_sum_exact(Tuple(a for pair in family for a in pair if a !== nothing))
        return (
            time = model.inner.time[], step = model.inner.step[],
            field_energy = em_field_energy(model.inner), field_energy_exact = energy,
            electric_memory_square_sum = memory_norm(model._memory.E),
            magnetic_memory_square_sum = memory_norm(model._memory.B),
            initial_energy = model.initial_energy,
            last_exchange = model.last_exchange, total_exchange = model.total_exchange,
            last_field_defect = model.last_field_defect, total_field_defect = model.total_field_defect,
            bookkeeping_residual = energy + model.total_exchange - model.initial_energy - model.total_field_defect,
        )
    finally
        model.busy[] = false
    end
end
