# Fixed-bit rational intervals for actual SBP quadrature. No local rounded
# reservoir totals or same-odd-sign theorem is used on signed differences.
function _open_hall_electron_terms(st::OpenHallHybrid{D,T}, f) where {D,T}
    terms = Dict{UInt64,_OpenHallQ}()
    temperature = _exact_rational(st.model.closure.Te)
    iszero(temperature) && return terms
    counts = Dict{Tuple{T,Bool},Int}()
    for I in CartesianIndices(f.n)
        key = (f.n[I], I[1] == 1 || I[1] == st.op.n[1])
        counts[key] = Base.checked_add(get(counts, key, 0), 1)
    end
    for ((value, endpoint), count) in counts
        density = _exact_rational(value)
        volume = _open_hall_volume(st.op, endpoint ? 1 : 2)
        _isothermal_add_density_terms!(terms, density, temperature * density * volume * count)
    end
    return terms
end

function _open_hall_energy_interval(st, f, species)
    terms = _open_hall_electron_terms(st, f)
    lower, upper = _isothermal_log_bounds(terms, st.controls.log_bits)
    offset = _kinetic_energy_rational(species) + _open_hall_magnetic(f.B, st.op)
    return lower + offset, upper + offset
end

# Merge signed density coefficients before enclosure, so identical endpoint
# logarithms cancel exactly rather than leaving two independent interval widths.
function _open_hall_energy_change_interval(st, old, endpoint, old_species, new_species)
    terms = Dict{UInt64,_OpenHallQ}()
    for (fields, sign) in ((old, -1), (endpoint, 1))
        for (base, coefficient) in _open_hall_electron_terms(st, fields)
            _isothermal_add_log!(terms, BigInt(base), sign * coefficient)
        end
    end
    lower, upper = _isothermal_log_bounds(terms, st.controls.log_bits)
    offset =
        _kinetic_energy_rational(new_species) - _kinetic_energy_rational(old_species) +
        _open_hall_magnetic(endpoint.B, st.op) - _open_hall_magnetic(old.B, st.op)
    return lower + offset, upper + offset
end

# Endpoint electron enthalpy estimate, using the declared represented reciprocal
# velocity expression ui_x-J_x*ninv as a complete rational diagnostic source.
# This does not introduce a second Ohm or particle-force kernel.
function _open_hall_enthalpy_terms(st, f)
    terms = Dict{UInt64,_OpenHallQ}()
    offset = zero(_OpenHallQ)
    Q = _exact_rational
    area = prod(Q, Base.tail(st.op.dx))
    temperature = Q(st.model.closure.Te)
    for side = 1:2, I in CartesianIndices(Base.tail(st.op.n))
        ix = side == 1 ? 1 : st.op.n[1]
        K = CartesianIndex((ix, Tuple(I)...))
        normal = side == 1 ? -1 : 1
        electron_velocity = Q(f.ui[1][K]) - Q(f.J[1][K]) * Q(f.ninv[K])
        flow = normal * area * electron_velocity
        density = Q(f.n[K])
        _isothermal_add_density_terms!(terms, density, flow * temperature * density)
        offset += flow * Q(f.pe[K])
    end
    return terms, offset
end
function _open_hall_enthalpy_interval(st, old, endpoint, h)
    first, a = _open_hall_enthalpy_terms(st, old)
    last, b = _open_hall_enthalpy_terms(st, endpoint)
    terms = Dict{UInt64,_OpenHallQ}()
    half = _exact_rational(h) / 2
    for group in (first, last), (base, coefficient) in group
        _isothermal_add_log!(terms, BigInt(base), half * coefficient)
    end
    lower, upper = _isothermal_log_bounds(terms, st.controls.log_bits)
    offset = half * (a + b)
    return lower + offset, upper + offset
end

function _open_hall_interval_value(::Type{T}, interval) where {T}
    lower, upper = map(x -> _from_exact(T, x), interval)
    return isfinite(lower) && isequal(lower, upper) ? lower : nothing
end

"""
    open_hall_diagnostics(model, species=model.owners)

Detached exact particle/source/work scalars and fixed-log-bit energy intervals.
The electron boundary enthalpy is an endpoint trapezoid estimate; the total
residual includes unresolved coupling/spatial/time/reservoir terms and is not
an exact conservation claim. Stored-norm and declared geometry are separate.
Private transforms may change while direct B constraints are checked.
"""
function open_hall_diagnostics(st::OpenHallHybrid{D,T}, species = st.owners) where {D,T}
    _acquire_open_hall!(st)
    try
        _open_hall_boundary(st, species)
        totals = _open_hall_particle_totals(species)
        energy = _open_hall_energy_interval(st, st.fields, species)
        grid_charge = _open_hall_grid_charge(st, st.fields)
        physical_volume = prod(_exact_rational, st.op.L)
        physical_face = prod(_exact_rational, Base.tail(st.op.L))
        quadrature_face = prod(BigInt(st.op.n[d]) * _exact_rational(st.op.dx[d]) for d = 2:D)
        norm_volume = sum(_exact_rational, st.op.sbp.H) * quadrature_face
        return (
            time = st.time[],
            step = st.step[],
            phase = :synchronous_integer,
            kinetic = totals.kinetic,
            magnetic = _open_hall_magnetic(st.fields.B, st.op),
            total_energy_interval = energy,
            certified_total = _open_hall_interval_value(T, energy),
            particle_totals = totals,
            represented_grid_charge = grid_charge,
            physical_particle_charge = totals.charge,
            grid_charge_defect = grid_charge - totals.charge,
            declared_volume = physical_volume,
            stored_norm_volume = norm_volume,
            declared_face_area = physical_face,
            quadrature_face_area = quadrature_face,
            volume_discrepancy = norm_volume - physical_volume,
            face_area_discrepancy = quadrature_face - physical_face,
            last_kick_defect = st.last_report === nothing ? zero(_OpenHallQ) :
                               st.last_report.kick_energy - st.last_report.kick_work,
            cumulative_kick_defect = st.cumulative.kick_energy - st.cumulative.kick_work,
            source = deepcopy(_open_hall_source_state(st)),
            last_report = deepcopy(st.last_report),
            log_bits = st.controls.log_bits,
            scope = :qualified_perpendicular_numerical_reservoir,
        )
    finally
        st.busy = false
    end
end
