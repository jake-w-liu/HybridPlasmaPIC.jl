struct ToroidalAnnulusStageReport{T}
    dt::T
    divergence::_CompositeQ
    divergence_tolerance::_CompositeQ
    nu::NTuple{2,_CompositeQ}
    kinetic_power::_CompositeQ
    particle_work::_CompositeQ
    grid_work::_CompositeQ
    current_work::_CompositeQ
    magnetic_power::_CompositeQ
    thermal_power::_CompositeQ
    pressure_work::_CompositeQ
    poynting::NTuple{2,_CompositeQ}
    enthalpy::NTuple{2,_CompositeQ}
    electric_controller::NTuple{2,_CompositeQ}
    thermal_controller::NTuple{2,_CompositeQ}
    lorentz_defect::_CompositeQ
    shape_defect::_CompositeQ
    ohm_defect::_CompositeQ
    magnetic_defect::_CompositeQ
    thermal_defect::_CompositeQ
    balance_defect::_CompositeQ
    logical_curl_defect::_CompositeQ
    pressure_assembly_defect::_CompositeQ
end

Base.isequal(a::ToroidalAnnulusStageReport, b::ToroidalAnnulusStageReport) =
    typeof(a) === typeof(b) &&
    all(name -> isequal(getfield(a, name), getfield(b, name)), fieldnames(typeof(a)))
Base.hash(report::ToroidalAnnulusStageReport, seed::UInt) =
    hash((typeof(report), ntuple(i -> getfield(report, i), fieldcount(typeof(report)))), seed)

_tah_dot(a, b) = sum((x * y for (x, y) in zip(a, b)); init = zero(_CompositeQ))
_tah_row(a, I) = ntuple(c -> _exact_rational(a[c][I]), 3)
function _tah_cross(a, b)
    return ntuple(
        c -> a[mod1(c + 1, 3)] * b[mod1(c + 2, 3)] - a[mod1(c + 2, 3)] * b[mod1(c + 1, 3)],
        3,
    )
end

function _tah_fields(::Type{T}, n) where {T}
    scalar() = zeros(T, n)
    vector() = ntuple(_ -> scalar(), 3)
    return (
        n = scalar(),
        Ji = vector(),
        ui = vector(),
        B = vector(),
        Eohm = vector(),
        E = vector(),
        Ecov = vector(),
        J = vector(),
        pe = scalar(),
        gradp = vector(),
        ninv = scalar(),
        ue = vector(),
    )
end
_tah_field_arrays(f) = (
    f.n,
    f.Ji...,
    f.ui...,
    f.B...,
    f.Eohm...,
    f.E...,
    f.Ecov...,
    f.J...,
    f.pe,
    f.gradp...,
    f.ninv,
    f.ue...,
)

function _tah_divergence(st, beta)
    mesh, op, work = st.mesh, st.mesh.op, st.work
    mixed_divergence!(work.divbeta, beta, op)
    Q = _exact_rational
    bounds = (
        2 / Q(op.dx[1]),
        sum(abs ∘ Q, op.transverse[1].ky; init = zero(_CompositeQ)),
        sum(abs ∘ Q, op.transverse[2].ky; init = zero(_CompositeQ)),
    )
    magnitude = ntuple(c -> maximum(x -> abs(Q(x)), beta[c]), 3)
    residual = scale = zero(_CompositeQ)
    for I in CartesianIndices(work.divbeta)
        _, J = _tah_metric(mesh, I)
        residual = max(residual, abs(Q(work.divbeta[I])) / J)
        scale = max(scale, _tah_dot(bounds, magnitude) / J)
    end
    tolerance = Q(st.div_atol) + Q(st.div_rtol) * scale
    residual <= tolerance ||
        throw(ArgumentError("annular magnetic flux failed its fresh physical-divergence gate"))
    return residual, tolerance
end

function _tah_pressure_source(st, pressure, I; flux = st.work.U[1][I])
    mesh = st.mesh
    side = I[1] == 1 ? 1 : I[1] == mesh.n[1] ? 2 : 0
    side == 0 && return zero(_CompositeQ), zero(_CompositeQ), 0
    _, J = _tah_metric(mesh, I)
    normal = side == 1 ? -1 : 1
    incoming = max(-normal * _exact_rational(flux), zero(_CompositeQ))
    nu = incoming / (J * _exact_rational(mesh.op.sbp.H[I[1]]))
    source = -nu * (_exact_rational(pressure) - _exact_rational(st.pressure_reservoir[side]))
    return source, nu, side
end

function _tah_volume_flux(mesh::ToroidalAnnulusGrid{T}, I, electron) where {T}
    basis, (metric, J) = _tah_basis(mesh, I), _tah_metric(mesh, I)
    return ntuple(
        c -> _tah_round(
            T,
            J * _tah_dot(basis[c], electron) / metric[c],
            "annular electron volume flux";
            preserve_nonzero = true,
        ),
        3,
    )
end

function _tah_contact_observation(st, f)
    mesh, Q = st.mesh, _exact_rational
    gamma = Q(st.closure.γ)
    nu, thermal, electric, poynting, enthalpy =
        ntuple(_ -> [zero(_CompositeQ), zero(_CompositeQ)], 5)
    for (face, side) in enumerate((1, mesh.n[1])), index in CartesianIndices(Base.tail(mesh.n))
        I = CartesianIndex((side, Tuple(index)...))
        basis, (metric, _) = _tah_basis(mesh, I), _tah_metric(mesh, I)
        volume = _tah_volume(mesh, I)
        E, Ji, j, ue = _tah_row(f.E, I), _tah_row(f.Ji, I), _tah_row(f.J, I), _tah_row(f.ue, I)
        flux = _tah_volume_flux(mesh, I, ue)
        source, frequency, _ = _tah_pressure_source(st, f.pe[I], I; flux = flux[1])
        nu[face] = max(nu[face], frequency)
        thermal[face] += volume * source / (gamma - 1)
        electric[face] +=
            volume * _tah_dot(ntuple(c -> Ji[c] - j[c], 3), ntuple(c -> E[c] - Q(f.Eohm[c][I]), 3))
        normal = face == 1 ? -1 : 1
        area = metric[2] * metric[3] * Q(mesh.dx[2]) * Q(mesh.dx[3])
        poynting[face] += normal * area * _tah_dot(_tah_cross(E, _tah_row(f.B, I)), basis[1])
        enthalpy[face] += normal * area * gamma * Q(f.pe[I]) * _tah_dot(ue, basis[1]) / (gamma - 1)
    end
    return (
        nu = Tuple(nu),
        thermal_controller = Tuple(thermal),
        electric_controller = Tuple(electric),
        poynting = Tuple(poynting),
        enthalpy = Tuple(enthalpy),
    )
end

function _tah_contact_fields(st, f)
    mesh, Q = st.mesh, _exact_rational
    for side in (1, mesh.n[1]), index in CartesianIndices(Base.tail(mesh.n))
        I = CartesianIndex((side, Tuple(index)...))
        _particle_stored_zero(f.Ecov[2][I]) && _particle_stored_zero(f.Ecov[3][I]) || throw(
            ArgumentError("annular contact tangential electric covectors must be stored zero"),
        )
        basis, (metric, _) = _tah_basis(mesh, I), _tah_metric(mesh, I)
        actual = _tah_row(f.E, I)
        for c = 1:3
            intended = Q(f.Ecov[c][I])
            observed = metric[c] * _tah_dot(basis[c], actual)
            scale = sum(abs(metric[c] * basis[c][d] * actual[d]) for d = 1:3) + abs(intended)
            abs(observed - intended) <= Q(st.refresh_rtol) * scale || throw(
                ArgumentError(
                    "annular controlled electric reconstruction failed its relative gate",
                ),
            )
        end
    end
    return nothing
end

const _TAH_ENDPOINT_RAW_POWER_KEYS = (
    :kinetic_power,
    :particle_work,
    :grid_work,
    :current_work,
    :magnetic_power,
    :thermal_power,
    :pressure_work,
)
const _TAH_ENDPOINT_FACE_POWER_KEYS =
    (:poynting, :enthalpy, :electric_controller, :thermal_controller)
const _TAH_ENDPOINT_DEFECT_KEYS = (
    :lorentz_defect,
    :shape_defect,
    :ohm_defect,
    :magnetic_defect,
    :thermal_defect,
    :balance_defect,
    :logical_curl_defect,
    :pressure_assembly_defect,
)

function _tah_endpoint_receipt(saved, refreshed, fields, st)
    _tah_report_valid(saved)
    _tah_report_valid(refreshed)
    contact = _tah_contact_observation(st, fields)
    all(
        name -> getfield(saved, name) == getproperty(contact, name),
        (:nu, _TAH_ENDPOINT_FACE_POWER_KEYS...),
    ) || throw(
        ArgumentError("annular endpoint contact powers/relaxation differ from its actual fields"),
    )
    saved.divergence_tolerance == refreshed.divergence_tolerance ||
        throw(ArgumentError("annular endpoint divergence source allowance differs"))
    # Native derivative refresh is qualified by source work, while the saved
    # endpoint's own exact algebra/contact data remain unchanged. Signed defect
    # cancellation never supplies its own relative-error denominator.
    scale = zero(_CompositeQ)
    for receipt in (saved, refreshed)
        for name in _TAH_ENDPOINT_RAW_POWER_KEYS
            scale += abs(getfield(receipt, name))
        end
        for name in _TAH_ENDPOINT_FACE_POWER_KEYS
            scale += sum(abs, getfield(receipt, name))
        end
    end
    allowance = _exact_rational(st.refresh_rtol) * scale
    for name in (_TAH_ENDPOINT_RAW_POWER_KEYS..., _TAH_ENDPOINT_DEFECT_KEYS...)
        abs(getfield(saved, name) - getfield(refreshed, name)) <= allowance || throw(
            ArgumentError("annular endpoint power/defect failed its source-work refresh gate"),
        )
    end
    for name in _TAH_ENDPOINT_FACE_POWER_KEYS
        all(
            pair -> abs(pair[1] - pair[2]) <= allowance,
            zip(getfield(saved, name), getfield(refreshed, name)),
        ) || throw(ArgumentError("annular endpoint face power failed its source-work refresh gate"))
    end
    return nothing
end

function _tah_prepare_fields!(st, stage, rate, h::T) where {T}
    mesh, w, Q = st.mesh, st.work, _exact_rational
    f = w.fields
    foreach(p -> _tah_particle_values(p, mesh), stage.species)
    all(a -> all(isfinite, a), stage.beta) && all(p -> isfinite(p) && p >= 0, stage.pe) ||
        throw(ArgumentError("annular primary field/pressure stage is invalid"))
    _tah_normal_flux(st, stage.beta)
    divergence, tolerance = _tah_divergence(st, stage.beta)
    _tah_moments!(st, stage.species)
    copyto!(f.pe, stage.pe)
    for I in CartesianIndices(f.n)
        basis, (metric, J) = _tah_basis(mesh, I), _tah_metric(mesh, I)
        physical = ntuple(c -> metric[c] * Q(stage.beta[c][I]) / J, 3)
        for c = 1:3
            f.B[c][I] = _tah_round(
                T,
                sum(basis[d][c] * physical[d] for d = 1:3),
                "annular physical magnetic field",
            )
            w.Bcov[c][I] =
                _tah_round(T, metric[c]^2 * Q(stage.beta[c][I]) / J, "annular magnetic covector")
        end
    end
    mixed_curl!(w.curlB, w.Bcov, mesh.op)
    mixed_gradient!(w.gradient, stage.pe, mesh.op)
    for I in CartesianIndices(f.n)
        basis, (metric, J) = _tah_basis(mesh, I), _tah_metric(mesh, I)
        for c = 1:3
            # sigma=-1: actual current flux is minus the logical curl.
            w.alpha[c][I] = -w.curlB[c][I]
        end
        for c = 1:3
            f.J[c][I] = _tah_round(
                T,
                sum(basis[d][c] * metric[d] * Q(w.alpha[d][I]) / J for d = 1:3),
                "annular physical current",
            )
            f.gradp[c][I] = _tah_round(
                T,
                sum(basis[d][c] * Q(w.gradient[d][I]) / metric[d] for d = 1:3),
                "annular physical pressure gradient",
            )
        end
    end
    risk =
        _ohm_field_rows!(f.Eohm, f.ui, f.B, f.J, w.zero_current, f.gradp, f.ninv, zero(T), zero(T))
    risk && throw(ArgumentError("annular original-product Ohm source-risk rejects before Faraday"))
    all(a -> all(isfinite, a), f.Eohm) || throw(ArgumentError("annular Ohm field is not finite"))
    electron_velocity!(f.ue, f.ui, f.J, f.n; nfloor = minimum(f.n))
    all(a -> all(isfinite, a), f.ue) ||
        throw(ArgumentError("annular electron velocity is not finite"))
    for I in CartesianIndices(f.n)
        basis, (metric, J) = _tah_basis(mesh, I), _tah_metric(mesh, I)
        electric, electron = _tah_row(f.Eohm, I), _tah_row(f.ue, I)
        flux = _tah_volume_flux(mesh, I, electron)
        face = I[1] == 1 || I[1] == mesh.n[1]
        for c = 1:3
            covariant = metric[c] * _tah_dot(basis[c], electric)
            f.Ecov[c][I] =
                face && c != 1 ? zero(T) : _tah_round(T, covariant, "annular electric covector")
            w.U[c][I] = flux[c]
        end
        for c = 1:3
            f.E[c][I] =
                face ?
                _tah_round(
                    T,
                    basis[1][c] * Q(f.Ecov[1][I]) / metric[1],
                    "annular controlled physical electric field",
                ) : f.Eohm[c][I]
        end
    end
    _tah_contact_fields(st, f)
    mixed_curl!(w.curlE, f.Ecov, mesh.op)
    for c = 1:3
        # beta_dot=-sigma*curl(Ecov)=+curl(Ecov).
        copyto!(rate.beta[c], w.curlE[c])
    end
    for side in (1, mesh.n[1]), index in CartesianIndices(Base.tail(mesh.n))
        I = CartesianIndex((side, Tuple(index)...))
        _particle_stored_zero(rate.beta[1][I]) ||
            throw(ArgumentError("annular electrical contact changed its frozen normal flux rate"))
    end
    mixed_divergence!(w.divue, w.U, mesh.op)
    nu = [zero(_CompositeQ), zero(_CompositeQ)]
    thermal = [zero(_CompositeQ), zero(_CompositeQ)]
    assembly = zero(_CompositeQ)
    gamma = Q(st.closure.γ)
    for I in CartesianIndices(f.n)
        _, J = _tah_metric(mesh, I)
        w.divue[I] = _tah_round(T, Q(w.divue[I]) / J, "annular electron divergence")
        # Every operand, including gamma, enters the canonical primitive as Q.
        raw = _electron_pressure_rate(
            (Q(f.pe[I]),),
            ntuple(c -> (Q(f.ue[c][I]),), 3),
            ntuple(c -> (Q(f.gradp[c][I]),), 3),
            (Q(w.divue[I]),),
            gamma,
            1,
            Val(3),
        )
        source, frequency, side = _tah_pressure_source(st, f.pe[I], I)
        Q(h) * frequency <= 1 ||
            throw(ArgumentError("annular h*nu exceeds the SAT-only relaxation bound"))
        rate.pe[I] =
            _tah_round(T, raw + source, "annular complete pressure rate"; preserve_nonzero = true)
        volume = _tah_volume(mesh, I)
        assembly += volume * (Q(rate.pe[I]) - raw - source) / (gamma - 1)
        if side != 0
            nu[side] = max(nu[side], frequency)
            thermal[side] += volume * source / (gamma - 1)
        end
    end
    grid = current = magnetic = pressure = A = zero(_CompositeQ)
    for I in CartesianIndices(f.n)
        basis, (metric, J) = _tah_basis(mesh, I), _tah_metric(mesh, I)
        volume = _tah_volume(mesh, I)
        E, Ji, j, ue = _tah_row(f.E, I), _tah_row(f.Ji, I), _tah_row(f.J, I), _tah_row(f.ue, I)
        grid += volume * _tah_dot(Ji, E)
        current += volume * _tah_dot(j, E)
        A += volume * _tah_dot(ue, _tah_row(f.gradp, I))
        pressure += volume * Q(rate.pe[I]) / (gamma - 1)
        logical_volume = volume / J
        magnetic +=
            logical_volume *
            sum(metric[c]^2 * Q(stage.beta[c][I]) * Q(rate.beta[c][I]) / J for c = 1:3)
    end
    contact = _tah_contact_observation(st, f)
    Tuple(nu) == contact.nu && Tuple(thermal) == contact.thermal_controller || throw(
        ArgumentError("annular consumed pressure source differs from its contact observation"),
    )
    budget = mixed_magnetic_rate_budget(w.Bcov, f.Ecov, w.alpha, rate.beta, mesh.op)
    # The canonical ledger's face orientation is positive; our order has sigma=-1.
    curl_defect = budget.magnetic_power + budget.current_work - budget.outward_poynting
    return (
        divergence = divergence,
        tolerance = tolerance,
        nu = Tuple(nu),
        grid = grid,
        current = current,
        magnetic = magnetic,
        pressure = pressure,
        A = A,
        poynting = contact.poynting,
        enthalpy = contact.enthalpy,
        electric = contact.electric_controller,
        thermal = Tuple(thermal),
        curl_defect = curl_defect,
        assembly = assembly,
    )
end

function _tah_stage_report(::Type{T}, h, data, kinetic, particle) where {T}
    lorentz = kinetic - particle
    shape = particle - data.grid
    ohm = data.grid - data.current + data.A - sum(data.electric)
    magnetic = data.magnetic + data.current + sum(data.poynting)
    thermal = data.pressure - sum(data.thermal) - data.A + sum(data.enthalpy)
    return ToroidalAnnulusStageReport{T}(
        h,
        data.divergence,
        data.tolerance,
        data.nu,
        kinetic,
        particle,
        data.grid,
        data.current,
        data.magnetic,
        data.pressure,
        data.A,
        data.poynting,
        data.enthalpy,
        data.electric,
        data.thermal,
        lorentz,
        shape,
        ohm,
        magnetic,
        thermal,
        lorentz + shape + ohm + magnetic + thermal,
        data.curl_defect,
        data.assembly,
    )
end
