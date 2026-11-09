"Receipt of a cold composite Q1 solve, independently checked in original charge units."
struct AdaptiveCompositeESSolveReport{T<:Union{Float32,Float64}}
    iterations::Int
    status::String
    solved::Bool
    max_residual::_CompositeQ
    source_norm::_CompositeQ
    tolerance::_CompositeQ
    scale::_CompositeQ
end
struct _AdaptiveQ1Data{M}
    mesh::M
    elements::Vector{Matrix{_CompositeQ}}
    diagonal::Vector{_CompositeQ}
    scale::_CompositeQ
end
_acq_max(a) = maximum(abs, a; init = zero(_CompositeQ))
function _acq_element(lower::NTuple{D}, upper) where {D}
    bits = _camr_bits(Val(D))
    h = ntuple(d -> upper[d] - lower[d], D)
    K = Matrix{_CompositeQ}(undef, length(bits), length(bits))
    for (a, aa) in enumerate(bits), (b, bb) in enumerate(bits)
        value = zero(_CompositeQ)
        for d = 1:D
            term = ((iszero(aa[d]) ? -1 : 1) * (iszero(bb[d]) ? -1 : 1)) / h[d]
            for j = 1:D
                j == d && continue
                term *= h[j] * (aa[j] == bb[j] ? 2 : 1) / 6
            end
            value += term
        end
        K[a,b] = value
    end
    return K
end
function _acq_data(mesh)
    elements = [_acq_element(mesh.lower[k], mesh.upper[k]) for k in eachindex(mesh.addresses)]
    diagonal = fill(zero(_CompositeQ), length(mesh.free))
    for k in eachindex(elements), (a, raw_a) in enumerate(mesh.corners[k]), (b, raw_b) in enumerate(mesh.corners[k])
        coefficient = elements[k][a,b]
        for (master_a, pa) in mesh.rows[raw_a], (master_b, pb) in mesh.rows[raw_b]
            master_a == master_b || continue
            f = mesh.free_index[master_a]
            iszero(f) || (diagonal[f] += pa * coefficient * pb)
        end
    end
    all(>(zero(_CompositeQ)), diagonal) || throw(ArgumentError("composite grounded Q1 diagonal is not positive"))
    return _AdaptiveQ1Data(mesh, elements, diagonal, maximum(diagonal))
end
function _acq_validate(data, mesh)
    data.mesh === mesh || throw(ArgumentError("composite Q1 operator mesh binding changed"))
    expected = _acq_data(mesh)
    isequal(data.elements, expected.elements) && isequal(data.diagonal, expected.diagonal) && data.scale == expected.scale ||
        throw(ArgumentError("composite original Q1 operator cache changed"))
    return nothing
end
function _acq_master_potential(mesh, phi)
    length(phi) == length(mesh.free) || throw(DimensionMismatch("composite free potential length differs"))
    out = fill(zero(_CompositeQ), length(mesh.masters))
    for (f, master) in enumerate(mesh.free)
        out[master] = _exact_rational(phi[f])
    end
    return out
end
function _acq_raw_potential(mesh, phi)
    u = _acq_master_potential(mesh, phi)
    return [sum(pair -> pair.second * u[pair.first], row; init = zero(_CompositeQ)) for row in mesh.rows]
end
function _acq_apply(data, u)
    mesh = data.mesh
    length(u) == length(mesh.masters) || throw(DimensionMismatch("composite master operand length differs"))
    raw = [sum(pair -> pair.second * u[pair.first], row; init = zero(_CompositeQ)) for row in mesh.rows]
    out = fill(zero(_CompositeQ), length(mesh.masters))
    for k in eachindex(data.elements), (a, v) in enumerate(mesh.corners[k])
        value = sum(b -> data.elements[k][a,b] * raw[mesh.corners[k][b]], eachindex(mesh.corners[k]); init = zero(_CompositeQ))
        for (master, coefficient) in mesh.rows[v]
            out[master] += coefficient * value
        end
    end
    return out
end
_acq_free_apply(data, phi) = _acq_apply(data, _acq_master_potential(data.mesh, phi))[data.mesh.free]
struct _AdaptiveQ1Operator{T,O}
    data::O
end
struct _AdaptiveQ1Jacobi{T,O}
    data::O
end
Base.size(A::_AdaptiveQ1Operator) = (length(A.data.diagonal), length(A.data.diagonal))
Base.size(A::_AdaptiveQ1Jacobi) = (length(A.data.diagonal), length(A.data.diagonal))
Base.size(A::Union{_AdaptiveQ1Operator,_AdaptiveQ1Jacobi}, d::Integer) = d <= 2 ? length(A.data.diagonal) : 1
Base.eltype(::_AdaptiveQ1Operator{T}) where {T} = T
Base.eltype(::_AdaptiveQ1Jacobi{T}) where {T} = T
function LinearAlgebra.mul!(out::Vector{T}, A::_AdaptiveQ1Operator{T}, input::Vector{T}) where {T}
    value = _acq_free_apply(A.data, input)
    for i in eachindex(out)
        out[i] = _from_exact(T, value[i] / A.data.scale)
        isfinite(out[i]) || throw(ArgumentError("composite Q1 matvec exceeds represented range"))
    end
    return out
end
function LinearAlgebra.mul!(out::Vector{T}, A::_AdaptiveQ1Jacobi{T}, input::Vector{T}) where {T}
    for i in eachindex(out)
        out[i] = _from_exact(T, _exact_rational(input[i]) * A.data.scale / A.data.diagonal[i])
        isfinite(out[i]) || throw(ArgumentError("composite Q1 preconditioner exceeds represented range"))
    end
    return out
end
function _acq_residual(data, phi, source)
    all = _acq_apply(data, _acq_master_potential(data.mesh, phi)) - source
    return all, all[data.mesh.free]
end
function _acq_solve(data, source, rtol::T, atol::T, maxiter::Int) where {T}
    data.mesh.periodic && !iszero(sum(source;init = zero(_CompositeQ))) && throw(ArgumentError("periodic Q1 solve requires exact original neutrality"))
    n = length(data.diagonal)
    iterations = iszero(maxiter) ? Base.checked_mul(2, n) : maxiter
    observed = data.mesh.periodic ? source : source[data.mesh.free]
    target = _exact_rational(atol) + _exact_rational(rtol) * _acq_max(observed)
    divisor = data.scale * (data.mesh.periodic ? n : 1)
    scaled = _from_exact(T, target / divisor)
    isfinite(scaled) || throw(ArgumentError("composite Q1 scaled tolerance exceeds represented range"))
    _exact_rational(scaled) > target / divisor && (scaled = prevfloat(scaled))
    rhs = [_from_exact(T, source[master] / data.scale) for master in data.mesh.free]
    all(isfinite, rhs) || throw(ArgumentError("composite Q1 scaled source exceeds represented range"))
    solver = Krylov.CgWorkspace(n, n, Vector{T})
    operator = _AdaptiveQ1Operator{T,typeof(data)}(data)
    jacobi = _AdaptiveQ1Jacobi{T,typeof(data)}(data)
    Krylov.cg!(solver, operator, rhs; M = jacobi, ldiv = false, radius = zero(T), linesearch = false, history = false, rtol = zero(T), atol = scaled, itmax = iterations)
    phi = copy(Krylov.solution(solver))
    all(isfinite, phi) || throw(ArgumentError("composite Q1 potential is nonfinite"))
    residual, free = _acq_residual(data, phi, source)
    admitted = data.mesh.periodic ? residual : free
    valid = solver.stats.solved && !solver.stats.inconsistent && !solver.stats.indefinite && _acq_max(admitted) <= target
    valid || throw(ArgumentError("composite Q1 solve failed original residual/status admission: $(solver.stats.status), residual=$(_acq_max(admitted)), target=$target"))
    report = AdaptiveCompositeESSolveReport{T}(solver.stats.niter, solver.stats.status, true, _acq_max(admitted), _acq_max(observed), target, data.scale)
    return phi, report
end
function _acq_report(data, phi, source, report, rtol::T, atol::T, maxiter) where {T}
    data.mesh.periodic && !iszero(sum(source;init = zero(_CompositeQ))) && throw(ArgumentError("periodic Q1 receipt requires exact original neutrality"))
    report isa AdaptiveCompositeESSolveReport{T} || throw(ArgumentError("composite solve receipt type differs"))
    all_residual, free_residual = _acq_residual(data, phi, source)
    residual = data.mesh.periodic ? all_residual : free_residual
    norm = _acq_max(data.mesh.periodic ? source : source[data.mesh.free])
    target = _exact_rational(atol) + _exact_rational(rtol) * norm
    limit = iszero(maxiter) ? Base.checked_mul(2, length(phi)) : maxiter
    report.solved && 0 <= report.iterations <= limit && !isempty(report.status) &&
    all(_sce_q_finite, (report.max_residual, report.source_norm, report.tolerance, report.scale)) &&
    report.max_residual == _acq_max(residual) <= target && report.source_norm == norm && report.tolerance == target && report.scale == data.scale ||
        throw(ArgumentError("composite solve receipt does not match original represented field/source"))
    return nothing
end
function _acq_energy(data, phi)
    u = _acq_master_potential(data.mesh, phi)
    image = _acq_apply(data, u)
    return sum(i -> u[i] * image[i], eachindex(u); init = zero(_CompositeQ)) / 2
end
function _acq_field(mesh::CartesianAMRHierarchy{D}, raw, k::Int, x) where {D}
    bits = _camr_bits(Val(D))
    return ntuple(3) do d
        d > D && return zero(_CompositeQ)
        value = zero(_CompositeQ)
        for (a, vertex) in enumerate(mesh.corners[k])
            coefficient = (iszero(bits[a][d]) ? -one(_CompositeQ) : one(_CompositeQ)) / (mesh.upper[k][d] - mesh.lower[k][d])
            for j = 1:D
                j == d && continue
                coefficient *= iszero(bits[a][j]) ? (mesh.upper[k][j] - x[j]) / (mesh.upper[k][j] - mesh.lower[k][j]) : (x[j] - mesh.lower[k][j]) / (mesh.upper[k][j] - mesh.lower[k][j])
            end
            value -= raw[vertex] * coefficient
        end
        return value
    end
end
function _acq_box_field(mesh::CartesianAMRHierarchy{D}, raw, lo, hi) where {D}
    total = fill(zero(_CompositeQ), 3)
    covered = zero(_CompositeQ)
    volume = prod(d -> hi[d] - lo[d], 1:D)
    volume > 0 || throw(ArgumentError("composite field averaging box must have positive volume"))
    for k in _camr_query(mesh, lo, hi)
        l = ntuple(d -> max(lo[d], mesh.lower[k][d]), D)
        u = ntuple(d -> min(hi[d], mesh.upper[k][d]), D)
        v = prod(d -> u[d] - l[d], 1:D)
        covered += v
        point = ntuple(d -> (l[d] + u[d]) / 2, D)
        E = _acq_field(mesh, raw, k, point)
        for c = 1:3
            total[c] += v * E[c]
        end
    end
    covered == volume || throw(ArgumentError("composite field averaging box is not covered exactly once"))
    return ntuple(c -> total[c] / volume, 3)
end
function _acq_means(mesh::CartesianAMRHierarchy{D,T}, phi) where {D,T}
    raw = _acq_raw_potential(mesh, phi)
    out = ntuple(_ -> Vector{T}(undef, length(mesh.addresses)), 3)
    for k in eachindex(mesh.addresses)
        mean = _acq_field(mesh, raw, k, ntuple(d -> (mesh.lower[k][d] + mesh.upper[k][d]) / 2, D))
        for c = 1:3
            out[c][k] = _from_exact(T, mean[c])
            isfinite(out[c][k]) || throw(ArgumentError("composite mean electric field exceeds represented range"))
        end
    end
    return out
end
