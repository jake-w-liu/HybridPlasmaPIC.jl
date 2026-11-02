# Specific SBP-x / periodic-transverse operators. No advancing shock state.

"""
    MixedShockOperators(n, L)

Cached D2/3 CPU Float32/64 operators: SBP-(2,1) in x, registered
SpectralOperators Fourier derivatives transversely. All vector fields have three
components. Work/metadata are private, mutable and single-caller; use a separate
operator per concurrent task. Public operations stage outputs before publication.
Storage checks use Base.mightalias; callers must supply Arrays whose overlap it
reports. See `docs/MIXED_HALL_OPERATORS.md` for the qualified range/FFT contract.
"""
struct MixedShockOperators{D,T,F,P}
    n::NTuple{D,Int}
    L::NTuple{D,T}
    dx::NTuple{D,T}
    sbp::SBP1D{T}
    coordinates::NTuple{D,Vector{T}}
    transverse::F
    packed::P
    differentiated::P
    stage::NTuple{3,Array{T,D}}
    temp::Array{T,D}
    weighted::Array{T,D}
end

function _mixed_mode(::Type{T}) where {T}
    T in (Float32, Float64) || throw(ArgumentError("mixed operators require Float32 or Float64"))
    rounding(T) === RoundNearest && !get_zero_subnormals() || throw(
        ArgumentError("mixed operators require nearest rounding and gradual underflow"),
    )
    return nothing
end

function _mixed_count(n)
    count = 1
    for size in n
        count = Base.Checked.checked_mul(count, size)
    end
    return count
end

function MixedShockOperators(n::NTuple{D,Integer}, L::NTuple{D,T}) where {D,T<:Union{Float32,Float64}}
    D in (2, 3) || throw(ArgumentError("mixed operators require two or three spatial dimensions"))
    _mixed_mode(T)
    nn = ntuple(d -> _require_positive_intlike("n[$d]", n[d]), Val(D))
    nn[1] >= 3 || throw(ArgumentError("mixed SBP x needs at least three nodes"))
    lengths = ntuple(d -> _require_finite_positive_real("L[$d]", L[d], T), Val(D))
    count = _mixed_count(nn)
    sbp = SBP1D(nn[1], lengths[1])
    dx = ntuple(d -> d == 1 ? sbp.dx : lengths[d] / nn[d], Val(D))
    all(h -> isfinite(h) && h > zero(T), dx) ||
        throw(ArgumentError("mixed spacings must be positive and finite"))
    sbp.dx >= floatmin(T) && sbp.dx / T(2) >= floatmin(T) && isfinite(T(2) * sbp.dx) ||
        throw(ArgumentError("mixed SBP dx and half dx must be normal, and twice dx finite"))
    coordinates = ntuple(Val(D)) do d
        d == 1 ? collect(range(zero(T), lengths[d]; length = nn[d])) :
        [(i - 1) * lengths[d] / nn[d] for i = 1:nn[d]]
    end
    all(a -> all(isfinite, a), coordinates) ||
        throw(ArgumentError("mixed coordinates must be finite"))
    transverse = ntuple(i -> FourierDerivYWorkspace(count ÷ nn[i+1], nn[i+1], lengths[i+1]), Val(D - 1))
    packed = ntuple(i -> zeros(T, count ÷ nn[i+1], nn[i+1]), Val(D - 1))
    differentiated = map(a -> similar(a), packed)
    A() = zeros(T, nn)
    op = MixedShockOperators(nn, lengths, dx, sbp, coordinates, transverse,
        packed, differentiated, (A(), A(), A()), A(), A())
    _mixed_storage(op)
    return op
end

_mixed_operator_arrays(op) = (op.coordinates..., op.sbp.H,
    (w.ky for w in op.transverse)..., (w.cbuf for w in op.transverse)...,
    op.packed..., op.differentiated..., op.stage..., op.temp, op.weighted)

function _mixed_noalias(arrays, context)
    for i = 1:length(arrays), j = (i+1):length(arrays)
        (arrays[i] === arrays[j] || Base.mightalias(arrays[i], arrays[j])) &&
            throw(ArgumentError("$context storage roles must be independent"))
    end
    return nothing
end

function _mixed_array(a, op::MixedShockOperators{D,T}, context; finite = false) where {D,T}
    a isa Array{T,D} || throw(ArgumentError("$context requires ordinary $D-D Array{$T} storage"))
    size(a) == op.n || throw(DimensionMismatch("$context shape must be $(op.n)"))
    axes(a) == map(Base.OneTo, op.n) || throw(DimensionMismatch("$context requires one-based axes"))
    finite && _mixed_finite(a, context)
    return nothing
end

function _mixed_finite(a, context)
    all(isfinite, a) || throw(ArgumentError(
        "$context crosses a finite materialized mixed-operator boundary; later scaling is unsupported"))
    return nothing
end

function _mixed_storage(op::MixedShockOperators{D,T}) where {D,T}
    _mixed_mode(T)
    D in (2, 3) && op.n[1] >= 3 && all(>(0), op.n) ||
        throw(ArgumentError("invalid mixed geometry"))
    count = _mixed_count(op.n)
    all(x -> isfinite(x) && x > zero(T), op.L) || throw(ArgumentError("invalid mixed lengths"))
    expected_dx = ntuple(d -> op.L[d] / (d == 1 ? op.n[1] - 1 : op.n[d]), Val(D))
    all(h -> isfinite(h) && h > zero(T), op.dx) ||
        throw(ArgumentError("mixed spacings must be positive and finite"))
    op.dx == expected_dx && op.sbp.n == op.n[1] && op.sbp.dx == op.dx[1] ||
        throw(ArgumentError("mixed spacing cache does not match geometry"))
    dx = op.dx[1]
    dx >= floatmin(T) && dx / T(2) >= floatmin(T) && isfinite(T(2) * dx) && isfinite(one(T) / dx) ||
        throw(ArgumentError("invalid mixed SBP coefficients"))
    length(op.sbp.H) == op.n[1] || throw(DimensionMismatch("mixed SBP norm shape changed"))
    for i = 1:op.n[1]
        op.sbp.H[i] == (i == 1 || i == op.n[1] ? dx / T(2) : dx) ||
            throw(ArgumentError("mixed SBP norm cache changed"))
    end
    for d = 1:D
        length(op.coordinates[d]) == op.n[d] || throw(DimensionMismatch("mixed coordinate shape changed"))
        expected = d == 1 ? range(zero(T), op.L[d]; length = op.n[d]) : nothing
        for i = 1:op.n[d]
            x = d == 1 ? expected[i] : (i - 1) * op.L[d] / op.n[d]
            isfinite(x) && op.coordinates[d][i] == x || throw(ArgumentError("mixed coordinates changed"))
        end
    end
    length(op.transverse) == D - 1 && length(op.packed) == D - 1 && length(op.differentiated) == D - 1 ||
        throw(DimensionMismatch("mixed transverse workspace count changed"))
    for d = 2:D
        w = op.transverse[d-1]
        dims = (count ÷ op.n[d], op.n[d])
        (w.nx, w.ny) == dims && w.Ly == op.L[d] && size(w.cbuf) == dims &&
            size(op.packed[d-1]) == dims && size(op.differentiated[d-1]) == dims ||
            throw(DimensionMismatch("mixed transverse workspace shape/config changed"))
        length(w.ky) == op.n[d] || throw(DimensionMismatch("mixed Fourier coefficient count changed"))
        for m = 0:(op.n[d]-1)
            mp = m <= op.n[d] ÷ 2 ? m : m - op.n[d]
            k = iseven(op.n[d]) && m == op.n[d] ÷ 2 ? zero(T) : T(2π) * mp / op.L[d]
            isfinite(k) && w.ky[m+1] == k || throw(ArgumentError("mixed Fourier coefficients changed"))
        end
    end
    for a in (op.stage..., op.temp, op.weighted)
        _mixed_array(a, op, "mixed private field")
    end
    _mixed_noalias(_mixed_operator_arrays(op), "mixed operator")
    return nothing
end

function _mixed_preflight(op, outputs, inputs; extra = ())
    _mixed_storage(op)
    for a in outputs
        _mixed_array(a, op, "mixed output")
    end
    for a in inputs
        _mixed_array(a, op, "mixed input"; finite = true)
    end
    _mixed_noalias((outputs..., inputs..., extra..., _mixed_operator_arrays(op)...), "mixed call")
    return nothing
end

# Raw owners use only private validated roles. No public input is written.
function _mixed_derivative!(out, f, op::MixedShockOperators{D,T}, axis) where {D,T}
    if axis > D
        fill!(out, zero(T))
    elseif axis == 1
        for I in CartesianIndices(Base.tail(op.n))
            tail = Tuple(I)
            sbp_deriv!(view(out, :, tail...), view(f, :, tail...), op.sbp)
        end
    else
        pack, diff = op.packed[axis-1], op.differentiated[axis-1]
        for (row, I) in enumerate(CartesianIndices(ntuple(d -> d == axis ? 1 : op.n[d], Val(D))))
            for j = 1:op.n[axis]
                K = CartesianIndex(ntuple(d -> d == axis ? j : I[d], Val(D)))
                pack[row, j] = f[K]
            end
        end
        fourier_deriv_y!(diff, pack, op.transverse[axis-1])
        for (row, I) in enumerate(CartesianIndices(ntuple(d -> d == axis ? 1 : op.n[d], Val(D))))
            for j = 1:op.n[axis]
                K = CartesianIndex(ntuple(d -> d == axis ? j : I[d], Val(D)))
                out[K] = diff[row, j]
            end
        end
    end
    _mixed_finite(out, "mixed derivative axis $axis")
    return out
end

function _mixed_transpose!(out, f, op::MixedShockOperators{D,T}, axis) where {D,T}
    if axis == 1
        for I in CartesianIndices(f)
            op.weighted[I] = (I[1] == 1 || I[1] == op.n[1]) ? T(2) * f[I] : f[I]
        end
        _mixed_finite(op.weighted, "mixed transpose Q inverse input")
        _mixed_derivative!(out, op.weighted, op, 1)
        for I in CartesianIndices(out)
            q = I[1] == 1 || I[1] == op.n[1] ? T(1) / T(2) : one(T)
            boundary = I[1] == 1 ? -one(T) : I[1] == op.n[1] ? one(T) : zero(T)
            out[I] = -q * out[I] + boundary / op.dx[1] * op.weighted[I]
        end
    else
        _mixed_derivative!(out, f, op, axis)
        out .*= -one(T)
    end
    _mixed_finite(out, "mixed transpose axis $axis")
    return out
end

function _mixed_curl!(out, vector, op)
    for (c, plus, plus_axis, minus, minus_axis) in _MIXED_CURL_CHANNELS
        _mixed_derivative!(out[c], vector[plus], op, plus_axis)
        _mixed_derivative!(op.temp, vector[minus], op, minus_axis)
        out[c] .-= op.temp
        _mixed_finite(out[c], "mixed curl component $c")
    end
    return out
end

function _mixed_divergence!(out, vector, op::MixedShockOperators{D}) where {D}
    _mixed_derivative!(out, vector[1], op, 1)
    for d = 2:D
        _mixed_derivative!(op.temp, vector[d], op, d)
        out .+= op.temp
        _mixed_finite(out, "mixed divergence")
    end
    return out
end

"Stage the D-component mixed gradient before copying any public output."
function mixed_gradient!(out::NTuple{D,Array{T,D}}, f::Array{T,D}, op::MixedShockOperators{D,T}) where {D,T}
    _mixed_preflight(op, out, (f,))
    for d = 1:D
        _mixed_derivative!(op.stage[d], f, op, d)
    end
    for d = 1:D
        copyto!(out[d], op.stage[d])
    end
    return out
end

"Stage all three curl components; absent spatial derivatives vanish in D2."
function mixed_curl!(out::NTuple{3,Array{T,D}}, vector::NTuple{3,Array{T,D}}, op::MixedShockOperators{D,T}) where {D,T}
    _mixed_preflight(op, out, vector)
    _mixed_curl!(op.stage, vector, op)
    for c = 1:3
        copyto!(out[c], op.stage[c])
    end
    return out
end

"Stage mixed divergence; all three supplied components must be finite, including inactive ones."
function mixed_divergence!(out::Array{T,D}, vector::NTuple{3,Array{T,D}}, op::MixedShockOperators{D,T}) where {D,T}
    _mixed_preflight(op, (out,), vector)
    _mixed_divergence!(op.stage[1], vector, op)
    copyto!(out, op.stage[1])
    return out
end

"Stage the represented ordinary Euclidean derivative transpose (not the weighted adjoint)."
function mixed_derivative_transpose!(out::Array{T,D}, f::Array{T,D}, op::MixedShockOperators{D,T}, axis::Integer) where {D,T}
    1 <= axis <= D || throw(ArgumentError("mixed transpose axis must be between 1 and $D"))
    _mixed_preflight(op, (out,), (f,))
    _mixed_transpose!(op.stage[1], f, op, axis)
    copyto!(out, op.stage[1])
    return out
end

"Private derived fields for mixed_hall_rhs!; not physical state or checkpoint data."
struct MixedHallWorkspace{D,T}
    n::NTuple{D,Int}
    pe::Array{T,D}
    gradp::NTuple{D,Array{T,D}}
    ninv::Array{T,D}
    E::NTuple{3,Array{T,D}}
    J::NTuple{3,Array{T,D}}
    lapJ::NTuple{3,Array{T,D}}
    rate::NTuple{3,Array{T,D}}
    floor_count::Base.RefValue{Int}
end

function MixedHallWorkspace(op::MixedShockOperators{D,T}) where {D,T}
    _mixed_storage(op)
    A() = zeros(T, op.n)
    V() = (A(), A(), A())
    return MixedHallWorkspace(op.n, A(), ntuple(_ -> A(), Val(D)), A(), V(), V(), V(), V(), Ref(0))
end

_mixed_hall_arrays(w) = (w.pe, w.gradp..., w.ninv, w.E..., w.J..., w.lapJ..., w.rate...)

"""
    mixed_hall_rhs!(out, B, n, ui, physics, op, work)

Raw full Hall/isothermal induction, with strict density and η=ηH=de2=0. Reuses
canonical pressure, reciprocal, completed Ohm rows and Faraday sign. Original-
product differentiated-source recovery is not supplied: source-risk events
reject before publication. Return private borrowed E/J/pressure/reciprocal/rate
and an exact stored-source face-work budget, valid until workspace reuse.
No SAT, projection or time step is implicitly performed.
"""
function mixed_hall_rhs!(out::NTuple{3,Array{T,D}}, B::NTuple{3,Array{T,D}}, n::Array{T,D},
    ui::NTuple{3,Array{T,D}}, physics::HybridModel, op::MixedShockOperators{D,T}, work::MixedHallWorkspace{D,T}) where {D,T}
    physics.closure isa IsothermalElectrons{T} ||
        throw(ArgumentError("mixed Hall requires isothermal electrons at grid precision"))
    _require_finite_nonnegative_real("Te", physics.closure.Te, T)
    all(_particle_stored_zero, (physics.η, physics.ηH, physics.de2)) ||
        throw(ArgumentError("mixed Hall requires original η=ηH=de2=0"))
    floor = _require_finite_positive_real("nfloor", physics.nfloor, T)
    work.n == op.n || throw(DimensionMismatch("mixed Hall workspace geometry changed"))
    for a in _mixed_hall_arrays(work)
        _mixed_array(a, op, "mixed Hall workspace")
    end
    _mixed_preflight(op, out, (B..., n, ui...); extra = _mixed_hall_arrays(work))
    all(x -> x >= floor, n) || throw(ArgumentError("mixed Hall density must be at least nfloor; clipping is unsupported"))
    electron_pressure!(work.pe, n, physics.closure)
    _mixed_finite(work.pe, "mixed Hall pressure")
    for d = 1:D
        _mixed_derivative!(work.gradp[d], work.pe, op, d)
    end
    _ohm_ninv!(work.ninv, n, floor, work.floor_count)
    iszero(work.floor_count[]) || throw(ArgumentError("mixed Hall density clipping is unsupported"))
    _mixed_finite(work.ninv, "mixed Hall reciprocal density")
    _mixed_curl!(work.J, B, op)
    for a in work.lapJ
        fill!(a, zero(T))
    end
    risk = _ohm_field_rows!(work.E, ui, B, work.J, work.lapJ, work.gradp, work.ninv, zero(T), zero(T))
    foreach(a -> _mixed_finite(a, "mixed Hall electric field"), work.E)
    risk && throw(ArgumentError(
        "mixed Hall original-product Faraday range recovery is unsupported; a finite local Ohm field does not settle its differentiated source"))
    _mixed_curl!(work.rate, work.E, op)
    _faraday_negate!(work.rate)
    foreach(a -> _mixed_finite(a, "mixed Hall magnetic rate"), work.rate)
    budget = _mixed_magnetic_budget(B, work.E, work.J, work.rate, op)
    report = (E = work.E, J = work.J, pe = work.pe, gradp = work.gradp,
        ninv = work.ninv, rate = work.rate, source_risk = risk, budget = budget)
    for c = 1:3
        copyto!(out[c], work.rate[c])
    end
    return report
end

function _mixed_qdot(a, b, op)
    value = zero(Rational{BigInt})
    for c = 1:3, I in CartesianIndices(a[c])
        weight = I[1] == 1 || I[1] == op.n[1] ? BigInt(1) // BigInt(2) : one(Rational{BigInt})
        value += weight * _exact_rational(a[c][I]) * _exact_rational(b[c][I])
    end
    return value
end

function _mixed_magnetic_budget(B, E, J, rate, op::MixedShockOperators{D}) where {D}
    volume = prod(map(_exact_rational, op.dx))
    area = prod(map(_exact_rational, Base.tail(op.dx)))
    poynting = zero(Rational{BigInt})
    for side in (1, op.n[1]), I in CartesianIndices(Base.tail(op.n))
        K = CartesianIndex((side, Tuple(I)...))
        normal = side == 1 ? -1 : 1
        poynting += normal * area * (_exact_rational(E[2][K]) * _exact_rational(B[3][K]) -
            _exact_rational(E[3][K]) * _exact_rational(B[2][K]))
    end
    power = volume * _mixed_qdot(B, rate, op)
    current_work = volume * _mixed_qdot(J, E, op)
    return (magnetic_power = power, current_work = current_work, outward_poynting = poynting,
        identity_residual = power + current_work + poynting)
end

"Exact rational stored-source mixed magnetic power, current work and two-face Poynting ledger."
function mixed_magnetic_rate_budget(B::NTuple{3,Array{T,D}}, E::NTuple{3,Array{T,D}},
    J::NTuple{3,Array{T,D}}, rate::NTuple{3,Array{T,D}}, op::MixedShockOperators{D,T}) where {D,T}
    _mixed_preflight(op, (), (B..., E..., J..., rate...))
    return _mixed_magnetic_budget(B, E, J, rate, op)
end
