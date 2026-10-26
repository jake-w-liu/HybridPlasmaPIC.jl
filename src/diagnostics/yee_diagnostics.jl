# Yee field diagnostics. Float32 fields are evaluated in Float64: every binary32
# square, product, interpolant and grid factor below is a normal Float64 value,
# so that path has no range loss before its final rounding. Float64 fields use
# native arithmetic with cheap range checks; detected overflow or a harmful
# underflow (a nonzero intermediate leaving the normal range before later
# scaling) selects the exact rational expression, rounded once.

@inline _yee_wide(::Type{Float32}) = Float64
@inline _yee_wide(::Type{Float64}) = Float64

# Nonnegative pairwise sum with its underflow certificate: n squares below the
# normal range contribute at most n*2^-1075 < u*total once total >= n*floatmin.
@inline function _yee_square_sum(::Type{W}, arrays) where {W}
    total = zero(W)
    count = 0
    for a in arrays
        total += sum(x -> W(x)^2, a)
        count += length(a)
    end
    return total, count
end

@inline _yee_tiny_sum(total, count, arrays) =
    iszero(total) ? !all(a -> all(iszero, a), arrays) : total < count * floatmin(total)

"""
    em_field_energy(model::MaxwellVacuum)

Physical field-energy quadrature at synchronized time, `Σ(E²+c²B²)*prod(dx)/2`.
PEC endpoint trapezoid weights multiply constrained zero entries, so all nonzero
degrees of freedom have unit weights. This energy has second-order temporal
oscillations; it is not the exactly conserved modified leapfrog energy.
Sums use Julia's pairwise summation (Float64 for Float32 fields); range loss is
completed exactly, and an exact value above the floating-point range is `Inf`.
"""
function em_field_energy(model::MaxwellVacuum{D,T}) where {D,T}
    _yee_validate_fields(model.fields, model.g)
    _yee_validate_pec(model.fields, model.g)
    E, B = model.fields.E, model.fields.B
    W = _yee_wide(T)
    electric, electric_count = _yee_square_sum(W, E)
    magnetic, magnetic_count = _yee_square_sum(W, B)
    speed = W(model.c)
    squared = speed * speed
    weighted = squared * magnetic
    value = W(0.5) * (electric + weighted)
    # Partial products must stay zero or normal before each spacing factor.
    lost = false
    for d = 1:D
        lost |= !iszero(value) & !_yee_normal(value)
        value *= W(model.g.dx[d])
    end
    T === Float32 && return T(value)
    lost |= !isfinite(value)
    lost |= !iszero(magnetic) & (!_yee_normal(squared) | !_yee_normal(weighted))
    lost |= _yee_tiny_sum(electric, electric_count, E) || _yee_tiny_sum(magnetic, magnetic_count, B)
    lost || return value
    return _em_field_energy_exact(E, B, model.c, model.g.dx)
end

# ------------------------------------------------------------ Poynting vector

# Sum of the entries at I and at the upper neighbor U along each listed axis.
@inline _yee_corner_sum(::Type{W}, a, I, U, ::Tuple{}) where {W} = W(a[I...]), 1
@inline function _yee_corner_sum(::Type{W}, a, I, U, axes::Tuple{Int}) where {W}
    d = axes[1]
    return W(a[I...]) + W(a[Base.setindex(I, U[d], d)...]), 2
end
@inline function _yee_corner_sum(::Type{W}, a, I, U, axes::Tuple{Int,Int}) where {W}
    d, e = axes
    Id = Base.setindex(I, U[d], d)
    left = W(a[I...]) + W(a[Id...])
    right = W(a[Base.setindex(I, U[e], e)...]) + W(a[Base.setindex(Id, U[e], e)...])
    return left + right, 4
end

# Exact interpolant: the same corners summed in rational arithmetic.
function _yee_interpolant_exact(a, I, U, axes)
    total, count = _yee_corner_sum(Rational{BigInt}, a, I, U, axes)
    return total / count
end

# Linear interpolation: the corner sum times an exact power-of-two weight. A
# checked weight is exact whenever the scaled interpolant stays normal.
@inline function _yee_interpolant(::Type{W}, a, I, U, axes, checked::Bool) where {W}
    total, count = _yee_corner_sum(W, a, I, U, axes)
    lost = checked & !iszero(total) & (abs(total) < count * floatmin(W))
    return total * W(inv(count)), lost
end

@inline _yee_product_underflow(x, y, product) =
    !iszero(x) & !iszero(y) & (abs(product) < floatmin(product))

# c²(xj*yk - xk*yj) for one Poynting component, with its checked product range.
@inline function _yee_cross(xj, yk, xk, yj, checked::Bool)
    first, second = xj * yk, xk * yj
    lost =
        checked & (_yee_product_underflow(xj, yk, first) | _yee_product_underflow(xk, yj, second))
    return first - second, lost
end

@inline function _yee_nodal_axes(electric::Bool, ::Val{c}, ::Val{D}) where {c,D}
    electric || return c <= D ? (c,) : ()
    j, k = mod1(c + 1, 3), mod1(c + 2, 3)
    return j <= D ? (k <= D ? (j, k) : (j,)) : (k <= D ? (k,) : ())
end

@inline _yee_upper(g::YeeGrid{D}, I::NTuple{D,Int}) where {D} =
    ntuple(d -> g.periodic[d] && I[d] == g.n[d] ? 1 : I[d] + 1, D)

function _yee_poynting_cell_exact(model::MaxwellVacuum{D}, I, U) where {D}
    E = ntuple(
        c -> _yee_interpolant_exact(model.fields.E[c], I, U, _yee_nodal_axes(true, Val(c), Val(D))),
        3,
    )
    B = ntuple(
        c -> _yee_interpolant_exact(
            model.fields.B[c],
            I,
            U,
            _yee_nodal_axes(false, Val(c), Val(D)),
        ),
        3,
    )
    squared = _exact_rational(model.c)^2
    return ntuple(3) do c
        j, k = mod1(c + 1, 3), mod1(c + 2, 3)
        squared * (E[j] * B[k] - E[k] * B[j])
    end
end

function _yee_poynting_cells!(result, model::MaxwellVacuum{D,T}) where {D,T}
    E, B, g = model.fields.E, model.fields.B, model.g
    W = _yee_wide(T)
    checked = W === T
    speed = W(model.c)
    squared = speed * speed
    lost_squared = checked & !_yee_normal(squared)
    for J in CartesianIndices(g.n)
        I = Tuple(J)
        U = _yee_upper(g, I)
        e1, l1 = _yee_interpolant(W, E[1], I, U, _yee_nodal_axes(true, Val(1), Val(D)), checked)
        e2, l2 = _yee_interpolant(W, E[2], I, U, _yee_nodal_axes(true, Val(2), Val(D)), checked)
        e3, l3 = _yee_interpolant(W, E[3], I, U, _yee_nodal_axes(true, Val(3), Val(D)), checked)
        b1, m1 = _yee_interpolant(W, B[1], I, U, _yee_nodal_axes(false, Val(1), Val(D)), checked)
        b2, m2 = _yee_interpolant(W, B[2], I, U, _yee_nodal_axes(false, Val(2), Val(D)), checked)
        b3, m3 = _yee_interpolant(W, B[3], I, U, _yee_nodal_axes(false, Val(3), Val(D)), checked)
        x1, n1 = _yee_cross(e2, b3, e3, b2, checked)
        x2, n2 = _yee_cross(e3, b1, e1, b3, checked)
        x3, n3 = _yee_cross(e1, b2, e2, b1, checked)
        s1, s2, s3 = squared * x1, squared * x2, squared * x3
        lost = lost_squared | l1 | l2 | l3 | m1 | m2 | m3 | n1 | n2 | n3
        lost |= checked & !(isfinite(s1) & isfinite(s2) & isfinite(s3))
        if lost
            exact = _yee_poynting_cell_exact(model, I, U)
            for c = 1:3
                result[c][J] = _from_exact(T, exact[c])
            end
        else
            result[1][J], result[2][J], result[3][J] = T(s1), T(s2), T(s3)
        end
    end
    return result
end

"""
    yee_poynting(model)

Return three cell-centered arrays of instantaneous `c² E×B`, using tensor-product
linear interpolation of synchronized Yee fields. Normwise contract: each entry
agrees with the complete expression from the original samples to a few `eps(T)`
times `c²(|E_j||B_k| + |E_k||B_j|)` (interpolants of absolute samples) plus a few
subnormal quanta; overflowing or underflowing intermediates are evaluated
exactly. This local visualization diagnostic is not a discrete face-flux
balance for the leapfrog modified energy. Closed PEC walls have zero normal
physical Poynting flux; paired periodic faces cancel. Accordingly this vacuum
model exchanges no energy with its boundaries.
"""
function yee_poynting(model::MaxwellVacuum{D,T}) where {D,T}
    _yee_validate_fields(model.fields, model.g)
    _yee_validate_pec(model.fields, model.g)
    return _yee_poynting_cells!(ntuple(_ -> zeros(T, model.g.n), 3), model)
end

# --------------------------------------------------------- face Poynting flux
#
# Face samples: tangential E sits on the face node (PEC: wall node, zero by
# constraint; periodic: node 1 for both faces). Tangential B uses the even
# extension at a PEC wall (adjacent center) or the average across a periodic
# seam. Tangential nodal axes are averaged as in the cell interpolation.

@inline function _yee_face_axes(tangent::Int, normal::Int, across::Bool)
    tangent == 0 && return across ? (normal,) : ()
    return across ? (tangent, normal) : (tangent,)
end

@inline function _yee_face_interpolant(
    ::Type{W},
    a,
    I,
    U,
    tangent,
    normal,
    across,
    checked,
) where {W}
    if tangent == 0
        return across ? _yee_interpolant(W, a, I, U, (normal,), checked) :
               _yee_interpolant(W, a, I, U, (), checked)
    end
    return across ? _yee_interpolant(W, a, I, U, (tangent, normal), checked) :
           _yee_interpolant(W, a, I, U, (tangent,), checked)
end

# Index tuples of one tangential face cell: E on the face node, B on the
# adjacent (or seam-straddling) center, and upper neighbors per axis.
@inline function _yee_face_indices(g::YeeGrid{D}, cell, normal::Int, upper::Bool) where {D}
    I = Tuple(cell)
    n, periodic = g.n[normal], g.periodic[normal]
    U = ntuple(d -> d == normal ? 1 : (g.periodic[d] && I[d] == g.n[d] ? 1 : I[d] + 1), D)
    node = periodic || !upper ? 1 : n + 1
    center = periodic || upper ? n : 1
    return Base.setindex(I, node, normal), Base.setindex(I, center, normal), U
end

function _yee_face_exact(model::MaxwellVacuum{D,T}, normal::Int, upper::Bool) where {D,T}
    g = model.g
    j, k = mod1(normal + 1, 3), mod1(normal + 2, 3)
    across = g.periodic[normal]
    total = zero(Rational{BigInt})
    for cell in CartesianIndices(ntuple(d -> d == normal ? 1 : g.n[d], D))
        IE, IB, U = _yee_face_indices(g, cell, normal, upper)
        ej = _yee_interpolant_exact(
            model.fields.E[j],
            IE,
            U,
            _yee_face_axes(k <= D ? k : 0, normal, false),
        )
        ek = _yee_interpolant_exact(
            model.fields.E[k],
            IE,
            U,
            _yee_face_axes(j <= D ? j : 0, normal, false),
        )
        bj = _yee_interpolant_exact(
            model.fields.B[j],
            IB,
            U,
            _yee_face_axes(j <= D ? j : 0, normal, across),
        )
        bk = _yee_interpolant_exact(
            model.fields.B[k],
            IB,
            U,
            _yee_face_axes(k <= D ? k : 0, normal, across),
        )
        total += ej * bk - ek * bj
    end
    total *= _exact_rational(model.c)^2
    for d = 1:D
        d == normal || (total *= _exact_rational(g.dx[d]))
    end
    return _from_exact(T, upper ? total : -total)
end

function _yee_face_flux(model::MaxwellVacuum{D,T}, normal::Int, upper::Bool) where {D,T}
    E, B, g = model.fields.E, model.fields.B, model.g
    W = _yee_wide(T)
    checked = W === T
    j, k = mod1(normal + 1, 3), mod1(normal + 2, 3)
    tj, tk = j <= D ? j : 0, k <= D ? k : 0
    across = g.periodic[normal]
    cells = CartesianIndices(ntuple(d -> d == normal ? 1 : g.n[d], D))
    term = index -> begin
        IE, IB, U = _yee_face_indices(g, cells[index], normal, upper)
        ej, lj = _yee_face_interpolant(W, E[j], IE, U, tk, normal, false, checked)
        ek, lk = _yee_face_interpolant(W, E[k], IE, U, tj, normal, false, checked)
        bj, mj = _yee_face_interpolant(W, B[j], IB, U, tj, normal, across, checked)
        bk, mk = _yee_face_interpolant(W, B[k], IB, U, tk, normal, across, checked)
        local cross, cross_lost = _yee_cross(ej, bk, ek, bj, checked)
        (cross, cross_lost | lj | lk | mj | mk)
    end
    # Pairwise reduction over the face cells.
    total, lost = mapreduce(term, (x, y) -> (x[1] + y[1], x[2] | y[2]), 1:length(cells))
    speed = W(model.c)
    squared = speed * speed
    # Each partial product of total*c²*area must be zero or normal before the
    # next factor; only the final rounding may underflow.
    lost |= checked & (!_yee_normal(squared) | (!iszero(total) & !_yee_normal(total)))
    value = total * squared
    for d = 1:D
        d == normal && continue
        lost |= checked & !iszero(value) & !_yee_normal(value)
        value *= W(g.dx[d])
    end
    lost |= checked & !isfinite(value)
    lost && return _yee_face_exact(model, normal, upper)
    return T(upper ? value : -value)
end

"""
    yee_poynting_flux(model)

Instantaneous outward physical Poynting flux integrated over each pair of domain
faces, as `((lower₁, upper₁), ...)`. Values use tangential cell-center quadrature
and centered field interpolation. PEC faces use the even tangential-B extension
and zero tangential E; periodic faces use the same wrapped interpolation with
opposite outward normals. Individual periodic-face fluxes can be nonzero, while
their sum cancels exactly. Face sums use pairwise summation; the normwise
contract of `yee_poynting` applies to the complete face integral, with scale
`c²·area·Σ(|E_j||B_k| + |E_k||B_j|)`. This is not a modified-leapfrog-energy
balance diagnostic.
"""
function yee_poynting_flux(model::MaxwellVacuum{D,T}) where {D,T}
    _yee_validate_fields(model.fields, model.g)
    _yee_validate_pec(model.fields, model.g)
    return ntuple(D) do normal
        ntuple(side -> _yee_face_flux(model, normal, side == 2), 2)
    end
end

# ------------------------------------------------ prescribed-current diagnostics

function _yee_source_work_exact(old, new, impulse, dx::T) where {T}
    total = zero(Rational{BigInt})
    for c = 1:3, i in eachindex(old[c])
        total +=
            _exact_rational(impulse[c][i]) *
            (_exact_rational(old[c][i]) + _exact_rational(new[c][i]))
    end
    return _from_exact(T, _exact_rational(dx) * total / 2)
end

# Whether some product impulse*(old+new) has two nonzero factors.
_yee_any_source_product(impulse, old, new) =
    any(c -> any(i -> !iszero(impulse[c][i]) && old[c][i] != -new[c][i], eachindex(old[c])), 1:3)

"""
    yee_source_work(Eold, Enew, impulse, g::YeeGrid{1})

Return `dx/2 * sum(impulse*(Eold+Enew))` in epsilon0=1 units with unit transverse
area. Positive work transfers energy out of the fields. Supply the electric
arrays immediately before and after the sourced step and its integrated current;
no history is retained by the model. Inputs are read-only, finite, correctly
shaped and consistent with homogeneous tangential PEC constraints.

The sum is pairwise (Float64 for Float32 data). Normwise contract: the result
agrees with the exact value to a few `eps(T)` per summed term times
`dx/2 * sum(|impulse|*(|Eold|+|Enew|))`; overflow, or underflow that is not
negligible at that scale, is evaluated exactly, and an exact value outside the
range is signed infinity. Exact-arithmetic B-half/E-full/B-half evolution obeys
`Hnew-Hold = -work` for `yee_modified_energy` at the same fixed dt;
floating-point staging introduces roundoff in that identity.
"""
function yee_source_work(old, new, impulse, g::YeeGrid{1,T}) where {T}
    _yee_validate_electric_data(old, g)
    _yee_validate_electric_data(new, g)
    _yee_validate_impulse(impulse, g)
    W = _yee_wide(T)
    total = zero(W)
    magnitude = zero(W)
    count = 0
    for c = 1:3
        p, a, b = impulse[c], old[c], new[c]
        total += sum(i -> W(p[i]) * (W(a[i]) + W(b[i])), eachindex(p))
        magnitude += sum(i -> abs(W(p[i])) * (abs(W(a[i])) + abs(W(b[i]))), eachindex(p))
        count += length(p)
    end
    value = (total * W(g.dx[1])) * W(0.5)
    T === Float32 && return T(value)
    # n underflowed products contribute at most n*2^-1075 < u*magnitude once
    # magnitude >= n*floatmin, independent of the later spacing scale.
    lost = !isfinite(value) | !isfinite(magnitude)
    lost = lost || (magnitude < count * floatmin(W) && _yee_any_source_product(impulse, old, new))
    lost || return value
    return _yee_source_work_exact(old, new, impulse, g.dx[1])
end

function _yee_modified_energy_exact(model::MaxwellVacuum{1,T}, dt::T) where {T}
    e2 = sum(_exact_rational(x)^2 for a in model.fields.E for x in a)
    b2 = sum(_exact_rational(x)^2 for a in model.fields.B for x in a)
    curl2 = zero(Rational{BigInt})
    for c = 1:3, I in CartesianIndices(model.fields.B[c])
        curl2 += _yee_curl_exact(model.fields.E, I, c, model.g, true)^2
    end
    result = e2 + _exact_rational(model.c)^2 * (b2 - _exact_rational(dt)^2 * curl2 / 4)
    return _from_exact(T, _exact_rational(model.g.dx[1]) * result / 2)
end

# Σ((dt*c/(2dx)) * ΔE)² over the 1D curl(E) faces. Since dt*c < dx, each scaled
# sample is below |ΔE|/2, so its square is dominated by the E² sum.
function _yee_scaled_curl_square_sum(factor::W, E, g::YeeGrid{1}) where {W}
    total = zero(W)
    count = 0
    for c = 2:3
        a = E[c]
        n = length(a)
        faces = g.periodic[1] ? n : n - 1
        total += sum(i -> (factor * (W(a[i == n ? 1 : i + 1]) - W(a[i])))^2, 1:faces)
        count += faces
    end
    return total, count
end

"""
    yee_modified_energy(model::MaxwellVacuum{1}, dt)

Return the fixed-step Yee quadratic invariant
`dx/2 * (sum(E²) + c²*sum(B²) - dt²*c²/4*sum(curl(E)²))`.
The supplied dt must be positive and strictly below this model's CFL bound.
The form is positive definite at that bound's strict interior and is conserved
by exact-arithmetic vacuum stepping at this fixed dt; prescribed current changes
it by `-yee_source_work(...)`. Changing dt changes the quadratic form itself.

This is distinct from synchronized physical `em_field_energy` and physical
Poynting flux. Sums are pairwise (Float64 for Float32 fields). Normwise contract:
the result agrees with the exact value to a few `eps(T)` per summed term times
the sum of the three nonnegative contributions; overflow or non-negligible
underflow is evaluated exactly, and an exact value outside the range is `Inf`.
Inputs are read-only and validated, and no prior state is retained.
"""
function yee_modified_energy(model::MaxwellVacuum{1,T}, dt::Real) where {T}
    dtT = _require_finite_positive_real("dt", dt, T)
    dtT < maxwell_cfl(model) ||
        throw(ArgumentError("dt must be strictly below the Maxwell CFL limit"))
    _yee_validate_fields(model.fields, model.g)
    _yee_validate_pec(model.fields, model.g)
    E, B, g = model.fields.E, model.fields.B, model.g
    W = _yee_wide(T)
    electric, electric_count = _yee_square_sum(W, E)
    magnetic, magnetic_count = _yee_square_sum(W, B)
    speed = W(model.c)
    squared = speed * speed
    weighted = squared * magnetic
    # An underflowed factor or scaled sample is negligible against the E² sum.
    factor = ((W(dtT) * speed) * W(0.5)) / W(g.dx[1])
    curl, _ = _yee_scaled_curl_square_sum(factor, E, g)
    value = (electric + weighted) - curl
    lost = !iszero(value) & !_yee_normal(value)
    value = (value * W(g.dx[1])) * W(0.5)
    T === Float32 && return T(value)
    lost |= !isfinite(value)
    lost |= !iszero(magnetic) & (!_yee_normal(squared) | !_yee_normal(weighted))
    lost |= _yee_tiny_sum(electric, electric_count, E) || _yee_tiny_sum(magnetic, magnetic_count, B)
    lost || return value
    return _yee_modified_energy_exact(model, dtT)
end
