# full_pic.jl — full electromagnetic PIC models, SEPARATE from the hybrid loop.
# `EMPIC1D` is the specialized 1D3V leapfrog/subcycling solver below; `EMPIC`
# is the dimension-parametric periodic 1D/2D/3D solver later in this file. Both
# use mobile kinetic electrons (q=−1, m=1). The positive species can be EITHER a
# uniform immobile ion background (+n0, the default) OR a second mobile kinetic
# ion species (q=+1, mass mi, default mi=1836·me). ε0=1, so the electron plasma
# frequency is ω_pe = √(n0) and ω_pe = 1 at n0 = 1. Speed of light c is a free
# parameter (default 5).
#
# Fields live on the collocated periodic Fourier grid; spatial derivatives are
# spectral (∂x ↔ ik with the Nyquist mode zeroed, supplied by FourierGrid{1}).
#
# Field components carried:
#   transverse EM  : Ey, Bz   with   ∂tEy = −c²∂xBz − Jy/ε0 ,  ∂tBz = −∂xEy
#   longitudinal   : Ex        with   ∂tEx = −Jx/ε0           (Ampère, ε0=1)
# (Bx, By, Ez are identically zero for a 1D problem with this polarization split
# and are never stored.)
#
# ---- Time levels (Yee-style leapfrog; E and B staggered by dt/2) -------------
#   particle position x : integer  n        (→ n+1)
#   particle velocity v : half     n−1/2    (→ n+1/2)
#   electric field  Ex,Ey : integer n       (→ n+1)
#   magnetic field  Bz    : half   n−1/2    (→ n+1/2)
# For changing timesteps, the half levels lie dt_previous/2 behind and dt/2
# ahead of the integer time. Velocity and B kicks span their sum, while position
# drifts, continuity, and the electric-field update use the new dt.
#
# A single field/electron substep advances all of these by one dt (the cycle is
# documented in _substep_em!):
#   1. B half-advance to n:  Bz^n = Bz^{n−1/2} − (dt/2) ∂xEy^n  (stored separately
#      as a time-centered B for the particle push; the canonical Bz stays at the
#      half level and is advanced a full dt below).
#   2. Gather (Ex^n,Ey^n,Bz^n) to particles; Boris push v^{n−1/2}→v^{n+1/2} and
#      drift x^n→x^{n+1}; record the midpoint x^{n+1/2}; wrap to the periodic box.
#   3. Charge-conserving longitudinal current Jx^{n+1/2} from ρ^n (x^n) and ρ^{n+1}
#      (x^{n+1}) via the exact spectral continuity solve  ik·Ĵx = −(ρ̂^{n+1}−ρ̂^n)/dt
#      (so ∂tρ + ∂xJx = 0 to roundoff). Transverse Jy^{n+1/2} = Σ q w v_y S at the
#      midpoint positions x^{n+1/2} (CIC). Both species contribute when mobile.
#   4. Bz: full leapfrog advance  Bz^{n+1/2} = Bz^{n−1/2} − dt ∂xEy^n.
#   5. Ey: Ey^{n+1} = Ey^n − dt(c²∂xBz^{n+1/2} + Jy^{n+1/2}).
#      Ex: Ex^{n+1} = Ex^n − dt·Jx^{n+1/2}.
#
# 2nd-order accuracy (verified 2026-07, deep-debug): the seeded initial state is physical (v^0 and
# Bz^0), but the leapfrog carries v and Bz at the HALF level, so the FIRST step primes both once
# (see the priming block in step_empic!): (a) velocity v^{-1/2} = v^0 − h·a^0 for BOTH species
# (h = dt_e/2 electrons, dt_ion/2 ions — electrons alone is not enough when ions are mobile), and
# (b) field Bz^{-1/2} = Bz^0 + (dt_e/2)∂xEy^0. Without either, the integer-level fields (Ex, Ey)
# are only 1st-order in dt; with both, cold-deterministic self-convergence gives rate → 2.0 for
# the longitudinal Langmuir Ex and the transverse Ey. (The stored Bz stays at the half level, so a
# diagnostic reading es.Bz as "Bz at the integer time" sees the inherent O(dt) Yee stagger — use
# the time-centered Bz from step 1 if an integer-level Bz is needed.)
#
# ---- Mobile ions -------------------------------------------------------------
# When `mobile=true`, a second kinetic species (the ions passed to init/step)
# is pushed alongside the electrons. Ions are heavy (mi≫me) so they respond
# slowly; both species deposit charge ρ = ρ_ion + ρ_electron = (q_i n_i + q_e n_e)
# and current J = J_ion + J_electron. The longitudinal Jx is still built from the
# TOTAL ρ snapshots via the exact continuity solve, so charge conservation holds
# for the combined system to roundoff.
#
# ---- Electron subcycling -----------------------------------------------------
# `n_sub` (default 1) splits each ion step `dt` into n_sub electron substeps of
# dt_e = dt/n_sub. The fields and electrons are advanced on the FAST substep
# (electrons resolve ω_pe and the EM CFL). The heavy ion VELOCITY is kicked once
# per ion step (over dt_ion = n_sub·dt_e) on the FIRST substep, but its POSITION
# is drifted dt_e on EVERY substep with the post-kick velocity, so the ion's
# charge/current is deposited at the time-consistent position on each substep.
# n_sub=1 reproduces the single-rate scheme exactly (and bit-for-bit reproduces
# the legacy immobile-ion run). This drift-splitting makes the mobile scheme
# 2nd-order in dt for ANY n_sub (verified 2026-07, deep-debug: cold self-
# convergence rate → 2.0 for n_sub=1..4, both light and heavy ions) while keeping
# exact Esirkepov charge conservation (per-substep residual ~1e-13). [Kicking on
# the middle substep with the whole-dt_ion drift lumped the ion current into one
# substep and left the others time-inconsistent → the earlier 1st-order behavior.]
#
# Verified oracles (test_empic1d.jl):
#   (1) transverse EM-wave dispersion ω² = ω_pe² + c²k² (≥2 k values, <3%);
#   (2) charge conservation max|ρ^{n+1}−ρ^n + dt ∂xJx|/scale < 1e-10;
#   (3) total energy bounded over the run;
#   (4) mobile ions leave the high-frequency EM dispersion unchanged (<3%);
#   (5) n_sub=1 vs n_sub=2 agree on the dispersion (convergence, few %);
#   (6) a v≈0.9c relativistic electron beam does not blow up the field energy.

"""
    EMPIC1D(g, Nparticles; n0=1.0, c=5.0, shape=CIC(), relativistic=false,
            mobile=false, mi=1836.0, n_sub=1)

State for a full electromagnetic 1D3V PIC on the periodic grid `g`
(`FourierGrid{1}`). Mobile electrons (charge `−1`, mass `1`) move on a positive
species that is, by default, a uniform immobile ion background of density `n0`
(background charge `+n0`); set `mobile=true` to instead push a second kinetic
ion species (charge `+1`, mass `mi`, default `mi=1836`) supplied to
[`init_empic!`](@ref)/[`step_empic!`](@ref). `ε0 = 1`, so `ω_pe = √n0`. `c` is
the speed of light; `shape` the deposition/gather kernel. Set
`relativistic=true` to use the relativistic Boris push (γ factor) for both
species. `n_sub` (≥1) is the electron subcycling factor: each ion step `dt` is
split into `n_sub` electron/field substeps of `dt/n_sub`.

Fields: `Ex, Ey, Bz` (length-`n` grid arrays), plus internal scratch. The
canonical `Bz` is held at the half time level `n−1/2`.
"""
mutable struct EMPIC1D{T,G,SH<:ShapeFunction}
    g::G
    n0::T
    c::T
    shape::SH
    relativistic::Bool
    mobile::Bool           # true ⇒ ions are a second kinetic species
    mi::T                  # ion mass (electron mass = 1)
    n_sub::Int             # electron subcycling factor (≥1)
    # fields
    Ex::Vector{T}          # longitudinal E, integer level n
    Ey::Vector{T}          # transverse E,  integer level n
    Bz::Vector{T}          # transverse B,  half level n−1/2
    Bzc::Vector{T}         # time-centered Bz^n for the particle push (scratch)
    # currents at n+1/2 (TOTAL over species)
    Jx::Vector{T}
    Jy::Vector{T}
    # density snapshots (TOTAL charge density)
    rho_n::Vector{T}       # ρ^n   (charge density at x^n)
    rho_np1::Vector{T}     # ρ^{n+1}
    ne::Vector{T}          # per-species number-density scratch
    # real-space derivative scratch
    dEy::Vector{T}
    dBz::Vector{T}
    # electron gather buffers
    Exp::Vector{T}
    Eyp::Vector{T}
    Bzp::Vector{T}
    work::Vector{T}
    # ion gather / scratch buffers (sized lazily on first mobile step)
    Exi::Vector{T}
    Eyi::Vector{T}
    Bzi::Vector{T}
    worki::Vector{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    last_dt::Base.RefValue{T}       # preceding full ion drift interval
    last_sub_dt::Base.RefValue{T}   # preceding electron/field drift interval
end

function EMPIC1D(
    g::FourierGrid{1,T},
    Nparticles::Integer;
    n0 = 1.0,
    c = 5.0,
    shape::ShapeFunction = CIC(),
    relativistic::Bool = false,
    mobile::Bool = false,
    mi = 1836.0,
    n_sub::Integer = 1,
) where {T}
    n_sub >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    Np = _particle_length(Nparticles)
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    cT = _require_finite_positive_real("c", c, T)
    miT = _require_finite_positive_real("mi", mi, T)
    n = g.n[1]
    z() = zeros(T, n)
    vz() = zeros(T, Np)
    EMPIC1D{T,typeof(g),typeof(shape)}(
        g,
        n0T,
        cT,
        shape,
        relativistic,
        mobile,
        miT,
        Int(n_sub),
        z(),       # Ex
        z(),       # Ey
        z(),       # Bz
        z(),       # Bzc
        z(),       # Jx
        z(),       # Jy
        z(),       # rho_n
        z(),       # rho_np1
        z(),       # ne
        z(),       # dEy
        z(),       # dBz
        vz(),      # Exp
        vz(),      # Eyp
        vz(),      # Bzp
        vz(),      # work
        T[],       # Exi  (sized lazily)
        T[],       # Eyi
        T[],       # Bzi
        T[],       # worki
        Ref(zero(T)),
        Ref(0),
        Ref(zero(T)),
        Ref(zero(T)),
    )
end

# Ensure the ion gather buffers are sized to the ion particle count.
@inline function _ensure_ion_buffers!(es::EMPIC1D{T}, ni::Integer) where {T}
    if length(es.Exi) != ni
        resize!(es.Exi, ni)
        resize!(es.Eyi, ni)
        resize!(es.Bzi, ni)
        resize!(es.worki, ni)
    end
    return es
end

# Electron gather buffers (also reused as position scratch in _substep_em!) are sized to Np at
# construction; re-size to the current electron count each step so a growing electron population
# (e.g. ionize_mcc! secondaries) composes with the push, symmetric with the ion buffers above.
@inline function _ensure_electron_buffers!(es::EMPIC1D{T}, ne::Integer) where {T}
    if length(es.Exp) != ne
        resize!(es.Exp, ne)
        resize!(es.Eyp, ne)
        resize!(es.Bzp, ne)
        resize!(es.work, ne)
    end
    return es
end

@inline function _require_empic_electrons(e::ParticleSet{D,T}) where {D,T}
    q = _require_finite_real("electron charge q", e.q, T)
    m = _require_finite_real("electron mass m", e.m, T)
    q == -one(T) ||
        throw(ArgumentError("electromagnetic PIC requires electron ParticleSet with q = -1"))
    m == one(T) ||
        throw(ArgumentError("electromagnetic PIC requires electron ParticleSet with m = 1"))
    return nothing
end

@inline function _require_empic_ions(ions::ParticleSet{D,T}, mi::T) where {D,T}
    q = _require_finite_real("ion charge q", ions.q, T)
    m = _require_finite_real("ion mass m", ions.m, T)
    q == one(T) || throw(ArgumentError("electromagnetic PIC requires ion ParticleSet with q = 1"))
    m > zero(T) || throw(ArgumentError("electromagnetic PIC requires ion ParticleSet with m > 0"))
    # the dynamics use the ParticleSet's per-particle mass; require it to match the constructor
    # `mi` so a user following the docstring cannot silently get a different ion mass.
    m == mi || throw(
        ArgumentError(
            "ion ParticleSet mass m=$m must equal the EMPIC constructor mi=$mi (the ion mass)",
        ),
    )
    return nothing
end

# ---------------------------------------------------------------- charge density

# Range status must survive per-species normalization: even a correctly rounded
# density can lose a contribution that the total charge can represent.
function _empic_raw_charge_density!(
    rho::Array{T,D},
    ne::Array{T,D},
    es,
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}},
) where {D,T}
    mobile = es.mobile && ions !== nothing
    eligible = true
    lost = false
    if mobile
        fill!(rho, zero(T))
        for species in (ions, e)
            state = _density_with_range!(ne, species, es.g, es.shape)
            eligible &= state.eligible && isfinite(species.q)
            lost |= state.lost
            @inbounds for I in eachindex(rho)
                contribution = species.q * ne[I]
                value = rho[I] + contribution
                lost |=
                    _deposition_product_lost(species.q, ne[I], contribution) ||
                    _deposition_range_lost(value)
                rho[I] = value
            end
        end
    else
        state = _density_with_range!(ne, e, es.g, es.shape)
        eligible = state.eligible && isfinite(es.n0)
        lost = state.lost
        @inbounds for I in eachindex(rho)
            value = es.n0 - ne[I]
            lost |= _deposition_range_lost(value)
            rho[I] = value
        end
    end
    if eligible && lost
        # Preserve the ordinary electron-density scratch. The complete charge
        # replay reads original particle operands and combines before conversion.
        exact = _exact_moment(e, (e.weight,), es.g, es.shape, mobile ? e.q : one(T))
        if mobile
            ion_exact = _exact_moment(ions, (ions.weight,), es.g, es.shape, ions.q)
            for I in eachindex(exact)
                exact[I] += ion_exact[I]
            end
        else
            background = _exact_rational(es.n0)
            for I in eachindex(exact)
                exact[I] = background - exact[I]
            end
        end
        for I in eachindex(rho)
            rho[I] = _from_exact(T, exact[I])
        end
    end
    return (; eligible, lost)
end

# Exact totals use bounded integer accumulation when the represented exponents
# permit it; the wider path retains every finite stored operand.
function _empic_exact_sum(values::AbstractArray{T}) where {T<:Union{Float16,Float32,Float64}}
    Q = Rational{BigInt}
    isempty(values) && return Q(0)
    reference = first(values)
    isfinite(reference) || throw(ArgumentError("exact dyadic sum requires finite values"))
    minimum_exponent = iszero(reference) ? typemax(Int) : last(frexp(reference))
    maximum_exponent = iszero(reference) ? typemin(Int) : minimum_exponent
    all_equal = true
    for value in values
        if value != reference
            all_equal = false
            isfinite(value) || throw(ArgumentError("exact dyadic sum requires finite values"))
            if !iszero(value)
                _, e = frexp(value)
                minimum_exponent = min(minimum_exponent, e)
                maximum_exponent = max(maximum_exponent, e)
            end
        end
    end
    if all_equal
        return length(values) == 1 ? Q(reference) : Q(reference) * length(values)
    end
    minimum_exponent == typemax(Int) && return Q(0)
    p = precision(T)
    # n < 2^ndigits(n), so this bounds every absolute partial sum strictly
    # below 2^127 in the selected quantum, including arbitrary sign order.
    if p + maximum_exponent - minimum_exponent + ndigits(length(values); base = 2) > 127
        return _empic_exact_repeated_sum(values)
    end
    total = zero(Int128)
    for value in values
        if !iszero(value)
            fraction, e = frexp(value)
            significand = Int128(ldexp(fraction, p))
            total += significand << (e - minimum_exponent)
        end
    end
    quantum = minimum_exponent - p
    integer = BigInt(total)
    return quantum >= 0 ? Q(integer << quantum) : integer // (BigInt(1) << -quantum)
end

function _empic_exact_repeated_sum(values)
    Q = Rational{BigInt}
    total = Q(0)
    previous = zero(eltype(values))
    count = 0
    for value in values
        if count > 0 && value != previous
            total += count == 1 ? Q(previous) : Q(previous) * count
            count = 0
        end
        previous = value
        count += 1
    end
    count > 0 && (total += count == 1 ? Q(previous) : Q(previous) * count)
    return total
end

function _empic_charge_shift_bounds(
    ::Type{T},
    shift::Rational{BigInt},
) where {T<:Union{Float32,Float64}}
    denominator(shift) > 0 || throw(ArgumentError("uniform shift must be finite"))
    rounded = _from_exact(T, shift)
    if isinf(rounded)
        return rounded > 0 ? (floatmax(T), rounded) : (rounded, -floatmax(T))
    end
    represented = _exact_rational(rounded)
    represented == shift && return (rounded, rounded)
    return represented < shift ? (rounded, nextfloat(rounded)) : (prevfloat(rounded), rounded)
end

# Monotonicity certifies a cell when both bracketing shifts round to the same
# result. Otherwise round the complete difference, including overflowing shifts.
function _empic_complete_charge_shift!(
    cells::AbstractArray{T},
    shift::Rational{BigInt},
) where {T<:Union{Float32,Float64}}
    lower, upper = _empic_charge_shift_bounds(T, shift)
    all(isfinite, cells) || throw(ArgumentError("uniform shift requires finite input cells"))
    for I in eachindex(cells)
        value = cells[I]
        lower_result, upper_result = value - upper, value - lower
        if isequal(lower_result, upper_result)
            cells[I] = lower_result
        else
            cells[I] = _from_exact(T, _exact_rational(value) - shift)
        end
    end
    return cells
end

# The target uses actual stored particle weights and the complete width product.
# After the uniform correction, the total error is exactly the sum of the final
# cell-rounding errors; a non-dyadic target need not be attained exactly.
function _empic_reconcile_charge!(rho, es, e, ions)
    volume = prod(_exact_rational.(es.g.dx))
    if es.mobile && ions !== nothing
        particles =
            _exact_rational(e.q) * _empic_exact_sum(e.weight) +
            _exact_rational(ions.q) * _empic_exact_sum(ions.weight)
        target = particles / volume
    else
        target = length(rho) * _exact_rational(es.n0) - _empic_exact_sum(e.weight) / volume
    end
    shift = (_empic_exact_sum(rho) - target) / length(rho)
    _empic_complete_charge_shift!(rho, shift)
    return rho
end

# Capture reconciles built-in partition shapes after the raw arithmetic succeeds.
# Custom stencils retain their represented coefficients and nonfinite data retain
# the raw propagation. Neither the density scratch nor particle storage is changed.
function _empic_charge_density!(
    rho::Array{T,D},
    ne::Array{T,D},
    es,
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}},
) where {D,T}
    state = _empic_raw_charge_density!(rho, ne, es, e, ions)
    if state.eligible && es.shape isa Union{NGP,CIC,TSC} && all(isfinite, rho)
        _empic_reconcile_charge!(rho, es, e, ions)
    end
    return rho
end

# Total charge density at current positions: ions plus electrons when mobile,
# or the uniform positive background minus electron number density otherwise.
function _charge_density!(
    rho::Vector{T},
    ne::Vector{T},
    es::EMPIC1D{T},
    e::ParticleSet{1,T},
    ions::Union{Nothing,ParticleSet{1,T}},
) where {T}
    return _empic_charge_density!(rho, ne, es, e, ions)
end

# ---------------------------------------------------------------- Esirkepov-exact
# 1-D charge-conserving longitudinal current. Continuity ∂tρ + ∂xJx = 0 in
# spectral form is  (ρ̂^{n+1}−ρ̂^n)/dt + ik_m Ĵx_m = 0, so
#   Ĵx_m = −(ρ̂^{n+1}−ρ̂^n)/(dt·ik_m)   for k_m ≠ 0.
# The k=0 mode of continuity is automatically satisfied (total charge conserved
# ⇒ ρ̂_0 constant), and is not constrained by it; the physical net current
# (k=0) is set from the directly deposited current mean (summed over species).
# Because the SAME ik as every other spectral operator is used (Nyquist zeroed),
# the discrete identity ∂tρ + ∂xJx = 0 holds to roundoff for every represented
# mode of the TOTAL charge density.
function _esirkepov_Jx!(
    es::EMPIC1D{T},
    e::ParticleSet{1,T},
    ions::Union{Nothing,ParticleSet{1,T}},
    dt::T,
) where {T}
    g = es.g
    n = g.n[1]
    cb = g.cbuf
    tb = g.tbuf
    # direct-deposit current mean (k=0 physical net current): Jx0 = Σ_species q Σ w v_x / L
    Jx0 = _combined_Jx0(es, e, ions)
    # ρ̂^{n+1} − ρ̂^n
    @inbounds for i = 1:n
        cb[i] = Complex{T}(es.rho_np1[i] - es.rho_n[i])
    end
    g.plan * cb
    ik = g.ik[1]
    @inbounds for m = 1:n
        ikm = ik[m]
        if ikm == 0
            tb[m] = zero(Complex{T})         # k=0 handled separately (set after ifft)
        else
            tb[m] = -(cb[m] / dt) / ikm      # Ĵx_m
        end
    end
    g.iplan * tb
    @inbounds for i = 1:n
        es.Jx[i] = real(tb[i]) + Jx0          # add back the k=0 net current
    end
    return es.Jx
end

@inline function _species_Jx0_with_range(s::ParticleSet{1,T}, L::T) where {T}
    vx = s.v[1]
    w = s.weight
    acc = zero(T)
    finite = isfinite(s.q) && isfinite(L) && L > zero(T)
    lost = false
    @inbounds for p in eachindex(w)
        product = w[p] * vx[p]
        finite &= isfinite(w[p]) && isfinite(vx[p])
        lost |= _deposition_product_lost(w[p], vx[p], product)
        acc += product
        lost |= _deposition_range_lost(acc)
    end
    numerator = s.q * acc
    value = numerator / L
    lost |=
        _deposition_product_lost(s.q, acc, numerator) ||
        _deposition_range_lost(value) ||
        (iszero(value) && !iszero(numerator))
    if finite && lost
        value = _from_exact(T, _exact_species_Jx0(s, L))
    end
    return (; value, eligible = finite, lost)
end

function _exact_species_Jx0(s, L)
    exact = zero(Rational{BigInt})
    for p in eachindex(s.weight)
        exact += _exact_rational(s.weight[p]) * _exact_rational(s.v[1][p])
    end
    return exact * _exact_rational(s.q) / _exact_rational(L)
end

@inline _species_Jx0(s::ParticleSet{1,T}, L::T) where {T} = _species_Jx0_with_range(s, L).value

function _combined_Jx0(es::EMPIC1D{T}, e, ions) where {T}
    electron = _species_Jx0_with_range(e, es.g.L[1])
    if es.mobile && ions !== nothing
        ion = _species_Jx0_with_range(ions, es.g.L[1])
        value = electron.value + ion.value
        if electron.eligible &&
           ion.eligible &&
           (electron.lost || ion.lost || _deposition_range_lost(value))
            value = _from_exact(
                T,
                _exact_species_Jx0(e, es.g.L[1]) + _exact_species_Jx0(ions, es.g.L[1]),
            )
        end
        return value
    end
    return electron.value
end

# The ordinary species passes restore their positions before this cold replay.
# Borrow each species' own midpoint buffer only while its exact stencil is read.
function _exact_midpoint_current(s, mid, g, shape, component)
    _swap_positions_with_midpoints!(s, mid)
    try
        return _exact_moment(s, (s.weight, s.v[component]), g, shape, s.q)
    finally
        _swap_positions_with_midpoints!(s, mid)
    end
end

function _combine_midpoint_current!(
    out::Array{T},
    es,
    e,
    ions,
    emid,
    imid,
    component,
    electron,
    ion,
) where {T}
    if electron.eligible &&
       ion.eligible &&
       (electron.lost || ion.lost || any(_deposition_range_lost, out))
        exact_e = _exact_midpoint_current(e, emid, es.g, es.shape, component)
        exact_i = _exact_midpoint_current(ions, imid, es.g, es.shape, component)
        for I in eachindex(out)
            out[I] = _from_exact(T, exact_e[I] + exact_i[I])
        end
    end
    return out
end

# Transverse current Jy^{n+1/2} = Σ_species q Σ w v_y S(x^{n+1/2}) (CIC). For each
# species the deposit must happen at its OWN midpoint positions; the caller parks
# the midpoint positions in the species' x array before calling and restores
# x^{n+1} afterwards. `accumulate` selects whether to zero Jy first.
function _add_transverse_Jy!(
    es::EMPIC1D{T},
    s::ParticleSet{1,T},
    work::Vector{T},
    scratch::Vector{T},
    accumulate::Bool,
) where {T}
    @inbounds @. work = s.weight * s.v[2]
    _validate_deposition_inputs(scratch, s, work, es.g)
    factors = (s.weight, s.v[2])
    state = _deposit_scalar_kernel!(scratch, s, factors, es.g, es.shape, Val(true))
    return _finish_moment!(
        es.Jy,
        scratch,
        s,
        factors,
        es.g,
        es.shape,
        s.q,
        state,
        Val(:prescale),
        accumulate,
    )
end

# Remove the Nyquist DFT mode of a real grid vector in place (no-op for odd n):
# f_i -= (-1)^(i-1) · (1/n) Σ_j (-1)^(j-1) f_j. The spectral ∂x zeroes the
# Nyquist multiplier, so ∂xBz has no Nyquist content and Ey's Nyquist mode would
# otherwise integrate the raw deposition noise of Jy with no restoring term —
# the same undamped random walk _esirkepov_Jx! excludes from the longitudinal
# channel, and the transverse twin of EMPIC's continuity correction, which
# zeroes every current component on the pure-Nyquist modes.
function _zero_nyquist_mode!(f::Vector{T}) where {T}
    n = length(f)
    iseven(n) || return f
    acc = zero(T)
    s = one(T)
    eligible = T <: Union{Float32,Float64}
    lost = false
    largest = zero(T)
    @inbounds for i = 1:n
        value = f[i]
        if T <: Union{Float32,Float64}
            eligible &= isfinite(value)
            largest = max(largest, abs(value))
        end
        acc += s * value
        lost |= _deposition_range_lost(acc)
        s = -s
    end
    total = acc
    acc /= n
    lost |= _deposition_range_lost(acc) || (iszero(acc) && !iszero(total))
    if eligible && n > 0 && (lost || largest > floatmax(T) - abs(acc))
        # Keep the original alternating sum and its divided mean complete until
        # each final cell is rounded. The two views are disjoint.
        odd = view(f, 1:2:n)
        even = view(f, 2:2:n)
        mean = (_empic_exact_sum(odd) - _empic_exact_sum(even)) / n
        _empic_complete_charge_shift!(odd, mean)
        _empic_complete_charge_shift!(even, -mean)
        return f
    end
    s = one(T)
    @inbounds for i = 1:n
        f[i] -= s * acc
        s = -s
    end
    return f
end

# ---------------------------------------------------------------- relativistic push

# Relativistic Boris push for a single particle (Birdsall & Langdon / Vay-free
# classic Boris): works in momentum u = γv. Returns updated (vx,vy,vz). c is the
# speed of light. Reduces to the non-relativistic boris_kick as c→∞.
@inline function _boris_rel(
    vx::T,
    vy::T,
    vz::T,
    Ex::T,
    Ey::T,
    Ez::T,
    Bx::T,
    By::T,
    Bz::T,
    qm::T,
    dt::T,
    c::T,
) where {T}
    c2 = c * c
    βx = vx / c
    βy = vy / c
    βz = vz / c
    β2 = βx * βx + βy * βy + βz * βz
    (isfinite(β2) && β2 < one(T)) ||
        throw(ArgumentError("_boris_rel: input velocity must be finite and satisfy |v| < c"))
    ginit = one(T) / sqrt(one(T) - β2)
    # to momentum u = γ v
    ux = ginit * vx
    uy = ginit * vy
    uz = ginit * vz
    h = qm * dt / 2
    umx = ux + h * Ex
    umy = uy + h * Ey
    umz = uz + h * Ez
    γm = sqrt(one(T) + (umx * umx + umy * umy + umz * umz) / c2)
    tx = h * Bx / γm
    ty = h * By / γm
    tz = h * Bz / γm
    t2 = tx * tx + ty * ty + tz * tz
    f = 2 / (1 + t2)
    sx = f * tx
    sy = f * ty
    sz = f * tz
    upx = umx + (umy * tz - umz * ty)
    upy = umy + (umz * tx - umx * tz)
    upz = umz + (umx * ty - umy * tx)
    unx = umx + (upy * sz - upz * sy)
    uny = umy + (upz * sx - upx * sz)
    unz = umz + (upx * sy - upy * sx)
    ux2 = unx + h * Ex
    uy2 = uny + h * Ey
    uz2 = unz + h * Ez
    γf = sqrt(one(T) + (ux2 * ux2 + uy2 * uy2 + uz2 * uz2) / c2)
    return (ux2 / γf, uy2 / γf, uz2 / γf)
end

# Push one species' velocities (in place) given per-particle gathered fields.
# Ez,Bx,By ≡ 0 for the 1D polarization split.
@inline function _push_velocities!(
    es::EMPIC1D{T},
    s::ParticleSet{1,T},
    Exp::Vector{T},
    Eyp::Vector{T},
    Bzp::Vector{T},
    dt::T,
) where {T}
    qm = s.q / s.m
    vx = s.v[1]
    vy = s.v[2]
    vz = s.v[3]
    z = zero(T)
    if es.relativistic
        cc = es.c
        @inbounds for p in eachindex(s.weight)
            nx, ny, nz =
                _boris_rel(vx[p], vy[p], vz[p], Exp[p], Eyp[p], z, z, z, Bzp[p], qm, dt, cc)
            vx[p] = nx
            vy[p] = ny
            vz[p] = nz
        end
    else
        @inbounds for p in eachindex(s.weight)
            nx, ny, nz = boris_kick(vx[p], vy[p], vz[p], Exp[p], Eyp[p], z, z, z, Bzp[p], qm, dt)
            vx[p] = nx
            vy[p] = ny
            vz[p] = nz
        end
    end
    return s
end

# Drift one species x^n → x^{n+1}; record the wrapped midpoint x^{n+1/2} into
# `mid`. Caller wraps x afterwards with apply_periodic!.
@inline function _drift_record_mid!(s::ParticleSet{1,T}, mid::Vector{T}, dt::T, L::T) where {T}
    xx = s.x[1]
    vx = s.v[1]
    @inbounds for p in eachindex(xx)
        xmid = xx[p] + (dt / 2) * vx[p]
        xx[p] += dt * vx[p]
        mid[p] = mod(xmid, L)
    end
    return s
end

# ---------------------------------------------------------------- init

"""
    init_empic!(es, electrons[, ions])

Initialize the run: deposit ρ^n from the loaded species and solve the initial
longitudinal field Ex from Poisson (−∂²φ = ρ/ε0, E = −∂φ; spectrally
Êx_m = −iρ̂_m/(ε0 k_m), Êx_0 = 0). The transverse fields Ey, Bz keep whatever
seed the caller has placed in `es.Ey`, `es.Bz`. Pass `ions` (a second
`ParticleSet{1}` with `q=+1`) when `es.mobile`; omit it for the immobile-ion
background. Call once after loading.
"""
function init_empic!(
    es::EMPIC1D{T},
    e::ParticleSet{1,T},
    ions::Union{Nothing,ParticleSet{1,T}} = nothing,
) where {T}
    _require_empic_electrons(e)
    _ensure_electron_buffers!(es, nparticles(e))
    if es.mobile
        ions === nothing &&
            throw(ArgumentError("mobile=true requires an ion ParticleSet in init_empic!"))
        _require_empic_ions(ions, es.mi)
        _ensure_ion_buffers!(es, nparticles(ions))
    end
    g = es.g
    n = g.n[1]
    cb = g.cbuf
    _charge_density!(es.rho_n, es.ne, es, e, ions)
    @inbounds for i = 1:n
        cb[i] = Complex{T}(es.rho_n[i])
    end
    g.plan * cb
    k = g.kvec[1]
    @inbounds for m = 1:n
        km = k[m]
        cb[m] = km == 0 ? zero(Complex{T}) : (-im * cb[m]) / km   # Êx = −iρ̂/k
    end
    g.iplan * cb
    @inbounds for i = 1:n
        es.Ex[i] = real(cb[i])
    end
    es.time[] = zero(T)                 # (re)start at t=0; step==0 triggers leapfrog priming
    es.step[] = 0
    es.last_dt[] = zero(T)
    es.last_sub_dt[] = zero(T)
    return es
end

# ---------------------------------------------------------------- substep

# Keep the ordinary Maxwell operation order, but detect intermediate range loss
# before accepting its result. A subnormal product can lose bits which a later
# large factor restores to a normal, representable electric-field update.
@inline function _em_increment_lost_range(forcing::T, increment::T, candidate::T) where {T}
    return !isfinite(forcing) ||
           !isfinite(increment) ||
           (!iszero(forcing) && abs(increment) < floatmin(T)) ||
           !isfinite(candidate)
end

@inline function _em_update_lost_range(
    c2::T,
    derivative::T,
    magnetic::T,
    forcing::T,
    increment::T,
    candidate::T,
) where {T}
    smallest_normal = floatmin(T)
    return !isfinite(c2) ||
           c2 < smallest_normal ||
           !isfinite(magnetic) ||
           (!iszero(derivative) && abs(magnetic) < smallest_normal) ||
           _em_increment_lost_range(forcing, increment, candidate)
end

function _em_exact_update_coefficients(c, dt)
    exact_dt = _exact_rational(dt)
    exact_c = _exact_rational(c)
    return exact_dt, exact_dt * exact_c * exact_c
end

function _em_exact_electric_update(
    E::T,
    curlB::Union{T,Rational{BigInt}},
    J::T,
    coefficients,
) where {T}
    exact_dt, exact_dtc2 = coefficients
    return _from_exact(
        T,
        _exact_rational(E) + exact_dtc2 * _exact_rational(curlB) - exact_dt * _exact_rational(J),
    )
end

# The public spectral operators retain their input validation and ordinary FFT
# arithmetic. Only a range-sensitive result is recomputed as separately scaled
# scalar derivatives. These are still FFTW-precision derivatives, not exact DFTs.
function _em_field_exponents(field::AbstractArray{T}) where {T}
    smallest, largest = T(Inf), zero(T)
    @inbounds @simd for i in eachindex(field)
        magnitude = abs(field[i])
        largest = max(largest, magnitude)
        smallest = min(smallest, ifelse(iszero(magnitude), T(Inf), magnitude))
    end
    isfinite(largest) || throw(ArgumentError("spectral operator input values must be finite"))
    iszero(largest) && return typemax(Int), typemin(Int)
    return frexp(smallest)[2], frexp(largest)[2]
end

@inline function _em_scale_bounds(::Type{T}, lo, hi, klo, khi, n) where {T}
    (lo == typemax(Int) || klo == typemax(Int)) && return 0, 0, false
    emin = frexp(floatmin(T))[2]
    emax = frexp(floatmax(T))[2]
    ceil_log_n = 8sizeof(Int) - leading_zeros(n - 1)
    # Match the canonical operator's 8N multiplier headroom with two additional
    # binary places, and keep both the input and its smallest k product normal.
    lower = hi + max(0, khi) - (emax - ceil_log_n - 5)
    upper = lo + min(0, klo) - (emin + 2)
    return lower, upper, !(lower <= 0 <= upper)
end

function _em_power_of_two(shift::Int)
    return shift >= 0 ? (BigInt(1) << shift) // BigInt(1) : BigInt(1) // (BigInt(1) << -shift)
end

function _em_dyadic_bit_range(value::Rational{BigInt})
    iszero(value) && return nothing
    n, d = abs(numerator(value)), denominator(value)
    ispow2(d) || throw(ArgumentError("normalized EM transform requires dyadic inputs"))
    shift = ndigits(d; base = 2) - 1
    return trailing_zeros(n) - shift, ndigits(n; base = 2) - 1 - shift
end

# Exact signed digits preserve low bits such as those in floatmax(T)-tiny.
# FFTW still rounds each transform; this is range protection, not an exact DFT.
# Stream the occupied bands instead of retaining one array for every exponent.
function _em_dyadic_transform(
    values::Array{Complex{Rational{BigInt}},D},
    g::FourierGrid{D,T},
    inverse::Bool = false,
) where {D,T}
    ceil_log_n = 8sizeof(Int) - leading_zeros(length(values) - 1)
    # Reserve integer-sum headroom as well as overflow headroom. When the grid
    # exceeds the mantissa's integer capacity even a one-bit band can round.
    width = max(1, precision(T) - ceil_log_n - 2)
    lowest, highest = typemax(Int), typemin(Int)
    for value in values, component in (real(value), imag(value))
        bits = _em_dyadic_bit_range(component)
        bits === nothing && continue
        lowest = min(lowest, bits[1])
        highest = max(highest, bits[2])
    end
    result = fill(zero(Complex{Rational{BigInt}}), size(values))
    lowest == typemax(Int) && return result
    work = similar(g.cbuf)
    mask = (BigInt(1) << width) - 1
    # The finite input/product bit ranges bound this loop and the integer
    # widths. Empty bands consume no FFT and no additional retained array.
    for band = fld(lowest, width):fld(highest, width)
        scale = _em_power_of_two(band * width + width)
        occupied = false
        @inbounds for I in eachindex(values)
            digits = map((real(values[I]), imag(values[I]))) do component
                # The input scan certified a power-of-two denominator. Shift
                # the magnitude before masking, preserving floor division for
                # negative bands and keeping the original signed-digit order.
                offset = trailing_zeros(denominator(component)) + band * width
                digit = (abs(numerator(component)) >> offset) & mask
                normalized = ldexp(T(digit), -width)
                signbit(component) ? -normalized : normalized
            end
            work[I] = complex(digits...)
            occupied |= !iszero(work[I])
        end
        occupied || continue
        inverse ? (g.iplan * work) : (g.plan * work)
        @inbounds for I in eachindex(work)
            coefficient = _spectral_exact_complex(work[I])::Complex{Rational{BigInt}}
            result[I] += scale * coefficient
        end
    end
    return result
end

@inline function _em_density_transform_bounds(::Type{T}, n::Int) where {T}
    ceil_log_n = 8sizeof(Int) - leading_zeros(n - 1)
    # Leave headroom for forward coefficients and the unnormalized inverse
    # butterfly sums, before the inverse plan applies its factor 1/N.
    headroom = 2ceil_log_n + 5
    return ldexp(floatmin(T), headroom), ldexp(floatmax(T), -headroom)
end

@inline function _em_density_transform_lost_range(value, bounds)
    magnitude = abs(value)
    return !isfinite(value) || (!iszero(value) && !(bounds[1] <= magnitude <= bounds[2]))
end

function _em_continuity_scaled(g::FourierGrid{D,T}, rho_n, rho_np1, J, dt::T) where {D,T}
    delta = map(rho_np1, rho_n) do following, preceding
        complex(_exact_rational(following) - _exact_rational(preceding), zero(Rational{BigInt}))
    end
    density_hat = _em_dyadic_transform(delta, g)
    current_hat = [
        _em_dyadic_transform(map(value -> complex(_exact_rational(value)), field), g) for
        field in Iterators.take(J, D)
    ]
    residual_hat = similar(density_hat)
    exact_dt = _exact_rational(dt)
    dc = first(CartesianIndices(density_hat))
    @inbounds for I in CartesianIndices(density_hat)
        k = _spectral_wavevector(g, I)
        I != dc && _spectral_mode_is_null(k) && (density_hat[I] = zero(eltype(density_hat)))
        divergence = zero(eltype(density_hat))
        for d = 1:D
            divergence += im * _exact_rational(k[d]) * current_hat[d][I]
        end
        residual_hat[I] = density_hat[I] + exact_dt * divergence
    end
    density = _em_dyadic_transform(density_hat, g, true)
    residual = _em_dyadic_transform(residual_hat, g, true)
    rmax, scale = zero(Rational{BigInt}), zero(Rational{BigInt})
    @inbounds for I in eachindex(density)
        rmax = max(rmax, abs(real(residual[I])))
        scale = max(scale, abs(real(density[I])))
    end
    return _from_exact(T, rmax / (iszero(scale) ? one(scale) : scale))
end

@inline function _em_continuity_finish(
    g::FourierGrid{D,T},
    rho_n,
    rho_np1,
    J,
    dt::T,
    density,
    divergence,
    lost::Bool,
) where {D,T}
    rmax, scale = zero(T), zero(T)
    @inbounds for I in eachindex(density)
        delta = real(density[I])
        increment = dt * divergence[I]
        residual = delta + increment
        lost |= _spectral_product_lost_range(dt, divergence[I], increment) || !isfinite(residual)
        rmax = max(rmax, abs(residual))
        scale = max(scale, abs(delta))
    end
    if !lost
        # Canonical divergence uses a common FFT scale across components.
        # Include every field/axis pair when checking that shared scaling.
        wave_exponents = map(_em_field_exponents, g.kvec)
        for field in Iterators.take(J, D)
            field_exponents = _em_field_exponents(field)
            for wave_exponent in wave_exponents
                lost |=
                    _em_scale_bounds(T, field_exponents..., wave_exponent..., length(density))[3]
            end
        end
    end
    if lost && isfinite(dt) && all(isfinite, rho_n) && all(isfinite, rho_np1)
        return _em_continuity_scaled(g, rho_n, rho_np1, J, dt)
    end
    return rmax / (scale > 0 ? scale : one(T))
end

function _em_normalized_derivative_parts(
    field::Array{T,D},
    g,
    axis,
    exponents,
    kexponents,
) where {T,D}
    parts = Tuple{Array{T,D},Rational{BigInt}}[]
    lo, hi = exponents
    klo, khi = kexponents
    (lo == typemax(Int) || klo == typemax(Int)) && return parts
    lower, upper, needed = _em_scale_bounds(T, lo, hi, klo, khi, length(field))
    if !needed
        part = similar(field)
        deriv!(part, field, g, axis)
        push!(parts, (part, BigInt(1) // BigInt(1)))
        return parts
    end

    base_lower, base_upper, _ = _em_scale_bounds(T, 0, 0, klo, khi, length(field))
    width = base_upper - base_lower
    width >= 0 || throw(ArgumentError("spectral multiplier range leaves no normal input band"))
    work = similar(field)
    band_hi = hi
    # Each nonzero datum belongs to exactly one disjoint band. At most
    # cld(hi - lo + 1, width + 1) arrays survive this loop; hi - lo is bounded
    # by the hardware exponent range. Empty bands allocate no output array.
    while band_hi >= lo
        band_lo = max(lo, band_hi - width)
        lower, upper, _ = _em_scale_bounds(T, band_lo, band_hi, klo, khi, length(field))
        shift = clamp(band_hi + fld(khi, 2), lower, upper)
        nonzero = false
        @inbounds for i in eachindex(field)
            value = field[i]
            exponent = iszero(value) ? typemin(Int) : frexp(abs(value))[2]
            selected = band_lo <= exponent <= band_hi
            work[i] = selected ? ldexp(value, -shift) : zero(T)
            nonzero |= selected
        end
        if nonzero
            part = similar(field)
            deriv!(part, work, g, axis)
            push!(parts, (part, _em_power_of_two(shift)))
        end
        band_hi = band_lo - 1
    end
    return parts
end

function _em_derivative_parts!(out, field, g, axis)
    deriv!(out, field, g, axis)
    exponents = _em_field_exponents(field)
    kexponents = _em_field_exponents(g.kvec[axis])
    _em_scale_bounds(eltype(field), exponents..., kexponents..., length(field))[3] || return nothing
    return _em_normalized_derivative_parts(field, g, axis, exponents, kexponents)
end

function _em_curl_parts!(out, fields, g::FourierGrid{D,T}) where {D,T}
    curl!(out, fields, g)
    # Every component contributes to the canonical curl's common FFT scale,
    # including longitudinal components and diagonal field/axis pairs.
    exponents = map(_em_field_exponents, fields)
    kexponents = map(_em_field_exponents, g.kvec)
    needed = false
    for c = 1:3, axis = 1:D
        needed |= _em_scale_bounds(T, exponents[c]..., kexponents[axis]..., length(fields[c]))[3]
    end
    needed || return nothing
    return ntuple(3) do c
        plus, minus = mod1(c + 1, 3), mod1(c + 2, 3)
        parts = Tuple{Array{T,D},Rational{BigInt}}[]
        if plus <= D
            append!(
                parts,
                _em_normalized_derivative_parts(
                    fields[minus],
                    g,
                    plus,
                    exponents[minus],
                    kexponents[plus],
                ),
            )
        end
        if minus <= D
            negative = _em_normalized_derivative_parts(
                fields[plus],
                g,
                minus,
                exponents[plus],
                kexponents[minus],
            )
            for (part, scale) in negative
                push!(parts, (part, -scale))
            end
        end
        parts
    end
end

function _em_exact_derivative(parts, i)
    result = BigInt(0) // BigInt(1)
    for (part, scale) in parts
        result += scale * _exact_rational(part[i])
    end
    return result
end

function _em_magnetic_update!(
    out,
    old,
    derivative,
    parts,
    preceding::T,
    following::T,
    backward::Bool = false,
) where {T}
    before, after = preceding / 2, following / 2
    interval = before + after
    half_lost_range =
        (!iszero(preceding) && abs(before) < floatmin(T)) ||
        (!iszero(following) && abs(after) < floatmin(T))
    exact_interval = nothing
    @inbounds for i in eachindex(out)
        value, slope = old[i], derivative[i]
        increment = interval * slope
        candidate = backward ? value + increment : value - increment
        if (
               parts !== nothing ||
               half_lost_range ||
               _em_increment_lost_range(slope, increment, candidate)
           ) &&
           isfinite(value) &&
           (parts !== nothing || isfinite(slope))
            if exact_interval === nothing
                exact_interval = (_exact_rational(preceding) + _exact_rational(following)) / 2
            end
            exact_slope =
                parts === nothing ? _exact_rational(slope) : _em_exact_derivative(parts, i)
            exact_increment = exact_interval * exact_slope
            candidate = _from_exact(
                T,
                backward ? _exact_rational(value) + exact_increment :
                _exact_rational(value) - exact_increment,
            )
        end
        out[i] = candidate
    end
    return out
end

# One field+electron substep of size dt. Pushes the electrons `e` always; pushes
# the ions `ions` only when `push_ions` is true (heavy ions are pushed once per
# full ion step, on its first substep). When mobile, the ions
# always contribute their charge (rho snapshots) and current to the field solve.
function _substep_em!(
    es::EMPIC1D{T},
    e::ParticleSet{1,T},
    ions::Union{Nothing,ParticleSet{1,T}},
    dt::T,
    push_ions::Bool,
    dt_previous::T,
    dt_ion_push::T,
) where {T}
    g = es.g
    n = g.n[1]
    L = g.L[1]
    c2 = es.c * es.c
    z = zero(T)
    kick_dt = dt_previous / 2 + dt / 2

    # snapshot ρ^n (from x^n) for the charge-conserving current
    _charge_density!(es.rho_n, es.ne, es, e, ions)

    # (1) center B using the preceding substep's half interval, which is where
    # the stored field lies relative to the current integer time.
    electric_parts = _em_derivative_parts!(es.dEy, es.Ey, g, 1)
    _em_magnetic_update!(es.Bzc, es.Bz, es.dEy, electric_parts, dt_previous, z)

    # (2) gather fields at time n and push electrons
    gather_scalar!(es.Exp, es.Ex, e, g, es.shape)
    gather_scalar!(es.Eyp, es.Ey, e, g, es.shape)
    gather_scalar!(es.Bzp, es.Bzc, e, g, es.shape)
    _push_velocities!(es, e, es.Exp, es.Eyp, es.Bzp, kick_dt)
    _drift_record_mid!(e, es.work, dt, L)
    apply_periodic!(e, (z,), (L,))

    # The ion kick spans its preceding and next half levels, once on the first
    # substep; each electron substep below drifts with that new ion velocity.
    if es.mobile && ions !== nothing && push_ions
        gather_scalar!(es.Exi, es.Ex, ions, g, es.shape)
        gather_scalar!(es.Eyi, es.Ey, ions, g, es.shape)
        gather_scalar!(es.Bzi, es.Bzc, ions, g, es.shape)
        _push_velocities!(es, ions, es.Exi, es.Eyi, es.Bzi, dt_ion_push)
    end
    # ion POSITION drift: EVERY substep by dt (with the current ion velocity), so the ion charge and
    # current deposited on each substep sit at the time-consistent position. Drifting the whole
    # dt_ion at once on the straddling substep lumped the ion current into one substep → 1st-order;
    # splitting it makes the n_sub≥2 mobile scheme 2nd-order (n_sub=1 is unchanged: dt = dt_ion).
    if es.mobile && ions !== nothing
        _drift_record_mid!(ions, es.worki, dt, L)
        apply_periodic!(ions, (z,), (L,))
    end

    # (3) currents at n+1/2
    #   transverse Jy from midpoint positions: swap x^{n+1} ↔ x^{n+1/2} per species
    _deposit_Jy_at_midpoints!(es, e, ions)
    #   longitudinal Jx from total ρ^n, ρ^{n+1} (charge-conserving)
    _charge_density!(es.rho_np1, es.ne, es, e, ions)
    _esirkepov_Jx!(es, e, ions, dt)

    # (4) Advance B between its preceding and next half levels.
    _em_magnetic_update!(es.Bz, es.Bz, es.dEy, electric_parts, dt_previous, dt)

    # (5) E advance:
    #   Ey^{n+1} = Ey^n − dt(c²∂xBz^{n+1/2} + Jy)
    magnetic_parts = _em_derivative_parts!(es.dBz, es.Bz, g, 1)
    exact_coefficients = nothing
    @inbounds for i = 1:n
        old_E, derivative, current = es.Ey[i], es.dBz[i], es.Jy[i]
        magnetic = c2 * derivative
        forcing = magnetic + current
        increment = dt * forcing
        candidate = old_E - increment
        if (
               magnetic_parts !== nothing ||
               _em_update_lost_range(c2, derivative, magnetic, forcing, increment, candidate)
           ) &&
           isfinite(old_E) &&
           (magnetic_parts !== nothing || isfinite(derivative)) &&
           isfinite(current)
            if exact_coefficients === nothing
                exact_coefficients = _em_exact_update_coefficients(es.c, dt)
            end
            exact_derivative =
                magnetic_parts === nothing ? _exact_rational(derivative) :
                _em_exact_derivative(magnetic_parts, i)
            candidate =
                _em_exact_electric_update(old_E, -exact_derivative, current, exact_coefficients)
        end
        es.Ey[i] = candidate
    end
    #   Ex^{n+1} = Ex^n − dt·Jx
    @inbounds for i = 1:n
        old_E, current = es.Ex[i], es.Jx[i]
        increment = dt * current
        candidate = old_E - increment
        if _em_increment_lost_range(current, increment, candidate) &&
           isfinite(old_E) &&
           isfinite(current)
            if exact_coefficients === nothing
                exact_coefficients = _em_exact_update_coefficients(es.c, dt)
            end
            candidate = _em_exact_electric_update(old_E, zero(T), current, exact_coefficients)
        end
        es.Ex[i] = candidate
    end

    return es
end

# Deposit total transverse Jy at each species' midpoint positions. The midpoints
# were parked in es.work (electrons) / es.worki (ions) by the drift; x currently
# holds x^{n+1}. We temporarily swap x ↔ midpoint for the deposit, then restore.
function _deposit_Jy_at_midpoints!(
    es::EMPIC1D{T},
    e::ParticleSet{1,T},
    ions::Union{Nothing,ParticleSet{1,T}},
) where {T}
    # electrons: park x^{n+1} in Bzp (unused now), use midpoint for deposit
    xe = e.x[1]
    @inbounds for p in eachindex(xe)
        es.Bzp[p] = xe[p]
        xe[p] = es.work[p]
    end
    # Eyp (per-particle) = work buffer; dBz (grid-sized, free until the Ey update
    # below) = deposit scratch. (Using the per-particle Exp here was a latent OOB
    # when the species count < grid size.)
    electron = try
        _add_transverse_Jy!(es, e, es.Eyp, es.dBz, false)
    finally
        @inbounds for p in eachindex(xe)
            xe[p] = es.Bzp[p]
        end
    end
    if es.mobile && ions !== nothing
        xi = ions.x[1]
        @inbounds for p in eachindex(xi)
            es.Bzi[p] = xi[p]
            xi[p] = es.worki[p]
        end
        ion = try
            # dBz is free grid scratch; accumulate onto the electron current.
            _add_transverse_Jy!(es, ions, es.Eyi, es.dBz, true)
        finally
            @inbounds for p in eachindex(xi)
                xi[p] = es.Bzi[p]
            end
        end
        _combine_midpoint_current!(es.Jy, es, e, ions, (es.work,), (es.worki,), 2, electron, ion)
    end
    _zero_nyquist_mode!(es.Jy)     # drop the unrepresentable Nyquist mode (see helper)
    return es.Jy
end

# ---------------------------------------------------------------- step

# One-time leapfrog priming for EMPIC1D. The loaded velocity is the physical v^0, but the Boris
# push expects v^{-1/2}; back it up a half step v^{-1/2} = v^0 − h·a^0, a^0 = qm(E + v×B) with
# B = (0,0,Bz). `h` is the species' half-step (dt_e/2 for electrons, dt_ion/2 for ions). Restores
# 2nd-order temporal accuracy (verified: rate 1 → 2 for both species).
@inline function _prime_empic1d!(v, Exp, Eyp, Bzp, qm, h, np, rel::Bool, c)
    vx, vy, vz = v
    if rel
        # the relativistic Boris push works in momentum u=γv, so back up MOMENTUM:
        # u^{-1/2}=γ^0 v^0 − h·qm(E+v^0×B^0), then v^{-1/2}=u^{-1/2}/γ(u^{-1/2}). Backing up v
        # directly (below) leaves an O(dt) momentum error ∝ the bulk drift → 1st order.
        c2 = c * c
        o = one(c)
        @inbounds for p = 1:np
            bz = Bzp[p]
            v0x = vx[p]
            v0y = vy[p]
            v0z = vz[p]
            βx = v0x / c
            βy = v0y / c
            βz = v0z / c
            β2 = βx * βx + βy * βy + βz * βz
            (isfinite(β2) && β2 < o) || throw(
                ArgumentError("_prime_empic1d!: input velocity must be finite and satisfy |v| < c"),
            )
            g0 = o / sqrt(o - β2)
            ax = qm * (Exp[p] + v0y * bz)                 # a = qm(E + v×B), B=(0,0,B_z)
            ay = qm * (Eyp[p] - v0x * bz)
            umx = g0 * v0x - h * ax
            umy = g0 * v0y - h * ay
            umz = g0 * v0z
            gm = sqrt(o + (umx * umx + umy * umy + umz * umz) / c2)
            vx[p] = umx / gm
            vy[p] = umy / gm
            vz[p] = umz / gm
        end
    else
        @inbounds for p = 1:np
            bz = Bzp[p]
            ux, uy = vx[p], vy[p]
            vx[p] = ux - h * qm * (Exp[p] + uy * bz)      # a_x = qm(Ex + v_y B_z)
            vy[p] = uy - h * qm * (Eyp[p] - ux * bz)      # a_y = qm(Ey − v_x B_z)
        end
    end
    return v
end

"""
    step_empic!(es, electrons[, ions], dt)

Advance the electromagnetic 1D PIC one full ion step `dt`. With `es.n_sub > 1`
the electrons and fields are advanced in `n_sub` substeps of `dt/n_sub`; the ions
(when `es.mobile`) receive one kick between their ion half levels. Time levels follow the
header of `full_pic.jl`. Changes in `dt` or `n_sub` retain the electron, ion,
and magnetic-field half levels from the preceding step; zero timesteps leave
this staggering unchanged. Returns `es`.
"""
function step_empic!(
    es::EMPIC1D{T},
    e::ParticleSet{1,T},
    ions::Union{Nothing,ParticleSet{1,T}},
    dt::Real,
) where {T}
    _require_empic_electrons(e)
    if es.mobile
        ions === nothing &&
            throw(ArgumentError("mobile=true requires an ion ParticleSet in step_empic!"))
        _require_empic_ions(ions, es.mi)
    end
    dtT = _validated_nonnegative_dt(T, dt; name = "step_empic!")
    ns = es.n_sub
    ns >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    iszero(dtT) && return es            # dt=0 no-op: do not consume the one-time priming
    dt_e = dtT / ns
    dt_e > zero(T) || throw(ArgumentError("dt/n_sub must be representable and positive"))
    next_time, next_step = _validated_step_clock(es.time[], es.step[], dt_e, ns)
    if es.relativistic
        _require_subluminal_velocities(e, es.c, "step_empic!")
        es.mobile && _require_subluminal_velocities(ions, es.c, "step_empic! ion species")
    end
    _ensure_electron_buffers!(es, nparticles(e))
    if es.mobile
        _ensure_ion_buffers!(es, nparticles(ions))
    end
    dt_previous = es.step[] == 0 ? dt_e : es.last_sub_dt[]
    dt_ion_push = es.step[] == 0 ? dtT : es.last_dt[] / 2 + dtT / 2
    # prime the leapfrog once (both species): loaded v is physical v^0 → v^{-1/2} for 2nd order.
    if es.step[] == 0
        g0 = es.g
        gather_scalar!(es.Exp, es.Ex, e, g0, es.shape)
        gather_scalar!(es.Eyp, es.Ey, e, g0, es.shape)
        gather_scalar!(es.Bzp, es.Bz, e, g0, es.shape)
        _prime_empic1d!(
            e.v,
            es.Exp,
            es.Eyp,
            es.Bzp,
            -one(T),
            dt_e / 2,
            nparticles(e),
            es.relativistic,
            es.c,
        )
        if es.mobile && ions !== nothing
            gather_scalar!(es.Exi, es.Ex, ions, g0, es.shape)
            gather_scalar!(es.Eyi, es.Ey, ions, g0, es.shape)
            gather_scalar!(es.Bzi, es.Bz, ions, g0, es.shape)
            # ions leapfrog at the full ion step dt_ion = ns·dt_e
            _prime_empic1d!(
                ions.v,
                es.Exi,
                es.Eyi,
                es.Bzi,
                ions.q / ions.m,
                dtT / 2,
                nparticles(ions),
                es.relativistic,
                es.c,
            )
        end
        # field priming: the seeded Bz is the physical Bz^0, but the Yee leapfrog carries Bz at
        # the half level Bz^{-1/2}. Back it up a half step: Bz^{-1/2} = Bz^0 + (dt_e/2) ∂xEy^0
        # (∂tBz = −∂xEy). Without this the integer-level fields (Ex, Ey) are only 1st-order in dt.
        electric_parts = _em_derivative_parts!(es.dEy, es.Ey, g0, 1)
        _em_magnetic_update!(es.Bz, es.Bz, es.dEy, electric_parts, dt_e, zero(T), true)
    end
    # The heavy ion is kicked ONCE (on the first substep) so that every substep drifts it with the
    # SAME post-kick velocity v^{n+1/2}; combined with the per-substep ion drift in _substep_em!,
    # its charge/current is then time-consistent on every substep — the whole mobile scheme is
    # 2nd-order for any n_sub. (Kicking on the middle substep instead left the pre-kick substeps
    # drifting with the stale v^{n−1/2}, an O(dt) current error → 1st-order for n_sub≥3.) The ion
    # drift now fills es.worki on every substep before the Jy deposit reads it, so no seed is needed.
    push_sub = 1
    @inbounds for s = 1:ns
        _substep_em!(es, e, ions, dt_e, s == push_sub, dt_previous, dt_ion_push)
        dt_previous = dt_e
    end
    es.time[] = next_time
    es.step[] = next_step
    es.last_dt[] = dtT
    es.last_sub_dt[] = dt_e
    return es
end

# Convenience 3-arg method: immobile-ion background (no ion species). Preserves
# the legacy call signature step_empic!(es, e, dt) bit-for-bit when n_sub=1.
step_empic!(es::EMPIC1D{T}, e::ParticleSet{1,T}, dt::Real) where {T} =
    step_empic!(es, e, nothing, dt)

# ---------------------------------------------------------------- diagnostics

function _em_field_energy_exact(E, B, c::T, dx) where {T}
    electric = zero(Rational{BigInt})
    magnetic = zero(Rational{BigInt})
    for component in E, value in component
        electric += _exact_rational(value)^2
    end
    for component in B, value in component
        magnetic += _exact_rational(value)^2
    end
    volume = prod(_exact_rational, dx)
    return _from_exact(T, (electric + _exact_rational(c)^2 * magnetic) * volume / 2)
end

# After normal-product checks, each exact contribution is at most eight times
# its rounded value (three operations for c²B²). The nonnegative rounded sum
# bounds every rounded contribution, so the exact sum is at most 8*count*sum.
# Powers of two bound that sum, the cell volume and the final factor of 1/2.
function _field_energy_upper_exponent(E, B, sum_squares, dx)
    count = Int128(0)
    for component in E
        count += length(component)
    end
    for component in B
        count += length(component)
    end
    count_exponent = 128 - leading_zeros(count - 1)
    upper = Int128(exponent(sum_squares)) + 3 + count_exponent
    for spacing in dx
        upper += Int128(exponent(abs(spacing))) + 1
    end
    return upper
end

function _em_field_energy_finish(E, B, c::T, dx, sum_squares, minimum_e2, minimum_b2) where {T}
    c2 = c * c
    lost = _spectral_product_lost_range(c, c, c2)
    # Contributions are nonnegative: the smallest square of a nonzero field
    # detects square underflow, and c² times the smallest magnetic square detects
    # magnetic-product underflow. Overflow remains visible in the total below.
    # Zero fields use an Inf sentinel and do not require exact arithmetic.
    lost |= minimum_e2 < floatmin(T) || minimum_b2 < floatmin(T)
    if isfinite(minimum_b2)
        lost |= _spectral_product_lost_range(c2, minimum_b2, c2 * minimum_b2)
    end
    volume = prod(dx)
    half = T(0.5) * sum_squares
    result = half * volume
    lost |=
        _spectral_product_lost_range(T(0.5), sum_squares, half) ||
        _spectral_product_lost_range(half, volume, result)
    # A spacing product can lose range and later return to normal range. Check
    # every contiguous subproduct: these cover either association in D <= 3.
    for first = 1:length(dx)
        partial = one(T)
        for last = first:length(dx)
            next = partial * dx[last]
            lost |= _spectral_product_lost_range(partial, dx[last], next)
            partial = next
        end
    end
    if T <: Union{Float16,Float32,Float64} && !lost && !iszero(sum_squares) && !iszero(volume)
        lost = _field_energy_upper_exponent(E, B, sum_squares, dx) > exponent(floatmax(T))
    end
    if lost &&
       isfinite(c) &&
       all(isfinite, dx) &&
       all(a -> all(isfinite, a), E) &&
       all(a -> all(isfinite, a), B)

        return _em_field_energy_exact(E, B, c, dx)
    end
    return result
end

"""
    em_field_energy(es)

Total electromagnetic field energy ∫ ½(ε0|E|² + |B|²/μ0) dx in normalized units
(ε0 = 1, μ0 = 1/c²): ½ Σ_i (Ex² + Ey² + c²Bz²) dx. The c² on Bz is the 1/μ0
factor (μ0 = 1/c² so that the wave speed is c).
Range-sensitive finite inputs retain the complete integral before rounding;
ordinary inputs keep their existing summation order.
"""
function em_field_energy(es::EMPIC1D{T}) where {T}
    c2 = es.c * es.c
    minimum_e2 = minimum_b2 = T(Inf)
    s = zero(T)
    @inbounds for i in eachindex(es.Ex)
        ex, ey, bz = es.Ex[i], es.Ey[i], es.Bz[i]
        ex2, ey2, bz2 = ex^2, ey^2, bz^2
        magnetic = c2 * bz2
        minimum_e2 =
            min(minimum_e2, ifelse(iszero(ex), T(Inf), ex2), ifelse(iszero(ey), T(Inf), ey2))
        minimum_b2 = min(minimum_b2, ifelse(iszero(bz), T(Inf), bz2))
        s += ex2 + ey2 + magnetic
    end
    return _em_field_energy_finish(
        (es.Ex, es.Ey),
        (es.Bz,),
        es.c,
        es.g.dx,
        s,
        minimum_e2,
        minimum_b2,
    )
end

"""
    charge_conservation_residual(es, dt) -> Real

Discrete continuity residual `max|ρ^{n+1} − ρ^n + dt ∂xJx| / scale` from the most
recent substep, evaluated on the spectrum the field solver actually evolves. Here
ρ is the TOTAL charge density (ions + electrons) and `dt` is the SUBSTEP size
`dt_full/n_sub` used by the last `step_empic!`.

The longitudinal current is built so that `ik·Ĵx = −(ρ̂^{n+1}−ρ̂^n)/dt` exactly,
which makes `∂tρ + ∂xJx = 0` hold to roundoff for every represented Fourier mode.
On an even grid the spectral first-derivative multiplier `ik` zeros the Nyquist
mode by construction (it has no representable sine partner — see `spectral.jl`),
so `∂xJx` carries no Nyquist content and that single mode of `Δρ` is structurally
outside the operator's range. The residual is therefore measured on the
representable (non-Nyquist) spectrum — the exact space on which the discrete
conservation law is defined — by forming `Δρ + dt ∂xJx` after dropping the
Nyquist mode of `Δρ` (`∂xJx` already has none). `scale = max(|Δρ_repr|)` or 1.
Call after `step_empic!`.

Range-sensitive inputs retain the complete density and timestep-weighted
current expression. Spectral input validation and FFT accuracy limits still
apply, including loss of relative accuracy for ill-conditioned cancellations.
"""
function charge_conservation_residual(es::EMPIC1D{T}, dt::Real) where {T}
    g = es.g
    n = g.n[1]
    dtT = T(dt)
    # Δρ projected onto the representable spectrum: fft, zero Nyquist, ifft.
    cb = similar(es.Jx, Complex{T})
    bounds = _em_density_transform_bounds(T, n)
    lost = false
    @inbounds for i = 1:n
        delta = es.rho_np1[i] - es.rho_n[i]
        cb[i] = Complex{T}(delta)
        lost |= _em_density_transform_lost_range(delta, bounds)
    end
    g.plan * cb
    if iseven(n)
        cb[n÷2+1] = zero(Complex{T})        # drop the unrepresentable Nyquist mode
    end
    g.iplan * cb
    dJx = similar(es.Jx)
    deriv!(dJx, es.Jx, g, 1)                     # already Nyquist-free
    return _em_continuity_finish(g, es.rho_n, es.rho_np1, (es.Jx,), dtT, cb, dJx, lost)
end

# ---------------------------------------------------------------- dimension-parametric EM PIC

"""
    EMPIC(g, Nparticles; n0=1.0, c=5.0, shape=CIC(), relativistic=false,
          mobile=false, mi=1836.0, n_sub=1)

Dimension-parametric periodic electromagnetic PIC on `FourierGrid{D}` for
`D = 1, 2, 3`. Electrons are kinetic (`q=-1`, `m=1`) with either a uniform
immobile positive background of density `n0` or, when `mobile=true`, a second
kinetic ion species (`q=+1`) supplied to [`init_empic!`](@ref) and
[`step_empic!`](@ref). Fields are collocated spectral arrays with `E` and `B`
stored as 3-component tuples. `n_sub` (≥1) splits each full step into electron
and Maxwell substeps of `dt/n_sub`. Mobile ions receive one velocity kick per
full step and drift during every electron substep, contributing charge and
current at each substep's actual positions and midpoints.

The mover uses a leapfrog Maxwell update (here `dt` is the electron substep):
`B^{n} = B^{n-1/2} - (dt/2) curl(E^n)` for the particle push,
`B^{n+1/2} = B^{n-1/2} - dt curl(E^n)`, and
`E^{n+1} = E^n + dt(c^2 curl(B^{n+1/2}) - J^{n+1/2})`.
The midpoint particle current is spectrally corrected only in its longitudinal
component so that `Δρ/dt + div(J) = 0` holds to roundoff on the represented
Fourier modes; the pure-Nyquist mode combinations (`k² = 0` but not DC), which
are outside the spectral operators' range, are zeroed in all components.
"""
mutable struct EMPIC{D,T,G,SH<:ShapeFunction}
    g::G
    n0::T
    c::T
    shape::SH
    relativistic::Bool
    mobile::Bool
    mi::T
    n_sub::Int
    E::NTuple{3,Array{T,D}}
    B::NTuple{3,Array{T,D}}
    Bc::NTuple{3,Array{T,D}}
    J::NTuple{3,Array{T,D}}
    Jraw::NTuple{3,Array{T,D}}
    rho_n::Array{T,D}
    rho_np1::Array{T,D}
    ne::Array{T,D}
    curlE::NTuple{3,Array{T,D}}
    curlB::NTuple{3,Array{T,D}}
    Ep::NTuple{3,Vector{T}}
    Bp::NTuple{3,Vector{T}}
    mide::NTuple{D,Vector{T}}
    worke::Vector{T}
    Epi::NTuple{3,Vector{T}}
    Bpi::NTuple{3,Vector{T}}
    midi::NTuple{D,Vector{T}}
    worki::Vector{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    last_dt::Base.RefValue{T}       # preceding full ion drift interval
    last_sub_dt::Base.RefValue{T}   # preceding electron/field drift interval
end

function EMPIC(
    g::FourierGrid{D,T},
    Nparticles::Integer;
    n0 = 1.0,
    c = 5.0,
    shape::ShapeFunction = CIC(),
    relativistic::Bool = false,
    mobile::Bool = false,
    mi = 1836.0,
    n_sub::Integer = 1,
) where {D,T}
    1 <= D <= 3 || throw(ArgumentError("EMPIC supports D = 1, 2, or 3"))
    n_sub >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    Np = _particle_length(Nparticles)
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    cT = _require_finite_positive_real("c", c, T)
    miT = _require_finite_positive_real("mi", mi, T)
    grid_tuple() = ntuple(_ -> zeros(T, g.n), 3)
    particle_tuple() = ntuple(_ -> zeros(T, Np), 3)
    midpoint_tuple() = ntuple(_ -> zeros(T, Np), D)
    return EMPIC{D,T,typeof(g),typeof(shape)}(
        g,
        n0T,
        cT,
        shape,
        relativistic,
        mobile,
        miT,
        Int(n_sub),
        grid_tuple(),
        grid_tuple(),
        grid_tuple(),
        grid_tuple(),
        grid_tuple(),
        zeros(T, g.n),
        zeros(T, g.n),
        zeros(T, g.n),
        grid_tuple(),
        grid_tuple(),
        particle_tuple(),
        particle_tuple(),
        midpoint_tuple(),
        zeros(T, Np),
        ntuple(_ -> T[], 3),
        ntuple(_ -> T[], 3),
        ntuple(_ -> T[], D),
        T[],
        Ref(zero(T)),
        Ref(0),
        Ref(zero(T)),
        Ref(zero(T)),
    )
end

@inline function _ensure_empic_ion_buffers!(es::EMPIC{D,T}, ni::Integer) where {D,T}
    ni >= 0 || throw(ArgumentError("ion particle count must be nonnegative"))
    for c = 1:3
        resize!(es.Epi[c], ni)
        resize!(es.Bpi[c], ni)
    end
    for d = 1:D
        resize!(es.midi[d], ni)
    end
    resize!(es.worki, ni)
    return es
end

# Electron gather buffers are sized to Np at construction; re-size them to the current electron
# count each step so a growing electron population (e.g. ionize_mcc! secondaries) composes with
# the push, symmetric with the ion buffers above. resize! to the same length is a no-op.
@inline function _ensure_empic_electron_buffers!(es::EMPIC{D,T}, ne::Integer) where {D,T}
    ne >= 0 || throw(ArgumentError("electron particle count must be nonnegative"))
    for c = 1:3
        resize!(es.Ep[c], ne)
        resize!(es.Bp[c], ne)
    end
    for d = 1:D
        resize!(es.mide[d], ne)
    end
    resize!(es.worke, ne)
    return es
end

@inline _empic_component_buffer(g, c::Int) = c == 1 ? g.cbuf : c == 2 ? g.tbuf : g.abuf

function _enforce_gauss!(es::EMPIC{D,T}) where {D,T}
    g = es.g
    rhohat = g.sbuf
    rhohat .= es.rho_n
    g.plan * rhohat
    buffers = (g.cbuf, g.tbuf, g.abuf)
    for c = 1:3
        buffers[c] .= es.E[c]
        g.plan * buffers[c]
    end
    @inbounds for I in CartesianIndices(rhohat)
        k = _spectral_wavevector(g, I)
        if !_spectral_mode_is_null(k)
            original = ntuple(d -> buffers[d][I], Val(D))
            corrected = _spectral_project_mode(original, k, rhohat[I], one(T))
            for d = 1:D
                buffers[d][I] = corrected[d]
            end
        end
    end
    for c = 1:3
        g.iplan * buffers[c]
        es.E[c] .= real.(buffers[c])
    end
    return es.E
end

function _charge_density!(
    rho::Array{T,D},
    ne::Array{T,D},
    es::EMPIC{D,T},
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}},
) where {D,T}
    return _empic_charge_density!(rho, ne, es, e, ions)
end

function _add_species_current!(
    J::NTuple{3,Array{T,D}},
    es::EMPIC{D,T},
    s::ParticleSet{D,T},
    work::Vector{T},
    scratch::Array{T,D},
    accumulate::Bool,
) where {D,T}
    length(work) == nparticles(s) ||
        throw(DimensionMismatch("current work length must equal particle count"))
    return ntuple(Val(3)) do c
        @. work = s.weight * s.v[c]
        _validate_deposition_inputs(scratch, s, work, es.g)
        factors = (s.weight, s.v[c])
        state = _deposit_scalar_kernel!(scratch, s, factors, es.g, es.shape, Val(true))
        return _finish_moment!(
            J[c],
            scratch,
            s,
            factors,
            es.g,
            es.shape,
            s.q,
            state,
            Val(:prescale),
            accumulate,
        )
    end
end

@inline function _swap_positions_with_midpoints!(
    s::ParticleSet{D,T},
    mid::NTuple{D,Vector{T}},
) where {D,T}
    for d = 1:D
        length(mid[d]) == length(s.x[d]) ||
            throw(DimensionMismatch("midpoint buffer length must equal particle count"))
    end
    @inbounds for d = 1:D
        x = s.x[d]
        xm = mid[d]
        for p in eachindex(x)
            tmp = x[p]
            x[p] = xm[p]
            xm[p] = tmp
        end
    end
    return s
end

function _deposit_current_at_midpoints!(
    es::EMPIC{D,T},
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}},
) where {D,T}
    _swap_positions_with_midpoints!(e, es.mide)
    electron = try
        _add_species_current!(es.Jraw, es, e, es.worke, es.ne, false)
    finally
        _swap_positions_with_midpoints!(e, es.mide)
    end
    if es.mobile && ions !== nothing
        _swap_positions_with_midpoints!(ions, es.midi)
        ion = try
            _add_species_current!(es.Jraw, es, ions, es.worki, es.ne, true)
        finally
            _swap_positions_with_midpoints!(ions, es.midi)
        end
        for c = 1:3
            _combine_midpoint_current!(
                es.Jraw[c],
                es,
                e,
                ions,
                es.mide,
                es.midi,
                c,
                electron[c],
                ion[c],
            )
        end
    end
    return es.Jraw
end

function _correct_current_continuity!(es::EMPIC{D,T}, dt::T) where {D,T}
    g = es.g
    buffers = (g.cbuf, g.tbuf, g.abuf)
    for c = 1:3
        buffers[c] .= es.Jraw[c]
        g.plan * buffers[c]
    end
    delta = g.sbuf
    @inbounds for I in CartesianIndices(delta)
        delta[I] = Complex{T}(es.rho_np1[I] - es.rho_n[I])
    end
    g.plan * delta
    dc = first(CartesianIndices(delta))          # the true k=0 mode (all indices 1)
    @inbounds for I in CartesianIndices(delta)
        k = _spectral_wavevector(g, I)
        if !_spectral_mode_is_null(k)
            original = ntuple(d -> buffers[d][I], Val(D))
            corrected = _spectral_project_mode(original, k, -delta[I], dt)
            for d = 1:D
                buffers[d][I] = corrected[d]
            end
        elseif I != dc
            # All-zero k covers the true DC mode AND every pure DC/Nyquist index
            # combination (kvec zeros the per-axis Nyquist entry). Only DC is
            # physical (net current, kept from the raw deposit). The Nyquist
            # combinations are structurally outside the spectral operators'
            # range: every derivative vanishes there, so raw deposition noise
            # would integrate into E as an undamped random walk (a growing
            # grid-Nyquist sawtooth felt by particles through the gather).
            # Zero ALL current components there, exactly as _esirkepov_Jx!
            # (longitudinal) and _deposit_Jy_at_midpoints! (transverse) do in
            # EMPIC1D. Mixed modes (e.g. Nyquist in x with k_y ≠ 0)
            # are deliberately left alone: the axis-Nyquist J component there
            # is transverse to the effective wavevector — raw, like every
            # transverse current — and the curl coupling on the remaining axes
            # gives E a restoring term, so it stays bounded (verified in 2D:
            # no secular (Nyq_x,k_y) growth over 500 steps, unlike the
            # pre-fix pure-Nyquist modes, which grew without bound).
            for c = 1:3
                buffers[c][I] = zero(Complex{T})
            end
        end
    end
    for c = 1:3
        g.iplan * buffers[c]
        es.J[c] .= real.(buffers[c])
    end
    return es.J
end

@inline function _push_velocities!(
    es::EMPIC{D,T},
    s::ParticleSet{D,T},
    Ep::NTuple{3,Vector{T}},
    Bp::NTuple{3,Vector{T}},
    dt::T,
) where {D,T}
    qm = s.q / s.m
    vx, vy, vz = s.v
    if es.relativistic
        cc = es.c
        @inbounds for p in eachindex(s.weight)
            nx, ny, nz = _boris_rel(
                vx[p],
                vy[p],
                vz[p],
                Ep[1][p],
                Ep[2][p],
                Ep[3][p],
                Bp[1][p],
                Bp[2][p],
                Bp[3][p],
                qm,
                dt,
                cc,
            )
            vx[p] = nx
            vy[p] = ny
            vz[p] = nz
        end
    else
        @inbounds for p in eachindex(s.weight)
            nx, ny, nz = boris_kick(
                vx[p],
                vy[p],
                vz[p],
                Ep[1][p],
                Ep[2][p],
                Ep[3][p],
                Bp[1][p],
                Bp[2][p],
                Bp[3][p],
                qm,
                dt,
            )
            vx[p] = nx
            vy[p] = ny
            vz[p] = nz
        end
    end
    return s
end

@inline function _drift_record_mid!(
    s::ParticleSet{D,T},
    mid::NTuple{D,Vector{T}},
    dt::T,
    L::NTuple{D,T},
) where {D,T}
    @inbounds for d = 1:D
        x = s.x[d]
        v = s.v[d]
        xm = mid[d]
        length(xm) == length(x) ||
            throw(DimensionMismatch("midpoint buffer length must equal particle count"))
        for p in eachindex(x)
            xmid = x[p] + (dt / 2) * v[p]
            x[p] += dt * v[p]
            xm[p] = mod(xmid, L[d])
        end
    end
    return s
end

"""
    init_empic!(es::EMPIC, electrons[, ions])

Initialize a dimension-parametric EM PIC state. Existing transverse electric
field content is preserved; only the longitudinal component is corrected so
`div(E) = rho` for the loaded particles. `B` keeps the caller-provided initial
values.
"""
function init_empic!(
    es::EMPIC{D,T},
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}} = nothing,
) where {D,T}
    _require_empic_electrons(e)
    _ensure_empic_electron_buffers!(es, nparticles(e))
    if es.mobile
        ions === nothing &&
            throw(ArgumentError("mobile=true requires an ion ParticleSet in init_empic!"))
        _require_empic_ions(ions, es.mi)
        _ensure_empic_ion_buffers!(es, nparticles(ions))
    end
    _charge_density!(es.rho_n, es.ne, es, e, ions)
    _enforce_gauss!(es)
    es.time[] = zero(T)                 # (re)start at t=0; step==0 triggers leapfrog priming
    es.step[] = 0
    es.last_dt[] = zero(T)
    es.last_sub_dt[] = zero(T)
    return es
end

function _substep_empic!(
    es::EMPIC{D,T},
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}},
    dt::T,
    dt_previous::T,
    push_ions::Bool,
    dt_ion_push::T,
) where {D,T}
    g = es.g
    c2 = es.c * es.c
    kick_dt = dt_previous / 2 + dt / 2
    _charge_density!(es.rho_n, es.ne, es, e, ions)

    electric_parts = _em_curl_parts!(es.curlE, es.E, g)
    for c = 1:3
        parts = electric_parts === nothing ? nothing : electric_parts[c]
        _em_magnetic_update!(es.Bc[c], es.B[c], es.curlE[c], parts, dt_previous, zero(T))
    end

    gather_vector!(es.Ep, es.E, e, g, es.shape)
    gather_vector!(es.Bp, es.Bc, e, g, es.shape)
    _push_velocities!(es, e, es.Ep, es.Bp, kick_dt)
    _drift_record_mid!(e, es.mide, dt, g.L)
    apply_periodic!(e, ntuple(_ -> zero(T), D), g.L)

    if es.mobile && ions !== nothing
        if push_ions
            gather_vector!(es.Epi, es.E, ions, g, es.shape)
            gather_vector!(es.Bpi, es.Bc, ions, g, es.shape)
            _push_velocities!(es, ions, es.Epi, es.Bpi, dt_ion_push)
        end
        # Both species must drift and deposit at every fast substep, even when
        # the ion velocity stays at its full-step half level.
        _drift_record_mid!(ions, es.midi, dt, g.L)
        apply_periodic!(ions, ntuple(_ -> zero(T), D), g.L)
    end

    _deposit_current_at_midpoints!(es, e, ions)
    _charge_density!(es.rho_np1, es.ne, es, e, ions)
    _correct_current_continuity!(es, dt)

    for c = 1:3
        parts = electric_parts === nothing ? nothing : electric_parts[c]
        _em_magnetic_update!(es.B[c], es.B[c], es.curlE[c], parts, dt_previous, dt)
    end
    magnetic_parts = _em_curl_parts!(es.curlB, es.B, g)
    exact_coefficients = nothing
    @inbounds for c = 1:3, I in eachindex(es.E[c])
        old_E, derivative, current = es.E[c][I], es.curlB[c][I], es.J[c][I]
        magnetic = c2 * derivative
        forcing = magnetic - current
        increment = dt * forcing
        candidate = old_E + increment
        if (
               magnetic_parts !== nothing ||
               _em_update_lost_range(c2, derivative, magnetic, forcing, increment, candidate)
           ) &&
           isfinite(old_E) &&
           (magnetic_parts !== nothing || isfinite(derivative)) &&
           isfinite(current)
            if exact_coefficients === nothing
                exact_coefficients = _em_exact_update_coefficients(es.c, dt)
            end
            exact_derivative =
                magnetic_parts === nothing ? _exact_rational(derivative) :
                _em_exact_derivative(magnetic_parts[c], I)
            candidate =
                _em_exact_electric_update(old_E, exact_derivative, current, exact_coefficients)
        end
        es.E[c][I] = candidate
    end
    return es
end

"""
    step_empic!(es::EMPIC, electrons[, ions], dt)

Advance the dimension-parametric EM PIC state by one full step, split into
`es.n_sub` electron and Maxwell substeps of `dt/es.n_sub`. The optional ion
species is required when `mobile=true` and ignored otherwise. Ions receive one
velocity kick per full step and drift during every substep. Changes in `dt` or
`n_sub` retain the electron, ion, and magnetic-field half levels from the
preceding step; zero timesteps leave this staggering unchanged.
"""
function step_empic!(
    es::EMPIC{D,T},
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}},
    dt::Real,
) where {D,T}
    _require_empic_electrons(e)
    if es.mobile
        ions === nothing &&
            throw(ArgumentError("mobile=true requires an ion ParticleSet in step_empic!"))
        _require_empic_ions(ions, es.mi)
    end
    dtT = _validated_nonnegative_dt(T, dt; name = "step_empic!")
    ns = es.n_sub
    ns >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    iszero(dtT) && return es            # dt=0 no-op: do not consume the one-time priming
    dt_e = dtT / ns
    dt_e > zero(T) || throw(ArgumentError("dt/n_sub must be representable and positive"))
    next_time, next_step = _validated_step_clock(es.time[], es.step[], dt_e, ns)
    if es.relativistic
        _require_subluminal_velocities(e, es.c, "step_empic!")
        es.mobile && _require_subluminal_velocities(ions, es.c, "step_empic! ion species")
    end
    _ensure_empic_electron_buffers!(es, nparticles(e))
    if es.mobile
        _ensure_empic_ion_buffers!(es, nparticles(ions))
    end
    # Prime physical initial v and B onto their half levels: dt_e/2 for the
    # electrons and B, and dtT/2 for ions. See the EMPIC1D header note.
    if es.step[] == 0
        g0 = es.g
        h = dt_e / 2
        gather_vector!(es.Ep, es.E, e, g0, es.shape)
        gather_vector!(es.Bp, es.B, e, g0, es.shape)
        _prime_leapfrog!(e.v, es.Ep, es.Bp, -one(T), h, nparticles(e), es.relativistic, es.c)
        if es.mobile && ions !== nothing
            gather_vector!(es.Epi, es.E, ions, g0, es.shape)
            gather_vector!(es.Bpi, es.B, ions, g0, es.shape)
            _prime_leapfrog!(
                ions.v,
                es.Epi,
                es.Bpi,
                ions.q / ions.m,
                dtT / 2,
                nparticles(ions),
                es.relativistic,
                es.c,
            )
        end
        electric_parts = _em_curl_parts!(es.curlE, es.E, g0)
        for c = 1:3
            parts = electric_parts === nothing ? nothing : electric_parts[c]
            _em_magnetic_update!(es.B[c], es.B[c], es.curlE[c], parts, dt_e, zero(T), true)
        end
    end
    dt_previous = es.step[] == 0 ? dt_e : es.last_sub_dt[]
    dt_ion_push = es.step[] == 0 ? dtT : es.last_dt[] / 2 + dtT / 2
    for substep = 1:ns
        _substep_empic!(es, e, ions, dt_e, dt_previous, substep == 1, dt_ion_push)
        dt_previous = dt_e
    end
    es.time[] = next_time
    es.step[] = next_step
    es.last_dt[] = dtT
    es.last_sub_dt[] = dt_e
    return es
end

step_empic!(es::EMPIC{D,T}, e::ParticleSet{D,T}, dt::Real) where {D,T} =
    step_empic!(es, e, nothing, dt)

function em_field_energy(es::EMPIC{D,T}) where {D,T}
    c2 = es.c * es.c
    minimum_e2 = minimum_b2 = T(Inf)
    s = zero(T)
    @inbounds for c = 1:3, I in eachindex(es.E[c])
        electric, magnetic = es.E[c][I], es.B[c][I]
        e2, b2 = electric^2, magnetic^2
        weighted = c2 * b2
        minimum_e2 = min(minimum_e2, ifelse(iszero(electric), T(Inf), e2))
        minimum_b2 = min(minimum_b2, ifelse(iszero(magnetic), T(Inf), b2))
        s += e2 + weighted
    end
    return _em_field_energy_finish(es.E, es.B, es.c, es.g.dx, s, minimum_e2, minimum_b2)
end

"""
    charge_conservation_residual(es::EMPIC, dt)

Continuity residual on the represented Fourier modes from the most recent
electron substep. Pass its duration `dt_full/es.n_sub`, including when ions
are mobile and receive one kick over the full step. The true zero-frequency
mode is retained so a change in total charge contributes to the residual;
only non-DC modes annihilated by every spectral derivative are excluded.
Range-sensitive inputs retain the complete density and timestep-weighted
current expression; spectral input validation and FFT accuracy limits still apply.
"""
function charge_conservation_residual(es::EMPIC{D,T}, dt::Real) where {D,T}
    g = es.g
    dtT = T(dt)
    dρhat = g.cbuf
    bounds = _em_density_transform_bounds(T, length(dρhat))
    lost = false
    @inbounds for I in CartesianIndices(dρhat)
        delta = es.rho_np1[I] - es.rho_n[I]
        dρhat[I] = Complex{T}(delta)
        lost |= _em_density_transform_lost_range(delta, bounds)
    end
    g.plan * dρhat
    dc = first(CartesianIndices(dρhat))
    @inbounds for I in CartesianIndices(dρhat)
        I != dc &&
            _spectral_mode_is_null(_spectral_wavevector(g, I)) &&
            (dρhat[I] = zero(Complex{T}))
    end
    g.iplan * dρhat
    dρrepr = similar(es.rho_n)
    @inbounds for I in eachindex(dρrepr)
        dρrepr[I] = real(dρhat[I])
    end

    dJ = similar(es.rho_n)
    divergence!(dJ, es.J, g)
    return _em_continuity_finish(g, es.rho_n, es.rho_np1, es.J, dtT, dρrepr, dJ, lost)
end
