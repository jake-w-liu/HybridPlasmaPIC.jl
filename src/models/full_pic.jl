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

# Number-density-weighted charge density of ONE species accumulated into `rho`
# (rho += q · n_species). Electrons (q=−1) subtract, ions (q=+1) add.
@inline function _add_species_charge!(
    rho::Vector{T},
    ne::Vector{T},
    es::EMPIC1D{T},
    s::ParticleSet{1,T},
) where {T}
    density!(ne, s, es.g, es.shape)
    q = s.q
    @inbounds for i in eachindex(rho)
        rho[i] += q * ne[i]
    end
    return rho
end

# Total charge density ρ at the CURRENT particle positions.
#   immobile-ion mode : ρ = n0 + q_e·n_e = n0 − n_e   (uniform +n0 background)
#   mobile-ion mode   : ρ = q_i·n_i + q_e·n_e = n_i − n_e
function _charge_density!(
    rho::Vector{T},
    ne::Vector{T},
    es::EMPIC1D{T},
    e::ParticleSet{1,T},
    ions::Union{Nothing,ParticleSet{1,T}},
) where {T}
    if es.mobile && ions !== nothing
        fill!(rho, zero(T))
        _add_species_charge!(rho, ne, es, ions)   # + q_i n_i
        _add_species_charge!(rho, ne, es, e)       # + q_e n_e
    else
        density!(ne, e, es.g, es.shape)            # electron number density n_e
        @inbounds for i in eachindex(rho)
            rho[i] = es.n0 - ne[i]                 # n0 − n_e (immobile ion bg)
        end
    end
    return rho
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
    Jx0 = _species_Jx0(e, g.L[1])
    if es.mobile && ions !== nothing
        Jx0 += _species_Jx0(ions, g.L[1])
    end
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

@inline function _species_Jx0(s::ParticleSet{1,T}, L::T) where {T}
    vx = s.v[1]
    w = s.weight
    acc = zero(T)
    @inbounds for p in eachindex(w)
        acc += w[p] * vx[p]
    end
    return s.q * acc / L
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
    ΔV = prod(es.g.dx)
    @inbounds @. work = s.weight * s.v[2]
    deposit_scalar!(scratch, s, work, es.g, es.shape)
    f = s.q / ΔV
    if accumulate
        @inbounds for i in eachindex(es.Jy)
            es.Jy[i] += f * scratch[i]
        end
    else
        @inbounds for i in eachindex(es.Jy)
            es.Jy[i] = f * scratch[i]
        end
    end
    return es.Jy
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
    @inbounds for i = 1:n
        acc += s * f[i]
        s = -s
    end
    acc /= n
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
    # guard against superluminal input (β²≥1 would give NaN γ and poison the
    # fields): clamp β² to just below 1 so γ stays finite.
    β2 = (vx * vx + vy * vy + vz * vz) / c2
    β2 >= one(T) && (β2 = one(T) - eps(T))
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
    return es
end

# ---------------------------------------------------------------- substep

# One field+electron substep of size dt. Pushes the electrons `e` always; pushes
# the ions `ions` only when `push_ions` is true (heavy ions are pushed once per
# full ion step, on the substep straddling their midpoint). When mobile, the ions
# always contribute their charge (rho snapshots) and current to the field solve.
function _substep_em!(
    es::EMPIC1D{T},
    e::ParticleSet{1,T},
    ions::Union{Nothing,ParticleSet{1,T}},
    dt::T,
    push_ions::Bool,
) where {T}
    g = es.g
    n = g.n[1]
    L = g.L[1]
    c2 = es.c * es.c
    z = zero(T)

    # snapshot ρ^n (from x^n) for the charge-conserving current
    _charge_density!(es.rho_n, es.ne, es, e, ions)

    # (1) time-centered Bz^n = Bz^{n−1/2} − (dt/2)∂xEy^n  (for the particle push)
    deriv!(es.dEy, es.Ey, g, 1)
    @inbounds for i = 1:n
        es.Bzc[i] = es.Bz[i] - (dt / 2) * es.dEy[i]
    end

    # (2) gather fields at time n and push electrons
    gather_scalar!(es.Exp, es.Ex, e, g, es.shape)
    gather_scalar!(es.Eyp, es.Ey, e, g, es.shape)
    gather_scalar!(es.Bzp, es.Bzc, e, g, es.shape)
    _push_velocities!(es, e, es.Exp, es.Eyp, es.Bzp, dt)
    _drift_record_mid!(e, es.work, dt, L)
    apply_periodic!(e, (z,), (L,))

    # ion VELOCITY kick: once per full ion step, on the straddling substep so the field is
    # time-centred; the heavy ion resolves the slow scale, so one kick per dt_ion = n_sub·dt is enough.
    if es.mobile && ions !== nothing && push_ions
        gather_scalar!(es.Exi, es.Ex, ions, g, es.shape)
        gather_scalar!(es.Eyi, es.Ey, ions, g, es.shape)
        gather_scalar!(es.Bzi, es.Bzc, ions, g, es.shape)
        _push_velocities!(es, ions, es.Exi, es.Eyi, es.Bzi, T(es.n_sub) * dt)
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

    # (4) Bz full leapfrog advance: Bz^{n+1/2} = Bz^{n−1/2} − dt ∂xEy^n
    @inbounds for i = 1:n
        es.Bz[i] -= dt * es.dEy[i]
    end

    # (5) E advance:
    #   Ey^{n+1} = Ey^n − dt(c²∂xBz^{n+1/2} + Jy)
    deriv!(es.dBz, es.Bz, g, 1)
    @inbounds for i = 1:n
        es.Ey[i] -= dt * (c2 * es.dBz[i] + es.Jy[i])
    end
    #   Ex^{n+1} = Ex^n − dt·Jx
    @inbounds for i = 1:n
        es.Ex[i] -= dt * es.Jx[i]
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
    _add_transverse_Jy!(es, e, es.Eyp, es.dBz, false)
    @inbounds for p in eachindex(xe)
        xe[p] = es.Bzp[p]
    end
    if es.mobile && ions !== nothing
        xi = ions.x[1]
        @inbounds for p in eachindex(xi)
            es.Bzi[p] = xi[p]
            xi[p] = es.worki[p]
        end
        _add_transverse_Jy!(es, ions, es.Eyi, es.dBz, true)  # dBz = free grid scratch; accumulate onto Jy
        @inbounds for p in eachindex(xi)
            xi[p] = es.Bzi[p]
        end
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
            β2 = (v0x * v0x + v0y * v0y + v0z * v0z) / c2
            β2 >= o && (β2 = o - eps(c))
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
(when `es.mobile`) are pushed once over the full `dt`. Time levels follow the
header of `empic1d.jl`. Returns `es`.
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
    iszero(dtT) && return es            # dt=0 no-op: do not consume the one-time priming
    _ensure_electron_buffers!(es, nparticles(e))
    if es.mobile
        _ensure_ion_buffers!(es, nparticles(ions))
    end
    ns = es.n_sub
    dt_e = dtT / ns
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
                (T(ns) * dt_e) / 2,
                nparticles(ions),
                es.relativistic,
                es.c,
            )
        end
        # field priming: the seeded Bz is the physical Bz^0, but the Yee leapfrog carries Bz at
        # the half level Bz^{-1/2}. Back it up a half step: Bz^{-1/2} = Bz^0 + (dt_e/2) ∂xEy^0
        # (∂tBz = −∂xEy). Without this the integer-level fields (Ex, Ey) are only 1st-order in dt.
        deriv!(es.dEy, es.Ey, g0, 1)
        @inbounds for i = 1:g0.n[1]
            es.Bz[i] += (dt_e / 2) * es.dEy[i]
        end
    end
    # The heavy ion is kicked ONCE (on the first substep) so that every substep drifts it with the
    # SAME post-kick velocity v^{n+1/2}; combined with the per-substep ion drift in _substep_em!,
    # its charge/current is then time-consistent on every substep — the whole mobile scheme is
    # 2nd-order for any n_sub. (Kicking on the middle substep instead left the pre-kick substeps
    # drifting with the stale v^{n−1/2}, an O(dt) current error → 1st-order for n_sub≥3.) The ion
    # drift now fills es.worki on every substep before the Jy deposit reads it, so no seed is needed.
    push_sub = 1
    @inbounds for s = 1:ns
        _substep_em!(es, e, ions, dt_e, s == push_sub)
        es.time[] += dt_e
    end
    es.step[] += 1
    return es
end

# Convenience 3-arg method: immobile-ion background (no ion species). Preserves
# the legacy call signature step_empic!(es, e, dt) bit-for-bit when n_sub=1.
step_empic!(es::EMPIC1D{T}, e::ParticleSet{1,T}, dt::Real) where {T} =
    step_empic!(es, e, nothing, dt)

# ---------------------------------------------------------------- diagnostics

"""
    em_field_energy(es)

Total electromagnetic field energy ∫ ½(ε0|E|² + |B|²/μ0) dx in normalized units
(ε0 = 1, μ0 = 1/c²): ½ Σ_i (Ex² + Ey² + c²Bz²) dx. The c² on Bz is the 1/μ0
factor (μ0 = 1/c² so that the wave speed is c).
"""
function em_field_energy(es::EMPIC1D{T}) where {T}
    dx = es.g.dx[1]
    c2 = es.c * es.c
    s = zero(T)
    @inbounds for i in eachindex(es.Ex)
        s += es.Ex[i]^2 + es.Ey[i]^2 + c2 * es.Bz[i]^2
    end
    return T(0.5) * s * dx
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
"""
function charge_conservation_residual(es::EMPIC1D{T}, dt::Real) where {T}
    g = es.g
    n = g.n[1]
    dtT = T(dt)
    # Δρ projected onto the representable spectrum: fft, zero Nyquist, ifft.
    cb = similar(es.Jx, Complex{T})
    @inbounds for i = 1:n
        cb[i] = Complex{T}(es.rho_np1[i] - es.rho_n[i])
    end
    fft!(cb)
    if iseven(n)
        cb[n÷2+1] = zero(Complex{T})        # drop the unrepresentable Nyquist mode
    end
    ifft!(cb)
    dJx = similar(es.Jx)
    deriv!(dJx, es.Jx, g, 1)                     # already Nyquist-free
    rmax = zero(T)
    scale = zero(T)
    @inbounds for i = 1:n
        dρ = real(cb[i])
        r = dρ + dtT * dJx[i]
        rmax = max(rmax, abs(r))
        scale = max(scale, abs(dρ))
    end
    return rmax / (scale > 0 ? scale : one(T))
end

# ---------------------------------------------------------------- dimension-parametric EM PIC

"""
    EMPIC(g, Nparticles; n0=1.0, c=5.0, shape=CIC(), relativistic=false,
          mobile=false, mi=1836.0)

Dimension-parametric periodic electromagnetic PIC on `FourierGrid{D}` for
`D = 1, 2, 3`. Electrons are kinetic (`q=-1`, `m=1`) with either a uniform
immobile positive background of density `n0` or, when `mobile=true`, a second
kinetic ion species (`q=+1`) supplied to [`init_empic!`](@ref) and
[`step_empic!`](@ref). Fields are collocated spectral arrays with `E` and `B`
stored as 3-component tuples.

The mover uses a leapfrog Maxwell update:
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
) where {D,T}
    1 <= D <= 3 || throw(ArgumentError("EMPIC supports D = 1, 2, or 3"))
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

@inline function _empic_divergence_hat(
    buffers,
    g::FourierGrid{D,T},
    I::CartesianIndex{D},
) where {D,T}
    s = zero(Complex{T})
    idx = Tuple(I)
    @inbounds for d = 1:D
        s += Complex{T}(0, g.kvec[d][idx[d]]) * buffers[d][I]
    end
    return s
end

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
        k2 = _spectral_k2(g, I)
        if k2 != 0
            delta = rhohat[I] - _empic_divergence_hat(buffers, g, I)
            idx = Tuple(I)
            for d = 1:D
                buffers[d][I] += (-im * g.kvec[d][idx[d]] / k2) * delta
            end
        end
    end
    for c = 1:3
        g.iplan * buffers[c]
        es.E[c] .= real.(buffers[c])
    end
    return es.E
end

@inline function _add_species_charge!(
    rho::Array{T,D},
    ne::Array{T,D},
    es::EMPIC{D,T},
    s::ParticleSet{D,T},
) where {D,T}
    density!(ne, s, es.g, es.shape)
    q = s.q
    @inbounds for I in eachindex(rho)
        rho[I] += q * ne[I]
    end
    return rho
end

function _charge_density!(
    rho::Array{T,D},
    ne::Array{T,D},
    es::EMPIC{D,T},
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}},
) where {D,T}
    if es.mobile && ions !== nothing
        fill!(rho, zero(T))
        _add_species_charge!(rho, ne, es, ions)
        _add_species_charge!(rho, ne, es, e)
    else
        density!(ne, e, es.g, es.shape)
        @inbounds for I in eachindex(rho)
            rho[I] = es.n0 - ne[I]
        end
    end
    return rho
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
    ΔV = prod(es.g.dx)
    @inbounds for c = 1:3
        @. work = s.weight * s.v[c]
        deposit_scalar!(scratch, s, work, es.g, es.shape)
        f = s.q / ΔV
        if accumulate
            for I in eachindex(J[c])
                J[c][I] += f * scratch[I]
            end
        else
            for I in eachindex(J[c])
                J[c][I] = f * scratch[I]
            end
        end
    end
    return J
end

@inline function _swap_positions_with_midpoints!(
    s::ParticleSet{D,T},
    mid::NTuple{D,Vector{T}},
) where {D,T}
    @inbounds for d = 1:D
        x = s.x[d]
        xm = mid[d]
        length(xm) == length(x) ||
            throw(DimensionMismatch("midpoint buffer length must equal particle count"))
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
    _add_species_current!(es.Jraw, es, e, es.worke, es.ne, false)
    _swap_positions_with_midpoints!(e, es.mide)
    if es.mobile && ions !== nothing
        _swap_positions_with_midpoints!(ions, es.midi)
        _add_species_current!(es.Jraw, es, ions, es.worki, es.ne, true)
        _swap_positions_with_midpoints!(ions, es.midi)
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
        k2 = _spectral_k2(g, I)
        if k2 != 0
            target = -delta[I] / dt
            defect = target - _empic_divergence_hat(buffers, g, I)
            idx = Tuple(I)
            for d = 1:D
                buffers[d][I] += (-im * g.kvec[d][idx[d]] / k2) * defect
            end
        elseif I != dc
            # k2 == 0 covers the true DC mode AND every pure DC/Nyquist index
            # combination (kvec zeros the per-axis Nyquist entry). Only DC is
            # physical (net current, kept from the raw deposit). The Nyquist
            # combinations are structurally outside the spectral operators'
            # range: every derivative vanishes there, so raw deposition noise
            # would integrate into E as an undamped random walk (a growing
            # grid-Nyquist sawtooth felt by particles through the gather).
            # Zero ALL current components there, exactly as _esirkepov_Jx!
            # (longitudinal) and _deposit_Jy_at_midpoints! (transverse) do in
            # EMPIC1D. Mixed modes (e.g. Nyquist in x with k_y ≠ 0, so k2 ≠ 0)
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
    return es
end

function _substep_empic!(
    es::EMPIC{D,T},
    e::ParticleSet{D,T},
    ions::Union{Nothing,ParticleSet{D,T}},
    dt::T,
) where {D,T}
    g = es.g
    c2 = es.c * es.c
    _charge_density!(es.rho_n, es.ne, es, e, ions)

    curl!(es.curlE, es.E, g)
    @inbounds for c = 1:3, I in eachindex(es.B[c])
        es.Bc[c][I] = es.B[c][I] - (dt / 2) * es.curlE[c][I]
    end

    gather_vector!(es.Ep, es.E, e, g, es.shape)
    gather_vector!(es.Bp, es.Bc, e, g, es.shape)
    _push_velocities!(es, e, es.Ep, es.Bp, dt)
    _drift_record_mid!(e, es.mide, dt, g.L)
    apply_periodic!(e, ntuple(_ -> zero(T), D), g.L)

    if es.mobile && ions !== nothing
        gather_vector!(es.Epi, es.E, ions, g, es.shape)
        gather_vector!(es.Bpi, es.Bc, ions, g, es.shape)
        _push_velocities!(es, ions, es.Epi, es.Bpi, dt)
        _drift_record_mid!(ions, es.midi, dt, g.L)
        apply_periodic!(ions, ntuple(_ -> zero(T), D), g.L)
    end

    _deposit_current_at_midpoints!(es, e, ions)
    _charge_density!(es.rho_np1, es.ne, es, e, ions)
    _correct_current_continuity!(es, dt)

    @inbounds for c = 1:3, I in eachindex(es.B[c])
        es.B[c][I] -= dt * es.curlE[c][I]
    end
    curl!(es.curlB, es.B, g)
    @inbounds for c = 1:3, I in eachindex(es.E[c])
        es.E[c][I] += dt * (c2 * es.curlB[c][I] - es.J[c][I])
    end
    return es
end

"""
    step_empic!(es::EMPIC, electrons[, ions], dt)

Advance the dimension-parametric EM PIC state by one step. The optional ion
species is required when `mobile=true` and ignored otherwise.
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
    iszero(dtT) && return es            # dt=0 no-op: do not consume the one-time priming
    _ensure_empic_electron_buffers!(es, nparticles(e))
    if es.mobile
        _ensure_empic_ion_buffers!(es, nparticles(ions))
    end
    # prime the leapfrog once (2nd order): the seeded v and B are physical (v^0, B^0) but the
    # leapfrog carries them at the half level. Velocity: v^{-1/2}=v^0−(dt/2)a^0 for both species;
    # field: B^{-1/2}=B^0+(dt/2)∇×E^0 (∂tB=−∇×E). See the EMPIC1D header note.
    if es.step[] == 0
        g0 = es.g
        h = dtT / 2
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
                h,
                nparticles(ions),
                es.relativistic,
                es.c,
            )
        end
        curl!(es.curlE, es.E, g0)
        @inbounds for c = 1:3, I in eachindex(es.B[c])
            es.B[c][I] += h * es.curlE[c][I]
        end
    end
    _substep_empic!(es, e, ions, dtT)
    es.time[] += dtT
    es.step[] += 1
    return es
end

step_empic!(es::EMPIC{D,T}, e::ParticleSet{D,T}, dt::Real) where {D,T} =
    step_empic!(es, e, nothing, dt)

function em_field_energy(es::EMPIC{D,T}) where {D,T}
    c2 = es.c * es.c
    s = zero(T)
    @inbounds for c = 1:3, I in eachindex(es.E[c])
        s += es.E[c][I]^2 + c2 * es.B[c][I]^2
    end
    return T(0.5) * s * prod(es.g.dx)
end

function charge_conservation_residual(es::EMPIC{D,T}, dt::Real) where {D,T}
    g = es.g
    dtT = T(dt)
    dρhat = g.cbuf
    @inbounds for I in CartesianIndices(dρhat)
        dρhat[I] = Complex{T}(es.rho_np1[I] - es.rho_n[I])
    end
    g.plan * dρhat
    @inbounds for I in CartesianIndices(dρhat)
        _spectral_k2(g, I) == 0 && (dρhat[I] = zero(Complex{T}))
    end
    g.iplan * dρhat
    dρrepr = similar(es.rho_n)
    @inbounds for I in eachindex(dρrepr)
        dρrepr[I] = real(dρhat[I])
    end

    dJ = similar(es.rho_n)
    divergence!(dJ, es.J, g)
    rmax = zero(T)
    scale = zero(T)
    @inbounds for I in eachindex(dJ)
        dρ = dρrepr[I]
        r = dρ + dtT * dJ[I]
        rmax = max(rmax, abs(r))
        scale = max(scale, abs(dρ))
    end
    return rmax / (scale > 0 ? scale : one(T))
end
