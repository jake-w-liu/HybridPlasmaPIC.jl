# empic1d.jl — full electromagnetic 1D PIC (1D space, 3V velocity), SEPARATE from
# the hybrid loop. Mobile kinetic electrons (q=−1, m=1). The positive species can
# be EITHER a uniform immobile ion background (+n0, the default) OR a second
# mobile kinetic ion species (q=+1, mass mi, default mi=1836·me). ε0=1, so the
# electron plasma frequency is ω_pe = √(n0) and ω_pe = 1 at n0 = 1. Speed of light
# c is a free parameter (default 5).
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
# (electrons resolve ω_pe and the EM CFL); the heavy ions are pushed ONCE over
# the full dt using the field gathered at the substep that straddles their
# half-level midpoint, and they contribute their (slowly varying) charge/current
# to the field solve on every electron substep. n_sub=1 reproduces the single-
# rate scheme exactly (and bit-for-bit reproduces the legacy immobile-ion run).
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
    mi > 0 || throw(ArgumentError("ion mass mi must be > 0"))
    n = g.n[1]
    z() = zeros(T, n)
    vz() = zeros(T, Nparticles)
    EMPIC1D{T,typeof(g),typeof(shape)}(
        g,
        T(n0),
        T(c),
        shape,
        relativistic,
        mobile,
        T(mi),
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
    if es.mobile
        ions === nothing &&
            throw(ArgumentError("mobile=true requires an ion ParticleSet in init_empic!"))
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

    # push ions on the straddling substep (heavy: pushed once per full ion step)
    if es.mobile && ions !== nothing && push_ions
        gather_scalar!(es.Exi, es.Ex, ions, g, es.shape)
        gather_scalar!(es.Eyi, es.Ey, ions, g, es.shape)
        gather_scalar!(es.Bzi, es.Bzc, ions, g, es.shape)
        # ions advance over the FULL ion step dt_ion = n_sub·dt
        dt_ion = T(es.n_sub) * dt
        _push_velocities!(es, ions, es.Exi, es.Eyi, es.Bzi, dt_ion)
        _drift_record_mid!(ions, es.worki, dt_ion, L)
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
    return es.Jy
end

# ---------------------------------------------------------------- step

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
    if es.mobile
        ions === nothing &&
            throw(ArgumentError("mobile=true requires an ion ParticleSet in step_empic!"))
        _ensure_ion_buffers!(es, nparticles(ions))
    end
    dtT = T(dt)
    ns = es.n_sub
    dt_e = dtT / ns
    # ions pushed on the substep straddling their half-level midpoint (middle one)
    push_sub = (ns + 1) ÷ 2
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
