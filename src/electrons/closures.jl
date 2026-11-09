# hybrid.jl — the massless-electron hybrid field model in Ω_ci-normalized units
# (§6–7 of the checklist). Proton-only quasineutrality n_e = n_i = n; the
# electric field is algebraic (generalized Ohm's law); scalar pressure can
# optionally evolve with B.
#
# Normalized equations:
#   J = ∇×B
#   u_e = u_i − J/n
#   E = −u_i×B + (J×B)/n − ∇p_e/n + η J          (hyperresistivity optional, off here)
#   ∂B/∂t = −∇×E
# Boris with q/m = 1 advances protons in these E, B.

# ---------------------------------------------------------------- electron closures

abstract type ElectronClosure end

"""
    EvolvingElectrons(γ = 5/3)

Adiabatic scalar electron pressure, evolved as a persistent `HybridFields.pe`
state by `HybridStepper` or `ImplicitHybridStepper`. Set the initial pressure
before `init!`; zero is a valid cold state. Requires finite `γ > 1` and an ideal massless `HybridModel`
(`η = ηH = de2 = 0`). CAM-CL and Hall-MHD reject this closure.

The equation is `∂pe/∂t = −ue⋅∇pe − γ*pe*∇⋅ue`, with
`ue=ui−J/max(n,nfloor)`. Pressure and B share each RK4 stage and are stored
at integer time; explicit particle velocity retains its half-time level.
`ImplicitHybridStepper` instead solves ion velocity, B and pressure jointly at
one synchronous integer-time boundary, with its own residual/work contract.
Pressure must be finite and nonnegative at every admitted stage or implicit
trial. Neither scheme is proved positivity preserving. Rejected numerical steps
roll back public state and never clip pressure. Checkpoints preserve
`pe`, and electron-energy diagnostics require it explicitly through `pe=`.
No heating, heat flux, or discrete total-energy conservation is implied;
[`TransportElectrons`](@ref) adds paired constant Ohmic and hyperresistive heating.
The magnetic [`recommended_dt`](@ref) does not bound pressure advection or
compression. The explicit spectral pressure scheme is not positivity preserving.
"""
struct EvolvingElectrons{T<:AbstractFloat} <: ElectronClosure
    γ::T
    function EvolvingElectrons{T}(γ::Real) where {T<:AbstractFloat}
        value = _require_finite_real("γ", γ, T)
        value > one(T) || throw(ArgumentError("evolving electron γ must be greater than one"))
        return new{T}(value)
    end
end

EvolvingElectrons(γ::Real = 5 / 3) = EvolvingElectrons{float(typeof(γ))}(γ)
closure_gamma(c::EvolvingElectrons) = c.γ
_evolves_pressure(::ElectronClosure) = false
_evolves_pressure(::EvolvingElectrons) = true

electron_pressure!(::Any, ::Any, ::EvolvingElectrons) = throw(
    ArgumentError(
        "EvolvingElectrons requires stored pressure; initialize fields.pe and use ohms_law!",
    ),
)

"""
    TransportElectrons(γ = 5/3)

Stored scalar electron pressure with constant Ohmic and hyperresistive transport for a serial,
periodic CPU `HybridStepper`. Use `HybridModel(TransportElectrons(γ); η, ηH)` with
finite constant `η ≥ 0`, `ηH ≥ 0` and `de2 = 0`; `γ > 1`. Generalized Ohm's law is
split into the ion force and the resistive field,

    E_force = −ui×B + (J×B)/n − ∇pe/n,      E = E_force + ηJ − ηH∇²J,
    ∂B/∂t = −∇×E,   ∂pe/∂t = −ue⋅∇pe − γ*pe*∇⋅ue + (γ−1)*(η*|J|² + H_H),

with `ue = ui − J/n` and positive hyperresistive heat `H_H`. Current-gradient
squares plus the full-Laplacian Nyquist channel pair the periodic magnetic
loss; see `docs/TRANSPORT_HYPERRESISTIVITY_PLAN.md`. Faraday's law uses the full `E`, while the first-step
priming, every Boris kick and the endpoint predictor use `E_force`. The
resistive friction is thus paired with its ion recoil, and its heat `η|J|²` enters
the single electron reservoir `∫pe/(γ−1) dV`. This is a prescribed reduced
friction model, not a species-resolved Coulomb collision or electron–ion
temperature-equilibration model.

Heating uses the current of the same RK4 stage as the magnetic field. Density
must be finite and at least `nfloor` at initialization and at every midpoint,
endpoint and predictor deposition; floor activation is rejected rather than
used as an untracked source. Scalar populations need unit charge; use a
one-element collection otherwise. Pressure validation, rollback and the absence
of clipping follow [`EvolvingElectrons`](@ref). The stepper owns the carried
`E_force`, which checkpoints and archives restore verbatim. With `η = ηH = 0` the
closure reproduces `EvolvingElectrons` bitwise. MPI, CAM-CL and Hall-MHD reject
this closure.
"""
struct TransportElectrons{T<:AbstractFloat} <: ElectronClosure
    γ::T
    function TransportElectrons{T}(γ::Real) where {T<:AbstractFloat}
        value = _require_finite_real("γ", γ, T)
        value > one(T) || throw(ArgumentError("transport electron γ must be greater than one"))
        return new{T}(value)
    end
end

TransportElectrons(γ::Real = 5 / 3) = TransportElectrons{float(typeof(γ))}(γ)
closure_gamma(c::TransportElectrons) = c.γ
_evolves_pressure(::TransportElectrons) = true
# Scalar transport closures split the ion force from the full electric field.
_transports_electrons(::ElectronClosure) = false
_transports_electrons(::TransportElectrons) = true

electron_pressure!(::Any, ::Any, ::TransportElectrons) = throw(
    ArgumentError(
        "TransportElectrons requires stored pressure; initialize fields.pe and use ohms_law!",
    ),
)

"""
    ConductingElectrons(γ=5/3; κ)

Stored scalar pressure with constant isotropic Fourier heat flux
`q_e = -κ*grad(Te)`, prepared from the exact quotient of the original stored
pressure and actual deposited density, `Te = Q(pe)/Q(n)`.
Use with a serial periodic CPU Float32/64 `HybridStepper` and
`HybridModel(...; η, ηH, de2=0)`. It shares the ion-force/field split and
Ohmic/hyperresistive heating of [`TransportElectrons`](@ref).

Requires finite `γ>1`, `κ>=0`; nonzero κ must remain finite and nonzero at grid
precision. Active conduction requires nearest hardware rounding and gradual
underflow. Each axis removes an exact constant-line background, then prepares
two normalized IEEE terms with residual bounded by `u²*abs(residual_input)`.
The full stored-wavenumber Laplacian retains Nyquist and projects each axis's
computed output to exact zero line mean before complete signed rate accumulation.
FFTW accuracy remains qualified; arbitrary quotient differentiation is not exact.
The existing grid-precision `γ-1` coefficient is retained. Only the complete
computed rate is rounded and range-checked, with full public-step rollback.
Pressure must remain nonnegative at every RK stage; spectral diffusion is not
positivity preserving. Stored rates and stages must fit grid precision, even
when a duration-scaled update could fit.

Set stored pressure before `init!`. `κ=0` retains the transport arithmetic and
performs no quotient preparation, transform or projection. Checkpoints retain
pressure, carried ion force, κ and required conduction policy 2, also for κ=0;
policy-1 payloads are rejected. The fixed-step
[`CollisionalHybridStepper`](@ref) composes it with ion-ion scattering, which
leaves pressure fixed. No extra reservoir, wall flux, electron-ion exchange,
MPI/device, CAM-CL or Hall-MHD support is implied.
"""
struct ConductingElectrons{T<:AbstractFloat} <: ElectronClosure
    γ::T
    κ::T
    function ConductingElectrons{T}(γ::Real; κ::Real) where {T<:AbstractFloat}
        gamma = _require_finite_real("γ", γ, T)
        gamma > one(T) || throw(ArgumentError("conducting electron γ must be greater than one"))
        conductivity = _require_finite_real("κ", κ, T)
        _validate_conduction_mode(κ)
        conductivity >= zero(T) || throw(ArgumentError("κ must be nonnegative"))
        (_conduction_zero(κ) || !_conduction_zero(conductivity)) ||
            throw(ArgumentError("ConductingElectrons κ underflows to zero at closure precision"))
        return new{T}(gamma, conductivity)
    end
end

function ConductingElectrons(γ::Real = 5 / 3; κ::Real)
    T = float(promote_type(typeof(γ), typeof(κ)))
    return ConductingElectrons{T}(γ; κ)
end

"""
    FieldAlignedConductingElectrons(γ=5/3; κ_parallel, κ_perpendicular=0)

Stored scalar electron pressure with field-aligned Fourier conduction on a
serial periodic CPU `HybridStepper`. The conductivity tensor is
`κ_perpendicular*I + (κ_parallel-κ_perpendicular)*B*B'/dot(B,B)` using the actual
magnetic field at each joint B/pressure RK stage and all three field components.
Unequal grid-resolved conductivities require nonzero B at every grid point.
Equal conductivities use the isotropic operator of `ConductingElectrons`.

Both coefficients must be finite, nonnegative and survive grid conversion.
Active conduction uses original stored pressure/density quotients, exact
projectors and two-term preparations of gradient and flux inputs. The real
odd-gradient tensor has separate diagonal Nyquist complements. FFTW error and
the coupled explicit RK map have no unconditional positivity or stability
guarantee. Only the complete pressure source is rounded and range-checked.
Policy-3 checkpoints retain the operator convention and carried ion force.
Ion-ion collisional composition leaves pressure fixed during scattering.
Electron-ion exchange, CGL pressure, implicit and distributed/device evolution
are outside this closure's contract.
"""
struct FieldAlignedConductingElectrons{T<:AbstractFloat} <: ElectronClosure
    γ::T
    κ_parallel::T
    κ_perpendicular::T
    function FieldAlignedConductingElectrons{T}(
        γ::Real;
        κ_parallel::Real,
        κ_perpendicular::Real = 0,
    ) where {T<:AbstractFloat}
        gamma = _require_finite_real("γ", γ, T)
        gamma > one(T) || throw(ArgumentError("conducting electron γ must be greater than one"))
        parallel = _field_aligned_coefficient("κ_parallel", κ_parallel, T)
        perpendicular = _field_aligned_coefficient("κ_perpendicular", κ_perpendicular, T)
        return new{T}(gamma, parallel, perpendicular)
    end
end

function FieldAlignedConductingElectrons(
    γ::Real = 5 / 3;
    κ_parallel::Real,
    κ_perpendicular::Real = 0,
)
    T = float(promote_type(typeof(γ), typeof(κ_parallel), typeof(κ_perpendicular)))
    return FieldAlignedConductingElectrons{T}(γ; κ_parallel, κ_perpendicular)
end

function _field_aligned_coefficient(name, value, ::Type{T}) where {T}
    converted = _require_finite_real(name, value, T)
    _validate_conduction_mode(value)
    value >= zero(value) && converted >= zero(T) ||
        throw(ArgumentError("$name must be nonnegative"))
    (_conduction_zero(value) || !_conduction_zero(converted)) ||
        throw(ArgumentError("$name underflows to zero at conduction precision"))
    return converted
end

const _ScalarTransportClosure =
    Union{TransportElectrons,ConductingElectrons,FieldAlignedConductingElectrons}
closure_gamma(c::FieldAlignedConductingElectrons) = c.γ
_evolves_pressure(::FieldAlignedConductingElectrons) = true
_transports_electrons(::FieldAlignedConductingElectrons) = true
electron_pressure!(::Any, ::Any, ::FieldAlignedConductingElectrons) =
    throw(ArgumentError("FieldAlignedConductingElectrons requires stored pressure; use ohms_law!"))
closure_gamma(c::ConductingElectrons) = c.γ
_evolves_pressure(::ConductingElectrons) = true
_transports_electrons(::ConductingElectrons) = true

electron_pressure!(::Any, ::Any, ::ConductingElectrons) = throw(
    ArgumentError(
        "ConductingElectrons requires stored pressure; initialize fields.pe and use ohms_law!",
    ),
)

_conduction_zero(value::Union{Float16,Float32,Float64,BigFloat}) = iszero(_exact_rational(value))
_conduction_zero(value::Real) = iszero(value)
function _validate_conduction_mode(κ)
    (_conduction_zero(κ) || _boris_accuracy_mode()) || throw(
        ArgumentError("active electron conduction requires nearest rounding and gradual underflow"),
    )
    return nothing
end

_electron_conductivity(::ElectronClosure, ::Type{T}) where {T} = zero(T)
function _electron_conductivity(closure::ConductingElectrons, ::Type{T}) where {T}
    conductivity = _require_finite_real("κ", closure.κ, T)
    _validate_conduction_mode(closure.κ)
    conductivity >= zero(T) || throw(ArgumentError("κ must be nonnegative"))
    (_conduction_zero(closure.κ) || !_conduction_zero(conductivity)) ||
        throw(ArgumentError("ConductingElectrons κ underflows to zero at grid precision"))
    return conductivity
end

function _field_aligned_conductivities(c::FieldAlignedConductingElectrons, ::Type{T}) where {T}
    # Validate both original operands, including the smaller coefficient.
    return (
        _field_aligned_coefficient("κ_parallel", c.κ_parallel, T),
        _field_aligned_coefficient("κ_perpendicular", c.κ_perpendicular, T),
    )
end

function _electron_conductivity(c::FieldAlignedConductingElectrons, ::Type{T}) where {T}
    parallel, perpendicular = _field_aligned_conductivities(c, T)
    return max(parallel, perpendicular)
end

"Isothermal electrons: p_e = T_e · n  (T_e the normalized electron temperature)."
struct IsothermalElectrons{T} <: ElectronClosure
    Te::T
end

"""
    PolytropicElectrons(pe0, n0, γ)

Scalar pressure `pe = pe0 * (max(n, 0) / n0)^γ`, with finite `pe0 >= 0`,
`n0 > 0` and `γ > 0`. The density floor supports transient Hall-MHD predictors;
the committed state is separately validated against its density floor.

For Float32/64 data, coefficients and destinations, pressure evaluation in
`Array` storage retains the stored input factors. Supported wrappers are views,
reshapes, reinterpretations, permutations, transposes and adjoints backed by
`Array`. Consecutive equal finite densities in compatible arrays reuse their
pressure within that call; no cache survives the call.
Ordinary powers use a division-residual correction; sensitive range or large-
correction cases use scaled intervals and exact midpoint checks. Other element
types and storage backends retain their native broadcast arithmetic.
This pointwise pressure contract does not imply exact integrated electron energy.
"""
struct PolytropicElectrons{T} <: ElectronClosure
    pe0::T
    n0::T
    γ::T
end

function IsothermalElectrons(Te::Real)
    T = float(typeof(Te))
    return IsothermalElectrons{T}(_require_finite_nonnegative_real("Te", Te, T))
end

function PolytropicElectrons(pe0::Real, n0::Real, γ::Real)
    T = float(promote_type(typeof(pe0), typeof(n0), typeof(γ)))
    pe0T = _require_finite_nonnegative_real("pe0", pe0, T)
    n0T = _require_finite_positive_real("n0", n0, T)
    γT = _require_finite_positive_real("γ", γ, T)
    return PolytropicElectrons{T}(pe0T, n0T, γT)
end

@inline _isothermal_pressure_value(Te, n, ops::_ScalarOps = _DefaultScalarOps()) =
    _scalar_mul(ops, Te, n)
electron_pressure!(pe, n, c::IsothermalElectrons) =
    (@. pe = _isothermal_pressure_value(c.Te, n); pe)
# Floor the base at 0: particle-deposited densities are ≥0 (so this is a no-op for hybrid/full-PIC),
# but a Hall-MHD RK4 predictor stage can transiently undershoot a cell below 0, and (neg)^γ with a
# non-integer γ throws an uncontrolled DomainError. pe→0 there (continuous, since 0^γ=0 for γ>0);
# the committed state is separately validated against nfloor.
_legacy_polytropic_pressure!(pe, n, c::PolytropicElectrons) =
    (@. pe = c.pe0 * (max(n, zero(eltype(n))) / c.n0)^c.γ; pe)
electron_pressure!(pe, n, c::PolytropicElectrons) = _legacy_polytropic_pressure!(pe, n, c)

"Adiabatic index of a closure (for the electron internal-energy budget)."
closure_gamma(::IsothermalElectrons) = 1.0
closure_gamma(c::PolytropicElectrons) = float(c.γ)

"""
    CGLElectrons(p_perp0, p_par0, n0, B0)

Anisotropic **double-adiabatic (Chew–Goldberger–Low)** electron closure. The gyrotropic
electron pressures follow the two CGL adiabatic invariants `p_⊥/(nB)=const` (frozen `μ`)
and `p_∥ B²/n³=const` (frozen `J_∥`), so from reference values `p_⊥0, p_∥0` at `n0, B0`:

    p_⊥(n,B) = C_⊥ · n · B,      C_⊥ = p_⊥0 / (n0 B0)
    p_∥(n,B) = C_∥ · n³ / B²,    C_∥ = p_∥0 B0² / n0³

The electron pressure is then the gyrotropic tensor `P_e = p_⊥ I + (p_∥ − p_⊥) bb`,
`b = B/|B|`, so the generalized Ohm's law carries `−(∇·P_e)/n` in place of the isotropic
`−∇p_e/n`. Equal reference pressures specify an isotropic reference in the ideal
parameterization; they do not imply an isotropic law away from it. A tensor field
with `p_∥ = p_⊥` throughout reduces to scalar pressure. For positive stored
coefficients and `n > 0`, `|B| > 0`, constitutive isotropy requires
`B³/n² = C_∥/C_⊥`. In the exact equal-reference parameterization this ratio is
`B0³/n0²`; the numerical model retains the constructor's rounded coefficients.
Unlike the scalar closures it depends on `B` (both magnitude and
direction), so the force `∇·P_e` is recomputed each B-subcycle stage rather than
frozen once per step. A step-frozen direction generally changes the stress and
does not implement the same constitutive law.

**Stability.** Pressure anisotropy can drive firehose or mirror instabilities.
Stability depends on the full model, equilibrium and resolved scales;
`p_∥ > p_⊥` alone does not establish stability. CGL and kinetic mirror thresholds
differ; see [Snyder, Hammett & Dorland (1997)](https://w3.pppl.gov/~hammett/gyrofluid/papers/snyder_pop97.pdf).
This electron closure supplies no FLR/Landau stabilization or universal numerical
stability guarantee. Unstable growth does not imply that every calculation must
blow up. MHD thresholds are not admission boundaries for this kinetic-ion/Hall
implementation; a limiter or kinetic closure would require its own model contract.
"""
struct CGLElectrons{T} <: ElectronClosure
    Cperp::T   # p_⊥ / (n B)
    Cpar::T    # p_∥ B² / n³
end

@inline function _cgl_perp_coefficient(pressure::T, density::T, field::T) where {T}
    iszero(pressure) && return zero(T)
    pf, pe = frexp(pressure)
    nf, ne = frexp(density)
    bf, be = frexp(field)
    return ldexp(pf / (nf * bf), pe - ne - be)
end

@inline function _cgl_parallel_coefficient(pressure::T, density::T, field::T) where {T}
    iszero(pressure) && return zero(T)
    pf, pe = frexp(pressure)
    nf, ne = frexp(density)
    bf, be = frexp(field)
    return ldexp((pf * bf * bf) / (nf * nf * nf), pe + 2be - 3ne)
end

@inline function _cgl_perpendicular_pressure(coefficient::T, density::T, field::T) where {T}
    (iszero(coefficient) || iszero(density) || iszero(field)) && return zero(T)
    value = coefficient * density * field
    isfinite(value) && !iszero(value) && return value
    return _finite_product3(coefficient, density, field)
end

@inline function _cgl_parallel_pressure(coefficient::T, density::T, field::T) where {T}
    (iszero(coefficient) || iszero(density)) && return zero(T)
    numerator = coefficient * density * density * density
    denominator = field * field
    value = numerator / denominator
    isfinite(value) && !iszero(value) && return value
    cf, ce = frexp(coefficient)
    nf, ne = frexp(density)
    bf, be = frexp(field)
    return ldexp((cf * nf * nf * nf) / (bf * bf), ce + 3ne - 2be)
end

function CGLElectrons(p_perp0::Real, p_par0::Real, n0::Real, B0::Real)
    T = float(promote_type(typeof(p_perp0), typeof(p_par0), typeof(n0), typeof(B0)))
    pp0 = _require_finite_nonnegative_real("p_perp0", p_perp0, T)
    ppa0 = _require_finite_nonnegative_real("p_par0", p_par0, T)
    n0T = _require_finite_positive_real("n0", n0, T)
    B0T = _require_finite_positive_real("B0", B0, T)
    Cperp = _cgl_perp_coefficient(pp0, n0T, B0T)
    Cpar = _cgl_parallel_coefficient(ppa0, n0T, B0T)
    isfinite(Cperp) && isfinite(Cpar) ||
        throw(ArgumentError("CGL invariant coefficient exceeds numeric range"))
    (iszero(pp0) || !iszero(Cperp)) ||
        throw(ArgumentError("CGL perpendicular invariant coefficient is below numeric range"))
    (iszero(ppa0) || !iszero(Cpar)) ||
        throw(ArgumentError("CGL parallel invariant coefficient is below numeric range"))
    return CGLElectrons{T}(Cperp, Cpar)
end

"Perpendicular CGL electron pressure `p_⊥ = C_⊥ n |B|` at density `n`, field magnitude `Bmag`."
@inline cgl_pperp(c::CGLElectrons{T}, n, Bmag) where {T} =
    _cgl_perpendicular_pressure(c.Cperp, T(n), T(Bmag))
"Parallel CGL electron pressure `p_∥ = C_∥ n³ / |B|²`."
@inline cgl_ppar(c::CGLElectrons{T}, n, Bmag) where {T} =
    _cgl_parallel_pressure(c.Cpar, T(n), T(Bmag))

closure_gamma(::CGLElectrons) = 5 / 3      # effective γ; the anisotropic budget is separate

# CGL has no scalar pressure — the Ohm's-law path dispatches to anisotropic_pressure_force!
# instead of _ohm_prep!/electron_pressure!. This method makes accidental scalar use explicit.
electron_pressure!(::Any, ::Any, ::CGLElectrons) = throw(
    ArgumentError(
        "CGLElectrons is gyrotropic and has no scalar electron pressure; use anisotropic_pressure_force!",
    ),
)

"Whether a closure produces a gyrotropic (anisotropic) pressure needing `∇·P_e` in Ohm's law."
is_anisotropic(::ElectronClosure) = false
is_anisotropic(::CGLElectrons) = true

# ---------------------------------------------------------------- model + state
