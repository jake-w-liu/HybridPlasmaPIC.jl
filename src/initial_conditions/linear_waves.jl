# kdv.jl — KdV soliton verification (KDV-001). Not a hybrid-PIC test: it
# exercises nonlinear spectral differentiation, 2/3 dealiasing, and dispersive
# time integration (the same machinery the hybrid Ohm's-law products need).
#
#   u_t + c0 u_x + α u u_x + β u_xxx = 0
#
# Integrating-factor RK4: the stiff dispersive linear part L(k)=i(βk³−c0k) is
# integrated exactly via the factor e^{Ldt}; classical RK4 advances the
# nonlinear term, whose quadratic product is 2/3-dealiased.

"""
    kdv_soliton(A, c0, α, β, x, t, x0, Ld)

Analytic KdV soliton `A·sech²((x−Vt−x0)/L)` on a periodic domain of length `Ld`,
with `V = c0 + αA/3` and `L = √(12β/(αA))`.
"""
function kdv_soliton(A, c0, α, β, x, t, x0, Ld)
    Ls = sqrt(12β / (α * A))
    V = c0 + α * A / 3
    xc = mod(x0 + V * t, Ld)
    return A .* sech.((mod.(x .- xc .+ Ld / 2, Ld) .- Ld / 2) ./ Ls) .^ 2
end

function _complex_copy(v::AbstractVector{T}) where {T<:AbstractFloat}
    out = Vector{Complex{T}}(undef, length(v))
    z = zero(T)
    j = 0
    @inbounds for i in eachindex(v)
        j += 1
        out[j] = Complex{T}(v[i], z)
    end
    return out
end

function _real_copy(v::AbstractVector{T}) where {T<:AbstractFloat}
    out = Vector{T}(undef, length(v))
    j = 0
    @inbounds for i in eachindex(v)
        j += 1
        out[j] = v[i]
    end
    return out
end

"""
    kdv_solve(u0, Ld, c0, α, β, dt, nsteps; dealias=true)

Advance the KdV equation `nsteps` steps on a periodic domain of length `Ld`.
`dealias` applies the 2/3 rule to the quadratic nonlinearity (an explicit,
documented spectral filter; disable to study its effect).
"""
function kdv_solve(
    u0::AbstractVector{T},
    Ld::Real,
    c0::Real,
    α::Real,
    β::Real,
    dt::Real,
    nsteps::Integer;
    dealias::Bool = true,
) where {T<:AbstractFloat}
    n = length(u0)
    n > 0 || throw(ArgumentError("u0 must be non-empty"))
    nsteps >= 0 || throw(ArgumentError("nsteps must be non-negative"))
    nsteps <= typemax(Int) || throw(ArgumentError("nsteps must fit in Int"))
    LdT = T(Ld)
    c0T = T(c0)
    αT = T(α)
    βT = T(β)
    dtT = T(dt)
    isfinite(LdT) && LdT > zero(T) || throw(ArgumentError("Ld must be finite and positive"))
    isfinite(c0T) || throw(ArgumentError("c0 must be finite"))
    isfinite(αT) || throw(ArgumentError("α must be finite"))
    isfinite(βT) || throw(ArgumentError("β must be finite"))
    isfinite(dtT) || throw(ArgumentError("dt must be finite"))
    nsteps_i = Int(nsteps)
    nsteps_i == 0 && return _real_copy(u0)
    k = collect(T(2π) .* FFTW.fftfreq(n, T(n) / LdT))   # angular wavenumber
    if iseven(n)
        k[n÷2+1] = zero(T)                            # zero Nyquist (odd derivatives)
    end
    Lop = im .* (βT .* k .^ 3 .- c0T .* k)
    E = exp.(dtT .* Lop)
    E2 = exp.(dtT / 2 .* Lop)
    ik = im .* k
    kc = T(2 / 3) * (T(π) * n / LdT)                   # 2/3 dealias cutoff
    mask = abs.(k) .<= kc
    u0c = _complex_copy(u0)
    P = plan_fft!(u0c)
    Pinv = plan_ifft!(u0c)
    uh = copy(u0c)
    P * uh
    u = similar(u0c)
    u2 = similar(u0c)
    tmp = similar(u0c)
    k1 = similar(u0c)
    k2 = similar(u0c)
    k3 = similar(u0c)
    k4 = similar(u0c)
    coeff = -T(α) / 2
    function Nhat!(out, uhat)
        copyto!(u, uhat)
        Pinv * u
        @inbounds for i = 1:n
            ui = real(u[i])
            u2[i] = Complex{T}(ui * ui, zero(T))
        end
        P * u2
        @inbounds for i = 1:n
            val = dealias && !mask[i] ? zero(Complex{T}) : u2[i]
            out[i] = coeff * ik[i] * val                  # −α (u²/2)_x
        end
        return out
    end
    half = T(1) / 2
    sixth = T(1) / 6
    for _ = 1:nsteps_i
        Nhat!(k1, uh)
        @inbounds for i = 1:n
            k1[i] *= dtT
            tmp[i] = E2[i] * (uh[i] + half * k1[i])
        end
        Nhat!(k2, tmp)
        @inbounds for i = 1:n
            k2[i] *= dtT
            tmp[i] = E2[i] * uh[i] + half * k2[i]
        end
        Nhat!(k3, tmp)
        @inbounds for i = 1:n
            k3[i] *= dtT
            tmp[i] = E[i] * uh[i] + E2[i] * k3[i]
        end
        Nhat!(k4, tmp)
        @inbounds for i = 1:n
            k4[i] *= dtT
            uh[i] = E[i] * uh[i] + sixth * (E[i] * k1[i] + 2 * E2[i] * (k2[i] + k3[i]) + k4[i])
        end
    end
    copyto!(u, uh)
    Pinv * u
    out = Vector{T}(undef, n)
    @inbounds for i = 1:n
        out[i] = real(u[i])
    end
    return out
end

function kdv_solve(
    u0::AbstractVector{T},
    Ld::Real,
    c0::Real,
    α::Real,
    β::Real,
    dt::Real,
    nsteps;
    dealias::Bool = true,
) where {T<:AbstractFloat}
    throw(ArgumentError("nsteps must be an integer"))
end
