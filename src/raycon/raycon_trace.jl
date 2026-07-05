# raycon_trace.jl — adaptive Dormand–Prince 4(5) ray integration with the
# RAYCON conversion-monitor event logic (port of trajectory.m + the ode45 call
# in ray.m 'propagate').
#
# The upstream code integrates dz/dσ = J∇U with MATLAB ode45 (RelTol 1e-6,
# AbsTol 1e-7 for 'Con', InitialStep 1e-7·span) and records the conversion
# monitor mon2 = |tr DD| along the ray; a terminal event fires where the
# quadratic fit of the last 6 monitor samples has zero time-derivative (the
# closest approach of the avoided crossing), after a 15-sample warmup. This
# port records monitors at accepted steps and applies the same fit/eventing;
# the event resolves at step resolution (the conversion analysis only needs a
# point near the monitor minimum, as upstream notes).

# Dormand–Prince 4(5) coefficients
const _DP_A = (
    (1 / 5,),
    (3 / 40, 9 / 40),
    (44 / 45, -56 / 15, 32 / 9),
    (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
    (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
    (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84),
)
const _DP_E = (71 / 57600, 0.0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40)

"""
    RayconTrace

One integrated ray segment: ray parameter values `sigma`, states `y` (4×N:
r, z, kr, kz), the recorded monitors, the stop `status`
(`:end_of_span`, `:conversion_event`, `:left_domain`, `:step_underflow`,
`:max_steps`), and — when a conversion event fired — the detection state with
the divided-difference velocity and acceleration needed by
[`analyze_conversion`](@ref).
"""
struct RayconTrace
    sigma::Vector{Float64}
    y::Matrix{Float64}
    mon2::Vector{Float64}
    status::Symbol
    zdot::Vector{Float64}      # empty unless status == :conversion_event
    zddot::Vector{Float64}
end

# quadratic least-squares fit of v(t) over the given samples; returns the
# coefficients (c0, c1, c2) of c0 + c1·τ + c2·τ² with τ = t/t[1] (upstream
# rescaling). Times must be positive.
function _quadfit(ts::AbstractVector{Float64}, vs::AbstractVector{Float64})
    sca = ts[1]
    τ = ts ./ sca
    A = hcat(ones(length(τ)), τ, τ .^ 2)
    cf = A \ vs
    return cf[1], cf[2], cf[3], sca
end

"""
    integrate_ray(prob, y0, sigma0, sigma_end; rtol=1e-6, atol=1e-7,
                  initial_step=1e-7*(sigma_end-sigma0), detect_conversion=true,
                  warmup=15, max_steps=200_000) -> RayconTrace

Integrate one ray with adaptive DP4(5) from `sigma0` to `sigma_end`, recording
the conversion monitor and stopping at a conversion event (see file header).
Leaving the plasma (`s ≥ 1` profile domain) terminates cleanly with
`:left_domain`.
"""
function integrate_ray(
    prob::RayconProblem,
    y0::AbstractVector{<:Real},
    sigma0::Real,
    sigma_end::Real;
    rtol::Real = 1e-6,
    atol::Real = 1e-7,
    initial_step::Real = 1e-7 * (Float64(sigma_end) - Float64(sigma0)),
    detect_conversion::Bool = true,
    warmup::Integer = 15,
    max_steps::Integer = 200_000,
)
    length(y0) == 4 || throw(ArgumentError("state must be (r, z, kr, kz)"))
    all(isfinite, y0) || throw(ArgumentError("state must be finite"))
    σ0 = Float64(sigma0)
    σe = Float64(sigma_end)
    σe > σ0 || throw(ArgumentError("sigma_end must exceed sigma0"))
    (isfinite(rtol) && rtol > 0 && isfinite(atol) && atol > 0) ||
        throw(ArgumentError("rtol and atol must be positive"))
    h = Float64(initial_step)
    (isfinite(h) && h > 0) || throw(ArgumentError("initial_step must be positive"))

    y = Float64.(collect(y0))
    σ = σ0
    sigmas = [σ]
    ys = reshape(copy(y), 4, 1)
    mons = Float64[]
    status = :end_of_span
    zdot = Float64[]
    zddot = Float64[]

    f = yy -> trajectory_rhs(prob, yy)
    local k1
    try
        k1 = f(y)
        push!(mons, conversion_monitors(prob, y).mon2)
    catch e
        e isa DomainError && return RayconTrace(sigmas, ys, mons, :left_domain, zdot, zddot)
        rethrow()
    end

    adjcnv = NaN                    # event normalization (upstream adjinit)
    lastev = NaN                    # previous event value
    nacc = 1
    ks = Vector{Vector{Float64}}(undef, 7)
    # step-underflow scale follows the σ-SPAN (a span of 1e-8 uses steps far
    # below eps(1.0); an absolute threshold would abort tiny spans immediately)
    hmin = 16 * eps(max(abs(σ0), abs(σe)))
    status = :max_steps             # overwritten on every regular exit path
    for step = 1:max_steps
        h = min(h, σe - σ)
        if h <= hmin
            status = σe - σ <= hmin ? :end_of_span : :step_underflow
            break
        end
        # DP stages
        ks[1] = k1
        ynew = y                 # overwritten by stage 6; definite assignment
        try
            for i = 1:6
                acc = copy(y)
                a = _DP_A[i]
                for j = 1:i
                    acc .+= (h * a[j]) .* ks[j]
                end
                ks[i+1] = f(acc)
                i == 6 && (ynew = acc)   # stage 7 point IS the 5th-order solution (FSAL)
            end
        catch e
            if e isa DomainError
                status = :left_domain
                break
            end
            rethrow()
        end
        # error estimate
        errv = zeros(4)
        for j = 1:7
            _DP_E[j] == 0 && continue
            errv .+= (h * _DP_E[j]) .* ks[j]
        end
        sc = atol .+ rtol .* max.(abs.(y), abs.(ynew))
        err = sqrt(sum((errv ./ sc) .^ 2) / 4)
        if err > 1.0
            h *= max(0.2, 0.9 * err^(-0.2))
            continue
        end
        # accept
        σ += h
        y = ynew
        k1 = ks[7]
        h *= min(5.0, max(0.2, err > 0 ? 0.9 * err^(-0.2) : 5.0))
        nacc += 1
        push!(sigmas, σ)
        ys = hcat(ys, y)
        try
            push!(mons, conversion_monitors(prob, y).mon2)
        catch e
            if e isa DomainError
                status = :left_domain
                break
            end
            rethrow()
        end

        # ---- conversion event (trajectory.m 'events' + normalization) ----
        if detect_conversion && nacc > warmup && length(mons) >= 6 && sigmas[end-5] > 0
            ts = sigmas[end-5:end]
            vs = mons[end-5:end]
            c0, c1, c2, sca = _quadfit(ts, vs)
            τnow = σ / sca
            moncnv = c1 + 2 * c2 * τnow          # fitted d(mon2)/dτ at current σ
            if isnan(adjcnv)
                adjcnv = moncnv != 0 ? 1 / moncnv : 1.0
            end
            ev = adjcnv * moncnv
            if !isnan(lastev) && isfinite(ev) && sign(ev) != sign(lastev) && ev != 0
                # confirmation fit on |mon2| (ray.m 'convert_which'): a minimum
                # requires positive curvature and a derivative sign change
                n5 = min(5, length(mons) - 1) + 1
                t5 = sigmas[end-n5+1:end]
                v5 = abs.(mons[end-n5+1:end])
                d0, d1, d2, s5 = _quadfit(t5, v5)
                τs = t5 ./ s5
                z1a = d1 + 2 * d2 * τs[1]
                z1b = d1 + 2 * d2 * τs[end]
                if d2 > 0 && sign(z1a) != sign(z1b)
                    # detection point: velocity/acceleration by 3-point divided
                    # differences of the trajectory (ray.m 'convert_list')
                    t0 = sigmas[end]
                    tm1 = sigmas[end-1]
                    tm2 = sigmas[end-2]
                    zv0 = ys[:, end]
                    zm1 = ys[:, end-1]
                    zm2 = ys[:, end-2]
                    dm2 = (tm2 - tm1) * (tm2 - t0)
                    dm1 = (tm1 - tm2) * (tm1 - t0)
                    d0d = (t0 - tm2) * (t0 - tm1)
                    zdot =
                        zm2 ./ dm2 .* (t0 - tm1) .+ zm1 ./ dm1 .* (t0 - tm2) .+
                        zv0 ./ d0d .* (2 * t0 - tm1 - tm2)
                    zddot = 2 .* (zm2 ./ dm2 .+ zm1 ./ dm1 .+ zv0 ./ d0d)
                    status = :conversion_event
                    break
                end
            end
            lastev = ev
        end
        if σ >= σe
            status = :end_of_span
            break
        end
    end
    return RayconTrace(sigmas, ys, mons, status, zdot, zddot)
end
