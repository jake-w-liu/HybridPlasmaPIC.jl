# raycon_driver.jl — programmatic equivalent of RAYCON's main.m driver loop:
# antenna launch, propagation, conversion detection, ray splitting, bookkeeping.

"""
    launch_ray(prob; s, theta, kr, kz, m=0.0) -> Vector{Float64}

Antenna launch (main.m): map the flux-coordinate antenna position `(s, θ)` to
`(R, Z)`, then adjust the wavevector onto the dispersion surface holding the
poloidal mode contribution `kθ = m/ρ` fixed ([`adjust_to_dispersion`](@ref)).
`kr`, `kz` are the antenna guesses (`kant[1]`, `kant[3]`).
"""
function launch_ray(prob::RayconProblem; s::Real, theta::Real, kr::Real, kz::Real, m::Real = 0.0)
    pos = map_flux(prob.eq, s, theta)
    return adjust_to_dispersion(prob, (pos.r, pos.z, kr, kz); m)
end

"""
    trace_rays(prob; s, theta, kr, kz, m=0.0, sigma_span=5e-2,
               max_conversions=3, rtol=1e-6, atol=1e-7) -> NamedTuple

Trace a ray with mode-conversion splitting (the main.m 'Con' loop): launch
from the antenna, integrate until a conversion event or the end of the ray-
parameter span; at each valid conversion record `(τ, β, saddle)`, queue the
TRANSMITTED ray as a new initial condition and continue the incident ray as
the CONVERTED one (upstream convention). After `max_conversions` splits
(upstream limit 3; `0` disables splitting) the ray keeps tracing without
further conversion detection — upstream stops the ray there instead, losing
the remaining trajectory. Conversion analysis requires `:cld2x2`; with
`:cld3x3` rays are traced without splitting.

Returns `(; rays, conversions)`:
  * `rays` — vector of `(; parent, sigma, y, status)` per traced ray
    (`parent = 0` for the antenna ray, otherwise the index of the ray whose
    conversion spawned it);
  * `conversions` — vector of `(; ray, sigma, conversion::RayconConversion)`.
"""
function trace_rays(
    prob::RayconProblem;
    s::Real,
    theta::Real,
    kr::Real,
    kz::Real,
    m::Real = 0.0,
    sigma_span::Real = 5e-2,
    max_conversions::Integer = 3,
    rtol::Real = 1e-6,
    atol::Real = 1e-7,
)
    (isfinite(sigma_span) && sigma_span > 0) || throw(ArgumentError("sigma_span must be positive"))
    max_conversions >= 0 || throw(ArgumentError("max_conversions must be ≥ 0"))
    y0 = launch_ray(prob; s, theta, kr, kz, m)
    span = Float64(sigma_span)

    queue = [(parent = 0, y = y0, sigma0 = 0.0)]
    rays = @NamedTuple{parent::Int, sigma::Vector{Float64}, y::Matrix{Float64}, status::Symbol}[]
    conversions = @NamedTuple{ray::Int, sigma::Float64, conversion::RayconConversion}[]

    while !isempty(queue)
        item = popfirst!(queue)
        segs_sigma = Float64[]
        segs_y = Matrix{Float64}(undef, 4, 0)
        y = collect(item.y)
        σ = item.sigma0
        nconv = 0
        nseg = 0
        status = :end_of_span
        while σ < span
            (nseg += 1) <= 50 || (status = :segment_limit; break)
            # conversion detection only where the analysis is available
            # (:cld2x2) and while the per-ray cap has headroom — with
            # max_conversions = 0 or :cld3x3 the ray traces the full span
            det = prob.model === :cld2x2 && nconv < max_conversions
            tr = integrate_ray(prob, y, σ, span; rtol, atol, detect_conversion = det)
            segs_sigma = vcat(segs_sigma, tr.sigma)
            segs_y = hcat(segs_y, tr.y)
            status = tr.status
            if tr.status !== :conversion_event
                break
            end
            σ = tr.sigma[end]
            y = tr.y[:, end]
            # a failed analysis (saddle outside the plasma, degenerate normal
            # form, non-converging polarization) skips the split and lets the
            # ray continue — upstream's abort-conversion path
            conv = try
                analyze_conversion(prob, y, tr.zdot, tr.zddot)
            catch e
                (e isa DomainError || e isa ArgumentError || e isa ErrorException) || rethrow()
                nothing
            end
            rayidx = length(rays) + 1        # index this ray will get below
            if conv !== nothing && is_valid(conv)
                nconv += 1
                push!(conversions, (; ray = rayidx, sigma = σ, conversion = conv))
                push!(queue, (parent = rayidx, y = collect(conv.transmitted), sigma0 = σ))
                # incident ray continues as the converted ray from the stop point
            end
            # resume from the stop point with a fresh monitor history
            # (upstream propagate restart)
            σ >= span && break
        end
        push!(rays, (; parent = item.parent, sigma = segs_sigma, y = segs_y, status))
    end
    return (; rays, conversions)
end
