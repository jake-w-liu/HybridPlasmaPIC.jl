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
the remaining trajectory. Conversion analysis runs for `:cld2x2`
(MATLAB-pinned) and `:cld3x3` (the corrected extension of upstream's
sign-broken 3×3 layer; see the port notes).

Returns `(; rays, conversions)`:
  * `rays` — vector of `(; parent, sigma, y, status)` per traced ray
    (`parent = 0` for the antenna ray, otherwise the index of the ray whose
    conversion spawned it);
  * `conversions` — vector of `(; ray, sigma, conversion::RayconConversion)`.

With `amplitude = true` the result also carries `amplitude` — one entry per
ray with the WKB transport along it (`sigma`, `y`, focusing tensor `W`,
`lnE2`, `phase`, per-species deposition `dep`, `inkspace`, `nmaslov`,
`status`): the antenna ray starts from `W0` (default: the flux-surface-aligned
[`antenna_focusing`](@ref)), `lnE20`, `phase0`; at each conversion the
TRANSMITTED child inherits the focusing tensor with `lnE² − 2πη²` and a fresh
deposition ledger, while the CONVERTED continuation gains `ln|β|²`, `arg β`,
and a re-matched focusing tensor (upstream 'Trs'/'Cnv' Amp blocks).
`damping = false` disables collisional-less absorption/deposition.
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
    amplitude::Bool = false,
    W0::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    lnE20::Real = 0.0,
    phase0::Real = 0.0,
    damping::Bool = true,
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
            # conversion detection while the per-ray cap has headroom (with
            # max_conversions = 0 the ray traces the full span); supported for
            # :cld2x2 (MATLAB-pinned) and :cld3x3 (corrected extension)
            det = prob.model in (:cld2x2, :cld3x3) && nconv < max_conversions
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
    if !amplitude
        return (; rays, conversions)
    end
    amp = _amplitude_pass(prob, rays, conversions, y0; W0, lnE20, phase0, damping, rtol, atol)
    return (; rays, conversions, amplitude = amp)
end

# Amplitude reconstruction along an already-traced ray tree: re-integrate the
# amplitude transport over each ray's σ range, applying the conversion
# bookkeeping at every split (transmitted child: lnE² − 2πη², deposition
# reset, focusing tensor carried over — upstream 'Trs'; converted
# continuation: lnE² + ln|β|², phase + arg β, focusing tensor re-matched via
# Slam — upstream 'Cnv'). The j-th recorded conversion spawned ray j+1, and
# parents always precede children, so a single forward pass suffices.
function _amplitude_pass(
    prob::RayconProblem,
    rays,
    conversions,
    y0::Vector{Float64};
    W0::Union{Nothing,AbstractMatrix{<:Real}},
    lnE20::Real,
    phase0::Real,
    damping::Bool,
    rtol::Real,
    atol::Real,
)
    W0v = W0 === nothing ? antenna_focusing(prob, y0) : Matrix{Float64}(W0)
    nrays = length(rays)
    ns = length(prob.amass)
    # spawn packets: (W, lnE2, phase) at each ray's birth; a child stays
    # `nothing` if its parent's amplitude transport failed before the split
    spawnW = Vector{Union{Nothing,Matrix{Float64}}}(nothing, nrays)
    spawnln = fill(Float64(lnE20), nrays)
    spawnph = fill(Float64(phase0), nrays)
    spawnW[1] = W0v
    out = Vector{
        @NamedTuple{
            sigma::Vector{Float64},
            y::Matrix{Float64},
            W::Matrix{Float64},
            lnE2::Vector{Float64},
            phase::Vector{Float64},
            dep::Matrix{Float64},
            inkspace::Vector{Bool},
            nmaslov::Int,
            status::Symbol,
        }
    }(
        undef,
        nrays,
    )
    for i = 1:nrays
        if spawnW[i] === nothing
            out[i] = (;
                sigma = Float64[],
                y = zeros(4, 0),
                W = zeros(3, 0),
                lnE2 = Float64[],
                phase = Float64[],
                dep = zeros(ns, 0),
                inkspace = Bool[],
                nmaslov = 0,
                status = :parent_amplitude_unavailable,
            )
            continue
        end
        r = rays[i]
        σ0 = r.sigma[1]
        σend = r.sigma[end]
        # this ray's own conversion events, in σ order (already sorted)
        events = [(j, conversions[j]) for j = 1:length(conversions) if conversions[j].ray == i]
        Wcur = spawnW[i]::Matrix{Float64}
        lncur = spawnln[i]
        phcur = spawnph[i]
        segs = AmplitudeTrace[]
        σa = σ0
        ya = r.y[:, 1]
        ok = true
        for (j, ev) in events
            tr = integrate_ray_amplitude(
                prob,
                ya,
                Wcur,
                σa,
                ev.sigma;
                lnE20 = lncur,
                phase0 = phcur,
                damping,
                rtol,
                atol,
            )
            push!(segs, tr)
            if tr.status !== :end_of_span
                ok = false
                break
            end
            uend = vcat(tr.y[:, end], tr.W[:, end], tr.lnE2[end], tr.phase[end], tr.dep[:, end])
            if tr.inkspace[end]
                # conversion bookkeeping lives in x-space: transform back
                # (upstream 'convert_list' does the same temporary k→x switch)
                _, mok = _maslov!(uend, true)
                mok || (ok = false; break)
            end
            Wcur = [uend[5] uend[6]; uend[6] uend[7]]
            lncur = uend[8]
            phcur = uend[9]
            cv = ev.conversion
            # transmitted child (ray j+1): same W, amplitude reduced by the
            # transmission factor, fresh deposition ledger
            spawnW[j+1] = copy(Wcur)
            spawnln[j+1] = lncur - 2π * cv.eta2
            spawnph[j+1] = phcur
            # converted continuation of THIS ray
            lncur += log(abs2(cv.beta))
            phcur += angle(cv.beta)
            Wcur = _slam_matching(Wcur, cv.gdalf, cv.gdlam)
            σa = ev.sigma
            ya = collect(cv.converted)
        end
        if ok && σend > σa
            tr = integrate_ray_amplitude(
                prob,
                ya,
                Wcur,
                σa,
                σend;
                lnE20 = lncur,
                phase0 = phcur,
                damping,
                rtol,
                atol,
            )
            push!(segs, tr)
        end
        # concatenate segments
        sigma = reduce(vcat, [s.sigma for s in segs])
        y = reduce(hcat, [s.y for s in segs])
        W = reduce(hcat, [s.W for s in segs])
        lnE2 = reduce(vcat, [s.lnE2 for s in segs])
        phase = reduce(vcat, [s.phase for s in segs])
        # deposition ledgers accumulate across segments (each segment starts
        # at zero; offset by the running total at its start)
        dep = similar(segs[1].dep, length(prob.amass), 0)
        run = zeros(length(prob.amass))
        for s in segs
            dep = hcat(dep, s.dep .+ run)
            run = dep[:, end]
        end
        inkspace = reduce(vcat, [s.inkspace for s in segs])
        out[i] = (;
            sigma,
            y,
            W,
            lnE2,
            phase,
            dep,
            inkspace,
            nmaslov = sum(s.nmaslov for s in segs),
            status = segs[end].status,
        )
    end
    return out
end
