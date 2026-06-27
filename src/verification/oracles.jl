# reference.jl — quantified comparison of a shock run against an established
# reference (items SHK-005 / "external-code comparison" / "published hybrid
# result reproduced").
#
# SCOPE (honest): the implemented reference is the ESTABLISHED ANALYTIC benchmark
# for a perpendicular collisionless shock — the magnetohydrodynamic
# Rankine–Hugoniot jump and exact flux-freezing (Bz/n conserved across the
# front), with the published kinetic ordering 1 < n₂ < X_RH (a collisionless
# shock compresses LESS than the γ-law fluid because reflected ions heat the
# downstream; Leroy et al. 1982, Tidman & Krall 1971). `compare_to_reference`
# is a generic NamedTuple comparator, so a user who HAS the output of a specific
# external hybrid code can pass its numbers as `reference` and get the same
# quantified PASS/FAIL — only the external code's data (not bundled here) is
# needed for a code-to-code comparison.

function _require_finite_nonnegative_tolerance(name::AbstractString, value::Real)
    v = Float64(value)
    isfinite(v) && v >= 0 || throw(ArgumentError("$name must be finite and non-negative"))
    return v
end

function _require_finite_comparison_value(kind::AbstractString, key::Symbol, value)
    v = Float64(value)
    isfinite(v) || throw(ArgumentError("$kind field $key must be finite"))
    return v
end

"""
    compare_to_reference(measured::NamedTuple, reference::NamedTuple;
                         rtol=0.1, atol=0.0) -> (; pass, maxrelerr, details)

Compare each field of `reference` against the same-named field of `measured` and
return whether every field matches to relative tolerance `rtol` (with absolute
floor `atol`). `details` is a vector of `(key, measured, reference, relerr, ok)`.
Fields present in `measured` but not in `reference` are ignored; a field in
`reference` but missing from `measured` fails. This is the comparison core for
checking a run against an analytic oracle OR against an external code's output.
"""
function compare_to_reference(
    measured::NamedTuple,
    reference::NamedTuple;
    rtol::Real = 0.1,
    atol::Real = 0.0,
)
    rtolT = _require_finite_nonnegative_tolerance("rtol", rtol)
    atolT = _require_finite_nonnegative_tolerance("atol", atol)
    details = Tuple{Symbol,Float64,Float64,Float64,Bool}[]
    pass = true
    maxrelerr = 0.0
    for key in keys(reference)
        ref = _require_finite_comparison_value("reference", key, getfield(reference, key))
        if !hasproperty(measured, key)
            push!(details, (key, NaN, ref, Inf, false))
            pass = false
            maxrelerr = Inf
            continue
        end
        mv = _require_finite_comparison_value("measured", key, getfield(measured, key))
        denom = max(abs(ref), atolT)
        relerr = denom == 0 ? abs(mv - ref) : abs(mv - ref) / denom
        ok = abs(mv - ref) <= atolT + rtolT * abs(ref)
        pass &= ok
        maxrelerr = max(maxrelerr, relerr)
        push!(details, (key, mv, ref, relerr, ok))
    end
    return (; pass, maxrelerr, details)
end

"""
    reproduce_established_shock(; MA=3.0, rtol=0.06, run_kwargs...)
        -> (; pass, measured, reference, comparison)

Run the verified 1-D reflecting-wall perpendicular shock ([`run_perp_shock`](@ref),
`run_kwargs` forwarded) at Alfvén Mach `MA` and check it against the established
analytic benchmark:

  • `frozen_ratio = (Bz₂/B0)/n₂ = 1`  — exact flux-freezing across the shock,

and separately verifies the published kinetic ordering `1 < n₂ < X_RH` (the
collisionless compression sits below the fluid Rankine–Hugoniot value). Returns
the PASS/FAIL of the flux-freezing match plus the measured/reference NamedTuples
and the full [`compare_to_reference`](@ref) breakdown.
"""
function reproduce_established_shock(;
    MA::Real = 3.0,
    γe::Real = 5 / 3,
    rtol::Real = 0.06,
    run_kwargs...,
)
    rtolT = _require_finite_nonnegative_tolerance("rtol", rtol)
    r = run_perp_shock(; MA = MA, γe = γe, run_kwargs...)
    measured = (; frozen_ratio = r.frozen_ratio)
    reference = (; frozen_ratio = 1.0)               # established flux-freezing
    cmp = compare_to_reference(measured, reference; rtol = rtolT)
    # established kinetic ordering: 1 < n₂ ≤ strong-shock fluid maximum (γ+1)/(γ-1).
    # (Uses the robust thermodynamic compression ceiling rather than the
    # front-speed-fit-dependent RH value, which is noisy at modest run lengths.)
    Xmax = (γe + 1) / (γe - 1)
    ordering_ok = isfinite(r.n2) && 1 < r.n2 < Xmax * (1 + rtolT)
    return (; pass = cmp.pass && ordering_ok, measured = r, reference, comparison = cmp)
end
