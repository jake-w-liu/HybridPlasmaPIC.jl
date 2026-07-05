# reference.jl — quantified comparison of a shock run against an established
# reference (items SHK-005 / "external-code comparison" / "published hybrid
# result reproduced").
#
# SCOPE (honest): the implemented references are:
#   * the established analytic perpendicular-shock benchmark: MHD Rankine-
#     Hugoniot jump, exact flux-freezing, and the kinetic ordering 1 < n₂ < X_RH;
#   * a compact summary of one published external hybrid-code dataset from
#     Preisser et al. 2020, Zenodo DOI 10.5281/zenodo.3697360.
#
# The published-data summary is deliberately small and offline: it stores source
# metadata, the upstream file checksum, and scalar observables derived from the
# HDF5 figure-data file. It does not bundle the full HDF5 profile.

const _PREISSER2020_SHOCK_REFERENCE_ID = :preisser2020_65deg_Bavg_y

const _PREISSER2020_SHOCK_REFERENCE_META = (
    id = _PREISSER2020_SHOCK_REFERENCE_ID,
    title = "Influence of He++ and shock geometry on interplanetary shocks in the solar wind: 2D Hybrid simulations",
    creators = (
        "Luis Preisser",
        "Xochitl Blanco-Cano",
        "David Burgess",
        "Domenico Trotta",
        "Primoz Kajdic",
    ),
    doi = "10.5281/zenodo.3697360",
    doi_url = "https://doi.org/10.5281/zenodo.3697360",
    publication_date = "2020-03-05",
    license = "CC-BY-4.0",
    resource_type = "dataset",
    notes = "HDF5 data corresponding to publication figure panels; 2D local hybrid simulations with particle ions and massless-fluid electrons.",
    file = "Fig2_65deg_1perc_5perc_10perc_Bavg_y.h5",
    file_url = "https://zenodo.org/api/records/3697360/files/Fig2_65deg_1perc_5perc_10perc_Bavg_y.h5/content",
    file_checksum = "md5:2ea4f239d7221cd52607705f383161a1",
    derived_with = "HDF5.jl 0.17.3 read of the Zenodo file; scalar min/max/mean over each 1000-sample Bavg_y dataset.",
)

const _PREISSER2020_BAVG_REFERENCES = (
    (
        dataset = "65deg_1perc_Bavg_y",
        thetaBn_deg = 65.0,
        alpha_fraction = 0.01,
        nsamples = 1000,
        Bavg_y_min = 0.9968850596311658,
        Bavg_y_max = 3.698209909020219,
        Bavg_y_mean = 1.924506582194626,
    ),
    (
        dataset = "65deg_5perc_Bavg_y",
        thetaBn_deg = 65.0,
        alpha_fraction = 0.05,
        nsamples = 1000,
        Bavg_y_min = 0.997214995499391,
        Bavg_y_max = 3.500705042674052,
        Bavg_y_mean = 1.9245920883357377,
    ),
    (
        dataset = "65deg_10perc_Bavg_y",
        thetaBn_deg = 65.0,
        alpha_fraction = 0.10,
        nsamples = 1000,
        Bavg_y_min = 0.9969985157435545,
        Bavg_y_max = 3.5037309774874927,
        Bavg_y_mean = 1.9378039653740708,
    ),
)

_unknown_published_reference_error(id) = ArgumentError(
    "unknown published hybrid reference $(id); expected $(_PREISSER2020_SHOCK_REFERENCE_ID)",
)

"""
    published_hybrid_reference_ids() -> Tuple{Symbol}

Return the published external hybrid-code reference identifiers bundled as
compact scalar summaries.
"""
published_hybrid_reference_ids() = (_PREISSER2020_SHOCK_REFERENCE_ID,)

"""
    published_hybrid_reference_metadata(id=:preisser2020_65deg_Bavg_y)

Return provenance for a bundled published external hybrid-code reference:
publication DOI, license, source HDF5 filename, source checksum, and derivation
notes. The full upstream data are not bundled.
"""
function published_hybrid_reference_metadata(id::Symbol = _PREISSER2020_SHOCK_REFERENCE_ID)
    id == _PREISSER2020_SHOCK_REFERENCE_ID || throw(_unknown_published_reference_error(id))
    return _PREISSER2020_SHOCK_REFERENCE_META
end

function _require_known_alpha_fraction(alpha_fraction::Real)
    a = Float64(alpha_fraction)
    isfinite(a) || throw(ArgumentError("alpha_fraction must be finite"))
    for ref in _PREISSER2020_BAVG_REFERENCES
        if a == ref.alpha_fraction
            return ref
        end
    end
    throw(
        ArgumentError(
            "unsupported alpha_fraction $(alpha_fraction); expected one of 0.01, 0.05, 0.10",
        ),
    )
end

"""
    published_hybrid_reference(; id=:preisser2020_65deg_Bavg_y, alpha_fraction=0.01)

Return scalar observables from the bundled published external hybrid-code
reference. For the Preisser et al. Zenodo reference, `alpha_fraction` selects the
1%, 5%, or 10% He++ 65-degree shock `Bavg_y` profile and returns the profile
sample count plus min/max/mean magnetic-field magnitude summary.
"""
function published_hybrid_reference(;
    id::Symbol = _PREISSER2020_SHOCK_REFERENCE_ID,
    alpha_fraction::Real = 0.01,
)
    id == _PREISSER2020_SHOCK_REFERENCE_ID || throw(_unknown_published_reference_error(id))
    return _require_known_alpha_fraction(alpha_fraction)
end

"""
    compare_to_published_hybrid_reference(measured; id=:preisser2020_65deg_Bavg_y,
                                          alpha_fraction=0.01, rtol=0.1, atol=0.0)

Compare `measured` scalar observables against a bundled published external
hybrid-code reference using [`compare_to_reference`](@ref). Extra fields in
`measured` are ignored; required fields are `thetaBn_deg`, `alpha_fraction`,
`nsamples`, `Bavg_y_min`, `Bavg_y_max`, and `Bavg_y_mean`.
"""
function compare_to_published_hybrid_reference(
    measured::NamedTuple;
    id::Symbol = _PREISSER2020_SHOCK_REFERENCE_ID,
    alpha_fraction::Real = 0.01,
    rtol::Real = 0.1,
    atol::Real = 0.0,
)
    reference = published_hybrid_reference(; id, alpha_fraction)
    numeric_reference = (;
        thetaBn_deg = reference.thetaBn_deg,
        alpha_fraction = reference.alpha_fraction,
        nsamples = reference.nsamples,
        Bavg_y_min = reference.Bavg_y_min,
        Bavg_y_max = reference.Bavg_y_max,
        Bavg_y_mean = reference.Bavg_y_mean,
    )
    comparison = compare_to_reference(measured, numeric_reference; rtol, atol)
    return (;
        pass = comparison.pass,
        measured,
        reference,
        metadata = published_hybrid_reference_metadata(id),
        comparison,
    )
end

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
    # established kinetic ordering: 1 < n₂ < X_RH — the collisionless compression
    # sits below the run's fluid Rankine–Hugoniot value. r.X_rh is evaluated at
    # the realized Mach from the mass-conservation shock speed (not a front fit),
    # so it is a robust per-run ceiling; gating on the M→∞ limit (γe+1)/(γe−1)
    # instead would leave compressions in (X_RH, (γe+1)/(γe−1)) — physically
    # impossible at the run's Mach number — undetected.
    ordering_ok = isfinite(r.n2) && isfinite(r.X_rh) && 1 < r.n2 < r.X_rh * (1 + rtolT)
    return (; pass = cmp.pass && ordering_ok, measured = r, reference, comparison = cmp)
end
