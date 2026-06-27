# validation.jl - runtime guards that must remain active in optimized runs.

function _require_all_finite(name::AbstractString, a, context::AbstractString)
    all(isfinite, a) || error("$(name) went non-finite ($(context))")
    return nothing
end
