# Adaptive selection uses the existing base-field curvature oracle. Only the
# mesh transaction owns applying these tags to a physical PIC state.
function _dce_tagged_refinements(field::Vector{T}, threshold::T, buffer::Integer) where {T}
    isfinite(threshold) && _stored_nonnegative(threshold) ||
        throw(ArgumentError("adaptive composite threshold must be finite and nonnegative"))
    0 <= buffer <= typemax(Int) ||
        throw(ArgumentError("adaptive composite buffer must be a nonnegative Int-sized integer"))
    all(isfinite, field) || throw(ArgumentError("adaptive base field must be finite"))
    flags = refine_flags(field, threshold)
    expanded = copy(flags)
    n, margin = length(flags), Int(buffer)
    for i in eachindex(flags)
        flags[i] || continue
        # Bound each distance before addition/subtraction, including huge margins.
        first = i - min(margin, i - 1)
        last = i + min(margin, n - i)
        fill!(view(expanded, first:last), true)
    end
    specs = Tuple{Int,Int,Int}[]
    first = 0
    for i in eachindex(expanded)
        if expanded[i]
            iszero(first) && (first = i)
        elseif !iszero(first)
            push!(specs, (1, first, i - 1))
            first = 0
        end
    end
    iszero(first) || push!(specs, (1, first, n))
    return Tuple(specs)
end
