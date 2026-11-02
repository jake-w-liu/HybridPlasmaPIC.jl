# Concrete unpreconditioned restarted GMRES storage shared by the shock and
# coupled implicit owners. No physical model or convergence policy lives here.
_gmres_workspace_vectors(s::Krylov.GmresWorkspace) = (s.x, s.w, s.Δx, s.V...)

function _gmres_workspace_arrays(s::Krylov.GmresWorkspace)
    return (
        _gmres_workspace_vectors(s)...,
        s.c,
        s.s,
        s.z,
        s.R,
        s.stats.residuals,
        s.stats.Aresiduals,
        s.stats.qvals,
        s.stats.Acond,
    )
end

function _validate_gmres_workspace_layout(s::Krylov.GmresWorkspace, n::Int; name = "GMRES")
    k = length(s.V)
    s.m == n &&
        s.n == n &&
        1 <= k <= n &&
        length(s.c) == k &&
        length(s.s) == k &&
        length(s.z) == k &&
        length(s.R) == k * (k + 1) ÷ 2 &&
        isempty(s.p) &&
        isempty(s.q) || throw(DimensionMismatch("$name workspace dimensions differ"))
    return nothing
end
