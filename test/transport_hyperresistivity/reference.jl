# Direct real DFT matrices; no production spectral operators or mode metadata.
function th_matrix(n, L, kind)
    modes = [m <= n ÷ 2 ? m : m - n for m = 0:(n-1)]
    return [
        sum(modes) do m
            k = 2pi * m / Float64(L)
            phase = 2pi * m * (i - j) / n
            kind === :lap ? -k^2 * cos(phase) / n :
            kind === :nyquist ? (iseven(n) && m == n ÷ 2 ? k * cos(phase) / n : 0.0) :
            (iseven(n) && m == n ÷ 2 ? 0.0 : -k * sin(phase) / n)
        end for i = 0:(n-1), j = 0:(n-1)
    ]
end
function th_apply(M, a, d)
    out = zeros(Float64, size(a))
    for I in CartesianIndices(a), j in axes(a, d)
        K = CartesianIndex(ntuple(c -> c == d ? j : I[c], ndims(a)))
        out[I] += M[I[d], j] * a[K]
    end
    return out
end
function th_reference(B, p, n, ui, L, γ, η, ηH)
    D = ndims(p)
    d(a, j) = j <= D ? th_apply(th_matrix(size(p, j), L[j], :first), a, j) : zeros(size(p))
    lap(a) = sum(th_apply(th_matrix(size(p, j), L[j], :lap), a, j) for j = 1:D)
    curl(A) = (d(A[3], 2) - d(A[2], 3), d(A[1], 3) - d(A[3], 1), d(A[2], 1) - d(A[1], 2))
    cross(A, B) =
        (A[2] .* B[3] - A[3] .* B[2], A[3] .* B[1] - A[1] .* B[3], A[1] .* B[2] - A[2] .* B[1])
    J = curl(B)
    LJ = map(lap, J)
    gp = ntuple(j -> d(p, j), 3)
    ub = cross(ui, B)
    jb = cross(J, B)
    force = ntuple(c -> -ub[c] + (jb[c] - gp[c]) ./ n, 3)
    E = ntuple(c -> force[c] + η .* J[c] - ηH .* LJ[c], 3)
    ue = ntuple(c -> ui[c] - J[c] ./ n, 3)
    div = sum(d(ue[j], j) for j = 1:D)
    advection = -sum(ue[j] .* gp[j] for j = 1:D) - γ .* p .* div
    ordinary = sum(d(J[c], j) .^ 2 for c = 1:3, j = 1:D)
    nyquist =
        sum(th_apply(th_matrix(size(p, j), L[j], :nyquist), J[c], j) .^ 2 for c = 1:3, j = 1:D)
    heat = η .* sum(a .^ 2 for a in J) + ηH .* (ordinary + nyquist)
    return (;
        J,
        LJ,
        force,
        E,
        outB = map(a -> -a, curl(E)),
        outp = advection + (γ - 1) .* heat,
        heat,
        ordinary,
        nyquist,
    )
end
