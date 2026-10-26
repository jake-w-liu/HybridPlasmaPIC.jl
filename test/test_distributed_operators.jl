using HybridPlasmaPIC, MPI, Test
import PencilArrays
import PencilFFTs

# Compilation dominates this suite. Helpers therefore neither specialize nor
# infer on function or argument-tuple values, vector fields are explicit three-tuples
# (never `ntuple` with a runtime length), loops never iterate over tuples of
# heterogeneous argument tuples, and component reductions are written out.

_triple(f) = (f(1), f(2), f(3))

# A recoverable rank-local failure must reach every rank before any transpose:
# the origin keeps its own exception and peers receive the stage ErrorException.
Base.@nospecializeinfer function _rejected_call(
    @nospecialize(operation),
    @nospecialize(arguments),
    @nospecialize(keywords),
    origin,
    comm,
)
    rank = MPI.Comm_rank(comm)
    caught = try
        operation(arguments...; keywords...)
        nothing
    catch err
        err
    end
    @test caught isa (rank == origin ? ArgumentError : ErrorException)
    @test MPI.Allreduce(1, +, comm) == MPI.Comm_size(comm)
    return nothing
end

# Invalid arguments are used only on `origin`; no output may be published.
Base.@nospecializeinfer function _check_rejected!(
    @nospecialize(operation),
    @nospecialize(valid),
    @nospecialize(invalid),
    @nospecialize(keywords),
    @nospecialize(outputs),
    sentinel,
    origin,
    comm,
)
    for output in outputs
        fill!(parent(output), sentinel)
    end
    arguments = MPI.Comm_rank(comm) == origin ? invalid : valid
    _rejected_call(operation, arguments, keywords, origin, comm)
    @test all(output -> all(==(sentinel), parent(output)), outputs)
end

# In-place variant: the operated fields must keep their exact contents.
Base.@nospecializeinfer function _check_rejected_in_place(
    @nospecialize(operation),
    @nospecialize(valid),
    @nospecialize(invalid),
    @nospecialize(keywords),
    @nospecialize(fields),
    origin,
    comm,
)
    snapshot = Any[copy(parent(field)) for field in fields]
    arguments = MPI.Comm_rank(comm) == origin ? invalid : valid
    _rejected_call(operation, arguments, keywords, origin, comm)
    @test all(i -> isequal(parent(fields[i]), snapshot[i]), eachindex(snapshot))
end

# Fill a rank-owned field from a function of zero-based global logical indices.
Base.@nospecializeinfer function _fill_global!(field, @nospecialize(value))
    view = PencilArrays.global_view(field)
    for I in CartesianIndices(view)
        view[I] = value(Tuple(I) .- 1)
    end
    return field
end

# Communicator-wide maximum of |field - expected| against Float64 expectations.
Base.@nospecializeinfer function _global_error(field, @nospecialize(expected), comm)
    view = PencilArrays.global_view(field)
    local_error = 0.0
    for I in CartesianIndices(view)
        local_error = max(local_error, abs(Float64(view[I]) - expected(Tuple(I) .- 1)))
    end
    return MPI.Allreduce(local_error, MPI.MAX, comm)
end

function _global_difference(a, b, comm)
    local_error = 0.0
    for I in eachindex(parent(a), parent(b))
        local_error = max(local_error, abs(Float64(parent(a)[I]) - Float64(parent(b)[I])))
    end
    return MPI.Allreduce(local_error, MPI.MAX, comm)
end

_global_magnitude(field, comm) = MPI.Allreduce(Float64(maximum(abs, parent(field))), MPI.MAX, comm)

_global_magnitude(v::Tuple, comm) =
    max(_global_magnitude(v[1], comm), _global_magnitude(v[2], comm), _global_magnitude(v[3], comm))

_vector_difference(a, b, comm) = max(
    _global_difference(a[1], b[1], comm),
    _global_difference(a[2], b[2], comm),
    _global_difference(a[3], b[3], comm),
)

_snapshot(v) = (copy(parent(v[1])), copy(parent(v[2])), copy(parent(v[3])))

_unchanged(v, snapshot) =
    isequal(parent(v[1]), snapshot[1]) &&
    isequal(parent(v[2]), snapshot[2]) &&
    isequal(parent(v[3]), snapshot[3])

_filled(v, value) =
    all(==(value), parent(v[1])) && all(==(value), parent(v[2])) && all(==(value), parent(v[3]))

# Deterministic broadband data containing every mode, including DC and Nyquist,
# defined on global indices so it is independent of the rank layout.
_broadband(seed) =
    i -> sin(0.37 * (i[1] + 1) * seed + 1.1 * (i[2] + 2)^2 + 0.7 * i[3] * (i[1] + 3) + seed)

# Largest Nyquist-zeroed derivative wavenumber per axis, computed independently
# from the signed modes -((n - 1) ÷ 2):(n ÷ 2).
function _derivative_gains(n, L)
    gain(d) = 2pi * maximum(m -> 2m == n[d] ? 0 : abs(m), (-((n[d] - 1) ÷ 2)):(n[d]÷2)) / L[d]
    return (gain(1), gain(2), gain(3))
end

# Local sums for <phi, div v> + <grad phi, v> = 0 and its normwise scale.
function _integration_by_parts_sums(phi, div, grad, v)
    left = right = scale = 0.0
    for I in eachindex(parent(phi))
        a = Float64(parent(phi)[I]) * Float64(parent(div)[I])
        b = (
            Float64(parent(grad[1])[I]) * Float64(parent(v[1])[I]),
            Float64(parent(grad[2])[I]) * Float64(parent(v[2])[I]),
            Float64(parent(grad[3])[I]) * Float64(parent(v[3])[I]),
        )
        left += a
        right += b[1] + b[2] + b[3]
        scale += abs(a) + abs(b[1]) + abs(b[2]) + abs(b[3])
    end
    return [left, right, scale]
end

@testset "Distributed gradient, divergence, curl and Laplacian modes" begin
    ensure_mpi_initialized!()
    comm = MPI.COMM_WORLD
    for T in (Float32, Float64), n in ((8, 6, 4), (5, 7, 9), (5, 4, 3))
        L = (T(2), T(3), T(5))
        g = distributed_fourier_grid(n, L; comm)
        k = (2pi / Float64(L[1]), 2pi / Float64(L[2]), 2pi / Float64(L[3]))
        kmax, ksum = maximum(k), sum(k)
        θ(i, d) = 2pi * i[d] / n[d]
        parity(i, d) = iseven(i[d]) ? 1.0 : -1.0
        tolerance(scale) = 64 * eps(T) * scale
        field(value) = _fill_global!(distributed_field(g), value)
        fields() = (distributed_field(g), distributed_field(g), distributed_field(g))

        f = field(i -> sin(θ(i, 1)) * cos(θ(i, 2)) + cos(θ(i, 3)))
        original = copy(parent(f))
        grad = fields()
        @test gradient!(grad, f, g) === grad
        @test parent(f) == original
        @test _global_error(grad[1], i -> k[1] * cos(θ(i, 1)) * cos(θ(i, 2)), comm) <=
              tolerance(2kmax)
        @test _global_error(grad[2], i -> -k[2] * sin(θ(i, 1)) * sin(θ(i, 2)), comm) <=
              tolerance(2kmax)
        @test _global_error(grad[3], i -> -k[3] * sin(θ(i, 3)), comm) <= tolerance(2kmax)
        # Unscaled inputs use the derivative's exact multiplier sequence.
        for j = 1:3
            @test parent(grad[j]) == parent(deriv(f, g, j))
        end
        # One output may alias the input; all components are staged first.
        aliased = (copy(f), distributed_field(g), distributed_field(g))
        gradient!(aliased, aliased[1], g)
        @test _unchanged(aliased, _snapshot(grad))

        v = (
            field(i -> sin(θ(i, 1)) * cos(θ(i, 2))),
            field(i -> cos(θ(i, 2)) * sin(θ(i, 3))),
            field(i -> sin(θ(i, 3)) * cos(θ(i, 1))),
        )
        div = distributed_field(g)
        @test divergence!(div, v, g) === div
        expected_div(i) =
            k[1] * cos(θ(i, 1)) * cos(θ(i, 2)) - k[2] * sin(θ(i, 2)) * sin(θ(i, 3)) +
            k[3] * cos(θ(i, 3)) * cos(θ(i, 1))
        @test _global_error(div, expected_div, comm) <= tolerance(ksum)
        # Only the first three components are used, and `out` may alias one.
        extra = distributed_field(g)
        divergence!(extra, (v[1], v[2], v[3], f), g)
        @test parent(extra) == parent(div)
        w = (copy(v[1]), copy(v[2]), copy(v[3]))
        divergence!(w[2], w, g)
        @test parent(w[2]) == parent(div)

        A = (
            field(i -> sin(θ(i, 2)) + cos(θ(i, 3))),
            field(i -> cos(θ(i, 1)) + sin(θ(i, 3))),
            field(i -> sin(θ(i, 1)) + cos(θ(i, 2))),
        )
        A_original = _snapshot(A)
        C = fields()
        @test curl!(C, A, g) === C
        @test _unchanged(A, A_original)
        @test _global_error(C[1], i -> -k[2] * sin(θ(i, 2)) - k[3] * cos(θ(i, 3)), comm) <=
              tolerance(2kmax)
        @test _global_error(C[2], i -> -k[3] * sin(θ(i, 3)) - k[1] * cos(θ(i, 1)), comm) <=
              tolerance(2kmax)
        @test _global_error(C[3], i -> -k[1] * sin(θ(i, 1)) - k[2] * cos(θ(i, 2)), comm) <=
              tolerance(2kmax)

        lap = distributed_field(g)
        @test laplacian!(lap, f, g) === lap
        expected_lap(i) = -(k[1]^2 + k[2]^2) * sin(θ(i, 1)) * cos(θ(i, 2)) - k[3]^2 * cos(θ(i, 3))
        @test _global_error(lap, expected_lap, comm) <= tolerance(2 * sum(abs2, k))
        in_place = copy(f)
        @test laplacian!(in_place, in_place, g) === in_place
        @test parent(in_place) == parent(lap)

        zero_field = distributed_field(g)
        for d = 1:3
            iseven(n[d]) || continue
            e = mod1(d + 1, 3)
            nyquist = field(i -> parity(i, d))
            knyquist = pi * n[d] / Float64(L[d])
            # The Laplacian retains Nyquist modes.
            laplacian!(lap, nyquist, g)
            @test _global_error(lap, i -> -knyquist^2 * parity(i, d), comm) <= tolerance(knyquist^2)
            # Odd derivatives drop the Nyquist axis only (mixed Nyquist modes).
            gradient!(grad, field(i -> parity(i, d) * sin(θ(i, e))), g)
            @test _global_error(grad[d], i -> 0.0, comm) <= tolerance(kmax)
            @test _global_error(grad[e], i -> parity(i, d) * k[e] * cos(θ(i, e)), comm) <=
                  tolerance(kmax)
            @test _global_error(grad[6-d-e], i -> 0.0, comm) <= tolerance(kmax)
            components = _triple(c -> c == d ? nyquist : zero_field)
            divergence!(div, components, g)
            @test _global_magnitude(div, comm) <= tolerance(ksum)
            curl!(C, _triple(c -> c == e ? nyquist : zero_field), g)
            @test _global_magnitude(C, comm) <= tolerance(ksum)
        end

        # Discrete identities on broadband data, normwise at each operation scale.
        gains = _derivative_gains(n, L)
        broadband = _triple(c -> field(_broadband(c)))
        curl!(C, broadband, g)
        divergence!(div, C, g)
        @test _global_magnitude(div, comm) <= tolerance(sum(gains) * _global_magnitude(C, comm))
        phi = field(_broadband(4))
        gradient!(grad, phi, g)
        curl!(C, grad, g)
        @test _global_magnitude(C, comm) <= tolerance(sum(gains) * _global_magnitude(grad, comm))
        # Integration by parts: <phi, div v> = -<grad phi, v>.
        divergence!(div, broadband, g)
        sums = MPI.Allreduce(_integration_by_parts_sums(phi, div, grad, broadband), +, comm)
        @test abs(sums[1] + sums[2]) <= tolerance(sums[3])
        # Laplacian equals div(grad) for band-limited fields (no Nyquist content).
        band = dealias_two_thirds!(field(_broadband(5)), g)
        gradient!(grad, band, g)
        divergence!(div, grad, g)
        laplacian!(lap, band, g)
        lap_scale =
            sum(gains) * _global_magnitude(grad, comm) +
            sum(abs2, gains) * _global_magnitude(band, comm)
        @test _global_difference(div, lap, comm) <= tolerance(lap_scale)
    end
end

@testset "Distributed divergence-free projection" begin
    comm = MPI.COMM_WORLD
    for T in (Float32, Float64), n in ((8, 6, 4), (5, 7, 9))
        L = (T(2), T(3), T(5))
        g = distributed_fourier_grid(n, L; comm)
        k = (2pi / Float64(L[1]), 2pi / Float64(L[2]), 2pi / Float64(L[3]))
        θ(i, d) = 2pi * i[d] / n[d]
        parity(i, d) = iseven(i[d]) ? 1.0 : -1.0
        tolerance(scale) = 64 * eps(T) * scale
        field(value) = _fill_global!(distributed_field(g), value)
        # Solenoidal modes, a potential gradient, DC, and on even axis 1 a
        # derivative-null Nyquist mode plus a mixed Nyquist/k_y mode whose x
        # part is orthogonal to its derivative-resolved wavenumber.
        solenoidal = (i -> sin(θ(i, 2)), i -> sin(θ(i, 3)), i -> sin(θ(i, 1)))
        potential = (
            i -> -k[1] * sin(θ(i, 1)) * cos(θ(i, 2)) * cos(θ(i, 3)),
            i -> -k[2] * cos(θ(i, 1)) * sin(θ(i, 2)) * cos(θ(i, 3)),
            i -> -k[3] * cos(θ(i, 1)) * cos(θ(i, 2)) * sin(θ(i, 3)),
        )
        mean = (0.5, -0.25, 0.75)
        even = iseven(n[1])
        null(i, c) = even && c != 3 ? 0.3 * parity(i, 1) : 0.0
        mixed(i, c) = even && c != 3 ? parity(i, 1) * sin(θ(i, 2)) : 0.0
        kept(i, c) = solenoidal[c](i) + mean[c] + null(i, c) + (c == 1 ? mixed(i, 1) : 0.0)
        B = _triple(
            c -> field(
                i -> solenoidal[c](i) + potential[c](i) + mean[c] + null(i, c) + mixed(i, c),
            ),
        )
        @test project_divfree!(B, g) === B
        for c = 1:3
            @test _global_error(B[c], i -> kept(i, c), comm) <= tolerance(4 * maximum(k))
        end

        gains = _derivative_gains(n, L)
        R = _triple(c -> field(_broadband(c)))
        project_divfree!(R, g)
        div = divergence!(distributed_field(g), R, g)
        scale = _global_magnitude(R, comm)
        @test _global_magnitude(div, comm) <= tolerance(sum(gains) * scale)
        projected = (copy(R[1]), copy(R[2]), copy(R[3]))
        project_divfree!(R, g)
        @test _vector_difference(R, projected, comm) <= tolerance(scale)

        # The projector is dimensionless: isotropically rescaled domains agree,
        # including lengths whose squared wavenumbers are not representable.
        extremes = T === Float32 ? (T(1e-30), T(1e30)) : (T(1e-250), T(1e250))
        for factor in extremes
            h = distributed_fourier_grid(n, (factor * L[1], factor * L[2], factor * L[3]); comm)
            S = _triple(c -> _fill_global!(distributed_field(h), _broadband(c)))
            project_divfree!(S, h)
            @test _vector_difference(S, projected, comm) <= tolerance(scale)
        end
        # Finite near-overflow amplitudes are rescaled before the transform.
        amplitude = floatmax(T) / T(64)
        H = _triple(c -> field(i -> amplitude * T(_broadband(c)(i))))
        project_divfree!(H, g)
        for c = 1:3
            @test all(isfinite, parent(H[c]))
            local_error = 0.0
            for I in eachindex(parent(H[c]))
                expected = Float64(amplitude) * Float64(parent(projected[c])[I])
                local_error = max(local_error, abs(Float64(parent(H[c])[I]) - expected))
            end
            @test MPI.Allreduce(local_error, MPI.MAX, comm) <= tolerance(Float64(amplitude) * scale)
        end
    end
end

@testset "Distributed spectral filters" begin
    comm = MPI.COMM_WORLD
    for T in (Float32, Float64)
        n = (6, 9, 4)
        g = distributed_fourier_grid(n, (T(2), T(3), T(5)); comm)
        c(i, d, m) = cos(2pi * m * i[d] / n[d])
        s(i, d, m) = sin(2pi * m * i[d] / n[d])
        tolerance(scale) = 64 * eps(T) * scale
        field(value) = _fill_global!(distributed_field(g), value)
        # Strict 2/3 boundary: 3|m| = n is removed on axes divisible by three.
        retained(i) = 0.75 + c(i, 1, 1) + s(i, 2, 2) + c(i, 3, 1) + c(i, 1, 1) * s(i, 2, 2)
        removed(i) = c(i, 1, 2) + c(i, 2, 3) + c(i, 3, 2) + c(i, 1, 1) * c(i, 2, 3) + s(i, 2, 4)
        h = field(i -> retained(i) + removed(i))
        @test dealias_two_thirds!(h, g) === h
        @test _global_error(h, retained, comm) <= tolerance(8)

        # Mode-ratio filter against an independent transfer-function evaluation.
        modes = ((0, 0, 0), (1, 0, 0), (0, 2, 0), (0, 0, 1), (1, 2, 0), (3, 0, 0), (0, 4, 2))
        values(i) = (
            0.75,
            c(i, 1, 1),
            s(i, 2, 2),
            c(i, 3, 1),
            c(i, 1, 1) * s(i, 2, 2),
            c(i, 1, 3),
            s(i, 2, 4) * c(i, 3, 2),
        )
        transfer(m, α, p) = exp(-α * max((m[1] / 3)^2, (m[2] / 4)^2, (m[3] / 2)^2)^p)
        signal(i) = sum(values(i))
        filtered(i, α, p) = sum(transfer.(modes, α, p) .* values(i))
        for (α, p) in ((2, 1), (36, 8), (0.5, 3))
            h = field(signal)
            @test exp_filter!(h, g; α, p) === h
            @test _global_error(h, i -> filtered(i, α, p), comm) <= tolerance(8)
        end
        h = field(signal)
        @test exp_filter!(h, g) === h
        @test _global_error(h, i -> filtered(i, 36, 8), comm) <= tolerance(8)
        # Equivalent parameter types convert to the same multiplier.
        reference = exp_filter!(field(signal), g; α = 2, p = 1)
        @test parent(exp_filter!(field(signal), g; α = 2 // 1, p = Int32(1))) == parent(reference)
        # α = Inf keeps only the mean; α = 0 is the identity up to roundoff.
        h = field(signal)
        exp_filter!(h, g; α = Inf, p = 8)
        @test _global_error(h, i -> 0.75, comm) <= tolerance(8)
        h = field(signal)
        exp_filter!(h, g; α = 0, p = 8)
        @test _global_error(h, signal, comm) <= tolerance(8)
    end
end

@testset "Distributed operators agree with serial FourierGrid" begin
    comm = MPI.COMM_WORLD
    rank, ranks = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    dims = MPI.Dims_create(ranks, (Cint(0), Cint(0)))
    # PencilFFTs rejects empty partitions, so singleton axes exist only where
    # the process grid is one-wide in the decomposed directions.
    shapes = [(5, 4, 3), (8, 6, 4)]
    dims[2] == 1 && push!(shapes, (6, 4, 1))
    dims[1] == 1 && append!(shapes, ((1, 4, 6), (1, 1, 1)))
    for T in (Float32, Float64), n in shapes
        L = (T(2), T(3), T(5))
        g = distributed_fourier_grid(n, L; comm)
        f = _fill_global!(distributed_field(g), _broadband(1))
        A = _triple(c -> _fill_global!(distributed_field(g), _broadband(c + 1)))
        out = distributed_field(g)
        O = (distributed_field(g), distributed_field(g), distributed_field(g))
        full_f = PencilArrays.gather(f)
        full_A = map(PencilArrays.gather, A)
        serial = rank == 0 ? FourierGrid(n, L) : nothing
        # Test-only gathered comparison: FFT orderings differ, so agreement is
        # normwise at the serial output scale rather than bitwise. Every rank
        # calls the collective gather; results are combined without short-circuit.
        function agree(distributed, reference)
            gathered = PencilArrays.gather(distributed)
            rank == 0 || return true
            scale = max(maximum(abs, reference), one(T))
            return maximum(abs, gathered .- reference) <= 32 * eps(T) * scale
        end
        agree(distributed::Tuple, reference) =
            agree(distributed[1], rank == 0 ? reference[1] : nothing) &
            agree(distributed[2], rank == 0 ? reference[2] : nothing) &
            agree(distributed[3], rank == 0 ? reference[3] : nothing)
        serial_vector() = (similar(full_f), similar(full_f), similar(full_f))
        gradient!(O, f, g)
        @test agree(O, rank == 0 ? gradient!(serial_vector(), full_f, serial) : nothing)
        divergence!(out, A, g)
        @test agree(out, rank == 0 ? divergence!(similar(full_f), full_A, serial) : nothing)
        curl!(O, A, g)
        @test agree(O, rank == 0 ? curl!(serial_vector(), full_A, serial) : nothing)
        laplacian!(out, f, g)
        @test agree(out, rank == 0 ? laplacian!(similar(full_f), full_f, serial) : nothing)
        B = (copy(A[1]), copy(A[2]), copy(A[3]))
        project_divfree!(B, g)
        @test agree(B, rank == 0 ? project_divfree!(map(copy, full_A), serial) : nothing)
        for (α, p) in ((2, 1), (36, 8))
            h = exp_filter!(copy(f), g; α, p)
            @test agree(h, rank == 0 ? exp_filter!(copy(full_f), serial; α, p) : nothing)
        end
        h = dealias_two_thirds!(copy(f), g)
        @test agree(h, rank == 0 ? dealias_two_thirds!(copy(full_f), serial) : nothing)
    end
end

@testset "Distributed operator range boundaries" begin
    comm = MPI.COMM_WORLD
    rank, ranks = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    for T in (Float32, Float64)
        # Power-of-two axes transform constants to exactly zero non-DC modes.
        n = (8, 4, 4)
        g = distributed_fourier_grid(n, (T(2), T(3), T(5)); comm)
        huge = floatmax(T) / T(2)
        sentinel = T(5)
        constant() = fill!(distributed_field(g), huge)
        marked() = fill!(distributed_field(g), sentinel)
        near_huge(x) = all(y -> isapprox(y, huge; rtol = 4eps(T)), parent(x))
        @test _filled(gradient!((marked(), marked(), marked()), constant(), g), zero(T))
        @test all(iszero, parent(divergence!(marked(), (constant(), constant(), constant()), g)))
        vector = (constant(), constant(), constant())
        @test _filled(curl!((marked(), marked(), marked()), vector, g), zero(T))
        @test all(iszero, parent(laplacian!(marked(), constant(), g)))
        B = project_divfree!(vector, g)
        @test near_huge(B[1]) && near_huge(B[2]) && near_huge(B[3])
        @test near_huge(exp_filter!(constant(), g))
        @test near_huge(dealias_two_thirds!(constant(), g))

        # A nonzero lost to the common scale fails before the next transpose.
        mixed = fill!(distributed_field(g), one(T))
        rank == 0 && (parent(mixed)[1] = nextfloat(zero(T)))
        rank == ranks - 1 && (parent(mixed)[end] = huge)
        unit = fill!(distributed_field(g), one(T))
        out = marked()
        O = (marked(), marked(), marked())
        _rejected_call(gradient!, (O, mixed, g), (;), 0, comm)
        _rejected_call(divergence!, (out, (unit, mixed, unit), g), (;), 0, comm)
        _rejected_call(curl!, (O, (unit, mixed, unit), g), (;), 0, comm)
        _rejected_call(laplacian!, (out, mixed, g), (;), 0, comm)
        @test all(==(sentinel), parent(out)) && _filled(O, sentinel)
        B = (copy(unit), copy(mixed), copy(unit))
        B_original = _snapshot(B)
        _rejected_call(project_divfree!, (B, g), (;), 0, comm)
        @test _unchanged(B, B_original)
        h = copy(mixed)
        _rejected_call(exp_filter!, (h, g), (;), 0, comm)
        _rejected_call(dealias_two_thirds!, (h, g), (;), 0, comm)
        @test isequal(parent(h), parent(mixed))

        # Squared wavenumbers overflow: Laplacian rejects on every rank, while
        # first-derivative operators still rescale huge constants to zero.
        tiny = T === Float32 ? T(1e-30) : T(1e-250)
        small = distributed_fourier_grid(n, (tiny, tiny, tiny); comm)
        x = fill!(distributed_field(small), huge)
        y = fill!(distributed_field(small), sentinel)
        Y = (distributed_field(small), distributed_field(small), distributed_field(small))
        mark!(v) = (foreach(z -> fill!(parent(z), sentinel), v); v)
        @test_throws ArgumentError laplacian!(y, x, small)
        @test all(==(sentinel), parent(y))
        @test MPI.Allreduce(1, +, comm) == ranks
        @test _filled(gradient!(mark!(Y), x, small), zero(T))
        @test all(iszero, parent(divergence!(y, (x, x, x), small)))
        @test _filled(curl!(mark!(Y), (x, x, x), small), zero(T))

        # The summed divergence/curl bound overflows before the gradient bound.
        edge = T === Float32 ? T(4.5e-38) : T(1e-307)
        narrow = distributed_fourier_grid((4, 4, 4), (edge, edge, edge); comm)
        x = fill!(distributed_field(narrow), one(T))
        y = fill!(distributed_field(narrow), sentinel)
        Y = mark!((distributed_field(narrow), distributed_field(narrow), distributed_field(narrow)))
        @test isinf(sum(map(T, _derivative_gains((4, 4, 4), (edge, edge, edge)))))
        @test_throws ArgumentError divergence!(y, (x, x, x), narrow)
        @test_throws ArgumentError curl!(Y, (x, x, x), narrow)
        @test_throws ArgumentError laplacian!(y, x, narrow)
        @test all(==(sentinel), parent(y)) && _filled(Y, sentinel)
        @test MPI.Allreduce(1, +, comm) == ranks
        @test _filled(gradient!(Y, x, narrow), zero(T))

        # Squared wavenumbers underflow, but the ordered product stays representable.
        wide, amplitude = T === Float32 ? (T(1e30), T(1e30)) : (T(1e200), T(1e200))
        large = distributed_fourier_grid(n, (wide, wide, wide); comm)
        x = _fill_global!(distributed_field(large), i -> amplitude * T(cos(2pi * i[1] / n[1])))
        y = laplacian!(distributed_field(large), x, large)
        kwide = 2 * big(pi) / big(wide)
        peak = Float64(kwide * big(amplitude) * kwide)
        @test peak > 0 && iszero(abs2(T(kwide)))
        @test _global_error(y, i -> -peak * cos(2pi * i[1] / n[1]), comm) <= 64 * eps(T) * peak
    end
end

@testset "Distributed operator failure and publication contract" begin
    comm = MPI.COMM_WORLD
    rank, ranks = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    origins = unique((0, ranks ÷ 2, ranks - 1))
    n = (6, 4, 4)
    for T in (Float32, Float64)
        L = (T(2), T(3), T(5))
        g = distributed_fourier_grid(n, L; comm)
        other = distributed_fourier_grid(n, L; comm)
        f = _fill_global!(distributed_field(g), _broadband(1))
        A = _triple(c -> _fill_global!(distributed_field(g), _broadband(c + 1)))
        sentinel = T(7)
        out = distributed_field(g)
        O = (distributed_field(g), distributed_field(g), distributed_field(g))
        outputs = (out, O[1], O[2], O[3])
        # Operand replacements used only on the origin rank.
        pencil = distributed_field(other)
        precision = fill!(similar(f, T === Float32 ? Float64 : Float32), 1)
        workspace = g.staging_work[1]
        poisoned(x) = (y = copy(x); parent(y)[end] = T(NaN); y)
        filter_keywords = (α = 2, p = 1)
        # Valid references; the reused workspace reproduces them bitwise.
        grad_reference = _snapshot(gradient!(O, f, g))
        div_reference = copy(parent(divergence!(out, A, g)))
        curl_reference = _snapshot(curl!(O, A, g))
        lap_reference = copy(parent(laplacian!(out, f, g)))
        projection_reference = _snapshot(project_divfree!(map(copy, A), g))
        filter_reference = copy(parent(exp_filter!(copy(f), g; filter_keywords...)))
        dealias_reference = copy(parent(dealias_two_thirds!(copy(f), g)))
        for origin in origins
            reject!(operation, valid, invalid) =
                _check_rejected!(operation, valid, invalid, (;), outputs, sentinel, origin, comm)

            valid = (O, f, g)
            reject!(gradient!, valid, (O, poisoned(f), g))
            reject!(gradient!, valid, (O, pencil, g))
            reject!(gradient!, valid, ((O[1], precision, O[3]), f, g))
            reject!(gradient!, valid, ((O[1], O[2], workspace), f, g))
            reject!(gradient!, valid, ((O[1], O[1], O[3]), f, g))
            reject!(gradient!, valid, ((O[1], O[2]), f, g))
            @test _unchanged(gradient!(O, f, g), grad_reference)

            valid = (out, A, g)
            reject!(divergence!, valid, (out, (A[1], poisoned(A[2]), A[3]), g))
            reject!(divergence!, valid, (out, (A[1], A[2], precision), g))
            reject!(divergence!, valid, (pencil, A, g))
            reject!(divergence!, valid, (workspace, A, g))
            reject!(divergence!, valid, (out, (A[1], A[2]), g))
            @test parent(divergence!(out, A, g)) == div_reference

            valid = (O, A, g)
            reject!(curl!, valid, (O, (poisoned(A[1]), A[2], A[3]), g))
            reject!(curl!, valid, (O, (A[1], workspace, A[3]), g))
            reject!(curl!, valid, ((O[1], O[2], pencil), A, g))
            reject!(curl!, valid, ((precision, O[2], O[3]), A, g))
            reject!(curl!, valid, ((O[1], O[1], O[3]), A, g))
            reject!(curl!, valid, ((O[1], A[3], O[3]), A, g))
            reject!(curl!, valid, ([O[1], O[2], O[3]], A, g))
            @test _unchanged(curl!(O, A, g), curl_reference)

            valid = (out, f, g)
            reject!(laplacian!, valid, (out, poisoned(f), g))
            reject!(laplacian!, valid, (out, pencil, g))
            reject!(laplacian!, valid, (precision, f, g))
            reject!(laplacian!, valid, (workspace, f, g))
            @test parent(laplacian!(out, f, g)) == lap_reference

            for kind in (:nan, :pencil, :type, :workspace, :alias, :shape)
                B = map(copy, A)
                replacement =
                    kind === :nan ? poisoned(B[2]) :
                    kind === :pencil ? pencil : kind === :type ? precision : workspace
                invalid =
                    kind === :alias ? (B[1], B[1], B[3]) :
                    kind === :shape ? (B[1], B[2]) : (B[1], replacement, B[3])
                _check_rejected_in_place(
                    project_divfree!,
                    (B, g),
                    (invalid, g),
                    (;),
                    B,
                    origin,
                    comm,
                )
            end
            @test _unchanged(project_divfree!(map(copy, A), g), projection_reference)

            for replacement in Any[poisoned(f), pencil, precision, workspace]
                h = copy(f)
                for (operation, keywords) in
                    ((exp_filter!, filter_keywords), (dealias_two_thirds!, (;)))
                    valid, invalid = (h, g), (replacement, g)
                    _check_rejected_in_place(
                        operation,
                        valid,
                        invalid,
                        keywords,
                        (h,),
                        origin,
                        comm,
                    )
                end
            end
            @test parent(exp_filter!(copy(f), g; filter_keywords...)) == filter_reference
            @test parent(dealias_two_thirds!(copy(f), g)) == dealias_reference
            for keywords in Any[
                (α = -1, p = 1),
                (α = NaN, p = 1),
                (α = "2", p = 1),
                (α = 2, p = 0),
                (α = 2, p = 1.5),
                (α = 2, p = big(typemax(Int)) + 1),
            ]
                h = copy(f)
                selected = rank == origin ? keywords : filter_keywords
                _rejected_call(exp_filter!, (h, g), selected, origin, comm)
                @test parent(h) == parent(f)
            end
        end
        # Individually valid filter parameters must agree on every rank.
        if ranks > 1
            for keywords in Any[(α = 3, p = 1), (α = 2, p = 2)]
                h = copy(f)
                selected = rank == ranks - 1 ? keywords : filter_keywords
                _rejected_call(exp_filter!, (h, g), selected, ranks - 1, comm)
                @test parent(h) == parent(f)
            end
        end
    end
end

@testset "Distributed operator workspace allocations" begin
    comm = MPI.COMM_WORLD
    for T in (Float32, Float64)
        g = distributed_fourier_grid((8, 6, 4), (T(2), T(3), T(5)); comm)
        f = _fill_global!(distributed_field(g), _broadband(1))
        A = _triple(c -> _fill_global!(distributed_field(g), _broadband(c + 1)))
        B = map(copy, A)
        out = copy(f)
        O = (distributed_field(g), distributed_field(g), distributed_field(g))
        # Warm calls allocate only inside PencilFFTs transforms and in the
        # scalar collectives of each stage; workspaces belong to the grid.
        HybridPlasmaPIC.distributed_fft_forward!(g.spectral, g.fft, g.physical)
        forward = @allocated HybridPlasmaPIC.distributed_fft_forward!(g.spectral, g.fft, g.physical)
        HybridPlasmaPIC.distributed_fft_inverse!(g.physical, g.fft, g.spectral)
        inverse = @allocated HybridPlasmaPIC.distributed_fft_inverse!(g.physical, g.fft, g.spectral)
        overhead(bytes, forwards, inverses) = bytes - forwards * forward - inverses * inverse
        gradient!(O, f, g)
        @test overhead((@allocated gradient!(O, f, g)), 1, 3) <= 1024
        divergence!(out, A, g)
        @test overhead((@allocated divergence!(out, A, g)), 3, 1) <= 1024
        curl!(O, A, g)
        @test overhead((@allocated curl!(O, A, g)), 3, 3) <= 1024
        laplacian!(out, f, g)
        @test overhead((@allocated laplacian!(out, f, g)), 1, 1) <= 1024
        project_divfree!(B, g)
        @test overhead((@allocated project_divfree!(B, g)), 3, 3) <= 1024
        exp_filter!(out, g; α = 2, p = 1)
        @test overhead((@allocated exp_filter!(out, g; α = 2, p = 1)), 1, 1) <= 1024
        dealias_two_thirds!(out, g)
        @test overhead((@allocated dealias_two_thirds!(out, g)), 1, 1) <= 1024
    end
end
