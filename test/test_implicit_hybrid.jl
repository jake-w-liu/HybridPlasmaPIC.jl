module ImplicitHybridAcceptanceTests
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC
Q(x) = Rational{BigInt}(x)

# Uniform lattice on the grid nodes: CIC deposits a constant density from any
# translated copy of it, so n, ui, pe, J=curl(B) and E stay spatially uniform
# and every particle of a species sees the same fields. Zero net transverse ui
# is not required; pe = Te*n is uniform so grad(pe) = 0 and J = 0 for a uniform
# seeded B. The represented endpoint is then the fixed point
#   v1 = Boris_h(v0; Emid, B0),  Emid = -ubar x B0,
#   ubar = sum_s q_s n_s (v0+v1)/2 / sum_s q_s n_s,
# solved here by exact rational Picard iteration, which shares no code with the
# production Newton/GMRES solver.
function fixture(::Type{T}; Te = T(1 / 64), masses = (2, 5), N = 8) where {T}
    g = FourierGrid((N,), (one(T),))
    species = [ParticleSet{1,T}(N; q = one(T), m = T(m)) for m in masses]
    dx = g.L[1] / g.n[1]
    for (k, ps) in enumerate(species)
        ps.x[1] .= [T((p - 1) * dx) for p = 1:N]
        fill!(ps.weight, one(T))
        for c = 1:3
            fill!(ps.v[c], T(0.03 * c * k + 0.01 * (c == 3 ? -1 : 1)))
        end
    end
    model = HybridModel(IsothermalElectrons(Te); nfloor = T(1e-6))
    st = ImplicitHybridStepper(g, model, CIC(), species)
    fill!(st.fields.B[3], T(1 / 2))
    return st, species
end

function snapshot(st, species)
    io = IOBuffer()
    serialize(
        io,
        (
            fields = deepcopy((B = st.fields.B, E = st.fields.E, n = st.fields.n)),
            phase = (st.time[], st.step[], st.last_dt[]),
            species = [
                deepcopy((
                    x = ps.x,
                    v = ps.v,
                    weight = ps.weight,
                    id = ps.id,
                    tag = ps.tag,
                    q = ps.q,
                    m = ps.m,
                )) for ps in species
            ],
            initialized = st.initialized[],
            report = st.last_report,
        ),
    )
    return take!(io)
end

cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

function boris_reference(v, alpha, duration, E, B)
    h = alpha * duration / 2
    u = v .+ h .* E
    t = h .* B
    s = (2 / (1 + sum(abs2, t))) .* t
    vp = u .+ cross3(u, t)
    return u .+ cross3(vp, s) .+ h .* E
end

# Exact fixed point of the represented coupled equations on a uniform lattice:
# densities n_s = w*N_s/L are constant, so nbar = sum(q_s n_s), ubar is the
# charge-weighted mean of vbar = (v0+v1)/2, J and grad(pe) vanish, and
# B1 = B0 - h*curl(Emid) = B0. Picard iteration contracts for the small
# gyro-angle steps used here; 200 iterations drive its defect far below the
# solver tolerance, so the remaining difference is the stepper's own residual.
function coupled_oracle(v0, charge_mass, ns, B0, h)
    alpha = [Q(q) / Q(m) for (q, m) in charge_mass]
    weights = [Q(q) * Q(n) for ((q, m), n) in zip(charge_mass, ns)]
    nbar = sum(weights)
    vbar = copy(v0)
    B = Q.(B0)
    for _ = 1:200
        ubar = ntuple(c -> sum(weights[s] * vbar[s][c] for s in eachindex(vbar)) / nbar, 3)
        Emid = .-cross3(ubar, B)
        for s in eachindex(vbar)
            v1 = boris_reference(v0[s], alpha[s], h, Emid, B)
            vbar[s] = (v0[s] .+ v1) ./ 2
        end
    end
    ubar = ntuple(c -> sum(weights[s] * vbar[s][c] for s in eachindex(vbar)) / nbar, 3)
    v1 = [boris_reference(v0[s], alpha[s], h, .-cross3(ubar, B), B) for s in eachindex(vbar)]
    u1 = ntuple(c -> sum(weights[s] * v1[s][c] for s in eachindex(vbar)) / nbar, 3)
    return v1, .-cross3(u1, B)
end

@testset "Direct construction, coupled solve and rational oracle" begin
    for T in (Float32, Float64)
        st, species = fixture(T)
        @test_throws ArgumentError step!(st, species, T(1 / 64))
        init!(st, species)
        @test st.initialized[] && iszero(st.time[]) && st.last_report === nothing
        B0 = deepcopy.(st.fields.B)
        v0 = [Tuple(Q(ps.v[c][1]) for c = 1:3) for ps in species]
        x0 = [Q(ps.x[1][1]) for ps in species]
        h = T(1 / 64)
        step!(st, species, h)
        r = st.last_report
        @test r !== nothing
        @test r.velocity_residual <= r.tolerance && r.magnetic_residual <= r.tolerance
        @test r.newton_iterations >= 0 && r.linear_iterations >= 0 && r.backtracks >= 0
        @test isfinite(r.initial_residual) && r.initial_residual >= 0
        qm = [(ps.q, ps.m) for ps in species]
        ns = [nparticles(ps) for ps in species]
        v1, E1 = coupled_oracle(v0, qm, ns, (B0[1][1], B0[2][1], B0[3][1]), Q(h))
        # Bound: the accepted normalized residual is at most `tolerance`; the
        # Picard oracle converges far below it. 16x covers the residual-to-
        # unknown conditioning of these small systems plus Float rounding.
        tol = 16 * Float64(r.tolerance) + 256eps(Float64)
        for (i, ps) in enumerate(species), c = 1:3
            @test abs(ps.v[c][1] - Float64(v1[i][c])) <= tol
        end
        for (i, ps) in enumerate(species)
            @test abs(
                mod(ps.x[1][1], one(T)) -
                Float64(mod(x0[i] + Q(h) * (v0[i][1] + v1[i][1]) / 2, Q(1))),
            ) <= tol
        end
        # Uniform Emid => curl = 0 => B unchanged to field rounding.
        for c = 1:3
            @test all(isapprox(st.fields.B[c], B0[c]; rtol = 512eps(T), atol = 512eps(T)))
        end
        # Endpoint Ohm field at a node matches the independent -u1 x B1.
        for c = 1:3
            @test abs(st.fields.E[c][1] - Float64(E1[c])) <= tol
        end
        d = implicit_hybrid_diagnostics(st, species)
        @test d.phase === :synchronous_integer && d.last_solve === st.last_report
        @test d.step == st.step[] == 1 && d.time == st.time[] == h
        before = snapshot(st, species)
        step!(st, species, zero(T))
        @test snapshot(st, species) == before
        @test_throws MethodError recommended_dt(st)
    end
end

@testset "Failure leaves public state and prior report intact" begin
    for T in (Float32, Float64)
        st, species = fixture(T)
        fill!(st.fields.B[3], T(NaN))
        @test_throws ArgumentError init!(st, species)
        fill!(st.fields.B[3], T(1 / 2))
        init!(st, species)
        step!(st, species, T(1 / 64))
        for bad_dt in (T(-1), T(NaN), T(Inf))
            before = snapshot(st, species)
            @test_throws ArgumentError step!(st, species, bad_dt)
            @test snapshot(st, species) == before
        end
        # Public-field corruption is rejected at the step boundary without
        # touching any state, including the corrupted field itself.
        fill!(st.fields.B[1], T(Inf))
        corrupted = snapshot(st, species)
        @test_throws ArgumentError step!(st, species, T(1 / 64))
        @test snapshot(st, species) == corrupted
        fill!(st.fields.B[1], zero(T))
        step!(st, species, T(1 / 64))
        @test st.step[] == 2
        # A nonzero residual with max_newton = 0 exhausts the bounded solver:
        # the attempt fails as ConvergenceError with state and (empty) report
        # intact.
        st0, species0 = fixture(T)
        st0 = ImplicitHybridStepper(st0.g, st0.model, CIC(), species0; max_newton = 0)
        fill!(st0.fields.B[3], T(1 / 2))
        init!(st0, species0)
        before0 = snapshot(st0, species0)
        @test_throws ImplicitHybridConvergenceError step!(st0, species0, T(1 / 64))
        @test snapshot(st0, species0) == before0
        @test st0.last_report === nothing && iszero(st0.step[])
    end
end

@testset "Checkpoint/archive round-trip continues bitwise" begin
    for T in (Float32, Float64)
        st, species = fixture(T)
        fill!(st.fields.B[3], T(1 / 4))
        init!(st, species)
        step!(st, species, T(1 / 64))
        mktempdir() do dir
            path = joinpath(dir, "ih.chk")
            save_checkpoint(path, st, species)
            dest_species = [ParticleSet{1,T}(nparticles(ps); q = ps.q, m = ps.m) for ps in species]
            for (dst, src) in zip(dest_species, species)
                dst.x[1] .= src.x[1]
                fill!(dst.weight, one(T))
            end
            restored = ImplicitHybridStepper(
                FourierGrid(st.g.n, st.g.L),
                HybridModel(IsothermalElectrons(st.model.closure.Te); nfloor = st.model.nfloor),
                CIC(),
                dest_species,
            )
            load_checkpoint!(restored, dest_species, path)
            @test snapshot(restored, dest_species) == snapshot(st, species)
            for h in (T(1 / 128), zero(T))
                step!(st, species, h)
                step!(restored, dest_species, h)
                @test snapshot(restored, dest_species) == snapshot(st, species)
            end
        end
    end
end
end # module
