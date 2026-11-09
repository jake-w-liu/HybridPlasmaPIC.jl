module ImplicitHybridPressureAcceptanceTests
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

record_values(r) =
    NamedTuple{fieldnames(typeof(r))}(ntuple(i -> getfield(r, i), fieldcount(typeof(r))))
report_values(::Nothing) = nothing
report_values(r::ImplicitHybridSolveReport) = record_values(r)
report_values(r::NamedTuple) =
    (solver = record_values(r.solver), pressure = record_values(r.pressure))

function snapshot(st, species)
    return deepcopy((
        fields = H._hybrid_field_state(st.fields),
        floor = st.fields.floor_count[],
        populations = [
            (x = ps.x, v = ps.v, weight = ps.weight, id = ps.id, tag = ps.tag, q = ps.q, m = ps.m)
            for ps in species
        ],
        counts = st._population_counts,
        phase = (st.initialized[], st.time[], st.step[], st.last_dt[]),
        report = report_values(st.last_report),
    ))
end

function fixture(::Type{T}, ::Val{D}; cold = false, gamma = T(2), kwargs...) where {T,D}
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
    species = [
        ParticleSet{D,T}(prod(g.n); q = one(T), m = T(2)),
        ParticleSet{D,T}(prod(g.n); q = T(2), m = T(5)),
    ]
    sine, cosine = T[0, 1, 0, -1], T[1, 0, -1, 0]
    for (s, ps) in enumerate(species), (p, I) in enumerate(CartesianIndices(g.n))
        for d = 1:D
            ps.x[d][p] = T(I[d] - 1)
        end
        ps.weight[p] = s == 1 ? T(1 / 2) : T(1 / 4)
        for c = 1:3
            ps.v[c][p] = (s == 1 ? T(1 / 16) : -T(1 / 32)) * sine[I[mod1(c, D)]]
        end
    end
    st = ImplicitHybridStepper(g, HybridModel(EvolvingElectrons(gamma)), CIC(), species; kwargs...)
    for I in CartesianIndices(g.n)
        st.fields.pe[I] =
            cold ? zero(T) : one(T) + sum((cosine[I[d]] + sine[I[d]] / 2) / T(32d) for d = 1:D)
        st.fields.B[1][I] = D == 1 ? T(1 / 4) : sine[I[2]] / 16
        st.fields.B[2][I] = D == 3 ? sine[I[3]] / 16 : cosine[I[1]] / 16
        st.fields.B[3][I] = T(1 / 2) + sine[I[1]] / 16
    end
    return st, species
end

# On four periodic nodes, the odd Fourier derivative (zero Nyquist) is exactly
# k1*(a_next-a_previous)/2. This reference uses stored k1 and exact Rational
# rows, with no production FFT, pressure-rate or derivative helper.
function derivative4(a, I::CartesianIndex{D}, d, g) where {D}
    next = CartesianIndex(ntuple(c -> c == d ? mod1(I[c] + 1, 4) : I[c], D))
    previous = CartesianIndex(ntuple(c -> c == d ? mod1(I[c] - 1, 4) : I[c], D))
    return Q(g.kfull[d][2]) * (Q(a[next]) - Q(a[previous])) / 2
end

function pressure_reference(st::ImplicitHybridStepper{D,T}) where {D,T}
    w, g = st.work, st.g
    reference = zeros(Q, g.n)
    for I in CartesianIndices(g.n)
        advection =
            sum(Q(w.electron_velocity[d][I]) * derivative4(w.midpoint.pe, I, d, g) for d = 1:D)
        divergence = sum(derivative4(w.electron_velocity[d], I, d, g) for d = 1:D)
        reference[I] = -advection - Q(T(st.model.closure.γ)) * Q(w.midpoint.pe[I]) * divergence
    end
    return reference
end

@testset "Joint evolving pressure / independent four-node Fourier rows" begin
    for (T, D) in ((Float32, 1), (Float64, 1), (Float64, 2), (Float64, 3))
        st, species = fixture(T, Val(D))
        init!(st, species)
        old_pe, old_B = copy(st.fields.pe), map(copy, st.fields.B)
        old_v = [map(copy, ps.v) for ps in species]
        arrays = H._implicit_hybrid_field_arrays(st.fields)
        owners = copy(species)
        h = T(1 / 64)
        step!(st, species, h)
        r, w = st.last_report, st.work
        @test r.solver.initial_residual > r.solver.tolerance
        @test r.solver.velocity_residual <= r.solver.tolerance
        @test r.solver.magnetic_residual <= r.solver.tolerance
        @test r.pressure.exact_residual <= Q(r.solver.tolerance)
        @test r.pressure.residual <= r.solver.tolerance
        @test !isequal(old_pe, st.fields.pe)
        @test any(c -> !isequal(old_B[c], st.fields.B[c]), 1:3)
        @test any(s -> !isequal(old_v[s], species[s].v), eachindex(species))
        @test all(x -> x >= zero(T), st.fields.pe)
        @test all(
            i -> arrays[i] === H._implicit_hybrid_field_arrays(st.fields)[i],
            eachindex(arrays),
        )
        @test all(s -> species[s] === owners[s], eachindex(species))
        reference = pressure_reference(st)
        # Engineering budget for the small represented FFT vs exact four-row
        # matrix; the nonlinear residual bound below is a separate exact gate.
        budget = Q(256) * Q(eps(T)) * (1 + maximum(abs, reference))
        @test maximum(abs, Q.(w.pressure_rate) - reference) <= budget
        independent_row =
            (Q.(st.fields.pe) - Q.(w.old.pe) - Q(h) * reference) ./ Q(w.pressure_scale)
        @test maximum(abs, independent_row) <=
              Q(r.solver.tolerance) + Q(h) * budget / Q(w.pressure_scale)
        for d = 1:D
            @test any(
                I -> !iszero(derivative4(w.midpoint.pe, I, d, st.g)),
                CartesianIndices(st.g.n),
            )
        end
        p = r.pressure
        @test p.delta_electron - p.pressure_rate_work == p.pressure_residual_energy
        @test p.energy_change ==
              p.delta_kinetic + p.delta_magnetic + p.delta_electron ==
              sum(p.defects)
        diagnostics = implicit_hybrid_diagnostics(st, species)
        @test diagnostics.phase === :synchronous_integer
        @test isequal(
            diagnostics.energy.electron_internal,
            electron_internal_energy(st.fields.n, st.model.closure, st.g; pe = st.fields.pe),
        )
        before, report = snapshot(st, species), st.last_report
        step!(st, species, zero(T))
        @test isequal(snapshot(st, species), before) && st.last_report === report
    end
end

@testset "Both directions of the coupled pressure dependency" begin
    st, species = fixture(Float64, Val(2))
    init!(st, species)
    before = snapshot(st, species)
    H._implicit_hybrid_stage!(st, species, 1 / 64)
    w = st.work
    z, out = copy(w.z), similar(w.z)
    nv, cells = 3sum(nparticles, species), length(st.fields.n)
    H._implicit_hybrid_residual!(out, st, z)
    base_velocity, base_E, base_pressure =
        copy(out[1:nv]), map(copy, w.midpoint.E), copy(w.pressure_rate)
    trial = copy(z)
    trial[nv+3cells+1] += 1 / 16
    H._implicit_hybrid_residual!(out, st, trial)
    @test !isequal(out[1:nv], base_velocity)
    @test !isequal(w.midpoint.E, base_E)
    trial = copy(z)
    trial[1] += 1 / 16
    H._implicit_hybrid_residual!(out, st, trial)
    @test !isequal(w.pressure_rate, base_pressure)
    trial = copy(z)
    trial[nv+2cells+1] += 1 / 16
    H._implicit_hybrid_residual!(out, st, trial)
    @test !isequal(w.pressure_rate, base_pressure)
    trial[nv+3cells+1] = -1
    @test_throws H._ImplicitHybridTrialError H._implicit_hybrid_residual!(out, st, trial)
    @test isequal(snapshot(st, species), before)
end

@testset "Cold pressure and rejected operations preserve the whole state" begin
    for T in (Float32, Float64)
        cold, particles = fixture(T, Val(1); cold = true)
        init!(cold, particles)
        step!(cold, particles, T(1 / 64))
        @test all(iszero, cold.fields.pe) && iszero(cold.last_report.pressure.exact_residual)
        for invalid in (T(-1), T(NaN), T(Inf))
            st, species = fixture(T, Val(1))
            st.fields.pe[1] = invalid
            before = snapshot(st, species)
            @test_throws ArgumentError init!(st, species)
            @test isequal(snapshot(st, species), before)
        end
        st, species = fixture(T, Val(1); max_newton = 0)
        # B=0 and uniform pe make the velocity/Faraday rows exactly zero at
        # the initial trial. Compression leaves only the pressure row active.
        foreach(a -> fill!(a, zero(T)), st.fields.B)
        fill!(st.fields.pe, one(T))
        init!(st, species)
        before = snapshot(st, species)
        failure = try
            step!(st, species, T(1 / 64))
            nothing
        catch err
            err
        end
        @test failure isa ImplicitHybridConvergenceError
        @test failure.reason === :nonlinear_limit
        @test iszero(failure.velocity_residual) && iszero(failure.magnetic_residual)
        @test failure.pressure_residual > failure.tolerance
        @test isequal(snapshot(st, species), before)
        st.fields.pe[1] = -one(T)
        before = snapshot(st, species)
        @test_throws ArgumentError step!(st, species, T(1 / 64))
        @test isequal(snapshot(st, species), before)
        empty_density, populations = fixture(T, Val(1))
        foreach(ps -> fill!(ps.weight, zero(T)), populations)
        before = snapshot(empty_density, populations)
        @test_throws ImplicitHybridConvergenceError init!(empty_density, populations)
        @test isequal(snapshot(empty_density, populations), before)
    end
    grid = FourierGrid((4,), (Float32(4),))
    particles = [ParticleSet{1,Float32}(0)]
    @test_throws ArgumentError ImplicitHybridStepper(
        grid,
        HybridModel(EvolvingElectrons(nextfloat(1.0))),
        CIC(),
        particles,
    )
end

@testset "Schema2 complete checkpoint/archive and variable-step continuation" begin
    for T in (Float32, Float64)
        # Float64 exercises an original MPFR-backed coefficient through both
        # formats; the evolved equation still uses the ordinary grid gamma = 2.
        gamma, other_gamma = setprecision(BigFloat, 128) do
            value = BigFloat(2)
            value, nextfloat(value)
        end
        st, species = fixture(T, Val(1); gamma = T === Float64 ? gamma : T(2))
        init!(st, species)
        step!(st, species, T(1 / 64))
        step!(st, species, T(1 / 128))
        mktempdir() do directory
            checkpoint, archive =
                joinpath(directory, "pressure.chk"), joinpath(directory, "pressure.run")
            save_checkpoint(checkpoint, st, species)
            saved = deserialize(checkpoint)
            @test saved.schema === 2 && saved.integrator === :implicit_hybrid_midpoint
            if T === Float64
                @test typeof(st.model.closure.γ) === BigFloat
                @test saved.configuration.model.gamma_type === BigFloat
                @test saved.configuration.model.gamma_numerator == "2"
                @test saved.configuration.model.gamma_denominator == "1"
                @test T(other_gamma) == T(gamma) && Q(other_gamma) != Q(gamma)
            end
            dest_species = [ParticleSet{1,T}(0; q = ps.q, m = ps.m) for ps in species]
            restored = ImplicitHybridStepper(
                FourierGrid(st.g.n, st.g.L),
                deepcopy(st.model),
                CIC(),
                dest_species,
            )
            for a in H._implicit_hybrid_field_arrays(restored.fields)
                fill!(a, T(NaN))
            end
            restored.time[], restored.step[], restored.last_dt[] = T(NaN), -1, -one(T)
            restored.fields.floor_count[] = -1
            load_checkpoint!(restored, dest_species, checkpoint)
            @test isequal(snapshot(restored, dest_species), snapshot(st, species))
            for h in (T(3 / 256), zero(T))
                step!(st, species, h)
                step!(restored, dest_species, h)
                @test isequal(snapshot(restored, dest_species), snapshot(st, species))
            end
            archive_run(archive, st, species; rng_seed = 0x517)
            step!(restored, dest_species, T(1 / 128))
            load_archive!(restored, dest_species, archive)
            @test isequal(snapshot(restored, dest_species), snapshot(st, species))
            if T === Float64
                mismatch_species = [ParticleSet{1,T}(0; q = ps.q, m = ps.m) for ps in species]
                mismatch = ImplicitHybridStepper(
                    FourierGrid(st.g.n, st.g.L),
                    HybridModel(EvolvingElectrons(other_gamma)),
                    CIC(),
                    mismatch_species,
                )
                before = snapshot(mismatch, mismatch_species)
                @test_throws ArgumentError load_checkpoint!(mismatch, mismatch_species, checkpoint)
                @test isequal(snapshot(mismatch, mismatch_species), before)
                @test_throws ArgumentError load_archive!(mismatch, mismatch_species, archive)
                @test isequal(snapshot(mismatch, mismatch_species), before)
            end
            for h in (T(1 / 256),)
                step!(st, species, h)
                step!(restored, dest_species, h)
                @test isequal(snapshot(restored, dest_species), snapshot(st, species))
            end
            bad_pressure = deepcopy(saved)
            bad_pressure.fields.pe[1] = -one(T)
            values = record_values(saved.last_report.pressure)
            bad_receipt = typeof(saved.last_report.pressure)(
                Base.values(merge(values, (energy_change = values.energy_change + Q(1),)))...,
            )
            bad_report = merge(saved.last_report, (pressure = bad_receipt,))
            for bad in (
                merge(saved, (schema = 1,)),
                bad_pressure,
                merge(saved, (last_report = bad_report,)),
            )
                serialize(checkpoint, bad)
                before, report = snapshot(restored, dest_species), restored.last_report
                @test_throws ArgumentError load_checkpoint!(restored, dest_species, checkpoint)
                @test isequal(snapshot(restored, dest_species), before) &&
                      restored.last_report === report
            end
        end
    end
end
end # module
