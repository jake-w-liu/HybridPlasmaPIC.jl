module CollisionalCGLTests
using HybridPlasmaPIC, Test, Random, Serialization
const H = HybridPlasmaPIC
include("collisional_fixtures.jl")

function fixture(
    ::Type{T};
    D = 1,
    shape = CIC(),
    logs = T(0.05),
    R = Xoshiro,
    closure = CGLElectrons(T(6e-4), T(3e-4), T(3), T(1 / 4)),
) where {T}
    st, ps = _composition_fixture(Val(D), T, R; shape, logs, closure, magnetic = T(0.2))
    for I in CartesianIndices(st.grid.n)
        angle = T(2pi) * T(I[1] - 1) / T(st.grid.n[1])
        st.inner.fields.B[2][I] = T(0.03) * sin(angle)
        st.inner.fields.B[3][I] = T(0.15) + T(0.02) * cos(angle)
    end
    return st, ps
end
function fresh_density(st, species)
    n, part = zeros(eltype(st.inner.fields.n), st.grid.n), similar(st.inner.fields.n)
    for p in species
        density!(part, p, st.grid, st.inner.shape)
        n .+= p.q .* part
    end
    return n
end

@testset "CGL charged composition, full anisotropic phase and source invariants" begin
    for (T, D, shape) in
        ((Float32, 1, NGP()), (Float64, 1, CIC()), (Float32, 2, TSC()), (Float64, 3, CIC())),
        R in (Xoshiro, MersenneTwister)

        st, ps = fixture(T; D, shape, R)
        reference, rps = fixture(T; D, shape, R)
        @test st.transaction.Eforce === nothing
        @test all(
            a -> size(a) == st.grid.n,
            (st.transaction.fields.pforce..., st.transaction.fields.cgl_dp),
        )
        init!(st, ps)
        init!(reference, rps)
        for _ = 1:2
            step!(st, ps, st.dt)
            _composition_component!(reference, rps)
            @test isequal(_composition_physical(st, ps), _composition_physical(reference, rps))
        end
        # Inspect the real ion-only subflow on the already advanced reference;
        # its charge and constitutive energy use freshly deposited positions.
        positions = [map(copy, p.x) for p in rps]
        B, n = map(copy, reference.inner.fields.B), fresh_density(reference, rps)
        energy = electron_internal_energy(n, B, reference.inner.model.closure, reference.grid)
        collide_coulomb!(
            rps,
            reference.grid,
            reference.dt;
            coulomb_log = reference.coulomb_log,
            epsilon0 = reference.epsilon0,
            weight_to_number = reference.weight_to_number,
            u_floor = reference.u_floor,
            rng = reference.rng,
            workspace = reference.workspace,
        )
        @test [p.x for p in rps] == positions
        @test reference.inner.fields.B == B && fresh_density(reference, rps) == n
        @test electron_internal_energy(
            n,
            reference.inner.fields.B,
            reference.inner.model.closure,
            reference.grid,
        ) == energy
        before = _composition_snapshot(st, ps)
        step!(st, ps, zero(T))
        @test isequal(before, _composition_snapshot(st, ps))
        stream = _composition_rng_state(st.rng)
        init!(st, ps)
        @test isequal(stream, _composition_rng_state(st.rng)) && st.source_count == 0
    end
    # Closure precision belongs to its original coefficients, independently
    # of the field type. Zero invariants still allocate all anisotropic roles.
    for closure in (
        CGLElectrons(6e-4, 3e-4, 3.0, 0.25),
        CGLElectrons(Float16(6e-4), Float16(3e-4), Float16(3), Float16(1 / 4)),
        CGLElectrons(Float32(0), Float32(0), Float32(1), Float32(1)),
    )
        st, ps = fixture(Float32; closure, logs = 0.0f0)
        ref, rps = fixture(Float32; closure, logs = 0.0f0)
        init!(st, ps)
        init!(ref, rps)
        stream = _composition_rng_state(st.rng)
        step!(st, ps, st.dt)
        step!(ref.inner, rps, ref.dt)
        @test isequal(_composition_physical(st, ps), _composition_physical(ref, rps))
        @test isequal(stream, _composition_rng_state(st.rng))
    end
    empty_st, ep = _composition_fixture(
        Val(1),
        Float64,
        Xoshiro;
        count = 0,
        logs = 0.0,
        closure = CGLElectrons(0.0, 0.0, 1.0, 1.0),
    )
    init!(empty_st, ep)
    stream = _composition_rng_state(empty_st.rng)
    step!(empty_st, ep, empty_st.dt)
    @test isempty(empty_st.transaction.x[1]) && empty_st.source_count == 1
    @test isequal(stream, _composition_rng_state(empty_st.rng))
end

@testset "CGL admission, anisotropic rollback roles and real field rejection" begin
    st, ps = fixture(Float64)
    function wrap(inner; electron_ion = nothing)
        CollisionalHybridStepper(
            inner,
            ps;
            dt = st.dt,
            rng = Xoshiro(1),
            coulomb_log = 0.05,
            epsilon0 = 1.0,
            weight_to_number = 1.0,
            electron_ion,
        )
    end
    request = ElectronIonThermalExchange(; electron_mass = 1 / 16, coulomb_log = 1.0)
    @test_throws ArgumentError wrap(st.inner; electron_ion = request)
    for c in (
        CGLElectrons{Float64}(-1.0, 0.1),
        CGLElectrons{Float64}(0.1, Inf),
        CGLElectrons{Float64}(NaN, 0.1),
        CGLElectrons{Union{Float32,Float64}}(0.1f0, 0.1),
    )
        @test_throws ArgumentError wrap(
            HybridStepper(st.grid, HybridModel(c), CIC(), sum(nparticles, ps)),
        )
    end
    for name in (:η, :ηH, :de2)
        model =
            name === :η ? HybridModel(st.model.closure; η = eps()) :
            name === :ηH ? HybridModel(st.model.closure; ηH = eps()) :
            HybridModel(st.model.closure; de2 = eps())
        @test_throws ArgumentError wrap(HybridStepper(st.grid, model, CIC(), sum(nparticles, ps)))
    end
    init!(st, ps)
    saved = st.transaction.fields
    original = saved.pforce
    saved.pforce = (st.inner.fields.pforce[1], original[2], original[3])
    entry = _composition_snapshot(st, ps)
    @test_throws ArgumentError step!(st, ps, st.dt)
    @test isequal(entry, _composition_snapshot(st, ps))
    saved.pforce = original
    # Uniform n and B give a valid initial state. Raising B makes actual
    # pperp=n*B unrepresentable; the canonical CGL force must reject it.
    bad, bp = _composition_fixture(
        Val(1),
        Float64,
        Xoshiro;
        shape = NGP(),
        closure = CGLElectrons(1.0, 1.0, 1.0, 1.0),
    )
    for p in bp, v in p.v
        fill!(v, 0.0)
    end
    init!(bad, bp)
    fill!(bad.inner.fields.B[1], floatmax(Float64))
    entry = _composition_snapshot(bad, bp)
    failure = try
        step!(bad, bp, bad.dt)
        nothing
    catch err
        err
    end
    @test failure isa ArgumentError
    @test occursin("CGL pressure exceeds numeric range", sprint(showerror, failure))
    @test isequal(entry, _composition_snapshot(bad, bp)) && !bad.busy[]
end

@testset "CGL source failure restores advanced fields and RNG for exact retry" begin
    rejected = 0
    for seed = 1:8
        closure = CGLElectrons(0.0, 0.0, 1.0, 1.0)
        st, ps = _composition_sampled_failure_fixture(Float64, Xoshiro, seed; closure)
        init!(st, ps)
        before = _composition_snapshot(st, ps)
        failure = try
            step!(st, ps, st.dt)
            nothing
        catch err
            err
        end
        failure === nothing && continue
        rejected += 1
        @test failure isa ArgumentError
        @test isequal(before, _composition_snapshot(st, ps))
        retry = try
            step!(st, ps, st.dt)
            nothing
        catch err
            err
        end
        @test typeof(retry) === typeof(failure) &&
              sprint(showerror, retry) == sprint(showerror, failure)
        @test isequal(before, _composition_snapshot(st, ps))
        direct, dp = _composition_sampled_failure_fixture(Float64, Xoshiro, seed; closure)
        init!(direct, dp)
        stream = _composition_rng_state(direct.rng)
        cause = try
            _composition_component!(direct, dp)
            nothing
        catch err
            err
        end
        @test typeof(cause) === typeof(failure)
        @test direct.inner.step[] == 1 && !isequal(stream, _composition_rng_state(direct.rng))
    end
    @test rejected > 0
end

@testset "CGL schema4 raw/archive carries every anisotropic field" begin
    for R in (Xoshiro, MersenneTwister)
        st, ps = fixture(Float64; R)
        init!(st, ps)
        step!(st, ps, st.dt)
        payload = deepcopy(H._checkpoint_state(st, ps))
        @test payload.schema == 4 && payload.cgl_policy_version === 1
        @test payload.cgl_force_floor === Float64(1e-12)
        @test payload.cgl_force_operator === :global_invariants_trial_B_tensor
        mktempdir() do dir
            raw, archive = joinpath(dir, "cgl.chk"), joinpath(dir, "cgl.run")
            save_checkpoint(raw, st, ps)
            archive_run(archive, st, ps; rng_seed = 712)
            target, tp = fixture(Float64; R)
            init!(target, tp)
            step!(target, tp, target.dt)
            step!(target, tp, target.dt)
            ids = map(
                objectid,
                H._collisional_field_arrays(values(H._hybrid_field_state(target.inner.fields))),
            )
            load_checkpoint!(target, tp, raw)
            @test isequal(_composition_snapshot(target, tp), _composition_snapshot(st, ps))
            @test map(
                objectid,
                H._collisional_field_arrays(values(H._hybrid_field_state(target.inner.fields))),
            ) == ids
            step!(st, ps, st.dt)
            step!(target, tp, target.dt)
            @test isequal(_composition_snapshot(target, tp), _composition_snapshot(st, ps))
            load_archive!(target, tp, archive)
            for damage in (:policy, :floor, :operator, :schema, :model, :layout)
                incoming = deepcopy(payload)
                incoming =
                    damage === :policy ? merge(incoming, (cgl_policy_version = true,)) :
                    damage === :floor ? merge(incoming, (cgl_force_floor = 2e-12,)) :
                    damage === :operator ? merge(incoming, (cgl_force_operator = :scalar,)) :
                    damage === :schema ? merge(incoming, (schema = 1,)) :
                    damage === :model ?
                    merge(
                        incoming,
                        (
                            inner = merge(
                                incoming.inner,
                                (
                                    model = HybridModel(
                                        CGLElectrons{Float64}(
                                            2st.model.closure.Cperp,
                                            st.model.closure.Cpar,
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ) : merge(incoming, (fields = merge(incoming.fields, (cgl_dp = zeros(1),)),))
                before = _composition_snapshot(target, tp)
                @test_throws ArgumentError H._restore_checkpoint_state!(target, tp, incoming)
                @test isequal(before, _composition_snapshot(target, tp))
            end
        end
    end
end

@testset "CGL original wide coefficient and fixed-context continuation" begin
    setprecision(BigFloat, 128) do
        # Cpar underflows in Float64, while Cpar*n^3/B^2 = 2^-920 is finite.
        closure = CGLElectrons{BigFloat}(BigFloat(0), BigFloat(2)^(-4000))
        function make()
            st, ps = _composition_fixture(
                Val(1),
                Float64,
                Xoshiro;
                ns = 1,
                shape = NGP(),
                logs = 0.0,
                closure,
                magnetic = 2.0^(-40),
            )
            for p in ps
                fill!(p.weight, 2.0^998)
                for i in eachindex(p.weight)
                    p.x[1][i] = (mod(i - 1, 4) + 0.25) * 0.5
                end
                foreach(v -> fill!(v, 0.0), p.v)
            end
            init!(st, ps)
            return st, ps
        end
        st, ps = make()
        @test Float64(st.model.closure.Cpar) == 0 && st.model.closure.Cpar > 0
        n = fresh_density(st, ps)
        @test all(==(2.0^1000), n)
        @test electron_internal_energy(n, st.inner.fields.B, st.model.closure, st.grid) ==
              BigFloat(2)^(-920)
        step!(st, ps, st.dt)
        mktempdir() do dir
            raw, archive = joinpath(dir, "wide.chk"), joinpath(dir, "wide.run")
            save_checkpoint(raw, st, ps)
            archive_run(archive, st, ps; rng_seed = 712)
            for (path, loader) in ((raw, load_checkpoint!), (archive, load_archive!))
                target, tp = make()
                loader(target, tp, path)
                @test isequal(_composition_snapshot(target, tp), _composition_snapshot(st, ps))
                step!(target, tp, target.dt)
                reference, rp = make()
                H._restore_checkpoint_state!(reference, rp, H._checkpoint_state(st, ps))
                @test reference.inner.time[] == st.inner.time[]
                @test reference.inner.step[] == st.inner.step[]
                @test reference.dt == st.dt
                # Continue through the public canonical loader as well. This is
                # the supported restart path and keeps the target/reference
                # comparison on identical restored state.
                loader(reference, rp, path)
                step!(reference, rp, reference.dt)
                @test isequal(
                    _composition_snapshot(target, tp),
                    _composition_snapshot(reference, rp),
                )
            end
        end
    end
end
end # module
