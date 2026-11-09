# Run write/check in separate Julia processes with the pinned old/new projects.
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC

length(ARGS) == 3 || throw(ArgumentError("expected write|check, output directory, package root"))
mode, directory, expected_root = ARGS
mode in ("write", "check") || throw(ArgumentError("unsupported compatibility operation"))
realpath(pkgdir(H)) == realpath(expected_root) ||
    error("compatibility process loaded the wrong package source")

function fixture(::Type{T}; count = 4) where {T}
    g = FourierGrid((4,), (T(4),))
    species = [ParticleSet{1,T}(count; q = 1, m = 2), ParticleSet{1,T}(count; q = 2, m = 5)]
    for (index, p) in enumerate(species)
        for i in eachindex(p.weight)
            p.x[1][i] = T(i - 1)
        end
        fill!(p.weight, T(index == 1 ? 1 / 2 : 1 / 4))
        fill!(p.v[1], T(index == 1 ? 1 / 16 : -1 / 32))
    end
    st = ImplicitHybridStepper(g, HybridModel(IsothermalElectrons(T(1 / 8))), CIC(), species)
    fill!(st.fields.B[1], one(T))
    return st, species
end

function snapshot(st, species)
    report = st.last_report
    values =
        report === nothing ? nothing :
        NamedTuple{fieldnames(typeof(report))}(
            ntuple(i -> getfield(report, i), fieldcount(typeof(report))),
        )
    return deepcopy((
        fields = H._hybrid_field_state(st.fields),
        floor = st.fields.floor_count[],
        species = [
            (x = p.x, v = p.v, weight = p.weight, id = p.id, tag = p.tag, q = p.q, m = p.m) for
            p in species
        ],
        counts = st._population_counts,
        phase = (st.initialized[], st.time[], st.step[], st.last_dt[]),
        report = values,
    ))
end

@testset "Historical advanced implicit-hybrid schema1 / $mode" begin
    for T in (Float32, Float64)
        prefix = joinpath(directory, string(T))
        if mode == "write"
            mkpath(directory)
            st, species = fixture(T)
            init!(st, species)
            step!(st, species, T(1 / 64))
            step!(st, species, T(1 / 128))
            @test st.step[] == 2 && fieldcount(typeof(st.last_report)) == 19
            save_checkpoint(prefix * ".chk", st, species)
            archive_run(prefix * ".run", st, species; rng_seed = 0)
            @test deserialize(prefix * ".chk").schema === 1
            before = snapshot(st, species)
            step!(st, species, T(3 / 256))
            serialize(prefix * ".expected", (before = before, after = snapshot(st, species)))
        else
            expected = deserialize(prefix * ".expected")
            for format in (:raw, :archive)
                st, species = fixture(T; count = 0)
                owners, fields = copy(species), H._implicit_hybrid_field_arrays(st.fields)
                foreach(a -> fill!(a, T(NaN)), fields)
                st.time[], st.step[], st.last_dt[] = T(NaN), -1, -one(T)
                if format === :raw
                    load_checkpoint!(st, species, prefix * ".chk")
                else
                    load_archive!(st, species, prefix * ".run")
                end
                @test st.last_report isa ImplicitHybridSolveReport{T}
                @test isequal(snapshot(st, species), expected.before)
                @test all(i -> species[i] === owners[i], eachindex(species))
                @test all(
                    i -> H._implicit_hybrid_field_arrays(st.fields)[i] === fields[i],
                    eachindex(fields),
                )
                step!(st, species, T(3 / 256))
                @test isequal(snapshot(st, species), expected.after)
            end
        end
        println((julia = VERSION, mode = mode, precision = T, package_source = pkgdir(H)))
    end
end
