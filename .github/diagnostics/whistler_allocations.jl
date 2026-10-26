using HybridPlasmaPIC, InteractiveUtils, Profile, SHA, Test, TOML

length(ARGS) == 2 || error("Expected source checkout and evidence directory paths")
const source = realpath(ARGS[1])
const evidence = abspath(ARGS[2])
realpath(pkgdir(HybridPlasmaPIC)) == source || error("Wrong HybridPlasmaPIC source")
mkpath(evidence)
const source_sha = readchomp(Cmd(["git", "-C", source, "rev-parse", "HEAD"]))
source_sha == ENV["HPIC_DIAGNOSTIC_SOURCE_SHA"] || error("Wrong source revision")
const source_files = (
    "Project.toml",
    "src/integrators/semi_implicit.jl",
    "src/utils/exact_arithmetic.jl",
    "test/test_exact_arithmetic.jl",
    "test/test_linear_cn_range.jl",
    "test/test_whistler_packet_allocation.jl",
)
const provider = TOML.parsefile(joinpath(pkgdir(HybridPlasmaPIC.SpectralOperators), "Project.toml"))
provider["version"] == "0.1.6" || error("Wrong SpectralOperators version")
open(joinpath(evidence, "identity.toml"), "w") do io
    TOML.print(
        io,
        Dict(
            "source_sha" => source_sha,
            "runtime" => string(VERSION),
            "machine" => Sys.MACHINE,
            "cpu" => Sys.CPU_NAME,
            "threads" => Threads.nthreads(),
            "source_hashes" =>
                Dict(f => bytes2hex(sha256(read(joinpath(source, f)))) for f in source_files),
            "provider_project" => provider,
            "provider_root" => pkgdir(HybridPlasmaPIC.SpectralOperators),
            "julia_options" => sprint(show, Base.JLOptions()),
        ),
    )
end
cp(Base.active_project(), joinpath(evidence, "Project.toml"); force = false)
cp(
    joinpath(dirname(Base.active_project()), "Manifest.toml"),
    joinpath(evidence, "Manifest.toml");
    force = false,
)

@testset "Unmodified CI prefix and allocation diagnosis" begin
    # Retain the real CI prefix, all assertions and their evaluation order.
    include(joinpath(source, "test/test_exact_arithmetic.jl"))
    include(joinpath(source, "test/test_linear_cn_range.jl"))
    include(joinpath(source, "test/test_whistler_packet_allocation.jl"))
    bytes = TestWhistlerPacketAllocation.packet_bytes
    rows = []
    for steps in (0, 1, 100), n in (8, 128)
        count = BigInt(steps)
        for _ = 1:2
            bytes(n, count, 1e-6, 1.0, 1.0, Float64(2π))
            bytes(n, steps, 1e-6, 1.0, 1.0, Float64(2π))
        end
        samples = []
        for sample = 1:3
            big_bytes = bytes(n, count, 1e-6, 1.0, 1.0, Float64(2π))
            fixed_bytes = bytes(n, steps, 1e-6, 1.0, 1.0, Float64(2π))
            push!(
                samples,
                Dict(
                    "sample" => sample,
                    "big_bytes" => big_bytes,
                    "fixed_bytes" => fixed_bytes,
                    "extra_bytes" => big_bytes - fixed_bytes,
                ),
            )
        end
        push!(rows, Dict("steps" => steps, "modes" => n, "samples" => samples))
        println("RAW_PACKET_PROFILE ", last(rows))
        flush(stdout)
    end
    open(joinpath(evidence, "samples.toml"), "w") do io
        TOML.print(
            io,
            Dict(
                "profiles" => rows,
                "scope" => "All three declared samples; no minimum or retry selection",
            ),
        )
    end

    # Profiling occurs after the original gate and raw samples. Preserve every
    # allocation and frame; profiler-instrumented bytes are separate evidence.
    for steps in (0, 1, 100), n in (8, 128), big_counter in (false, true)
        count = big_counter ? BigInt(steps) : steps
        Profile.Allocs.clear()
        Profile.Allocs.@profile sample_rate = 1.0 bytes(n, count, 1e-6, 1.0, 1.0, Float64(2π))
        allocations = Profile.Allocs.fetch().allocs
        records = [
            Dict(
                "type" => string(a.type),
                "size" => a.size,
                "frames" => [
                    Dict(
                        "function" => string(f.func),
                        "file" => string(f.file),
                        "line" => f.line,
                        "inlined" => f.inlined,
                        "from_c" => f.from_c,
                    ) for f in a.stacktrace
                ],
            ) for a in allocations
        ]
        open(joinpath(evidence, "trace-$(steps)-$(n)-$(big_counter).toml"), "w") do io
            TOML.print(
                io,
                Dict(
                    "steps" => steps,
                    "modes" => n,
                    "big_counter" => big_counter,
                    "allocations" => records,
                ),
            )
        end
    end
    Profile.Allocs.clear()
    open(joinpath(evidence, "typed.txt"), "w") do io
        for I in (Int, BigInt)
            code_warntype(
                io,
                bytes,
                Tuple{Int,I,Float64,Float64,Float64,Float64};
                debuginfo = :source,
            )
            show(
                io,
                MIME"text/plain"(),
                code_typed(bytes, Tuple{Int,I,Float64,Float64,Float64,Float64}; optimize = true),
            )
            println(io)
        end
    end
end
