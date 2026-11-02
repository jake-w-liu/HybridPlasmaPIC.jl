# Remote-only auxiliary observation. Run in a fresh Julia process AFTER the
# unchanged strict suite; this does not replace or accept any test/JET gate.
using HybridPlasmaPIC
using InteractiveUtils
using SHA

const H = HybridPlasmaPIC
const OBSERVATION_OUTPUT = get(ENV, "YEE_CONTEXT_OUTPUT", "yee-stored-zero-context")
mkpath(OBSERVATION_OUTPUT)

# These reproduce the indexed expression without hiding or changing inference.
# Their lowering is auxiliary; they are not the unchanged Test macro caller.
indexed_shared(a::Vector{T}) where {T} = H._particle_stored_zero(a[1])
indexed_delegate(a::Vector{T}) where {T} = H._yl_stored_zero(a[1])

versioninfo()
println((
    purpose = :stored_zero_context_observation,
    acceptance = false,
    julia = VERSION,
    build_commit = Base.GIT_VERSION_INFO.commit,
    arch = Sys.ARCH,
    cpu = Sys.CPU_NAME,
    threads = Threads.nthreads(),
    julia_cpu_target = get(ENV, "JULIA_CPU_TARGET", "<unset>"),
))
for source in (
    "src/utils/validation.jl",
    "src/models/yee_longitudinal_pic.jl",
    "test/test_yee_longitudinal_pic.jl",
    "test/test_yee_particle_charge.jl",
)
    println((source = source, sha256 = bytes2hex(sha256(read(source)))))
end

# Retain the actual runner's Base source alongside generated operations; local
# source inspection does not establish the exact deployed runner build.
base_output = joinpath(OBSERVATION_OUTPUT, "base-source")
mkpath(base_output)
for name in (
    "essentials.jl",
    "int.jl",
    "promotion.jl",
    "rounding.jl",
    "reinterpretarray.jl",
    "subarray.jl",
    "float.jl",
)
    source = Base.find_source_file(name)
    source === nothing && error("runner Base source unavailable: $name")
    destination = joinpath(base_output, name)
    cp(source, destination)
    println((
        runner_base_source = source,
        retained = destination,
        sha256 = bytes2hex(sha256(read(destination))),
    ))
end

# The same four precision/context rows as the unchanged longitudinal test.
# Construct and retain the word view before changing the arithmetic context.
# No floating display or logging is performed inside the changed context.
for T in (Float32, Float64)
    U = T === Float32 ? UInt32 : UInt64
    subnormal = T[reinterpret(T, one(U))]
    words = reinterpret(U, subnormal)
    for (roundmode, flush) in ((RoundNearest, true), (RoundDown, false))
        saved = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
        record = nothing
        try
            requested_fenv = Base.Rounding.to_fenv(roundmode)
            rounding_setter = Base.Rounding.setrounding_raw(Float64, requested_fenv)
            flush_setter = Base.set_zero_subnormals(flush)
            active_before = (
                fenv_f32 = Base.Rounding.rounding_raw(Float32),
                fenv_f64 = Base.Rounding.rounding_raw(Float64),
                rounding_f32 = rounding(Float32),
                rounding_f64 = rounding(Float64),
                flush = get_zero_subnormals(),
                thread = Threads.threadid(),
            )
            memory_word_before = words[1]
            loaded = subnormal[1]
            loaded_word = reinterpret(U, loaded)
            masked_word = loaded_word & (typemax(U) >> 1)
            masked_iszero = iszero(masked_word)
            shared = H._particle_stored_zero(loaded)
            delegate = H._yl_stored_zero(loaded)
            shared_indexed = indexed_shared(subnormal)
            delegate_indexed = indexed_delegate(subnormal)
            memory_word_after = words[1]
            active_after = (
                fenv_f32 = Base.Rounding.rounding_raw(Float32),
                fenv_f64 = Base.Rounding.rounding_raw(Float64),
                rounding_f32 = rounding(Float32),
                rounding_f64 = rounding(Float64),
                flush = get_zero_subnormals(),
                thread = Threads.threadid(),
            )
            record = (
                precision = T,
                requested_rounding = roundmode,
                requested_flush = flush,
                requested_fenv = requested_fenv,
                rounding_setter_return = rounding_setter,
                flush_setter_return = flush_setter,
                active_before = active_before,
                memory_word_before = memory_word_before,
                loaded_word = loaded_word,
                masked_word = masked_word,
                masked_iszero = masked_iszero,
                shared = shared,
                delegate = delegate,
                shared_indexed = shared_indexed,
                delegate_indexed = delegate_indexed,
                memory_word_after = memory_word_after,
                active_after = active_after,
            )
        finally
            Base.Rounding.setrounding_raw(Float64, saved[1])
            Base.set_zero_subnormals(saved[2])
        end
        println(record)
        println((
            restored_fenv = Base.Rounding.rounding_raw(Float64),
            restored_flush = get_zero_subnormals(),
            saved = saved,
            current_thread = Threads.threadid(),
        ))
    end
end

# Capture complete unfiltered code for the unchanged shared helper/delegate,
# and for the auxiliary indexed callers if their behavior differs. Capture is
# after the mode rows, with the original context restored by finally above.
for T in (Float32, Float64)
    for (label, operation, signature) in (
        (:shared, H._particle_stored_zero, Tuple{T}),
        (:delegate, H._yl_stored_zero, Tuple{T}),
        (:indexed_shared, indexed_shared, Tuple{Vector{T}}),
        (:indexed_delegate, indexed_delegate, Tuple{Vector{T}}),
    )
        prefix = joinpath(OBSERVATION_OUTPUT, string(label, "-", T))
        for optimize in (false, true)
            path = string(prefix, optimize ? "-optimized.ll" : "-unoptimized.ll")
            open(path, "w") do io
                code_llvm(
                    io,
                    operation,
                    signature;
                    raw = true,
                    dump_module = true,
                    optimize = optimize,
                    debuginfo = :source,
                )
            end
            println((generated = path, sha256 = bytes2hex(sha256(read(path)))))
        end
        path = string(prefix, ".native")
        open(path, "w") do io
            code_native(
                io,
                operation,
                signature;
                raw = true,
                dump_module = true,
                binary = true,
                debuginfo = :source,
            )
        end
        println((generated = path, sha256 = bytes2hex(sha256(read(path)))))
    end
end
println((
    observation_completed = true,
    acceptance = false,
    note = "Original unmodified numerical and JET gates remain authoritative.",
))
