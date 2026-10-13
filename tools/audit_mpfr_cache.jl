# Run in a fresh process: julia --startup-file=no --threads=4 tools/audit_mpfr_cache.jl
# This detects gross constant-cache corruption; a passing sample is not a
# certification of MPFR thread safety. It does not load HybridPlasmaPIC.
using Libdl

function audit_mpfr_constant(bits, which, ready, go)
    # Reach the gate before allocating the test value, so an allocation failure
    # cannot strand the coordinator waiting for this worker's readiness signal.
    put!(ready, nothing)
    take!(go)
    value = BigFloat(; precision = bits)
    if which === :pi
        ccall(
            (:mpfr_const_pi, Base.MPFR.libmpfr),
            Cint,
            (Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
            value,
            Base.MPFR.RoundNearest,
        )
    else
        ccall(
            (:mpfr_const_log2, Base.MPFR.libmpfr),
            Cint,
            (Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
            value,
            Base.MPFR.RoundNearest,
        )
    end
    # These hardware-precision references do not use MPFR's shared caches.
    expected = which === :pi ? Float64(pi) : log(2.0)
    observed = Float64(value)
    valid = isfinite(observed) && isapprox(observed, expected; rtol = 1e-12, atol = 0)
    return (; bits, observed, expected, valid)
end

function audit_mpfr_cache()
    Threads.nthreads(:default) >= 2 ||
        error("The MPFR concurrency probe requires at least two default threads; use --threads=4")
    library = Base.MPFR.libmpfr
    tls = ccall((:mpfr_buildopt_tls_p, Base.MPFR.libmpfr), Cint, ())
    shared = ccall((:mpfr_buildopt_sharedcache_p, Base.MPFR.libmpfr), Cint, ())
    patches = unsafe_string(ccall((:mpfr_get_patches, Base.MPFR.libmpfr), Cstring, ()))
    println((;
        julia = VERSION,
        mpfr = Base.MPFR.version(),
        library = Libdl.dlpath(library),
        threads = Threads.nthreads(:default),
        tls,
        shared,
        patches,
    ))
    flush(stdout)
    tls != 0 || error("MPFR was built without thread-local state; concurrent calls are unsupported")
    # Ascending requests force cache growth. Keep the workload bounded at 64
    # values per constant and 8192 bits, matching the original reproducer.
    precisions = 128:128:8192
    failures = 0
    for which in (:pi, :log2)
        ready = Channel{Nothing}(length(precisions))
        go = Channel{Nothing}(length(precisions))
        results = @sync begin
            jobs = map(precisions) do bits
                Threads.@spawn audit_mpfr_constant(bits, which, ready, go)
            end
            foreach(_ -> take!(ready), precisions)
            foreach(_ -> put!(go, nothing), precisions)
            fetch.(jobs)
        end
        incorrect = filter(result -> !result.valid, results)
        failures += length(incorrect)
        println((; constant = which, checks = length(results), incorrect))
        flush(stdout)
    end
    failures == 0 ||
        error("MPFR returned $failures incorrect constants; see docs/MPFR_CONCURRENCY_PLAN.md")
    println(
        "No incorrect constants observed in this fresh process; this is not a thread-safety certification.",
    )
    return nothing
end

audit_mpfr_cache()
