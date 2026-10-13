# Concurrent BigFloat runtime investigation

Status: a native MPFR shared-cache defect is reproduced on the installed
Julia 1.12.7 and 1.13.1 runtimes. A temporary native repair passes the focused
C reproducer. Installed Julia libraries are unchanged. The radiation arithmetic
candidate remains uncommitted while its concurrent-runtime contract is resolved.

## Reproduction and cause

VERIFIED: `/tmp/hybrid-radiation-base-concurrency.jl` reproduces an incorrect
`log(BigFloat(0.6875))` using only Julia Base. Its failing run returns about
-45.87345 rather than -0.37469345:
`/tmp/hybrid-radiation-base-concurrency112-24.log`.

Direct `mpfr_const_pi` calls at concurrently increasing precisions return
incorrect constants on both installed runtimes:
`/tmp/hybrid-mpfr-cache-concurrency112-1.log` and
`/tmp/hybrid-mpfr-cache-concurrency113-1.log`. Both load MPFR 4.2.2 with TLS
and shared caches enabled; `mpfr_get_patches()` returns an empty string.
This failure does not require loading HybridPlasmaPIC.

The [official MPFR 4.2.2 source](https://www.mpfr.org/mpfr-4.2.2/)
archive has SHA256
`b67ba0383ef7e8a8563734e2e889ef5ec3c3b898a01d00fa0a6869ad81c6ce01`.
Its `src/cache.c:107` downgrades the cache lock after updating local `cprec`.
The POSIX macro in `src/mpfr-thread.h:137` releases the write lock and then
acquires a read lock. Another writer can increase the cache precision between
those operations. The subsequent rounding then uses the old local precision
with the new mantissa storage.

VERIFIED: a standalone native C build reproduces wrong pi in 16 of 24 fresh
runs (64 constant calls per run). Instrumentation at that lock transition
records stale precision alongside wrong results, for example old=1152,
current=1408 followed by pi=1.27014. Re-reading `cprec` after reacquiring the
read lock changes this matched build to 24 successful runs, with all 1,536
constant checks passing. This establishes the observed cache race; it is not
an audit of all MPFR behavior.

The temporary build, original/instrumented/repaired sources, C reproducer and
logs are at:

`/var/folders/93/r1t98m8d3014c24llb8j119w0000gn/T/hybrid-mpfr-race-mkgvgt31`

Key logs: `native-original.log`, `native-instrumented.log`, `native-fixed.log`.
The [experimental native patch](patches/mpfr-4.2.2-shared-cache-precision.patch)
is retained for review against MPFR 4.2.2 (LGPL-3.0-or-later). It is not applied
by this package and has not been submitted upstream.
The repair changes only the precision refresh in the native library. No
package-local mutex or alternative logarithm has been added.

## Integration and remaining work

A temporary shared library is linked against Julia's existing GMP and selected
only for isolated test processes. Loaded-library paths are checked before any
such result is attributed to it. These checks do not establish that an ordinary
Julia launch uses the repair.

The unchanged radiation candidate passes 12 fresh four-thread processes per
Julia version with the temporary library (386 assertions per process). Logs:
`julia-patched-radiation112-1.log` through `-12.log` and the corresponding
`113` files under the build root (4,632 assertions per version).

The selected upstream `tversion`, `tconst_pi`, `tconst_log2`, `tlog`, `tlog1p`,
`texp` and `texpm1` executables all exit successfully against the repaired
library. Their logs are `upstream-focused-v2-<name>.log`. The first upstream
check invocation ran 198 programs and reported 194 lock-abort failures, one
pass and three skips. Its test allocator declared `mpfr_lock_memory` without
initializing it. The temporary test support now initializes that lock in
`tests_memory_start` and destroys it in `tests_memory_end`; the numeric and
memory assertions remain intact. The original failed logs are retained in
`upstream-transcendental-tests.log` and `mpfr-4.2.2/tests/test-suite.log`.
Only the seven selected programs were rerun successfully; a complete upstream
suite pass is not claimed.

Retain the native patch and reproducer
for review. A distributable runtime correction or an explicit supported
concurrency contract is still required before closing this item. Do not hide
the failure with warm-up, retries, a package-only lock, or weakened assertions.
The unrelated [Julia MPFR/GC deadlock report](https://github.com/JuliaLang/julia/issues/62753)
describes a different locking failure; its proposed repair has not been used
as evidence that this wrong-result race is fixed.
