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

## Standalone diagnostic

Run the retained [Julia probe](../tools/audit_mpfr_cache.jl) in a fresh process:

```sh
julia --startup-file=no --threads=4 tools/audit_mpfr_cache.jl
```

It loads no plasma package, reports the actual MPFR library and build flags,
and checks 64 growing-precision values of each of pi and log(2). Incorrect
constants produce a nonzero exit with their precision, observed value and
hardware-precision reference. The tolerance detects gross cache corruption;
it is not a test of every rounding bit. A passing process does not establish
thread safety because the race is intermittent. Do not retry until green.

The final tool reproduces incorrect constants on both installed Julia versions
in `/tmp/hybrid-mpfr-repo-probe112-2.log` and `113-2.log` (exit 1).
The preceding fresh process passed on each version. With the temporary repaired
library, all 128 checks pass on each version (exit 0):
`/tmp/hybrid-mpfr-repo-probe-fixed112.log` and `fixed113.log`.
A one-thread invocation rejects the invalid probe setup before spawning tasks
on both versions: `/tmp/hybrid-mpfr-repo-probe-single112.log` and `single113.log`.

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
initializing it. The first temporary test correction initializes that lock in
`tests_memory_start` and destroys it in `tests_memory_end`; the numeric and
memory assertions remain intact. The original failed run is retained in
`upstream-transcendental-tests.log`; the generated `tests/test-suite.log` is
replaced by subsequent runs.
A subsequent full run passes 194 programs, skips three unsupported format tests,
and fails `tasin` with a lock abort (`upstream-full-v3.log`). The temporary support
incorrectly destroyed its lock at `tests_memory_end`; `tests/tasin.c:272–279`
continues using the installed allocator after that check. The support now
initializes the lock once and destroys it at process exit. The complete rerun
then passes 195 programs, with three unsupported-format skips and no failures
(`upstream-full-v4.log`, terminal exit 0). Numeric and memory assertions remain
unchanged; the allocator-support corrections are distinct from the native
cache-precision patch.

Retain the native patch and reproducer
for review. A distributable runtime correction or an explicit supported
concurrency contract is still required before closing this item. Do not hide
the failure with warm-up, retries, a package-only lock, or weakened assertions.
The unrelated [Julia MPFR/GC deadlock report](https://github.com/JuliaLang/julia/issues/62753)
describes a different locking failure; its proposed repair has not been used
as evidence that this wrong-result race is fixed.

## Runtime configuration comparison

A second temporary build uses the unmodified official MPFR 4.2.2 source with
`--disable-shared-cache --enable-thread-safe`. It is linked against the GMP
already loaded by each Julia runtime and selected only for isolated subprocesses.
The probe confirms `tls=1`, `shared=0` and the actual loaded library paths.
All 128 constant checks pass on each version. The unchanged radiation candidate
also passes six fresh concurrent processes per version, each with 386 assertions
(2,316 per version). Build and logs:

`/var/folders/93/r1t98m8d3014c24llb8j119w0000gn/T/hybrid-mpfr-unshared-cot3scwu`

Files are `julia-load112.log`, `julia-constants112.log`,
`julia-radiation112-1.log` through `-6.log`, and their `113` counterparts.
The pristine build's complete upstream suite passes **195 programs**, with three
unsupported-format skips and no failures or errors (`upstream-tests.log`,
terminal success). Unlike the shared-cache build, this run requires no test
allocator-lock support changes.
No installed Julia or package dependency has been replaced.
This tests MPFR's supported cache configuration as a possible runtime correction;
it does not establish that the distributed Julia binaries use that configuration.

## Standalone native reproducer and upstream handoff plan

- Mode: retain a native diagnostic and prepare an upstream report.
- Requirement: reproduce incorrect shared-cache constants without Julia or plasma
  code; reject nonfinite constants too, report build flags and return failure on
  any incorrect result. Keep process startup and thread failures observable.
- Baseline: the retained native reproducer fails in 16/24 original-library runs;
  instrumented lock transitions prove stale precision, while the precision-refresh
  patch passes 24/24. The existing Julia-only diagnostic remains available.
- Cause: the cache lock downgrade permits a new writer before rounding resumes,
  so local precision can disagree with the newly shared mantissa.
- Contract: four pthread workers request 64 growing-precision pi values per fresh
  process. This detects gross corruption; a passing run cannot certify thread
  safety or every rounding bit. No solver or installed runtime is changed.
- Owner/surfaces: a standalone `tools/audit_mpfr_cache.c`, this evidence plan,
  the existing experimental native patch and a reviewable upstream report.
- Compatibility: no package dependency, API, build hook or automatic runtime
  replacement. The C diagnostic requires an MPFR build with thread-local state.
- Strategy: retain the demonstrated native workload, validate finite results,
  release thread-local caches in workers and shared caches after joins, check
  pthread operations, and record compiler/library versions. Prepare the report
  according to MPFR's official reporting instructions; sending remains separate.
- Verification: strict C compiler warnings; native original, patched-shared and
  pristine-unshared controls; actual loaded library paths; bounded fresh-process
  samples with all exit statuses retained; inspect the final report against logs.

The final [C diagnostic](../tools/audit_mpfr_cache.c) builds with
`-std=c11 -O2 -Wall -Wextra -Werror -pedantic -pthread` against all four library
variants. Eight fresh processes per variant give one failure with Julia 1.12.7's
library, six with Julia 1.13.1's, and none with either temporary control
(512 constant checks per control). Separate `DYLD_PRINT_LIBRARIES=1` runs verify
the actual MPFR/GMP paths. The linker reports harmless duplicate rpath entries
for the two stock variants; compiler diagnostics are clean.
Driver: `/tmp/hybrid-mpfr-native-handoff-build.py`; build log:
`/tmp/hybrid-mpfr-native-handoff-build.log` (terminal exit 0). Per-process output,
exit statuses, loaded paths and library manifest are retained under:

`/var/folders/93/r1t98m8d3014c24llb8j119w0000gn/T/hybrid-mpfr-native-handoff-1jflmbhb`

The [prepared upstream report](MPFR_UPSTREAM_REPORT.md) includes the standalone
source, measured environment, source-level cause, experimental patch and bounded
control results. It has not been sent. The runtime integration gate stays open.

The report's shell command was extracted and executed against the unshared
control with the final C source; it compiles and passes all 64 checks (exit 0).
Log: `/tmp/hybrid-mpfr-report-command.log`.
