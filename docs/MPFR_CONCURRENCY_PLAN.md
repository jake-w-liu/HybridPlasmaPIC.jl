# Concurrent BigFloat runtime investigation

Status: the native MPFR shared-cache defect is reproduced on Julia 1.12.7,
1.13.1, 1.13.2 and 1.13.3. Fresh private-cache native controls pass their bounded
suite, but the copied-runtime experiment stops at its fresh-worker exponent
range gate. The canonical Julia integration source is sealed and reviewed.
The unchanged 198-program native suite passes 195 tests plus three original
format skips against its exact private GMP. Private MPFR installation, default
dependencies, release/debug core compilation and strict link/signature checks
pass. Corrected prebootstrap fixtures pass the stock negative control and both
actual cores: 12 serial thread contexts, 36 state snapshots and 517 exact symbol
owner records. The helper widens only the positive controls' exponent ranges,
preserves flags/precision/rounding and makes zero measured GMP allocations.
Original loader and test-fixture failures remain retained. Static-analysis
dependency installation passes. The actual MPFR analysis stops on two
GCChecker diagnostics for external calls from helpers annotated NOTSAFEPOINT.
The reviewed repair retains the early initializer's contract with eight local
constant/TLS metadata declarations and removes the cleanup annotation because
cleanup calls the current GMP free callback. Its existing caller enters a live
GC-unsafe region before teardown. This two-file repair is being prepared in a
new isolated source copy; the failed analysis and original source are preserved.
No checker exemption or runtime-safety conclusion follows. Threading
analysis, Julia bootstrap, numerical consumers, GC/lifecycle, P4 child/depot and
platform gates remain unaccepted.
Installed Julia libraries are unchanged and NUM-001
remains open.

The earlier experiments below are historical records. Their temporary build
paths are absent in the takeover environment and are not fresh acceptance
evidence. The [fresh native/runtime tracker](#fresh-native-build-and-integration-tracker-2026-10-27)
and [sealed integration design](#sealed-canonical-julia-integration-design) distinguish
retained native results, failed runtime gates and unimplemented work.

The radiation-candidate references below describe historical experiments.
The landed fixed-width particle radiation path uses native/scaled arithmetic
and exact rational invariants without BigFloat. That change does not close
the runtime defect. The [BigFloat timestep repair](BIGFLOAT_TIMESTEP_CONTEXT_PLAN.md)
assumes correctly rounded MPFR primitives; task-scoped numerical settings do
not repair shared native caches.

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

## Julia 1.13.2 refresh

**VERIFIED:** the freshly installed official Julia 1.13.2 still reproduces the
constant-cache failure. One fresh four-thread invocation of
`julia +1.13.2 --startup-file=no --threads=4 tools/audit_mpfr_cache.jl` exits 1,
with two incorrect pi and two incorrect log(2) values among 128 checks. The
reported library is the installed 1.13.2 `libmpfr.6.dylib`, MPFR 4.2.2, TLS=1,
shared=1, with an empty patch string. For example, a pi request at 6,912 bits
returns about 2.72846. Full values, precision, library path and exit evidence:
`/tmp/hybrid-mpfr-repo-probe1132-1.log`. No retry, runtime-library replacement or
package-local masking was performed. Runtime integration remains open.

## Julia 1.13.3 refresh

**VERIFIED:** one predeclared fresh four-thread invocation of
`julia +1.13.3 --startup-file=no --threads=4 tools/audit_mpfr_cache.jl`
exits 1 with two incorrect constants among 128 checks. The pi request at
5,760 bits returns approximately 3.768553263367864; log(2) at 6,656 bits
returns approximately 0.7999924657811566. The loaded stock Julia 1.13.3
`libmpfr.6.dylib` reports MPFR 4.2.2, TLS=1, shared=1 and an empty patch string.
The source, command, exit status, loaded path and complete log are retained
under `HybridPlasmaPIC-audit-evidence/2026-10-27/takeover/mpfr-latest-refresh/`.
No retry or runtime-library replacement was used. This reproduces the existing
wrong-result behavior on the latest tested runtime; it does not re-prove the
native cause independently of the earlier instrumented controls.

## Fresh native build and integration tracker (2026-10-27)

Earlier temporary paths above record historical experiments. The three named
temporary build/handoff directories are absent in this takeover environment;
those runs are not reused as fresh acceptance gates. A new isolated build uses
official Julia 1.13.3 source `67ac47b9` and unmodified MPFR 4.2.2 source/tests.
Its sole tracked recipe change selects `--disable-shared-cache` while retaining
thread-local state and the canonical default GMP 6.3.0+2 dependency.

**VERIFIED native results:** `make -C <isolated-source>/deps -j1 check-mpfr`
exits zero: 195 programs pass, three disabled-format tests skip, and none fail
or error. The 572 original source/test files remain unchanged. Sixteen
predeclared serial C processes, eight per variant with no retries, retain all
logs: stock fails in five processes with 21 incorrect pi values out of 512;
the private-cache build returns 512 correct values across eight processes.
Every process records the expected MPFR and GMP loader paths. Identity probes
verify version, TLS/cache mode and public layout, not complete semantic ABI
equivalence. Root independently checks all raw control logs and 1,974 preserved
source/evidence pins.

Evidence: `mpfr-runtime-build-r1-r2-v1/` and
`takeover/mpfr-native-root-review-v2.json` in the dated external audit archive.
The first configure run's unresolved native GMP rpath and a root pin-checker
format error remain retained; neither is numerical acceptance. The successful
build uses a native CPU target and absolute isolated GMP rpath. It is a native
experiment, not a portable artifact or a corrected Julia runtime. Installed
Julia libraries remain unchanged.

| Stage | State | Acceptance still needed |
| --- | --- | --- |
| R1 source/recipe | Native build and source pins checked | Portable target and relocatable dependency configuration |
| R2 native suite/controls | Complete upstream suite and fixed native sample pass | Recheck against the exact GMP selected for any new runtime artifact |
| R3 runtime selection | Copied official prefix selects the exact corrected MPFR/GMP pair; H stops at fresh-worker range gate (exit 3) | Canonical Julia-owned thread initialization and input admission; no numerical consumers after H |
| R4 ABI/GC | Bounded loaded Julia/native layout, symbol owners and allocator identities checked before the H stop | Actual allocation/finalization/GC and thread-exit behavior of a newly implemented canonical runtime |
| R5 concurrency | Stock Julia still fails; corrected native sample passes | Fresh-thread Julia constants, precision/rounding/exponent contexts and all declared processes |
| R6 consumers/resources | Pending corrected Julia runtime | Timestep and remaining BigFloat consumers; bounded storage/cost gates |
| R7 artifact/launcher | Copied-prefix provenance/signing and scoped MPFR/GMP relocation measured; full-prefix G remains unmet | Broad dependency closure and SuiteSparse operations, clean caches and distributable immutable selection |
| R8 child/MPI | Unimplemented | Every child/rank uses the same verified artifact |
| R9 platforms/CI | Unimplemented | Minimum/latest Julia and every claimed target/platform |
| R10 closure | Open | Fresh independent audit; no claim that ordinary stock launches are repaired |

## Sealed canonical Julia integration design

The source/contract phase is frozen in the dated external evidence directory
`mpfr-runtime-julia-integration-contract-v1/`. Root independently checked the
1,983 evidence pins, 1,818 unchanged Julia source pins, 43 source files / 85 cited
spans and the review hashes. The whole-evidence manifest SHA-256 is
`9b7330865a6c001ce1dc9e92fd8a937acdbd6e2da2acd6b9d65845f8fc24389e`.
This records the design review. The subsequent static implementation and
native build gates below do not establish a loaded runtime correction.

The checked call inventory covers 132 public MPFR symbols at 114 owned Julia
sites. The static candidate implements 111 typed same-native-frame forwarders, one
operand-role owner, intrinsic TLS/private-cache mode rejection during initialization, and
Julia-owned initial-thread/worker context setup. Adopted and existing foreign
callback contexts remain caller-owned. Sign-only copysign donors, metadata
queries, raw limb/IEEE/identity paths and arbitrary-format internal helpers keep
their explicit exceptions. No arithmetic lock, warmup, retry or per-call context
reset is proposed.

Bounded native probes establish zero GMP allocation for the measured scalar
getter/setter hook phase, nonempty local-cache cleanup without new allocation,
and valid sign-donor/alias/format-role behavior. They do not establish Julia GC,
foreign destructor ordering, permanent-worker lifetime bounds, other platforms
or concurrent numerical safety. Cache cleanup belongs at actual registered
thread exit while allocator accounting is live; suspended workers and
`jl_atexit_hook` do not imply thread exit or zero retained caches. Pre-init host
GMP allocator ownership remains an explicit external-client contract.

The isolated static candidate is sealed in `mpfr-runtime-julia-candidate-v1/`;
root checks 4,107 evidence pins and all 1,796 source entries. Its actual native
preparation/suite/install/dependency and release/debug core gates pass with
canonical matched GMP/native MPFR and an isolated depot. Root independently
checks all 755 continuation pins and the complete three-process prebootstrap
state/symbol evidence in `takeover/mpfr-prebootstrap-v2-actual-root-review-v1.json`.
The stock negative helper rejects shared caches; both actual cores accept the
private-cache configuration. These probes create no MPFR numerical objects and
do not initialize Julia. Actual Julia GC/lifecycle and numerical behavior remain
unaccepted. Direct GC-thread snapshots require a separately reviewed diagnostic
core variant, not user tasks scheduled onto GC workers. Fresh raw thread state
must pass before numerical consumers. The P0–P10 implementation/loader/context/
resource/consumer/platform gates are preserved in that sealed handoff; R3–R10
above remain open. No package runtime-distribution API, installed/default
runtime replacement, upstream submission or public artifact is supplied.
