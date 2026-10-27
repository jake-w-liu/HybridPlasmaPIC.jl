# MPFR 4.2.2 shared-cache precision race

Historical prepared report, not submitted. Its fixed-sample counts refer to the
recorded builds; some temporary paths are no longer available. Use the
[current investigation and sealed design](MPFR_CONCURRENCY_PLAN.md) for fresh
evidence and remaining runtime work.

Prepared report scope: The reproducer uses C and pthreads; it does not
load Julia or HybridPlasmaPIC. The source inspected is the official MPFR 4.2.2
release, not current development HEAD.

## Observed behavior

Concurrent calls to `mpfr_const_pi` at increasing precisions intermittently return
grossly incorrect values with thread-local state and shared caches enabled.
For example, one fresh process using Julia 1.12.7's MPFR library reports:

```text
mpfr=4.2.2 header=4.2.2 gmp=6.3.0 tls=1 shared=1 patches=
BAD_PI precision=128 observed=1.2701436221184967 expected=3.1415926535897931
```

The [standalone reproducer](../tools/audit_mpfr_cache.c) starts four workers
together and makes 64 requests, covering 128 through 8192 bits. Each worker owns
its destination. A result that is nonfinite or differs from the independent
`acos(-1.0)` reference by more than `1e-12` fails the process. This tolerance
detects gross corruption; it does not verify every rounding bit. A passing
process cannot establish thread safety because the failure is intermittent.

## Environment and fixed sample results

- `uname -srm`: `Darwin 25.5.0 arm64`.
- Compiler: Apple clang 21.0.0 (`clang-2100.1.1.101`), POSIX thread model.
- Strict build flags: `-std=c11 -O2 -Wall -Wextra -Werror -pedantic -pthread`.
- MPFR runtime and header: 4.2.2; GMP: 6.3.0; runtime patch string: empty.
- The original native source build used `--enable-shared-cache
  --enable-thread-safe --disable-shared --enable-static`, with GMP selected
  through `--with-gmp`.

The final reproducer was compiled separately against each library below.
Actual MPFR and GMP loaded-library paths were checked in a separate diagnostic
process before eight fresh sample processes per variant. No failed sample was
replaced by a retry.

| MPFR library | TLS / shared cache | Processes with incorrect constants |
| --- | --- | --- |
| Bundled with Julia 1.12.7 | 1 / 1 | 1 of 8 |
| Bundled with Julia 1.13.1 | 1 / 1 | 6 of 8 |
| Temporary 4.2.2 build with precision-refresh patch | 1 / 1 | 0 of 8 |
| Unmodified 4.2.2 built without shared caches | 1 / 0 | 0 of 8 |

Each control therefore passes 512 constant checks in this fixed sample. These
counts describe the observed schedules, not an estimated failure probability.
The installed Julia libraries have not been replaced.

## Reproduction command

The command below was checked with explicit MPFR/GMP headers and library files
on the system above. Set the include/library variables to the installation
being tested and `probe_binary` to a new output filename. Keep MPFR and GMP
ABI-compatible. The library arguments are full file paths.

```sh
cc -std=c11 -O2 -Wall -Wextra -Werror -pedantic -pthread \
  -I"$mpfr_include" -I"$gmp_include" tools/audit_mpfr_cache.c \
  "$mpfr_library" "$gmp_library" \
  -Wl,-rpath,"$(dirname "$mpfr_library")" \
  -Wl,-rpath,"$(dirname "$gmp_library")" \
  -lm -o "$probe_binary"
"$probe_binary"
```

Check loaded paths as well as reported versions when multiple installations
exist. On the tested macOS system, `DYLD_PRINT_LIBRARIES=1` records those paths.
Preserve the exit status and output of every fresh process. Do not use repeated
execution to turn an observed failure into a passing result.

## Verified cause and experimental repair

**VERIFIED:** in the official 4.2.2 source, `src/cache.c:107` downgrades the
write lock through `MPFR_LOCK_WRITE2READ`. The POSIX definition at
`src/mpfr-thread.h:137` unlocks and then acquires a read lock. Another writer can
change the cached precision and mantissa during that interval. The local
`cprec` then remains stale when rounding resumes.

Instrumentation immediately after this lock transition recorded:

```text
CACHE_PREC_RACE old=1152 current=1408
BAD_PI precision=1024 value=1.2701436221184967
```

The [experimental patch](patches/mpfr-4.2.2-shared-cache-precision.patch) reads
`cprec = MPFR_PREC (cache->x)` again after the read lock is acquired. In the
earlier matched native build, the original code failed in 16 of 24 fresh
processes; the patched code passed all 24 processes (1,536 constant checks).
The final reproducer's independent samples are reported above.

The unmodified build with `--disable-shared-cache --enable-thread-safe` also
passes MPFR's complete upstream suite: 195 programs passed, three unsupported
format tests skipped, no failures or errors. The patched shared-cache build
reaches the same result only after a separate correction to the test allocator's
lock initialization and lifetime; it is not an unmodified-test-suite result.
The [local investigation record](MPFR_CONCURRENCY_PLAN.md) preserves those
earlier test failures and the support changes separately from the cache patch.

This report establishes a wrong-result race in the inspected release and
builds. It does not claim that current development HEAD is affected, that the
issue is unknown upstream, or that the experimental patch audits all MPFR
concurrency. MPFR's [reporting instructions](https://www.mpfr.org/mpfr-4.2.2/mpfr.html#Reporting-Bugs)
identify the information and standalone testcase requested for an upstream
report.
