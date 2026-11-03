# Whistler allocation CI investigation

The ±16/32-byte packet variation on Ubuntu now has a measured cause; see
[usable-size accounting](#usable-size-accounting-2026-11-03). The larger
historical failures below (4,768 and 6,352 bytes) are not attributed to it, and
their cause remains unresolved. The later
[bounded-counter optimization](WHISTLER_BIGINT_COUNTER_PLAN.md) is accepted at
`cd8b94d`, with six green jobs and complete compiler review; that result does
not establish the cause of earlier failures. Exact revision `a7fb3c8` fails both package
jobs in [run 40393658443](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40393658443):
`test_whistler_packet_allocation.jl:158` measures mode-count differences of
4,768 bytes on Julia 1.13.3 and 6,352 on 1.12.7 against the unchanged 4,096-byte
limit. Formatting and all three MPI jobs pass. No billing block was observed.
This does not establish the cause or invalidate the numerical RK4 checks.

The packet implementation and permanent test are byte-identical at the prior
green `f830b2a`, failed `a7fb3c8` and `0f78536`. Twelve predeclared local commands
all exit zero, including eight observations on both source revisions and both
runtimes, with and without coverage. All three samples per fixture are kept;
BigInt overhead is mode independent locally. These packet-only observations
do not reproduce Ubuntu or the preceding exact-arithmetic and linear-CN tests.
Evidence is in `takeover/whistler-iteration-ci` under the dated audit archive.

## Diagnostic contract

The manual-only [workflow](../.github/workflows/WhistlerAllocationDiagnostic.yml)
compares the two frozen source revisions on Ubuntu with exact Julia 1.12.7 and
1.13.3. It resolves the declared test environment under each source Project's
compatibility and pins registered SpectralOperators 0.1.6, as both failed jobs
did. Each artifact preserves the complete resolved manifest, source hashes,
runtime/options, original test output, all three paired allocation samples per
fixture, every sampled allocation and stack frame, and typed code.

The diagnostic runs the original first three test files in order. It retains
all assertions and the 4,096-byte gate; an outer testset allows the subsequent
observations to be written before a failed assertion gives the job a nonzero
exit. Allocation tracing occurs after the original tests and raw samples.
Trace-instrumented allocation data is separate from the ordinary measurements.
The six-job main CI workflow is unchanged and remains required for every push.

## Tracker

- [x] Preserve exact failed CI logs and source identity.
- [x] Check every local command, log hash and raw paired sample.
- [x] Predeclare Ubuntu source/runtime matrix, original prefix and trace scope.
- [x] Verify the diagnostic scripts locally: 71,638 assertions per runtime;
  all raw samples and complete trace packets checked. Wrong-source and
  forced-assertion pipeline controls preserve visible nonzero failure. Both
  Julia files pass formatting/syntax/idempotence checks; the original formatter
  harness line-location mismatch is retained.
- [x] Push diagnostic `f93e17d`; its exact-SHA run retains a latest-Julia package failure and five passing jobs in the [audit tracker](DEEP_DEBUG_TRACKER.md).
- [x] Dispatch and review all four Ubuntu artifacts; preserve raw byte variation and the unsuccessful historical gate.
- [x] Establish the cause of the ±16/32-byte packet variation with measured
  evidence: Linux usable-size accounting, identical requested allocations.
- [ ] Attribute the historical 4,768/6,352-byte failures. Usable-size
  accounting has not been measured on that BigInt-counter path, so it is not
  assumed to be their cause.
- [x] Implement and reverify the separately measured bounded-counter optimization with all existing gates retained; this does not attribute the historical variation.
- [x] Push `cd8b94d`; [all six exact-SHA jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40423035764) and complete compiler review pass.

No tolerance increase, disabled coverage, test skip, selected minimum, retry to
green or production optimization is part of this investigation slice. The
original failed revision remains failed regardless of later results.

## Usable-size accounting (2026-11-03)

`test_linear_cn_range.jl` compared one `@allocated run_whistler(...)` sample
for each ordinary parameter type with one Float64 sample. It failed
intermittently on Ubuntu by exactly 32 bytes on both runtimes (for example
`23360 == 23328` in [run 41655715723](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41655715723)).
Taking the minimum of three samples (`12586fc`) did not remove the failure.

Julia v1.12.7 and v1.13.3 count malloc-backed array buffers and counted mallocs
(GMP limbs) by `memory_block_usable_size`: `malloc_usable_size` on Linux and
`malloc_size` on macOS (`src/gc-common.c`, `src/gc-stock.c`). On Linux, the
counted size is the usable size of the chunk glibc returns, which can exceed
the request.

The Ubuntu [diagnostic run 41671275834](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41671275834)
used the CI test flags on both runtimes:

- The test's single-shot comparison mismatched 22 of 480 times on 1.12.7 and
  43 of 480 on 1.13.3, by ±16 or 32 bytes; some differences were negative.
- The profiled requested allocations (type and size, sample rate 1) had one
  distinct inventory. The requested bytes (23,000 / 23,128) and allocation
  counts (27 / 30) were identical for the Float64 and integer parameters.
- Measured bytes minus the usable-size slack of the reachable malloc-backed
  buffers were constant (for example, 8,136 bytes = 8,080 + 56 for
  `Vector{Float64}(undef, 1002)`).

HPIC therefore requests identical allocations for these parameters; the byte
difference lies in what the counter adds for them. The packet check now
compares the requested allocation records and allocation counts. Byte counts
of malloc-backed buffers are not compared exactly. Twenty repeats of the
revised file in one Ubuntu process ([run 41672612207](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41672612207))
passed 335,360 of 335,360 assertions on each runtime, including all 480 packet
comparisons.

Julia pools a `Memory` only when its data plus 16-byte header fit in
`GC_MAX_SZCLASS` (2,024 bytes; `src/genericmemory.c`). Larger buffers, such as
the packet's 2,048-byte field and its RNG state, are malloc-backed. The only
other equality between nonzero `@allocated` counts in the suite,
`frequency_bytes` in the same file, measures pooled frequency vectors of at most
1,024 bytes. Every remaining `@allocated` assertion checks zero or an upper
limit.
