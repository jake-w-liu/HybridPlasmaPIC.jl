# Whistler allocation CI investigation

The allocation failure remains open. Exact revision `a7fb3c8` fails both package
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
- [ ] Push the diagnostic and current status to main; require all six CI jobs.
- [ ] Dispatch once and review all four Ubuntu artifacts, including failures.
- [ ] Establish the allocation owner and cause with measured evidence.
- [ ] Repair that owner, retain all existing gates and reverify both runtimes.
- [ ] Push the repair separately and require exact-SHA green CI.

No tolerance increase, disabled coverage, test skip, selected minimum, retry to
green or production optimization is part of this investigation slice. The
original failed revision remains failed regardless of later results.
