# Whistler packet allocation plan

Current follow-up (2026-10-27): [bounded BigInt counters](WHISTLER_BIGINT_COUNTER_PLAN.md)
extend multiplier reuse to positive BigInt counts fitting Int. Generic
arithmetic and other counters retain the original path. The sections below
preserve the earlier V2 campaign and its then-current compatibility decision;
its blanket BigInt exclusion is superseded only within that exact boundary.
Original integration is included in green `f830b2a`. Subsequent packet-allocation
CI failures and raw Ubuntu observations remain in the
[CI investigation](WHISTLER_ALLOCATION_CI_DIAGNOSIS.md).

This is a separate optimization of the frozen linear-CN V5c candidate. It does
not alter scalar CN arithmetic, frequency construction, physical tolerances or
the original pending latency gates. Evidence lives in
`linear-cn-range/packet-hoist-v2`; the failed V1 candidate and all `repair-v5`
sources/harnesses remain unchanged. V2 functional, allocation, retention, compiler
and formatting checks pass on Julia 1.12.7 and 1.13.2. Landed on `main` after the
linear CN repair `03e4eb3` (2026-10-26): `test_whistler_packet_allocation.jl`
(5,793 assertions), the linear CN/phase suites, Aqua and JET pass on both
runtimes; the largest measured n=128/100-step call drops from about 692 MB to
9.4 MB. The hoist landed as `cf2c4c5`; the containing `4724947` passes
[all six jobs and complete compiler review](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40325153253).
That containing-revision result does not establish the outcome of the original
hoist revision or close the separate latency work.

## Pre-edit brief

- **Mode:** measured allocation optimization at the canonical packet driver.
- **Baseline:** isolated worktree starts at local snapshot `fbc7bed` of the nine
  hash-verified V5c files (original base `c5e8472`). The V5c resource campaign
  measured complete n=128/100-step exceptional packet allocations up to
  691,889,568 bytes on Julia 1.12.7 and 691,889,696 on 1.13.2. Both raw campaigns
  and their independent resource audits are retained. No latency claim follows.
- **Verified owner:** `run_whistler` in
  [semi_implicit.jl](../src/integrators/semi_implicit.jl) constructs a private
  `Vector{Float64}` of frequencies, a private seeded `ComplexF64` packet and one
  fixed timestep. Its step-major loop recomputes the same multiplier for every
  mode on every timestep. The only public caller is
  `compare_integrators_whistler`; validation case 19 and the CN/phase tests call
  the packet entry point directly. The shock CN solver is independent.
- **Invariants:** there is no cross-mode coupling in this loop. Each amplitude
  must undergo the identical ordered sequence of `bhat[m] *= multiplier`
  operations and Float64 stores; replacing the sequence with a power is not
  equivalent. Frequency construction, RNG draws, band filtering, both energy
  reductions, returned fields and validation order stay unchanged. Zero steps
  must not evaluate a multiplier. Nonrepresentable frequencies continue to
  raise their existing `ArgumentError` before packet construction.
- **Compatibility decision:** hoist only for concrete built-in Float16/32/64
  timesteps with fixed-width integer step counts (including Bool). Their
  validated finite values and the local Float64 frequency
  determine a fixed multiplier for this call. Keep the original step-major
  path for all other timesteps or step-count types, preserving precision-sensitive
  arithmetic, custom method order/count and arbitrary-precision iteration costs.
  No global precision, rounding or exponent settings are modified. Both CN and
  Euler use the same structural optimization; Euler's existing nonfinite results
  remain visible.
- **Strategy:** guard zero steps, interchange the two loops for hardware
  timesteps and retain one scalar multiplier per mode. Allocate no multiplier
  array, add no public API and change no scalar/frequency formula. The generic
  loop is a deliberate compatibility boundary, not a retry or approximation.
- **Predeclared allocation targets:** compare complete warmed public calls
  against frozen V5c on both runtimes. Every original ordinary default/override
  packet must use no more bytes. For all nine retained exceptional frequency
  cases at n=4/128, allocation at 100 steps must equal allocation at one step
  for hardware timesteps, and one-step allocation must not exceed V5c. For the
  seven frequency cases whose retained 100-step/one-step allocation ratio
  exceeds 20, require at least 95% reduction at both n values. Both remaining
  cases must not regress. No timing acceptance is authorized in this slice.
- **Verification:** compare exact complete returned values at every prefix of
  ordinary and exceptional histories, including negative/zero/subnormal/huge
  timesteps, both methods, signed zero, odd/even grids, zero steps, band limits,
  seeds and integer parameter overrides. Independently derive rational CN and
  frequency references for finite cases. Compare exception type/message and
  input/settings before/after; preserve unsupported/nonrepresentable cases.
  Exercise BigFloat stored/current precision, rounding modes, rational/integer
  timesteps and a custom `Real` with observable arithmetic call order. Run
  permanent CN tests, public comparison/case-19 callers, concrete JET and
  configured formatting on 1.12.7/1.13.2. Keep sources, manifests, commands,
  outputs and failed evidence; root owns integration and main push/CI.
- **Wrong patches rejected:** multiplier exponentiation changes rounded
  histories; eager multiplier construction changes zero-step/error behavior;
  caching arbitrary custom `Real` changes arithmetic evaluation counts/order;
  adding a multiplier vector violates ordinary allocation gates.

## Tracker

- [x] Read canonical source, public callers, tests, validation case and CI.
- [x] Verify and snapshot V5c without changing its source or frozen harness.
- [x] Freeze compatibility and allocation targets before production edits.
- [x] Implement the bounded loop change and permanent history/error tests.
- [x] Preserve V1's confirmed BigInt step-count allocation regression and repair
  the admission boundary; retain a mode-count-independent iteration-cost test.
- [x] Reverify complete histories and all allocation targets on both runtimes.
- [x] Reverify generic inputs, settings and nonrepresentable/error paths.
- [x] Run neighboring caller, compiler and formatting checks; review final diff.
- [x] Freeze exact patch and evidence for root review in `packet-hoist-v2`.
- [x] Root integration and individual main push (`cf2c4c5`).
- [x] Containing `4724947`: all six jobs and complete compiler review, as linked above.
- [x] Recheck the original `cf2c4c5` [run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40275325475): both package jobs failed; formatter and MPI 2/4/8 passed. This historical failure remains failed.

Latency acceptance remains open under the unchanged V5c contract. The complete
plasma-simulator roadmap remains separate from this linear reference optimization.

## Retained V1 evidence and regression

V1 passed its original hardware-timestep inventory but was not eligible for
promotion: extending the final resource review to BigInt step counts found
eight regressions. For n=128 and `nsteps=BigInt(100)`, V5c uses 43,040/43,168 bytes
on 1.12.7/1.13.2; V1 used 2,463,152/2,463,280. Subtracting each matching
fixed-width packet proves that iterator allocation was multiplied by the mode
count. Julia's ordinal-range iteration adds a step to its BigInt counter; V1
repeated that allocating iteration once per mode. The failed source/diff,
commands and both probe records remain in `packet-hoist-v1`.

V2 requires a built-in fixed-width counter before interchanging the loops.
BigInt/custom counters keep the original order without narrowing their values.
The extended acceptance target is no increased complete-call allocation for
all tested counter types at zero, one and 100 steps, plus the unchanged original
hardware packet targets. V2's source-supported pre-edit brief is retained in
`packet-hoist-v2/PRE_EDIT_BRIEF.md`.

## Verified V2 results

Both Julia 1.12.7 and 1.13.2 pass 5,805 assertions in the registered packet-test
context, including Aqua and package JET after the custom `Real` extension.
Independent V5c comparison passes 2,568 assertions per runtime, retains 1,562
output rows, compares all prefixes 0–100 of nine exceptional n=128 packets and
preserves all 15 recorded exception types/messages. BigFloat rounding/precision
and custom `Real` arithmetic order/count remain unchanged. An explicit control
rejects replacing ordered products with a multiplier power.

The complete-call allocation campaign passes all 76 profiles and 882 assertions
on each runtime after eight warmups and three alternating paired samples. All
20 original ordinary/default/override profiles have unchanged bytes; all 18
exceptional one/100-step pairs now allocate equally. Every one of the 14
expensive profiles exceeds the predeclared 95% reduction target; the minimum
reduction is 98.4896% on 1.12.7 and 98.4888% on 1.13.2. The largest measured
V5c call drops from 691,889,392 to 9,400,816 bytes on 1.12.7, and from
691,889,520 to 9,400,944 on 1.13.2. These measurements use the new direct-call
comparison boundary; the earlier V5c resource-boundary figures remain retained.
An independent record audit rejects 12 mutations per runtime.

The extended counter campaign passes 22 profiles and 176 assertions per runtime,
covering Int8, UInt128, BigInt and Bool counters. All paired allocations equal
V5c, including 43,040/43,168 bytes for n=128 and `BigInt(100)` steps. Measured
samples have zero compilation/recompilation. A separate audit rejects ten
record mutations per runtime. The permanent regression test requires BigInt
iterator overhead to remain independent of mode count.

The unchanged pinned GMP observer passes 20 profiles and 40 windows per runtime,
covering 18 exceptional packets and two BigInt-counter packets, with one and
eight complete calls. All requested/usable native deltas are zero; the maximum
corrected Julia live delta is 64 bytes under the original 64 KiB gate. The
maximum raw live delta is 8,644,608 bytes; raw values and exact accounting drift
remain recorded. Held-buffer controls reject retained memory and verify
release; an independent resource audit rejects 11 mutations per runtime. This
measures allocation and retention, not exact transient memory peaks or execution
latency.

Each runtime also passes 27 concrete JET checks, nine Aqua checks, 16,782 existing
CN/phase assertions, and all eight metrics from the actual validation-case-19
caller. Complete compiler trees match V5c after checked unchanged-line mappings:
29 reports on 1.12.7 and 27 on 1.13.2; no report was suppressed or left unmapped.
Final formatting checks pass on both runtimes for all three affected Julia
files, with parsed syntax equality checked. The six-file final diff contains
only the canonical packet loop, its contract, tests and documentation.

The first V1 formatter attempt used an existing auxiliary environment whose
CSTParser 3.5.0 source was unavailable. Its failure is retained as
`packet-hoist-v1/format-v1.log`. The copied, pinned evidence formatter environment
uses JuliaFormatter 1.0.62 with the installed CSTParser 3.5.1. No global package
or tool installation was performed. V5c and the rejected V1 source/harness
hashes are checked again in the final handoff freeze.

Generic timesteps and arbitrary-precision/custom counters intentionally retain
the original evaluation path. Nonrepresentable frequencies still fail visibly;
Euler's existing overflow results remain observable. Latency acceptance,
exact transient native peaks, root integration and GitHub CI are not claimed
by these local allocation and functional checks.
