# Complete legacy transverse-current filtering

Status: final native acceptance passes on Julia 1.12.7 and 1.13.2: 1,090
focused/resource/concurrency, 14,100 quality/numerical and 29,473 complete
EM/EM1D/checkpoint assertions per runtime. All concrete/package compiler reports
match baseline. The tested repair is pushed as `b73abfe`; its exact-revision
[GitHub run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38777494822)
is fully green: 369,711/369,712 assertions on actual Julia 1.12.7/1.13.2,
plus formatting and 2-/4-/8-rank MPI. [CI evidence](CI_RELIABILITY_PLAN.md#acceptance-run-at-b73abfe).

## Pre-edit brief

- Mode: numerical range repair in the existing legacy EM current filter.
- Acceptance: for range-sensitive finite Float32/64 vectors, round each complete
  Nyquist-removed cell once. Preserve ordinary operation order, odd-length no-op,
  empty vectors, nonfinite propagation, particle restoration and public step phase.
- Baseline: 68 defect-witness assertions pass on 1.13.2, terminal exit 0, in
  `/tmp/hybrid-nyquist-range-baseline1132.log`; the same 68 witnesses also pass on
  1.12.7, terminal exit 0, in `...baseline112.log`. At lengths 4/6/8, alternating
  `±floatmax(T)/2` should filter to zero but produces alternating infinities.
  Alternating least-subnormal/zero values should each round to even zero but
  remain unchanged. Public four-cell `step_empic!` with stationary NGP electrons
  and alternating transverse velocities produces infinite Jy and Ey instead of
  zero for both types; particles remain unchanged.
- Cause — VERIFIED: `src/models/full_pic.jl::_zero_nyquist_mode!` (lines 614–629
  at `e9c4fa7`) accumulates the alternating sum and divides it in T before cell
  correction. Its caller `_deposit_Jy_at_midpoints!` filters at line 1302 after
  restoring particle positions; the public step consumes the returned current.
  The independent baseline uses an exact rational alternating-sum oracle.
- Contract: with even positive n, the complete result is nearest-even rounding
  of `f[i] - (-1)^(i-1) * sum((-1)^(j-1)*f[j])/n` from original stored values.
  Select that expression after intermediate range loss or a conservative final
  subtraction overflow bound. Ordinary cancellation accuracy remains separate.
  True final overflow stays visible; no clamping, changed mode policy or tolerance.
- Owner: only `_zero_nyquist_mode!`; reuse existing `_empic_exact_sum` on disjoint
  odd/even views and `_empic_complete_charge_shift!` with opposite signs. Add a
  dedicated permanent test file and its include. No second rounding algorithm.
- Compatibility: preserve ordinary and nonfinite results, including the checked
  Float16/BigFloat private controls. Hardware cold results intentionally change. No
  public API, dependency, schema, cache or particle-storage changes.
- Strategy: collect eligibility, accumulator range status and largest magnitude
  before mutation. Check divided mean including loss to zero; detect possible
  final overflow using `largest > floatmax(T)-abs(mean)`. For selected finite
  Float32/64 inputs, form the exact alternating mean and apply the shared complete
  shift owner to disjoint odd/even cells. Otherwise execute the original loop.
- Verification: exact independent rounding cells; positive/negative overflow,
  subnormal ties, nonzero DC, mixed signals, odd/empty/nonfinite and generic types;
  original ordinary bits and allocations; public steps, particle/error boundaries,
  affected EM and restart tests, concrete/package compiler/Aqua and formatting on
  1.12.7/1.13.2, then exact-revision CI. Reject zeroing the vector (nonzero DC),
  checking only infinities (subnormal ties), or rounding the mean first (half-ulp
  controls). Preserve the original failing evidence and unchanged tolerances.

The formatter-idempotent candidate is `/tmp/hybrid-nyquist-v1`; all 181 Julia
source/test hashes are frozen in `/tmp/hybrid-nyquist-frozen-hashes.json`.
Full-PIC SHA256: `fc992cbbec802c4d3cc0a7f3d0ff7acec3a88aac40d0ed22752778895e1dabce`.
Its first permanent test run passes 596 exact range/tie checks, 148 compatibility
checks and 120 public-step checks, terminal exit 0. Both native final combined
numerical/compiler/quality gates and the complete EM/EM1D/checkpoint suites
use that frozen source. Logs: `/tmp/hybrid-nyquist-{quality,integration}{112,1132}.log`.

## Resource and additional boundary checks

Both versions pass 126 complete-call resource/compatibility assertions and four
concurrency assertions, terminal exit 0. Ten ordinary cases retain bitwise
outputs and **0 B** allocation. Eight alternating `±floatmax(T)/2` overflow
cases allocate **2,912 B**
(Float32) or **5,280 B** (Float64), independent of the tested lengths 4–4,096.
The concurrent controls use 22 independent vectors, with 20 calls per task;
results match serial capture and original inputs/ambient precision are unchanged.
Logs and parsed measurements: `/tmp/hybrid-nyquist-resources{112,1132}.{log,json}`.
Ordinary shared-host median ratios span 1.072–2.450 on 1.12.7 and 1.175–4.186 on
1.13.2; retaining zero allocation does not remove this measured kernel overhead.
No portable timing bound or full-step overhead claim follows.

Mixed normal/subnormal vectors expose a larger cold cost. Both versions pass
64 independent rounding/repeat/state assertions over eight cases (lengths
8/32/128/1,024), terminal exit 0. At 1,024 cells, complete calls allocate
**1,794,064 B** for Float32 and **4,243,792 B** for Float64. These aggregate
allocations are not peak memory or a general bound; the equal-amplitude cases
above do not characterize all exceptional inputs. Logs and driver:
`/tmp/hybrid-nyquist-wide-resources{112,1132}.log` and `...wide-resources.jl`.

A further 16-check control on 1.13.2 rejects removing the intermediate range
flag: inputs `[u,1,u,-1]` and `[3u,1,3u,-1]`, where u is the least subnormal,
have a tiny prefix that later normal arithmetic loses. The complete final small
cells lie at subnormal midpoints. The candidate rounds them correctly; a mutant
checking only the final mean fails all eight output comparisons. Evidence:
`/tmp/hybrid-nyquist-prefix-control{.jl,1132.log}`. The two corresponding vectors
are added to the final permanent-test artifact, while the running frozen source
and integration snapshot remain untouched. The final permanent file passes
**880 assertions per version**, terminal exit 0: 612 exact range/tie, 148
compatibility and 120 public-step checks. Its SHA256 is
`eb64cf229b4de4940d0b072ad68ed6c480fb86cb75bd82323f5ab3186503656c`;
logs: `/tmp/hybrid-nyquist-final-test{112,1132}.log`.

## Compiler gate investigation

The initial combined gate passes all numerical testsets but exits 1 on each
version: 22 of 28 concrete calls have zero reports; the six no-ion public-step
calls each report `_require_subluminal_velocities(::Nothing, ...)` and
`nparticles(::Nothing)` in existing `es.mobile` branches. The initial logs remain
at `/tmp/hybrid-nyquist-quality{112,1132}.log`. These were initially unclassified;
the completed comparison below establishes baseline equivalence for all six calls. `/tmp/hybrid-nyquist-compiler-v2.jl` analyzes
and executes the same concrete calls against unchanged `e9c4fa7` source and the
candidate, records complete reports, and runs package/Aqua checks. Report
equivalence and valid/invalid-ion execution must be observed before promotion;
no package ceiling, permanent assertion, or diagnostic is suppressed.

The first baseline execution harness also failed once on 1.12: it called the
private midpoint deposition phase without filling its midpoint work arrays.
The public step fills them before that call; initialization does not promise
those scratch values. The retained v2 log reports a nonfinite position read from
that scratch (`/tmp/hybrid-nyquist-compiler-baseline112.log`). The v3 harness
explicitly copies the stationary particle positions into both midpoint arrays
before executing the private phase. Neither production source nor tests change.

The v3 minimum-version rerun additionally exposes uninitialized particle
positions in the harness before `init_empic!` (19 passes/one error). VERIFIED:
`ParticleSet{D,T}(N)` allocates x with `undef` in
`src/particles/particle_set.jl:57`; the v4 harness now explicitly sets both
species' positions, velocities, weights and midpoint inputs before execution.
Earlier harness logs remain retained and do not count as successful gates.
Final driver: `/tmp/hybrid-nyquist-compiler-v4.jl`.

## Compatibility refinement before promotion

Pre-edit brief extension: the final caller/type review found that the existing
private filter also accepts complex vectors. VERIFIED: on 1.13.2, the original
kernel returns finite results for four ComplexF32/ComplexF64 fixtures while the
isolated candidate throws `MethodError` from its new unconditional `max` scan
(`/tmp/hybrid-nyquist-generic-probe1132.log`, eight witness assertions).
The range selector is intended only for Float32/64. Restrict its finite/magnitude
scan by the compile-time element-type predicate, leaving the original generic
arithmetic intact. Add complex empty/odd/even/nonfinite compatibility controls.
This changes no hardware range expression or public signature. Freeze a v2
snapshot; verify the final permanent tests and compiler/quality gates on both
versions, compare generated hardware code, and run affected integration before
promotion. The running v1 snapshot remains untouched.

The v1 full EM/EM1D/electromagnetic-checkpoint suites finish successfully:
**29,113 assertions on each runtime**, terminal exit 0, in the recorded
`/tmp/hybrid-nyquist-integration{112,1132}.log`. They predate the complex guard.
The next frozen acceptance snapshot composes that guard with the separately
verified dyadic-extraction optimization (PERF-005), whose hunk is disjoint in
`full_pic.jl`. Run final tests, package/concrete diagnostics and the complete
integration suite on this combined source. Each repair retains a separate
commit and push, with its own GitHub run; no held radiation changes are included.

## Completed compiler comparison and final composed snapshot

The corrected v4 baseline and v1-candidate compiler processes all finish with
exit 0. Each executes 34 valid/private/error controls and passes package JET
(three assertions) and Aqua (nine). All 28 full concrete report renderings match
unchanged `e9c4fa7`, including the six no-ion calls with two reports each; all
22/21 package terminal reports and 14/13 complete trees match on 1.12.7/1.13.2.
This establishes baseline equivalence for those findings, not a claim that all
compiler findings are harmless. Logs/artifacts:
`/tmp/hybrid-nyquist-compiler-v4-{baseline,candidate}{112,1132}.{log,toml}`;
comparisons: `...-comparison{112,1132}.log`.

The final composed source is frozen at `/tmp/hybrid-nyquist-dyadic-v1` with
181 source/test hashes in `...-v1-hashes.json`. Full-PIC SHA256:
`eeb61b0a42802dcbc5d98a6e0d6931065fb359e40e192273ab0f1c29251220db`.
The final permanent file adds complex compatibility cases and contains 896
assertions. Together with complete-call resource, wide-range and four-thread
checks, it passes **1,090 assertions per version**, terminal exit 0. Logs: `/tmp/hybrid-nyquist-dyadic-focused-resources{112,1132}.log`.
Final quality gates use the same frozen source and finish with exit 0: **14,100
assertions per version** (14,046 numerical, 34 concrete-call execution, eight
additional concrete compiler, three package JET and nine Aqua). All 28 concrete
report renderings and 14/13 complete package trees match unchanged source;
terminal report counts remain 22/21. Comparison logs:
`/tmp/hybrid-nyquist-dyadic-compiler-comparison{112,1132}.log`. The final complete
EM/EM1D/checkpoint integration gates pass **29,473 assertions per version**,
terminal exit 0 (23m23.1s/23m32.6s). Logs:
`/tmp/hybrid-nyquist-dyadic-integration{112,1132}.log`. All 181 frozen source/test
hashes remain unchanged after these gates; the final patch replays exactly on
`afc324f`, including the four-file formatter-idempotence check.

Final composed resource results retain 0 B for all ten ordinary cases and
2,912/5,280 B for the eight pure-alternating Float32/64 cases. Observed ordinary
median timing ratios are 1.245–2.933 on 1.12.7 and 1.477–2.112 on 1.13.2;
shared-host timings are descriptive. Parsed complete-call records:
`/tmp/hybrid-nyquist-dyadic-resources{112,1132}.json`.

Generated LLVM was captured for the v1/v2 Float32/64 kernels on both versions.
It is not text-identical after function-name normalization (block labels and,
on 1.13.2, GC frame assignments differ), so it is not used as an equivalence
proof. Acceptance uses the completed numerical/integration runs above. Raw IR
and differences: `/tmp/hybrid-nyquist-ir-*.ll`, `...diff`; comparison JSON:
`...-comparison.json`. The first capture harness hit Julia's method-world guard;
the corrected driver evaluates definitions before a new top-level capture loop.
Both capture logs are retained under `/tmp/hybrid-nyquist-hardware-ir*`.

An isolated IEEE-bit range-predicate optimization was measured and not adopted.
Each version passes 69 controls: all 65,536 Float16 encodings, 203,584 Float32
and 228,672 Float64 sampled/boundary encodings, full filter comparisons, zero
ordinary allocations and consumed-output checks. Seven-pair timing medians span
0.850–1.311 on 1.12.7 and 0.953–1.147 on 1.13.2 relative to the frozen kernel,
without consistent improvement. The shared predicate remains unchanged. Evidence:
`/tmp/hybrid-range-flag-{brief.md,probe.jl,probe112.log,probe1132.log}`; both
processes exit 0. The measured ordinary filtering overhead remains open.
