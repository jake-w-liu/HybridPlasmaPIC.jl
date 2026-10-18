# CI reliability and current Julia compatibility

GitHub CI must pass for each pushed change. Billing blocks are reported only
when supported by the run's diagnostics; they do not excuse code failures.

## Verified failures at fadc03f

[Run 38710456393](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38710456393)
passes formatting and the 2-/4-rank jobs. The latest-Julia, minimum-Julia and
8-rank jobs fail. No billing block is evidenced. The minimum-Julia suite finishes
with 352,282 passing assertions and five failures in 65m57.7s; JET and Aqua pass.
Its complete log is `/tmp/hybrid-ci-38710456393-julia112.log`.

- **JET compatibility:** the fresh Julia 1.13.1 test resolver permits only
  JET 0.12, while `Project.toml` permits 0.9–0.11. The resolver diagnostic is
  retained in `/tmp/hybrid-ci-38710456393-julia-latest.log`. JET's
  [v0.12.3 manifest](https://github.com/aviatesk/JET.jl/blob/v0.12.3/Project.toml)
  supports Julia 1.12 onward; its current documentation explicitly covers 1.13.
- **MPI fixture ordering:** the publication test reads preserved rank files,
  then rank zero renames their directory and replaces it with a symlink. Its
  intervening broadcast does not ensure all ranks have finished the reads.
  GitHub's artifact records missing rank files 5–8. An isolated 8-rank replay
  delays the last reader and widens the rename/symlink interval; it reproduces
  missing files and a missing-generation directory listing. A 2-rank replay
  without the widened interval passes; that pass does not disprove the race.
  Logs: `/tmp/hybrid-ci-38710456393-mpi8-artifact/mpi-smoke.log`,
  `/tmp/hybrid-ci-mpi-race-baseline112.log`, and
  `/tmp/hybrid-ci-mpi-race-baseline8-112.log`.

- **RayCon normalization comparisons:** minimum-Julia CI fails four comparisons
  in `test/test_raycon.jl`: conversion tau/eta2, a saddle coordinate and beam
  amplitude. Native runs do not reproduce those exact comparisons, but an
  independent one-ulp launch probe reproduces sensitivity at the magnetic
  derivative owner. Existing tolerances remain unchanged; details follow below.
- **KdV error ordering:** filtered and unfiltered soliton errors differ by about
  4.03e-16, reversing an asserted ordering. This is separate from the confirmed
  nonlinear filtering defects below. No physical accuracy bound is relaxed.

## Pre-edit brief

- **Mode:** repair required CI dependency compatibility and fixture ordering.
- **Acceptance:** fresh resolution and actual JET analysis on Julia 1.12/1.13;
  all rank readers precede the fixture rename; GitHub jobs pass.
- **Baseline:** resolver failure and missing-file artifact/replay above.
- **Root causes:** excluded compatible JET series; missing read-to-rename barrier.
- **Contract:** retain the JET ceiling of 30, both concrete loader analyses,
  all MPI assertions, and the full rank matrix. Required JET import failures
  must fail the gate rather than skip it.
- **Owners:** `Project.toml`, `test/test_jet.jl`, and the publication fixture in
  `test/test_mpi_checkpoint_collections.jl`.
- **Compatibility:** add JET 0.12 to the existing range; no public API, serialized
  format, solver behavior or production MPI synchronization change.
- **Strategy:** require JET import and place an explicit barrier immediately
  before rank zero renames the test generation directory.
- **Verification:** separately resolved native environments; JET, Aqua and
  archive regressions; delayed 8-rank replay; full GitHub package/MPI/format jobs.

## Checks against misleading repairs

Raising the report ceiling would hide the dependency/inference requirement;
keep it unchanged and exercise both concrete calls. Catching an import error
would report a skipped gate; direct import now propagates it. A sleep or another
broadcast would leave MPI correctness dependent on timing; retain the widened
rename interval and delayed-reader replay to challenge the explicit barrier.
No assertions are removed or relaxed.

## Verification results

The MPI fixture repair is committed as `fadfbf1`.

The baseline delayed 8-rank replay exits 1. With the barrier, the same replay
passes all **48 assertions on each of eight ranks** on each native Julia version
and exits 0 (`/tmp/hybrid-ci-mpi-race-fixed8-{112,113}.log`). Both delays remain in the
isolated replay; neither is added to the permanent test. The latest-version
MPI run uses the earlier immutable `6341c4b` runtime environment; all 85 Julia
source files match the candidate source exactly. Each captured log contains eight
passing testset summaries and no test/error markers. The final two changed
Julia test files are unchanged by JuliaFormatter 1.0.62 at margin 100.

Fresh manifest-free environments resolve JET **0.12.3** on both native Julia
versions. The current General registry tree is
`f182e4d8e8597200dd5e99cbb1ae4037d8e9c4fb`; JET's tree is
`2c350a2325c19490c0e18ece0137366eafbc4733`.

Each native gate exits 0 with **3 JET assertions**, **9 Aqua checks** and
**80 archive assertions**. Package reports are **22 / 21** on Julia 1.12.7 /
1.13.1, below the unchanged ceiling of 30; both concrete loader signatures
report zero findings. Evidence: `/tmp/hybrid-ci-jet012-gate{112,113}.log`.
Source/test/environment hashes are retained in `/tmp/hybrid-ci-repair-hashes.json`.

The direct-import negative check exits 1 when JET is absent from the active
test environment (`/tmp/hybrid-ci-jet-required-negative112.log`); it cannot turn
a missing analyzer into a skipped success.

## Status

- [x] Read current CI logs and artifact; reproduce MPI ordering failure.
- [x] Implement compatibility/import and fixture-barrier changes.
- [x] Observe final native JET/Aqua/archive gates on both versions.
- [x] Observe fixed delayed MPI replay including all publication failure cases on both native versions.
- [x] Review exact diff, matching snapshot hashes and formatting.
- [x] Commit and push the verified fixes.
- [ ] Observe all GitHub jobs green for the pushed revision.

## KdV filtering repair — pre-edit brief

- **Mode:** repair the public KdV nonlinear operator and its acceptance tests.
- **Acceptance:** a strict retained band has no aliases from excluded inputs or
  cutoff products; linear modes and the unfiltered operator retain their behavior.
- **Baseline:** `/tmp/hybrid-kdv-filter-baseline113.log` shows changes of
  0.00150–0.00175 after one 0.001 step for excluded cosine modes that should have
  zero projected nonlinear RHS. Nyquist plus mode 7 changes by 0.00500.
- **Cause — VERIFIED:** `Nhat!` projects only the product, so high input products
  alias into retained modes. The inclusive cutoff also admits self-aliasing
  boundary modes. The derivative wave-number array zeros Nyquist and therefore
  cannot define the input projection mask.
- **Contract:** evaluate `P[(P u)^2]`, where retained integer modes satisfy
  `3abs(m) < n`. Keep the complete state for linear propagation. No signature,
  validation, zero-step or unfiltered-operator change.
- **Owner:** `src/initial_conditions/linear_waves.jl`, its docstring and
  `test/test_kdv.jl`; the only other direct consumer is validation case 06.
- **Compatibility:** intentionally correct filtered nonlinear results; retain
  all soliton shape, amplitude, phase, mass and momentum accuracy thresholds.
- **Strategy:** construct the strict mask from integer mode indices and project
  inputs before squaring as well as outputs afterward. Replace the unsupported
  error-ordering assertion with an analytic bound for the unfiltered soliton and
  direct Fourier-mode alias tests.
- **Verification:** Float32/64, even/odd/small grids, strict cutoff, Nyquist cross
  products, nonzero retained products, both flag values and high-mode linear
  propagation; both native Julia versions, downstream validation, compiler and
  formatting gates, then GitHub CI.

The soliton error-ordering assertion is not a valid filter acceptance criterion:
ignoring `dealias` would satisfy it by equality. The baseline soliton's excluded
spectral norm is 1.50e-14 of its total; native filtered/unfiltered errors are
1.216270138268942e-7 / 1.2162701385078242e-7, while GitHub reverses their ordering.
Evidence: `/tmp/hybrid-kdv-soliton-baseline113.log` and the full CI log. The new
mode tests instead reject an ignored flag, output-only filtering, an inclusive
cutoff, a Nyquist mask derived from derivative wave numbers, and a zero nonlinear
operator. Physical accuracy and conservation gates remain unchanged.

- [x] Reproduce excluded-mode aliasing and inspect all KdV callers.
- [x] Implement strict operand/product projection and direct alias tests.
- [x] Observe both native test results and red baseline tests.
- [x] Verify downstream case, concrete compiler analysis, format and exact diff.
- [ ] Commit, push and observe GitHub acceptance.

The final KdV suite passes **61 assertions** on each native version, including
**44** direct mode/flag checks. Against the old production source those mode
checks give **30 passes and 14 failures**. Logs:
`/tmp/hybrid-kdv-repair{112,113}.log` and
`/tmp/hybrid-kdv-oracle-red113.log`. The final source and tests match the tested
snapshots byte-for-byte; hashes are in `/tmp/hybrid-kdv-repair-hashes.json`.

Both versions pass two concrete JET analyses with zero reports and all **eight**
validation-case-06 metrics. Its soliton relative error is **4.31395e-8** against
the unchanged **0.002** threshold. Evidence:
`/tmp/hybrid-kdv-gates-final{112,113}.log`. An earlier gate driver omitted FFTW
from its temporary environment; its compiler checks passed, but validation
could not load. Adding that existing package dependency to the temporary gate
environment resolves the harness error; no production workaround was added.

## RayCon magnetic derivative repair

- **Mode:** repair field derivatives used by conversion and amplitude transport.
- **Acceptance:** analytic derivatives agree with independent high-precision
  differences across equilibrium scales; unchanged normalized/SI and MATLAB
  reference gates pass; one-ulp launch changes respect the existing 1e-9
  amplitude tolerance.
- **Baseline:** all 781 existing RayCon assertions pass natively on Julia 1.12.7,
  so the exact Linux failures remain unverified until CI. An independent field
  oracle yields 18 passes and 162 failures against the old derivative owner.
  A launch-rounding block yields seven passes and three failures.
- **Cause — VERIFIED:** the fixed 1e-5 finite-difference stencil amplifies
  last-bit coordinate changes in second derivatives. Changing only `_dmagnetic`
  reduces the largest measured relative focusing perturbation from 4.87e-7 to
  3.76e-13 (`/tmp/raycon-roundoff-diagnostic112.log`).
- **Contract:** derive the gradient and Hessian of the existing Solovev field
  magnitude, retain geometry chaining, public APIs, error behavior and every
  existing parity/reference tolerance.
- **Owner:** `_dmagnetic` in `src/raycon/raycon_magnetic.jl`; geometry is its only
  caller. The stencil-only `_bfield_at` helper has no remaining callers and is
  removed. Tests and port notes describe the final owner.
- **Compatibility:** intentional numerical correction, no interface migration.
- **Strategy:** differentiate the three field components and their norm
  analytically; do not round normalized inputs or widen solver tolerances.
- **Verification:** 256-bit independent differences at three length scales,
  three elongations and four locations; one-ulp launch perturbations; complete
  native RayCon suites, exact diff/format checks, combined static gates and CI.

Rounding states would leave the derivative oracle failing. Increasing a fixed
step would still fail its scale sweep. Correcting only first derivatives would
fail the second and mixed derivative checks. No fixture branches, fallbacks or
suppressed diagnostics are introduced.

Each native version passes the complete **961-assertion** RayCon suite, including
180 added oracle checks. The final eight launch assertions were then run against
the exact final block on both versions: **10/10** checks pass (two overlap the
complete suite). Existing MATLAB-reference tests contribute **376** assertions
per version. Logs: `/tmp/raycon-repair{112,113}.log`,
`/tmp/raycon-amplitude-roundoff-green{112,113}.log`. Red oracle logs are
`/tmp/raycon-magnetic-oracle-red.log` and
`/tmp/raycon-amplitude-roundoff-red112.log`.

The final production source matches the native test snapshots byte-for-byte.
JuliaFormatter 1 and whitespace checks pass. Combined compiler checks pass
as recorded below. Linux GitHub acceptance remains pending.

## Combined native acceptance

KdV is committed as `c5d7e94`; RayCon as `6a6b316`. All **264 tracked Julia
files** in both combined snapshots match committed source `4124225` byte-for-byte
(`/tmp/hybrid-ci-combined-hashes.json`). Held radiation changes are excluded.

Both versions pass **three package JET assertions** and **nine Aqua checks**.
Report counts remain **22 / 21** on Julia 1.12.7 / 1.13.1; the existing abstract
FFTW finding at `kdv_solve` remains, while concrete Float32/64 calls are clean.
Logs: `/tmp/hybrid-ci-combined-gates{112,113}.log`.

The temporary combined driver then used an unqualified RayCon type name, so that
driver exited 1 after the passing package/Aqua gates. Qualifying the existing
type fixes the harness. The corrected final concrete driver exits 0 with all
**three analyses clean** on each version (Float32/64 KdV and the magnetic
derivative owner): `/tmp/hybrid-ci-concrete-final{112,113}.log`. No production
change or repeated whole-package analysis was needed for this driver correction.

The four changed Julia files are formatter-idempotent at margin 100
(`/tmp/hybrid-ci-combined-format.log`). All verified fixes are pushed through `55d996b`.
[GitHub run 38722319490](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38722319490)
has passed formatting and the 2-/4-rank MPI jobs. The 8-rank job timed out;
both package jobs are running at the last observed check. GitHub acceptance remains pending.

## Eight-rank execution budget — CI-003

**VERIFIED:** job 116379636722 ends with GitHub's explicit 20-minute timeout
diagnostic at `/tmp/hybrid-ci-38722319490-mpi8-job.log:326`. Its captured artifact
contains passing results through collective moment-buffer validation, including
all **48 checkpoint-publication checks on each of eight ranks**. The repaired
reader barrier therefore passes this GitHub replay. The log contains no failed
assertion; it does not establish that the unfinished suite would pass.

- **Mode:** repair CI execution configuration.
- **Acceptance:** complete successful 2-, 4- and 8-rank jobs on GitHub.
- **Baseline:** the 8-rank step is interrupted after 20 minutes; the 4-rank
  step completes in 7m50s. In those jobs, ordered collection continuation takes
  11m09.8s and 3m22.1s respectively; the 4-rank final ion collection suite takes
  1m15.3s. Evidence: the two retained `mpi{4,8}-job.log` files and the 8-rank
  artifact under `/tmp/hybrid-ci-38722319490-mpi8-artifact/`.
- **Cause:** the explicit step deadline terminates the 8-rank run. The cause
  of its longer runtime is unmeasured; no deadlock or hardware-limit claim is made.
- **Contract:** retain every assertion, all three ranks, failure propagation,
  and always-uploaded logs. Continue to bound hung jobs.
- **Owner:** `jobs.mpi-smoke` in `.github/workflows/HybridPlasmaPIC-CI.yml`.
- **Compatibility:** no production code, launcher, test threshold or API change.
- **Strategy:** set the 8-rank step budget to 40 minutes, twice the observed
  interrupted duration; retain 20 minutes for 2/4 ranks. This is a bounded
  trial budget, not a measured full-runtime prediction.
- **Verification:** workflow lint and parsed structural comparison; require
  the complete GitHub suite to pass before calling the change accepted.

Dropping ranks or tests would hide missing coverage. Suppressing a failed exit
would hide defects. An unlimited timeout would lose the termination bound.
The final diff must reject all three alternatives.

- [x] Capture explicit timeout diagnostic and passing partial artifact.
- [x] Implement the rank-specific bounded budget.
- [x] Validate workflow syntax, matrix and unchanged execution semantics.
- [ ] Push and observe complete GitHub acceptance.

`actionlint` 1.7.12 exits 0 on the final workflow
(`/tmp/hybrid-ci-timeout-actionlint.log`). A Ruby YAML structural comparison
against committed `a2154cf` passes: after restoring the timeout scalar in memory,
the parsed documents are identical, and the rank matrix remains `[2, 4, 8]`.
Whitespace checks pass. GitHub acceptance is still pending; the current run is
allowed to finish before another push cancels its package jobs.

The Julia 1.13.1 job in run `38722319490` completes successfully with **352,519
assertions** in 66m46.4s and 21 JET reports. This accepts the latest-Julia
dependency repair and the included RayCon/KdV changes on the Linux runner.
The Julia 1.12.7 job also completes successfully with **352,519 assertions** in
137m59.4s and 22 JET reports. Its source log is
`/tmp/hybrid-ci-38722319490-julia112.log`, lines 125 (version), 993–994 (JET) and
1195–1197 (terminal result). The eight-rank timeout is the completed run's only
failed job; no billing failure is reported. Latest-version source log:
`/tmp/hybrid-ci-38722319490-julia-latest.log`, lines
125 (version), 974–975 (JET) and 1151–1153 (terminal package result).

## Required Aqua gate — CI-004

- Mode: repair required test-dependency error handling.
- Acceptance: missing or broken Aqua imports fail; installed Aqua executes all
  existing checks with unchanged options.
- Baseline: VERIFIED, `JULIA_LOAD_PATH='@:@stdlib' julia +1.13.1
  --project=/tmp/hybrid-threaded-alias113` reports `aqua = nothing`, one broken
  test and `QUALITY_GATE_RETURNED_SUCCESS`, exit 0. Log:
  `/tmp/hybrid-aqua-missing-baseline113.log`.
- Cause: `test/test_quality.jl` catches every import error and uses `@test_skip`.
  `Project.toml` lists Aqua in the required test target; contribution guidance
  also says this file requires Aqua.
- Invariants: retain `Aqua.test_all` and its existing exclusions; no production,
  dependency, numerical or quality threshold changes.
- Owner: `test/test_quality.jl`, included by the package test runner.
- Compatibility: standalone quality testing needs the declared test environment.
  Missing-dependency runs intentionally change from success to failure.
- Strategy: unconditional import, matching the mandatory JET gate; remove the
  catch, availability flag and skipped-success branch.
- Verification: native installed-Aqua runs on 1.12/1.13, an actual missing-package
  failure, an injected import-time failure, formatting and final diff review.
  Local results pass; final GitHub acceptance remains pending.

Keeping the broad catch, replacing it with a warning-only branch, or omitting
`Aqua.test_all` would preserve false success; the negative probes and installed
check counts must reject those alternatives.

The repaired file passes all **9 Aqua checks per native version** with Aqua
0.8.18 (both exit 0). Both negative probes exit 1 and preserve their original
errors: missing package, and `AQUA_IMPORT_FAILURE_PROBE` from an isolated test
package that throws during import. No quality check is skipped. Logs:
`/tmp/hybrid-aqua-repair-installed{112,113}.log`,
`/tmp/hybrid-aqua-repair-missing113.log`, and
`/tmp/hybrid-aqua-repair-broken113.log`.
The tested copy is byte-identical to `test/test_quality.jl`; formatter
idempotence and whitespace checks pass. The diff retains all Aqua options.

## Acceptance run at 240f77c

**VERIFIED:** at `240f77c`, `.github/workflows/HybridPlasmaPIC-CI.yml:18–20`
uses one cancellable group per branch. Pushes were held while this run was
active. CI-005 below changes the group boundary so each revision can finish.

[Run 38741178196](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38741178196)
passes formatting and all 2-/4-/8-rank MPI jobs. The
[eight-rank job](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38741178196/job/116442965228)
finishes its smoke step in **19m53s**, within the configured 40-minute limit,
and reaches the final ion-collection summaries. Retained log:
`/tmp/hybrid-ci-38741178196-mpi8-job.log`, lines 648–663. This run measures
successful completion; it does not establish a runtime ceiling or the cause of
variation between runs. The minimum-version package job now passes on actual
Julia **1.12.7**: **357,258 assertions** in **114m37.1s**, with 22 JET reports.
Its job log records the runtime at line 125 and summary at lines 1168–1172:
`/tmp/hybrid-ci-38741178196-julia112-job.log`.

The latest-version job also completes successfully on actual **Julia 1.13.2**:
**357,259 assertions** in **140m52.1s**, with 21 JET reports. Its runtime and
terminal summary are recorded at lines 125 and 1125–1126 of
`/tmp/hybrid-ci-38741178196-julia-latest-job.log`. The GitHub API reports the
whole run `completed/success`, with all six jobs successful. This exact SHA is
fully green, including the latest stable Julia; no billing exception was used.
Later revisions retain their own acceptance evidence below.

## Acceptance run at 4cd642c

At `4cd642c`, the latest-version package job in
[run 38750764319](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38750764319)
passes **358,979 assertions** in **97m39.2s** on actual **Julia 1.13.2**, with
22 JET reports. Its minimum-version package job also passes **358,978 assertions**
in **109m26.7s** on actual **1.12.7**, with 23 JET reports. Formatting and all
MPI jobs pass; the whole run is `completed/success`. This exact revision is
fully green without a billing exception. Minimum-version log:
`/tmp/hybrid-ci-38750764319-julia112-job.log`, lines 125, 1040–1041 and 1258–1260.
Source log: `/tmp/hybrid-ci-38750764319-julia-latest-job.log`, lines 125,
1012–1014 and 1205–1207.

## Acceptance run at 4c03c87

The latest-version job at scheduling repair `4c03c87` passes **357,259 assertions**
in **141m36.0s** on actual **Julia 1.13.2**, with 22 JET reports. The minimum job
passes **357,258 assertions** in **149m40.9s** on actual **1.12.7**, with 23 JET
reports. Formatting and all MPI ranks pass; the whole run is `completed/success`.
Logs: `/tmp/hybrid-ci-38750474902-julia-latest-job.log`, lines 125, 1014 and
1206–1208; `/tmp/hybrid-ci-38750474902-julia112-job.log`, lines 125, 917 and
1134–1135.
[Run 38750474902](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38750474902)
retains this revision independently of later pushes.

## Acceptance run at 6ad7276

At `6ad7276`, the latest-version package job in
[run 38754771205](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38754771205)
passes **359,871 assertions** in **97m34.8s** on actual **Julia 1.13.2**, with
22 JET reports. Its minimum-version job passes **359,870 assertions** in
**111m18.1s** on actual **Julia 1.12.7**, with 23 JET reports. Formatting and all
MPI ranks pass; the complete run is `completed/success`, without a billing
exception. Logs: `/tmp/hybrid-ci-38754771205-julia-latest-job.log`, lines 125,
1014–1015 and 1207–1209; `/tmp/hybrid-ci-38754771205-julia112-job.log`, lines
125, 1034 and 1251–1252.

The independent frozen local full-package run also completes on 1.13.2:
359,871 assertions in 155m28.1s, terminal exit 0, with 21 package JET reports.
All 179 source/test files still match `6ad7276` after the run. Evidence:
`/tmp/hybrid-full-6ad7276-julia1132{.log,-verification.json}`.

## Acceptance run at e9c4fa7

Charge-conservation repair `e9c4fa7` is pushed to main with its exact locally
verified source and permanent tests. [Run 38767093527](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38767093527)
is **completed/success**, with formatting and all three MPI jobs green. Actual
Julia 1.13.2 passes 368,456 assertions in 99m44.3s (22 JET reports); actual
1.12.7 passes 368,455 in 109m39.8s (23 reports). Logs:
`/tmp/hybrid-ci-38767093527-julia-latest-job.log`, lines 125, 877 and 1070–1071;
`...-julia112-job.log`, lines 125, 917 and 1135–1136. The original watch process
exited after an API read timeout; a fresh run-status query and final watch both
confirm success. No billing exception is used.

## Acceptance run at afc324f

The dyadic-allocation optimization is pushed as `afc324f` with all 180
source/test blobs matching native acceptance. [Run 38774781309](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38774781309)
is **completed/success**, including formatting and all three MPI ranks. Actual
Julia 1.12.7 passes 368,815 assertions in 112m01.6s; actual 1.13.2 passes
368,816 in 84m09.1s. Package JET counts remain 23/22. No billing exception is
used. Native commands and complete diagnostic comparison
are recorded in [PERF-005](CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md#dyadic-digit-extraction-allocation--perf-005).

## Acceptance run at b73abfe

The transverse Nyquist repair is pushed as `b73abfe`; all 181 source/test index
hashes match the final frozen native acceptance snapshot. [Run 38777494822](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38777494822)
is **completed/success**, with formatting and all three MPI ranks green. Actual
Julia 1.12.7 passes 369,711 assertions in 146m07.2s; actual 1.13.2 passes
369,712 in 107m10.0s. Package JET counts are 23/22. Native focused/resource/concurrency (1,090), quality/numerical
(14,100) and complete solver/checkpoint (29,473) assertions pass on each of
1.12.7/1.13.2. Full compiler reports match baseline; exact GitHub acceptance
is complete. [Contract and evidence](EMPIC_NYQUIST_RANGE_PLAN.md).

## Acceptance run at 2330379

The energy-documentation correction and progress tracker are pushed as
`2330379`. Removing docstrings leaves identical executable text in
`src/diagnostics/conservation.jl`; JuliaFormatter 1.0.62 is idempotent. The
energy balance is checked by 17 exact symbolic assertions. Its independent
[run 38781510377](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38781510377)
is **completed/success**. Actual Julia 1.12.7/1.13.2 pass 369,711/369,712
assertions in 114m31.1s/100m59.8s, respectively, with 23/22 package JET reports.
Formatting and all MPI ranks also pass.

## Acceptance run at 2bc7085

The README energy qualification and progress update are pushed as `2bc7085`.
[Run 38785224021](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38785224021)
is **completed/success**, including formatting and 2-/4-/8-rank MPI. Actual
Julia 1.12.7 passes 369,711 assertions in 94m16.2s; actual 1.13.2 passes
369,712 in 141m16.9s. Package JET counts are 23/22. No billing exception is used.

All four runs above were rechecked on October 18. Each exact SHA has six
successful jobs in `/tmp/hybrid-ci-<run>-final-status.json`. Full package logs
are `/tmp/hybrid-ci-<run>-complete-julia{112,1132}.log`; the checked versions,
assertion counts and timings are in `/tmp/hybrid-ci-complete-20261018-summary.json`.

## Preserve acceptance for every revision — CI-005

Pre-edit brief:

- Mode: CI scheduling repair for the user's requirement that every push pass.
- Acceptance: a new revision's workflow does not cancel an earlier revision;
  retain the existing jobs, dependencies, failure behavior, rank/version matrix
  and artifact publication. Observe actual GitHub completion for each push.
- Baseline: at `240f77c`, the group uses only workflow and ref with
  `cancel-in-progress: true`; earlier runs were superseded by subsequent pushes.
  Run `38741178196` still has two active package jobs, so another push with the
  same group would prevent their acceptance.
- Cause: **VERIFIED** in workflow lines 18–20 and GitHub's
  [concurrency documentation](https://docs.github.com/en/actions/how-tos/write-workflows/choose-when-workflows-run/control-workflow-concurrency):
  cancellation is scoped to equal group keys. A branch-wide key groups distinct
  commits together.
- Invariants: every distinct pushed SHA retains its own test run. The
  [GitHub context reference](https://docs.github.com/en/actions/reference/workflows-and-actions/contexts#github-context)
  defines `github.sha` as the triggering commit and permits `github` in workflow
  concurrency expressions. Duplicate runs for the same workflow/ref/SHA retain
  the existing cancellation policy.
- Owner: the single workflow-level `concurrency.group`; no application or test
  source, job, dependency, timeout or assertion changes.
- Compatibility: append SHA to the existing workflow/ref key. Independent
  revisions may execute concurrently, subject to GitHub's runner capacity;
  acceptance still requires observing all jobs, not just dispatch.
- Strategy: one group per workflow/ref/SHA. Retain `cancel-in-progress: true`
  within that exact group; do not cancel or rerun any existing workflow manually.
- Verification: actionlint, parsed YAML equality after reverting only the group
  value, distinct-key controls for current/next commits and workflow/ref changes,
  exact diff, then observe the old and new GitHub runs remain active and finish.

Simply disabling cancellation with the old default single-pending queue can
still replace a pending revision; retaining the branch key with cancellation
leaves the observed problem. Dropping tests would not satisfy acceptance. The
selected change keeps every job and separates distinct revisions explicitly.

- [x] Pass workflow lint and structural/key checks.
- [x] Commit and push only this scheduling repair (`4c03c87`).
- [x] Verify earlier and later revisions remain running independently.
- [ ] Observe every pushed revision finish green (or record actual billing evidence).

Local CI-005 checks pass: actionlint 1.7.12 exits 0
(`/tmp/hybrid-ci-per-commit-actionlint.log`), and Ruby's parsed YAML comparison
confirms the only semantic workflow change is `concurrency.group`. Explicit key
controls separate distinct SHA/workflow/ref values and the active old group,
while retaining equal keys for duplicate SHA runs. Evidence:
`/tmp/hybrid-ci-per-commit-structure.log`. GitHub execution remains pending.

Scheduling repair `4c03c87` is pushed to main. Its
[run 38750474902](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38750474902)
and the earlier `240f77c` run both finish green. They remained independently
active across later pushes, verifying dispatch, retention and completion for
distinct revisions. Later pushes retain their own acceptance requirement.

## Latest stable refresh — October 15

The [official manual downloads page](https://julialang.org/downloads/manual-downloads/)
now lists Julia **1.13.2**, released October 14, 2026. Existing 1.13.1 measurements
remain valid for that tested version; they do not establish 1.13.2 acceptance.
The explicit 1.13.2 local channel installs successfully; `julia +1.13.2 --version`
reports 1.13.2 and `juliaup status` retains the original default. The full native
`Pkg.test()` run against immutable pushed `6ad7276` subsequently passes 359,871
assertions in 155m28.1s, exit 0. Its 179 Julia source/test hashes still match
the revision. Evidence: `/tmp/hybrid-full-6ad7276-julia1132-source.json`,
`...-verification.json` and `/tmp/hybrid-full-6ad7276-julia1132.log`.
Final new-candidate checks must include 1.13.2 and minimum supported 1.12.7.
CI already requests version `1`; record the actual selected version from each
completed job log rather than inferring it from the matrix label.


The first complete latest-stable GitHub acceptance is now recorded above at
`240f77c` (actual 1.13.2, 357,259 assertions). The local watcher for that run
exited on a GitHub API read timeout; its retained watch log is not a CI failure.
Subsequent direct API inspection reports `completed/success` for the whole run
and all six jobs, and the completed job logs provide the version and test evidence.
Native `6ad7276` is complete; later pushed revisions retain their own final gates.

## October 18 pressure-feature acceptance

Documentation revision `e33935024aff59f0dbcf1da467c5e1c353fb9676` has fully green
[run 39091128686](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39091128686).
All six jobs pass. Actual Julia 1.12.7 passes 369,711 assertions in 111m46.9s;
1.13.2 passes 369,712 in 98m51.0s. Package JET counts are 23/22. Status and full
job logs are retained in `/tmp/hybrid-ci-followup-20261018-0440`, with parsed
version/count/hash evidence in `39091128686-evidence.json`.

Documentation revision `1150ed3a5d53843902d75869bf1f1ce70c6958da` also has all six
jobs green in [run 39093601976](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39093601976).
Actual Julia 1.12.7 passes 369,711 assertions in 113m02.0s; 1.13.2 passes 369,712
in 98m57.8s. Package JET counts remain 23/22. The same evidence directory retains
all six job logs and `39093601976-evidence.json` with their parsed version/count
and package-log hashes. The pressure
implementation is pushed as `013a70b9a6703769fda6e7069c8a00ce6192a8b2`;
[run 39097756477](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39097756477)
has green formatting and all MPI ranks, with both package jobs pending. API requests
briefly returned rate-limit errors; subsequent reads succeeded. This is not
billing evidence or a CI completion result. Each exact SHA must finish green.
