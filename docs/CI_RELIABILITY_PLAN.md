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
has passed formatting and the 2-/4-rank MPI jobs; both package jobs and the
8-rank MPI job are running at the last observed check. GitHub acceptance remains pending.
