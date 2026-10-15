# CI reliability and current Julia compatibility

GitHub CI must pass for each pushed change. Billing blocks are reported only
when supported by the run's diagnostics; they do not excuse code failures.

## Verified failures at fadc03f

[Run 38710456393](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38710456393)
passes formatting and the 2-/4-rank jobs. The latest-Julia job and 8-rank job
fail; minimum-Julia testing is still running. Neither failure is a billing block.

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
- [ ] Commit and push the verified fixes.
- [ ] Observe all GitHub jobs green for the pushed revision.
