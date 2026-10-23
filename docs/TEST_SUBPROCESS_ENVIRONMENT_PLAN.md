# Metadata subprocess test environment

Status: 2026-10-23; final local validation passes. Separate push and exact-SHA
GitHub CI remain required.

## Pre-edit brief

- Mode: test-harness repair; production serialization is unchanged.
- Requirement: both checksum subprocesses must use the resolved parent test
  environment and load the same package checkout, including when that checkout
  has no root Manifest.toml.
- Baseline: `../HybridPlasmaPIC-audit-evidence/2026-10-23/full1132.log`
  records 599,361 passing assertions, zero failed assertions and one error on
  Julia 1.13.2. This full-suite run failed; it is not a green run.
- **VERIFIED cause:** `test/test_metadata.jl:118,129,144` selects
  `pkgdir(HybridPlasmaPIC)` as the child project. The recorded child command
  switches from the resolved external environment to the manifest-free
  ionization checkout and fails to load FFTW before reaching the checksum test.
- Contract: retain independent writer/reader processes, both-direction state
  and metadata assertions, native child failure propagation and temporary-file
  cleanup. `src/io/metadata.jl:167,191` remains the persistence owner.
- Canonical surface: the two subprocess launches in `test/test_metadata.jl`.
  The endpoint fixture in `test/test_exact_square_root.jl:228` already uses the
  active project for this same environment boundary.
- Compatibility: no public API, file format, dependency or production change.
- Strategy: select `Base.active_project()`, require a project, and pass the
  expected package path to each child for an explicit source-identity check.
  Disable startup files, matching the isolated endpoint fixture.
- Verification: retain the original failure; run the complete metadata file
  from resolved external environments on Julia 1.12.7 and 1.13.2 with no
  checkout manifest. Also exercise a resolved checkout environment on each
  runtime, check formatting and the complete diff, then require all six
  exact-commit GitHub jobs.

## Progress

- [x] Trace both child commands and the actual full-suite failure.
- [x] Repair the launch environment without changing checksum assertions.
- [x] Pass external-environment and checkout-environment metadata runs on both runtimes.
- [x] Verify formatting and unchanged production files; held hashes are checked at promotion.
- [ ] Commit and push separately; require exact-SHA GitHub CI.

Skipping the child, ignoring its exit status, and selecting an arbitrary
installed checkout would each invalidate this test. The final diff must retain
both `run` calls and all persistence assertions; child source-identity checks
must reject a mismatched checkout.

## Final local evidence

Both native versions reproduce the original FFTW dependency-load failure in a
manifest-free checkout selected by the old child command. The repaired complete
metadata suite passes 80 assertions from each external and checkout environment
on both Julia 1.12.7 and 1.13.2. Each command also passes six source-identity
checks: both actual child scripts exit unsuccessfully with the expected error
when given a mismatched expected package source. Those intentional child errors
remain visible in the logs.

After composition onto the field-energy commit `5c2796b`, both external-environment
commands pass the same 86 assertions again. All six final command/log hashes,
the frozen source manifest and 356 unchanged committed Julia files are verified.
Formatting and parsed-semantics checks pass. The only changed Julia file is the
metadata test; both persistence `run` calls, state checks, metadata checks and
temporary-file cleanup remain present. No runtime dependency or production
serialization owner changed.

Persistent evidence:
`../HybridPlasmaPIC-audit-evidence/2026-10-23/metadata-environment`, including
`baseline-commands.json`, `commands.json`, `composed-commands.json`,
`final-verification.json` and the full logs. The earlier full-suite error remains
recorded in `../full1132.log`; focused replay does not turn that failed run green.
