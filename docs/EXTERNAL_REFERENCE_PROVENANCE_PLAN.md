# External reference provenance — CAP-007A

Status: local implementation and validation gates passed; exact-SHA GitHub CI
remains pending. Production dispersion and full shock-profile comparisons remain
CAP-007B/C.

## Pre-edit brief

- Mode: repair the reference build/cache boundary and implement provenance.
- Requirement: cases 27/28 consume only output tied to pinned source, actual
  input, tools, build commands and output checksums. Changes invalidate reuse;
  failures retain their cause and artifacts. Never install global tools or
  remove existing caches.
- Verified baseline: `build_nhds.sh` clones moving HEAD, removes `.git`, and
  reuses an executable by existence. `build_vpic.sh` clones a moving branch,
  removes `.git`, can install OpenMPI, and deletes build/run directories.
  Both `run.jl` files reuse output by existence and have empty build catches.
- Canonical owner: a validation-local Julia module, using SHA/TOML stdlibs.
  Both shell entry points and Julia cases delegate to it; there is no package
  runtime dependency or second cache policy.
- Contract: detached exact-revision sources, out-of-source builds, preserved
  licenses, explicit command/environment/tool identities, complete artifact
  checksums, exclusive writer lock and atomic publication. Invalid cached
  artifacts fail visibly; a changed input/tool identity creates a new attempt.
  Failed attempts remain available and never replace a successful receipt.
- Compatibility: keep both documented shell entry points and all physics
  gates. Legacy cache directories remain untouched and are not trusted.
  New output paths come from the verified receipt instead of hard-coded paths.
- Implementation: pin both fetched commits in one lock file; build NHDS with
  its configure/Makefile and Hybrid-VPIC through upstream CMake's `USER_DECKS`
  hook (CMakeLists.txt:194), copying the deck and sibling `injection.cxx` into
  the attempt. Preserve unmodified upstream trees. No source patch is needed.
- Verification: real builds/runs with installed toolchains, cache replay and
  mutations of source/input/tool/output/receipt, interrupted/failed command
  publication, lock contention, path escapes, spaces in paths, and both Julia
  versions. Register offline cache-boundary tests in the package suite.

## Guardrail exceptions

The shell wrappers retain documented commands and delegate all behavior to one
owner. A typed failure at the external-process boundary carries the command,
exit cause and retained log; validation converts only that typed failure into a
failed metric. Interruptions and programming errors propagate.

## Independent checks

Tests must reject plausible incorrect implementations that (1) check only file
existence, (2) fingerprint inputs but ignore changed output/source, (3) publish a
partial failed run, or (4) bypass verification in one of the Julia callers.
Hashes establish local reproducibility metadata, not external authenticity or
bitwise reproducibility across machines. Physical acceptance remains owned by
the unchanged case metrics.

## Progress

- [x] Inspect both builders, runners, documentation, validation driver, CI and
  upstream build systems. Fetch exact source revisions while retaining `.git`.
- [x] Locate actual licenses: NHDS `LICENSE` (BSD-2-Clause), Hybrid-VPIC
  `README.md` copyright section (BSD-3-Clause).
- [x] Implement shared owner and migrate both shell/Julia callers.
- [x] Run offline mutation/failure checks on Julia 1.12.7 and 1.13.2.
- [x] Build and execute both pinned external references; verify cache replay.
- [x] Review physics metrics, formatting, final diff and documentation.
- [ ] Push accepted slice to main and review all six exact-SHA CI jobs.

## Verification ledger

- Offline V3: 42 checks pass on each supported runtime, including source,
  output, manifest, revision, path, lock, process and tool-diagnostic failures.
- Pinned NHDS and Hybrid-VPIC both build and run in frozen candidate V2.
  System `mpirun --version` fails with PMIx runtime/build incompatibility and
  exit 213. An isolated relocation of already-installed MPICH artifacts passes
  compiler/version and real two-rank Allreduce checks, then runs the unchanged
  Hybrid-VPIC deck successfully. No global tools were installed or modified.
- Independent V1 boundary witness: four failed assertions on each Julia.
  A receipt could select an existing source file as output, and `mv(force=true)`
  replaced an existing directory. Repair binds output to the source lock and
  uses the package publication owner's same-directory `Base.Filesystem.rename`
  pattern, with no delete/copy fallback. Original logs remain retained.

Evidence directory:
`../HybridPlasmaPIC-audit-evidence/2026-10-23/reference-provenance`.
Command and snapshot records include `offline-v3-commands.json`,
`real-build-v2-commands.json`, `build-candidate-v2-freeze.json`,
`isolated-mpich-toolchain-v2.json` and `independent-boundaries-v1-commands.json`.
Final V3 adds 50 registered and four independent boundary checks per runtime.
Both real reference builds complete, and all ten unchanged physics metrics pass
on Julia 1.12.7 and 1.13.2 (`final-v3-commands.json` and complete case logs).
NHDS maximum relative differences are 0.0122201550 at high k and 0.0494183000
over the full range; shock overshoot/compression relative differences are
0.181943324 and 0.0543385109. These remain scoped oracle/scalar comparisons.

The final V3 runner reused one artifact directory per runtime; the single-case
driver intentionally cleans CSV files, so case 28 replaced case 27's CSV files.
All complete logs and reference outputs remain retained. Final replay will use
separate per-case artifact directories. The path repairs below still require
final-source integration checks before promotion; those checks are recorded below.

## Build-path repair brief

- Mode: repair the canonical NHDS build command.
- Reproduction: `real-cache-v1-112.log` fails before replay at configure;
  the retained `attempt-kNXCOH/build.log` reports `unsafe srcdir value` for
  an absolute cache path containing spaces. Pinned `configure:1951–1961`
  permits spaces in the working directory but rejects them in `srcdir`.
- Contract: caller-selected cache directories may contain spaces, and source
  remains unmodified. Build failures retain their attempts and never publish.
- Change: invoke the configure script by `../source/configure` from the owned
  build directory, so upstream receives a safe relative source path.
- Verification: rebuild and replay real references in cache paths containing
  spaces on both supported Julia versions; preserve V1 failure. Run tool/deck
  mutation, command-failure publication, caller and offline boundary checks.
- Harness correction: establish a runtime-specific baseline before replay
  assertions. The recorded `Sys.MACHINE` and PATH differ between the original
  builds and this invocation; changed identities correctly request a build.
- V2 confirms that the relative configure path repairs the cache-directory
  failure: ten replay/corruption checks pass. A compiler executable whose name
  contains spaces then fails. Its retained `build/config.log:89–103` shows
  Autoconf splitting the FC command at the first space.
- Guardrail exception: the NHDS adapter needs one owned compiler launcher on
  its build PATH because upstream interprets FC as command text. The launcher
  executes the recorded compiler path through a quoted environment variable;
  it contains no compiler policy, fallback or retry. Record its bytes and the
  effective environment in the hashed build log. Preserve unchanged source,
  original failures and actual-tool hashing. Real compiler paths containing
  spaces and shell metacharacters must execute literally in integration checks.

## Final V5 acceptance evidence

- All 50 registered and four independent boundary assertions pass per runtime
  (`final-v5-offline-commands.json`). Formatter/AST/idempotence: 10 checks pass;
  both shell scripts pass `bash -n`.
- Real-cache integration passes 44 assertions per runtime: replay, corruption
  rejection, retained previous receipts, actual compiler/deck invalidation and
  failed compilation publication. Cache paths contain spaces; the successful
  compiler filename also contains a semicolon, dollar sign and apostrophe.
  Both builds preserve the unmodified upstream source checkout.
- These integration commands retain their subsequent harness failures: the
  first caller probe used a bare module without `include`; its replacement
  treated the documented relative artifact basename as an absolute path.
  The corrected probe uses the real shared-module loader and resolves artifacts
  relative to their output directory. All 16 checks pass per runtime with
  `--depwarn=error` (`caller-quality-v2-commands.json`); Aqua passes 9/9 each.
  No production assertion or physics tolerance was weakened.
- Final frozen-source physics replay builds both references and passes all ten
  original metrics on each runtime. Each case now retains separate CSV outputs,
  raw reference outputs and its provenance manifest (`physics-v5-commands.json`,
  `physics-v5-source.json`). These metrics match the values recorded above.
- Both actual shell entry points replay their verified caches without creating
  an attempt or changing the current receipt (`shell-replay-v5-commands.json`).
- `final-review-v5.json` records all source hashes, 20 passing physics metrics,
  167 checked relative documentation links and unchanged held root files. The
  final diff changes validation/test dependencies, builders, callers and docs;
  production solver source and physical acceptance thresholds are unchanged.

The local system MPI failure and all product/harness failure records remain
retained. No global tool installation, automatic retry or cache deletion was
introduced. Full package tests and complete compiler-log review remain required
on the pushed GitHub revision.
