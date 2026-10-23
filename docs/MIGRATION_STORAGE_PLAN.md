# Particle migration storage integrity

Status (2026-10-23): **GROWTH-002 final V3 passes all local acceptance gates**
on Julia 1.12.7 and 1.13.2 and is pushed separately as `7c6e1e4`.
[Exact-commit CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39892881923)
has passed formatting; the other five jobs remain running. GROWTH-001's append validation is a separate
accepted fix (`16c9776`) and does not repair this finding.

Recovery note (2026-10-23): the uncommitted `/tmp` candidate and evidence
locations below are unavailable. No prior running handle remains accessible.
The reconstructed source and new acceptance evidence are recorded below;
historical results are not current acceptance evidence. The separate GROWTH-001 repair is now
pushed as `16c9776` with all six GitHub CI jobs green.

## Verified public behavior

The pinned-main probe at
`/tmp/hybrid-particle-migration-storage-review-20261020` finishes with 36
passing control assertions and 24 failed integrity assertions on each runtime.
It covers Float32/64 in 1D–3D. `results.json` retains commands and log hashes.

- With a view-backed tag in either logical rank, migration throws while
  assigning a vector subset to the tag field. Earlier position, velocity,
  weight and id fields have already been replaced. For the failing 3D rank,
  those lengths are zero while tag length remains one. A later-rank failure
  can also leave an earlier rank changed.
- When position x₁ and weight share storage within one population, a periodic
  wrap from 1.25 to 0.25 also changes weight from 1.25 to 0.25. The call returns
  successfully. The same weight change occurs in a second population when its
  weight shares the first population's x₁. Transport must preserve these
  weights or reject the unsupported ownership before mutation.
- Independent vector-backed controls move both particles to their expected
  ranks, preserve weights and produce positions 0.25 and 0.75.

**VERIFIED mechanism:** `src/particles/migration.jl` wraps original position
storage before classification, and `_replace_particles!` assigns component
subsets one at a time. The public probes reach the late conversion failure and
show the aliased-weight changes. These measured/source paths establish the
serial finding. Real MPI also calls these owners, but no invalid-input baseline
MPI run has been launched and no MPI hang is claimed. The repaired candidates'
rejection/recovery runs are recorded below.

## Pre-edit contract

Prepare complete migration results in private storage before publishing any
particle fields. Reuse the same subset/replacement owner in local and MPI
transport. Preserve charge/mass checks, finite-position checks, periodic and
nonperiodic classification, particle ordering, all non-position fields, public
return values and the original `ParticleSet` object identities.

The intended canonical split is preparation and publication. Preparation
indexes and converts every component to its destination field type before
constructing a replacement of that same `ParticleSet` type. A storage-type
failure must precede writes to the original population. Publication assigns
already-compatible fields; it must not perform a second indexing/conversion
pass. Inspect all private callers before replacing the existing helper.

For local migration, classify using the same wrapped-coordinate operations as
the existing owner, build kept and outbound subsets from original inputs, and
wrap only those private subsets. Append arrivals to the private kept
populations. Publish all ranks only after every result is ready. This also
allows read-only shared input arrays to be copied consistently, instead of
corrupting another field through in-place wrapping. Keep compatible custom
storage; do not impose a blanket `Vector` requirement merely to simplify the
implementation.

For real MPI, audit collective ordering before editing. Prepare local kept and
outbound state without publishing it, communicate predictable preparation
errors collectively before entering transport, then receive and assemble the
private local result. Communicate assembly rejection before publication and
final statistics. Preserve useful local exceptions and give other ranks an
explicit failure; do not swallow errors or permit a rank to skip a collective
that other ranks enter. Resource exhaustion or broken MPI infrastructure is
not a promised recoverable transaction without separate evidence.

## Required checks

- Re-run the original public witnesses and success controls on both supported
  Julia versions. Check whole-population rollback, not only array lengths.
- Exercise early and late incompatible fields, overlapping components and
  populations, compatible abstract vector fields and disjoint views, empty
  ranks, all-lost and all-kept ranks, periodic/nonperiodic layouts and
  transverse coordinates. Compare complete non-position payloads by particle
  id and preserve ordinary baseline ordering/positions.
- Inspect alias ownership after successful migration. Private preparation
  must not write through retained source parents or another population.
- Verify the MPI collective stages first, then test one-rank rejection and a
  subsequent valid collective on real 2/4/8-rank runs. Keep partial results
  unpublished when a rank rejects. Do not launch an intentionally hanging
  baseline reproduction.
- Measure preparation and transport storage/cost. Retained and outbound
  particles already require indexed copies; avoid another full population
  copy or quadratic concatenation. Report measured overhead and allocations.
- Run neighboring particle, injection, migration, hybrid/MPI, checkpoint and
  compiler/quality gates, with full compiler-tree comparisons and formatting.
  Push the verified repair separately to main and require all six exact-SHA
  GitHub jobs.

The implementation remains open. Finish GROWTH-001's acceptance independently;
its success cannot be used as evidence that migration storage is repaired.


## First private candidate and follow-up

The serial preparation prototype at
`/tmp/hybrid-migration-storage-v1-20261020` is based on the unpromoted GROWTH-001
V3 source. It prepares compatible private subsets, wraps those subsets, appends
arrivals privately and only then publishes particle fields. On both Julia versions,
the initial public witnesses and existing migration suite pass 216 assertions.
MPI collective integration is deliberately still unimplemented in this first
prototype; it is not a promotable repair.

A fresh publication-owner probe catches another case within this scope:
repeating the same mutable `ParticleSet` object in two rank slots cannot retain
independent final rank states. On each Julia version, both baseline and first
candidate retain 18 control passes and six rank-ownership failures. Distinct
populations sharing input position arrays remain valid controls. The next
candidate must reject repeated mutable owners before preparation, retain
compatible shared input arrays, and check the rank container's one-based axes.
All four owner-probe commands are complete; failures remain recorded.

For MPI integration, factor the existing checkpoint local-stage protocol into
a common collective stage with lazy context/phase error text. Preserve its
original local exception/stack, peer checkpoint message, and fatal-error
exclusions. Apply it only to local storage validation, local preparation and
local assembly; do not wrap a function that itself enters rank-dependent
collectives. Keep the existing position/species collective validation and
publish compatible prepared fields after all transport/statistics stages.
Both checkpoint and migration regressions must cover this shared-owner change.

Private preparation retains original arrays longer. Measure allocation and
peak live storage rather than assuming that unchanged indexed-copy counts imply
unchanged peak memory. Record any cost before promotion.

The second private candidate at
`/tmp/hybrid-migration-storage-v2-20261020` rejects repeated population owners,
checks the rank container axes and uses the common collective stage for local
MPI storage validation, preparation and assembly. It publishes only after
successful transport, assembly and statistics. Its initial serial suite passes
240 assertions on each runtime. Expanded serial rollback/custom-storage tests
pass 870 assertions per runtime; focused real 2/4/8-rank runs also pass on both
versions. The later transport preparation follow-up below prevents accepting
V2 as the complete repair. The incomplete snapshot-creation script failed before any native
verification; its corrected snapshot is frozen by `run.py` before testing.

## MPI byte transport follow-up: count and displacement preparation

**VERIFIED on both Julia versions:** the allocation-free count witness at
`/tmp/hybrid-mpi-count-review-20261020/probe_v2.jl` constructs the installed
MPI 0.20.28 `VBuffer` with counts `[2147483647, 1, 1, 0]`; its displacements
are `[0, 2147483647, -2147483648, -2147483647]`. No MPI collective or large
payload allocation is used. The installed `MPI/pa7Yi/src/buffers.jl:329`
accumulates offsets in `Cint`, matching that overflow. The first witness
fixture lacked the required stride method and errored before these assertions;
that setup error is retained, and V2 adds only the unit-stride contract.

The package's `_mpi_alltoallv_bytes` checks each chunk's length locally before
`Alltoall`, then relies on that implicit displacement construction before
`Alltoallv!`. This source path leaves local count/copy failures uncoordinated and
does not validate prefix offsets. This is a transport preparation gap within
GROWTH-002; no baseline hang or invalid memory access has been attempted.

Pre-edit decision: retain the existing Cint count capacity, validate nonempty
chunk displacements with checked native-Int prefix arithmetic, and pass explicit
validated displacement vectors to MPI. Empty chunks may use zero displacement;
a final nonempty chunk may end beyond Cint capacity when its count and starting
displacement each fit. Do not impose an unnecessary total-byte limit. Share the
existing local-stage failure protocol at communicator level. Perform send
count/packing preparation collectively before exchanging counts, and receive
displacement/allocation preparation collectively before exchanging payloads.
Every recoverable local exception remains an error, with the original local
exception preserved. Fatal resource/MPI failures retain their stated boundary.

Verification: pure metadata boundary tests require no huge allocations; real
MPI tests reject one-rank invalid chunk count, excessive count, unsupported
prefix displacement (where rank count permits), and packing failure before a
subsequent valid exchange. Retain the ordinary destination-routed transport
suite, checkpoint rejection tests, both Julia versions and full compiler
comparisons. The existing V2 snapshot remains unchanged while the permanent
candidate integrates this additional preparation repair.


## Permanent candidate and resource observations

V3 integrates the permanent serial/MPI tests and the checked byte preparation.
Formatting passes 350 assertions and the source-owner audit passes 92. Its
minimum-Julia focused run retains 1,260 passes and two oracle errors:
`sum(Int, empty_counts)` has no inferred additive identity. Both errors occur
in the test's expected-value expression, while the production prefix loop
handles empty counts. One-rank migration and byte preparation checks pass.
The latest and real-MPI lanes are not run after that failure.

V4 adds `init = 0` to those two expected sums, retaining all assertions.
All 88 production files are byte-identical to V3, recorded in
`/tmp/hybrid-migration-storage-v4-20261020/production-identity.json`.
Its formatting and 92 owner checks pass, followed by 1,262 focused checks on
each runtime. All six focused real-MPI commands (2/4/8 ranks on both runtimes)
also pass. The full two-rank suites on both Julia versions pass, including
checkpoint restore/publication, ordered-collection continuation, moment
validation, ion collection and evolving-pressure checks. All eight MPI command
exits, complete log hashes and frozen sources are verified. The full-suite
elapsed times are 3,304.79 and 3,297.16 seconds. V3's failed logs are preserved.

V4's first concrete compiler command retains 14 passes and six failed
zero-report assertions. Each failing MPI-migration specialization reports two
`matches_to_methods(::Nothing, ::Core.TypeName, ::Nothing)` diagnostics inside
Julia Serialization. The full log is
`/tmp/hybrid-migration-storage-v4-20261020/candidate-concrete112.log`.
No production cause or pre-existing status is inferred from that location.
Matched baseline/candidate captures of migration and the unchanged serializer
are running at `/tmp/hybrid-migration-compiler-review-20261020`. The original
failed gate remains failed; later compiler gates did not run, and the candidate
is not accepted.

The first matched capture driver itself failed before producing reports: it
used the package-report field on `JETCallResult`. V2 uses the installed JET
call-report printing interface and retains that failed driver/log separately.
All four baseline/candidate captures then complete on both runtimes.
`/tmp/hybrid-migration-compiler-review-v2-20261020/comparison.json` verifies all
12 concrete cases: serial migration has no reports, and complete serializer
trees match after source-path/equal-line normalization. Each MPI tree retains
the same serializer subtree and adds only the collective-stage/preparation-
closure frames. All typed inputs and those source locations were checked;
`prefix-review.json` records the review. The original zero-report gate remains
failed because its expectation also fails on the baseline. Latest-Julia checks
of the new helper owners pass all eight assertions; full unloaded/fixture-loaded
package captures at `/tmp/hybrid-migration-compiler-followup-20261020` now
complete all seven commands. All command/log and frozen-file hashes pass.
Both complete fixture-loaded trees match baseline: 22/21 reports in 14/13
groups, with no added or removed group. Unloaded comparisons await the original
growth captures. Final combined-source MPI regressions and transport cost
measurements remain required.

The initial resource matrix at `/tmp/hybrid-migration-resources-20261020`
passes 953,733 value/resource-observation assertions in each of four commands
(actual GROWTH-001 V3 owner versus migration V2, both Julia versions). All 54
output fingerprints per runtime match exactly. For the 1,048,576-particle 3D
case, migration allocation changes from 1,279,928,320 to 1,279,929,360 bytes on
minimum Julia and from 877,275,136 to 877,276,176 bytes on latest Julia. Whole-
process peak RSS is 788,316,160 versus 822,591,488 bytes, and 624,230,400 versus
674,283,520 bytes, respectively. These are measured process peaks including
compilation/allocator state; they do not isolate the live particle arrays.

One small-profile allocation outlier in the initial matrix requires follow-up.
The first repeated protocol retains one compilation-contaminated sample per
minimum-Julia owner (359 passes/one failed zero-compilation assertion each),
with an extra 820,416 bytes in that sample. Both latest-Julia commands pass 360
checks. This does not establish the cause of the earlier outlier, whose compile
metadata was not recorded. A second protocol collects garbage after input
construction, outside every measured call, and keeps all twelve profiles, ten
samples and value/zero-compilation assertions.
`/tmp/hybrid-migration-allocation-review-v2-20261020` completes all four
commands with 360 passes each. All 120 samples per owner/runtime retain zero
compilation/recompilation and constant allocation within each profile. The
candidate adds 352, 544, 608 or 1,040 bytes per call, depending on the profile.
This compares the actual GROWTH-001 V3 owner with migration V2; final MPI
transport costs still require measurement. No elapsed-time acceptance claim
is made.

## Composition with the revised growth owner

V4 inherits GROWTH-001 V3's newly reproduced cross-species birth alias defect
and predates the pushed rank-classification optimization. It cannot be
promoted as a complete source snapshot. Its completed checks remain evidence
for the migration changes.

The isolated `/tmp/hybrid-migration-storage-v5-20261020` snapshot applies only
V4's seven changed source/test paths to the revised growth/rank candidate.
Six paths have an unchanged parent and are copied exactly; `test/runtests.jl`
adds both migration includes while retaining the append, birth-alias and rank
tests. `integration-manifest.json` verifies this union against both frozen
parents, preserving the revised ionization guard. Dependency TOML files
preserve their resolved contents apart from the package source path. Native
checks on this combined snapshot have not run.

No migration algorithm change is intended in this composition. Before
promotion, require the revised growth parent's acceptance, combined focused
storage/rank checks, full compiler comparisons with identical fixture methods,
real MPI rejection/recovery and integration checks, and resource measurements
against the combined parent. The older cost queue still references the failed
zero-report expectation and must not provide acceptance for the new source.

V5's controller is queued behind the revised growth initial, integration,
resource and complete compiler gates, plus V4's full MPI and compiler
comparisons. Its dependency TOML structures match the revised growth parent
after normalizing only the package source path. All 11 Julia drivers parse on
minimum Julia; Python syntax, driver paths, comparison filenames and all 268
source/test/extension hashes are checked. These are preparation checks, not
native acceptance.

The queued matrix combines the original migration witnesses, permanent
migration/growth/birth-alias/rank tests and quality checks. It requires 66
direct zero-report compiler calls per runtime, complete baseline/candidate
MPI/serializer trees, full unloaded/fixture-loaded package comparisons, six
focused real-MPI runs and both full two-rank suites. The known MPI Serialization
reports are compared in full rather than assigned a false zero-report
expectation. Source review of the retained wrapper frames, final resource
measurements and exact-commit GitHub acceptance still remain separate gates.

The final resource protocol will compare the revised growth/rank parent with
V5, preserving the 54 serial Float32/64, 1D–3D, kept/exchanged/lost profiles and
the million-particle case. Exchange profiles with 32 and 4,096 particles per
rank retain ten measured samples. Every measured call follows garbage
collection outside its timing window and must report zero compilation and
recompilation. Independent field/id invariants and complete output hashes
cover every sample. Process peak RSS remains a whole-process observation.

The two-rank MPI matrix retains 0/16/1,024 particles per rank, three warmups and
16 samples. It will check complete positions, velocities, weights, ids, tags,
species properties and routing, and gather each rank's full output hash outside
the timed call. These measurements record the added collective cost; unpaired
elapsed observations cannot establish a speedup. Run them only after combined
native/compiler acceptance and review of the complete typed MPI trees.

V5's first controller stopped before native work because the growth parent's comparison
had a verified source-path spelling defect. The replacement controller
requires the corrected full comparisons, retains all original semantic gates,
and additionally requires the four direct comparisons with committed main.
Production source and Julia test drivers are unchanged.

The final measurement harness at
`/tmp/hybrid-migration-final-resources-v1-20261020` passes a minimum-Julia
preflight: seven drivers parse, 864 serial protocol assertions pass, and the
real two-rank byte-digest gather control passes. This verifies the harness,
not the final migration costs. Its eight measurement commands are queued
behind V5's 32 native commands and explicit review of the complete typed
MPI trees.

The corrected parent comparisons have completed, and V5's native matrix is
running. Its four matched typed captures and complete comparison pass on both
runtimes. All twelve cases retain the same two serializer diagnostics as the
baseline. The complete serializer subtrees match, and source review at
`parallel/mpi.jl:548–562` and `1508–1518` accounts for the two additional
collective-stage/preparation-closure frames. Every public signature and
complete captured closure type matches its inputs. Local recoverable
exceptions are rethrown after notifying peers; peers receive an explicit error.
The original zero-report failure remains retained.

`/tmp/hybrid-migration-storage-v5-20261020/typed-review-accepted.json` pins
the comparison, MPI source and all 72 capture hashes. All four command logs
and 723 frozen files were checked. This closes the manual typed-tree review;
full package comparisons, real-MPI integration, final resource measurements
and exact-commit GitHub acceptance remain required.


## Persistent reconstruction (2026-10-23)

Fresh serial probes against `16c9776` retain 30 passes and 30 failed integrity
assertions on both runtimes. They cover whole-rank rollback, within/cross-rank
aliases and repeated owners. Both native MPI metadata probes reproduce negative
Cint displacements without allocating or communicating a huge buffer.
`../HybridPlasmaPIC-audit-evidence/2026-10-23/migration/pre-edit-brief.md` records
the source-backed contract and collective error policy before production edits.

The reconstructed prepare/publish repair passes 8,432 focused assertions per
runtime and first real two-rank rejection/recovery runs on both runtimes.
V2 adds axes, disjoint-view and empty-byte controls without a production algorithm
change; its four files are frozen under `migration/v2-files`. Broader native,
six real-MPI lanes, compiler and resource checks remain in progress.

All six V2 focused real-MPI commands now exit successfully (2/4/8 ranks on both
runtimes); `migration/v2-mpi-results.json` records commands and log hashes.
This does not close the final review. The safe `position-stage-witness.jl`
singleton probe retains four passes and one failure per runtime: a custom
position getter throws inside the stage in the direct control, but outside it
through public migration. Source `parallel/mpi.jl:1177` reads positions before
the following Allreduce. This verifies an uncoordinated local failure; no
multirank hang was attempted or observed.

V3's planned repair moves the existing storage-validation stage into the
canonical migration input validator, covering both axes validation and the
finite-position scan. It preserves the following collective position/species
checks and error text without adding a collective. The permanent MPI regression
must reject one rank's read error, preserve all original populations, then
complete a subsequent valid collective. V2 snapshots remain frozen, and V3
requires its own source-bound acceptance checks.

V3 is frozen in `../HybridPlasmaPIC-migration-v3-20261023`, based on `b251cb7`.
Its position-read witness passes five assertions on each runtime; the original
storage witness passes 60. The first real two-rank run passes on minimum Julia.
A separate combined V3 snapshot includes both permanent tests in the package
and real-MPI runners. Its six formatting checks pass. All twelve fully typed
MPI compiler captures retain the same two serializer diagnostics as the parent;
the complete serializer subtrees match. The two added frames are the reviewed
local collective stage and preparation closure, with exact captured input types.
No diagnostic was suppressed. `migration-combined-v3/typed-comparison.json`
pins the source and all capture hashes.

A separate combined worktree at `../HybridPlasmaPIC-migration-combined-20261023`
uses `b251cb7` as parent, composes the four frozen files, and adds both test-runner
includes so GitHub's package and real-MPI jobs exercise the new regressions.
Its 268-file source/test snapshot is frozen; six formatting checks pass.
Its combined suites have completed; this older V2 snapshot remains separate
from the final V3 acceptance below.

### Current V3 acceptance

Both final direct suites pass 8,495 assertions, and all six focused real-MPI
commands pass (2/4/8 ranks on Julia 1.12.7 and 1.13.2). The final combined V3
worktree is `../HybridPlasmaPIC-migration-combined-v3-20261023`; its neighboring
suite passes 99,191 assertions per runtime. Both full two-rank integration
commands pass, including checkpoint continuation, collective rejection and
rollback, multispecies and ion-collection integration. The older V2 combined snapshot is retained separately
and does not substitute for V3 acceptance.

Fresh complete single-thread compiler captures match the parent in all groups
and group order (30/27 reports). The earlier four-thread equality failures are
retained: each runtime selected one different entry caller into an otherwise
identical existing diagnostic subtree. Source review verifies JET's concurrent
signature analysis and deduplication by final error frame, with no entry caller
in that key. The affected project sources are byte-identical. Both the original
differences and the complete single-thread comparisons are retained under
`migration-combined-v3`; no diagnostic is filtered or declared harmless.

Final serial resource commands pass 11,672 assertions each. All 541 complete
outputs match the parent per runtime across 55 profiles, with no measured
compilation or recompilation. The 1,048,576-particle 3D exchange allocates
270,764,976 bytes on minimum Julia and 268,929,968 on latest, versus
267,094,576 baseline bytes: increases of 3,670,400 and 1,835,392 bytes.
Smaller profiles have both increases and reductions. The million-particle case
has one measured sample after three warmups; it does not establish a latency
distribution. Whole-process peak RSS is retained as an observation.

The real two-rank resource matrix passes 496 assertions per rank per command.
All 96 complete rank samples match per runtime, including gathered output
hashes. Each profile adds 240 allocated bytes per rank per call, with zero
measured compilation/recompilation. The unpaired elapsed measurements do not
establish a speedup. `migration/final-evidence-check.json` records all 22 final
command logs and both source manifests; comparison files preserve every sample.

The final root check in `migration-combined-v3/final-verification.json` verifies
all 34 accepted command logs, the three frozen source manifests, the six owned
files and 347 unchanged committed Julia files. Complete typed/package compiler
and resource comparisons are rerun from the retained captures. Git diff checks
pass. No held radiation changes belong to this repair.
