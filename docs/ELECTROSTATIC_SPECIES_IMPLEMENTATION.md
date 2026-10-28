# CAP-008 periodic electrostatic species collections

Status: landed on `main` (2026-10-26) as the first CAP-008 slice. On `cf2c4c5`
the species, ES/EM, charge, I/O, restart, variable-step, clock, dimension,
moment, MPI-collection-checkpoint, Aqua and JET suites pass on Julia 1.12.7 and
1.13.2. V7 functional, storage, allocation and compiler results are recorded
below; coordinated timing acceptance remains open. The feature landed as
`d46029e`; its original CI run failed, while containing `4724947` passes
[all six jobs and complete compiler review](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40325153253).
That containing result does not relabel the original failure. Arbitrary EM
collections and the remaining CAP-008 integrations stay open.

## Pre-edit brief (2026-10-25)

- Mode: extend the existing periodic electrostatic solver.
- Requirement: ordered CPU Float32/Float64 populations in 1D–3D, arbitrary finite
  charge and positive mass, shared Poisson source and field energy, explicit
  background policy, synchronous leapfrog phase, and deterministic restart.
- Baseline: `electrostatic.jl` validates every scalar population as q=-1,m=1;
  `ElectrostaticPIC` stores only one gather buffer and electron density.
  `checkpoint.jl` stores a scalar electron payload; `restart.jl` delegates to it.
- Design gap: collection charge, storage, phase, neutrality and restore contracts
  do not exist. This is an extension, not a bug in the documented scalar model.
- Invariants: every population sees the same old electric field; charge sources
  sum q_s n_s; every velocity shares one half-time level; independent storage;
  invalid input is rejected before physical state changes. A staged collection
  step also prevents arithmetic failures from partially publishing particles.
- Canonical owners: `electrostatic.jl` for charge, Poisson and leapfrog;
  `checkpoint.jl`/`restart.jl` for serialization; existing kinetic/momentum owners
  for particle totals. Collection overloads reuse these numerical owners.
- Compatibility: retain the integer-count electron constructors and their n0,
  ne, Ep and checkpoint contracts. Their historical uniform neutralization is
  documented. New collection construction defaults to requiring neutral physical
  charge, with a declared roundoff tolerance. Explicit `:neutralize` adds a
  uniform compensating background, reported by diagnostics. `ne` for a collection
  denotes electron-equivalent charge density -sum(q_s n_s), so n0-ne remains the
  sole Poisson source. No EM contract is changed.
- Strategy: collection workspace attached to ElectrostaticPIC; reusable staged
  particles/fields; shared scalar/collection kick and Poisson kernels; identity
  binding on initialization; validated collection schema beside legacy scalar
  schema. All populations advance with one dt; replacing/reordering populations
  after initialization requires reinitialization. Particle growth is permitted
  only with new velocities supplied at the existing shared half level.
- Verification: pair plasma and unequal q/m, direct charge deposition, scalar
  compatibility, permutations, neutral/empty cases, field energy counted once,
  variable dt, alias and invalid-final-population rejection, cold-mode frequency
  and refinement, checkpoint/archive replay, neighboring scalar/EM tests on
  Julia 1.12.7 and 1.13.2. Formatting/Aqua and full final diff review follow.

Guardrail exception: the collection constructor/overloads extend the canonical
ElectrostaticPIC owner; the integer constructor retains its released electron
contract. There is no duplicate Poisson or pusher physics. Collection policy is
explicit because periodic Poisson cannot support a nonzero mean source. The
legacy scalar schema remains readable; collection payloads cannot be loaded as
scalar payloads. No fallback catches or weakened test gates are planned.

## Physics and representation contract

The [Smilei PIC algorithm specification](https://smileipic.github.io/Smilei/Understand/algorithms.html)
defines total charge as a sum over all species and weighted particles, gathers
one common field at the old particle positions, and advances staggered velocities
before projecting the new charge. This implementation uses those same ownership
and time-level conventions, with the repository's periodic spectral Poisson and
CIC operators. The spectral zero mode is fixed to zero; integrating periodic
Poisson over the box requires zero mean charge. `:require` checks exact represented
physical charge to the declared tolerance; `:neutralize` declares a uniform
compensator. Ordinary CIC rounding can leave a small represented grid residual,
which the zero-mode solve removes. Diagnostics distinguish that grid integral
from exact physical and compensating charges.

[Fitzpatrick's plasma-frequency derivation](https://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node6.html)
gives each species contribution `ω_ps² = n_s q_s²/(ε0 m_s)`. For the cold collective
mode used in the tests, linearizing continuity and acceleration with a common
longitudinal field gives `ω² = Σ_s ω_ps²`. Here `ε0=1`. Pair-plasma and unequal
charge/mass tests compare this prediction and timestep refinement, and separately
check momentum and centered-energy convergence. This is an explicit leapfrog
integrator; exact discrete energy conservation is not claimed. Reported ordinary
energy uses stored half-time velocities and full-time fields.

Initialization binds ordered population identities. Replacement populations must
fit the constructor vector's element type. A new state is required to change a
concretely typed storage representation. The documented scalar/count constructors are unchanged; the concrete state
type gains a workspace parameter. Homogeneous strided views and explicitly
heterogeneous typed vectors are supported. Restored values may come from a
different supported source representation; destination storage is independently
validated before publication. Initialization resets phase without converting
velocities, so callers must supply physical full-time velocities when restarting
the phase manually. After manually changing positions, weights or charges, or
adding particles, call `poisson_E!` to refresh the self-field while preserving
phase. The step uses the supplied common old field; it does not silently replace
an imposed field. New particles still require velocities at the existing half
level. Population replacement/reordering follows the initialization contract.

## Progress

- [x] Map source, callers, diagnostics, checkpoint/archive and existing tests.
- [x] Record design and compatibility decisions before production edits.
- [x] Implement collection construction, validation, charge and synchronized step.
- [x] Implement shared-field diagnostics and collection checkpoint/archive.
- [x] Verify functional, physical, failure and compatibility cases on both Julia versions.
- [x] Format, static/architecture review and final source freeze.
- [x] Verify scalar bytes/allocations and collection storage/growth on both runtimes.
- [x] Preflight complete timing harness, including output retention on failure.
- [x] Revalidate V7 boundary-scan cleanup and freeze the candidate: all 14 final command/log receipts pass, as recorded below.
- [ ] Complete coordinated timing acceptance.
- [x] Root integration and main push as `d46029e`; original failed CI and containing `4724947` six-job/compiler acceptance are distinguished above.


## Audit matrix

| Surface | Required evidence | Status |
|---|---|---|
| Charge/Poisson | q-weighted reference, mean policy, permutations | Verified on Julia 1.12.7/1.13.2, including three active charges |
| Time integration | shared old field, variable dt, cold-mode refinement | Pair/unequal-q/m dispersion, energy refinement and growth phase pass |
| Ownership | aliases, invalid last species, phase binding | Atomic failure, concrete/heterogeneous views and owner-type changes pass |
| Diagnostics | species totals and field energy once, stated time levels | Exact charge, grid-integral range and field-energy range regressions pass |
| Persistence | schema validation, no partial restore, exact continuation | Archive, variable dt, malformed schema and Array/view transfers pass |
| Compatibility | existing scalar tests and shared kernels | Scalar tests pass; 96 scalar byte/allocation rows match parent per runtime; neighboring EM owners unchanged from passed dual-version tests |
| Resources | reusable workspace shape, measured full-call costs | 144 storage/growth rows pass per runtime; timing acceptance pending |

Plausible wrong patches and rejecting tests: electron-only charge summation is
rejected by unequal-charge tests; sequential field updates between species are
rejected by permutation and shared-force tests; field energy counted per species
is rejected by diagnostics decomposition; resetting the half-level on restart
is rejected by variable-step exact continuation; validating inside the species
update loop is rejected by invalid-final-species snapshots.

## Dependent CAP-008 work

EM arbitrary species/current summation, multirate and relativistic scheduling,
local Coulomb composition, bounded domains, distributed/device ownership and
portable cross-version archives remain separate, open implementation tracks.
This slice uses CPU periodic CIC and synchronized nonrelativistic velocities.


## Verification ledger

Evidence is retained outside the repository in the audit workspace under
`electrostatic-species/`; commands record the worktree/runtime and preserve
failed logs. No performance result is inferred from functional test durations.

- `collection-v1-1132.log`: 832 assertions passed on Julia 1.13.2.
- `scalar-regression-v1-1132.log`: existing ES and scalar checkpoint suites
  completed with exit 0. This precedes the final field-energy range repair.
- `collection-charge-v2-1132.log`: collection plus existing EM complete-charge
  combination/conservation suites completed with exit 0. This precedes the
  exact-neutrality owner move and final diagnostic additions.
- `field-energy-range-before.log`: verified finite-result overflow in the old
  ES `field_energy` for both Float32 and Float64. The new implementation reuses
  the existing field-energy range owner; final dual-version regressions pass.
- Formatter v1 retained 4 failed AST comparisons: its line-node removal did not
  recursively remove macro source locations. No functional gate was weakened.
  Corrected recursive normalization is used by formatter v2/v3; v3 passed 22
  syntax/idempotence checks. The final diff still receives manual review.
- `source-v3.json` and `commands-v3.json`: frozen final-candidate source hashes
  and sequential Julia 1.12.7/1.13.2 functional, Aqua and compiler checks.
  Both runs terminated; the failed ES view fixture and compiler findings are
  recorded below. Neighboring EM testsets passed, with unchanged-owner hashes
  recorded in `retained-neighbor-evidence-v6.json`.

The physical charge check now uses the shared exact dyadic-sum owner. It detects
nonzero charge whose products cancel only after floating-point rounding; an
explicit zero neutrality tolerance therefore has exact represented-input meaning.
Diagnostics retain exact rational charge totals. Density/Poisson and particle
updates still use their declared floating-point precision.


### Independent v3 findings (not accepted)

- The public collection `poisson_E!` source-mean check rejected a previously
  initialized neutral pair plasma: exact physical net charge 0, deposited
  density sum -1/281474976710656. Its scale was the small cancelled grid source.
  The staged `init_espic!`/`step_espic!` path passed the physical policy correctly.
  Repair adopted: public collection Poisson refreshes its bound populations
  through the same physical validation and staged field owner, preserving phase.
- Deposited-charge diagnostics returned -Inf for four finite cancelling cells
  whose exact integral is 0. Repair adopted: retain the complete grid sum and
  spacings through the existing exact owner.
- Julia 1.12.7 v3: 16,857 assertions passed, 0 failed, 1 fixture error. A shortened
  Vector record could not be assigned into the fixture's concrete vector of
  view-backed records. The production restore was never reached by that case.
  The correction uses a heterogeneous record container and shortens both saved
  populations, preserving neutrality so the restore-storage gate is exercised.
- Concrete JET 18/18 and Aqua 9/9 passed. Package JET reported three added paths
  (32 total versus the reviewed local 29 baseline), failing the unchanged ceiling:
  longitudinal field tuple indexing and scalar/collection schema cross-branches.
  Fixes make the tuple embedding and payload-layout dispatch explicit; no
  compiler suppression or ceiling change is planned.

The frozen v3 run and failed probes are retained. Source files are never replaced
beneath a running verification snapshot.

### v4 revalidation and v5 follow-up

Both v4 functional runs completed with 6,140 assertions passed and four failed.
All four failures were in the former floating `deposited_charge ≈ 4` fixture:
exact rational diagnostics use an exact default comparison against an integer,
whereas the stored CIC grid has ordinary floating deposition roundoff. The v5
oracle checks the exact represented grid integral independently and separately
bounds deposition error in the solver precision. The public refresh, field-range,
phase, malformed-schema and view-backed continuation regressions passed.

Concrete v4 compiler checks passed 24/30 on each runtime. The six failures traced
to public refresh reading owners through an erased `ParticleSet{D,T}` element
type. V5 retains the constructor collection element type in the state/workspace
and validates replacement/restore owner types before mutation. It does not impose
destination representation on saved source arrays. A concrete Julia 1.13.2 public
refresh probe now reports no JET errors; the frozen v5 dual-version matrix
subsequently completed. It adds three active charged populations with unequal sizes, direct
shared-force checks, Array/view checkpoint transfers, and concrete view compiler
checks. Formatting passed 24 syntax/idempotence assertions before the v5 freeze.


### Final independent review

The v5 minimum-runtime functional matrix reports 6,376 passes and 18 failures in
cross-storage comparisons. Array and view checkpoints necessarily serialize
different type/parent graphs; the fixture incorrectly used byte snapshots for
that comparison. The corrected fixture compares every physical field, phase and
particle record directly. Byte snapshots remain the failure-before-publication
oracle when the storage representation is unchanged. Both v5 quality runs pass
48 concrete compiler assertions, Aqua, and the unchanged package gate (29 and
27 reports); full report-tree comparison follows the final freeze.

Review also clarified the manual reinitialization contract: initialization resets
the phase and does not convert half-level velocities. New tests cover mixed
Array/view collections, replacing members within an explicitly abstract element
type, growth at the existing half level, and reuse of stage/mesh objects. The subsequent untimed resource campaign measures complete-step allocations and
retained workspace identity. Logical array lengths alone do not establish a heap
capacity bound; runtime timing acceptance remains open.


## V6 verified evidence and remaining acceptance

`source-v6.json` and `commands-v6.json` preserve the V6 Julia source hashes,
commands, process IDs, exit codes and log hashes. Every command exited 0:

| Gate | Julia 1.12.7 | Julia 1.13.2 |
|---|---:|---:|
| New collection, diversity, scalar ES and ES checkpoint assertions | 6,496 | 6,497 |
| Concrete Array/view compiler assertions | 48 | 48 |
| Aqua assertions | 9 | 9 |
| Raw package/loader compiler assertions | 4 | 4 |
| Scalar complete-step byte/allocation assertions | 216 | 216 |
| Collection storage/growth assertions | 720 | 720 |

The public pair-plasma example also executes on both runtimes. Formatting passes
24 AST/idempotence assertions. Complete package compiler trees match the parent
at 29/27 reports after checked source-line mapping, module/gensym normalization,
and removal of one explicitly retained interleaved precompile progress burst.
No report, signature or message is omitted; no source frame is unmapped and no
compiler ceiling changes. See `quality-v6-*-trees.json` and the empty tree diffs.

`resource-gates-v1-*.toml` contains 96 scalar and 144 collection rows per runtime.
Scalar particles, E/ne/Ep and phase match the parent byte-for-byte after complete
calls, including variable/zero dt; both versions allocate 0 bytes in those rows.
This compatibility gate does not compare the intentionally repaired field-energy
rounding/range behavior. Collection rows cover S=2/4/8, balance/imbalance and
warm/grow/shrink/empty-population transitions. They retain stage and mesh objects,
validate logical lengths, and report allocations (10,360–416,544 bytes across the
measured rows, including growth). Reachable Julia heap measurements include owner
references to caller populations and exclude native FFT plan storage.

The prepared cost harness has 150 workload pairs and 30 representative control
fixtures per runtime. It compares original scalar/candidate scalar and
scalar/one-population physical equivalence, then collection species-count and
imbalance costs, initialization, public refresh and diagnostics. Default mode is
untimed preflight. Timing mode uses repeated, order-alternated pairs plus identical
and actual doubled-call controls; root must coordinate an exclusive CPU window.
Raw wall time from tests is not performance evidence. No full optimization,
throughput, cross-backend or full CAP-008 completion claim is made here.


## V7 boundary ownership review

`apply_periodic!` in `src/boundaries/periodic.jl` validates all positions before
wrapping any coordinate. The shared ES pusher previously scanned those same
proposed positions immediately before calling it. V7 removes that duplicate
post-push scan, retains velocity validation and retains the separate collection
input preflight. This is a source-backed redundant-work removal; no measured
speedup or global optimization claim follows from it.

V6 sources are archived unchanged. `duplicate-position-scan-v7.diff` and its JSON
record retain the exact one-function delta and the canonical boundary source
hash. V7 proportional verification passed on both runtimes: the full ES functional matrix,
an independent positive/negative-stride checkpoint probe, scalar byte/allocation
and collection storage gates, 48 concrete compiler calls, and the unchanged cost
harness preflight on both runtimes. Whole-package/Aqua and unchanged neighboring
EM evidence remain bound to their V6/V3 snapshots.

The frozen V2 cost harness passed 1,464 preflight checks per runtime on V6. It
writes full raw rows to stdout and preserves TOML output in `finally` while
propagating failures. A deliberate one-assertion failure exited nonzero and left
`completed=false` plus the raw probe row (`failure-retention-v2.json`). No timing
acceptance has run. The final candidate must retain its source and harness hashes
through the separately coordinated timing window.


The final V7 results on Julia 1.12.7 / 1.13.2 pass: 6,496 / 6,497 functional assertions plus the public example, 120 independent strided-storage assertions per runtime, 216 scalar and 720 collection resource assertions, 48 concrete compiler calls, and 1,464 V2 preflight checks. Commands and hashes are retained in `commands-v7.json`. The stride probe uses positive and negative stride-two views, unequal charge/mass, variable and zero steps, public checkpoint restore, and unchanged sentinels outside each view.

Final harness review found that V2 covered fixed total particle count but omitted the resource plan's fixed particles per species sweep and general nonuniform-weight path. V3 adds those cases and explicit single-thread Julia/FFTW settings, retaining V2's unchanged sampling, controls, acceptance thresholds and failure-output path. Both runtimes pass 1,908 untimed preflight assertions; all 198 workload IDs match the independent manifest in `cost-v3-exact-coverage.json`. Timing mode has 72 additional control rows per runtime. Actual Array backing-capacity probes pass 864 assertions / 144 rows per runtime after growth, shrink and empty transitions; their internal storage access is pinned to both Julia runtime sources and confined to audit tooling. `commands-v7-supplement.json` retains command and log hashes. Process maximum RSS includes compilation and all prior fixtures; it is not a per-operation peak-memory measurement or native FFT bound.

## Pre-existing particle-kinematics range defect

An independent audit found that finite q/m and dt can overflow the shared Boris
coefficient `q/m*dt/2` and produce NaNs on an exactly neutral zero-field
equilibrium. The defect lives in the shared pusher and priming arithmetic used
by every solver, not in the collection code, so it is repaired separately under
the normwise particle-kinematics plan rather than blocking this slice.
