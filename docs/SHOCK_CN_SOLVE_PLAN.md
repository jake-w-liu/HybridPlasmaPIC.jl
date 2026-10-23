# Shock Crank–Nicolson solve repair

Owner: `/root/completeness_architecture`; integration and exact-SHA CI: `/root`.
Baseline: `df75656390ed416f7b0a37cb1173058cd9465e71`. This is CAP-005A:
repair the existing reduced 3-D shock induction solve. It does not complete
CAP-005's nonlinear particle–field coupling or add Hall induction to this model.

## Progress and evidence

| Stage | State | Evidence / next gate |
| --- | --- | --- |
| Map canonical solver, callers, persistence and dependencies | Checked | Source map below; no existing general nonsymmetric matrix-free solver dependency |
| Reproduce against independently assembled trapezoid | Confirmed defect | 28 profiles per native lane; both lanes: 66 pass, 50 fail, 0 errors |
| Review repair, compatibility and transaction boundary | Reviewed | Parent accepted maintained solver, full-step transaction and Vector-only commit scope |
| Implement residual-controlled canonical solve | Frozen candidate6 | Matrix-free GMRES, initial/final true residual and full-step transaction; all confirmed candidate gaps repaired |
| CRC numerical, failure, caller and restart verification | Pass | Candidate6: 666/666 broad and 24/24 independent public-boundary checks on each native lane |
| Measure storage, allocation, ordinary and error-matched costs | Pass with recorded repair costs | First fixed matrix and one prescribed four-profile confirmation complete; all raw samples retained |
| Matched fresh-process compiler comparison | Pass with retained abstract diagnostics | 28 full records match on minimum plus two removed/one reviewed generic KdV report; all 27 match on latest; concrete KdV/FFT differentials pass |
| Independent final audit and exact final diff | Local and root source review pass | Boundary audit passes; root reviewed six source/test diffs and full owner/tests without a new confirmed defect; final composition is verified; exact-SHA CI remains |
| Push main and exact-SHA GitHub CI | Pushed ce52836; CI running | All 363 committed Julia files match the candidate; all six exact-SHA jobs required |

Durable evidence directory, outside the source worktree:
`HybridPlasmaPIC-audit-evidence/2026-10-23/shock-cn-solve/`.
`baseline-source.json` records 279 source/test/project file identities;
`baseline-commands.json` records exact commands, versions, process ownership,
exit codes and hashes. `baseline-witness.jl` constructs SBP and Fourier
derivative matrices independently of production FFT routines. Its two native
logs and TOML metric files retain successful controls and failed expectations.
The baseline source identities were checked before and after both runs.
An initial multi-file patch failed its exact-context check and applied no changes;
the corrected patch used the inspected source context.

## Pre-edit brief

- **Mode:** confirmed numerical repair with a new explicit failure boundary.
- **Requirement:** a successful `field_method=:cn` step must solve its frozen
  moment trapezoid equation to a stated residual tolerance; failure must leave
  the entire caller-owned shock and particle state unchanged.
- **Baseline:** native Julia 1.12.7 and 1.13.2 each reproduce 50 failures in 116
  assertions. All independent-RHS and direct-solve residual controls pass.
  For Float64 public stepping with zero particles, pure transverse resistive
  decay, `dt=4`, `NB=1`, the current field reaches `max(abs(B))=3.125`; the
  independently assembled CN solution reaches approximately `0.533709` and
  the current normalized equation residual is approximately 8. At `dt=10000`,
  the current maximum is approximately `1.56219e14`, versus approximately
  `0.601087` for the direct reference. Exact lane values remain in the TOML.
- **Verified cause:** `_cn_b3d!` accepts four fixed Picard iterates without a
  convergence check. For scalar `z=h*lambda`, those four iterates have
  amplification `1+z+z^2/2+z^3/4+z^4/8`; this differs from the trapezoid
  amplification `(1+z/2)/(1-z/2)`. The negative-real diffusion witness above
  exercises this discrepancy through the public caller.
- **Contract:** retain the existing spatial operator, frozen moments, SAT
  target, particle push, boundaries and subcycling. Replace the inaccurate
  algebraic solve. A true residual, finite result and bounded iteration count
  determine acceptance. Solver failure is observable and never silently
  replaced by explicit stepping or a partial iterate.
- **Canonical owner:** `src/initial_conditions/reflecting_wall_shock3d.jl`;
  `_b_rhs3d!`, `_cn_b3d!`, `step_shock3d!`, `run_perp_shock3d`. Affected
  surfaces: dependency/import/export declarations, focused tests and docs.
- **Compatibility:** `:rk4` and driver defaults stay unchanged. No field is
  added to `PerpShock3D`; checkpoint serialization retains its existing shape.
  An optional reusable CN workspace and explicit solver controls are added.
  CN now rejects view/custom particle storage, including empty populations;
  this deliberate storage restriction makes nonthrowing commit possible.
  The ordinary RK4 storage dispatch is unchanged.
  Correct CN results intentionally differ from the erroneous four-iterate
  approximation. The internal `iters` keyword is not a public convergence API.
- **Strategy:** maintained restarted GMRES applied to the existing induction
  operator; explicit trial state stages the entire CN public step; commit only
  after every magnetic substep and final electric field succeed.
- **Verification:** independent assembled solve, analytic mode amplification,
  refinement, full-step rollback, exact restart, RK4 differential state,
  package/compiler checks on both native versions, then coordinated timings.

## Source map and bounded physics scope

Line references in this section refer to the frozen baseline.

| Surface | Checked source | Consequence |
| --- | --- | --- |
| Reduced induction | `reflecting_wall_shock3d.jl:1–20,337–359` | Preserve `E_ind=-u×B+eta*curl(B)`; full particle Ohm law is separate |
| State and constructor | `reflecting_wall_shock3d.jl:25–61,76–148` | CPU arrays and scalar parameters; no solver workspace field |
| Spatial derivatives | `reflecting_wall_shock3d.jl:150–209` | SBP x, Fourier y/z, zero even Nyquist derivative |
| Moment update | `reflecting_wall_shock3d.jl:234–275` | Post-push moments are frozen during all magnetic substeps |
| Particle electric field | `reflecting_wall_shock3d.jl:305–332` | Final E recomputation is part of transaction success |
| Existing CN approximation | `reflecting_wall_shock3d.jl:384–404` | Four Picard updates; no residual or failure result |
| Public step | `reflecting_wall_shock3d.jl:423–505` | Particle push, absorption/compaction and moments precede field advancement |
| Driver defaults | `reflecting_wall_shock3d.jl:638–671` | `dt=.02`, `field_method=:rk4`, `NB=2` in loop |
| Raw restart | `verification/shock_campaign.jl:12–34,62–69` | Serializes `(sh,ps)` directly; continuation compares every state field |
| Existing CN regression | `test/test_phase_final.jl:260–281` | Small-step macroscopic shock comparison does not establish linear solve accuracy |

The source filename in the first seven rows is under
`src/initial_conditions/`. The restart source is under `src/`.
Existing general `recommended_dt` guidance belongs to the periodic hybrid
integrator; this public shock driver does not call it. A converged CN magnetic
substep does not eliminate explicit particle-push, coupling or resolution
accuracy requirements. The full SBP/SAT system is not asserted to be globally
energy stable merely because the scalar trapezoid method is A-stable.

## Solver and failure contract

At fixed deposited velocity, write the existing induction RHS as `f(B)=L*B+b`,
where `b` is the upstream SAT target contribution. Solve

```text
M = I - h*L/2
rhs = (I + h*L/2)*B_old + h*b
M*B_new = rhs
```

The homogeneous matrix-vector product must call the same induction owner with
zero SAT target. No second curl, resistive operator or SAT convention is
implemented inside a solver adapter. The small assembled matrix is an
independent test oracle only; production must not allocate an O(grid^2) matrix.

User controls are `cn_rtol=100eps(T)`, `cn_atol=0`,
`cn_maxiter=200`, and optional `cn_workspace`. `ShockCNWorkspace(sh,ps;
memory=20)` owns bounded restarted Krylov storage and trial state. These
defaults pass the independent equation, analytic-mode and refinement checks.
Tolerance changes must be justified against the fixed accuracy tests, never
loosened to hide a failed solve. Tolerances must be finite and nonnegative, relative
tolerance below one, at least one tolerance positive, and iteration/restart
limits positive integers.

After GMRES, recompute the true unpreconditioned residual using the canonical
operator. Success requires finite input/output/residual working values and

```text
norm(M*B_new - rhs, 2) <= cn_atol + cn_rtol*norm(rhs, 2).
```

The library's estimated residual alone is not acceptance evidence. A typed
convergence exception should expose the true residual, acceptance threshold
and iteration count; nonfinite arithmetic must reject the step. This does not
promise arbitrary-range arithmetic or convergence for singular/ill-conditioned
operators. An admissible solve either meets the declared residual or fails
without a public partial update.

The transaction starts before the first particle mutation. The explicit
workspace owns a separate trial shock and trial particle storage. Canonical
push/deposit/subcycle/E operations run on the trial. Successful commit copies
into the caller's existing arrays, preserving their identities and accepted
absorption/compaction semantics. There is no exception-to-success catch.
On rejection, caller fields, scratch, moments, particle values, IDs/tags and
lengths are unchanged; workspace scratch may change and can be reused.
Storage must be validated before mutation or commit; unsupported aliases or
incompatible workspace dimensions must fail before touching caller state.

Workspace ownership is explicit and single-task; concurrent use of the same
workspace is unsupported. No global cache or hidden warm-start history is
introduced. A fresh workspace reconstructed after checkpoint restore must
produce identical continuation to a reused workspace. The driver creates
one reusable workspace outside its loop. Construction, steady-state cost and
storage growth with changing particle counts are measured separately.

## Dependency decision

`Project.toml` and source searches contain no existing generic nonsymmetric
matrix-free Krylov solver. `LinearAlgebra` provides the independent dense
oracle but does not supply the proposed production iterative solve. More
Picard iterations do not provide convergence for the reproduced stiff case.

Chosen dependency: Krylov.jl 0.10.11, with a compatible 0.10 version bound.
The initial source review used 0.10.10; resolution selected 0.10.11. The installed
project/license and complete GMRES change were checked before testing. The
dependency tree is `e0574e7e2cd5c40d9a6ae444511fbadd51afbf30`;
`dependency-resolved-source.json` records actual installed source hashes.
The new optional weighted inner product is unused (`W=I`); the updated
restart normalizes by the newly recomputed residual. The minimum bound names
the version actually tested, rather than claiming unrun 0.10.10 support.
Its pinned [project declaration](https://raw.githubusercontent.com/JuliaSmoothOptimizers/Krylov.jl/v0.10.10/Project.toml)
lists only LinearAlgebra, Printf and SparseArrays runtime dependencies and
Julia 1.10 compatibility. Its [license](https://raw.githubusercontent.com/JuliaSmoothOptimizers/Krylov.jl/v0.10.10/LICENSE.md)
declares MPL 2.0. No dependency source is copied into this repository.
Official [in-place workspace documentation](https://jso.dev/Krylov.jl/stable/inplace/)
and [matrix-free interface](https://jso.dev/Krylov.jl/stable/matrix_free/)
support reusable storage and a focused `mul!`/`size`/`eltype` operator adapter.
Pinned source snapshots and SHA256 hashes are retained in
`dependency-review.json`. Runtime, inference, load/storage and cost
measurements are recorded below. No security or performance guarantee is inferred
from a dependency's existence.

## Acceptance matrix

| Dimension | Required check | State |
| --- | --- | --- |
| Equation correctness | All 28 frozen baseline profiles pass direct CN and true residual checks on both precisions/versions | Pass — candidate6 on both native Julia lanes |
| Independent physics | Fourier advection has unit CN amplification magnitude; resistive transverse energy follows analytic multiplier; SAT-only affine relaxation reaches the declared target equation | Pass — candidate6 on both native Julia lanes |
| Accuracy | Small-step agreement plus second-order magnetic time refinement against matrix exponential/analytic solution at fixed moments | Pass — candidate6 on both native Julia lanes |
| Spatial conventions | Odd/even periodic axes, singleton transverse axes, nonuniform mixed fields, nonzero SAT forcing and off-diagonal component coupling | Pass — candidate6 on both native Julia lanes |
| Public composition | Initialized nonempty particles, changing moments, NB changes, particle reflection/absorption, and final E from accepted fields | Pass — candidate6 on both native Julia lanes |
| Failure | Forced iteration exhaustion, malformed controls/workspace, nonfinite arithmetic, later substep rejection, and finite final E rejection preserve all caller state | Pass — candidate6 on both native Julia lanes |
| State ownership | Caller array identities retained, trial aliases rejected, independent workspaces do not interfere, failed workspace reusable | Pass — candidate6 on both native Julia lanes |
| Restart | Checkpoint before/after CN, fresh versus reused workspace, varying dt/NB, full-state bitwise continuation | Pass — candidate6 on both native Julia lanes |
| Compatibility | Baseline RK4 full-state differential and existing 3-D shock/campaign tests, including invalid-input order | Pass — candidate6 on both native Julia lanes |
| Resources | Bounded Krylov basis, storage formula and measured summary size, repeated-step allocation, particle shrink/growth reuse | Pass — stored workspace layout unchanged; candidate6 warmed allocation/growth and RK4 differential gates |
| Static/integration | Concrete solver/caller JET, package JET baseline comparison, Aqua, normal test entry point, latest/minimum Julia | Pass — candidate6 on both native Julia lanes |
| Performance | Exclusive interleaved repeated-workspace ordinary profiles and separately stiff/error-matched profiles; frozen settings and raw samples | Pass — fixed matrix and one prescribed reverse-order confirmation; measured CN repair costs retained |
| Final independent audit | Review exact final diff through analytic-mode/public-restart route, including untested assumptions and all failure states | Pass — separate 24-check public boundary audit and root exact source/test diff review |

Wrong patches rejected by these gates include: merely increasing Picard
iterations (stiff direct oracle), solving the homogeneous problem but dropping
SAT forcing (affine equilibrium), trusting an estimated GMRES residual (true
residual check), committing each substep before all succeed (later-failure
rollback), and retaining workspace history that changes restart (fresh/reused
bitwise continuation). These are test obligations, not claims already verified.

Optimization acceptance will be fixed before timing: RK4 numerical states must
remain exact; solver acceptance and error targets stay identical across
candidates; steady-state storage is O(grid*restart + particles), with bounded
restart independent of iteration cap. A converged method cannot fairly claim
a speedup by comparing its cost against the baseline's inaccurate stiff step.
Report ordinary correctness-repair overhead separately from equal-error work.


## Numerical-reference correction and first acceptance

The original failing baseline and its original Float64 dense oracle remain
unchanged in `baseline-witness.jl`. Candidate1 passed 291 latest-lane checks,
including analytic amplification, full public rollback and exact restart, but
failed the stiff Float64 dense forward-error comparison: `3.48e-10` on 1.12
and `4.23e-9` on 1.13.2 against the unchanged `1000eps(Float64)` gate.

The isolated, single-BLAS-thread `oracle-conditioning.jl` checks both analytic
transverse diffusion and a separately assembled BigFloat matrix. On both
native versions, candidate error against the analytic solution and independent
BigFloat solution is `1.11e-16`; the BigFloat reference's analytic error is
`5.55e-17`. The Float64 matrix has condition number about `3.824e8` and a
constant-field operator defect `4.44e-16`; its Fourier derivative row-sum
error is `2.22e-16`. Independent BigFloat assembly reduces those defects to
about `5.18e-77` and `2.59e-77`. Merely solving the already-rounded Float64
matrix in BigFloat still has `5.39e-10` forward error. The accepted candidate
retains the constant mean within `5.55e-17`.

Accordingly, the stiff test now independently assembles its reference in
BigFloat, keeping the same forward-error threshold. No production branch
recognizes this test. Candidate2 passes 308/308 focused assertions on each
native version. The corrected oracle run against a separate, unmodified
`df75656` worktree retains six successful controls and four expected failures,
so it still rejects the original Picard implementation. Old test source and
all failed/passing logs are retained; this is an oracle correction justified
by independent controls, not a tolerance relaxation.

All 48 RK4 snapshots per version serialize identically between baseline and
candidate, including complete shock scratch and particle state. Both versions
produce SHA256 `89d09e485f113dd66ab8c12404526fae250205013520db34e06c52a271385c7a`.
This verifies numerical compatibility for that matrix; timing compatibility
is a separate pending measurement.

Concrete JET identified two paths in the stdlib generic norm's recursive
container check for each CN signature. The explicitly supported CPU
Float32/64 vectors can use the public scaled `BLAS.nrm2` primitive directly.
Before changing that call, `norm-probe.jl` passed 22 checks per version: empty,
zero, ordinary, finite extreme, overflowing finite-component norm, subnormal,
NaN and infinity cases, plus zero concrete JET reports for the primitive.
Both residual norms now use that primitive and retain explicit finite norm
and finite threshold checks. Additional public tests reject both an
unrepresentable norm with individually finite components and an overflowing
tolerance. No inference warning was suppressed.

One compiler harness attempt placed `const` inside a testset's local scope;
it failed before acceptance and was corrected by wrapping `include` from the
outer driver. Its log is retained as a harness error, not a package failure.
Candidate3 reruns focused checks, existing 3-D shock/final-phase regressions,
concrete compilation, Aqua and the unmodified package JET gate on both native
versions. This section records progress; it is not final completion evidence.


## Independent ownership audit and final gates

Candidate4 adds external particle growth/reuse and a later-substep overflow
witness. For admitted negative SAT damping with scalar lambda=1 and h=1, CN
amplifies by three. An initial upstream component floatmax(T)/8 therefore has
a representable first substep and an unrepresentable second. The measured
Float32/64 public witness throws, retains exactly the first accepted substep
in trial storage, and leaves the full caller state unchanged.

The subsequent dependency ownership audit found a candidate transaction defect:
Krylov 0.10.11 `reset!(SimpleStats)` empties `residuals`, `Aresiduals`, `qvals`
and `Acond`, even though GMRES history collection is disabled. Assigning the
first vector to caller `sh.x` reproduced a BoundsError and emptied caller
geometry (`stats-alias-probe.jl`, `.log`). All four arrays now participate in
validation, including workspace cross-aliases and identical empty arrays; all
40 alias configurations and unchanged-caller assertions pass in candidate5. This
candidate failure is retained rather than hidden by the earlier passing suite.

The integer-control probe also confirmed that NB=big(2) promotes h to BigFloat
and fails the concrete Float32/64 CN method, whereas Int8/Int128(2) succeed.
The CN boundary now uses the existing positive-count normalizer (fits Int)
before staging; exact Int8/Int128/BigInt differentials and overflow rejection pass; RK4 dispatch remains unchanged. Nondefault CN controls and any
nonnull CN workspace are rejected with RK4. Explicitly supplied default control
values are accepted and unused in RK4, since Julia keyword values do not retain
whether the caller typed them.

The final frozen revision must rerun focused physics/failure/restart checks,
concrete compiler signatures and package quality gates after these amendments.
The earlier fresh minimum-lane package JET comparison has baseline 30
and guarded candidate 29 reports. The complete reports are retained; their
unrelated AMR/KdV paths are not byte-identical, so counts are not treated as proof
of identical findings. Direct inspection finds no remaining shock CN/Krylov
report after enforcing concrete Float32/64 precision at the operator boundary.
No JET suppression or package ceiling change is used.

The existing Fourier derivative creates fresh buffers/plans. In the measured
1024-cell profile, Float64 allocates 17200 bytes per derivative, 137600 per RHS,
and 963200 per private solve; Float32 allocates 9008, 72064 and 288256 respectively.
The exact multiples attribute that measured solve allocation to its RHS calls.
State copy measures zero bytes; validation measures 3552 in candidate3. These
are allocation observations, not runtime speed claims. FFT workspace caching
is a separately reviewed follow-up after this repair is frozen, so it cannot
hold the confirmed numerical correction or bypass its unchanged accuracy gates.


### Stationary-state dependency defect and accepted repair

The independent boundary audit after candidate5's 602/602 passes found a valid
stationary nonzero field rejected with true residual 1 against tolerance
1.1920929e-5 (Float32). The caller was unchanged. Installed Krylov 0.10.11
`gmres.jl:165` adds the explicit initial guess when restart is enabled; its
zero-initial-residual branch adds that guess again at line 182. The source
witness and failed runtime probe (`final-contract-audit-v1.jl`,
`final-contract-audit-1132.log`) are retained. This is an observed failure in
a valid equilibrium, not a reason to weaken the true residual check.

The reviewed repair checks the initial guess against the same canonical
operator and finite true-residual threshold before GMRES. If it already
satisfies that equation, B remains unchanged. Otherwise the maintained solver
and the final true-residual gate run as before. No exception is swallowed and
no explicit-method fallback is added. Exact stationary/nonzero-field cases,
restart memory 1 and 3, zero/tiny/ordinary intervals, and nonstationary →
stationary → nonstationary reuse are added to permanent coverage. The extra
operator evaluation on nontrivial solves must be included in final costs.
Workspace solver statistics describe its most recent library call; they are
scratch and are not a separate public CN result contract.

Candidate5's untimed cost harness passes all 18 checks on each baseline/candidate
and Julia-version combination. Eight RK4 allocation profiles per lane match
baseline exactly, and all 48 complete-state snapshots per lane remain
byte-identical. These checks precede the stationary repair; source hashes and
revision snapshots preserve their exact scope. Candidate5's equal-error
transverse-mode absolute error is about 1.19235e-6 versus a fixed 5e-6 gate;
the baseline error is about 1.18875e-6. No elapsed-time samples have been taken.
The stationary repair must rerun these gates before final reporting.


The final checkpoint test also saves immediately before the absorbing step,
then compares uninterrupted/reused-workspace and restored/fresh-workspace
execution through that same particle removal. The existing checkpoint after
absorption and subsequent varying-dt/NB continuation remains. This explicitly
covers absorption during restored execution rather than inferring it from
an already-compacted checkpoint.


## Frozen candidate6 evidence

The final source identity manifest is `candidate6-source.json`; all six changed
source/test/project files matched before and after the native correctness runs.
`candidate6-code.patch` applies cleanly to the original baseline. The completed cost matrix and compiler classification are recorded below.
The final composition also updates the README, changelog and master trackers.

- Julia 1.12.7 and 1.13.2: 666/666 broad assertions each, including the complete
  permanent CN suite, existing 3-D shock/final-phase suites, eight concrete
  JET signatures, Aqua and the unchanged package JET gate.
- Independent public-boundary audit: 24/24 each, including singular matrix
  rejection, nonzero stationary fields with restart memory 1, invalid memory,
  nonfinite inputs with byte-preserved caller state and driver option handling.
- Full package compiler logs retain 29/27 abstract reports. Direct path review
  (`candidate6-compiler-path-review.json`) finds no shock CN, shock3d or Krylov
  path among those reports; every concrete signature has zero reports. These
  findings do not reclassify unrelated package reports as harmless.
- Formatting/parsed-semantics gate: 10/10; code patch/application and whitespace
  checks pass. No assertion, threshold, JET ceiling or dependency warning is
  suppressed.
- Four untimed baseline/candidate cost preflights: 18/18 each. All eight warmed
  RK4 allocation profiles and all 48 serialized state snapshots per Julia lane
  match baseline exactly. CN ordinary batch allocations are 1.015–1.842× the
  old method on 1.12 and 1.015–1.398× on 1.13.2; this includes the initial true
  residual application and the full transaction. Equal-error allocation is
 1.0085× in both lanes. These are measured repair costs, not a speedup claim.

### Final matched compiler and elapsed gates

The first fixed four-process cost matrix completed with 18/18 numerical gates
in every process. Four RK4 median ratios exceed the predeclared 1.10 flag:

| Julia | Precision / grid / particles | Candidate / baseline median |
| --- | --- | --- |
| 1.12.7 | Float32 / 3×5×3 / 64 | 1.1417 |
| 1.12.7 | Float64 / 8×8×4 / 64 | 1.1166 |
| 1.13.2 | Float64 / 3×5×3 / 64 | 1.2999 |
| 1.13.2 | Float64 / 8×8×4 / 64 | 1.1371 |

The one prescribed reverse-order confirmation is complete: ratios are
0.969113 / 0.883692 / 1.087395 / 1.030978 in the table order, all within 1.10;
allocations match for every profile. The driver preserves the original
profile-function bytes, nine samples, ten-step batches and three warm batches.
No threshold or sample count changed, and no further retry was run. Other task
native work was held during sampling. All first-run samples and comparisons are in
`candidate6-first-timing-comparison.json` and the four raw timing TOML/log pairs.

The baseline package-only compiler driver and candidate after broad tests have
different entry contexts. Their 30 versus 29/27 counts and absent CN paths do
not complete compiler acceptance. Static full-tree classification retains 22
identical normalized chains, six same-location/same-message chains with display
differences, two removed generic leaves, and one newly displayed KdV/FFTW tuple
leaf in the minimum lane. All 15 referenced repository source files are
byte-identical, which establishes source ownership but not runtime harmlessness.

The required replacement comparison uses one identical isolated driver in four
fresh processes, retaining every frame, full specialization/signature and leaf
message without depth limits. Within each runtime, resolved Project/Manifest
graphs match except the package worktree path. Both resolved test environments
include Krylov; the baseline production Project and source remain unchanged and
do not declare/import it. If the newly displayed KdV branch persists, concrete
compiler and runtime differentials in both worktrees must settle that branch.
`matched-environment-manifest.json` records the exact dependency identities;
`matched-package-compiler.jl` and `run-matched-final-gate.py` have now run in
all four fresh processes. Their full raw-record comparison matches 28 records
on minimum Julia, with two removed generic records and one added KdV/FFTW
record; latest matches all 27. Exact method-site checks justify two generated
keyword-function number mappings; no diagnostic data are discarded. The
original differing trees remain retained.

The added minimum-lane generic report loses the FFT region type before an
abstract sparse-array iterator branch. Source inspection verifies that the
actual KdV owner constructs an ordinary complex vector and FFTW
constructs its default region as `(1,)`. Matched concrete KdV/FFT probes pass
86 assertions and the complete existing 61-assertion KdV suite in all four
processes. Six concrete signatures per process have no reports; all 16 complete
output/FFT bit patterns match within each runtime.
`candidate6-matched-compiler-review.md` records this classification and its
limits. The abstract report remains visible; no package ceiling was changed.

Ordinary CN median runtime ratios are 0.925–1.828 on minimum Julia and
1.157–1.851 on latest; equal-error ratios are 1.168 / 1.517. These are measured
correctness-repair costs, not speedup claims. The solved method and transactional
state staging perform work absent from the inaccurate parent approximation.

The independent audit identified and repaired statistics-array alias mutation,
BigInt subcycle promotion and zero-residual restarted warm-start behavior after
earlier suites were green. Every original failure and superseded candidate is
retained in evidence. The final public-boundary audit found no additional
confirmed defect in this bounded CN repair surface. Nonlinear coupling, full
Hall induction, GPU/distributed implicit solves and large-cluster scalability
remain separate capabilities in the master plan.


## Final composition

The repair is composed onto `083082d` without modifying any validated Julia
source. All 363 Julia files match the accepted candidate byte for byte;
`Project.toml` also matches its frozen hash. README and changelog describe the
new solve/failure controls and intentional CPU/Vector storage contract. The
master trackers retain the broader implicit-coupling and model gaps. Pushed
separately as `ce52836`; all 363 committed Julia files match the
validated composition and held hashes remain unchanged.
[Exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39999012813)
passes all six jobs. Complete logs verify 624,178 / 624,179 assertions on Julia
1.12.7 / 1.13.2. Full compiler trees (24 / 21 reports) match the separately
reviewed `ba3694c` trees byte for byte; the two commits differ only in Markdown.
Every job log contains the exact checkout SHA. Evidence and report classification
are recorded in the [audit tracker](DEEP_DEBUG_TRACKER.md).

A follow-up optimization should profile and cache the existing Fourier
operator's repeated plans/buffers behind the same operator and residual
contract. The measured allocation evidence above identifies that work; its
implementation and complete-call cost verification remain open.
