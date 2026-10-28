# Relativistic output publication and rollback (C1)

Status: isolated implementation on `fd0a1686deea4bbffa1a45fb9688a909aad6ec4a`;
full acceptance and root integration are pending. The parent and earlier V5
sources/evidence remain immutable. The original evidence is in
`HybridPlasmaPIC-audit-evidence/2026-10-27/relativistic-output-cone-c1`;
V2 preserves the test-only migration and V3 the strengthened diagnostic capture.
V3 stops after six passed commands and one failed fullstep gate. The separately
versioned oracle-only V4 continuation completes 20 numerical/resource commands
and carries those six receipts after exact source/harness/environment proof.
It then stops at eight baseline zero-report compiler expectations. Those failures
remain failed; full matched diagnostic collection and actual Metal acceptance
remain pending.
The current-main root integration evidence is in `takeover/c1-root-integration`.

## Supported behavior and evidence

The supported stored-velocity domain is unchanged: finite components, a finite
native `(vx/c)^2+(vy/c)^2+(vz/c)^2 < 1`, and an exact original-component
norm comparison on its accepted side within eight epsilons of one. Call this
D_narrow. Exactly interior vectors with rounded beta² at least one remain
unsupported. Shared input validation and scalar publication use one predicate.

VERIFIED on the actual parent: with c=qm=h=1, dt=2 and a particle initially at
rest, an axial electric field of 2^12 (Float32) or 2^27 (Float64) produces a
rounded component equal to c in scalar, cold, native and both actual batch
owners. The actual batch certificate accepts these moderate witnesses. Both
priming owners similarly publish -c. The exact physical result is strictly
inside the cone but beyond the midpoint between c and its predecessor. Extra
arithmetic precision alone cannot move its correctly rounded component inside.
The earlier design also bounds the large gamma error of inward projection;
C1 performs no projection. Receipts: `parent-bounded-11333.log` and TOML,
38 assertions/24 rows; exact representability evidence is retained in the
sibling `relativistic-output-closure-design` directory.

The uniform four-particle public parent fixture confirms partial mutation on
a rejected first step and an invalid successful established step. Subsequent
stepping, energy and persistence reject that invalid state. A separate actual
parent BigFloat/RoundUp witness at precision 11 uses c=dt=1, K=2^-44,
E=(0,1,0), B=0, and velocity components 1181/2048. Radiation returns components
591/1024: their exact squared norm is 1047843/1048576<1 but their rounded beta²
is one. The next public call rejects them. Receipt: `generic-parent-v1.log`.
The failed earlier axial hypothesis and its harness/log remain retained.

## Canonical change and compatibility

Boris checks native outputs, including its batch bypasses. A failed IEEE native
publication is re-evaluated by the existing guarded/cold owner before deciding.
Priming uses the same policy. If the resulting tuple is outside D_narrow, an
explicit supported-range `ArgumentError` replaces the inadmissible return.
No clipping, tolerance change, cone expansion or momentum state is introduced.
Raw numerical repair helpers retain their stage role; private raw batch/prime
helpers do not gain an atomic whole-array guarantee.

Generic-precision radiation retains its ambient native arithmetic but stages
three velocity copies and publishes them after every candidate passes. IEEE
radiation factor arithmetic is unchanged; its component-nonincreasing output
preserves its input domain. The existing K=0/dt=0 no-op order is retained.
Initialization checks every required relativistic species before resizing,
deposition and clock reset. The unused ion argument remains ignored when
mobile=false. Public zero-dt step semantics are unchanged.

All original kinematics inputs remain, including the 100,000 seeded ordinary
inputs per IEEE type. Independent stored-component D_narrow classification is
performed before interpreting each old output. Valid historical publication
keeps its original bits; formerly invalid publication is compared with the
independent 4096-bit stage reference or rejected if that rounded reference is
also unsupported. Existing normwise tolerances and end-to-end coverage remain.
This replaces the old blanket-return gate; that old gate is not relabelled as
passed. Earlier V5 priming conditioning compatibility and its threshold remain.

## Public transaction and persistence

Only public relativistic `step_empic!` wraps the unchanged advance body. For
standard writable CPU arrays and views it saves, once per call:

- Every required species x/v component, weight, id, tag, q and m.
- All 11 legacy or 24 general grid arrays, every electron and ion gather/work
  vector with its original length, and all four complex FFT buffers.
- Time, step and both prior timestep values.

An advance error restores those values in place, preserving array identities
and gather lengths. Local borrowed-midpoint cleanup remains in force. Preflight
failures occur before the cache is touched. The lazy undo cache itself is
scratch and may be created or resized on an expected output error. It is
excluded from the enumerated physical-state comparison and explicit persistence
schema. Checkpoint/archive data and staggered-phase continuation are unchanged;
a fresh destination starts with no undo cache and creates it on its next advance.
The concrete solver struct has a new internal field; raw Serialization of solver
objects across code versions is not a portable persistence interface.

Read-only particle metadata is copied to make the undo payload complete. This
does not establish forward correctness for arbitrary aliases. Custom write or
restore failures, resource exhaustion, side-effecting callbacks and concurrent
use of a solver/grid are outside the atomicity claim. Grid metadata, native FFT
plans and configuration are not mutated by the canonical advance and are not
copied. No transaction is imposed on nonrelativistic public steps.

## Resource contract

The pre-edit inventory and budget are frozen in `PRE_EDIT_BRIEF_AND_BUDGET.md`
in the evidence directory. Let b=sizeof(T), G=grid cells and n=the total required
particle count. Count actual old gather lengths when counts differ. The element
payload at matching steady counts is b*(19G+9n)+12n bytes for EMPIC1D and
b*(32G+(11+2D)n)+12n for general EMPIC, plus q/m scalars. At fixed warmed counts,
all cache arrays and retained size must remain stable. The declared retained
summarysize ceiling is logical payload plus 512*(A+1), where A is the number of
owned element arrays. This allowance is a measurement gate, not a runtime
constant or an RSS guarantee.

Each successful public relativistic step copies that enumerated logical payload
once; a failed advance adds one restore. There are no per-substep snapshots,
plan copies or hidden momentum. Generic RR's three copied velocity vectors are
an intentional separate per-call cost. Complete warmed allocations are measured
against the actual parent; logical copy counts do not measure bus traffic.
No timing acceptance is inferred under shared-machine contention.

## Verification tracker

| Gate | Status |
| --- | --- |
| C1.1 actual-parent numerical, native-batch/public and generic-RR reproduction | Verified on Julia 1.13.3 with compile=min; evidence above |
| C1.2 shared publication, preflight and typed whole-advance undo | Implemented in isolation; final source review pending |
| C1.3 borrowed arrays/midpoint, old gather lengths, ignored ions, changing dt/subcycles, fresh restart | Both final-source focused runs pass all 2,615 assertions; remaining integration gates are below |
| C1.4 retained ordinary/range fixtures independently classified; valid-parent bits | Both full seeded/companion runs pass 46,725 assertions with every original 100,000 input/type retained |
| C1.5 generic scale/vector/stored-precision/mode classification and atomic rejection | Both full focused runs include 432 independently classified rows (368 valid native-bit outputs, 64 input rejections) and distinct atomic output rejection |
| C1.6 exact old advance and frozen IEEE RR/DD/persistence arithmetic | Final isolated source passes 48 expression/byte checks on each Julia 1.12.7/1.13.3, including reconstruction of the entire original07792 IEEE prefix under the sealed12-symbol map; formatter syntax/idempotence passes20 checks |
| C1.7 warm retained/copy budget, ordinary allocations, concrete/full compiler, Aqua | Both normal-compiler resource campaigns pass: all 32 complete-step profiles per runtime have zero allocation increase; all 12 ordinary profiles retain three 0 B samples. Count-change retention and copy/retained ceilings pass. Failed-call heap/latency and near-cone costs need additive measurement. Full compiler/Aqua acceptance remains pending after the baseline zero-report failure |
| C1.8 min/latest affected regressions, task contexts, actual Metal, formatting | Both runtimes pass neighbors 75,380, affected IEEE RR 54,877, continued fullstep 38 stage/eight public, tasks 512, independent RR 9,728 and public RR 1,719 assertions. Actual Metal and remaining quality gates are pending. No relativistic resident Metal path exists in this scope |
| C1.9 root review, main push and exact-SHA GitHub CI | Pending; root owns integration |

Wrong fixes rejected by these gates include scalar-only checks that miss batch
bypasses; inward clipping that conceals representability; particle-only rollback
that leaves fields/FFT scratch/history changed; snapshots after resizing; and
generic per-particle writes that mutate an earlier particle before rejection.

C1 establishes supported local publication or explicit rejection. It does not
claim trajectory-wide accuracy or accurate recovery of a large physical gamma
from rounded velocity storage. Exact-cone widening requires auditing every gamma,
diagnostic, initialization and persistence owner; hidden-momentum storage requires
a separate state/schema and method design. Their [C2/C3 work orders](RELATIVISTIC_VALIDATION_PLAN.md#planned-domain-and-state-extensions) remain unimplemented.

## Current-main integration

- [x] Preserve the original V1 failed legacy run; independently classify every
  retained input before changing legacy publication expectations.
- [x] Freeze the V2 test-only migration and additive full-bit replay.
- [x] Verify the V3 complete diagnostic capture on both actual Julia runtimes;
  root checks all 2,023 source/harness/environment pins.
- [x] Compose C1 on `cd8b94d` and compare every one of the 465 Julia files
  against exact frozen C1/current-main source or their complete three-way merge.
  The subsequent RR docstring clarification changes no production statements;
  forwarding the docs-only ancestor to `e72db6c` preserves all 465 Julia files
  and all 3,017 initial/resolved pins.
- [x] Independently classify the failed fullstep fixture on both runtimes: two
  Float32 prime outputs are outside unchanged D_narrow; 2,072 supported parent
  publications retain all bits and ten raw stages retain physical accuracy.
- [x] Review all 20 completed numerical/resource continuation receipts and six
  explicitly carried passes, including every actual receipt/log hash.
- [ ] Complete the remaining four matched full compiler/Aqua collections and two
  actual Metal calls; review every diagnostic difference. The original V3/V4
  failed gates remain failed; pure Test context probes establish collection
  behavior only.
- [x] Review all six root setup graphs/providers and both 48-check source proofs;
  the explicit single-registry provenance migration changes no dependency source,
  version, edge or numerical test.
- [x] Pass both current-main focused calls: 50,147 assertions per runtime, plus
  minimum complete legacy replay: 46,725 assertions.
- [ ] Pass latest legacy replay and remaining root numerical, resource, neighbor,
  full compiler/Aqua and actual Metal gates.
- [ ] Finish final source/docs review, push separately and pass all six exact-SHA
  GitHub CI jobs.
