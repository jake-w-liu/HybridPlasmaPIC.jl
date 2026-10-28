# Relativistic output publication and rollback (C1)

Status: implemented, locally accepted and pushed separately as `1bc618b`.
All six exact-SHA CI jobs pass. Root checks the complete minimum/latest compiler
trees against the accepted parent after verified unchanged-source line mapping;
all 24/21 reports and 243/181 frames match. The audit tracker retains the evidence.
This page owns the C1 contract and local gate status; the
[audit tracker](DEEP_DEBUG_TRACKER.md) owns main/CI acceptance.

The dated external evidence archive retains all earlier failures and frozen
plans. The accepted root campaign is `takeover/c1-root-integration`; additive
resource review is `takeover/c1-v7-numerical-root-review-v1.json`. A bounded
matched compiler comparison is accepted explicitly; the original zero-report
assertions remain failed. Earlier stopped runs are not relabelled as passes.

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

The additive near-cone corpus records extra candidate costs: scalar calls
allocate 16–6,384 bytes, forced-native controls 16–9,120 bytes, and prime/batch
calls 16–108,528 bytes; corresponding parent calls allocate 0–16 bytes. Both
native batch certificates are false throughout this corpus. The failure-call
fixtures allocate up to 1,017,624 warmed bytes, with retained undo storage
2,080–39,168 bytes. These actual Julia 1.12.7/1.13.3 observations are retained in
`relativistic-output-cone-c1-v7-capture/C1_E7_RESOURCE_HANDOFF_V1.md`; they do not
attribute allocation causes or establish universal heap/native-memory/latency
bounds. Near-cone cost reduction is a separate open optimization.

## Verification tracker

| Gate | Status |
| --- | --- |
| C1.1 actual-parent numerical, native-batch/public and generic-RR reproduction | Verified on Julia 1.13.3 with compile=min; evidence above |
| C1.2 shared publication, preflight and typed whole-advance undo | Current-main source reviewed; original advance bodies and the shared domain predicate are retained |
| C1.3 borrowed arrays/midpoint, old gather lengths, ignored ions, changing dt/subcycles, fresh restart | Both current-main focused runs pass 50,147 assertions; both neighbor runs pass 75,380 |
| C1.4 retained ordinary/range fixtures independently classified; valid-parent bits | Both current-main legacy replays pass 46,725 assertions; every original 100,000 input/type remains. Both-runtime 2,084-row classification retains 2,072 supported publications, ten raw physical stages and two explicit Float32 prime rejections |
| C1.5 generic scale/vector/stored-precision/mode classification and atomic rejection | Both focused runs include the original 432-row independent classification and atomic rejected publication |
| C1.6 old advance and frozen IEEE RR/DD/persistence arithmetic | Both source proofs pass; all 465 Julia files are composed against frozen/current-main owners. Whole tracked-file JuliaFormatter 1.0.62 passes with zero changed files |
| C1.7 warm retained/copy budget, ordinary allocations, complete compiler/Aqua | All eight root resource calls pass: 32 matched complete-step profiles add zero bytes and 12 ordinary profiles retain three zero-byte samples per runtime. Four complete quality collections retain eight inherited zero-report failures each; root reconstructs every frame/source/render difference and accepts the bounded comparison, including all 30 empty concrete captures |
| C1.8 min/latest regressions, task contexts and actual Metal | Actual Julia 1.12.7/1.13.3 numerical/neighbor/resource gates pass. Both actual functional Metal calls pass 3,543 assertions per runtime; this does not establish a resident relativistic Metal solver |
| C1.9 additive failure/near-cone observations | Independent replay checks 128 failure configurations, 664 stage rows, 52 complete profiles, 2,649 input pins and 196 terminal artifacts. Original physical/domain/bit/copy/retained gates pass; full raw observations retained |
| C1.10 main push and exact-SHA GitHub CI | `1bc618b` pushed; [all six jobs pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40694351554). Actual Julia 1.12.7/1.13.3 pass 1,242,177/1,242,178 assertions; both complete compiler trees match the accepted parent with checked source-line mapping |

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

- [x] Preserve every stopped source/oracle/quality/export attempt and raw failure.
- [x] Independently classify all retained inputs before changing publication
  expectations; preserve the physical references and normwise tolerances.
- [x] Compose all 465 Julia files on current main and retain all 3,017 original
  source/environment pins through the docs-only `e72db6c` / `430bf48` ancestry.
- [x] Review all six accepted setup graphs and both source proofs. The declared
  registry-provenance migration changes no dependency version, source or edge.
- [x] Complete both current-main focused, full legacy and neighbor campaigns,
  all eight resource calls and both actual functional Metal calls.
- [x] Review all four complete compiler/Aqua collections with every original
  diagnostic retained under the explicit bounded disposition.
- [x] Independently replay the additive failure/near-cone observations. Native
  binaries retain path, size, full streamed SHA and magic; text retains its full
  contents. Every exposed allocation/frame field remains present.
- [x] Pass whole tracked-file CI-compatible formatting with no source changes.
- [x] Finish source/docs review and push separately as `1bc618b`.
- [x] Pass all six exact-SHA CI jobs and review both complete compiler reports.

A separate inherited nonrelativistic Float16 arithmetic gap remains open:
four normal-ratio native/scalar/public witnesses on both runtimes exceed the
unchanged four-ulp bound against original-factor and explicit-qm exact references.
Root independently checks their rational 3×3 rotation oracle. C1 does not repair
that gap or the original charge/mass-factor range defect. See the
[factor and accuracy work order](PARTICLE_CHARGE_MASS_RANGE_PLAN.md).
