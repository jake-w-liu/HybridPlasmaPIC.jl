# Distributed FFT precision selection

## Pre-edit brief (2026-10-26)

- Mode: implement the precision prerequisite for CAP-003 periodic operators.
- Requirement: `distributed_fft_plan(...; real_type=Float64)` accepts exactly
  Float32 and Float64, creates correctly typed plans/storage and agrees precision
  across every communicator rank before topology creation. Default behavior stays
  Float64. The existing constructor is the only API owner.
- VERIFIED baseline: `distributed-fft-precision/baseline-v1-*` records the default
  ComplexF64 input and MethodError for either explicit precision on both Julia
  1.12.7 and 1.13.2. This is a missing capability, not a violation of the existing
  documented default. PencilFFTs `plans.jl:296` and `global_params.jl:50` accept a
  real type; Transforms defines its supported FFTReal as Float32 or Float64.
- Canonical owner: HPIC Pencil extension validates configuration, serializes its
  value and constructs the dependency plan. Dependency owns precision-specific
  transform mathematics, allocation and transposes. No new dependency or copied
  FFT implementation is required.
- Invariants: valid integer shape, periodic-only geometry, same transforms,
  precision and communicator on all ranks, no empty partition, collective
  recoverable validation/planning failure, exact local origin exception.
- Compatibility: add one optional keyword with the existing Float64 default.
  Unsupported types/values produce ArgumentError collectively. Rank mismatch
  produces the existing disagreement error. Core missing-extension fallback
  remains the owner when dependencies are absent.
- Strategy: validate precision inside `_checked_fft_configuration`; include the
  selected type in the canonical serialized configuration; pass it to both
  GlobalFFTParams validation and PencilFFTPlan creation. Retain typed broadcasts,
  transport order and wrappers. No fallback, cast or alternate public method.
- Verification: both Julia runtimes and 1/2/4/8 real ranks; all eleven existing
  transform families, odd/even uneven shapes, serial FFTW reference, exact input
  and output precision, default/explicit Float64 agreement, local unsupported
  values from every origin, mismatched supported types and subsequent collective
  liveness. Preserve the old Float64 bounds and use the same precision-relative
  rule for Float32. Recheck hot-path outputs/allocations, complete compiler trees,
  extension loading, docs and actual CI step. Integrate the independently accepted
  diagnostic repair and rerun the combined extension before promotion.
- Wrong patches distinguished: accept the keyword but ignore it; validate type
  only locally; omit it from rank agreement; validate Float32 metadata but create
  a Float64 plan; narrow input arrays after a Float64 transform. Exact storage
  type, independent Float32 FFTW output and mismatch/liveness tests reject these.

## Progress

- [x] Inspect actual dependency constructors and retain baseline probes.
- [x] Freeze requirement, compatibility and numerical test contract.
- [x] Implement canonical precision propagation and permanent regressions.
- [x] Pass combined numerical/resource/compiler gates and update docs.
- [ ] Push main and verify all six exact-SHA jobs and complete compiler logs.

The frozen V2 candidate passes 17,177 combined assertions and 768 independent
boundary assertions per runtime across the complete 1/2/4/8-rank matrix. Its
complete compiler review and actual CI launcher pass. The original setup
measurement found 4,640 additional warmed bytes per Float64 construction.

V7 retains the V2 constructor body and adds one inlining annotation. All eight
original setup profiles remove 4,432 of those added bytes (95.5%), meeting the
fixed 75% target; 208 additional setup bytes remain. The serialized configuration
grows by two bytes. Both V7 compiler reviews and boundary suites pass. Earlier
failed optimization experiments remain archived. No latency or complete distributed
operator claim follows from these checks.

Landing verification (2026-10-26, applied to `563888f`): the CI entry
`test/test_distributed_fft.jl` passes under MPICH at 1, 2, 4 and 8 ranks on Julia
1.12.7 and 1.13.2. The complete `Pkg.test()` passes 798,477 assertions on 1.12.7
and 798,478 on 1.13.2, resolving the registered SpectralOperators 0.1.3. JET 0.12.3
reports 24 / 21 findings with normalized signatures identical to the parent.

## Integrated V2 contract

V2 starts at accepted diagnostic repair `563888f`. Its canonical constructor
remains byte-identical to the passing V1 preflight and its diagnostic helper
remains byte-identical to that accepted repair. The permanent diagnostic matrix
now also exercises Float32 and checks the returned error's precision, preserving
all existing Float64 checks. The actual CI entry already runs both matrices.
The original V1 source and baseline API-gap probes remain frozen. Final gates
run on this complete source before promotion.

## V3 setup-allocation pre-edit brief

- Requirement: preserve all V2 numerical/collective contracts while removing the
  measured avoidable plan-construction overhead. This is setup work; no transform
  latency improvement is claimed.
- VERIFIED baseline: `startup-v1-*` records exactly 4,640 additional warmed bytes
  per Float64 construction in all eight shape/family/runtime comparisons. The
  serialized configuration grows by two bytes. Full allocation profiles retain
  additional boxed plan/wrapper objects at the constructor boundary.
- Hypothesis to test: retrieving precision from a heterogeneous serialized-value
  tuple loses the singleton type known by the original keyword. Passing the same
  already validated/agreed keyword directly to PencilFFTPlan preserves it.
- Canonical change: retain precision in configuration validation, serialization
  and rank equality, but pass the original validated `real_type` to the existing
  dependency constructor. No new public method, type assertion, local-only
  validation, omitted agreement or copied FFT implementation.
- Acceptance: at least 75% reduction of V2's added setup bytes in every original
  Float64 profile; all Float32 behavior, per-rank invalid-type/mismatch liveness,
  complete compiler trees, full transform/diagnostic matrix and exact warmed
  operation resource/output comparisons must remain. Keep V2 frozen and retain
  its measurements; changing any gate requires a new explicit contract.
- Wrong patches rejected: forcing Float64; removing the precision agreement;
  keyword type restrictions that throw before collective error handling. The
  existing both-precision/invalid-origin/mismatch matrix distinguishes these.

### V3 result and V4 specialization experiment

V3 fails the fixed setup reduction gate in all eight Float64 profiles: the
allocations equal V2 exactly. The bypass hypothesis alone is rejected. Complete
outer-call inference still returns a concrete plan; profile traces show boxed
plan/result copies. V4 tests an explicit singleton-type specialization barrier
at the already validated plan-construction boundary. The private helper owns
pencil/plan construction and returns the canonical wrapper; no numerical or
collective path is duplicated. All configuration validation and agreement still
finish before that helper is entered. Unsupported precision values therefore
retain collective error handling. The same 75% reduction and full V2 gates apply.

### V5 inlining pre-edit experiment

V4's first runtime reduces setup allocation by 2,336 bytes but remains above the
unchanged 75% reduction target. V5 adds an inlining annotation only to the private
type-specialized construction boundary, testing whether the boxed return copy
can be removed without duplicating constructors or changing validation. Keep
V4 and all measurements; the same full numerical/collective and resource gates
apply. This is an allocation experiment, not a latency improvement claim.

### V6 typed configuration experiment

V4 and V5 both reduce setup overhead from 4,640 to 2,304 bytes in every original
Float64 profile on both runtimes; both fail the unchanged 75% reduction target.
Inlining the construction barrier adds no improvement. V6 carries the validated
precision as `Val(Float32)` or `Val(Float64)` inside the private configuration
instead of an erased DataType value, and specializes the existing construction
barrier on that tag. It removes the ineffective inlining annotation. Validation
still accepts exactly the same public values before constructing the tag, and
serialization/equality still agree the precision on every rank. The internal
payload is transient and not a checkpoint or public wire schema; all ranks run
the same package revision. Recheck its actual byte size and errors. No numerical
formula, unsupported-value guard or collective is removed; the existing 75%
setup reduction and complete correctness/resource/compiler gates stay fixed.

### V7 caller-specialization experiment

V6 removes only 16 additional setup bytes and increases the transient payload;
it also fails the unchanged target. V7 returns to the simplest frozen V2 body
and tests one inlining annotation on the public constructor. This tests whether
call-site precision specialization eliminates boxed intermediate plan returns,
without another helper or a changed configuration representation. V2–V6 remain
frozen. The same 75% setup target applies; no causal claim is made before this
controlled allocation comparison and compiler/profile review.
