# Complete electromagnetic charge combination

Status: pushed as `4cd642c`; [GitHub acceptance run 38750764319](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38750764319)
is active. Native numerical, compiler, quality, concurrency, resource and solver
integration gates pass on Julia 1.12.7 and 1.13.1. Final public-step controls
also pass on both versions. GitHub completion remains required.

Candidate: `/tmp/hybrid-empic-charge-combination-v1`, based on `240f77c`.
Both EMPIC families share the complete charge owner. Combined current, hybrid
bulk-velocity ratios, pressure and MPI reductions remain separate boundaries.

## Verified baseline

The retained 22-assertion probe reproduces opposite overflowing species
densities/currents becoming NaN instead of exact zero, and background
subtraction becoming −Inf instead of finite `-floatmax(T)/2`. It runs against
the final MOMENT-001 source on Julia 1.12.7 and exits 0:
`/tmp/hybrid-moment-combination-baseline112.log`. These are defect witnesses.

Both native versions also pass 28 independent underflow/tie witnesses. With
smallest positive subnormal `u` and cell width 8, an ion of weight `3u` and an
electron of weight `-3u` should contribute total `3u/4`, rounding to `u`;
separate densities round to zero. Background `u` minus electron density from
weight `4u` should give `u/2`, rounding to even zero; separate density rounding
instead gives `u`. Float32/64 and both families are covered. Exact rational
identities establish the fractions independently of production rounding.
Evidence: `/tmp/hybrid-charge-combination-boundaries{112,113}.log`.

## Pre-edit brief

- Mode: repair complete charge reduction at the electromagnetic model boundary.
- Acceptance: when finite represented inputs lose range during deposition,
  normalization, species scaling or combination, round the complete cell charge
  expression once. Preserve ordinary order, visible nonfinite propagation,
  particle state and final electron-density scratch contents.
- Baseline: the overflow/cancellation and underflow/midpoint witnesses above run
  against the final isolated MOMENT-001 source, so correcting per-species
  normalization alone does not close this defect.
- Cause: **VERIFIED**, both `_charge_density!` implementations in
  `src/models/full_pic.jl` call `density!` before combining species or subtracting
  from `n0`. The exact contributions in the witnesses are lost at that boundary.
- Invariants: mobile charge is `q_i*n_i + q_e*n_e`; immobile charge is `n0-n_e`.
  Public EMPIC validators require electron/ion charges −1/+1. Retain ion-then-
  electron ordinary order and the existing immobile expression. Empty species,
  signed weights, deterministic shapes and writable particle views remain valid.
- Owners and surfaces: both model-family `_charge_density!` entry points and
  their species-add helpers; the shared private density/normalization owner
  carries range/eligibility status without changing `density!`'s public return.
  Initialization and each step already call these charge owners. Reuse the
  per-species exact moment accumulator and common final-rounding owner.
- Compatibility: no new public diagnostic, field, checkpoint schema or species
  model. Preserve ordinary arithmetic and nonfinite/unsupported-coefficient
  behavior. Finite exceptional combined outputs intentionally change.
- Strategy: retain normalization status through the private density path; one
  shared model charge owner combines those statuses with coefficient/addition
  checks. Eligible range loss evaluates all included species and the background
  in rational cell expressions before one conversion. Retain the normal scratch
  writes, and never infer eligibility from rounded output alone.
- Verification: exact rational rounding cells, overflow and subnormal/midpoint
  controls, nonfinite/unsupported coefficients, zero/unequal populations, all
  shapes and dimensions, ordinary bitwise/scratch comparisons, initialization
  and neighboring step tests, concrete/package JET, Aqua, allocations/timing,
  both native versions, formatter, exact diff and GitHub CI.

The cold path is bounded to one complete combination replay per call and reuses
the deposition arithmetic owner. Its work is proportional to particles/stencils
and grid cells, with operand-dependent BigInt storage. Measure its actual costs;
do not claim a fixed memory or latency bound.

Three incomplete repairs must be rejected: testing only Inf/NaN misses the
subnormal controls; adding rounded exact per-species values still loses the
midpoint case; making all calls exact changes ordinary rounding and cost.

## Verification

The final production source is unchanged across candidate test snapshots v1–v4.
The v3 manifest retains 178 source/test hashes; v4 changes only the existing
neutral-population test block to exercise public steps, variable timesteps and
one/two subcycles. Both final v4 runs pass all 1,720 assertions and exit 0:
`...permanent-v4-{112,113}.log`.

| Gate | Observed result on each native Julia version |
| --- | --- |
| Charge oracle and compatibility | 1,720 final assertions pass, including 256 public initialization/step/subcycle checks |
| Normalization/storage/threaded neighbors | 4,663 assertions pass |
| Concrete compiler analyses | All 64 calls have zero reports |
| Package JET / Aqua | Three JET assertions and nine Aqua checks pass; terminal diagnostics match baseline (22 / 21 reports) |
| Four-thread replay | 433 checks pass; outputs/scratch equal serial results and particle storage is unchanged |
| Ordinary actual-baseline comparison | 32 checks across 16 configurations pass; outputs/scratch match `240f77c` bitwise and both paths allocate 0 B |
| Cold resource controls | Six cancellation checks pass |
| Existing deposition and electromagnetic integration | 24,329 assertions pass on each version, both exit 0 |

Compiler/quality logs: `/tmp/hybrid-empic-charge-combination-v1-final-gates-v3-112.log`
and `...final-gates-v2-113.log`, both exit 0. Latest Julia separately passes the
final v3 distinct-ion-stencil fixture in `...permanent-v3-113.log`, exit 0.
The full neighboring solver driver includes `test_deposit.jl`, `test_empic.jl`
and `test_empic1d.jl`; logs are `...integration{112,113}.log`.

Diagnostic comparison retains every terminal report with multiplicity. Only
gensyms and exported `QRPivoted`, `Adjoint` and `Transpose` qualifications are
normalized; three identity assertions pass on each runtime. Raw qualification
mismatches remain in `...jet-comparison113.log` and `...jet-comparison-v2-112.log`.
Final comparisons `...jet-comparison-v3-{112,113}.log` both exit 0 with no
missing or added messages. Type-identity logs:
`/tmp/hybrid-empic-charge-jet-types{112,113}.log`.

Concurrency logs: `...concurrency{112,113}.log`, both exit 0. Worker calls use
native operands and do not change ambient BigFloat precision. Resource logs:
`...resources112.log` and `...resources-v2-113.log`, both exit 0. Shared-host
median candidate/baseline time ratios are 0.782–1.756 on 1.12.7 and 0.769–1.491
on 1.13.1; these measurements establish no speedup or portable ceiling. Cold
16/256-particle, 8/64/512-cell calls allocate 214,720–19,634,272 B per call;
these are allocation totals, not peak RSS.

## Retained harness failures

No production change or tolerance relaxation was made for these failures:

- The first expanded latest-version run passes 760 and fails eight background
  rounding assertions (`...permanent113.log`, exit 1). **VERIFIED:** three-
  particle 2D/3D TSC arithmetic stays in range; the expectation incorrectly
  requires once-rounded exact output. Independent `lost=false` and ordinary-
  rounding controls pass eight assertions (`...background-control113.log`).
  Those controls remain permanent. Exceptional fixtures use 16 particles.
- The first combined gate passes 1,953 normalization and 1,192 charge checks,
  then rejects an initialization fixture whose ion mass is 1 but constructor
  `mi` is 1836 (`...final-gates113.log`, exit 1). **VERIFIED:** existing
  `_require_empic_ions` enforces this agreement. The fixture sets `mi=ions.m`.
- The first latest resource launch cannot load SpectralOperators through the
  wrapper project (`...resources113.log`, exit 1); no measurements were taken.
  The final wrapper exposes the already-resolved production dependency UUIDs,
  retaining the native manifests and dependency versions.

The initial/middle test snapshots, final source/test manifests and driver/
environment hashes are retained alongside the candidate. The final v4 source
and environment hashes were rechecked after the completed gates.

- [x] Reproduce baseline and independent rounding boundaries.
- [x] Trace and repair canonical charge/status owners in isolation.
- [x] Pass native numerical, compiler, quality, concurrency and resource gates.
- [x] Complete final public-step controls and both native integrations.
- [x] Review final diff/format and promote the exact tested patch.
- [x] Commit and push to main (`4cd642c`) with its own GitHub acceptance run.
- [ ] Observe all GitHub jobs green for the pushed revision.

## Promotion evidence

Final patch SHA256: `27ba1f5c641f2f372e6da15a4d13afac01ec28dede02bb9dad782e3fee95d9ce`.
All 178 source/test hashes and 14 driver/environment hashes match their retained
manifests after the gates. Read-only JuliaFormatter 1.0.62 checks pass for all
five changed Julia files (`/tmp/hybrid-empic-charge-combination-v4-format.log`).
Both duplicate `_add_species_charge!` methods are removed; every remaining
`_finish_moment!` caller was checked against its private status return. Public
serial/threaded density calls still return their destination arrays. No new
catch, suppressed diagnostic, public field, API or checkpoint schema is added.
