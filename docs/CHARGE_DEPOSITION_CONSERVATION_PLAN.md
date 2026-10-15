# Charge-deposition conservation

Status: VALID-001 passes local acceptance on Julia 1.12.7 and 1.13.2 and is
integrated. Unchanged wave thresholds, resource, concurrency, coupled evolution,
restart, compiler and quality gates pass. Exact-revision GitHub acceptance is
green at `e9c4fa7`: both Julia package jobs, formatting and 2-/4-/8-rank MPI pass
in [run 38767093527](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38767093527).

## Verified failure and cause

Validation case 20 evolves 10,240 CIC electrons on a periodic 32×4 mesh for
1,000 steps. Its final continuity residual is **9.889422779356765e-7**, above
the unchanged **1e-8** tolerance on both native versions. The frequency error
passes. The case is defined in
[`validation/20_empic_transverse_dispersion_2d/run.jl`](../validation/20_empic_transverse_dispersion_2d/run.jl).

**VERIFIED:** particle weights and charge are unchanged across the last step,
and their exact total charge plus background is zero at both endpoints under
the implementation's rounded cell-volume convention. Replaying
`_charge_density!` from captured particle positions exactly reproduces both
stored density fields. Their total charge nevertheless changes by **115/2^53**.
The canonical trace is `deposit_scalar!` in
[`src/coupling/deposit.jl`](../src/coupling/deposit.jl), `density!` in
[`src/coupling/moments.jl`](../src/coupling/moments.jl), then `_charge_density!`
in [`src/models/full_pic.jl`](../src/models/full_pic.jl).

The exact difference in the charge-error budget is:

| Stage | Contribution to the change in summed charge density |
| --- | ---: |
| Rounded shape partition | 2.959871930885427e-18 |
| Rounded stencil products | 7.310759377774246e-17 |
| Cell accumulation | 1.2522712962035057e-14 |
| Cell-volume division | 1.6878435544561575e-16 |
| Background subtraction | 0 |

The script retains exact rational contributions and verifies their sum, rather
than relying on this rounded display. Cell accumulation accounts for about
98.08% of this measured drift; improving only that stage has not been shown
to satisfy the complete acceptance criteria.

**VERIFIED:** periodic spectral divergence has zero DC. After removing only
the existing unrepresentable non-DC null modes, the exact density-change mean
is `115/2^60` and its maximum magnitude is `116286195/2^60`. Thus the residual
ratio for these fixed density snapshots is at least
`23/23257239 ≈ 9.889394007603395e-7`, already above the tolerance. Changing
current alone cannot remove this DC term. This statement is confined to these
snapshots and the existing periodic divergence operator.

An independent 512/1,024-bit direct DFT agrees: normalized residual
9.889406491600726e-7, with a non-DC residual magnitude of about 1.259e-22.
The exact DC bound does not depend on that floating-point reference.

## Evidence and reproduction

All listed completed scripts exit successfully on both native versions against
the immutable diagnostic snapshot `/tmp/pic005-diagnostic-prototype-rpzn57rx`
(full-PIC source SHA256
`ba62bd1a5c060b847205a5329e8563089db5c2ea275ae82776bb5b1a55bc6a50`).

- `/tmp/hybrid-validation20-recon.jl`: unchanged evolution with final-state
  capture, four grouped checks/version; logs
  `/tmp/hybrid-validation20-recon{112,113}.log`.
- `/tmp/hybrid-validation20-dc-probe.jl`: nine independent-DFT checks/version;
  logs `/tmp/hybrid-validation20-dc{112,113}.log`.
- `/tmp/hybrid-validation20-exact-dc-bound.jl`: six exact rational/null-mode
  checks/version; logs `/tmp/hybrid-validation20-exact-dc{112,113}.log`.
- `/tmp/hybrid-validation20-deposition-recon-v2.jl`: 13 particle/deposition
  checks/version; logs `/tmp/hybrid-validation20-deposition-v2-{112,113}.log`.
  Captured particles, fields and exact stage errors are retained under the
  matching `/tmp/hybrid-validation20-deposition-v2-{112,113}/deposition.bin`.
- `/tmp/hybrid-validation20-deposition-budget.jl`: nine exact accounting
  checks/version; logs `/tmp/hybrid-validation20-deposition-budget{112,113}.log`.

Run each driver with its native Julia environment, snapshot path and capture
directory/file as its two positional arguments. The first deposition driver
failed because its harness used unqualified `width`; those exit-1 logs remain
in `/tmp/hybrid-validation20-deposition{112,113}.log`. The corrected driver uses
the package's actual `width` function; no solver code was changed.

## Repair requirements and next work

1. Establish the canonical conservation boundary. Compare compensated cell
   deposition with reconciliation against the actual particle charge and
   background. The fixed-state experiment below checks both hypotheses; neither
   is a verified complete repair. The full-wave results below narrow the choice.
2. Preserve genuine net charge, changes in particle weights/populations, mobile
   species and the diagnostic's DC response. A uniform zero-charge assumption,
   DC masking or relaxed threshold is not an acceptable repair.
3. Define rounding guarantees for shape partition, signed weights, density
   division and charge combination. Cover built-in NGP/CIC/TSC and inspect the
   custom-shape contract. Keep nonfinite inputs observable. Do not silently
   replace the existing cell-volume convention while fixing accumulation.
4. Audit the shared deposit/gather adjoint contract and other density consumers
   before changing their results. State compatibility explicitly, including
   ordinary trajectories, initialization, checkpoints and both EM interfaces.
5. Complete an evidence-supported nine-field implementation brief before
   production edits. Require independent total-charge and spatial-mode oracles,
   both Julia versions, unchanged validation 20, Float32/64 and 1D–3D, mobile
   ions, changing timesteps/populations, cancellation and resource checks.
   Measure allocation and runtime cost; do not infer them from operation count.

This finding is distinct from the remaining current-materialization range and
severe FFT-cancellation accuracy work in the
[continuity plan](CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md).

## Fixed-state design experiment

The first bounded experiment passes 19 checks on each native version, with
identical measured values. It replays the captured before/after particles and
reapplies the existing current correction to the original raw current. No
production code, diagnostic or tolerance is changed. Driver and analysis brief:
`/tmp/hybrid-validation20-deposition-methods.jl` and
`/tmp/hybrid-validation20-deposition-analysis-brief.md`; logs
`/tmp/hybrid-validation20-deposition-methods{112,113}.log`.

| Accumulation / charge treatment | Change in summed charge density | Continuity residual |
| --- | ---: | ---: |
| Existing serial cell sum | 115/2^53 | 9.889422779356765e-7 |
| Kahan cell sum | 1/2^53 | 8.601256684984105e-9 |
| Neumaier cell sum | 1/2^53 | 8.601256684984105e-9 |
| Exact sum of the same rounded stencil products, rounded once per cell | 1/2^53 | 8.601256684984105e-9 |
| Existing deposition with uniform reconciliation to actual particle charge | 0 | 2.877175015333387e-12 |

The last prototype computes its target from the actual weights, species charge
and background under the existing cell-volume convention. A changed-weight
control has target -0.0125 and represented sum -0.012499999999999999, so it
preserves nonzero charge instead of clamping it away. Original captures remain
unchanged. Exact output-rounding error accounting is retained in the driver.

These fixed-state results do not establish an implemented repair. Full-wave
results follow; other dimensions/precisions/shapes, mobile and signed populations,
nonfinite behavior, spatial/adjoint guarantees and measured cost remain open.
Compensated accumulation passes this one threshold narrowly; no general
sufficiency claim follows. Choose the production owner and rounding contract
only after those comparisons.

## Full-wave experiment brief

- **Mode:** isolated numerical design experiment, with no production edits.
- **Acceptance:** run the unchanged 1,000-step validation case on both native
  Julia versions; retain frequency and continuity thresholds and failing runs.
  Compare final fields and actual total charge, and record execution costs.
- **Baseline:** the public case fails on both versions; fixed-state comparison
  gives the measured results above. The next copies start from pushed `85ddb1b`.
- **Verified gap:** the fixed-state test reprojects the original raw current.
  Changing generic deposition also changes deposited current during evolution,
  so that result alone does not establish coupled behavior.
- **Contract:** one prototype adds Kahan compensation to canonical scalar
  deposition; another reconciles total EM charge density to actual particle
  weights/charges and background after deposition. Both EM interfaces use the
  latter helper. The prototype explicitly limits reconciliation to finite raw
  densities/weights and a finite positive represented cell volume; it is not a
  proposed final nonfinite or extreme-range contract.
- **Owners:** temporary copies of `deposit_scalar!` or `_charge_density!` only.
  Preserve the actual current projector, diagnostic, field solver and validation.
- **Compatibility:** no live source/API/schema/dependency change. Generic
  compensation's allocation and adjoint implications remain open; exact-rational
  target summation is an analysis backend whose production cost is unverified.
- **Strategy:** baseline, compensated-deposition and reconciled-charge copies;
  capture the same final physical state and metrics from each complete run.
  Reconciliation uses the existing rounded cell-volume convention and actual
  net charge, without DC masking or a neutrality assumption.
- **Verification:** native 1.12.7/1.13.1, exact source hashes, untouched validation
  thresholds, actual final charge totals and field differences. Treat first-run
  timings as compilation-inclusive observations under shared load. Follow with
  broader precision/shape/species/range and warmed-cost work before selecting a
  production implementation.

## Full-wave results

The three isolated copies of pushed `85ddb1b` completed the same 1,000-step
wave on both native versions. The frequency error remains
0.013636642520216225, below 0.03, in all six runs. The continuity threshold is
unchanged at 1e-8.

| Prototype | Final continuity residual, both versions | Final-step summed-density drift |
| --- | ---: | ---: |
| Baseline | 9.889422779356765e-7, fails | 115/2^53 |
| Kahan scalar deposition | 6.87978521988381e-8, fails | -1/2^50 |
| EM charge reconciliation | 2.877175015333387e-12, passes | 0 |

**VERIFIED:** compensated scalar deposition is insufficient for this complete
case, despite passing the fixed-state comparison. Both candidate failures are
retained as exit-1 runs; the experiment orchestrator also exits 1. Reconciliation
passes five grouped checks/version. Its final E, B, corrected/raw current,
positions and velocities are numerically identical to the baseline in the
32-assertion final-state comparison. This result applies to this wave only.

The exact-rational prototype has substantial cost. Compilation-inclusive
allocated bytes for the complete run are 13,172,242,280 on 1.12 and
12,636,754,200 on 1.13, versus 542,162,464 and 12,239,056 for their baselines.
Observed wall times are 42.08/35.52 s versus 13.54/5.77 s. These single runs
include compilation, capture and output under shared machine load; they are
not warmed kernel benchmarks or portable performance bounds. A cheaper charge
summation/reconciliation algorithm and the broader acceptance matrix remain
required before production adoption.

Evidence: `/tmp/hybrid-charge-wave-driver.jl`, isolated copies and SHA256
manifests under `/tmp/hybrid-charge-wave-v1`, per-mode/runtime logs
`/tmp/hybrid-charge-wave-v1-{baseline,kahan,reconcile}{112,113}.log`, and
`/tmp/hybrid-charge-wave-v1-gates.log`. Final captured states are in the matching
`artifacts-*` directories. Comparison driver/log:
`/tmp/hybrid-charge-wave-compare.jl`, `/tmp/hybrid-charge-wave-v1-compare.log`.
No live deposition or EM evolution code was changed by this experiment.

## Exact-sum resource experiment

An isolated second prototype groups adjacent equal finite values before exact
summation. Each run contributes its actual value times its integer count;
varied and signed weights keep their actual charge. Skipping multiplication
for singleton runs removes the extra allocation observed in the first version.
This is an exact arithmetic identity, not a neutrality or equal-weight assumption.
The nine-field experiment brief is `/tmp/hybrid-charge-repeated-sum-brief.md`.

Both versions pass **665 oracle/state checks** for each of the two prototypes:
empty/constant/alternating/random Float32/64 sequences, signed cancellation,
subnormals, finite extremes, reversed order, unchanged inputs and all 17 captured
state fields. The second full-wave prototype also passes its five checks/version,
with the same 2.877175015333387e-12 continuity residual and numerically identical
captured state to the original exact-rational reconciliation run.

Three warmed allocation samples per fixture agree on both Julia versions:

| 10,000-value sum | Direct rational sum | Second run-grouping prototype |
| --- | ---: | ---: |
| Constant Float32 | 6,080,104 B | 960 B |
| Constant Float64 | 6,119,152 B | 960 B |
| Alternating Float64 | 6,125,976 B | 6,125,976 B |
| Seeded random Float64 | 6,091,880 B | 6,091,880 B |

Alternating/random Float32 samples also retain the direct sum's 6,080,104 B.
The complete second wave run allocates 714,790,128 B on 1.12 and 181,964,928 B
on 1.13, including compilation/output. Observed times are 16.01/7.99 s under
shared load. These are separate measurements from the warmed sum microchecks;
there is no portable timing or peak-memory claim. Varied-weight cost remains
large, and the latest-version wave still allocates substantially more than the
12,239,056 B baseline.

Evidence: `/tmp/hybrid-charge-wave-runs-v2-{112,113}.log` and
`/tmp/hybrid-charge-repeated-sum-oracle-v2-{112,113}.log`, driven by
`/tmp/hybrid-charge-wave-runs-v2-gates.py` and
`/tmp/hybrid-charge-repeated-sum-oracle-v2.jl`. All exit 0. First-prototype logs
remain at `/tmp/hybrid-charge-wave-runs{112,113}.log` and
`/tmp/hybrid-charge-repeated-sum-oracle{112,113}.log`.
Second-prototype full-PIC SHA256: `29424eb42c712ef8da3edde2e97bcf138472ef4204f401dce70a2dbab029d606`.

### Bounded dyadic accumulator and finite spatial matrix

The next isolated prototype uses an Int128 sum only when the significand width,
observed exponent span and particle count prove that every partial sum fits.
Larger spans use the existing exact-rational path. Constant inputs retain an
exact value-times-count shortcut. No production deposition code has changed.

On native Julia 1.12.7 and 1.13.1, the final v4 prototype passes **70,100 exact
oracle checks per version**, including every finite Float16 value, random bit
patterns, signed cancellation, subnormals, extreme exponents and accumulator
boundaries. In the measured 10,000-value fixtures, varied Float32 sums allocate
312 B, alternating Float64 328 B and random Float64 352 B. These are warmed
allocation measurements, not peak-memory or portable timing claims.

An ambient-precision replay exposed a prototype defect: converting the rational
correction directly to a hardware float used ambient BigFloat precision. The v3
replay fails two of 30 checks per version. Reusing the existing `_from_exact`
conversion fixes all **30 checks per version** at precision 2, 4, 8, 32 and 256;
the caller's precision is preserved. The failing experiment remains recorded.

The v4 full wave passes its five checks per version at the unchanged 1e-8 limit:
residual **2.877175015333387e-12**. Complete-run allocations, including compilation
and output, are 557,835,216 B / 19,986,176 B on 1.12 / 1.13. Timings were taken
under shared load and do not establish a speed ratio. All 17 captured state
fields match the original exact-reconciliation wave on both versions (34 checks,
`/tmp/hybrid-charge-dyadic-state-v4.log`, exit 0).

The finite-geometry matrix passes **58,817 assertions per version** across 512
configurations: Float32/64; general EM dimensions 1–3 and legacy 1D; unit and
fractional volumes; NGP/CIC/TSC/custom CIC shapes; mobile/immobile populations;
empty, neutral, varied and signed weights; and a changed particle weight. Exact
physical charge, per-cell nearest rounding, total rounding budgets, non-DC
Fourier changes and unchanged particle/scratch storage are checked independently.
This does not yet cover extreme cell volumes, public coupled evolution, restart
or nonfinite-input semantics. Those remain production-adoption gates.

Evidence: `/tmp/hybrid-charge-v4-oracle{112,113}.log`,
`/tmp/hybrid-charge-v3-ambient{112,113}.log` (expected failures),
`/tmp/hybrid-charge-v4-ambient{112,113}.log`,
`/tmp/hybrid-charge-wave-dyadic-v4-{112,113}.log`, and
`/tmp/hybrid-charge-v4-matrix{112,113}.log`. Final passing processes exited 0.
The isolated full-PIC SHA256 is
`85ced08fc67446c7ff61201d22becfca22f34bda46c7849361e1787e06e3bfa6`.

### Confirmed prerequisite: cell-volume range

**VERIFIED:** `density!` materializes `prod(g.dx)` before dividing the deposited
weights (`src/coupling/moments.jl`). On a 2×2 NGP grid with one particle at the
origin, finite individual cell widths and weights produce these results on both
native Julia versions:

| Precision | Domain length in both axes | Weight | Observed occupied-cell density | Exact-expression rounded density |
| --- | ---: | ---: | ---: | ---: |
| Float32 | 1e25 | 1e30 | 0 | about 4e-20 |
| Float32 | 1e-25 | 1e-30 | Inf | about 4e20 |
| Float64 | 1e200 | 1e300 | 0 | about 4e-100 |
| Float64 | 1e-200 | 1e-300 | Inf | about 4e100 |

All four cells become nonfinite in the small-volume cases. The exact rational
oracle uses the stored widths and weight; expected occupied-cell results are
finite and nonzero. Each native version passes 12 witness checks, which confirm
the defect rather than certify the baseline density implementation. Evidence:
`/tmp/hybrid-density-volume-repro.jl`, `/tmp/hybrid-density-volume112.log`, and
`/tmp/hybrid-density-volume113.log` (exit 0).

Repair must address canonical moment normalization before the reconciliation
prototype is adopted. Inspect numerator deposition, stencil products, signed
accumulation and the momentum/charge callers as well as volume formation; fixing
only the final division cannot recover a numerator already lost to range.
Specify the finite-volume rounding contract, preserve ordinary-path behavior
where required, then test exact occupied/empty-cell results, variable weights,
shape and dimension variants, charge coupling, compiler reports and resources.
The broader numerator cases were subsequently measured by the probe below.
The inspected callers include `momentum!` and `current!` in
`src/coupling/moments.jl`, and both general/legacy `_charge_density!` paths in
`src/models/full_pic.jl`. EMPIC validates electron/ion charges as −1/+1, so its additional range checks
must cover cancellation between species or against the background. The public
`current!` also multiplies by the species charge after materializing momentum.
The public-current witness below confirms charge-scaling failure; EMPIC species/background cancellation remains unmeasured here.

### Additional normalization witnesses

A second probe confirms that a cell-volume-only fix is insufficient. Both native
versions pass **24 witness assertions** for six baseline defects:

- With 16 colocated NGP particles of weight `floatmax(T)/4` and cell width 1e8,
  the raw sum overflows before division. The complete Float32/64 densities are
  about 1.36e31 / 7.19e300; the implementation returns Inf.
- A single CIC particle in 2D has a cross-node contribution that underflows
  before normalization. The complete densities are about 1e-20 / 1e-250 for
  Float32/64; the implementation returns zero. Both cell volumes are finite
  and positive. The exact oracle uses the actual stored stencil coefficients.
- `current!` forms overflowing momentum before applying a small finite species
  charge. The complete Float32/64 currents are about 1.36e19 / 7.19e108;
  the implementation returns Inf in all three components. These are public
  moment calls, separate from EMPIC's unit-charge restriction.

All positions, particle weights, velocities and charges used here are finite.
Evidence: `/tmp/hybrid-moment-range-expanded.jl` and
`/tmp/hybrid-moment-range-expanded{112,113}.log` (both exit 0). These passing
witness checks demonstrate baseline defects; they do not validate a repair.
The required implementation must preserve the complete weighted stencil sum,
cell normalization and applicable species coefficient through range-sensitive
intermediates. Keep generic deposit/gather adjoint semantics and validate both
scratch/no-scratch moment calls, signed/zero weights and charges, ordinary
results/resources, and coupled density/current consumers before adoption.

### Control against premature weight normalization

**VERIFIED:** predividing particle weights by the cell volume would regress a
currently correct case. A 2×2 unit-domain CIC grid with one particle at (0.25,
0.25) and weight `floatmax(T)/2` deposits one quarter of that weight into each
cell. The cell volume is one quarter, so each final density equals the original
finite weight. The current implementation matches the exact rational result
in all cells for Float32 and Float64. Predividing that weight by the volume
instead yields Inf before applying the stencil.

The independent control passes **eight assertions on each native version**.
Evidence: `/tmp/hybrid-moment-prescale-control.jl` and
`/tmp/hybrid-moment-prescale-control{112,113}.log`, both exit 0. Preserve this
case alongside the overflow/underflow witnesses: changing multiplication order
alone cannot serve as the complete normalization repair. The subsequent
[per-species moment repair](MOMENT_NORMALIZATION_PLAN.md) records implementation
and acceptance; [complete EMPIC charge combination](EMPIC_CHARGE_COMBINATION_PLAN.md)
remains a separate follow-up.

## Post-combination baseline and complete correction

The frozen source matching pushed `4cd642c` reproduces validation 20's unchanged
failure on both native versions: residual **9.889422779356765e-7** against
**1e-8**, with frequency error **0.013636642520216225** passing **0.03**. The
exact summed-density drift remains **115/2^53**. Each diagnostic driver passes
four defect-witness assertions and exits 0; the validation metric itself fails.
Logs: `/tmp/hybrid-charge-wave-post-combination{112,113}.log`. Source SHA256:
`6aeb60babf6ce7e6f33ddbccbddd31c595ce9975d699669e4783bffb6808ce8e`.

A nine-check probe per version compares represented `prod(dx)` with the exact
product of stored widths. For the retained wave, the former gives target zero;
the latter gives a small nonzero target. Both fixed-state replays have residual
2.877175015333387e-12. Exact target, output-rounding and drift identities are in
`/tmp/hybrid-charge-volume-convention{112,113}.log`. Passing that one residual
cannot select the volume contract.

An isolated complete-subtraction kernel brackets the exact rational uniform
shift by adjacent hardware values. When the two endpoint subtractions agree,
monotonic nearest rounding certifies the complete result; otherwise that cell
uses the existing exact conversion owner. It handles a shift that overflows
while its final subtraction is finite. **24,458 independent rounding-cell,
boundary, random, view, invalid-input and ambient-precision checks pass per
version**, exit 0. Brief, source and logs:
`/tmp/hybrid-charge-uniform-shift-{brief.md,prototype.jl}` and
`/tmp/hybrid-charge-uniform-shift{112,113}.log`.

Resource probes pass 24 checks per version. For the captured 128-cell wave,
represented-volume shifts require no exact cell conversion and allocate
616–664 B; exact-width shifts require three exact cell conversions and allocate
3,936–4,000 B. Evaluating every cell exactly allocates 119,568–139,360 B.
Two deliberately uncertain 256-cell controls require all 256 exact conversions
and allocate 211,632–454,856 B, so this is not a universal fast path. These are
complete kernel allocation totals, not peak RSS. Timing samples remain specific
to the shared host. Logs: `/tmp/hybrid-charge-uniform-shift-resources{112,113}.log`.

### Coupled v5 experiment — pre-edit brief

- Mode: isolated conservation design experiment, with no live production edits.
- Acceptance: both volume conventions must run the unchanged full validation
  and retain its frequency/continuity thresholds. Record actual physical state,
  total-charge rounding and resources before selecting a production contract.
- Baseline: post-MOMENT-002 witnesses and complete-subtraction controls above.
  Start both copies from pushed `4cd642c`, preserving its tested charge owner.
- Design gap: the earlier prototype assumes a finite represented volume and
  rounds the correction before subtraction. Range repair now evaluates complete
  width products; the ordinary density path still uses the represented product.
  Selecting one global target is an explicit compatibility decision.
- Invariants: use actual signed species weights and the immobile background;
  never assume neutrality or suppress the diagnostic DC response. Compare one
  fixed geometry convention per run, independent of particle position/range
  flags. Preserve visible nonfinite results and all particle/scratch storage.
- Owners: the isolated canonical EMPIC charge owner after complete combination;
  reuse the previously verified bounded exact-sum prototype and complete-shift
  kernel. Generic deposition/gather and solver/current/diagnostic code stay fixed.
- Compatibility: only built-in NGP/CIC/TSC shapes enter this experiment's global
  normalization. Existing custom-shape fixtures can deliberately have nonunit
  stencil sums; their accepted coefficient semantics must remain unchanged.
  No public shape trait/API is introduced. This boundary needs explicit tests
  before production adoption. Nonfinite density/operand cases retain their
  existing arithmetic; finite complete normalization is evaluated separately.
- Strategy: two immutable copies. One targets the represented volume when it is
  finite and positive, otherwise the exact width product; the other always
  targets the exact width product. Eligible finite built-in results receive one
  complete uniform correction. Keep the two policies and source hashes distinct.
- Verification: both native full waves and independent state/target accounting;
  then exact spatial/rounding tests, extreme geometry, custom/nonfinite controls,
  public coupled evolution/restart, compiler/quality and cost gates before any
  production adoption. This experiment alone does not close VALID-001.

### V5 coupled results and rejected volume rule

Both policies complete the unchanged 1,000-step wave on native Julia 1.12.7 and
1.13.1: each run passes five assertions, with frequency error
0.013636642520216225 and continuity residual 2.877175015333387e-12 against the
original 0.03/1e-8 thresholds. Represented-volume totals and drift are exactly
zero. Complete-width totals differ by 1/5070602400912917605986812821504 between
the final two snapshots. Logs:
`/tmp/hybrid-charge-reconcile-v5-{represented,complete}{112,113}.log`.
These measurements include compilation/capture; they are not warmed benchmarks.

A strict bitwise physical-state hypothesis passes for the represented policy
but fails for the complete policy on both versions: 27 passes/five failures per
comparison driver. All positions, weights and geometry remain identical. The
complete policy's maximum final differences are 1.893e-23 in E, 2.426e-24 in B,
1.935e-24 in J, 1.117e-23 in Jraw, and 3.389e-21 in velocity. The logs retain the
full component measurements and failed assertions; no equality requirement was
silently relaxed. This disproves bitwise compatibility for that policy; its
physical acceptance still requires convergence and rounding-error evidence.
Sources: `/tmp/hybrid-charge-reconcile-v5-state.jl` and `...state{112,113}.log`.

**VERIFIED: finite positive `prod(dx)` is insufficient for choosing the
represented-volume target.** On stored Float64 widths
(1e-160, 1e-160, 1e200), the native product is about 0.999988867182683 times the
exact product. Per-species range recovery already accounts for the lost
intermediate. The represented reconciliation then injects about 1.392e114 into
each empty cell, equivalent to 9.708e9 ULPs of the occupied cell. The Float32
counterpart injects about 29,002 occupied-cell ULPs; a positive subnormal volume
also differs by a factor 4/3. This rejects v5's finite-positive decision rule.
The exact-width policy and any alternative geometry-only rule must pass these
controls before adoption. Probe and first completed evidence:
`/tmp/hybrid-valid001-contract-audit.jl` and
`/tmp/hybrid-valid001-contract-audit-represented113.log` (144 witness assertions,
exit 0; this is evidence of the defective policy, not its acceptance).

The independent expanded matrix compares against the actual archived charge
implementation, loaded in a separate module. It checks complete shifted-cell
rounding, exact global error accounting, non-DC Fourier rounding bounds, signed
and changed weights, empty/unequal populations, custom/nonfinite compatibility,
unrepresentable native volumes, and finite cells with totals beyond hardware
range. Driver: `/tmp/hybrid-charge-reconcile-v5-matrix.jl`. Initial execution
could not load transitive SpectralOperators from the test environment; dedicated
matrix environments expose the existing SpectralOperators/MPI dependencies.
The original failed log remains retained; production source is unchanged.
Matrix completion, contract choice and broader public/restart/compiler/resource
gates remain pending. VALID-001 is still open.

The complete-width matrix now passes **60,874 assertions on each native
version**, both exit 0: 44,125 ordinary spatial/charge, 15,177 extreme-volume,
1,040 custom/nonfinite, and 532 finite-cell/unrepresentable-total checks.
Logs: `/tmp/hybrid-charge-reconcile-v5-matrix-complete112.log` and
`...matrix-complete113-v2.log`. Only test-environment dependency exposure changed
after the initial failed launch. Both source hashes still match the recorded v5
complete candidate and archived MOMENT-002 baseline.

### Production candidate — pre-edit brief

- Mode: isolated repair candidate for periodic EMPIC charge conservation,
  retaining the raw charge arithmetic as an independently testable stage.
- Acceptance: both EMPIC families capture charge whose exact represented sum
  differs from actual particle/background charge only by the final cell-rounding
  errors. Pass the unchanged validation-20 tolerances, independent spatial/charge
  bounds, public variable-step/subcycled/restart tests, existing integrations,
  compiler/quality, concurrency and measured resource gates on both Julia versions.
- Baseline: pushed `6ad7276` includes the verified current/charge range repairs.
  The wave failure at `4cd642c`, rejected represented-volume witnesses and the
  completed v5 exact-width gates above establish the required correction. Current
  changes do not touch the charge owner; final integration must still be tested.
- Cause: VERIFIED in the exact wave accounting above: rounding during stencil
  evaluation, accumulation and charge subtraction changes the represented DC
  total between snapshots. The periodic spectral divergence cannot correct DC;
  the unchanged diagnostic observes the resulting residual. A final-cell finite
  check also cannot certify the accuracy of an earlier volume product.
- Contract: first compute the existing complete raw charge from represented
  coefficients and inputs. For eligible finite results from built-in NGP/CIC/TSC,
  set V to the exact product of stored widths and Qgrid to the actual signed
  species weight/charge total divided by V, plus Ncell*n0 for an immobile
  background. Set delta=(exact sum(raw)-Qgrid)/Ncell, then round each complete
  raw[I]-delta once. The exact final sum error is the sum of those rounding errors;
  do not promise exact equality to a non-dyadic target. Custom shapes retain raw
  coefficient semantics; nonfinite/ineligible data retain visible propagation.
- Canonical owner: split the existing common charge implementation into one raw
  deposition stage and its existing `_empic_charge_density!` capture owner.
  Keep both `_charge_density!` delegates, public initialization/step callers,
  moment/density APIs, particle/scratch storage, diagnostics and field solvers.
  No model field, archive schema, public option, or dependency is added.
- Compatibility: finite built-in charge snapshots intentionally receive the
  uniform correction, so their ordinary bitwise and per-cell nearest-stencil
  contract changes. Preserve and test the previous raw stage's exact range and
  ordinary-order contracts separately; retain the original assertions and add
  independent final-capture assertions. Custom coefficients, raw density/gather
  adjoint semantics, scratch contents and error propagation remain unchanged.
  Existing archives retain saved snapshots; exact continuation between runs of
  the new implementation is required, while cross-version bitwise trajectories
  are not promised. Inspect restore behavior before accepting this boundary.
- Strategy: port the verified bounded dyadic sum and monotonic complete-shift
  certificate into the isolated common owner. Reuse canonical exact conversion;
  no ambient BigFloat mutation, tolerance-dependent policy, special wave branch,
  suppressed diagnostic, or catch. Keep complete grid/cell arithmetic independent
  of particle positions and range flags. Avoid new cached model state.
- Verification: preserve both baseline and rejected prototype evidence; reject
  neutralization by checking nonzero/changed/signed/empty/unequal totals; reject
  rounded-volume normalization using recovered/subnormal controls; reject
  two-round correction using ties and overflowing shifts; reject custom-shape
  normalization with nonunit coefficients. Add coupled convergence, checkpoint
  continuation and failure/scratch controls. Compare full terminal JET diagnostics
  and measure complete ordinary/cold calls against actual `6ad7276` source before
  promotion. The matrix and wave results alone do not close these remaining gates.

### V5 complete-call resource characterization

The actual `6ad7276` charge owner and unchanged complete-width v5 prototype
pass **504 resource/oracle assertions** on Julia 1.13.1: 48 ordinary, six cold
and three fixed-grid population-size cases. Every case checks scratch/particle
preservation and independent final rounding/total bounds before timing.
Ordinary 5,000-particle cases allocate 3,528–6,088 B per corrected call versus
0 B in the baseline. Shared-host median ratios range 0.505–1.966; mixed-load
variation prevents a portable speed claim. Complete cold calls allocate
179,272–5,136,504 B versus 164,920–4,649,400 B, with ratios 0.710–2.262.
These are aggregate allocations, not peak resident memory.

For the same eight-cell grid and equal weights, 16/256/4,096 particles require
3,728–3,784 B per call. The 16-particle case takes 7.208 microseconds versus
0.875 microseconds in the baseline (shared-host medians). This scalar correction
cost remains an optimization target; no zero-allocation or maximum-performance
claim applies to corrected capture. The final production source still needs its
own gates on 1.12.7 and the newly released 1.13.2.
Driver: `/tmp/hybrid-charge-conservation-v6-resources.jl`; results:
`/tmp/hybrid-charge-reconcile-v5-root-resources113.{log,toml}`, exit 0.

### Final v6 wave and resource gates

The final candidate starts from pushed `6ad7276`; its full-PIC SHA256 is
`7390c5803abe7148e056d577aa1295ab64c757754d651c2d57cbe0b9fc50e6aa`.
Native Julia **1.12.7 and 1.13.2** each pass the unchanged 1,000-step wave's
five assertions. The frequency error is 0.013636642520216225 and continuity
residual is 2.877175015333387e-12, below the original 0.03/1e-8 thresholds.
The final-step exact summed-density drift is
1/5070602400912917605986812821504. Both current-revision baselines reproduce
the original failed residual of 9.889422779356765e-7 and drift 115/2^53.
Evidence: `/tmp/hybrid-charge-conservation-v6-{base-wave,wave}{112,1132}.log`.
The candidate wave commands succeeded before their sequential integration
commands began. Those complete integrations now pass **29,113 assertions per
version**, terminal exit 0, across the full `test_empic.jl`, `test_empic1d.jl`
and `test_io_electromagnetic.jl` files. The last includes exact checkpoint
continuation, changed timestep/subcycling and edited particle storage controls.
Logs: `/tmp/hybrid-charge-conservation-v6-integration{112,1132}.log`.

Each version also passes **504 complete-call resource/oracle assertions**,
with terminal exit 0. All 48 ordinary cases allocate 3,528–6,088 B versus 0 B
in actual `6ad7276`; six cold cases allocate 179,272–5,136,504 B versus
164,920–4,649,400 B. On the eight-cell equal-weight control, allocations remain
3,728–3,784 B at 16/256/4,096 particles. For 16 particles, the median call takes
6.709 microseconds versus 0.875 on 1.12.7, and 7.542 versus 1.000 on 1.13.2.
Ordinary-case median ratios span 0.597–2.000 and 0.185–3.650 respectively under
shared load; these observations support no portable timing claim. The small-call
cost remains an optimization item. Logs and structured measurements:
`/tmp/hybrid-charge-conservation-v6-resources{112,1132}.{log,toml}`.

Checkpoint inspection confirms `_empic_checkpoint_fields` includes both saved
charge arrays and `_restore_checkpoint_state!` restores the field bundle directly
through `_restore_checkpoint_fields!`; it does not recalculate charge on load.
Owners: `src/io/checkpoint.jl:722` and `:1025`. Exact new-version continuation
passes in the complete electromagnetic I/O suite above.


Final independent four-thread controls pass **2,273 assertions per version**,
with terminal exit 0 on Julia 1.12.7 and 1.13.2. Across 162 fixtures they cover
single-cell and odd grids through 3D, both EM families, Float32/64, all built-in
shapes, mobile/background charge, unequal populations, contiguous particle
views and extreme cancellation. Each fixture checks independent nearest-even
cell bounds, exact charge-error identity, non-DC cell differences, scratch and
particle preservation; 20 concurrent captures per fixture match the serial
result before and after a changed particle weight. Ambient precision is unchanged.
Logs: `/tmp/hybrid-charge-conservation-v6-concurrency{112,1132}.log`; driver:
`/tmp/hybrid-charge-conservation-v6-concurrency.jl`. The first minimum-version
launch used an unavailable `1.12.7` Juliaup channel and exited before Julia ran;
its log remains `...concurrency112-channel-failure.log`. The installed `1.12`
channel was checked to report 1.12.7 and used for the completed run. All 180
frozen Julia source/test hashes still match the candidate manifest.


The candidate's Julia 1.12.7 numerical/compiler/quality gate completes with
**13,268 assertions**, terminal exit 0: 1,953 normalization, 2,936 strengthened
raw/final charge, 892 current, 7,369 conservation, 106 concrete compiler,
three package-analysis and nine Aqua checks. All 22 package terminal diagnostics
match the byte-identical `6ad7276` baseline; the report ceiling is unchanged.
Logs: `/tmp/hybrid-charge-conservation-v6-quality112.log` and
`...diagnostics112.log`. The 1.13.2 candidate also passes **13,268 assertions**,
terminal exit 0. Every one of its 21 package terminal diagnostics matches the
same-runtime actual `6ad7276` baseline, which separately passes 4,577 assertions.
Logs: `/tmp/hybrid-charge-conservation-v6-quality1132.log`,
`...baseline-quality1132.log` and `...diagnostics1132.log`.

The four-file production/test patch was replayed on a fresh `6ad7276` copy;
all 180 resulting Julia source/test files match the frozen candidate byte for
byte. JuliaFormatter 1.0.62 is idempotent for all four changed Julia files.
Permanent tests retain the previous raw-stage assertions and add independent
capture rounding, exact target, geometry and state-preservation controls.
The README and changelog describe the changed capture contract and remaining
accuracy boundary. No diagnostic, tolerance, schema or dependency changes.


### Remaining local-accuracy boundary

**VERIFIED:** the DC contract does not certify local deposition accuracy under
strong ordinary cancellation. An independent 48-check witness on each of 1.12.7
and 1.13.2 uses three colocated NGP
electrons with weights `[2^precision(T), 1, -2^precision(T)]` for Float32/64 and
both EM families/background modes. The native raw charge is `[0, 0]` with no
range-loss flag; final capture restores the correct total -1 as `[-0.5, -0.5]`,
while exact local deposition is `[-1, 0]`. The uniform correction preserves the
raw non-DC error. This is an explicit accuracy limitation, not evidence that the
DC target or complete-shift rounding failed. Source and completed evidence:
`/tmp/hybrid-charge-conservation-v6-local-cancellation-boundary.jl` and
`...local-cancellation{112,1132}.log` (48 assertions each, exit 0).

Follow-up: evaluate compensated or binned local accumulation against independent
spatial moment oracles for signed and unequal weights. Define its changed rounding
contract separately from generic raw deposition, then require coupled accuracy,
compiler, threaded/reduction and measured complete-call cost gates. A global
correction alone cannot supply the missing non-DC information in this witness.
