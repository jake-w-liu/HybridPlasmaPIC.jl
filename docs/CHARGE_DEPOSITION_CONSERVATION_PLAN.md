# Charge-deposition conservation

Status: VALID-001 cause verified on Julia 1.12.7 and 1.13.1; repair design under
investigation. No deposition change or validation-tolerance change is implemented.

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
The broader numerator cases are not yet measured by this four-case probe.
The inspected callers include `momentum!` and `current!` in
`src/coupling/moments.jl`, and both general/legacy `_charge_density!` paths in
`src/models/full_pic.jl`. EMPIC validates electron/ion charges as −1/+1, so its additional range checks
must cover cancellation between species or against the background. The public
`current!` also multiplies by the species charge after materializing momentum.
Those specific extreme-input outcomes have not yet been measured here.
