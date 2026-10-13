# Gyro averaging correctness and cost

Status: geometry and streaming-mean repair verified on Julia 1.12.7 and 1.13.1.
The measured ordinary-path slowdown remains open.

## Verified defects and contract

The baseline Float16 mean of constant one is 0.5 at 8,192 samples and zero at
70,000. At 12,000 points the phase overflows and sampling throws after 10,428
callbacks. Tiny magnetic components and retained BigFloat cancellation also
lose representable coordinates or means. Evidence:
`/tmp/gyro_ring_baseline_probe112_fixed.log` and
`/tmp/gyro_ring_retained_geometry_probe112.log`.

Form each phase from a bounded fraction. Retain radius-weighted basis terms from
raw magnetic components through the complete center-plus-displacement sum.
Use a certified expansion for the running sample sum and switch once to exact
arithmetic when its certificate fails. The shared conversion owner rounds the
final mean. Samples are neither stored nor replayed.

Successful calls evaluate the callback exactly once per point in phase order,
including zero radius. Reject nonfinite coordinates before that callback.
Invalid samples and callback exceptions stop sampling and retain earlier callback
effects. Public drift and guiding-centre update arithmetic remain unchanged.

## Verification

Both versions pass **268,187 assertions across 19 testsets**, including 2,939
fusion assertions, 21,366 exact-sum/mean checks and 25,053 independent physical
geometry checks. Exact rounding-cell tests are independent of the production
converter. The physical basis oracle uses 8,192-bit arithmetic for extreme
planes. Callback boundaries and a quadratic-field analytic average are covered.
Logs: `/tmp/gyro_ring_frozen_regressions112.log` and
`/tmp/gyro_ring_frozen_regressions113.log`.

Each version passes 118 concrete JET calls with zero reports and 135 runtime
assertions. The Julia 1.12 package gate retains 24 reports, matching the baseline
after path/line normalization; the ceiling remains 30. Logs:
`/tmp/gyro_ring_concrete_jet112.log`, `/tmp/gyro_ring_concrete_jet113.log` and
`/tmp/gyro_ring_package_jet112.log`.

## Measured cost and next implementation

Seven alternating local batches compare the same callbacks against the committed
baseline. Ordinary Float32/64 calls at 16 points remain allocation-free, but
measured medians are 11.2–18.8 μs with paired ratios **3.78–20.0× slower**.
Variation is substantial on the shared machine; these are not portable limits.
Log: `/tmp/gyro_ring_performance112.log`.

Zero-radius constant means at 8,192 samples allocate 0 B; exact division at
12,000/70,000 allocates 3,400–3,752 B. Extreme geometry and uncertain means use
38,992–57,056 B and 4,424–8,696 B respectively. Bounded stored-BigFloat cases
use 282,760–299,352 B. Arbitrary exponent-span resource bounds remain unverified.

Next, reduce the cost of ordinary point construction and certified accumulation
while retaining all range, cancellation, callback and rounding tests. Measure
components with every output consumed; the existing geometry-only timing can
omit work through dead-code elimination and cannot establish complete setup
cost. Require paired end-to-end improvement, unchanged allocation behavior and
both Julia gates before claiming optimization. Public drift costs remain a
separate follow-up.
