# Gyro averaging correctness and cost

Status: geometry and streaming-mean repair pushed as `de6151d`. The ordinary
ring optimization below is verified on Julia 1.12.7 and 1.13.1 and included
in this commit. Further measured cost reduction remains open.

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

## Correctness-repair cost

Seven alternating local batches compare the same callbacks against the committed
baseline. Ordinary Float32/64 calls at 16 points remain allocation-free, but
measured medians are 11.2–18.8 μs with paired ratios **3.78–20.0× slower**.
Variation is substantial on the shared machine; these are not portable limits.
Log: `/tmp/gyro_ring_performance112.log`.

Zero-radius constant means at 8,192 samples allocate 0 B; exact division at
12,000/70,000 allocates 3,400–3,752 B. Extreme geometry and uncertain means use
38,992–57,056 B and 4,424–8,696 B respectively. Bounded stored-BigFloat cases
use 282,760–299,352 B. Arbitrary exponent-span resource bounds remain unverified.

## Verified ordinary ring optimization

Cache lossless binary alignment of the center and radius-weighted basis terms
once per ring. Eligible Float32/64 points use compensated phase products and
sum residuals. Uncertain cancellation and range boundaries retain the previous
complete-expression calculation; Float16 and BigFloat retain their prior path.
The canonical physical basis, callback order and exact mean of the floating
sample stream remain unchanged. This is not a universal bitwise-coordinate or
correctly-rounded-coordinate claim.

Both versions pass **269,015 assertions across 19 testsets**, including all
retained force, position, drift, ring and BigFloat oracles and 828 new ring
checks. Each version also passes 70,308 independent geometry/alignment checks,
124 concrete JET calls with zero reports, and 141 runtime assertions. Julia
1.12 package JET retains 24 reports; all 221 normalized report lines match the
committed baseline. Logs:
`/tmp/gyro_ring_opt_final_regressions_v3_112.log`,
`/tmp/gyro_ring_opt_final_regressions_v3_113.log`,
`/tmp/gyro_ring_opt_invariant_compiler112.log`,
`/tmp/gyro_ring_opt_final_compiler113.log`.
Final source SHA256: `4d5c733b8937324f0450af0d1f008a9308a62e3246cb719edd6853005da52f96`.

Seven alternating paired local batches compare against `de6151d`. Default
16-point constant/quadratic cases take **18–37% less time**, with medians of
3.15–4.51 μs and unchanged 0 B allocation. Evolving quadratic inputs take about
32% less time; 128-point paired time ratios are 0.565–0.623, also 0 B.
Three-point results are mixed (new/old ratios 0.777–1.211), retaining
2,184–2,840 B allocations. The shared machine was not isolated.

Fully consumed setup/point/callback measurements improve from 3.62 to 2.31 μs
for Float32 and 3.40 to 2.47 μs for Float64, with 0 B. First evaluated calls,
including compilation, cost 1.01/1.39 s and 58.9/189.7 MB for Float32/64 versus
0.65/0.87 s and 47.9/114.1 MB in the baseline. Candidate-first ordering prevents
an order-neutral startup comparison. Large-count and exceptional paths retain
their allocations, with mixed timings; uncertain-mean measurements include
slowdowns despite unchanged mean arithmetic. Retained BigFloat cases allocate
284,040–300,632 B. Raw measurements and interquartile ranges are in
`/tmp/gyro_ring_opt_final_resources112.log`.

An initial allocation test captured an enclosing mutable binding and measured
144 B on both baseline and candidate. Explicit local arguments correct the
fixture; the strict 0 B requirement remains. An explicit cache-presence check
inside the tuple closure removes one new abstract compiler finding without
changing arithmetic. Failed logs remain in the verification record:
`/tmp/gyro_ring_opt_final_handoff.md`.

These measurements establish improvement over the correctness repair for the
stated cases. Pre-repair speed, uniform speedup, arbitrary BigFloat exponent-span
resource bounds and optimal performance remain unverified. Public drift and
further ring cost reduction remain separate follow-ups.
