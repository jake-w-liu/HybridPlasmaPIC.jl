# Oblique Hall-MHD production validation

CAP-007B3 remains open beyond this bounded periodic CPU Float64 fluid slice.

Landing status (2026-10-26, on `cda63f1`): `test_oblique_hall_history.jl` and the
neighboring parallel-wave, mode-history and validation-frequency suites pass on
Julia 1.12.7 and 1.13.2. The earlier report that a fixed memory limit failed on
Julia 1.12.7 does not reproduce: the landed test and `validation/oblique_hall_history.jl`
contain no resource assertion, and the test passes on both runtimes. Case 35 is
a long campaign (60 profiles; 46 completed in about 25 minutes on a shared
Apple M5), so it is a manual validation run rather than part of `Pkg.test`.
The scientific matrix and gates were frozen before production measurements;
the tracker below records subsequent evidence. Base: `839c1098df6bcb444ffd1a958f4bc8ad689c8460`.
The independent worktree is `HybridPlasmaPIC-oblique-hall-20261025`; evidence is
`HybridPlasmaPIC-audit-evidence/2026-10-25/oblique-hall-v1`.

## Pre-edit brief and architecture

- Implement missing oblique production histories in new case 35, using the
  unchanged `HallMHDState`, `hall_mhd_rhs!`, `step_hall_mhd!`, Fourier diagnostic
  and `validation/mode_history.jl` fitter. Do not modify case 34, its helper,
  the canonical solver or shared roadmap documents.
- The equations in `src/models/hall_mhd.jl` evolve density, bulk velocity and B.
  The momentum RHS includes both `Ti*gradient(n)` and the electron pressure
  gradient prepared by `src/electrons/ohms_law.jl`. The actual
  `IsothermalElectrons` closure in `src/electrons/closures.jl` is `pe=Te*n`.
  Thus the linear total-pressure coefficient is exactly `cs²=Ti+Te` at n0=1.
  Scalar pressure has zero curl in the linear induction equation about uniform
  density. Evolving-electron, anisotropic and kinetic closures are excluded.
- Reuse `test/oracles/hybrid_dispersion_oracle.jl` as the independent eigenproblem
  owner. It contains no production spatial operators. Derive its density/vector
  mapping separately below and test it against a separately written reduced
  matrix and central differences of the actual production RHS. Do not duplicate
  a production solver or fit simulated theory in place of measured data.
- The initializer owns a deterministic physical eigenvector. The sampler owns
  actual seven-component complex Fourier histories from production fields.
  The fitter receives only a measured component and dt; eigenvalues never enter
  it. The case owns fixed profiles, acceptance metrics, raw CSVs and receipts.
- Preserve errors visibly. Invalid configurations, degenerate branch seeds,
  ambiguous tracking and incomplete runs fail. At the case boundary, record a
  thrown exception and existing raw-file hashes in the receipt, then rethrow.
  This catch documents failure and never returns a passing/default result.
  Stream each successful sample to its raw CSV so a failed evolution retains
  the completed prefix. Use fresh artifact directories for evidence campaigns.

## Independent equations and normalization

Use `exp(i*k*x-i*omega*t)`, n0=|B0|=1, u0=0,
`B0=(cos(theta),sin(theta),0)`, and k>0. The angle is between k and B, not merely
between a parallel wave and a coordinate axis. This first production matrix is
one spatial dimension with three vector components; it does not establish
multidimensional spatial-derivative convergence.

For `q=(delta_n,u_x,u_y,u_z,b_y,b_z)`, let c=cos(theta), s=sin(theta) and C=Ti+Te.
The independently reduced time-generator equations are:

```text
dn/dt   = -i*k*u_x
du_x/dt = -i*k*(C*delta_n + s*b_y)
du_y/dt =  i*k*c*b_y
du_z/dt =  i*k*c*b_z
db_y/dt =  i*k*(c*u_y-s*u_x) - k²*c*b_z
db_z/dt =  i*k*c*u_z          + k²*c*b_y
```

The oracle uses `(u_x,u_y,u_z,b_x,b_y,b_z,delta_p)` and eigenvalue omega.
For every nonzero propagating root, recover `delta_n=k*u_x/omega` and verify
`delta_p=C*delta_n`, `b_x=0`, and the reduced residual independently. This remains
valid at C=0 without dividing by pressure. Require real ideal frequencies and
reject zero/near-zero roots before the density mapping.

Normalize the physical seven-vector `(delta_n,u_x,u_y,u_z,b_x,b_y,b_z)` to unit
Euclidean norm and fix its largest component real positive. Initialize each
real field with `amplitude*real(v*exp(i*k*x))` plus its uniform background.
Record `2*mode_amplitude(field,g,(mode,))/n` for every component: the factor two
accounts for the conjugate spatial mode of a real field. Use the initially
largest component for the single-exponential fit; save every component so the
measurement and polarization can be independently replayed. The initial mode
norm is the declared amplitude; density remains positive at the declared
amplitudes. No static mode is called a propagating wave.

## Branch identity frozen before measurements

At theta=0, the positive cold electromagnetic seed frequencies are
`(sqrt(k^4+4*k²)+k²)/2` (whistler) and
`(sqrt(k^4+4*k²)-k²)/2` (ion-cyclotron). Warm acoustic seeds are `k*sqrt(C)`.
Use their signed counterparts for negative-frequency branches. Cold acoustic
and zero-frequency roots are excluded explicitly. Reject coincident seed
frequencies within `1e-6*max(1,maximum(abs,roots))`; never start tracking from an
arbitrary degenerate basis. This also rejects warm acoustic/electromagnetic
crossings at the seed instead of relabeling them silently.

Continue independently for each sign from 0 degrees to the requested angle in
fixed 0.1-degree increments (within the originally frozen maximum of 1 degree). At each point, choose the one-to-one
assignment of normalized physical eigenvectors that maximizes total absolute
inner-product overlap. Require each matched overlap >=0.8, assignment-score
gap >1e-3, and root separation above the seed bound. Fix phases by overlap,
retain every continuation node and its frequency/overlap/gap, and only then
select the named branch. Supported angles are finite and in [0,90) degrees;
near-degenerate inputs may fail these explicit admissibility checks.

## Frozen production matrix and gates

Run 60 profiles: angles 30 and 60 degrees; temperatures (Ti,Te)=(0,0) and
(0.3,0.2); two cold branches or three warm branches; both frequency signs;
three independent resolutions. Every profile uses mode 1, L=2pi, amplitude
1e-5, duration 80, eta=etaH=0, nfloor=1e-6 and final magnetic projection enabled.
Baseline is n16/dt0.01; grid refinement is n32/dt0.01; time refinement is
n16/dt0.005. These steps resolve the retained highest-grid Hall frequencies
within the imaginary-axis RK4 stability interval. This is not a general CFL
estimator or a latency benchmark.

Fit the final three quarters, retaining exact start/end indices and every raw
sample. Require, per profile:

- Signed relative frequency error <=1e-3 and absolute growth/damping error <=1e-4.
- Phase and log-magnitude RMS <=1e-3; half-window frequency/rate differences <=1e-4.
- At least two measured periods; zero/static fits cannot meet this gate.
- Maximum full seven-vector history error relative to its initial norm <=1e-3
  against the independently initialized linear eigenmode. Finite-amplitude
  compressive nonlinear effects are included in this error, not subtracted.
- Maximum polarization leakage relative to measured mode norm <=1e-3.
- Relative periodic mass error <=1e-12 and maximum varying longitudinal B <=1e-12.

For each branch/angle/temperature/sign, the grid- and time-refined signed
frequency differences from baseline must be <=1e-4; rate differences <=1e-5.
Report raw errors as well; do not claim an order from a spectral/roundoff floor.

Four additional timestep-order families use mode 2, angle 30 degrees, whistler,
both temperatures and signs, n16, amplitude1e-6, duration1, and dt0.04/0.02/0.01.
Compare the full measured vector at every time to independent linear evolution;
require both consecutive orders >=3.5. No post-measurement exclusions or relaxed
limits are permitted. These 12 raw histories are separate from the fit matrix.

## Verification and tracker

- [x] Read actual Hall/closure/operator owners, independent oracle, case 34,
  fitter, callers, test registration and CI; establish pre-edit contract.
- [x] Implement new validation helper/case/plotter and local case README.
- [x] Verify reduced matrix, cold/warm closed-form limits, nondegenerate seeds,
  branch identity, eigenvector normalization and independent production-RHS
  central differences, including held-out angles, mode numbers and temperature
  splits. Confirm no pressure term is omitted or inserted into induction.
- [x] Test actual short histories, both signs, raw prefix/failure receipts,
  invalid configuration, static-mode rejection and artifact replay.
- [x] Run the frozen 60-profile and 12-history matrix on Julia 1.12.7/1.13.2;
  retain all failures and reverify unchanged gates after any confirmed repair.
- [x] Check neighboring Hall, oracle, fitter and case-34 tests, reviewed compiler
  output, package sandbox imports, formatting, CLI discovery and final diff.
- [x] Generate all four PDFs on each runtime; check PDF structure/text/bounds.
  No visual image inspection was performed.
- [x] Verify complete-call allocation slope, stable samples and zero measured
  compilation; preserve original failed measurement setups.
- [ ] Pass the original retention gate on both runtimes: final Julia 1.12.7
  RSS growth remains above 16 MiB. Native diagnostic plateaus do not pass it.
- [x] Freeze source/dependency/run/artifact hashes and separate patch/handoff.
- [ ] Parent review, individual main push and exact-commit CI.

Wrong patches to reject: reversing Faraday/Hall or frequency signs; omitting
Te from momentum; dividing cold delta_p by zero; selecting a sorted eigenvalue
without continuation; choosing an arbitrary degenerate seed; using theoretical
rates as measurements; normalizing real Fourier modes without the factor two;
accepting static/multimode histories; or dropping failed raw data.

This slice supplies deterministic oblique fluid histories and numerical gamma
checks, not kinetic damping/growth, anisotropic/evolving electron pressure,
oblique hybrid PIC, mass-ratio scans, nonperiodic/GPU/distributed behavior or
production-to-NHDS convergence. Those CAP-007B3/B4 requirements remain open.

## Retained reference-path failure

The initial independent test fails visibly at held-out k=0.7, C=0.5: distinct
parallel acoustic/ion-cyclotron seeds are separated by 0.001661956684. A 1-degree
jump has matched overlap 0.7162344503, below the unchanged 0.8 gate. The retained
reference-only probe measures 0.9702401424 for 0.1 degrees. The continuation now
uses fixed 0.1-degree nodes; assignment, root-gap and overlap gates are unchanged.
The failed source and log remain under `initial-source` and `test-initial-112.log`;
no full production campaign preceded this repair.

The fixed 0.1-degree reference path passes all 2,905 permanent assertions on both
Julia 1.12.7 and 1.13.2 before the full campaign, including reduced-matrix,
production central-RHS, normalization, signed short-history and failure-receipt
checks. This is initialization resolution, not production-solver convergence.

Resource checks are separate from latency: after eight warmups and three
allocation samples, measure complete calls of 8, 1,000 and 2,000 steps at n16.
Require incremental allocation `(bytes2000-bytes8)/(bytes1000-bytes8) <= 2.2`
with stable measured bytes and no measured compilation/recompilation. After
20 additional short calls and full garbage collection, require Julia live-byte
growth <=1 MiB and process resident-memory growth <=16 MiB from the warmed
baseline. These are bounded retention checks, not exact transient native peaks.

## Final V5 evidence and unresolved resource acceptance

This candidate is **not ready for promotion**: the original RSS gate fails on
Julia 1.12.7. Scientific, inference and allocation passes remain separate from
that failed resource acceptance. No latency acceptance was run.

- `campaign-v5` runs the actual case CLI on Julia 1.12.7 and 1.13.2: all 60
  profiles, 12 order histories and 748 metrics pass. The independent NumPy
  reduced-matrix/raw-history replay verifies every metric policy, archived
  source/environment hash, and rejects 16 mutations per runtime. All 134
  scientific artifacts and initial eigenvectors per runtime are bitwise
  identical to the original V1 campaign after the inference/allocation repair.
- `sandbox-v5` runs the actual Pkg test target: all 2,905 permanent assertions
  pass on each runtime, including inferred reference returns. `quality-v7`
  passes the neighboring Hall, actual PIC dispersion, fitter/case-34, Aqua and
  package compiler gates. Full package JET retains 29/27 reports; these are not
  zero-report claims. Minimal stdlib controls reproduce the concrete norm
  reports and the Julia 1.12 view-alias fallback; full trees and exact checked
  correspondences remain in the logs and `COMPILER_REVIEW_V5.md`.
- `formatting-v5` checks all five changed Julia files with the pinned formatter
  on both runtimes. `plots-v1` contains four PDFs per runtime, with successful
  PDF/text/bounds checks. Final plotting inputs are bitwise unchanged from V1.
  There was no visual image inspection or latency measurement.
- `binding-repair-v5` compares 528 complete references and 20 actual history
  objects exactly on each runtime. Separating a captured seed binding, removing
  three transient phase-column copies with views, and reusing the checked
  branch index reduces ordinary reference allocation from 8,339,536 B to
  7,416,480 B on Julia 1.12.7 and to 7,387,552 B on Julia 1.13.2. This meets the
  predeclared 7.5 MB target. The binding-only 7,891,680 B miss and earlier failed
  harness/source versions remain retained. Norm arithmetic is unchanged.
- `resource-v6` retains all work and samples, with full GC after the original
  eight warmups so deferred FFTW finalizers compile during setup. The prior
  measured-finalizer-compilation failure remains under `resource-v5`, with its
  concrete methods traced in `resource-compile-trace-v1`. Final samples have
  stable bytes and zero compile/recompile; allocation slope is 2.0081104 and
  2.0081170, below the original 2.2 gate. Julia live deltas are -80 B and 384 B.
  RSS deltas are **20,398,080 B (FAIL)** and 15,286,272 B (pass), respectively.
- The original `resource-v1` RSS growth of 22,183,936 B also failed on Julia
  1.12.7. Its result file incorrectly says pass because the original harness
  wrote status before a failing testset returned. Its exit 1 and log are
  authoritative. Corrected receipts write pass only after success.
- `native-ownership-v3` records five fixed 20-call windows per runtime and
  passing retained 64 MiB allocation/free sensitivity controls. Malloc-live
  endpoints plateau from windows 2 through 5; RSS/compression vary. This is
  diagnostic evidence only, not proof of no leaks or a cause for the failed
  RSS gate. `RESOURCE_REVIEW_V5.md` keeps the bounded next investigation open.

The final seven-file patch and hashes are in the evidence directory's
`handoff-v5`. Numerical sources match the completed V5 runs; this final tracker
update records their results after source-pinned runs finished. Canonical solver,
case 34 and the shared fitter remain unchanged. CAP-007B3 remains open beyond the
bounded oblique fluid scope, and parent integration/push/CI remain unperformed.
