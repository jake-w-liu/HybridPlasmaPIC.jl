# Production wave-history validation

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

Recorded status: the first CAP-007B slice passes local correctness, resource, compiler,
artifact-replay and CLI/plot checks on Julia 1.12.7 and 1.13.2. Pushed as
`839c109`; [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40194825415)
passes five jobs but fails twelve Julia 1.12 log-amplitude residual assertions.
The [centering repair](MODE_HISTORY_CENTERING_PLAN.md) preserves those failures
and the unchanged tolerance. Corrected CI and full compiler review remain.
The full [completeness roadmap](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md)
remains open; this plan does not replace its other solver and parameter regimes.

## Pre-edit brief

- Mode: implement missing production time-series measurement and comparison.
- Requirement: measure signed frequency and growth/damping from actual evolved
  complex Fourier coefficients. Preserve raw histories, fit windows, diagnostics,
  refinement results and failures. Theory may initialize a mode and supply the
  comparison, but must not enter the fitted frequency or rate.
- Baseline: `27372b2`; case 27 compares an independent Hall-MHD eigenvalue oracle
  to NHDS. Cases 18/20 and `test_dispersion.jl` run production solvers but measure
  unsigned frequencies without this growth/fit-window assessment.
- Owners: new `validation/mode_history.jl` owns the single-exponential fit;
  case 34 will own the production run matrix and artifacts. Existing time
  steppers, Fourier-mode diagnostic and external-reference cache remain canonical.
- Contract: fit `a(t)=A exp((gamma-i*omega)t)` over a caller-selected, uniformly
  sampled complex history. Inputs must be finite, nonzero and representable in
  Float64; at least eight samples allow two disjoint half-window fits. Reject
  invalid cadence and adjacent wrapped phase increments at least pi/2. This
  is a sampling guard, not a proof of absence of aliasing. Refinement must
  resolve alias ambiguity. A constant nonzero history has omega=gamma=0.
- Fit: regress unwrapped phase and log magnitude independently on centered,
  normalized time. Report signed omega, gamma, residual RMS, conditional OLS
  slope standard errors, and half-window slope differences. Standard errors
  assume independent constant-variance residuals; correlated PIC histories
  require independent seeds/refinement and cannot inherit that confidence claim.
- Compatibility: validation-only addition; no package export, runtime solver,
  checkpoint or existing estimator contract changes. The common validation
  loader delegates to this owner. No extra dependencies or swallowed failures.
- Verification: synthetic exact exponentials with both signs, growth/damping,
  phases, magnitude ranges, views and shifted windows; an independent matrix
  least-squares oracle; invalid inputs; phase alias boundary and multimode/chirp
  controls. Preserve input histories, check concrete compiler calls, allocation
  scaling, actual Pkg test sandboxes and full repository formatting. Production
  runs and their fixed gates follow before feature acceptance.

OLS residual scaling follows [NIST's least-squares definition](https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd431.htm),
checked on 2026-10-25. This reference supplies the conditional regression
statistic, not a plasma-noise model or numerical-discretization uncertainty.

## Staged implementation and fixed gates

- [x] **007B1:** implement and independently verify the fit owner. Reject wrong
  patches that return theory, erase the frequency sign, discard rate information,
  hide zero/nonfinite samples, or report a multimode/chirped history as a clean
  single exponential. The caller gates residuals/window sensitivity explicitly.
- [x] **007B2:** case 34 evolves circularly polarized parallel waves with
  `step_hall_mhd!` and `HybridStepper`. Set n0=B0=1, zero ion/electron temperature,
  amplitude 1e-4, modes 1 and 2 in L=2pi. Initialize both eigenbranches of the
  independent two-component matrix `[[0,ik],[ik,ik^2-eta*k^2]]` in (u,B).
  For Hall-MHD use eta=0 and 0.03; hybrid uses eta=0 because its resistive force
  closure differs from the fluid momentum equation. Use the final three quarters
  of a fixed-duration history, save every sample and compare both half-windows.
  No expected eigenvalue is passed to the fit owner.
- [x] Freeze the exact run matrix before solver measurements: baseline n16,
  dt0.01, particles/cell16, duration 40; refine n32, dt0.005 and particles/cell32
  independently. Use eight magnetic subcycles for hybrid. Require relative
  frequency error <=0.03, absolute rate error <=0.01, phase/log-amplitude RMS
  <=0.03, and half-window frequency/rate differences <=0.01. Require at least
  two measured periods in the fit window. Record each refinement's error and
  noise diagnostics; no automatic convergence-order claim from a noise floor.
- [x] Add timestep order studies against actual full-history errors: use the
  mode-2 whistler, n16, particles/cell16, eta0, amplitude1e-6, duration1 and
  dt0.04/0.02/0.01 for both solvers. Compare every complex sample with the
  independent continuum (Hall) or discrete CIC (hybrid) matrix exponential.
  Require both consecutive observed orders >=3.5 (Hall RK4) or >=1.5 (hybrid).
  The higher-frequency branch keeps the RK4 truncation error measurable at
  these stable steps. There are no post-measurement roundoff exclusions.
- [x] Add ordinary test registration, CLI/discovery integration, README,
  per-case plotter, raw CSVs and a source/environment/run receipt. Both Julia
  versions pass local checks. This first matrix is periodic CPU Float64;
  unsupported backend/model claims remain open.
- [x] Push the locally accepted slice individually to main as `839c109`.
- [ ] Complete all six exact-SHA CI jobs plus full compiler-tree review.
- [ ] **007B3:** extend to oblique fluid and cold hybrid modes, electron/ion
  temperature and mass-ratio scans, production EM branches and physical
  damping/growth with independent kinetic references. Freeze initialization,
  fitting windows, seeds and branch tracking for each regime before running.
- [ ] **007B4:** compare converged production frequencies/rates with pinned NHDS
  outputs where model assumptions match. Propagate independent-seed statistics,
  finite-window variation and discretization bias separately. Keep case27's
  eigenvalue-only comparison labeled as such.

The first parallel circular wave has transverse velocity and constant magnetic
magnitude; its cold fluid equations reduce to the stated two-mode matrix. The
existing independent Hall test checks the matrix exponential against the full
solver (`test/test_hall_mhd.jl`, helical-wave test). This source fact supports the
chosen oracle; the completed campaign is checked against this oracle below.

## Initial production trial and discrete initialization

The unchanged baseline matrix was run first on Julia 1.12.7. All eight frequency
errors are below 3%, but the mode-2 hybrid ion-cyclotron history exceeds both
fixed residual limits: phase RMS 0.0341386504 and log-amplitude RMS 0.0343242399
versus 0.03. `pilot-v1-112.toml`, its raw histories, source hashes and log retain
this failed gate; limits remain unchanged. The first fit driver also had a
wrong checkout path; its two failed launches are retained separately from the
corrected 2,634-assertion-per-runtime fit tests.

The next initialization follows the discrete cold CIC equations, without
changing the production solver or fitted result. Let `G_p` be interpolation of
the unit grid mode at particle p, and `Q=mean(abs2(G_p))` for the uniform lattice.
Grid deposition is the adjoint interpolation, so for the grid velocity Fourier
coefficient U and grid magnetic coefficient b:

```text
dU/dt = i(Q-1)U + i*k*Q*b
db/dt = i*k*U + i*k^2*b
v_p(t=0) = U*G_p/Q
```

This follows by depositing `dv_p/dt = i*G_p*U + i*k*G_p*b - i*v_p`.
The continuum matrix is recovered at Q=1. Initialize a single eigenvector of
this discrete matrix and retain the continuum eigenvalue as the comparison
reference. Q comes from the actual lattice/CIC basis, never from a fitted
frequency. Verification must check adjoint normalization and the initial
linearized force/induction independently. The independent initializer regression now checks actual deposited moments and
production-gathered Lorentz forces: 6,748 assertions, plus 37 configuration and
history checks, pass on each Julia version.

The second source-pinned trial runs all four hybrid histories and meets the
unchanged residual limits: worst phase RMS 2.191e-5 and log-amplitude RMS 2.203e-5.
The worst frequency error remains 1.758%, since the continuum comparison was
not replaced by the discrete initialization frequency. Independent initial-force
checks and the full two-runtime campaign subsequently pass. Both trial source
snapshots are retained and verified against their original recorded hashes.

The full comparison matrix contains 40 profiles: two branches, two modes and
baseline/grid/time refinements for each Hall resistivity (24); the same two
branches/modes with baseline/grid/time/particle refinements for hybrid (16).
Six additional histories implement the fixed timestep-order studies above.

## Full campaign and fit refinement

The first frozen 40-profile/six-order campaign passes all 324 metrics on both
Julia 1.12.7 and 1.13.2. Source snapshots, commands, logs, raw histories and receipt
hashes are retained in `production-dispersion/campaign-v1-*` in the audit evidence.
The final integration and plot results are recorded below.

A separate 10,000-sample reversed history fails the original fitter's 128eps
frequency gate on both runtimes: -1.0000000000000868 versus -1. Phase accumulation
is now compensated; the original accuracy gate is unchanged and the case joins
the permanent suite. The original resource estimate 16n+2048 also fails four
profiles per runtime. VERIFIED: constructing only two bare Float64 vectors at
n=10,000 already allocates 163,968 bytes versus 160,000 logical payload bytes
(the original limit 162,048). The revised resource contract compares against
those same bare work arrays plus 512 bytes, preserving bounded linear workspace
while accounting for runtime allocation size. Both failed logs remain retained;
this correction is not described as a measured memory optimization.


The formatted candidate passes both actual Pkg test sandboxes: 3,153 existing
frequency checks, 2,642 fit checks, 6,785 initializer/configuration checks and
nine Aqua checks per runtime. Discovery, the full 324-metric matrix and all
four concrete compiler plus 64 resource/accuracy assertions also pass each.
An independent artifact review then finds all 40 fit CSV rows replace configured
run duration 40 with fit duration 30 during NamedTuple merging. The correction
preserves `duration` and adds `fit_duration`; no histories, fits or physics gates
change. Both prior campaigns and the 40-row-per-runtime witness are retained.

Both shared plotting environments resolve PlotlySupply 1.9.3, while the existing
suite requires 1.8.0; both default plotting commands fail that explicit version
gate. A separate pinned 1.8.0 environment now exports all four PDFs on Julia 1.12.7.
The corrected-metadata replay and pinned plotting subsequently pass on both runtimes.
No shared environment or version gate was changed.


## Final local acceptance

The corrected, formatted source passes the complete 324-metric campaign on
both runtimes. An independent Python replay checks every raw history with a
separate phase-unwrapping/least-squares implementation and a closed quadratic
reference. It verifies 144 source/environment/data hashes per runtime and rejects
200 changed-rate, changed-residual and duration-metadata records per runtime.
The observed Hall orders are 3.99959 and 3.99990; hybrid orders are 2.00163 and
1.99974. All limits are unchanged. Both pinned plotting commands export four
single-page PDFs. Rendered previews show readable titles, axes and traces;
the PDF text extractor's 0.056-point font-box overhang does not crop the visible
frequency label. Original failed checks and label probe remain retained.

Evidence under `HybridPlasmaPIC-audit-evidence/2026-10-23/production-dispersion/`:

| Check | Evidence | Result per runtime |
| --- | --- | --- |
| Actual package test sandbox | `final-v2-sandbox-*.log`, `final-v2-commands.json` | 3,153 existing frequency + 2,642 fit + 6,785 initialization + nine Aqua assertions |
| Concrete fit compiler/resource checks | `final-v2-quality-*.log/toml` | four concrete analyses; 64 bounded-storage/accuracy assertions |
| Final production matrix | `final-v3-campaign-*.log`, `final-v3-source.json` | 324 metrics; 40 profiles and six order histories |
| Independent raw-data replay | `independent-audit-final-v3-*.json`, `audit-campaign-v1.py` | all 46 histories, 144 hashes and 200 rejection controls |
| Plot export | `final-v3-plot-*.log` | four PDFs with pinned PlotlySupply 1.8.0 |
| Whole-tree formatting | `format-v2.log/toml` | 12 checks, including unchanged ASTs and two clean passes |

Package `src/`, extensions and dependency declarations are unchanged from
`27372b2`; the new validation helpers have separate concrete compiler checks.
This is local acceptance of the cold parallel validation slice. CAP-007B3/B4,
external matched production dispersion, complete solver/backend capabilities
and comparative performance evidence remain open.
