# PIC005 first slice — normalized continuity diagnostic

Status: diagnostic range repair verified on Julia 1.12.7 and 1.13.1. Current correction, Poisson/Gauss transform range and severe FFT cancellation accuracy remain open.

- **Mode:** Repair. Complete the existing public `charge_conservation_residual` expression on selected range-sensitive finite inputs, retaining its ordinary implementation and FFTW accuracy model.
- **Requirement and acceptance criteria:** Both EMPIC1D and general EMPIC D1–3 must evaluate the represented-mode `max(abs(Δrho + dt*div(J)))/max(abs(Δrho))` (or denominator1 when represented Δrho is zero) without first requiring the density difference, intermediate Fourier coefficients, inverse butterflies, standalone divergence, timestep product, or numerator/denominator maxima to fit in T. All original finite rho/J/dt operands survive until composition. True DC charge drift stays visible; only non-DC all-null effective wavevectors are excluded. Exact integer-root DFT witnesses cover both signs, finite and zero residuals, M±tiny bands, cancellation, all axes, coupled components, and finite ratios of unrepresentable maxima. Existing ordinary return values and allocations must remain measured-equivalent; exceptional costs and limits must be recorded. No claim of an exact arbitrary-grid DFT or complete PIC005 closure.
- **Baseline or reproduction evidence:** `/tmp/empic_boundary_audit_v2_112.log` is terminal0,256 witness assertions on Julia1.12.7: diagnostic inverse-density overflow, dt*div overflow followed by cancellation, divergence underflow restored by dt, and divergence overflow suppressed by dt, alongside the separate current cases. Driver and exact bounds: `/tmp/empic_boundary_audit_v2.jl`, `/tmp/empic_boundary_audit_handoff.md`. Current scalar-final source remains affected. `/tmp/pic005_raw_current_coupled_probe113.log` additionally verifies798 current-boundary assertions and8 exact-dyadic representation assertions on Julia1.13.1; these remain queued, not claimed repaired by this slice. Temporary complete diagnostic prototype `/tmp/pic005_dyadic_fft_prototype_v2.jl` passes1,211/1,211 on Julia1.12.7, including96 diagnostic cases and exact forward/inverse Gaussian-integer DFTs on D1–3. `/tmp/pic005_dyadic_fft_prototype_v2_112.log` prints the actual frozen module and EM hash. The identical prototype also passes1,211/1,211 on Julia1.13.1, terminal0: `/tmp/pic005_dyadic_fft_prototype_v2_113.log`.
- **Root cause or design gap (VERIFIED):** `full_pic.jl:1110–1137` and`:1731–1774` materialize Δrho in T, apply a raw forward/inverse transform, materialize a standalone derivative/divergence, then form the dt product, numerator and normalization in T. The installed AbstractFFTs inverse plan executes the unnormalized backward transform before multiplying by1/N (`definitions.jl:269,279–285`), so mathematically finite final data do not protect its butterflies. A one-mantissa normalizer also loses the tiny term of exact M−tiny; `/tmp/pic005_raw_current_coupled_probe113.log` verifies this separately. The existing finite-coefficient scalar projection owner cannot reconstruct prior transform loss.
- **Contract and invariants:** Keep original ordinary arithmetic and FFT calls when selection does not find a range-sensitive stage. Continue to call canonical `deriv!`/`divergence!` before the exceptional diagnostic path, preserving its finite-input, dynamic-range and multiplier-bound rejections. Keep existing acceptance of finite zero/negative diagnostic dt; no new public dt validation, and nonfinite dt/source operands stay on the existing path rather than entering Rational conversion. Diagnostic physical rho/J/E/B state remains unchanged; grid scratch is allowed to change and the grid is already non-thread-safe. No state/checkpoint/schema/cache/dependency/API addition. BigFloat is not used by production. This does not establish MPFR thread safety or support another hardware rounding/flush-to-zero mode. FFTW remains the transform: exact chunk decomposition and exact composition of its returned floats do not remove within-band FFT rounding/cancellation on general grids.
- **Canonical owner and affected surfaces:** Put one private normalized dyadic transform beside existing EM scaled spectral helpers in `full_pic.jl`; reuse `_em_power_of_two`, `_em_field_exponents`, `_em_scale_bounds`, shared `_exact_rational`/`_from_exact`, and `_spectral_mode_is_null`/`_spectral_wavevector`. Both public diagnostic overloads delegate their exceptional expression to one common owner over grid/rho/J/dt. Keep the previous ordinary wrapper loops, with range flags fused into existing density-copy and residual loops; reuse the existing exponent selector for canonical divergence restoration risk. Do not add a third TwoSum implementation or call gyro/AMR-private helpers. Root owns `em_field_energy` overloads and adjacent energy helpers in the same source; all PIC005 edits avoid them. Tests go in `test/test_empic.jl`, covering both families, so no include request is needed.
- **Compatibility and migration decision:** Public signatures, return type T, normalization, DC/Nyquist rules, ordinary operator order and errors remain. Selected exceptional finite-input results change from NaN/Inf/incorrect1 to the complete residual. Truly unrepresentable final normalized results retain floating Inf instead of being clamped; canonical rejection remains rejection. No changes to deposition, charge arrays, particle phase/clock, current outputs, ES/Gauss initialization or checkpoint data. All repository callers must be covered: `test_empic.jl`, `test_empic1d.jl`, electromagnetic I/O residual comparisons, and relevant validation examples (16 and20). The separate current transform failures remain open with their exact witnesses.
- **Implementation strategy:** Cold input delta is constructed exactly from original stored floats, then a finite dyadic radix expansion feeds the existing FFT plans. For every real/imaginary entry, denominator is a power of two; signed digits are grouped by binary weight and converted losslessly to T. Normalize each digit below1 and retain its exact power-of-two scale outside the transform. Choose an integer-derived width `max(1,p-ceil_log2(N)-2)`, using the existing integer ceil-log idiom (not floating log2). The two spare bits provide integer-sum headroom for the exact Gaussian-root acceptance fixtures when the bound is positive; no such exactness is promised when N is so large the width floors at1, or when twiddles are non-Gaussian. Stream one hardware work array and one exact accumulated output, skipping empty bands; do not retain all band arrays. After forward transforms, mask only the non-DC null density modes and compose `deltahat + dt*i*sum(k_d*Jhat_d)` exactly before normalized inverse transformation. Reconstruct represented density independently for the denominator, then compare exact absolute maxima and convert the final rational ratio once. All these modal operations are dyadic; non-dyadic current projection is a later slice with its own inverse-materialization contract. For lowest/highest nonzero bit positions L/H and width w, at most `floor(H/w)-floor(L/w)+1` bands exist. L/H are bounded by hardware input exponents and the finite number of k/dt products/FFT reconstructions; dense integer widths are likewise bounded for a fixed type and grid size. No data-dependent unbounded refinement/tolerance loop is introduced.
- **Verification strategy:** Before promotion, inspect the exact prototype arithmetic and run native-version focused oracles, selection checks, ordinary allocation/state comparison and a local consumed-output cost sample. Permanent tests must reject: (1) only checking nonfinite values (underflow-restored-by-dt witness); (2) one shared scale or one mantissa per exact coefficient (M±tiny and unequal-component cancellation); (3) materializing standalone divergence, numerator or denominator before the final ratio (four original boundaries); and (4) dropping DC to improve residuals (explicit charge-drift negative control). Add odd/non-power-of-two grids with independent high-precision DFT tolerances, preserving the stated floating-transform limit; check invalid input behavior, zero/negative dt, all-null/mixed Nyquist modes, scratch reuse and unchanged physical state. Native final-source focused/full EM-family affected-caller regressions, concrete and matched package JET, formatter1 idempotence and complete diff review are required. Final freezes must compose only owned full_pic hunks and test_empic against committed HEAD, excluding root's concurrent energy hunks and IO changes unless explicitly committed. No further broad suite until final source stabilizes.

## Prototype correction retained

The initial temporary width=p candidate passed the96 residual cases but failed12 strict arbitrary complex forward/inverse DFT checks: `/tmp/pic005_dyadic_fft_prototype112.log` (1,119pass/12fail). A p-bit digit leaves no spare mantissa capacity for butterfly sums. The independent child checked the actual v2 proof: Gaussian-root components sum N signed integer digits, and `w+ceil_log2(N)<=p−2` keeps those integer sums representable for the tested power-of-two shapes. v2 narrows the digit width, retains every original bit, and passes the same assertions without weakening their equality. This is finite-fixture and bit-capacity evidence, not a universal FFT accuracy proof.

## Earlier prototype observations

The prototype's maximum observed occupied band count is45 over340 transform samples, maximum represented bit span2,099; these are observations, not final upper bounds or production measurements. Final production measurements and gates are recorded below. Current correction and non-dyadic inverse policy remain an explicit subsequent slice.

The child independently rejected a fixed low-bit cutoff for later non-dyadic current coefficients: `/tmp/pic005_fixed_tail_counterexample.jl`, `/tmp/pic005_fixed_tail_counterexample112.log`,16/16 on Julia1.12.7. An omitted tail below qmin−p−2 changes both a nearest-neighbor midpoint and finite-versus-infinite overflow rounding in an exact DC inverse witness. This is a representation counterexample, not a newly reproduced production input or a proof that these denominators arise from admissible PIC state. It rules out a generic cutoff claim absent a reachable-value/separation argument. That is why this brief does not claim the non-dyadic current boundary solved.

## Final verification

The immutable `/tmp/pic005-diagnostic-prototype-rpzn57rx` snapshot is based on
`b725c71` and includes the pushed field-energy repair. Both native versions pass
**29,113 assertions** across complete EM, EM1D and electromagnetic I/O suites
(12m23.0s on 1.12.7; 11m52.0s on 1.13.1). The 900 new assertions cover independent
exact small-grid DFTs, range boundaries, error/state contracts and odd-grid
accuracy. Logs: `/tmp/pic005_final_v4_regression112.log` and `113.log`.

Final SHA256 values:

- `src/models/full_pic.jl`:
  `ba62bd1a5c060b847205a5329e8563089db5c2ea275ae82776bb5b1a55bc6a50`.
- `test/test_empic.jl`:
  `4fc4e2a5bfc12d346ad912e08820a19b394c543c911c10feb1ee784ae7617fdb`.

Both versions pass 33 concrete compiler assertions; package findings match the
baseline leaf messages exactly (26/25). Logs: `...v4_compiler112.log`, `113.log`
and `...v4_package_compare112.log`, `113.log`. The final iteration uses bounded
spatial-current iteration and an explicit exact-complex conversion type;
earlier unsuccessful type-constraint experiments were not promoted. An inline
shared finisher preserves ordinary allocation behavior.

Ordinary capture (96 checks) and comparison (65 checks) pass on each version:
32 configurations retain bitwise results and unchanged physical state and warmed
allocation counts. This means zero added allocation, not zero total allocation.
Logs: `...v4_ordinary{112,113}.log` and `...v4_ordinary_comparison{112,113}.log`.
Final singleton-grid boundaries pass 48 checks on 1.12. JuliaFormatter 1.0.62 is
idempotent at margin 100 (`...v4_formatter112.log`). Root reviewed the complete
production/test diff and checked the live bytes against the frozen hashes.
All abbreviated filenames above begin `/tmp/pic005_final`.

## Resource and timing measurements

The final cold-path sample passes 96 assertions per version across 32 fixtures,
with three equal allocation samples per fixture. Allocations are
91,800–621,008 B/call on 1.12 and 91,016–611,728 B/call on 1.13:
`/tmp/pic005_final_v4_cold{112,113}.log`.
The wider `floatmax(T) - smallest_subnormal(T)` density difference passes another
24 assertions over eight configurations, allocating 517,992–6,375,016 B/call.
Observed CPU times are 0.553–7.455 ms and wall times 0.546–19.557 ms on this
shared machine (`/tmp/pic005_final_v4_wide112.log`). These are allocated bytes,
not retained heap or peak RSS.

No plans, caches or new state are retained. A conservative Float64, 3D, 64-cell
bound is at most 468 transforms and 10,476-bit final ratio integers; this is a
work/representation bound, not expected latency. Derivation and executable check:
`/tmp/pic005_diagnostic_bound_review.md` and
`/tmp/pic005_diagnostic_band_bounds.py`.

One eight-pair timing window covers 16 configurations on Julia 1.12.7 with one
Julia and one FFTW thread. Every consumed checksum matches. Candidate/baseline
CPU median ratios span 0.978–1.110, with individual ratios 0.643–1.376; wall median
ratios span 0.905–1.380, with individual ratios 0.346–1.902. Other regression jobs
were active. This variation does not support an isolated or portable latency
claim. Raw pairs and load are retained in `/tmp/pic005_final_benchmark_pairs.json`,
`/tmp/pic005_final_benchmark_summary112.log` and `/tmp/pic005_benchmark_load112.log`.

## Preserved failures and remaining work

Validation cases 16 and 20 return the same eight metric values on the baseline
and candidate, on both Julia versions. Case 16 passes all six metrics in native
environments with HDF5 available. An earlier harness lacked a direct HDF5
installation and failed that metric; its failure remains in
`/tmp/pic005_final_regression112.log` and the dependency check in
`/tmp/pic005_original_environment_hdf5112.log`.

Case 20's charge residual is **9.889422779356765e-7**, exceeding its unchanged
**1e-8** tolerance on both baseline and candidate. Its frequency error remains
0.013636642520216225 against tolerance 0.03. This validation is not green and
requires a separate repair. The subsequent particle replay and exact DC bound
verify deposition-induced total-charge drift; see the
[charge-deposition plan](CHARGE_DEPOSITION_CONSERVATION_PLAN.md).
Logs: `/tmp/pic005_final_v4_validation{112,113}.log` and the corresponding
`validation_comparison` logs.

Severe arbitrary-grid cancellation accuracy is also unresolved. Stored nearly
opposite current components on a 5×5 grid produce residuals about 0.915 and
1.286 on the forced normalized path, against converged independent references
about 4.55e-8 and 6.44e-17 for Float32/64. The original 8-pass/2-fail probe remains
`/tmp/pic005_subnormal_divergence_probe112.log`.
**VERIFIED:** the independent decomposition in
`/tmp/pic005_conditioning_error_budget_probe112.log` (16 checks) measures
normalized current-forward error budgets 3.2875/2.7866, while inverse errors are
3.84e-8/2.15e-16. Significant error precedes exact modal composition. Merely
selecting this same path more often or refining only its inverse cannot recover
those lost coefficients. This finding is confined to the measured algorithm and
inputs; it does not establish a universal accuracy limit.

A stronger accuracy contract needs a separately verified transform or error
policy. The detailed assessment, formulas, source references and test requirements
are retained in `/tmp/pic005_conditioning_accuracy_assessment.md`. Source-scale
heuristics or agreement between two FFT band widths must not be labeled a
certified error bound. Keep the existing severe-case target and validation
thresholds unchanged while investigating that follow-up.


## Legacy diagnostic plan reuse — pre-edit brief

- **Mode:** measured performance repair in the existing EMPIC1D diagnostic.
- **Acceptance:** reuse owned FFT plans; preserve transform order, normalization,
  DC/Nyquist rules, valid grid-shaped state, errors and physical state. Require
  exact differential results and lower measured allocation; no accuracy claim.
- **Baseline:** AbstractFFTs `definitions.jl:61–68` constructs a new plan for each
  raw `fft!`/`ifft!`. The frozen two-line prototype passes 12,246 checks across
  Float32/64, 26 lengths, one/two FFTW threads, repeated buffers, public range and
  error cases, and changes to global thread settings after grid construction.
  Log: `/tmp/pic005_plan_reuse_probe112.log`.
- **Verified cause:** the legacy method plans per call despite grid-owned
  ESTIMATE plans in SpectralOperators `periodic_fourier.jl:103–104`. General EM
  and canonical spectral operators already use those plans.
- **Contract:** the same separate temporary density buffer survives canonical
  derivative scratch reuse. Existing grid thread-safety constraints apply;
  no cache, state, schema, dependency or public API is added.
- **Owner:** only the forward and inverse calls in
  `charge_conservation_residual(::EMPIC1D, dt)`; use `es.g.plan`/`es.g.iplan`.
- **Compatibility:** existing grids retain the plans selected at construction,
  matching the other spectral paths. A later global FFTW setting change does
  not re-plan this one diagnostic. No dynamic re-planning contract was documented.
- **Strategy:** replace the two raw calls with the existing plans. The measured
  source is frozen at `/tmp/pic005-plan-reuse-t516vhs3`, SHA256
  `699decbfa1eef1ff70cbb081ce6ddc656cc68199f067fbfbc8e422652863f9af`.
- **Verification:** latest-version differential probe; every complete testset
  directly calling the diagnostic in both EM test files; both-version concrete
  and package compiler gates, formatter and final diff. Existing full range
  assertions remain intact. Wrong normalization, overwritten scratch and wrong
  axes/null masking must fail independent DFT and nonzero-current controls.

One completed eight-pair benchmark over 24 configurations measures exactly
544 fewer allocated bytes/call. Candidate/baseline CPU median ratios are
0.025–0.033 at N=4/5, 0.087–0.123 at N=128, 0.154–0.186 at N=129,
0.439–0.582 at N=1024 and 0.734–0.884 at N=4096. All consumed checksums match;
both benchmark servers exit successfully. Large-grid individual ratios reach
2.159 under shared load, so this is a local measurement, not a portable speed
claim. Raw evidence: `/tmp/pic005_plan_reuse_benchmark112.log` and
`/tmp/pic005_plan_reuse_benchmark_pairs.json`.

Final gates pass on Julia 1.12.7 and 1.13.1: 19,980 runtime assertions per
version (12,246 differential/applicability, 4,456 complete EM caller testsets,
3,270 restart assertions and eight validation-16 checks). All six validation
metrics exactly match the committed baseline. All 33 concrete compiler calls
are clean; package findings retain the exact baseline leaf-message multisets
(26/25). JuliaFormatter 1.0.62 is idempotent and the final diff contains only
the two planned FFT-call replacements. Logs and executable commands are in
`/tmp/pic005_plan_reuse_handoff.md`; the runtime logs use the prefixes
`/tmp/pic005_plan_reuse_{probe,regression,integration}{112,113}.log`.
Validation 20 and severe cancellation accuracy remain separate open findings.


## Current-transform follow-up on Julia 1.13.2

The actual pushed `6ad7276` source still reproduces the separate current
materialization defects on Julia 1.13.2. The unchanged exact integer-root probe
passes **798 defect-witness assertions** plus eight dyadic representation
controls, terminal exit 0. Its finite expected currents become nonfinite after
raw forward FFT or inverse materialization; raw-current and charge arrays remain
unchanged. This verifies the remaining defect, not a repair. Source SHA256:
`07ffcf03d21f7572e9bf5068321dd7ab1c5df414654667d121a6a3e245130e49`.
Driver: `/tmp/pic005_raw_current_coupled_probe.jl`; refreshed log:
`/tmp/hybrid-current-materialization-6ad7276-baseline1132.log`.

An isolated rational-inverse experiment factors a common odd denominator from
all exact spectral coefficients, leaving dyadic scaled inputs for the existing
band transform. It restores that denominator after the transform and before
final hardware conversion. Unlike a fixed binary tail cutoff, this retains the
entire represented rational input. The least common multiple can grow with grid
size and denominator diversity; resource measurements remain required.

The first Julia 1.13.2 probe passes **508 assertions**, terminal exit 0, against
independent exact forward/inverse DFTs on 1D–3D grids with axes of length 1, 2 or
4. It covers signed complex inputs, widely separated exponents, different odd
denominators, dyadic equality, invalid rational inputs, and the retained
midpoint/overflow-tail counterexamples. This proves these controls only; FFTW
twiddle rounding on arbitrary grids and severe cancellation are still open.
No production code is edited. Brief and source:
`/tmp/hybrid-current-rational-transform-{brief.md,prototype.jl}`; log:
`/tmp/hybrid-current-rational-transform-prototype1132.log`.


Both refreshed Julia 1.12.7 runs also finish successfully: 806 baseline
witness/representation checks and 508 isolated rational-transform checks.
Logs use the same prefixes with suffix `112.log`. The actual-grid denominator
probe passes 36 divisibility/parity checks on 1.13.2. For anisotropic Float64
16³ geometry, 448 distinct odd denominator factors produce a 44,918-bit common
denominator although each individual factor has at most 112 bits. That width
alone spans 1,152 bands of the existing 39-bit band width. This is a measured
representation-growth warning, not an observed count of executed FFTs or full
call memory/time. Complete current performance must be measured before accepting
this representation. Evidence: `/tmp/hybrid-current-rational-denominators.jl`
and `...denominators1132.log`.


The subsequent complete-current prototype passes **884 assertions per version**
on 1.12.7 and 1.13.2 (876 raw-current/coupled-axis checks plus eight representation
controls), both terminal exit 0. It retains exact density/current coefficients
through projection and the rational inverse; all analytically finite expected
currents in the refreshed baseline witnesses become finite and match their
independent integer-root wave values. Raw current and charge snapshots remain
unchanged. This is an isolated general-EMPIC experiment, with no ordinary fast-path,
legacy-family, public-step or complete-call resource acceptance yet.
Source and logs: `/tmp/hybrid-current-rational-pipeline-prototype-v2.jl` and
`...prototype-v2-{112,1132}.log`. The original driver omitted FFTW from its imports
and failed with three passes/one error before evaluating the prototype current;
that original source/log is retained. The corrected driver adds only the missing
import, without altering numerical assertions or production code.


An independent density-source extension passes **9,912 assertions per version**
with terminal exit 0. It uses analytic cosine/sine modes and independent
nearest-even rounding cells to cover subtraction/forward/inverse range loss,
subnormal source and timestep factors, anisotropic derivative scales, and the
complete `floatmax - tiny` difference on every axis in 1D–3D. Source, raw-current
and physical field arrays remain unchanged. Driver:
`/tmp/hybrid-current-rational-density-prototype.jl`; logs:
`...density-prototype{112,1132}.log`. These are still isolated pipeline checks;
ordinary dispatch, arbitrary-grid conditioning, public steps, legacy current
and production resource acceptance remain open.

The first complete-call resource probe passes **736 exact spatial/repeat
assertions** on Julia 1.13.2, terminal exit 0. It uses origin impulses on eight
Float32/64 anisotropic grids with 4, 16, 32 or 64 cells and an independent
integer-root inverse DFT. All expected finite cells satisfy their nearest-even
rounding intervals; the baseline produces 1, 7 and 18 nonfinite component cells
on each type's 16-, 32- and 64-cell cases. The prototype allocates 177,016 to
8,694,008 bytes per complete call. Shared-host median times range 0.244 to
14.157 milliseconds; these are measured small-grid observations, not bounds or
ordinary-path speed comparisons. Driver: `/tmp/hybrid-current-rational-resource-prototype.jl`; log:
`...resource-prototype1132.log`. No production current path has been changed.

The larger-grid extension passes 40 finite-output, repeat and state-preservation
checks on 1.13.2, terminal exit 0. Complete 8³ calls allocate **810,838,152 B**
for Float32 and **811,238,592 B** for Float64; three warmed samples take
3.35–5.17 s and 2.28–2.61 s, respectively, under shared host load. These are
aggregate allocated bytes, not peak memory. The 8² controls allocate about
11.8 MB per call. This measures a substantial exceptional-path cost; it does
not establish exact spatial accuracy for those non-integer twiddles or a
portable speed bound. Evidence: `/tmp/hybrid-current-rational-resource-larger.jl`
and `...resource-larger1132.log`. A smaller representation or measured justified
cost policy is still required before production acceptance.
