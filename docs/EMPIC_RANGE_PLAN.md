# Electromagnetic field-update range repairs

Status: finite-curl repair pushed as `0f53361`. Scaled spectral updates are
verified on Julia 1.12.7 and 1.13.1 and pushed as `8e7d5be`. The finite-Fourier
longitudinal projection repair is verified on both versions. Complete continuity
pipeline range repairs remain open below.

## Verified gaps

- Squaring an accepted finite light speed can overflow. A zero-field vacuum
  then acquires nonfinite fields on both EM interfaces and Float32/64.
  Evidence: `/tmp/empic_light_speed_range_probe112.log`.
- Even when the squared speed is representable, forming the magnetic force or
  timestep increment first can overflow or underflow before the complete
  electric-field result is formed. A subnormal longitudinal increment can
  also lose the final rounding tie. Evidence:
  `/tmp/empic_maxwell_neighbor_probe112.log` and
  `/tmp/empic_longitudinal_range_probe112.log`.
- A spectral derivative can itself overflow before multiplication by a small
  timestep. Both EM interfaces produce infinity for a case with a finite
  expected field. Evidence: `/tmp/empic_spectral_range_probe112.log`.
  The finite-curl repair does not recover an already nonfinite derivative.

## Finite-curl repair

Keep the existing Maxwell operation order for ordinary inputs. Detect range
loss in the squared speed, magnetic product, force, increment and final result;
then evaluate the complete electric-field expression through the existing exact
arithmetic owner. Construct exact timestep coefficients lazily once per substep.
The specialized 1D and general 1D–3D solvers share this calculation, including the
specialized longitudinal field. Existing nonfinite-input rejection remains.

Independent high-precision expressions and discrete wave evolution check both
signs, all axes, extreme and neighboring speeds, changing timesteps/subcycles,
vacuum, cancellation and genuinely unrepresentable final results. Reordering
`dt*c*c` alone fails a separate underflow case; replacing nonfinite values with
zero would fail the nonzero-wave and input-rejection checks.

Frozen full EM regressions pass **9,501 assertions on each of Julia 1.12.7 and
1.13.1**. Logs: `/tmp/empic_maxwell_final_regression112.log` and
`/tmp/empic_maxwell_final_regression113.log`. The earlier fixture expected NaN
propagation despite the spectral input guard; the corrected test requires its
existing rejection, checked separately against the baseline (13 assertions).

All 64 ordinary configurations preserve captured state and allocations exactly
(129 comparison assertions; `/tmp/empic_maxwell_compare112.log`). Ordinary
steps allocate 1,280–7,680 B, unchanged. Sixteen measured extreme wave cases
allocate 10,560–403,360 B per step. Paired local timing has substantial variation;
no ordinary-path speedup is claimed. Logs:
`/tmp/empic_maxwell_cold_resources112.log` and
`/tmp/empic_maxwell_benchmark_pairs112.log`.

Final Julia 1.12 compiler checks pass 26 assertions: all 16 concrete helper and
ParticleSet-ion entrypoints have zero reports, with 10 runtime assertions.
Package JET retains 24 reports against the unchanged ceiling of 30. Eight
optional-Nothing entrypoints retain the same two reports as the baseline;
the exact report comparison passes 17 assertions. Logs:
`/tmp/empic_maxwell_jet_patched112_retry.log`,
`/tmp/empic_maxwell_jet_baseline112.log` and
`/tmp/empic_maxwell_jet_comparison112.log`. Two intervening import checks failed
after the snapshot manifest disappeared; restoring its identical baseline copy
allowed the retry to pass. The cause of that file removal is unverified.

## Scaled spectral repair

Trace derivative consumption in magnetic centering, Faraday advancement,
electric advancement and initial magnetic priming. Define a shared scaled
spectral operation that retains a representable complete update without first
materializing an unrepresentable curl. Preserve the dependency's normalization,
Nyquist rules, input rejection and ordinary results. Validate against an
independent discrete Fourier sum, then test both solvers, axes, precisions,
phase transitions, allocation and compiler behavior before committing.
The independent 2,048-bit discrete Fourier probe now reproduces both electric
advancement and initial magnetic-field failures in 48 assertions across
Float32/64 and both solver families in 1D–3D. Log:
`/tmp/empic_scaled_curl_ownership_probe112_retry.log`.

The canonical spectral call runs first, retaining input rejection and ordinary
FFT arithmetic. Range-sensitive inputs are differentiated in disjoint, exactly
scaled binary bands. Signed scales survive through the complete field update.
Selection includes longitudinal and diagonal components because the dependency
uses a common curl normalization. Magnetic priming, centering and advancement
share one owner that retains the original timesteps before halving them.
The state layout, public API and checkpoint schema remain unchanged.

These remain FFTW-precision derivatives. The repair does not recover modes lost
through ordinary FFT cancellation or bypass the spectral operator's dynamic
range limits. Exceptional scratch derivatives can be unrepresentable; physical
updates consume the scaled parts. Each derivative retains at most
`cld(hi-lo+1, width+1)` arrays, bounded by the hardware exponent range; empty
bands allocate no output array. No persistent cache is added.

Both Julia versions pass **19,765 assertions** across the complete EM families.
Independent real-space DFT checks cover finite complete updates, tiny half
intervals, signed curl cancellation, scale bands, phase transitions, changing
steps/subcycles and preserved rejection. Logs:
`/tmp/empic_scaled_final_v3_regression112.log`,
`/tmp/empic_scaled_final_v3_regression113.log`.
Source SHA256: `3e83ee7ddde7e448436a539ca427e17fc433d21b790b46e0a43ad7635d8861b7`.

All 64 ordinary configurations retain bitwise state/particle results and
1,280–7,680 B per step (129 comparison checks). Julia 1.12 JET reports zero
findings for all 54 concrete signatures; 62 assertions pass. Whole-package
analysis retains the matched `aba03a9` baseline's 26 reports against the
unchanged ceiling of 30. Optional-Nothing entrypoint reports also match exactly
(17 checks). Logs: `/tmp/empic_scaled_compare_v3_112.log`,
`/tmp/empic_scaled_jet_v3_112.log`,
`/tmp/empic_scaled_compare_compiler_v3_112.log`.

Exceptional public stepping passes 128 resource assertions across 32 cases on
four-cell-per-axis meshes. Measured allocations are 36,040–1,359,512 B per step;
case median times are 34–1,597 μs. These are fixture costs, not upper bounds.
Eight alternating paired ordinary batches give CPU time ratios of 0.887–1.182
and wall-time ratios of 0.896–1.142, with broad variation on the shared machine.
Specialized 1D vacuum CPU medians increase 18.2% for Float32 and 11.6% for
Float64. No uniform speedup or portable timing claim follows. Logs:
`/tmp/empic_scaled_final_cold_resources112.log`,
`/tmp/empic_scaled_final_benchmark_pairs112.log`.

Final formatter idempotence and diff checks pass. Superseded candidates,
corrected negative-control fixtures, complete resource protocol and raw timing
ranges are retained in `/tmp/empic_scaled_v3_evidence.md`. Integrated full
package runs remain a separate gate after the outstanding slices are ready.

## Finite-Fourier longitudinal projection repair

VERIFIED: general electrostatic and electromagnetic initialization produce zero
electric fields for a four-point nonuniform charge when a finite nonzero
wavenumber's square overflows or underflows. The specialized 1D solvers match
the finite inverse-DFT field. The public Float32/64 probe passes 28 witness
checks on each Julia version: `/tmp/pic_spectral_k2_probe112.log` and
`/tmp/pic_spectral_k2_probe113.log`. Further probes reproduce lost products even
with finite normal squared norms: `/tmp/pic_projection_complete_expression_probe112.log`.

Shared scalar kernels now classify null modes from the actual wavevector
components. Ordinary operations retain their order; range-sensitive modes
evaluate the complete longitudinal expression exactly from the original field,
source and divisor, rounding only the final components. Poisson, Gauss correction
and current correction share that owner. Current passes the original timestep
separately, preserving results whose intermediate charge-change/timestep ratio
is unrepresentable. True DC, pure Nyquist and mixed modes keep their distinct
contracts. The change begins with finite represented Fourier coefficients;
the seven pipeline boundaries below remain separate work.

The final immutable candidate passes **27,434 assertions on each Julia version**
across the complete ES, EM and specialized 1D suites. Independent exact real/imaginary
oracles check final rounding cells, signed underflow, overflow thresholds, ties,
all dimensions and axes, mixed scales, and null-mode behavior. Logs:
`/tmp/pic_projection_final_v2_full112.log` and `...full113.log`.
All 108 concrete compiler signatures per version have zero reports. Package
report contents match the baseline exactly: 26 on 1.12 and 25 on 1.13, with the
ceiling unchanged. Evidence: `/tmp/pic_projection_final_v2_package_comparison.log`.
An earlier tuple-cardinality inference regression was repaired before these gates.

Eight ordinary ES and 64 EM configurations preserve captured physical state and
particles bitwise, with unchanged step allocations (146 comparisons). Ordinary
helpers allocate zero bytes. Extreme public-operator samples allocate
18,464–3,941,136 B per call; these are fixture costs, not upper bounds. Logs:
`/tmp/pic_projection_final_v2_state_compare112.log` and `...cold_resources112.log`.

Eight balanced alternating timing pairs on the shared machine give particle-bearing
CPU median ratios of 0.901–1.032 and vacuum ratios of 0.930–1.223. Float32 2D/3D
vacuum medians rise 22.3%/18.9%; Float64 2D/3D rise 10.3%/6.5%. Paired ranges are
broad, including variation in the unchanged specialized 1D control. These measured
costs remain optimization targets; no speedup or isolated causal timing claim is
made. Raw evidence: `/tmp/pic_projection_final_v2_benchmark112.log`.
Formatter idempotence and final diff checks pass. Full-package integration remains
a separate gate after the outstanding I/O slice is ready.

## Continuity pipeline follow-up

VERIFIED: an independent four-point Fourier probe reproduces seven additional
boundaries on Julia 1.12.7 (256 witness assertions, terminal success): physical
density subtraction, forward FFT, inverse current FFT, diagnostic density
inverse FFT, timestep/divergence cancellation, and standalone divergence
underflow and overflow. Float32/64 and 1D/2D cases separate these stages;
expected complete currents or normalized residuals remain representable.
Log: `/tmp/empic_boundary_audit_v2_112.log`. The frozen source SHA256 is
`3e83ee7ddde7e448436a539ca427e17fc433d21b790b46e0a43ad7635d8861b7`.

For example, with M=floatmax(T), four-cell density delta=.1M*cos(x), dt=.25,
and k=1, corrected Fourier current coefficients approximately .8M remain
finite. The inverse transform nevertheless stores infinite current instead of
the representable -.4M*sin(x). Separately, a diagnostic can lose a tiny
divergence before multiplication by a large timestep, reporting one instead
of a residual near roundoff. These are failures of the complete pipeline;
standalone unrepresentable derivatives are not themselves dependency defects.

Plan a shared scaled density-difference/FFT path that preserves signed scales
through projection and final current materialization. Compose the diagnostic's
timestep-weighted divergence and density parts before forming its numerator
and normalization. Preserve ordinary arithmetic, DC/Nyquist behavior and input
rejection. Verify independent Fourier sums, all dimensions, precisions,
public step chronology, resources and both Julia versions. The scalar-mode
repair alone cannot close these boundaries. Full derivations, rejected
hypotheses and exact test limits are retained in
`/tmp/empic_boundary_audit_handoff.md`.
