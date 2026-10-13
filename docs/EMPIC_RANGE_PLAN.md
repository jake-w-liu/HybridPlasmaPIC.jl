# Electromagnetic field-update range repairs

Status: finite-curl electric-field repair verified on Julia 1.12.7 and 1.13.1.
Scaled spectral updates remain open.

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

## Next repair

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
