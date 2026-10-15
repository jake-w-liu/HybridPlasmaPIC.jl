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
   background. Both are **untested repair hypotheses**, not verified fixes.
   Prototype them against the captured states and complete validation case.
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

These are fixed-state results, not an implemented or end-to-end repair. Full
wave evolution, other dimensions/precisions/shapes, mobile and signed populations,
nonfinite behavior, spatial/adjoint guarantees and measured cost remain open.
Compensated accumulation passes this one threshold narrowly; no general
sufficiency claim follows. Choose the production owner and rounding contract
only after those comparisons.
