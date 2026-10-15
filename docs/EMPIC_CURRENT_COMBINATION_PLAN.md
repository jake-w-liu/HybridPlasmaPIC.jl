# Complete cross-species EMPIC current

Isolated baseline/candidate: `/tmp/hybrid-empic-current-combination-{base,v1}`, archived from `240f77c814d45cd966ac30ccbd38c5b96796af1d`. Charge combination is tracked separately in [MOMENT-002](EMPIC_CHARGE_COMBINATION_PLAN.md). Candidate prerequisite: `/tmp/hybrid-empic-charge-status-prerequisite-v1/manifest.json`, verified before copying only `src/coupling/moments.jl` and `src/parallel/threads.jl`. Their hashes are retained as `prototype-status-prerequisite.json`. Current-specific production edits belong only to the current sections of `src/models/full_pic.jl`.

## Pre-edit brief

- Mode: repair the current combination owner, with scoped cleanup of borrowed midpoint positions.
- Acceptance: for supported finite represented inputs with per-species or combined range loss, round the complete electron+ion current once. Preserve ordinary arithmetic, particle state, final scratch and exception visibility.
- Baseline: `/tmp/hybrid-empic-current-combination-baseline.jl` passes 42 range witnesses plus eight restoration-defect witnesses on both native versions, logs `...baseline{112,113}.log`, both exit 0. The actual midpoint owners return NaN for exact-zero overflow cancellation, zero for an exact 3/4 subnormal quantum, and one quantum for an exact half-quantum tie. Legacy DC mean has the same subnormal/tie defects. A deterministic ion-stencil exception leaves ions at midpoint coordinates in both families while electrons have already been restored.
- Cause: VERIFIED in base `full_pic.jl:307-310`, `:1049-1080`, `:1617-1630`: separately rounded species current contributions are added after their original information is lost. Temporary position borrowing lacks `finally`, so the witnessed stencil exception bypasses restoration.
- Invariants: public electron/ion charges are -1/+1; velocities are at the half level, and each species deposits at its own wrapped midpoint. An immobile background contributes no current. Ordinary order is electron then ion. Restore original positions and midpoint buffers on success and scoped deposition/replay exceptions; propagate the original exception without converting it into a result. No whole-step rollback promise.
- Canonical owners: `_deposit_current_at_midpoints!` (general three-component Jraw), `_deposit_Jy_at_midpoints!` (legacy transverse Jy), the mean combined before `_esirkepov_Jx!`, and their current helpers. Reuse the shared `_finish_moment!` status and the existing exact deposition/conversion owners. Do not edit charge sections.
- Compatibility: public interfaces and normal arrays/scratch remain unchanged. Preserve `_species_Jx0` scalar result through a delegate to its status-owning implementation; private deposition helpers carry the required status to their only production callers. Nonfinite/symbolic ineligible inputs retain legacy propagation. Finite exceptional totals and previously leaked temporary positions intentionally change.
- Strategy: retain per-component eligible/lost flags, check the native combination, and cold-replay both original species before one conversion. Revisit each species' own midpoint under scoped restoration; use no new particle copies or public model fields. Preserve normal electron-first processing and restore it before touching ions. Scalar mean uses the same pattern without position dependence. The cold replay is bounded per component, with operand-dependent exact storage.
- Verification: independent exact rational cells and tie controls, nonzero finite cancellation, distinct midpoint locations, views/unequal/empty populations, all native types/shapes/dimensions, nonfinite/symbolic/retained BigFloat coefficients, scratch and ordinary baseline equality, injected deposition/replay failures, selected existing continuity/equivalence/subcycling/error-order tests, both native JET/Aqua with full baseline diagnostics, formatter, exact diff, ordinary allocations/time and cold resource measurements.

## Ownership and ordering

General flow: drift fills `es.mide`/`es.midi`; `_deposit_current_at_midpoints!` temporarily swaps electron coordinates, fills `worke`/`ne` and `Jraw`, restores electrons, then repeats for ions with `worki`/`ne` and accumulation. `worke/worki` end at w*v[3]; `ne` ends as the last raw species/component deposit. Midpoint arrays must return to their entry values. `_correct_current_continuity!` then projects the longitudinal Fourier component and removes pure Nyquist modes; its independent FFT/projection range is unchanged.

Legacy flow: drift stores midpoints in `es.work`/`es.worki`; deposition parks x^(n+1) in `Bzp`/`Bzi`, uses `Eyp`/`Eyi` for w*v_y and `dBz` for raw grid scratch, then restores positions. These final scratch meanings remain. `_zero_nyquist_mode!` runs after complete species combination; its independent reduction/projection range is outside this slice. The longitudinal DC mean combines q*sum(w*v_x)/L in electron-first order and then feeds the unchanged continuity solve. Charge snapshots and the nonzero Fourier modes remain root/other-owner arithmetic.

Public `step_empic!` first validates electron/ion identity, timestep, n_sub/clock and relativistic velocities, then sizes buffers and advances. The slice does not move those checks. Temporary current positions are restored even if a shape rejects deposition; already-completed particle pushes or grid writes are not rolled back. General midpoint length validation will precede its first coordinate swap, preserving validation order while preventing partial-coordinate borrowing.

## Performance and scope

Ordinary native outputs and all scratch/positions must match 240f77c bitwise; no rational grid or new particle-sized buffer belongs to the ordinary path. Measure allocations and warm latency against the actual archived owner calls, not rewritten reference loops. Cold replay may revisit each species and allocate rational grids; measure complete three-component cost and retained storage without inventing a fixed ceiling. The previous 22/21 package JET terminal messages are the diagnostic baseline, not merely the unchanged ceiling of 30.

Open boundaries: generic raw deposition; charge/background combination (MOMENT-002); pressure; hybrid ratios; MPI reductions; ordinary in-range cancellation; stateful/unsupported custom-shape semantics; independent Fourier/Nyquist and field-update ranges. Wrong patches rejected by tests include final-Inf-only detection, adding already-rounded exact species, depositing both species at electron midpoints, retaining parked positions on error, or changing all ordinary calls to exact arithmetic.

## Native verification in progress

Each native Julia version passes 892 permanent assertions, 70 concrete JET
analyses with zero reports, three package JET assertions and nine Aqua checks.
The complete terminal diagnostic multisets match the actual baseline: 22 reports
on 1.12.7 and 21 on 1.13.1. Logs:
`/tmp/hybrid-empic-current-combination-quality{112,113}.log` and
`...diagnostics{112,113}.log`, all exit 0. Selected physical integrations and
public finite-charge/overflowing-current step controls remain required before
promotion. The current-only source remains isolated; charge/status prerequisites
are tracked separately.
