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

## Verification and promotion

Status: pushed as `6ad7276`; [GitHub run 38754771205](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38754771205)
is active. Final combined source passes both native gates; GitHub acceptance
remains pending. The merged checkout
`/tmp/hybrid-empic-current-combination-merged-v1` starts from `4cd642c` and
retains the complete charge repair. Exactly three production/test files change.

Both native Julia versions pass **4,805 checks** on the final merged source:
1,953 normalization + 1,720 charge + 892 current + 128 public current-step +
30 independent Nyquist/range proofs + 70 clean concrete JET + three package
JET + nine Aqua assertions. All 22/21 package terminal diagnostics match the
actual baseline; the ceiling remains unchanged. Logs and parsed totals:
`/tmp/hybrid-empic-current-combination-merged-quality{112,113}.log`,
`...merged-summary{112,113}.json`, `...merged-diagnostics{112,113}.log`, all exit 0.

Each version separately passes 1,862 affected existing integration assertions,
including legacy/general continuity, subcycling, equivalence, validation and
energy tests. The selection driver checks complete inclusion of its named sets;
it is selected coverage, not a full package suite. Logs:
`/tmp/hybrid-empic-current-combination-integration-quality{112,113}.log`.
The public-step control uses finite charge near floatmax/2, matched moving species
and individually overflowing currents; fields/current stay zero and particle
drift, weights, velocities and clocks retain their expected values for n_sub=1/2.

Ten ordinary archived-baseline comparisons preserve all output, particle and
seeded scratch arrays bitwise and allocate **0 B** on each version. Observed
minimum-time candidate/baseline ratios are 0.806–1.270 on 1.12 and 0.810–1.095
on 1.13; these shared-host observations do not establish a portable speed claim.
Six complete cold three-component cases allocate 790,680–61,238,136 B total,
not peak RSS. Exact replay cost depends on cells, particles and operand bit lengths.
Evidence: `...resources112.log`, `...supplement113.log` and candidate
`prototype-resource{112,113}.json`.

Final read-only formatting/idempotence passes for seven affected/prerequisite
files. All 179 merged production/test hashes, both prerequisite hashes and
exact current-region equality pass the final audit. The promoted three files
match the tested merged snapshot byte for byte. Merged `full_pic.jl` SHA256:
`07ffcf03d21f7572e9bf5068321dd7ab1c5df414654667d121a6a3e245130e49`.
Adjusted patch SHA256:
`ca2a0c849b09e9f156ff8646a151b35666d0f96b0831cc0615d5ba913922c064`.
Evidence: `...format-merged.log`, `...final-audit.log` and
`/tmp/hybrid-empic-current-combination-artifact-manifest.json`.

Initial failing fixture controls remain retained: an ordinary BigFloat component
had no range event; inactive ion scratch was uninitialized in a comparison;
and a public-step test referenced nonexistent `rho` instead of `rho_n`.
Independent controls corrected these fixtures without source or tolerance changes.
The full handoff records the failures and exact reproduction commands:
`/tmp/hybrid-empic-current-combination-handoff.md`.

- [x] Reproduce range and position-restoration defects on both native versions.
- [x] Pass numerical, error, compiler, quality, integration and resource gates.
- [x] Verify combined charge/current source on both native versions.
- [x] Review final diff/format and promote the exact tested patch.
- [x] Commit and push to main (`6ad7276`) with its own GitHub acceptance run.
- [ ] Observe every GitHub job green for that revision.
