# Resistive diagnostic range audit

## Pre-edit brief (2026-10-27)

- Mode: repair, with an independent transport/collision integration audit.
- Requirement: the public `resistive_dissipation(J, eta, g)` must evaluate its documented Ohmic quadratic integral for finite Float32/Float64 inputs without losing a representable answer through intermediate range loss. Retain ordinary operation order, validation and nonfinite-input behavior.
- Baseline: clean `cd558315956df0debb14b4a2fb7790eaef2821e6`; Julia 1.13.3, one cell of unit length, `J = (2^600, 0, 0)`. Actual eta=0 result is NaN (exact zero); eta=2^-600 result is Inf (exact finite 2^600). The raw `resistive-before-1133.log` is retained in the dated takeover evidence directory.
- Root cause: `src/diagnostics/conservation.jl` squares and sums the current before multiplying the coefficient and cell volume. Intermediate overflow or underflow may disappear in the complete expression.
- Contract: a read-only scalar diagnostic; no simulation state or RNG changes. The validated coefficient at grid precision and every represented cell spacing define the integral. Genuine final overflow remains Inf. Nonfinite currents retain the existing IEEE result; this repair does not silently sanitize them.
- Canonical owner: `resistive_dissipation` in `conservation.jl`, sharing `_spectral_product_lost_range`, `_field_energy_upper_exponent`, `_field_square_sum_exact` and `_from_exact` with existing numerical owners. No alternate diagnostic API or duplicated exact square reduction.
- Compatibility: retain the ordinary sum and product sequence bitwise whenever range checks certify it. Preserve dimensions, coefficient validation, return precision and read-only ownership.
- Strategy: track lost nonzero squares while performing the existing reduction; check coefficient/volume products and every contiguous spacing subproduct. On finite binary32/binary64 range loss, evaluate the complete integral through the existing exact sum and final converter. Keep nonfinite and other-precision arithmetic unchanged.
- Verification: original reproductions on minimum/latest Julia, exact independent rational reference for current/coefficient/volume boundaries, normal-state bitwise replay, zero/nonfinite behavior, malformed dimensions/coefficient rejection, unchanged input snapshots, warm allocation and concrete JET checks, relevant budget and transport/collision suites, formatting and full diff review.

## Audit matrix

| Surface | Domain | State / robustness | Resources | Status |
| --- | --- | --- | --- | --- |
| Transport stage Ohm/pressure and ion force | Paired stage current, eta=0 and common-time energy | Priming, Boris and predictor share Eforce; density/pressure rollback | Six force grids and warm step checks | Current min/latest suites pass |
| Collisional hybrid composition | Canonical collision placement and source interval | Full field/particle/RNG rollback, ledger, fixed config | Shared workspaces and bounded output ownership | Current min/latest suites pass |
| Raw/archive restart | Force policy and stochastic continuation | Incoming graph, runtime, phase and destination configuration | Detach incoming RNG metadata before commit | Current min/latest suites pass |
| Ohmic diagnostic | Complete nonnegative integral | Zero eta, finite range extremes, invalid input | Preserve normal allocation/operation order | Canonical repair and all local gates pass |
| Cross-feature closure boundary | No unproved combined collision/transport physics | Collisional wrapper explicitly rejects non-isothermal closure | Not applicable | Current bounded limitation, not a defect |

## Acceptance tracker

- [x] Rebuild canonical call and persistence paths from current sources.
- [x] Reproduce the diagnostic defect on latest Julia 1.13.3.
- [x] Implement the canonical repair with regression coverage.
- [x] Pass minimum/latest correctness, robustness and resource checks.
- [x] Review compiler reports and exact final diff.
- [ ] Root integration, main push and exact-SHA CI acceptance.

## Wrong-repair controls

- Returning early only for eta=0 misses the nonzero finite exchange after overflowing/underflowing current squares; the exact range matrix rejects it.
- Scaling current squares but using a rounded cell volume loses restored D=2/3 products; all spacing orderings are checked against the complete exact integral.
- Multiplying eta into each term or changing summation order can alter ordinary runs; array and view profiles compare every result bitwise with the original loop and require zero warm allocation.

## Verified local results

On Julia 1.12.7 and latest installed 1.13.3, using registry-resolved
SpectralOperators 0.1.3 (`036881a8627ae86e285825891bb427dd66508c08`):

- The unchanged parent fails all four original finite-current assertions per runtime: zero eta produces NaN and the small nonzero coefficient produces Inf. The repaired diagnostic passes the same assertions.
- The new regression module passes 4,218 assertions per runtime, including 108 ordinary array/view profiles with bitwise compatibility and zero warm allocation. Post-format runs repeat these results; budget, phase-6 diagnostic, transport and Aqua neighbors also pass (6,625 assertions per runtime including the new module).
- The fresh current-parent transport, transport I/O, collisional hybrid/checkpoint/physics, multispecies, evolving-energy and Aqua suites pass 7,292 assertions per runtime. These include bitwise eta=0 compatibility, common-time energy convergence, exact stochastic restart and malformed-state rejection.
- An independent dense-DFT audit of even/odd D1–3 grids passes 80 assertions per runtime: current, Faraday difference, pressure heating and integrated Ohmic exchange agree, and the unsupported collision/transport closure combination rejects before field/clock/RNG mutation.
- Current-owner failure/interrupt injection passes 544 assertions per runtime across eight boundaries and both RNGs. Independent serialized-snapshot replay passes 682 assertions per runtime, including particle, field and RNG mutation controls.
- Complete compiler report trees equal the parent after replacing only the exact worktree root (29 reports on minimum Julia, 27 on latest). Every frame and source line remains present; all 12 concrete array/view calls per runtime have zero reports. A changed or missing report is rejected by the comparison.
- JuliaFormatter AST/idempotence checks and `git diff --check` pass. No simulation state, checkpoint schema, collision operator, electron closure or RNG behavior is changed by the repair.

Evidence is retained in
`HybridPlasmaPIC-audit-evidence/2026-10-27/transport-collisions-takeover/`,
including raw compiler trees, exact command receipts, source/environment hashes,
original failures and corrected harness revisions. The initial formatter
comparison mistakenly compared line-number nodes; its failure is retained and
the normalized AST check and post-format tests pass. The independent DFT probe's
initial fixture failed only its fixed positive-signal gate; increasing the input
field amplitude preserved every assertion/tolerance. The collision harness was
updated to the actual current field/pair owners; its initial path and clock-
expectation failures are retained with source-backed corrections.

No latency acceptance claim is made under shared-machine contention. The
existing collision wrapper still supports only its documented ideal isothermal
closure; jointly coupled transport/collisions, electron transport S2–S8, MPI and
device sources remain separate implementation work. Root integration and all six
[exact-SHA CI jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40369103017) pass for `a356c6b`; complete compiler acceptance is recorded in the [audit tracker](DEEP_DEBUG_TRACKER.md#current-progress).
