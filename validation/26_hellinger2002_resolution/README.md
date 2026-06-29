# Hellinger et al. 2002 — Perpendicular-Shock Resolution Dependence

Validates against **Hellinger, Trávníček & Matsumoto (2002)**, *Reformation of
perpendicular shocks: Hybrid simulations*, Geophys. Res. Lett. **29**(24), 2234,
doi:`10.1029/2002GL015915` (paper in `ref/`).

## Reference result (Fig 1)

The maximum magnetic gradient `max|dB_y/dx|/B₀` in the shock front vs grid spacing
`dx`, for upstream β_p = 0.2, 0.5, 1.0 (β_e = 0.5, M_A ≈ 6.6):

- **β_p = 1.0 (hot):** gradient is **resolution-independent** (≈ const ~5) — proton
  reflection stops the steepening; the shock is quasi-stationary.
- **β_p = 0.2 (cold):** gradient **rises** as `dx` decreases (≈ 5 → 13 → 25 for
  dx = 0.16 → 0.08 → 0.04) — grid-determined, nonstationary. Hellinger's headline is
  that 1-D hybrid codes **cannot** describe the nonstationary (cold) case.

## What is tested

The physically meaningful **trend** (not the exact magnitudes): with sustained
injection, the model's front max-gradient is approximately resolution-independent for
β_p = 1.0 and rises with refinement for β_p = 0.2, and the cold case steepens more
strongly than the hot case. Only the resolvable `dx = 0.16, 0.32` are used — at
`dx ≲ 0.08` with a fixed `dt` the whistler CFL is violated and the field blows up
(`recommended_dt` computes the stable step; see §10.3).

## Result (honest)

**Gated (pass):** β_p=1.0 resolution-independence (relative change < 0.5); β_p=0.2
gradient rises with refinement; cold steepens more than hot. **Skip (informational):**
the absolute gradient magnitudes are recorded in the notes (β_p=1.0 ≈ 3.6→4.7 vs
Hellinger ~5; β_p=0.2 ≈ 5.4→11.3 rising vs Hellinger 5→13→25) — same trend, magnitudes
within a factor of order one.
