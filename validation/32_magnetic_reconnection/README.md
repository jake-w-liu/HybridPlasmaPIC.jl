# Case 32 — Collisionless magnetic reconnection (Harris-sheet tearing)

Phase-1 **kinetic phenomenon**, driven on the **2-D hybrid** engine. A periodic
double-Harris current sheet is perturbed and reconnects via the `m=1` tearing mode —
the iconic process behind magnetospheric substorms, solar flares, and sawtooth
crashes.

## Equilibrium
Reversing field with sheets at `y₁=L_y/4, y₂=3L_y/4`:

    B_x(y) = B₀[tanh((y−y₁)/λ) − tanh((y−y₂)/λ) − 1]

in pressure balance `n(y)(Tᵢ+Tₑ) + B_x²/2 = const` via the Harris density
`n(y) = n_b + n₀[sech²((y−y₁)/λ)+sech²((y−y₂)/λ)]`, `n₀ = B₀²/(2(Tᵢ+Tₑ))`, deposited
through per-particle weights. The electron fluid carries the sheet current `J = ∇×B`.
A **divergence-free** flux-function seed `δψ = δ cos(kₓx)[sech²₁−sech²₂]` perturbs both
`δB_x=∂_yδψ` and `δB_y=−∂_xδψ`, launching the island.

## Diagnostic (why the coherent mode)
The reconnected flux is the **coherent `m=1` power** of `B_y`, `Σ_y|B̂_y(kₓ,y)|²` via
`rfft` along x. This is essential: the *total* `B_y` energy is dominated by broadband
particle shot noise (~0.07 in both the sheet and a uniform control), which hides the
island; the `m=1` projection isolates the growing tearing mode.

## Threshold (analytic)
    tearing-unstable ⇔ sheet && kₓλ < 1   (Furth-Killeen-Rosenbluth 1963)

Here `kₓλ = (2π/L_x)·λ = 0.123 < 1` → unstable.

## What is tested (`reconnection_growth`)
Sheet (`sheet=true`): the `m=1` mode grows ~exponentially (measured ~13× the seed).
Uniform control (`sheet=false`, `B_x=B₀`, no free energy): the same seed stays flat
(~1.1×). **Gated (pass):** sheet growth ≥3×; uniform ≤2×; threshold flag matches.
~12× separation.

`default=false` (research validation; ~90 s for both runs at 64×128).
