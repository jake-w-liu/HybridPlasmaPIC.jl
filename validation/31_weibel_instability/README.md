# Case 31 — Weibel (current-filamentation) instability

Phase-1 **kinetic instability**, driven in the **full electromagnetic PIC** model
(`EMPIC`, kinetic electrons + immobile neutralizing ion background). Two
counter-streaming cold electron beams `±u₀ x̂` supply a velocity-space anisotropy;
the unstable `k ∥ ŷ` (⊥ the streaming) grows the out-of-plane field `B_z(y)` from
shot noise.

## Why full-PIC, not hybrid

The anisotropy instabilities (cases 29 firehose, 30 EMIC) run on the *hybrid* engine
because they are ion-scale and magnetized. The Weibel cannot: in the quasineutral
massless-electron hybrid model `B = 0` is an **exact fixed point** (with uniform
density the electric field is the curl-free `∇pₑ` term, so `−∇×E = 0`), and only the
*bulk* ion moments couple to the field — the two beams' transverse drift cancels and
a seeded `B` merely decays (verified: `B_z/seed` falls monotonically). The full-PIC
electron dynamics are essential, so this case uses `EMPIC`.

## Threshold (analytic)
With `A = (u₀/vth)²` (effective streaming anisotropy `⟨vₓ²⟩/⟨v_y²⟩ = 1 + A`):

    unstable ⇔ A > 1   (bimodal counter-streaming; Weibel 1959, Fried 1959)

## What is tested (`weibel_growth`)
Unstable `u₀=0.6, vth=0.1` → `A=36 ≫ 1`; stable `u₀=0` → single Maxwellian, `A=0`.
`B_z` grows exponentially (`½∫B_z²` from ~1e-4 to ~0.16 over `t≈20`), then saturates as
the beams isotropize.

**Gated (pass):** the counter-streaming case grows `B_z` energy above `0.02` (measured
~0.16); the no-stream case stays `≤0.005` (~5e-4, shot noise); the threshold flag
matches. ~300× separation.

`default=false` (research validation; ~10 s).
