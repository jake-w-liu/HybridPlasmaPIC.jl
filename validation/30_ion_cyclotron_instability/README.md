# Case 30 — Electromagnetic ion-cyclotron (EMIC) instability

Phase-1 **kinetic instability**, the `T_⊥ > T_∥` counterpart of the firehose (case 29),
validated by reusing the hybrid PIC engine. Anisotropic bi-Maxwellian ions in a uniform
`B₀ = x̂`; same parallel run + transverse-δB diagnostic, opposite anisotropy.

## Threshold (analytic)
With `A = T_⊥/T_∥ − 1` and `β_∥ = 2 T_∥`:

    unstable ⇔ A > 0.43 / β_∥^0.43   (Gary 1993 marginal-stability fit)

## What is tested (`ion_cyclotron_growth`)
`ion_cyclotron_growth` reuses `firehose_growth` for the run and applies the EMIC
threshold. Unstable `vth=(0.4,1.3)` → `A=9.6 ≫ threshold`; sub-threshold `vth=(0.8,0.9)`
→ `A=0.27 < threshold`.

**Gated (pass):** the unstable case grows `δB_⊥` above the noise floor (measured ~0.16 of
B₀ energy); the sub-threshold case stays ≤2% (~0.001, noise); the threshold flag matches
Gary's condition. ~110× separation.

`default=false` (research validation; ~10 s).
