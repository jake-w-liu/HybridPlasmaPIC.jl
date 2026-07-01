# Case 29 — Parallel firehose instability

Phase-1 **kinetic instability** validated by reusing the hybrid PIC engine (no new
solver). An anisotropic bi-Maxwellian ion distribution in a uniform `B₀ = x̂`.

## Threshold (analytic)
Fluid electrons are isotropic, so the firehose threshold is set by the ions alone:

    unstable ⇔ β_∥ − β_⊥ > 2 ⇔ T_∥ − T_⊥ > B₀² ⇔ vth_∥² − vth_⊥² > 1   (n₀=B₀=μ₀=1)

(Gary 1993, *Theory of Space Plasma Microinstabilities*.)

## What is tested (`firehose_growth`)
Two runs: **unstable** `vth=(1.5, 0.3)` → `β_∥−β_⊥ = 4.32`; **sub-threshold** `vth=(1.0, 0.8)`
→ `0.72`. The transverse fluctuation energy `δB_⊥ = (B_y,B_z)` as a fraction of the
background magnetic energy is measured.

**Gated (pass):** the unstable case grows `δB_⊥` to ≥10% of B₀ energy (measured ~0.54);
the sub-threshold case stays ≤2% (measured ~0.002, particle-noise floor); and the
theory-threshold flag matches the analytic condition. ~250× separation between them.

`default=false` (research validation; ~8 s). Run via `run_validation.jl` or directly.
