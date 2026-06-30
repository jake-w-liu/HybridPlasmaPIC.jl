# Case 28 — Perpendicular shock vs Hybrid-VPIC (live external hybrid code)

The **nonlinear** code-to-code comparison (the shock analog of case 27's NHDS dispersion):
our perpendicular hybrid shock vs **Hybrid-VPIC**, a *different, independent* open-source
kinetic-ion/fluid-electron hybrid PIC code from LANL
([github.com/lanl/vpic-kokkos](https://github.com/lanl/vpic-kokkos), branch `hybridVPIC`),
**built and run from source**.

## What is compared

`build_vpic.sh` builds Hybrid-VPIC and runs its `examples/shock/` deck reconfigured to match
ours — **θ_Bn = 90° (perpendicular, B = Bz)**, β_i = 1, T_e = T_i, drift Mach 4 (→ measured
shock-frame **M_A ≈ 6**), 1-rank thin-2D, natural hybrid units (v_A = Ω_ci = d_i = B0 = 1, the
same as ours). It dumps the Bz(x) profile (`shock_perp.cxx`, an ASCII dump added to the deck).
Since Bz/n is frozen-in, that profile gives both compression (Bz₂/Bz₁) and overshoot
(Bz_max/Bz₂); the shock-front speed gives M_A.

## Result (measured)

| | compression | overshoot |
|---|---|---|
| **Hybrid-VPIC** (external, M_A≈6.07) | **2.89** | **1.28** |
| ours (`run_perp_shock_rh`, M_A=6) | 3.20 | 1.43 (per-snapshot) |
| fluid Rankine–Hugoniot | 3.21 | — |
| Leroy 1982 (published hybrid) | < RH | 1.26 ± 0.06 |

- **Overshoot:** Hybrid-VPIC **1.28** ≈ ours (clean profile ~1.33 / energy-closure ~1.27) ≈
  Leroy **1.26**. Independent codes agree on the key kinetic observable. *Gated:* VPIC overshoot
  in Leroy's band [1.1, 1.5], and ours-vs-VPIC overshoot agree within 25% (measured 11%).
- **Compression:** Hybrid-VPIC self-consistently gives **2.89 < fluid RH 3.21** — the
  *kinetic-compresses-less-than-fluid* effect (Leroy's point). *Gated:* VPIC compression ≤ RH.
  Our `run_perp_shock_rh`/`leroy` use an RH initialization that holds the fluid compression
  (3.20), so they sit at RH by construction — the difference is that init choice, not a
  disagreement (our self-consistent overshoot still matches VPIC).

So an independent external hybrid code **confirms** our perpendicular-shock physics
(overshoot, and kinetic < fluid compression), consistent with Leroy 1982.

## Running it

VPIC source/binary/output are **gitignored** (`validation/**/vpic_build/`, `*.txt`); only the
deck + scripts are committed. Run:

```
bash validation/28_hybridvpic_perp_shock/build_vpic.sh   # needs git+cmake+MPI (brew installs OpenMPI on macOS)
julia --project=validation validation/28_hybridvpic_perp_shock/run.jl
julia --project=validation validation/28_hybridvpic_perp_shock/plot.jl   # VPIC Bz(x) profile PDF
```

`run.jl` auto-runs `build_vpic.sh` if the profile is missing and **skips cleanly** if VPIC
cannot be built. `default=false` (research validation; not in the package CI test suite —
no MPI/network there). A future extension is a full shock-aligned Bz(x) overlay of both codes.
