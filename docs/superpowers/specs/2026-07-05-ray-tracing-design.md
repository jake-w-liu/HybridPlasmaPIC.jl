# Ray Tracing (Geometric Optics) for Hybrid Wave Branches — Design

Date: 2026-07-05
Status: approved for implementation (autonomous session; assumptions stated inline)

## Purpose

Add WKB / geometric-optics **ray tracing** as a supporting method: trace wave-packet
trajectories `x(t)`, `k(t)` of the warm Hall-MHD (hybrid, massless-electron,
quasineutral) wave branches through a smoothly varying plasma — either an analytic
medium or a snapshot of the simulation's own grid fields `n(x)`, `B(x)`.

This is a diagnostic/supporting method, not a field solver: rays are traced through a
**frozen** snapshot (valid when the wave frequency is large compared to the background
evolution rate) on a single rank (post-processing on gathered global fields; MPI/GPU out
of scope).

## Physics (verified before design was finalized)

Units are the code's Ω_ci-normalization (length `d_i`, time `Ω_ci⁻¹`, velocity `v_A`,
`B` in `B0`, `n` in `n0`). The scalar dispersion relation for the warm Hall-MHD model —
the same physics as the HYB-006 eigenvalue oracle (`test/oracles/hybrid_dispersion_oracle.jl`)
— is, with `u = ω²`:

```
D(ω, k; n, B) = (u − A)(u² − P·u + Q) − R·u·(u − S)

A = k∥² v_A²          P = k²(v_A² + c_s²)      Q = k² k∥² v_A² c_s²
R = k² k∥² B²/n²      S = k² c_s²
v_A² = B²/n           k∥² = (k·B)²/B²          B² = |B|²
c_s²(n) = γ_e T_e n^{γ_e−1} + γ_i T_i n^{γ_i−1}   (combined sound speed, cf. HYB-005)
```

`D = 0` is a cubic in `u` whose three non-negative roots are the slow, intermediate
(shear-Alfvén / ion-cyclotron), and fast (magnetosonic / whistler) branches at arbitrary
angle, including the Hall term.

**Verification already performed** (scratchpad `verify_dispersion.jl`, run 2026-07-05):
roots of this cubic match `HybridDispersionOracle.dispersion_frequencies` over 499
random `(k, B, n, c_s)` draws to worst relative error `3.9e-11`, and reproduce the
closed-form parallel R/L modes `ω = ±k²/2 + k√(1+k²/4)` and perpendicular fast mode
`ω = k√(v_A²+c_s²)` to `1e-10`. The same check ships as test RAY-001.

## Ray equations

For time-stationary media, `ω` is constant along a ray and

```
dx/dt = −(∂D/∂k)/(∂D/∂ω) = v_g        dk/dt = +(∂D/∂x)/(∂D/∂ω)
```

- `∂D/∂ω`, `∂D/∂k_j`, `∂D/∂n`, `∂D/∂B_j`: **complex-step differentiation** of the single
  generic implementation of `D` (machine-precision, no step-size tuning; `D` is
  polynomial/rational plus `n^(γ−1)`, all complex-analytic for `Re n > 0`).
- `∂D/∂x_i = (∂D/∂n)·∂n/∂x_i + Σ_j (∂D/∂B_j)·∂B_j/∂x_i` — chain rule through the medium
  gradients supplied by the medium type.
- Integrator: classical RK4 with fixed user-chosen `dt` (matches the codebase's RK4
  use; no new dependencies).
- The relative residual `|D|/scale` is recorded each step as the standard ray-code
  accuracy diagnostic; the trace stops with a status flag on non-finite state,
  `∂D/∂ω = 0` (caustic/resonance), invalid medium (`n ≤ 0`, `|B| ≤ Bmin`), or residual
  drift beyond `residual_max`.

## Components (single file `src/diagnostics/ray_tracing.jl`)

1. **Dispersion core** (pointwise, medium-free):
   - `hybrid_wave_dispersion(ω, k, n, B; Te, γe, Ti, γi)` → `D` (generic over Complex).
   - `hybrid_wave_frequencies(k, n, B; ...)` → `NTuple{3}` of `ω ≥ 0` sorted ascending
     (`:slow`, `:intermediate`, `:fast`); cubic solved via 3×3 companion-matrix `eigvals`.
   - `hybrid_wavenumbers(ω, khat, n, B; ...)` → all positive `|k|` roots along direction
     `khat` at frequency `ω` (cubic in `|k|²`), each labeled with its branch — the
     standard way to launch rays at a prescribed frequency (e.g. from `antenna.jl`).
   - `wave_group_velocity(k, n, B; branch, ...)` → `(v_g, ω)` at a point.
2. **Media** (`RayMedium` abstract type):
   - `AnalyticRayMedium(nfun, Bfun; Te, γe, Ti, γi, h)` — user callables
     `nfun(x,y,z)::Real`, `Bfun(x,y,z)::NTuple{3}`; gradients by central finite
     difference with per-axis step `h·max(1,|xᵢ|)`, `h = cbrt(eps())` by default.
   - `GridRayMedium(g::FourierGrid{D}, n, B; Te, γe, Ti, γi)` for `D = 1,2,3` — takes a
     snapshot of grid `n` and `B = (Bx,By,Bz)`, precomputes **spectral** gradients
     (`deriv`) at construction, and CIC-interpolates values and gradients at ray
     positions (same node convention as `gather_at`/`deposit.jl`: node `i` at
     `(i−1)·dx`, periodic wrap). Rays carry 3-component positions and wavevectors
     always; the medium varies only along the grid's `D` axes (the package's
     "positions carry D coordinates, velocities carry 3" convention, applied to rays).
   - Sound-speed parameters map from the electron closures: `IsothermalElectrons(Te)` →
     `Te=Te, γe=1`; `PolytropicElectrons(pe0,n0,γ)` → `Te=pe0/n0^γ, γe=γ`. Warm ions via
     optional `Ti, γi`. CGL (anisotropic) is not supported by this scalar relation.
3. **Tracer**:
   - `trace_ray(med, x0, k0; branch=:fast, dt, nsteps, Bmin, residual_max)` →
     `(; t, x, k, vg, residual, ω, branch, status)` with `x, k, vg :: 3×(m+1)` matrices,
     truncated at early termination; `status ∈ (:ok, :nonfinite, :caustic,
     :invalid_medium, :residual)`.

Validation follows `utils/validation.jl` conventions (`_require_finite_real` etc.,
`ArgumentError` on bad input); construction validates array sizes against the grid,
finiteness, and `n > 0` everywhere.

## Alternatives considered

- **Full cold-plasma Stix dispersion** (electron-mass branches): rejected — inconsistent
  with the massless-electron model and its Ω_ci normalization; the hybrid relation covers
  every branch this code propagates (README HYB-002..006).
- **det(M − ωI) of the 7×7 oracle matrix as `D`**: rejected — complex-valued, spurious
  null modes, poor scaling for Hamiltonian derivatives; the closed-form cubic is exact
  and was verified against the oracle.
- **Adaptive ODE integration (DifferentialEquations.jl)**: rejected — new heavy
  dependency for a diagnostic; fixed-step RK4 matches codebase practice, and the
  recorded `|D|` residual exposes step-size inadequacy.
- **Analytic hand-derived Hamiltonian derivatives**: rejected in favor of complex-step
  on one generic `D` — one source of truth, machine precision, far smaller bug surface.

## Testing (test/test_ray_tracing.jl, RAY-001..009)

Oracle-based, matching repo culture: frequency roots vs HYB-006 oracle (random sweep);
factored/expanded consistency `D(ω_root) ≈ 0`; complex-step vs central-difference
derivatives; group velocity vs oracle `dω/dk`; uniform-medium straight-line rays
(analytic + grid); stratified perpendicular fast wave against the exact local relation
`k(x) = ω/√(v_A²(x)+c_s²(x))`; conservation of `k_y, k_z` in x-stratified media;
`hybrid_wavenumbers` round trip; validation-error paths; and 1D ≡ y-invariant-2D
dimension parametricity (DIM-001 style).

## Out of scope

MPI-distributed media, GPU, time-dependent media, mode conversion across branch
degeneracies (flagged via status, not followed), and amplitude transport along rays
(focusing/absorption) — the tracer computes trajectories, not wave energy budgets.
