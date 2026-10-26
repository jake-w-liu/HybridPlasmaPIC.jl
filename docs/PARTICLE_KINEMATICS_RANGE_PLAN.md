# Particle kinematics normwise range repair

Status: the normwise Boris/phase/coordinate range chain landed in `cd55831`.
Fresh near-light correctness work remains in the [audit tracker](DEEP_DEBUG_TRACKER.md).
The implemented chain uses the
[normwise standard](DEEP_DEBUG_TRACKER.md#numerical-accuracy-standard-user-decision-2026-10-26).

## Defects repaired

Finite inputs could produce NaN/Inf or lose kicks through intermediate range
loss: `h = qm*dt/2` overflowing in a zero field (`Inf*0 = NaN`), an overflowing
`|t|²` setting `f = 0` and dropping the rotation, an underflowed `qm*dt`
dropping a representable kick, a drift `x + dt*v` overflowing although the exact
coordinate is finite (shock steppers then absorbed the particle), the same
issues in Euler priming/predictor half-kicks, and subnormal `dt` (or `dt/2`)
erasing half-step kicks and midpoints.

## Contract

1. Well-conditioned ordinary inputs keep the native operation order: results are bitwise those
   of the previous formulas.
2. A guard over the native intermediates flags range loss; only flagged
   Float16/Float32/Float64 evaluations repeat the same expression wider and
   round once to `T`: Float16 → Float32, Float32 → Float64, Float64 → exact
   `Rational{BigInt}`. Relativistic expressions (square roots) use 256-bit
   BigFloat for unresolved flagged inputs with the exact Lorentz gap
   `1 - |v|²/c²`, because `(|u|/c)²` can exceed the next IEEE format for
   extreme `c`. Bounded IEEE priming first attempts a residual-certified
   double-word evaluation of that same expression.
   Other number types keep native semantics; non-finite inputs propagate.
3. Accuracy is a few ulps of the operation scale (`|v|`, `|h E|` for kicks,
   `|x| + |dt v|` for drifts); componentwise accuracy under cancellation is
   not certified in general. Relativistic priming whose native output has rounded normalized squared
   speed `>1/2` uses the repair/publication tier: a normwise one-ulp stage
   error there can amplify the next Lorentz gap. Every admitted fast repair
   is certified to lie in one rounding cell; unresolved cells use BigFloat.
4. A wide result that does not fit `T` is returned as `±Inf` for the caller's
   existing policy; a leapfrog/ES priming back-up that does not fit throws an
   `ArgumentError`.
5. Public steppers and pushers reject `0 < |dt| < 2*floatmin(T)` (and EM
   substeps `dt/n_sub` below that bound) with an `ArgumentError` before any
   mutation; `dt == 0` is unchanged.

The optional gathered-push midpoint buffers have a separate
[storage admission contract](GATHERED_MIDPOINT_STORAGE_PLAN.md): they must
not overlap particle arrays, gathered fields or each other. CPU and Metal use
the same preflight check before any particle update.

## Implementation

| Expression | Owner | Guard (beyond finite outputs) |
| --- | --- | --- |
| Boris kick | `boris_kick`, `src/particles/boris.jl` | `h` normal or legitimately zero, `|t|²` finite, `f ≥ floatmin`, `S = max(|v|,|hE|)` zero or `≥ 256 floatmin (1+τ)` |
| Relativistic Boris | `_boris_rel`, `src/models/full_pic.jl` | as above on momentum, plus `c²` normal, finite `γ⁻`, `γ⁺`, and a resolved Lorentz factor |
| Euler priming/predictor | `_euler_kick`, `src/particles/kinematics_range.jl` (used by `_prime_leapfrog!`, `_predict_half_kick!`, `_prime_empic1d!`) | `h*qm` normal or legitimately zero, finite `E + u×B`, scale `≥ 256 floatmin (1+|h qm|)` |
| Relativistic priming | `_rel_prime_in_range`/`_rel_prime_repair` | `c²` normal, `h` normal or zero, finite acceleration and `γ`, resolved Lorentz factor and output squared speed `≤ c²/2`, scale `≥ 256 floatmin (1+|h|)(1+|qm|)` |
| ES priming | `_es_prime_component` | `h*qm` normal or zero, scale `≥ 256 floatmin` |
| Drift/midpoint | `_particle_coordinate` (push_uniform!/push_gathered!, hybrid, EMPIC1D/EMPIC, ES, 7 shock updates) | finite result only (subnormal sums are exact) |

The wide-format claims (no overflow, harmless underflow) are derived in the
comments next to each repair. `dt` validation is `_require_resolved_dt` in
`src/utils/validation.jl`, called by `_validated_step_dt` (`step!`,
`step_camcl!`, `mpi_step!`, `step_shock!`, `step_leroy_shock!`,
`step_shock2d!`, `step_shock3d!`), `step_espic!` (three methods),
`step_empic!` (both models, also on `dt/n_sub`), `step_hall_mhd!`, and the CPU
and Metal `push_uniform!`/`push_gathered!`. It is not added to
`_validated_nonnegative_dt`, whose collision, ionization, radiation,
electron-pressure and face-injection callers have no half-step terms
(`test_inject.jl` exercises a subnormal injection `dt`).

Batch owners (`push_uniform!`, `push_gathered!`, the hybrid kick/drift, EMPIC
and EMPIC1D velocity pushes and drifts, ES priming/kick/drift, leapfrog and
EMPIC1D priming, the predictor and the shock steppers) first run one fused
reduction without calls or writes: max |v|, |E|, |B| (|x| for drifts) and the
smallest nonzero guard scale. A global certificate derived from these bounds
(comments at `_boris_batch_certificate` and siblings in `kinematics_range.jl`)
proves that every particle passes its guard; the batch then runs the native
historical loop with the same arithmetic order. Otherwise the guarded
loop runs on the untouched inputs. The scalar `boris_kick` keeps its per-call
guard. Shock loops gather per particle, so their certificate uses grid maxima
(nonnegative interpolation weights) and conservatively rejects particles at rest
in a nonzero electric field.

Metal kernels call the native-only `_boris_kick_native` and
`_particle_coordinate_native` (no guard or fallback on the device). Before a
uniform launch the host checks the velocity-independent part of the guard
(`h` normal or zero, finite `h E`, finite `|t|²`, `f ≥ floatmin`); velocity
dependent range loss on the device (which also flushes subnormals) is
unsupported.

## Verification

`test/test_particle_kinematics_range.jl` retains 10⁵ ordinary inputs per type.
Boris, relativistic Boris, nonrelativistic priming, predictor, ES priming and
drift compare bitwise against copies of the previous formulas. For both
relativistic priming owners, an independent 4096-bit reference classifies the
physical output: squared speed/c² ≤1/2 must retain historical bits; the faster
slice must equal the independently rounded stage. This is an explicit change
from the original blanket priming bit-equality contract; those original
failures are retained in the takeover evidence. Range
witnesses from the certified chain agree with exact Rational (Cayley) or
4096-bit BigFloat references within 4 ulps of the operation scale (Float64
exact fallbacks equal the correctly rounded result); unresolved `dt` is
rejected without state change on every covered entry point; zero-field
extreme-charge equilibria stay unchanged for hybrid, CAM-CL, ES, EMPIC1D and
EMPIC; shock steppers keep and reflect a particle whose rounded drift product
overflows.

## Performance (Apple M5, Julia 1.12.7, shared machine at load 27–70)

The original pre-takeover measurement recorded the parent's vector IR (120 `push_uniform!`,
62 EMPIC kicks, 34 priming, 34 predictor operations); the pre-passes vectorize.
The scalar `boris_kick` keeps its per-call guard, so a loop over it stays scalar
(`boris_bench.jl`: 14.2–14.7 against 3.6–3.8 / 7.7–7.9 ns per push, about
3.9× / 1.9×). `push_uniform!` costs about 1.3× the parent (one extra fused
read of velocities and positions against a very cheap native kernel);
`push_gathered!`, hybrid 2D and ES 2D steps are within the load noise.

## Lorentz-gap conditioning follow-up

The [takeover audit](KINEMATICS_AUDIT_20261027.md) records a reproduced error
in accepted near-light-speed states: finite native intermediates did not bound
the cancellation in `1-β²`. Relativistic Boris and both priming paths now
resolve the stored Lorentz gap with shared scaled error-free arithmetic when
rounded `β² > 1/2`; unresolved gaps and other range failures retain the
existing exact-gap replay. Their batch
certificates require every velocity component to be at most `c/4` before
bypassing scalar guards; other batches retain per-particle admission. The
accepted input domain and well-conditioned ordinary operation order are
unchanged. The scaled factor rejects gaps below its resolved arithmetic range; it never
clips a gap or velocity. Priming also checks the output conditioning: a
native published stage with rounded normalized squared speed `>1/2` enters
a bounded residual-certified DD tier.
Exact expansion residuals and outward interval propagation certify each
component's rounding cell, including the uncertainty of the final gamma.
Failed arithmetic/cell admission retains the existing BigFloat replay.
The old blanket ordinary-bit fixture is not claimed passed: high-output-speed
Float32/64 priming changes ordinary results by a few ulps. The revised test
keeps every original input and unchanged whole-step tolerances.
High-speed allocations and cost are measured
separately; no performance acceptance is claimed while the machine is shared.
