# Reference-unit normalization range repair

Status: **NORM-001 confirmed on Julia 1.12.7 and 1.13.2; repair not implemented.**
This is a prerequisite for accepting the new local collision normalization
helper, and also affects existing public unit conversion and RayCon adapters.

## Evidence and canonical ownership

The keyword `PlasmaUnits` constructor accepts finite positive density, field
and mass references. [normalization.jl](../src/verification/normalization.jl)
owns `alfven_speed`, `gyrofrequency`, `inertial_length` and `_scale`, used by
`to_SI` / `to_normalized`. The RayCon normalized equilibrium and ray adapters
call these accessors in [raycon_normalized.jl](../src/raycon/raycon_normalized.jl).
The isolated collision candidate also derives its setup scales from them.

**VERIFIED:** the unchanged source at `6e7051b` produces 45 passing assertions
and 25 failures per runtime in a four-profile witness. Its independent 512-bit
reference combines the stored physical constants and inputs before conversion
to Float64. All three accessors, eight conversion kinds and representable
inverse conversions are checked. Ordinary solar-wind values are controls.

| References `(n0, B0, mi)` | Public result | Independent reference |
| --- | --- | --- |
| `(1e300, 1e200, 1e100)` | Alfvén speed `0.0` | `892.0620578335771` |
| `(1e-300, 1e-200, 1e-100)` | Alfvén speed `Inf` | `892.0620578335771` |
| `(1e300, 1e-300, 1e100)` | Inertial length `NaN` | `5.567813428950414e-79` |
| `(1e-300, 1e-200, 1e-100)` | Current scale `0.0` | Approximately `1.429241e-316` |

The checked arithmetic forms `mu0*n0*mi` before its square root (line 61),
then divides independently rounded Alfvén and gyro scales (line 75). `_scale`
composes those results further. These operations are the verified failing
owners; the complete-scale oracle and recorded outputs establish the finite
results that they lose. True output overflow is recorded separately.

Evidence is retained under
`HybridPlasmaPIC-audit-evidence/2026-10-23/normalization-range/`:
`source.json`, `witness.jl`, `baseline-commands.json`, both full logs/TOMLs and
`baseline-review.json`. The source hash is checked by each executed witness.

## Repair and acceptance tracker

- [x] Trace accessors, public scalar/array conversions, RayCon callers and tests.
- [x] Reproduce complete-scale failures against independent references on both
  supported Julia versions; preserve the ordinary controls.
- [x] Characterize raw positional constructors, stored floating types and invalid
  constants: 112 records per runtime retain the mixed-positional `MethodError`
  and show that raw typed constructors accept zero, negative and nonfinite
  constants. Eight huge-exponent BigFloat profiles retain 20 passes / 16 failures.
  Precision, rounding and exponent limits are restored after characterization.
- [ ] Choose and document constructor/zero-sign compatibility before changing
  either constructor. Complete the [shared BigFloat materializer repair](EXACT_BIGFLOAT_CONVERSION_PLAN.md)
  before relying on exact complete-scale conversion.
- [ ] Repair the canonical complete scale, retaining ordinary outputs and types
  where the existing evaluation is safe. Reuse exact conversion/square-root
  ownership where applicable. Bound exceptional work; avoid expanding enormous
  BigFloat exponents into unbounded rational integers or changing global MPFR
  state. Do not patch only the collision normalization helper.
- [ ] Audit conversion of the caller's value separately: a finite result may
  involve a scale outside the destination range. Decide the full scalar/array
  contract before treating a scale-only repair as complete conversion coverage.
- [ ] Add independent boundary and exponent matrices, subnormal outputs,
  forward/inverse scalar and array checks, and retained-precision inputs.
  Reject zero-clamping, independent-factor rounding and round-trip-only patches
  with explicit negative controls.
- [ ] Run the original witness unchanged, existing normalization/RayCon tests,
  collision unit-equivalence tests on their final source, concrete/full compiler
  comparisons and threaded state checks on both Julia versions.
- [ ] Freeze ordinary and exceptional complete-call resource/cost profiles before
  measuring. Keep all failed measurements; no timing retries or relaxed gates.
- [ ] Update API docs, README and changelog; review the final diff; push this fix
  separately to main and require all six exact-SHA CI jobs.

This plan does not establish correctness of the unaccepted collision candidate
or close CAP-001. Its corrected unit-equivalence and coupled-channel tests still
require execution on the final source.
