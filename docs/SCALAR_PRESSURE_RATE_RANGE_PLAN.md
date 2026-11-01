# Coupled scalar pressure rate range

Updated: 2026-11-01. The implementation and essential regressions are prepared;
GitHub package, compiler, formatter and MPI validation remain pending.

The zero-hyperresistivity RK stage shares one signed pressure equation between
`EvolvingElectrons` and `TransportElectrons`:

```
rate = -ue·grad(pe) - gamma*pe*div(ue) + (gamma-1)*eta*|J|².
```

`EvolvingElectrons` requires `eta=0`. The stored `ue=ui-J/max(n,nfloor)` convention
and rounded `gamma-1` heat coefficient are retained. This update completes
pressure sources; it adds no heating to the ideal closure. Standalone
`advance_electron_pressure!` remains its separate Euler building block.

## Source-supported failures

A four-point centered derivative supplies independent exact dyadic witnesses.
For Float32/64, positive finite pressure, finite flow/current and finite required
physical E/magnetic rate can give separately overflowing ideal and Ohmic terms
whose complete pressure rate is exactly zero. Ideal advection and compression
can also cancel exactly without heat. The old stage rejected their nonfinite
intermediate result.

A separate half-quantum case has ideal and Ohmic contributions that round to
zero individually, while their complete sum is one representable subnormal.
Normal heat in other cells does not remove that local product-range event.
Another case loses a pressure derivative before a large finite flow rescues a
normal advective rate. The tests retain the exact four-point matrix and original
factors; Julia execution of these prepared cases is pending on GitHub.

## Canonical implementation

The ordinary stage retains its pressure formula, then its transport field/heat
addition. A range scan requests completion for nonfinite derivatives/rates,
nonzero product range loss and an active tiny whole rate. The existing S2
pressure-derivative routing is shared and checks original variation and the
stored gradient in each differentiated axis, alongside FFT magnitude checks.
A normal gradient in another axis cannot hide a consumed-axis range event.
These are conservative routing conditions, not proofs of
a nonzero source or a general FFT error bound.

Before exact extraction, finite original pressure, spatial electron-velocity
components, coefficients and active Ohmic current are required. Unused transverse
velocity and inactive current are not new operands. The existing `_complete_transport_pressure_rhs!`
then combines the signed sources using normalized derivative data and exact
represented coefficients with `etaH=0`. No current Laplacian is computed or
required for this branch. The S2 owner retains its default diagnostic context;
the shared ideal/Ohmic stage reports a complete electron pressure rate.

Complete zero and signed cancellation are valid. A genuinely overflowing
complete rate, or a nonzero complete norm that rounds wholly to zero, rejects
through the existing complete-rate owner. Required physical E/J/force checks and
public step rollback remain in their existing owners; raw stage outputs are
scratch and have no new atomic-publication guarantee. No persistent arrays,
checkpoint schema or model parameter is added.

## Progress and limits

- [x] Inspect canonical pressure, heat, stage, field, transaction and S2 owners.
- [x] Derive small exact overflow, half-quantum and derivative-rescue witnesses.
- [x] Implement shared cold routing and essential independent regressions.
- [ ] Pass exact-revision GitHub package/compiler, formatting and MPI checks.
- [ ] Measure complete public costs before allocation or performance claims.

The D2 axis-permutation regression also checks reconstruction of a tiny consumed
gradient alongside a normal transverse gradient. Its runtime result is pending:
the shared scaled derivative currently transforms all axes before applying the
directional symbol. Per-axis routing alone does not certify that reconstruction.

General Fourier conditioning, generation of `ue` and other derivative sources,
full hardware-mode pipelines and dt-scaled RK representations remain separate
boundaries. This source review does not certify every extreme coupled trajectory.
The [audit tracker](DEEP_DEBUG_TRACKER.md) owns integration and CI acceptance.
