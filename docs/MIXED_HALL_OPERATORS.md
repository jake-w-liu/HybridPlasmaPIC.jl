# Mixed Hall operators and constrained field projection

This source prerequisite supplies CPU Float32/64 operators and a raw induction
rate on an SBP x interval with periodic transverse axes. It owns no shock time
step, particles, reservoir, SAT, clock or RNG. A complete multidimensional
open-boundary full-Hall shock family remains open. Source/caller review is
separate from numerical acceptance: no new tests or execution were performed
under the current implementation-first instruction.

`MixedShockOperators(n, L)` supports D2/3 and always three vector components.
The x spacing is Lx/(nx−1), with nx≥3; transverse spacing is Ld/nd. The x
operator is registered SpectralOperators `SBP1D`; cached
`FourierDerivYWorkspace`s act on packed transverse fibres. Even Nyquist odd
symbols and singleton axes retain the dependency's zero derivative. Plans are
cached, but packing adds O(cells) copies and measured allocation/cost is pending.
No new dependency or SpectralOperators floor is required.

`mixed_gradient!`, `mixed_curl!`, `mixed_divergence!` and
`mixed_derivative_transpose!` stage their complete outputs. Curl geometry is
shared with the reduced 3D shock wrapper, preserving that wrapper's derivative
and subtraction order. The reduced induction equations, direct SAT, defaults
and checkpoint representation are unchanged by this prerequisite.

All supplied arrays must be ordinary one-based Arrays at grid precision. Roles
must be independent, including inactive components, operator metadata and live
workspace buffers. Storage checks use `Base.mightalias`; callers must supply
arrays whose overlap it reports. Metadata must remain unchanged. One operator
and workspace serve one caller at a time. Every public output is disjoint from
inputs/private work; failure before publication may change private work only.
Nearest rounding and gradual underflow are required. Materialized derivatives,
transposes and constitutive/rate intermediates must be finite. A later small
factor does not rescue a rejected materialized intermediate.
The coordinate expressions also retain their materialized range boundary;
complete original-geometry range closure is not implied.

The SBP norm is Hx=dx·Qx, Qx=diag(1/2,1,…,1,1/2). The ordinary Euclidean
transpose uses the existing forward derivative:

```text
Dxᵀ a = −Qx Dx(Qx⁻¹a) + (E/dx) Qx⁻¹a
E = diag(−1,0,…,0,1),    Dyᵀ=−Dy, Dzᵀ=−Dz.
```

This is an algebraic identity, not a guarantee of native FFT/row symmetry.
Intermediate Q-inverse and face-lift range is part of the admission boundary.

`MixedHallWorkspace(op)` and
`mixed_hall_rhs!(out, B, n, ui, physics, op, work)` supply raw isothermal
induction. The caller provides positive-ion collection-normalized n and ui.
`physics` must be `HybridModel(IsothermalElectrons{T}(...))` with original
η=ηH=de2=0. Require finite n≥nfloor without clipping. The method retains
canonical Te*n → mixed gradient → represented reciprocal, J=curl B,
completed canonical Ohm rows, and −curl E. A canonical original-product
source-risk event rejects: a finite completed local E does not establish a
representable differentiated original source. No periodic spectral recovery,
nonlinear pressure-chain rule or hidden projection is substituted.

The returned E/J/pe/gradp/ninv/rate are borrowed private fields, valid until
workspace reuse. `mixed_magnetic_rate_budget(B,E,J,rate,op)` and the RHS report
evaluate exact rational products of finite stored values and stored quadrature:
magnetic power, J·E, both outward Poynting faces and their signed identity
residual. The residual is reported rather than assumed zero. This does not
certify the supplied J equals curl B or the supplied rate equals −curl E;
the raw RHS actually constructs those fields. No full coupled energy theorem
or pressure/Hall boundary stability follows from this ledger.

`MixedHallProjectionWorkspace(op)` serves `project_mixed_field!` and
`project_mixed_rate!`. They target the restrictive perpendicular space

```text
V = {v: div v=0, vx(left/right)=0, (curl v)x(left/right)=0}.
```

It includes nontrivial transverse modes and arbitrary x-dependent tangential
fields in the transverse derivative nullspace. It is not an arbitrary-normal-B
or transparent Hall outflow law. All represented DC/Nyquist nulls and dependent
constraint rows remain; no inverse wavenumber or rank cutoff removes them.

Projection first eliminates normal face DOFs structurally (Rn). With
Gr=(dx·div, dx·left/right normal curl) and Q=I3⊗Qx, it solves the consistent
singular equation `Gr Rn Q⁻¹ Grᵀ λ = Gr Rn v` using cold, unshifted,
unpreconditioned MINRES-QLP. It builds the candidate and complete correction
privately. O(cells + constraint rows) storage includes six full solver vectors,
empty optional/history vectors, and fixed field/face/FFT packing buffers.
No assembled matrix, inverse, iteration history or checkpointed warm state is
used. Every matvec performs multiple mixed transforms. Convergence, native
adjoint conditioning, allocation and cost are unmeasured.

A solver success flag is required and insufficient: every solve is followed by
fresh exact stored-row true-residual comparison and direct divergence/normal
current checks. Normal faces must be stored zero. Source-scaled direct gates
use stored-symbol bounds `gdiv=dx*(2/dx+Σ|ky|+Σ|kz|)` and
`gj=dx*(Σ|ky|+Σ|kz|)`; these are scales, not FFT error theorems. Reports include
actual residuals and thresholds. Exact rational comparisons prevent overflowed
tolerances or underflowed scalar subtractions from silently accepting a solve.

Field projection also gates Q-orthogonality and energy nonincrease within its
explicit work budget, reporting the signed energy identity and subtraction
rounding defect. Rate projection separately admits B's constraints and gates
the actual represented change in dot_Q(B,rate). Exact power preservation holds
only for B truly in V and an exact projection; state and solve defects are not
silently treated as zero. Defaults are relative 100eps(T), absolute zero,
maxiter=200. Work/state controls are separately exposed; no measured adequacy
of these defaults is asserted. A failed gate raises `MixedHallProjectionError`
with actual residual/status context before any public output copy.

The global physical norm is the normalized Q norm times the common positive
stored cell volume. Projection reports use normalized Q units: divergence and
normal-current gates have dx-scaled field/rate units; work tolerances have
field² (field projection) or field-times-rate units (rate projection). Magnetic
face-work reports include the exact stored volume/area factors separately.
Only the correction arrays in a successful report are borrowed; scalar ledger
values are detached exact rationals. Rebuild workspaces from geometry after
restart; the actual evolving caller remains responsible for persisting its
physical fields and phase.

Package load, native Fourier/adjoint accuracy, mode/nullspace convergence,
rollback, compiler/JET/format/CI and resource acceptance remain deferred.
A full coupled Hall shock still needs an explicit physical boundary/SAT,
pressure/free-energy work, reservoir/newborn phase, stiffness/time controller
and complete atomic source/RNG/ledger persistence.
