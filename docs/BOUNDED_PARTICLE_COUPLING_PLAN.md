# Closed Yee particle charge, current and path work

`YeeChargeWorkspace`, `yee_particle_charge!` and `yee_longitudinal_impulse!`
provide the CPU Float32/Float64 charge/current prerequisite of CAP-002B on a
closed one-dimensional PEC interval. They do not select or advance a particle
trajectory. Coupled bounded electromagnetic PIC remains incomplete.

`yee_longitudinal_voltage(Ex,xold,xnew,g)` supplies the read-only scalar path
work prerequisite on the same closed interval. It does not choose a force or
advance particles.

An optional `YeeTriangularCloud(a)` shape now supplies static finite physical-width
charge/current/path work and a stationary work-adjoint field. Its source review is complete; essential regressions and exact-revision GitHub
acceptance remain pending.
It is a separate declared cloud model; point/CIC remains the default.

For N cells, nodal charge uses the CIC hats at the exact stored-input coordinate
`s=N*x/L`, including both wall nodes. It is extensive charge, without division
by dx. The longitudinal impulse integrates the oriented path overlap of the
same hats between caller-supplied endpoints, using the original species charge
and particle weight. It has N entries and requires no additional dt or dx factor.
In exact arithmetic, `Qnew-Qold+[I[1],diff(I)...,-I[end]]=0`. Independently rounded
charge and impulse have a normwise continuity residual, not bitwise continuity.

Accept one `ParticleSet`, a tuple or a one-based vector of species. Abstract
collection/component declarations and strided views remain supported; consumed
coordinates and weights have grid precision. Charge and weights must be finite,
weights nonnegative and endpoints in `[0,L]`. Mass, velocity and provenance
values are not operands, although all component shapes/storage are checked.
Outputs are ordinary vectors of grid precision. The output and both workspace
arrays must be independent of one another and all inputs; read-only inputs may
share storage. One workspace belongs to one call at a time.

All populations, endpoints, shapes and aliases are checked before scratch writes.
The default point/CIC native accumulation has bound `gamma*S`, where
`gamma=(P+10)*u/(1-(P+10)*u)`,
`u=2^-53`, P is particle count and S is the sum of absolute contributions to the
entry. Final conversion adds `eps(T)/2*(abs(exact)+gamma*S)` plus half a subnormal
quantum. Native arithmetic requires nearest hardware rounding, gradual underflow
and the count/coefficient/product/range checks. Otherwise the same kernel replays
exact rational original operands, including zero, order and geometry decisions.
Final exact conversion is nearest-even independently of hardware context. A complete
entry whose canonical rounding is nonfinite raises `ArgumentError` before output
publication. Scratch contents can change on arithmetic rejection and remain
reusable; invalid-input preflight preserves scratch too.

## Longitudinal path voltage

Let `C[p,e]` be the oriented normalized cell overlap already used by the impulse,
and let `V[p]=dx*sum(Ex[e]*C[p,e])`. The exact contraction identity is
`sum(q[p]*weight[p]*V[p])=dx*sum(Ex[e]*I[e])`. Supplying the exact temporal
mean `(Eold+Enew)/2` makes this the longitudinal part of the existing field
source-work expression. Constructing that mean, rounding the impulse and
rounding each voltage introduce their own errors; no bitwise equality between
independently rounded APIs is promised.

The Maxwell work metric uses stored `dx`, while endpoint geometry uses exact
`N*x/L`. For nonzero physical displacement, the derived work-adjoint average
is `V/(xnew-xold)`, including `N*dx/L`. A constant field therefore gives
`V=(N*dx/L)*Ex*(xnew-xold)` before rounding. The factor differs from one when
stored `L/N` is rounded; the operator preserves this existing discrete metric.
It includes longitudinal DC without imposing a fixed voltage or free plate
charge. Wall endpoints remain geometric endpoints, not absorption events.

The field is a CPU strided vector of grid precision with N finite entries.
Both scalar endpoints have that precision and lie in `[0,L]`. Validation
precedes calculation, including for a zero path. The complete expression uses
exact stored operands and rounds once through the shared nearest-even IEEE
owner, independently of hardware rounding/flush modes. A nonfinite final
rounding raises `ArgumentError`; inputs remain unchanged on every outcome.
A zero path returns zero voltage and does not select a stationary-node force.

Charge/impulse and voltage share the canonical CIC path stencil. Path voltage
scans N field entries for validation, then uses O(crossed cells) rational
operations and scalar temporaries, without a new workspace. A native fast path,
allocation and latency acceptance remain open.
This work contraction does not resolve the retained first-order CIC point-force
or discontinuous node-at-rest failures. A coupled force/pusher needs its own
shape, time-centering, stability and wall contracts.

The charge/impulse workspace payload is `(8+sizeof(T))*(N+1)` bytes plus object headers. Charge costs
O(N+P); impulse costs O(N+P+crossed cells), at worst O(N+P*N). Exact replay owns
O(N) temporary rational entries. No universal cold allocation or latency bound
is claimed.

The private native/exact stencil signatures express the admitted IEEE scalar
domain without converting positions or restricting public containers. Earlier
isolated physical/resource checks and the failed broad-input compiler diagnostics
remain retained in the external audit archive. New essential regressions and
concrete compiler checks run through GitHub CI; changed-source acceptance is
pending until those jobs pass.

## Optional physical-width triangular cloud

```julia
cloud = YeeTriangularCloud(0.05) # physical half-width, fixed under refinement
yee_particle_charge!(charge, particles, grid, workspace; shape=cloud)
yee_longitudinal_impulse!(impulse, particles, xold, xnew, grid, workspace; shape=cloud)
voltage = yee_longitudinal_voltage(Ex, xold[1], xnew[1], grid; shape=cloud)
stationary = yee_longitudinal_field(Ex, xold[1], grid; shape=cloud)
```

The scalar calls take one coordinate, whereas the impulse takes one endpoint
vector per supplied population. This example's scalar coordinates refer to a
single `ParticleSet`. The shape retains a finite positive Float32/Float64 width
at its original precision even if the grid has another precision. It is a
standalone Yee shape, with no periodic deposit/gather or GPU implementation.
The stationary helper requires explicit cloud selection and returns a field,
without q/m, acceleration or a trajectory.

The density is `S_a(t)=max(a-|t|,0)/a²`. Convolve the closed point/CIC hats with
this density, extending the left/right endpoint hat constantly outside its
nearest wall. Outside cloud tails remain **particle-bound endpoint support**;
they are not absorbed particles, free plate charge or an injection ledger.
All centers in `[0,L]`, coincident/overlapping clouds and every mesh node are
admitted. This static convention does not define a finite-radius mechanical
wall, reflection or absorption.

For exact `h=L/N`, let `P_a` be the density's cumulative and `R_a'=P_a`. Each edge
has one canonical cumulative and derivative:

```
A_e(x) = (R_a(x-e*h)-R_a(x-(e+1)*h))/h
g_e(x) = A'_e(x) = integral_cell_e S_a(r-x) dr/h
C(x) = [1-A_0, A_0-A_1, ..., A_last]
I_e = sum(q*w*(A_e(xnew)-A_e(xold)))
V = dx*sum(Ex*(A(xnew)-A(xold)))
G_E(x) = dx*sum(Ex*g(x))
```

`C` is nonnegative and sums to one. Exact continuity and work contraction follow
from these same cumulative differences. `g` is C1 through mesh nodes, so the
stationary field has one continuous value there. At a zero-field rest equilibrium,
the linear coupling matrix is `K_ef=dx*sum(q²*w*g_e*g_f/m)`, symmetric positive
semidefinite for positive masses and nonnegative weights. This is a spatial
coupling property, not nonlinear/time-step stability or an accepted mover.

Interior DC response retains `eta=N*dx/L`. Near a wall it is
`eta*Ex*integral_0^L S_a(r-x)dr`; at either wall with `a<=L` it is half interior
DC. The continuum `Ex=r` cloud force integral at the left wall is `a/6` when
`a<=L`; for `a>L` it is `L²/(2a)-L³/(3a²)`.
The discrete helper uses the actual edge samples and its declared work metric;
it does not return the unchanged point-wall target zero. Holding a fixed is
essential: the source-derived cell-average force bound is O(h²) with a constant
depending on a. Taking `a=O(h)` does not establish fixed-macrocharge point-force
convergence. Point-sheet/contact and point-wall rejection evidence stays open.

The initial cloud tier evaluates complete stored-input rational polynomials and
whole-entry sums, then rounds once using the shared nearest-even IEEE owner.
It includes original q and weight cancellation and is independent of hardware
rounding/FTZ. Invalid input preserves output and scratch; unrepresentable required
output rejects before publication, with reusable scratch. No point native error
bound is reused. Exact rational support bounds omit edges with identically
zero contributions. With k_p retained edges for particle p, charge/impulse use
O(N+sum(k_p)) rational operations, with worst case O(N+P*N). Fixed physical
width can still span an increasing number of cells under refinement. Scalar
field/voltage validate all N entries and contract only the retained support.
The workspace payload is unchanged. Native arithmetic, allocation and latency
acceptance remain open; no measured speedup is claimed.

The compact new regressions use independent polynomial volume/path quadrature,
whole original-source cancellation, mixed width/grid precision, broad/view
storage, continuity, node/overlap/wall/DC response, stationary Gram pairing,
fixed-width cell-average force bounds and rejected publication/context controls.
They have not run locally. GitHub runtime, formatting and full compiler review
must pass on the integrated revision before this prerequisite is accepted.

## Implementation tracker

The optional static-cloud API, surface-tail model and canonical owners above
were recorded before production edits and reviewed on 2026-11-01. The source
implementation preserves point defaults; exact-revision runtime/compiler/format
acceptance remains a separate gate.

- [x] Closed extensive charge and compatible longitudinal endpoint impulse.
- [x] Canonical longitudinal path voltage and exact discrete work contraction.
- [x] Optional explicitly named finite physical-width triangular-cloud static
  charge/current/voltage and stationary work-adjoint field source implementation.
- [ ] Accept the optional cloud's independent essential regressions, full compiler
  review and formatting on the exact integrated GitHub revision.
- [ ] GitHub validation of the path-voltage implementation and unchanged impulse.
- [ ] Choose a force/current pairing that passes retained point-sheet force,
  isolated-plate/DC, work and semidiscrete stability controls.
- [ ] Transverse current and a synchronized coupled particle/field method.
- [ ] Event-local reflection with crossing, grazing/contact and repeated-impact
  rules, atomic coupled failure and wall momentum/charge/work ledgers.
- [ ] Coupled energy/Gauss diagnostics, background and free plate-charge state.
- [ ] Versioned checkpoint/restart of fields, species, phases and wall/history state.
- [ ] Absorption/injection/circuits, relativistic coupling, higher dimensions,
  GPU/MPI and broader physical validation.

The existing reflecting utility and prescribed-source field step do not supply
these missing pieces. Retained force and impact-convergence failures remain open.
