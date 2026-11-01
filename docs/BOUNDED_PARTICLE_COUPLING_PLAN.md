# Closed Yee particle charge, current and path work

`YeeChargeWorkspace`, `yee_particle_charge!` and `yee_longitudinal_impulse!`
provide the CPU Float32/Float64 charge/current prerequisite of CAP-002B on a
closed one-dimensional PEC interval. They do not select or advance a particle
trajectory. Coupled bounded electromagnetic PIC remains incomplete.

`yee_longitudinal_voltage(Ex,xold,xnew,g)` supplies the read-only scalar path
work prerequisite on the same closed interval. It does not choose a force or
advance particles.

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
The native accumulation has bound `gamma*S`, where `gamma=(P+10)*u/(1-(P+10)*u)`,
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

## Remaining implementation

- [x] Closed extensive charge and compatible longitudinal endpoint impulse.
- [x] Canonical longitudinal path voltage and exact discrete work contraction.
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
