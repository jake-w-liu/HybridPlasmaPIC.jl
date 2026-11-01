# Closed Yee particle charge and current

`YeeChargeWorkspace`, `yee_particle_charge!` and `yee_longitudinal_impulse!`
provide the CPU Float32/Float64 charge/current prerequisite of CAP-002B on a
closed one-dimensional PEC interval. They do not select or advance a particle
trajectory. Coupled bounded electromagnetic PIC remains incomplete.

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

Workspace payload is `(8+sizeof(T))*(N+1)` bytes plus object headers. Charge costs
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
