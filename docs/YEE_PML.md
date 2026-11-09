# Directional ADE PML for serial Yee fields

The construction/advance smoke test passes for Float32/64 and D1–3 in group 10
of the exact-main [CI run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41813949753).
It advances a seeded pulse for twelve steps and checks finite fields and
diagnostics. Reflection/refinement, long-time stability, rollback/restart and
complete-step resource acceptance remain unverified. This is one field-only part of CAP-002C, not completion of bounded PIC.

`YeePML` implements source-free classical directional PML with causal auxiliary
differential equations (kappa=1, alpha=0). It supports serial CPU Float32/64 in
D1–3 and all three E/B components. It differs from the
[scalar-conductivity absorber](YEE_CONDUCTIVE_LAYER.md): each directional curl term
has its own memory. There is no general CPML tuning or material/sheath model.

The continuous family follows the directional coordinate stretch in
[Johnson §3.2 Eq7 and §3.4](https://arxiv.org/pdf/2108.05348).
[WarpX's PML equations and matching conditions](https://warpx.readthedocs.io/en/latest/theory/boundary_conditions.html)
illustrate auxiliary evolution and distinguish continuous matching from finite
discretization. The time map and storage/diagnostics below are this model's
explicit choices; no package-equivalence or discrete reflection result is inferred.

```julia
vacuum = MaxwellVacuum(g; c)
# Seed synchronized fields and explicitly supply valid homogeneous PEC values.
pml = YeePML(vacuum; dt=h, cells=((nlow, nhigh),),
    peak_rate=((rate_low, rate_high),), profile=:quadratic)
init!(pml)
step_maxwell!(pml, h)
d = pml_diagnostics(pml)
save_checkpoint(path, pml)
archive_run(archive_path, pml; rng_seed)
```

The pair tuples illustrate D1; supply one pair per axis in D2/3. The supplied
vacuum must be unadvanced and is exclusively borrowed, with one physical clock.
Do not mutate or advance it outside the wrapper. `init!` resets all memories,
clock/count and budget; it leaves the seeded fields unchanged. This is a declared
zero-prehistory initial value problem. Nonzero field seeds inside the PML do not
imply a reconstructed prior outgoing or steady-wave memory. Reinit resets history;
continuation uses checkpoint restore. Gauss-compatible physical source-free seeds
remain the caller's responsibility; no tolerance or projection is added.

## Directional profiles, state and boundaries

Use nonnegative lower/upper integer thicknesses retaining at least one interior
cell per axis. Enabled faces require finite positive rates after grid-T conversion;
disabled faces require resolved zero rates. Periodic axes must have both faces
disabled. Profiles are quadratic or cubic, with exact integer/half-index geometry
and canonical E/B node/center sampling. Rates and their stored T/type/bits define
the configuration; no wider-keyword original-value guarantee follows conversion.

For component c let j=c+1 and k=c+2 modulo3. The two curl slots differentiate
component k along j with sign+1, and component j along k with sign-1. Each memory
lives on its target E[c] or B[c] array. A suppressed/disabled coordinate has no
memory; an enabled coordinate has two E and two B memory arrays plus corresponding
private candidates. D1/2/3 maxima are4/8/12 live arrays. At a corner, each axis
retains its own rate/memory; rates are not summed into a scalar damping factor.
Derived axis tables require O(sum(n)) rational entries, without full rational fields.

Actual exterior boundaries remain homogeneous PEC/periodic. Tangential E, normal
B and each memory on that constrained target plane must be stored zero. The
canonical absent-wall/difference rules preserve these constraints without a clamp,
including face intersections. Memories must also be zero where their own rate is
zero. At initialized count0 every live memory must be globally stored zero,
including positive-rate unconstrained samples. Both signs of stored zero are valid.

D1 Ex/Bx have no derivative or memory and copy their original values, including
signed zeros. Both transverse polarizations evolve. Ordinary divergence is not
generally conserved inside the artificial layer; no magnetic/charge projection,
particle source/current law, injection or particle exit is provided. A prescribed
current argument rejects in every dimension; vacuum's existing source API remains
available separately. Forcing a PML requires a declared transformed-source history.

## Equations, midpoint map and admission

For each active unsigned directional derivative X=D_d F retain memory psi:

    psi_dot = -lambda_d*(psi + X)
    X_pml = X + psi.

E_dot is c² times the signed corrected curl of B; B_dot is minus the signed
corrected curl of E. Lambda has inverse-time units, psi_E has B/length units and
psi_B has E/length units. Zero-history Laplace elimination gives
X_pml=p*X/(p+lambda_d), the directional stretch. It is not scalar field damping.

A magnetic subflow freezes E and advances B/psi_B. An electric subflow freezes B
and advances E/psi_E. For duration delta, beta=-1 for magnetic and beta=c² for
electric. With r_d=delta*lambda_d/2, the joint midpoint expressions are

    psi1* = ((1-r_d)*psi0 - 2*r_d*X)/(1+r_d)
    F1* = F0 + beta*delta*sum(sign*(X+psi0)/(1+r_d)).

The full map is magnetic Q(h)/2, electric Q(h), magnetic Q(h)/2. Original stored
stage operands/dx and the complete rational expressions round only their required
final field/memory values. No rounded half-duration, derivative, factor, memory
endpoint or increment substitutes for those operands. Each retained stage state
must be finite; a required unrepresentable memory rejects even if a field-only
contraction would fit. Final underflow is ordinary conversion. This is a stage-input
representation contract, not an infinitely precise history or native-bit promise.

Require positive finite converted h equal to configured h, nearest rounding and
gradual underflow, strict canonical Maxwell CFL, and exact
`h*maximum_axis(maximum_face_rate)/2 < 1`. The latter keeps each midpoint denominator
in[1,2) and homogeneous memory multipliers in(0,1]; it excludes the stiff sign-flip
regime. It is a local conditioning policy, not proof of coupled PML stability.
The symmetric map is the proposed second-order discretization in ideal arithmetic;
no measured order is claimed. Variable h is unsupported. With all profiles zero,
there are no memory arrays and the unchanged canonical vacuum candidate acts
directly on the original fields, without ADE arithmetic.

## Signed exchange and publication

Use canonical H=prod(stored_dx)*sum(E²+c²B²)/2. In the artificial layer it is a
norm of stored coordinate fields, not a material energy. With exact averages of
actual accepted stored endpoints, record magnetic-stage exchange
P_A=c²*V*delta*sum(sign*Bbar*psibar_B) and electric-stage exchange
P_C=-c²*V*delta*sum(sign*Ebar*psibar_E). Sum both magnetic halves and the electric
stage in rational arithmetic. Exchange P can be signed; it is not positive heat.

The signed field defect is R=Hnew-Hold+P. The exact budget retains
Hnew-Hold+P=R and Hcurrent+totalP-initialH=totalR. It does not assert R=0 or a
positive augmented PML energy. The reported memory square sums have derivative
units and do not supply such an energy. Canonical energy/Poynting diagnostics
retain stored-field meanings; PEC outer flux remains zero and paired periodic
flux cancels. These outer fluxes are not a measured radiation-loss diagnostic.

Physical fields, live memories, clock and exact ledgers stage before fixed-size
publication. All inputs/configuration/aliases/current budgets validate first;
failed private candidates are overwritten on reuse. Arrays and clock Refs retain
identities. Expected-error rollback assumes ordinary writable arrays/exclusive
access; protected memory, resource exhaustion, custom setters and concurrent
mutation are excluded. No allocation or speed result has been measured.

## Complete restart

Distinct `:yee_ade_pml` schema1 saves the required causal kappa1/alpha0 midpoint-KDK
and staggered-profile policies, fixed configuration, nested synchronized vacuum
schema1, every live memory, ready/count state and exact exchange/defect budgets.
At count0 all memories/last-total ledgers must be zero and initial H=current H;
at count1 last equals total. Advanced memories restore verbatim. Only initialized
ready boundaries save; no waveform, RNG, plan, candidate or cache is invented.

Restore checks matching configuration/storage, validates and detaches all incoming
fields/memories together, checks zero-rate/PEC/count-zero constraints and budgets,
and prepares caches before publication. Old replaceable destination field/memory/
clock/ledger/cache values may be stale. Physical/memory arrays and Refs retain
identities; no reconstruction, field projection or silent vacuum/scalar migration
occurs. Atomic trusted Serialization and the existing checksummed archive4 remain
canonical. Runtime restart and broader source/material/backend/particle/PML
acceptance remain open beyond the construction/advance smoke test.
