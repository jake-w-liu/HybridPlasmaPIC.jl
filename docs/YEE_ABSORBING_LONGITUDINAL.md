# One isolated passive capture source candidate

`YeeAbsorbingLongitudinalPIC` is a distinct nonrelativistic CPU Float32/64
longitudinal source candidate. It admits zero captures or one uniquely isolated
capture per external step. Numerical/convergence, rollback/restart, compiler,
formatting and complete-step resource acceptance are deferred during the current
implementation-first pass. Existing authored controls remain unrun; new test
execution is paused. It does not complete [bounded PIC](BOUNDED_PARTICLE_COUPLING_PLAN.md).

```julia
state = YeeAbsorbingLongitudinalPIC(grid;
    shape = YeeTriangularCloud(0.0625), max_root_refinements = 128)
init!(state, particles)
step!(state, particles, dt)
diagnostics = yee_absorbing_longitudinal_diagnostics(state, particles)
```

Use a closed 1D YeeGrid, explicit fixed physical cloud width and one ParticleSet
or ordered one-based vector. All particle leaves must be ordinary Vectors.
Nearest hardware rounding, gradual underflow and stored-zero transverse inputs
are required. The first successful positive converted dt fixes external dt;
default root budget is precision(T), with a positive integer override. The cap
is a finite work budget, not a success tolerance. Public live centers are strict
interior; grazing/terminal contact, unresolved ordering/budget, a second event
or inconclusive coverage/passivity reject the whole candidate.

This declares an ideal **passive instantaneous rigid-cloud-to-isolated-plate**
law. When its center hits, the complete cloud retires into that plate, including
particle-bound endpoint tails. There is no gradual erosion, sheath, circuit,
emission, image force or transverse radiation model. Complete the event-duration
KDK final electric half-kick before removal. The wall receives the incoming
represented m*w*v momentum and m*w*v^2/2 kinetic energy. Zero-weight tracers may
be removed/count but carry no physical charge, momentum or energy.

For canonical cumulative A and Bleft=0/Bright=1, source death has exact impulse
`q*w*(Bside-A(wall))`, beyond the actual drift-to-wall impulse. This current
transfers the full signed q*w to the selected free plate and changes Ex through
the canonical E-I row; no voltage reset or Gauss projection occurs. Wide and
mixed-precision clouds retain their original physical support/type/bits. Both
plate charges are exact rational state: immutable `left_plate_charge` stores
initial configuration; `free_left_plate_charge[]` is current left charge, and
`right_plate_charge[]` is current right charge. Diagnostics use current plates.

The death-source work Wdeath is signed work out of the fields. Received wall
energy is Kcaptured+Wdeath. A negative represented received energy rejects.
Additionally its exact error lower bound must certify nonnegative received
energy for the original reference trajectory; a failed lower bound is an
inconclusive passive-capture certificate. This is an explicit reduced law, not
heat invented to cancel an observed conservation error. Captured kinetic,
electrical work and received energy remain separate exact last/cumulative
ledgers. The wall impulse ledger omits conductor/background support stress;
it does not imply conservation of total plasma+wall momentum.

Internal durations and halves, q/m and gathers remain original rational values
through canonical Boris. A constrained contact endpoint has the reviewed
four-ulp actual snap and separate predictor-time gate. Only one uniquely isolated
first hit is allowed. Other particles' premode prefixes exclude an earlier first
hit; they do not pretend forces remain smooth after capture. Survivor errors
align at the true hit, retain source/position/velocity uncertainty, apply the
same ideal death map, gather the new force and certify both walls over the
enlarged remainder. J bounds apply within each smooth epoch, never across the
finite force jump. There is no zero-error reset or reflected-wall departure rule.
These are conditional one-step bounds, not an entire-run accuracy/stability claim.

`last_transport_work[]` is the exact signed motion-source work;
`last_wall_work` holds exact per-wall death work. `last_work[]` rounds their
complete total once. `last_energy_defect[]` rounds complete DeltaKlive+Kcaptured
minus transport work. `last_physical_energy_change[]` rounds physical
Delta(Klive+U), and `last_total_energy_defect[]` rounds Delta(Klive+U)+received
wall energy. `last_field_energy_defect[]` retains exact DeltaU+total source work.
Independent current/field rounding contributes a measured row residual; no
diagnostic is forced to zero. Read-only energy/Gauss may overflow; required stored
IEEE histories must be finite before publication. Error bounds describe surviving
trajectories; empty survivors do not erase field/source or wall ledgers.

Stable survivor compaction preserves public particle objects and every component
array identity, q/m, surviving weights, ids and tags. Logical populations remain
ordered, including empty rows. Source/provenance/configuration and identities are
otherwise bound until explicit init/restore; external charge injection/resizing
is not supplied. Initial-minus-live counts/charge must equal cumulative captures.
The committed source binding is independent of shrinkable scratch, so a rejected
capture can retry. All physical fields, particles, metadata, clocks and histories
stage before publication. Custom setter/OS protection changes and concurrent
stepping/mutation are outside the ordinary-write guarantee; one instance cannot
be stepped concurrently. Matched views reject in this new tier; old models keep
their existing storage policies.

Distinct checkpoint tag `:yee_absorbing_longitudinal_pic`, schema1, saves original
config TYPE/bits, verbatim physical fields, live particles, both current free
plates, initial source records, complete capture/certificate/energy histories and
required versioned policies. No cross-model/schema migration is implicit.
The last-step capture record is absent on no-capture steps; otherwise its exact
local elapsed s and completed external step index are saved. It is not an absolute
timestamp recovered from rounded time-minus-dt. Restore validates/detaches all
incoming values, reserves Vector capacity and prepares transport/binding before
stores. Incoming dynamic plates derive from incoming sources, not destination
history. Archives reuse the atomic/checksummed container4.

## Conditional capture certificate

The reference starts from exact interpretations of this external step's stored
positions, velocities and fields. Let `Qj=qj*wj`, `alpha_j=qj/mj`, `h=L/N`,
`delta=stored_dx` and `a_cloud` be the original cloud half-width. Canonical cloud
bounds give `M=min(1/h,1/a_cloud)`, `sum(g)<=1/h` and
`sum(abs(g'))<=2/(a_cloud*h)`. Before and after death, separately, use

```text
E*=norminf(E0)+sum_initial abs(Qj)
Aj=abs(alpha_j)*delta*E*/h; Vj=abs(v0j)+Aj*tau
Jj=abs(alpha_j)*delta*((M/h)*sum_initial abs(Qk)*Vk
                              +2*E*/(a_cloud*h)*Vj).
```

Here initial means the live population at this external step's start. The field
bound survives capture because each original source's cumulative endpoint moves
within [0,1] and ends at Bside. The larger predeath J remains conservative within
each smooth epoch; no J estimate crosses the finite acceleration jump.

The root gates isolate one true incoming hit `t in [a,b]`, with
`0<a<=s<=b<tau` and `Delta=max(s-a,b-s)`. All omitted particles have both-wall
premode clearance through b, excluding an earlier or tied first hit. Measure the
actual predeath numerical map, including snap and kick/current/field rounding:

```text
dxj=xminus,j-(x0j+s*v0j+s^2*a0j/2)
kappa_j=vminus,j-v0j-s*(a0j+a1minus,j)/2
rho_pre=Eminus-E0+Ipre,*
Xj=abs(dxj)+Jj*s^3/6+Vj*Delta  (survivors); Xcaptured=0
WX=sum_survivor abs(Qj)*Xj; eEpre=norminf(rho_pre)+M*WX
Dminus,j=abs(alpha_j)*delta*(eEpre/h+2*E*/(a_cloud*h)*Xj)
Bj=abs(kappa_j)+Delta*abs(a0j+a1minus,j)/2+b*Dminus,j/2+Jj*b^2/4.
```

These errors align to the single true prehit state. If t>=s, the prehit Taylor
tube plus motion to t bounds X. If t<s, Taylor at t uses `t^3<=s^3`; shifting
only the predictor polynomial to s costs at most `Vj*Delta`. Thus neither case
extends physical predeath flow beyond its true hit. The captured centers agree
at the exact wall; `ev=Bcaptured` bounds its incoming velocity error. The
Lipschitz trapezoid error is at most `Jj*t^2/4`, giving B above.

With original represented mass `mu=m*w`, require `Ereceived>=eK+eW`, where

```text
eK=mu*(abs(vminus)*ev+ev^2/2)
eW=delta*sum(abs(Ideath,*))*eEpre
   +abs(delta*sum(Ideath,*(Eminus-Ideath,*/2))-Wdeath).
```

Both physical and numerical death maps use the same ideal current at the same
wall. The measured signed death-row residual adds to rho; removing the captured
source then preserves the surviving X/B errors. Fresh postdeath gathers give
surviving acceleration errors using total rho and surviving WX. The physical
remainder lies in `[tau-b,tau-a]`, while the numerical duration is `tau-s`.
Every survivor must pass numerical and carried both-wall clearance through the
enlarged horizon; terminal position/velocity bounds retain X/B, acceleration
error, time displacement and actual second-chunk residuals. There is no
departing-wall exemption. The canonical owners are the
[capture seed/death/compaction](../src/models/yee_absorbing_longitudinal_pic.jl)
and shared [remainder/terminal maps](../src/models/yee_impact_longitudinal_pic.jl).
The [reflection contract](BOUNDED_EVENT_COVERAGE_CONTRACT.md) has different
fixed-source premises and does not itself certify this capture law. All these
statements are conditional algebra; full numerical, scientific and resource
acceptance remains pending.

Only essential fixed independent death/plate/work/passivity, force-jump refinement,
late compaction rollback/retry, stable provenance/empty and strict detached restart
controls are included in the existing full/longitudinal suites. At `bb7da93`,
[minimum-Julia essentials](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41442148118/job/125266382813)
record 182 passes and 22 failures in these absorption groups on Julia 1.12.7:
18 snapshot comparisons and four JET assertions. The repair retains every
capture-record field in value-based snapshots and selects explicit history
scalars in the canonical validator. Numerical maps, IO and every assertion remain
unchanged; the repair requires exact-revision remote validation. Multiple
captures, general material absorption/emission, transverse walls, relativistic/
device/distributed execution and wider physical/resource acceptance remain
separate work.
