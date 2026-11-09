# Fixed-step collisional hybrid runs

`CollisionalHybridStepper` composes the existing hybrid field step with local
self/interspecies Coulomb collisions ([`collide_coulomb!`](LOCAL_COULOMB.md)).
Optional [electron–ion thermal moment exchange](ELECTRON_ION_THERMAL_EXCHANGE.md)
adds a reduced two-way scalar-pressure/ion-thermal source; its new numerical,
compiler/format, scientific and resource acceptance remains pending. The default
`electron_ion=nothing` retains the existing ion-only schema-1 behavior. This
slice is **first order** and CPU-only. It supports Float32/Float64 in one to three
spatial dimensions, periodic Fourier grids, positive-ion populations, standard
Vector particle storage, NGP/CIC/TSC, ideal massless isothermal or evolving
scalar electrons (`eta=etaH=de2=0`), and `TransportElectrons` or
`ConductingElectrons`/`FieldAlignedConductingElectrons` with constant
`eta>=0`, `etaH>=0`, `de2=0`. The ordered populations, counts, charges,
masses, model, normalization and positive timestep are fixed for the wrapper's
lifetime. The bound model is detached at construction; comparisons retain the
concrete model/closure type and original scalar values, including BigFloat
coefficients, rather than scalar backing-storage identity.

## Step order and time levels

Particle positions and fields are at integer time; stored velocity is at
`time-last_dt/2`, and the first loaded velocity is the physical initial
velocity (the `HybridStepper` convention). One nonzero step of size `h` from
`t` to `t+h`:

1. The canonical `HybridStepper` push: Boris kick to the half level `t+h/2`
   (priming on the first step), drift to `t+h`, magnetic field to `t+h` with
   frozen midpoint moments. With stored-pressure closures, pressure advances
   with B through the canonical RK4 stages; transport adds its existing Ohmic
   and applicable hyperresistive heat there.
2. `collide_coulomb!` for the interval `h`, at the new positions and the
   half-level velocities, with the wrapper's RNG and workspace. The collision
   interval is the field interval, and its velocities sit at that interval's
   midpoint. This ion-ion subflow leaves stored electron pressure fixed.
3. With `electron_ion` enabled, the NGP/grid-gamma-5/3 thermal source exchanges
   actual ion kinetic energy with stored electron pressure. The linked thermal
   guide owns its coefficients, rounding/entropy gates and limitations.
4. Integer-time moments, Ohm's law and the carried-E recentering from the
   collided velocities, giving E at `t+h`. Transport refreshes the independently
   carried ion force as well as the full field; its predictor uses that force.

With thermal exchange disabled, the result is bitwise identical to `step!(inner, ...)`, then
`collide_coulomb!(..., h; ...)`, then recomputing the moments and the closure's
carried field/force before carried-E recentering; the composition forms E once
instead of twice. This is a staggered
first-order (Lie-type) split, not a symmetric split or a synchronized exact
collision subflow. Skipping the final refresh (reusing pre-source E) is wrong: the
next kick then lacks `h (q/m) dE`, where `dE` is the change of `-u_i×B` from
the collisional change of the charge-weighted bulk velocity.

For `EvolvingElectrons`, `TransportElectrons` or `ConductingElectrons`, set
finite nonnegative `inner.fields.pe` before `init!`.
Initialization preserves it, and collision feedback enters the next pressure
interval through the refreshed moments. The default ion-only subflow adds no
electron heat; the opt-in thermal source transfers energy between the existing
pe and ion reservoirs. Pressure
retains the canonical explicit-stage positivity checks and is never clipped.

Transport retains `E = E_force + eta*J - etaH*laplacian(J)`, force-based ion
kicks and the canonical stage-based positive Ohmic and hyperresistive pressure
heat. S2's heat includes its stated Nyquist convention; see the
[model and discrete balance](TRANSPORT_HYPERRESISTIVITY_PLAN.md).
The ion-ion substep adds no electron heat source.
The bounded [conduction model](EVOLVING_ELECTRONS.md#constant-isotropic-conduction-source-candidate)
adds its original-quotient heat flux only in those field stages, using bounded
two-term preparation, full stored-k/Nyquist and exact per-axis line-mean
projection. Its FFTW accuracy remains qualified, and finite stored-rate,
strict pure-conduction and pressure admission still apply. κ is fixed as part
of the wrapper's bound model. Derived quotient scratch must be independent of
the outer force snapshot.

The wrapper combines the existing reduced transport model with ion-ion Coulomb scattering;
it does not reinterpret the transport terms as electron-ion Coulomb equilibration.
`eta`, `etaH` and the Coulomb logarithms remain independent model inputs; the wrapper
does not derive one from another.
Deposited density must remain finite and at least `nfloor` at initialization
and every midpoint, endpoint and predictor, as in the canonical transport model.
`eta=etaH=0` matches evolving electrons only on states admitted by that density
rule. S1's `etaH=0` arithmetic is retained. Inertial and anisotropic composition
remain unsupported. S2 inherits its required finite physical current Laplacian,
electric field and complete rates; see the [source-range contract](TRANSPORT_SOURCE_RANGE_PLAN.md).

## Timestep, RNG and failures

Zero `dt` is a true no-op, including RNG and priming. Any other `dt` must equal
the configured step at grid precision; a different value is rejected before any
state changes. Variable steps are not supported in this slice.

Supply an explicit `Xoshiro` or `MersenneTwister`. Construction deep-copies it,
so the caller's RNG is not consumed. The wrapper borrows its inner stepper and
populations exclusively: advance them through the wrapper while the composed
run is active. Independent wrappers need independent inner states and particle
storage. `init!` resets the physical clock and preserves the current owned RNG.

A failed step restores particle x/v, every public hybrid field including evolving
pressure, the transport carried force, clock, source ledger and RNG in place; no partially
collided velocities are published and a retry draws the same numbers. Private
scratch (inner workspaces, collision bins) may change. The standalone `collide_coulomb!`
keeps its weaker policy: it publishes no velocities on failure but may advance
its RNG.

The wrapper owns an independent three-grid force snapshot for transport, in
addition to the inner stepper's force transaction. Isothermal and ideal evolving closures
store `nothing` for this snapshot and allocate no added force arrays. New
complete-call allocation, memory and timing acceptance remains open.

Coulomb logarithms may be a scalar or a symmetric channel matrix; zero entries
disable channels. When thermal exchange is disabled, all-zero ion-ion logs run
exactly the plain `HybridStepper` step and leave the RNG unchanged. With thermal
exchange enabled, both sets of logarithms must be zero for that physical reduction. `epsilon0` and `weight_to_number` are explicit because
macro-weights are not automatically real-particle counts. Use
`coulomb_normalization` with your `PlasmaUnits` for a physical unit conversion;
do not infer units from the example values below. `u_floor` is a velocity and
must use the same units as particle velocity.

## Example and restart

The following dimensionless example exercises stepping and exact continuation
within the same Julia runtime and environment:

```julia
using HybridPlasmaPIC, Random, Test

function make_run()
    g = FourierGrid((8,), (2pi,))
    ions = [ParticleSet{1,Float64}(16; q=q, m=m) for (q,m) in ((1.0,1.0),(2.0,4.0))]
    for (a,p) in enumerate(ions), j in eachindex(p.weight)
        p.x[1][j] = (j - 0.5) * g.L[1] / length(p.weight)
        p.weight[j] = g.L[1] / length(p.weight)
        p.v[2][j] = 0.02 * sin(p.x[1][j] + a)
    end
    inner = HybridStepper(g, HybridModel(IsothermalElectrons(0.02)), CIC(), 32)
    fill!(inner.fields.B[1], 0.2)
    st = CollisionalHybridStepper(inner, ions;
        dt=0.001, rng=Xoshiro(42), coulomb_log=0.05,
        epsilon0=1.0, weight_to_number=1.0)
    return st, ions
end

st, ions = make_run()
init!(st, ions)
for _ in 1:10
    step!(st, ions, st.dt)
end
mktempdir() do directory
    checkpoint = joinpath(directory, "collisional-state.bin")
    save_checkpoint(checkpoint, st, ions)
    resumed, resumed_ions = make_run()
    load_checkpoint!(resumed, resumed_ions, checkpoint)
    # Loading restores initialization and the advanced RNG. Continue directly.
    for _ in 1:5
        step!(st, ions, st.dt)
        step!(resumed, resumed_ions, resumed.dt)
    end
    @test all(p.x == q.x && p.v == q.v for (p,q) in zip(ions,resumed_ions))
    @test st.inner.fields.B == resumed.inner.fields.B
    @test st.inner.fields.E == resumed.inner.fields.E
    @test st.inner.time[] == resumed.inner.time[]
    @test st.source_count == resumed.source_count
end
```

Only initialized, completed boundaries can be saved. The versioned checkpoint
stores the complete advanced RNG, all public hybrid fields, canonical particle
state, fixed collision configuration and source ledger. Construct a matching
fresh wrapper and call `load_checkpoint!`; do not initialize it afterward.
The destination must have the same grid, precision, model/shape values,
population order/counts/charges/masses and RNG algorithm. Arrays stay in place;
this slice does not resize particles on restore. Missing/malformed, mid-step or
incompatible state is rejected before publication. Evolving pressure is stored
in the existing schema-1 full-field and inner records; the reader requires those
copies to agree before restoring either state.
Transport uses the inner record's required `transport_policy_version=1` and
`Eforce`; neither is defaulted or reconstructed from the full E. Force arrays
must be finite, correctly shaped and independent of all saved fields and
particles. The force is restored in place, and outer failures restore its
entry value along with pressure and the RNG.
Conducting state also requires integer inner `conduction_policy_version=2`,
including for κ=0; missing, policy-1 or malformed versions reject before mutation.
Field-aligned conduction instead uses required inner policy 3, including its
zero-coefficient case. The outer schema remains 1 with thermal exchange disabled;
enabled thermal exchange requires schema 2 and its complete configuration/report.
The destination must match the saved conduction coefficients as
part of the fixed model, and quotient/transform scratch is not persisted.

`archive_run(path, st, ions; rng_seed=42)` uses the same state and existing
checksum/provenance archive. Its seed is metadata, not a replacement for advanced
RNG state. Stochastic restart requires the same Julia version/platform and
compatible environment; cross-version bitwise identity is not promised.

## Verification status

The ideal `EvolvingElectrons` and S1/S2 transport extensions have source-reviewed
admission, transaction and checkpoint paths, with essential component-sequence,
force/full-E, zero-log, eta-zero, pressure/force/RNG continuation and failure
regressions. An owner-level control restores the outer transaction after a
completed inner step and explicit test exception; it does not reproduce a stock
public post-completion failure.

The [focused run at `213ed8c`](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41278911757)
passes all 2,729 assertions on each actual Julia 1.12.7 and 1.13.3, covering
the three existing collisional test files listed below. This verifies those
component, persistence and production-invariant cases. Full exact-revision
package, formatting and compiler acceptance passed at `b5ce21a`
([run 41683122086](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41683122086),
Julia 1.12.7 and 1.13.3, all 26 jobs green); resource acceptance remains
pending. No new energy-conservation or splitting-order claim is made.

The `ConductingElectrons` composition extends these existing owners. Essential
component/zero-log, temperature/outer-force alias, guard rollback and advanced-RNG
checkpoint/archive controls are in the same files. Their exact-revision GitHub
runtime, compiler and formatter checks passed in the same run; resource
acceptance remains pending.

The following historical results cover **isothermal electrons only**
(2026-10-26, Julia 1.12.7 and 1.13.2):

- The composed step equals the explicit component sequence bitwise over D=1–3,
  Float32/64, both RNGs, scalar and matrix logs, NGP/CIC/TSC; zero logs equal a
  separate unwrapped `HybridStepper`, with interleaved zero steps.
- An uninterrupted run and a run restarted from a file at step 4 agree bitwise
  after 10 steps for both RNGs; resetting only the RNG to its seed diverges.
- Rejected steps (unrepresentable sampled velocities, field failures, invalid
  inputs, changed `dt`) restore the full entry state and RNG; retries repeat the
  same outcome. Injected failures and interrupts after the field advance, inside
  and after the collision, during E formation, after the final validation and
  during restore do the same.
- The next-kick effect of the carried-E recomputation equals `h (q/m) dE` within
  the Boris rotation bound in 48 measured profiles. The original Float32
  detection gate (`> 128 eps(Float32)` at `h = 0.01`) fails 32 assertions on
  both runtimes because that physical signal is 2.1e-6–2.3e-5, independent of
  precision; it is not an implementation defect.
- Save/load/archive restore exact state; 64 repeated cycles retain at most
  1,856 bytes of live heap. Their per-call allocated bytes vary by up to about
  9% between identical states because identity-hash tables (`IdDict` in
  `deepcopy`, Serialization and checkpoint validation) grow according to object
  addresses. The plain `HybridStepper` checkpoint shows the same variation.
- A composed step allocates no more than the plain step plus a separate
  `collide_coulomb!` call (288 small profiles per runtime), and exactly as much
  in 32k–131k-particle runs. Its process CPU time is 1.006–1.082 times their sum (minimum
  over 21 interleaved samples; 32k–131k particles, D=1–3, Float32/64, NB=2) on
  a heavily loaded Apple M5. This is not a latency acceptance campaign.

The tests are `test/test_collisional_hybrid.jl`, `test_collisional_checkpoint.jl`
and `test_collisional_physics.jl`. The [composition plan](COLLISION_COMPOSITION_RESTART_PLAN.md)
records the design history and evidence.

The reduced NGP thermal source has its own linked guide and pending acceptance.
CIC/TSC thermal overlap, consistent drag and electron inertia, variable timesteps,
higher-order splitting, anisotropic composition, CAM-CL, neutral/BGK sources,
ionization, radiation and GPU RNG/storage remain open. Separate bounded
[collisional ES](COLLISIONAL_ELECTROSTATIC.md) and
[EM](COLLISIONAL_EMPIC.md) families and their actual-pencil counterparts have
source implementations with their own deferred physical/runtime acceptance.
