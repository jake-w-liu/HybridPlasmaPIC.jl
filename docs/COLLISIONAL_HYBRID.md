# Fixed-step collisional hybrid runs

`CollisionalHybridStepper` composes the existing hybrid field step with local
self/interspecies Coulomb collisions ([`collide_coulomb!`](LOCAL_COULOMB.md)). This
slice is **first order** and CPU-only. It supports Float32/Float64 in one to three
spatial dimensions, periodic Fourier grids, positive-ion populations, standard
Vector particle storage, NGP/CIC/TSC, and ideal massless isothermal electrons
(`eta=etaH=de2=0`). The ordered populations, counts, charges, masses, model,
normalization and positive timestep are fixed for the wrapper's lifetime.

## Step order and time levels

Particle positions and fields are at integer time; stored velocity is at
`time-last_dt/2`, and the first loaded velocity is the physical initial
velocity (the `HybridStepper` convention). One nonzero step of size `h` from
`t` to `t+h`:

1. The canonical `HybridStepper` push: Boris kick to the half level `t+h/2`
   (priming on the first step), drift to `t+h`, magnetic field to `t+h` with
   frozen midpoint moments.
2. `collide_coulomb!` for the interval `h`, at the new positions and the
   half-level velocities, with the wrapper's RNG and workspace. The collision
   interval is the field interval, and its velocities sit at that interval's
   midpoint.
3. Integer-time moments, Ohm's law and the carried-E recentering from the
   collided velocities, giving E at `t+h`.

The result is bitwise identical to `step!(inner, ...)`, then
`collide_coulomb!(..., h; ...)`, then recomputing the moments, Ohm's law and the
carried E; the composition forms E once instead of twice. This is a staggered
first-order (Lie-type) split, not a symmetric split or a synchronized exact
collision subflow. Skipping step 3 (reusing the pre-collision E) is wrong: the
next kick then lacks `h (q/m) dE`, where `dE` is the change of `-u_i×B` from
the collisional change of the charge-weighted bulk velocity.

## Timestep, RNG and failures

Zero `dt` is a true no-op, including RNG and priming. Any other `dt` must equal
the configured step at grid precision; a different value is rejected before any
state changes. Variable steps are not supported in this slice.

Supply an explicit `Xoshiro` or `MersenneTwister`. Construction deep-copies it,
so the caller's RNG is not consumed. The wrapper borrows its inner stepper and
populations exclusively: advance them through the wrapper while the composed
run is active. Independent wrappers need independent inner states and particle
storage. `init!` resets the physical clock and preserves the current owned RNG.

A failed step restores particle x/v, every public hybrid field, the clock, the
source ledger and the RNG in place; no partially collided velocities are
published and a retry draws the same numbers. Private scratch (inner workspaces,
collision bins) may change. The standalone `collide_coulomb!` keeps its weaker
policy: it publishes no velocities on failure but may advance its RNG.

Coulomb logarithms may be a scalar or a symmetric channel matrix; zero entries
disable channels. All-zero logs run exactly the plain `HybridStepper` step and
leave the RNG unchanged. `epsilon0` and `weight_to_number` are explicit because
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
incompatible state is rejected before publication.

`archive_run(path, st, ions; rng_seed=42)` uses the same state and existing
checksum/provenance archive. Its seed is metadata, not a replacement for advanced
RNG state. Stochastic restart requires the same Julia version/platform and
compatible environment; cross-version bitwise identity is not promised.

## Verified state (2026-10-26, Julia 1.12.7 and 1.13.2)

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

Variable timesteps, second-order stochastic splitting, evolving/anisotropic
pressure, CAM-CL, electrostatic/electromagnetic PIC, neutral/BGK sources,
ionization, radiation, distributed cell ownership and GPU RNG/storage remain
separate planned extensions.
