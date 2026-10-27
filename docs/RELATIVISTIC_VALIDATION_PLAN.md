# Relativistic velocity validation

Status: shared boundary repair verified on Julia 1.12.7 and 1.13.1 and pushed
as `43881de`.

## Confirmed defect and ownership

VERIFIED: exact rational sums of the stored velocity components exceed `c²`
for Float16/32/64 examples that the divided/squared/summed guard accepts.
The shared `_require_subluminal_velocities` in `src/particles/boris.jl` serves
radiation, both relativistic electromagnetic PIC APIs, and electromagnetic
checkpoint validation for electrons and mobile ions.

The public baseline probe passes 28 witness assertions on each Julia version:
radiation accepts the three precisions in 1D–3D, and both PIC APIs accept the
Float32/64 examples. Logs:
`/tmp/hybrid-relativistic-public-baseline112.log` and
`/tmp/hybrid-relativistic-public-baseline113.log`.
The first probe stopped at a fixture error: default particle charge was +1,
while electromagnetic electrons require -1. The corrected probe sets -1 and
completed on both versions; its eight PIC results are explicitly `:accepted`.

## Repair and compatibility

Retain the finite-component and rounded-norm preflight, then compare the exact
stored squared components with exact `c²` when the accepted rounded norm lies
within eight machine epsilons of one. This repair rejects a truly invalid state
before the first particle or field is changed. It changes no public signature,
state representation, or checkpoint schema.

For hardware floats, the normal arithmetic has at most five accumulated
relative-rounding factors: the quotient contributes twice when squared, the
multiplication once, and the two positive additions once each. Bounding each
by one machine epsilon gives a lower bound of `(1-eps)^5 >= 1-5eps` for the
retained contribution. The eight-epsilon guard also leaves room for tiny
underflowed contributions; even Float16's minimum normal squared is less than
one epsilon by over five orders of magnitude. Six exact arithmetic checks of
the normal-error and conservative tiny-contribution bounds pass for the three
hardware types (`/tmp/hybrid-relativistic-roundoff-bound112.log`). This is a boundary predicate,
not a replacement numerical integrator.

Inputs whose rounded norm is already at least one remain rejected. Admitting
additional exactly subluminal states requires a separate repair of downstream
Lorentz-factor arithmetic, which still uses rounded squared norms. This slice
does not claim to solve those downstream range limitations or establish the separate
[radiation arithmetic contract](RADIATION_REACTION_PLAN.md).

## Acceptance checks

- Exact-rational boundary witnesses, both sides of the boundary, all axis/sign
  permutations, three power-of-two light-speed scales and 1D–3D particles.
- Seeded nearby random states classified by an independent exact norm.
- Radiation rejection before cooling the earlier valid particle, including
  zero-force calls; PIC rejection for electron and mobile-ion inputs.
- Checkpoint rejection preserving prior file bytes and destination state;
  continued acceptance of a valid checkpoint.
- Both native Julia versions, neighboring forcing/particle tests, concrete
  compiler calls, package analysis against the unchanged gate, ordinary and
  exceptional allocation measurements, and final formatter/diff review.

Wrong patches distinguished by these checks: rejecting every near-boundary
state loses valid neighbors; checking the rounded `v/c` tuple exactly still
loses division bits at nonunit scales; checking only the first particle or
checking after mutation breaks public rollback assertions. Non-power-of-two
scale witnesses exercise division rounding independently and reject that second
wrong patch on all three hardware precisions.

## Observed verification

The frozen v1 source passes 12,706 focused assertions on each Julia version:
12,486 exact-boundary/permutation/random checks, 36 radiation rollback checks,
and 184 PIC/checkpoint transaction checks. The final v2 adds only two testsets:
15 non-power-of-two division checks and 80 BigFloat stored-precision/rounding
checks, all passing on each version. Final source and retained test text match
v1 exactly; a byte comparison confirms the two added testsets are the entire
v2 test delta. The complete new regression file therefore has 12,801 verified
assertions. Logs:
`/tmp/hybrid-relativistic-v1-focused112.log`,
`/tmp/hybrid-relativistic-v1-focused113.log`,
`/tmp/hybrid-relativistic-v2-added-regressions112.log`,
`/tmp/hybrid-relativistic-v2-added-regressions113.log`.
The last two logs also contain all 1,810 neighboring forcing, particle and
long-orbit checks; their total is 1,905 passing assertions each.

All 18 concrete compiler calls per version have zero reports (36 analysis and
runtime assertions). Package JET retains 26/25 reports against the unchanged
ceiling of 30, matching the baseline; see the existing
[MPI report classification](MPI_CHECKPOINT_PLAN.md). Logs:
`/tmp/hybrid-relativistic-v1-jet112.log` and
`/tmp/hybrid-relativistic-v1-jet113.log`.

All 27 warmed ordinary Float16/32/64, 1D–3D, 0/1/256-particle guard cases remain
unchanged and allocate 0 B on Julia 1.12.7. The three valid boundary-neighbor
cases allocate 3,168/3,168/3,248 B respectively, compared with the old guard's
0 B; exact boundary discrimination has an explicit resource cost. Each resource
run passes 57 assertions. Logs:
`/tmp/hybrid-relativistic-baseline-resources112.log` and
`/tmp/hybrid-relativistic-v1-resources112.log`.
No optimization or whole-project completion claim follows from these checks.

## Output publication and public rollback (C1)

The isolated, unaccepted [C1 implementation](RELATIVISTIC_OUTPUT_PUBLICATION_PLAN.md)
factors this unchanged input predicate for scalar output admission. In that
candidate, Boris and priming retain a supported native result or use their
existing accurate cold stage; an unsupported rounded result raises
`ArgumentError`. Its initialization validates both required relativistic species
before changing gathers, deposition or clocks, and its public relativistic
EMPIC steps use a reusable undo workspace. Current-main source integration and
acceptance remain pending. These are a separate compatibility and storage change
from `43881de` above; the linked tracker owns C1 verification.

## Planned domain and state extensions

These are unimplemented follow-ons. C1 must first pass its isolated and
current-main acceptance gates. Root owns sequencing and integration; the
kinematics track owns the numerical contract. Neither follow-on may reinterpret
a velocity already stored by the existing API as hidden authoritative momentum.

| Stage | Canonical change | Required acceptance |
| --- | --- | --- |
| C2.1 exact stored domain | Inventory every gamma, velocity validator, diagnostic, prime, push, radiation and persistence caller; specify exact stored `sum(vᵢ²) < c²` admission, including vectors whose rounded beta² is at least one | Independent exact classification over all retained C1 inputs, both boundary sides, signed zeros, subnormals, nonunit c and mixed stored BigFloat precisions/modes; no caller left on a conflicting predicate |
| C2.2 conditioned arithmetic | Extend the existing original-component Lorentz-gap owner and exponent-safe stage arithmetic before changing public admission | Independent stage references at successively increased precision with stable target rounding; few-ulp normwise bounds fixed before measurement; supported C1 histories retain their original bits; explicitly reject an unrepresentable required result without mutation |
| C2.3 complete integration | Apply the agreed domain to initialization, every required species, radiation, energy/momentum diagnostics, checkpoint/archive validation and restart | First and established public steps, changed dt/subcycles, mobile/ignored ions, failure rollback, task precision/rounding restoration and restart continuation; no clipping or beta² tolerance |
| C3.1 state decision | Choose one explicit authoritative state and public mutation contract after comparing a dedicated momentum species type with a versioned extension of existing storage | Document which API owns x, momentum, charge/mass and any derived velocity; no two writable authoritative representations or silently stale cache; determine compatibility and migration before production edits |
| C3.2 momentum evolution | Integrate prime, Boris/field coupling, radiation and drift with that authority; derive velocity for consumers without reconstructing gamma from a rounded near-light velocity | Electric acceleration and magnetic rotation with independent momentum/energy references, radiation balance, extreme finite scale/cancellation and timestep convergence; define supported momentum/exponent bounds and distinguish physical growth from arithmetic range failure |
| C3.3 persistence and backends | Version the authoritative payload and staggered phase; route diagnostics, checkpoints, archives and accepted backend kernels through it | Fresh and changed-population restart, malformed/old payload migration or explicit rejection, corruption atomicity, actual CPU/device/MPI intersection gates; no resident/distributed capability claim from host-only coverage |

C2 expands which stored vectors can be evaluated. It cannot require every
physical interior velocity to have an admissible independently rounded tuple.
An axial tuple rounded to exactly c still requires rejection under a
velocity-authoritative contract. C3 addresses that representation limitation;
its domain and migration are decisions to settle before implementation.

Initial source inventory for C2 in the isolated C1 candidate (the full
cross-backend caller graph is still pending):

| Owner | C1 candidate contract to preserve or migrate |
| --- | --- |
| `kinematics_range.jl`: `_subluminal_velocity_supported`, `_conditioned_lorentz_factor`, `_exact_lorentz_gap` | Shared narrow admission and conditioned original-component gap; widening admission before extending its consumers is insufficient |
| `boris.jl`: `_require_subluminal_velocities` | Public particle/species validation delegates to the shared domain |
| `full_pic.jl`: `_boris_rel_parts`, `_boris_rel_native`, both prime/init/step owners | Native and repaired stages, batch bypass, publication and complete public-step rollback |
| `extrapolated_leapfrog.jl`: `_prime_leapfrog_loop!` | General priming uses the shared input/output policy |
| `radiation_reaction.jl`: `_radiation_reaction_generic!` | Generic native gamma currently forms a squared norm directly; IEEE factor arithmetic has its own exact-completion contract and must retain its validated bits |
| `extrapolated_leapfrog.jl`: `kinetic_energy_relativistic`, `_kinetic_energy_relativistic_exact` | Energy already resolves both uncertain sides from exact stored components and explicitly documents a broader domain than the pusher; C2 must not narrow this existing diagnostic |
| `checkpoint.jl`: `_validate_checkpoint_state` for electromagnetic steppers | Both required saved species use the shared validation; preserve validation-before-publication and restart phase |

Freeze resource budgets before either implementation. Measure ordinary,
boundary and rejected complete calls separately on the minimum and latest
Julia. For C3, record persistent bytes per particle, any derived-vector traffic,
per-step allocation, retained storage after population changes and exact
synchronization/transfer ownership. Require matched complete compiler trees,
actual frame fields and source-backed review of every difference; report counts
alone cannot accept either stage. Add latency/scaling acceptance only with a
fixed workload and a valid paired calibration. Push each accepted stage
separately and require all six exact-SHA normal CI jobs.

- [ ] C2 caller inventory and supported-range/rounding contract.
- [ ] C2 original-component numerical implementation and independent oracles.
- [ ] C2 end-to-end, resource, compiler, persistence and CI acceptance.
- [ ] C3 authoritative state/mutation/schema decision and pre-edit budget.
- [ ] C3 canonical evolution and diagnostics with independent trajectory tests.
- [ ] C3 migration/restart/backend integration, resource and exact-SHA CI gates.
