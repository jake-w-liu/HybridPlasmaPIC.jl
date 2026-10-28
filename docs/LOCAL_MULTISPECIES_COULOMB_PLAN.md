# Spatially local multispecies Coulomb collisions

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
and [local collision guide](LOCAL_COULOMB.md) for current scope.

Status: landed as standalone CPU operator `1b7b8b7` (2026-10-26);
[all six exact-SHA jobs and complete compiler review pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40272777081). See the
[allocation repair record](LOCAL_COULOMB_RATE_ALLOCATION_PLAN.md) for landing
verification. The fixed-step CPU isothermal [CAP-010 composition](COLLISIONAL_HYBRID.md)
has also landed; broader closures/models, MPI ownership and device kernels remain open.
Original base: `bf2a8bc615f514f8d2868934849bf44814140bbb`.
Historical V10 candidate assembly: `c5e8472084dff861ae92d8ff6466edfaaeeebe44`; production baseline is unchanged from `1097076` and its accepted normalization prefix is retained.
Evidence directory: `../HybridPlasmaPIC-audit-evidence/2026-10-23/local-multispecies-coulomb`.

## Verified starting point (original baseline)

`src/particles/collisions.jl:607` exposes one-population, whole-set collisions;
the permutation at line 656 does not use positions. Its pair scatter at line 389
contains the canonical angle sampling and relative-velocity rotation. The
midpoint at line 422 assumes equal physical masses. The public documentation
already identifies cell-local pairing as an upgrade (line 594).

`src/particles/sorting.jl:28` owns periodic cell indexing, including the rounded
upper-box seam. `particle_set.jl:136` owns storage-alias checks. The existing
collision owner stages velocities before publication. `src/io/checkpoint.jl:3`
explicitly assumes RNG-free stepping; `archive_run` records a seed, not live RNG
state. These are implementation and scope facts, not a runtime acceptance claim.

## Pre-edit brief

- Mode: design and implementation, preserving the documented 0-D operation.
- Requirement: local self and cross collisions for distinct physical species,
  with correct pairing rates, conservation semantics, reusable storage and
  observable validation failures.
- Canonical owner: `particles/collisions.jl`; extract its angle/rotation once,
  and use it for the old and new entry points. Extend the periodic indexing
  owner rather than copy its wrapping policy into collisions.
- Compatibility: retain `collide_coulomb!(ps, gcoeff, dt)`, seeded draw ordering,
  ordinary numerical results and allocation behavior. Its prescribed density
  coefficient remains distinct from the new operator's computed cell density.
- Invariants: particles in different cells never pair; positions, weights, IDs,
  tags, charges and masses stay unchanged; no particle self-pair; physical mass
  controls the center-of-mass transform. Equal macro-weights conserve momentum
  and energy to roundoff. Unequal macro-weights conserve them in expectation.
- Failure policy: validate every population before RNG draws or scratch writes.
  Stage all velocities and publish only after successful sampling. Invalid
  sampled states leave every population unchanged; RNG and derived scratch may
  advance. Resource exhaustion and user-provided storage/RNG callback side
  effects are not a stronger transactional guarantee.
- Strategy: periodic CPU operator first; explicit dependent stages below.
- Verification: independent pairing enumeration, analytic velocity-generator
  oracles, statistical and deterministic invariants, runtime/compiler/resource
  gates on Julia 1.12.7 and 1.13.2, and legacy differential controls.

## Public contract for the first implementation

Proposed interfaces:

```julia
workspace = CoulombCollisionWorkspace(species, grid)
collide_coulomb!(species, grid, dt;
    coulomb_log, epsilon0, weight_to_number, u_floor=1e-3,
    rng=Random.default_rng(), workspace=nothing)
```

`species` is a tuple or vector of independent `ParticleSet{D,T}` populations;
`T` is Float32 or Float64 and `D` is 1, 2 or 3. The grid is a periodic
`FourierGrid{D,T}`. The first implementation uses CPU vector storage; device
storage must fail explicitly. Inputs may contain finite periodic images outside
the principal box. Finite nonnegative weights are required; zero-weight
particles are excluded from pairing and remain unchanged.

`coulomb_log` is a required finite nonnegative scalar or symmetric species-by-
species matrix. A zero matrix entry disables that channel. Finite positive
`epsilon0`, `weight_to_number`, masses and `u_floor`, finite charges and
velocities, and finite nonnegative `dt` are required even for an inactive step.
A neutral species does not scatter. `weight_to_number` converts stored macro-
weight into the number of real particles it represents in the chosen unit
system. Both normalization keywords are required on this new API: physical-
number weights in a consistent unit system use `weight_to_number=1` explicitly.
There is no prior local API whose defaults need preserving.

For dimensional reference units `L0,T0,M0,Q0`, use positions in units `L0`,
velocities in `L0/T0`, masses in `M0`, charges in `Q0`, and
`epsilon0 = epsilon0_SI*M0*L0^3/(Q0^2*T0^2)`. If stored weights equal physical
particle counts divided by `N0`, set `weight_to_number=N0`. In particular,
`set_density_weight!` with density normalized to `n_ref` and positions normalized
to `L0` implies `N0=n_ref*L0^3`. For SI inputs whose weights already count real
particles, use `epsilon0_SI` and `weight_to_number=1`. A hybrid calculation
with dimensionless density, charge and mass cannot infer its physical collision
rate from `epsilon0=1`; it must supply these parameters from its physical
reference plasma. Field normalization alone does not determine the number of
physical particles represented in a normalized volume.
The regularizing `u_floor` uses the particle velocity units; changing units
also requires converting this floor. For `PlasmaUnits`, its SI value is
`to_SI(u_floor, :velocity, units)`.

`coulomb_normalization(units::PlasmaUnits)` in the existing
`src/verification/normalization.jl` owner supplies this pair of keywords from
the package's reference plasma. It evaluates the complete expressions from the stored physical references
with the shared bounded-significand normalization owner. For D<3, omitted coordinates have unit
normalized extent. The physical volume of a bin is `prod(g.dx)*L0^3`; its
resolved SI volume is `prod(g.dx)*L0^D` and the reference transverse measure is
`L0^(3-D)`. Thus normalized weights convert by `n_ref*L0^3` in every D. An SI
D-dimensional density weight instead counts particles per omitted transverse
measure, `w_SI=w_normalized*n_ref*L0^D`; its partner density is
`w_SI/prod(dx_SI)`. The complete SI↔normalized generator and seeded collision
updates are checked separately in each dimension.

The workspace owns only derived data: staged velocities, stable cell bins,
shuffled indices and offset/cursor arrays. It owns no RNG or simulation state.
Particle growth may resize its arrays before scattering; repeated calls with
unchanged dimensions reuse storage. No particle reordering is exposed. The
workspace must not be used concurrently; independent workspaces and RNGs can
serve independent simulations. Aliasing between populations, workspace arrays
or live particle arrays is rejected before mutation.

The operation visits cells and enabled unordered species pairs in a defined
order. An identical initial state, species order, explicit RNG state and Julia
environment reproduce the result. Reordering species need not produce identical
bits; exchanging species must preserve the physical distributions within
statistical error. A seed alone cannot restart an already advanced RNG.

## Pairing and rate derivation

The physical prefactor follows the two-body, nonrelativistic Takizuka–Abe
variance with reduced mass `mu = ma*mb/(ma+mb)`:

```
K_ab = (qa*qb)^2 * logLambda / (8*pi*epsilon0^2*mu^2)
variance(delta) = K_ab * effective_partner_density * dt / max(abs(g),u_floor)^3
delta = tan(theta/2)
```

This is also obtained by integrating the small-angle Rutherford deflection
`theta = 2*qa*qb/(4*pi*epsilon0*mu*g^2*b)` over impact annuli
`n*g*dt*2*pi*b*db`, then dividing `E(theta^2)` by four. This independent
derivation fixes the factor of two between angle-pathlength and tangent
variance. See [Takizuka and Abe (1977)](https://doi.org/10.1016/0021-9991(77)90099-7)
and [Higginson, Holod and Link (2020), equations 6 and 29–30](https://www.osti.gov/servlets/purl/1623345).

For one cross pair, use `weight_to_number*max(wi,wj)` in the rate and accept each physical
velocity update with probability `w_other/max(wi,wj)`. This cancels the maximum
in each particle's expected rate. The center-of-mass transform uses physical
masses, never macro-weighted masses. The downloaded HHL manuscript supports
these choices; its discussion also requires timestep convergence when weights
or counts are disparate. The isotropic branch retained from the existing
kernel is a strong-scattering model, not a substitute for convergence tests.

Within a cell, let `Na` and `Nb` count positive-weight particles. Shuffle both
cross lists, cycle the smaller list to cover the larger, and let `D_j` be the
actual reuse count of a smaller-list particle. Each pair uses normalization
`max(Na,Nb)/D_j`. Both effective weights are divided by the same `D_j`, so their
acceptance ratio is unchanged. For a fixed large-list particle, random pairing
selects small particle `j` with probability `D_j/max(Na,Nb)`; multiplication by
the normalization cancels that probability. The small particle's `D_j` visits
likewise cancel its denominator. This argument holds for non-integer count
ratios and weights correlated with velocity, at the first-order collision
generator. Finite-step evolution still requires convergence measurements.

For self collisions, the number of available partners is `N-1`, excluding the
particle itself. Even counts use disjoint random pairs. Odd counts use the
existing three-particle cycle at half variance plus disjoint pairs. For any
unordered pair, its expected graph coefficient is
`3/(N*(N-1)) + (N-3)/(N*(N-1)) = 1/(N-1)`; multiplying by `N-1` gives the
all-pairs generator. Empty and singleton cells are no-ops. Zero-weight entries
do not change these positive-particle counts.

This is an explicitly chosen empirical Landau generator: each represented
particle interacts with the other represented particles, excluding itself.
HHL equation 28 uses `N-1` for even self populations and a different odd-count
duplication scheme. The retained TA triplet needs the enumeration above,
not HHL's odd-duplication coefficient. A continuum finite-N prescription that
uses total density `N*w/V` instead of partner density `(N-1)*w/V` differs by
`N/(N-1)` for equal weights. The new API does not conceal this finite-N modeling
choice. A one-cell legacy comparison uses
`gcoeff=K_ab*weight_to_number*(N-1)*w/V`; it does not equate the two APIs for
arbitrary legacy `gcoeff`.

For a fixed unlike pair with relative vector `g` and mass fraction
`fa = mb/(ma+mb)`, the independent small-step oracle is

```
E[Delta va]/dt = -2*fa*K_ab*weight_to_number*wb*g/(V*abs(g)^3)
Cov(Delta va)/dt = 2*fa^2*K_ab*weight_to_number*wb/(V*abs(g)) * (I - g*g'/abs(g)^2)
```

Sum this expression over every allowed physical partner to test the stochastic
operator, independently of its pairing implementation. Exact pairing-graph
enumeration tests normalization without statistical noise; repeated seeded
collision trials test the sampled distribution and finite-step convergence.

## Small-angle and extreme-range arithmetic

The canonical sampler owns the angle, azimuth and rotation basis once. For a
Gaussian tangent sample `delta`, the local operator evaluates
`cos(theta)-1 = -2delta^2/(1+delta^2)` directly (using inverse-delta arithmetic
when needed). Applying the rotation matrix minus the identity gives the
relative increment without subtracting two nearly equal rotated/original
vectors. A zero-angle increment is exactly zero. The existing 0-D caller
continues consuming the old renormalized relative-vector representation;
its seeded numerical/RNG contract is compared separately against the parent.

The local physical updates are `va + mb/(ma+mb)*Delta_g` and
`vb - ma/(ma+mb)*Delta_g`. Range-sensitive completion combines stored values
as exact rationals before rounding. If the raw relative increment overflows,
the finite original and rotated relative components supply that exact
completion. Only accepted macro-particle proposals must be finite; rejected
proposals do not create a fictitious failure. An accepted unrepresentable
result aborts publication for every species.

Independent high-precision tests cover oblique zero/tiny angles, both
precisions, mass ratios, noncubic cells, extreme prefactors, and an
unrepresentable relative increment whose physical updates remain finite.
Pair and thermal conservation tests also constrain roundoff introduced by
using the direct increment rather than renormalization noise.

## Stages and acceptance gates

- [x] Inspect canonical collision, periodic grid, storage, model and persistence
  owners; preserve the primary paper and current source evidence.
- [x] Root reviewed the design before production edits; explicit unit conversion
  and finite-N generator choices are recorded above.
- [x] Prove pairing expectations by independent enumeration for even/odd self
  populations and unequal cross counts, including 2:3 and 3:7.
- [x] Extract the shared scatter and check legacy outputs, RNG state and warmed
  allocations against the frozen parent. Both runtimes pass 8,064 ordinary
  differential assertions and 672 Float16/BigFloat assertions. The allocation
  measurement controls and retained initial failures are recorded below;
  latency acceptance remains separate.
- [x] Implement reusable cell bins and transactional multispecies collisions;
  runtime and resource acceptance had not completed at that historical V10 checkpoint.
- [x] Test both precisions and all dimensions: locality, periodic seams,
  empty/odd/sparse cells, neutral/zero channels, mass and charge ratios,
  weight/count repartitioning, no self pairs, and untouched metadata.
- [x] Test validation and rollback: later-species invalid data, invalid samples,
  nonfinite parameters, numeric-range boundaries, overlapping storage,
  workspace reuse/growth, and recovery after a rejected call.
- [x] Check fixed-pair momentum/energy, weighted expectation conservation,
  drift/diffusion rates, thermal exchange and timestep refinement. Freeze seeds,
  sample sizes and error budgets before observing candidate outcomes.
- [x] Exercise continuation using restored particle and explicit RNG state,
  recreating scratch. Document that existing checkpoint files do not save RNG.
- [ ] Run neighboring collision, particle, sorting, collection, checkpoint and
  reproducibility tests; concrete JET, full package report comparison, Aqua,
  formatter AST/idempotence checks and exact final diff review.
- [ ] Measure warmed bytes and elapsed time over D1–3, Float32/64, 1/2/4 species,
  dense/sparse cells, equal/odd/unequal counts, zero/neutral channels and weight
  ratios. Reusable ordinary calls target zero allocations after warmup. Retain
  every failed gate; no same-source timing retries to obtain acceptance.
- [ ] For unchanged legacy work, require identical numerical results and bytes,
  and both paired-median and ratio-of-medians latency at most 1.5 times the
  parent under a coordinated exclusive timing window. New functionality has no
  working parent speed comparison; report its measured scaling and memory.
- [ ] Root composes the test runner, reviews, pushes the independently accepted
  stage, and checks all six exact-SHA GitHub CI jobs.

Plausible wrong implementations and their rejection gates: macro-weighted COM
passes conservation but fails physical equilibration; uncorrected duplication
passes equal-count tests but fails enumerated rates; cross-cell pairing can
isotropize globally but fails isolated-cell invariants; per-species publication
passes the happy path but fails later-species rollback; a new duplicated
rotation may pass broad tolerances but fails source ownership and legacy seeded
differential checks.

## Dependent work kept open

1. Collision-aware hybrid/EM integration needs a defined velocity time level,
   field/moment refresh policy and split-order convergence. A complete step
   must retain atomic state semantics across the field and collision substeps.
2. Stochastic restart needs a versioned RNG/configuration/split-phase payload
   owned by the simulation state and written through the existing atomic
   checkpoint publication path. An archive seed is insufficient.
3. Distributed collisions require all colliding particles in a cell on the
   owner rank after migration and collective error/rollback semantics; the
   first CPU operator alone does not establish this.
4. CUDA/Metal require device binning, device RNG ownership, paired updates and
   conservation/rate/resource tests. No CPU scalar-index fallback is implied.
5. **Confirmed legacy numerical follow-up:** the prescribed-rate 0-D overload
   loses tiny original velocities at an exact null angle.
   `legacy-null-v1-1132` fails all eight identity gates in both frozen parent
   and candidate. A separate high-precision witness for neighboring nonzero
   oblique angles is prepared but has not yet run. Migrating this existing
   overload to stable increments changes
   seeded floating-point results and requires a separately reviewed fix with
   its own correctness, compatibility and cost acceptance. This stage preserves
   that legacy behavior; the package-wide defect is not claimed fixed.
6. Relativistic and quantum/degenerate collisions, dynamic Coulomb logarithms,
   nonperiodic/AMR collision volumes, recombination and reaction networks remain
   separate physics implementations. These are not furnished by this stage.

The capability remains open until its accepted stages and dependent scope are
explicitly reconciled in the main completeness tracker. No package-wide or
maximum-performance claim follows from this plan.

## Iteration evidence

All paths below are under the evidence directory. JSON command records retain
the exact runtime executable, command, source hashes before/after, PID and exit
status. No latency gate has run.

- `pairing_oracle.py` / `.json`: exact rational enumeration passed self counts
  2–7 and cross counts 1–4 plus 3:7 and 7:3. Pair expectations equal one;
  additional rational-weight sums match every opposing partner exactly.
- `legacy-v1-{112,1132}.log`: 8,064 assertions per runtime over 576 profiles;
  exact legacy velocity/scratch/RNG results and unchanged measured bytes.
- `focused-v1-1132.log` and `focused-v2-{112,1132}.log` failed because the
  fixture RNG's parametric `randn` method was ambiguous with the stdlib method.
  Their source copies are retained. The fixture now uses an explicit union of
  supported type objects; production sampling was not bypassed.
- `focused-v3-{112,1132}.log`: 1,125 assertions per runtime passed, including
  pair physics, locality/conservation, unequal-weight expectations, rollback,
  explicit RNG continuation and complete dimensional normalization checks.
- `statistical-v1-{112,1132}.log`: 64 independent Landau drift/diffusion gates
  per runtime passed, using 32,768 independent cell samples per profile.
- `resources-v1-1132.log`: all 108 profiles failed the zero-allocation target.
  `allocation-diagnosis-v1-1132.log` attributes 2,304 B in its one-species case
  to mixed-type array tuples in validation and bin resizing. The checks were
  preserved while traversals were divided into homogeneous groups.
- `null-rotation-v1-{112,1132}.log`: eight zero-angle identity checks per runtime
  failed. Midpoint reconstruction discarded a small original velocity. The
  new local update now applies mass-weighted increments of the rotated relative
  velocity and uses exact cold arithmetic on range/cancellation boundaries.
  The existing 0-D reconstruction remains unchanged.
- `null-rotation-v2-1132.log`: all eight checks passed. `focused-v4` passes all
  1,125 assertions on each runtime after the local reconstruction change.
- `resources-v2-1132.log`: 36 profiles passed; D=2 reported 16 B/call and D=3 48 B/call (72 profiles).
  `allocation-diagnosis-v2-1132.log` reports zero allocations inside its concrete
  public call and all three previously allocating owners. The remaining
  initial concrete probe was D1; the later full-matrix attribution below
  identifies the remaining D2/D3 owner.
- `thermal-v1-1132.log`: six independent Maxwellian exchange/refinement gates
  and eight unequal-mass equilibration/conservation checks passed. The finest
  Float64 measured energy-exchange rate was 1.055295 versus the analytic
  1.074871, with measured standard error 0.013227. The minimum-Julia
  counterpart also passes all 14 thermal gates.

Remaining gates are the complete diagnostic call-tree comparison, final
neighboring regression checks, legacy latency and new scaling measurements,
composition onto current main, and exact-SHA CI. Earlier successful runs do not
replace those final-source checks.

- `resources-v3-{112,1132}.log`: a typed measurement barrier retained the same
  failures; changing the harness did not resolve them. Full-matrix allocation
  profiling in `resource-attribution-v4-1132.log` locates runtime prefix tuples
  in `_particle_cell_strides` (D2: 16 B; D3: 48 B, including empty populations).
  Replacing those tuples with a direct prefix-product loop passed the unchanged
  108-profile zero-allocation gate in `resources-v4-{112,1132}.log` on both
  runtimes. Earlier condensed descriptions incorrectly associated these bytes
  with nonempty populations; the raw logs and this corrected ledger are authoritative.

- `null-oblique-v1-1132.log`: four of eight null-angle checks failed from
  angle-independent normalization/renormalization noise. The independent
  `small-angle-v1-1132.log` then failed 48/60 increment-accuracy gates, showing
  a zero-only special case would not fix the neighboring numerical defect.
- `extreme-increment-v1-1132.log`: eight finite physical-update gates passed,
  but two rejected-proposal cases raised a fictitious range error. The local
  update now validates accepted proposals and tests accepted-overflow rollback.
- `expanded-v2-1132.log`: the stable-increment repair passes all expanded
  physical, small-angle, extreme, storage, generator and thermal gates.
  `legacy-v2-{112,1132}.log` retains all 8,064 exact differential assertions
  per runtime, including scratch, permutation, RNG state and warmed bytes.
- `neighbors-v1-112.log`: all 49,720 neighboring assertions passed before the
  stable-increment change. The later latest-Julia result is recorded below.
- `format-v1-1132.log`: all 20 AST-preservation/idempotence gates passed.
- `final-correctness-v1-{112,1132}.log`: all permanent local-collision tests,
  the unchanged 108-profile zero-allocation matrix, and the runnable normalized
  example passed on both runtimes after formatting. Later static review found
  that the original seeded SI↔normalized comparison used identical species
  velocities and therefore exercised a no-op; that limitation is detailed below.
- `quality-v3-{112,1132}.log`: 55 concrete JET paths and nine Aqua gates pass
  per runtime. Full-package diagnostic leaf/call-tree comparison remains under
  review; equal report counts alone are not acceptance.

- `legacy-other-v1-1132.log`: expanded legacy Float16/BigFloat coverage passes
  660/672 checks, exposing 12 active BigFloat allocation regressions from eager
  computation of the unused increment. Values and RNG continuation still match.
  The canonical sampler now accepts an internal requested-representation argument;
  legacy requests only its relative vector and skips increment arithmetic.
- `format-v2-1132.log`: the requested-representation repair passes all 20
  formatter AST/idempotence checks.
- `legacy-other-v2-112.log`: all 672 expanded precision checks pass. On Julia
  1.13.2 the same source passes 669/672: all values and RNG states match, while
  three exact-byte gates fail because candidate bytes are lower (15,744 versus
  15,752; 12,960 versus 12,968; 108,072 versus 108,344). The original 12 eager-
  arithmetic regressions are absent. The exact-byte contract is still pending;
  no gate has been relaxed. Full allocation attribution and identical-source
  controls are prepared to identify these remaining discrepancies.
- `neighbors-v2-1132.log`: all 49,720 neighboring assertions pass after the
  stable-increment repair, before the requested-representation optimization.
- `quality-v4-1132.log`: eight direct sampler/legacy compiler checks pass in
  addition to the 55 local paths and nine Aqua checks.
- `package-jet-v3-{112,1132}.log`: unrestricted raw signatures, complete frame
  lists and leaf messages are preserved for every parent/candidate diagnostic.
  Julia 1.12.7 reports 30/30 and Julia 1.13.2 reports 27/27; equal counts are
  not acceptance. Respectively 26 and 23 complete blocks match after stripping
  only package namespaces and source roots. Three additional blocks differ
  only in generated keyword-method indices in unchanged source. The remaining
  generic diagnostic reaches the same `Base.hypot`/`bitcast` leaf through a
  changed package entry site: the parent's relative-speed calculation versus
  the candidate's rotated-direction calculation. The complete unmatched
  blocks are retained in `package-jet-raw-diff-{112,1132}.txt`; this is not
  described as an identical source call. Concrete public and sampler checks
  pass in the earlier source state. The requested-representation repair still
  needs its own compiler recheck.
- `cost-screen-v1-112.log` passes all 216 legacy value/RNG/allocation checks;
  `local-cost-screen-v1-1132.log` passes all 216 new-operator allocation checks.
  These are untimed harness screens. The frozen latency controls, old/new
  comparisons and new-operator scaling measurements remain pending.
- Static coverage review added a four-species coupled-channel generator oracle
  to the permanent statistical tests. It sums independent analytic drift and
  diffusion rates at the initial velocities, with unequal masses, charges and
  weights, disabled matrix entries and reversed species order. Its 256 cells,
  128 seeds, timestep and seven-standard-error plus finite-step budgets match
  the frozen pairwise gates. The pre-measurement fixture is retained as
  `statistics-coupled-oracle-before-first-run.jl`; these 32 new checks have not
  yet run during the coordinated native timing hold.
- Static source review found both one-particle species in the original
  SI↔normalized fixture inherited the same loader seed and velocity vector.
  Its independent dimensional-generator algebra remains meaningful, but the
  seeded physical-update comparison does not establish scattering equivalence.
  The original file is retained as
  `local-tests-before-normalization-oracle-repair.jl`, with a pending failing
  witness in `normalization-null-fixture-witness.jl`. The replacement fixture
  uses distinct representable velocities, asserts nonzero updates for both
  species and both unit systems, checks a genuinely weak-angle rate rather than
  the saturated isotropic branch, and covers active/inactive velocity floors
  with explicit SI rescaling. These replacement checks remain untested during
  the native timing hold; prior green results are not their acceptance.

## Historical candidate composition and confirmed follow-up findings

- [x] Compose on accepted normalization without replacing its prefix.
- [x] Existing 1,716 correctness checks pass on Julia 1.12.7 and 1.13.2.
- [x] Normalization helper allocation witness repaired: 162 checks pass each, maximum 24,312 bytes versus 9,854,776 before.
- [x] Expanded complete-angle defect reproduced on both runtimes.
- [ ] Complete-angle repair, independent recheck and permanent regressions.
- [ ] Helper full range/rounding/rejection checks and permanent tests.
- [ ] Final resource, legacy, compiler, fixed-cost, integration and source review.
- [ ] Register new tests, update public docs/example, commit/push and exact-SHA CI.

The corrected non-null SI↔normalized fixture and all 32 coupled four-species
checks now run as part of the 1,716-check suite on both runtimes. Earlier
“untested during timing hold” entries above describe historical source states.
`helper-v2-commands.json` retains the successful composition/resource runs.
`complete-angle-v2-commands.json` retains 468 passes / 48 failures per runtime;
this implementation is still held. The original failing logs are retained.
The normalization helper now evaluates the complete stored-reference monomials
with bounded significands, rather than expanding rounded reference lengths.
Broader helper validation and all final acceptance gates remain pending.

### Latest completed angular iteration (V3)

- [x] Original extreme-angle matrix: 516/516 on both runtimes.
- [x] Existing local suites: 1,716/1,716 on both runtimes.
- [ ] Oblique physical oracle: 22,156 passes / 164 failures on each runtime.
- [ ] Independent center-of-mass reconstruction of actual sampled angles:
  22,181 passes / 139 failures on each runtime. All failing profiles use
  reciprocal-tangent backscatter with oblique Float64 velocities.
- [ ] Reused-workspace allocations: 96 passes / 12 failed profiles on each
  runtime. Failed profiles are Float32, 4,096 particles and 2 or 4 species;
  the recorded calls allocate nonzero bytes in these ordinary profiles.

The original physical oracle is retained. The sampled-angle reconstruction
separates integration loss from sensitivity to the rounded physical rate;
neither equal report counts nor a passing smaller matrix establishes acceptance.
V4 is repairing the reciprocal-tangent update. Numerical, resource and full
compatibility checks must all pass on the final source before code promotion.
The detailed command records are `complete-angle-v3-commands.json`,
`complete-angle-oblique-v2-commands.json` and
`complete-angle-sampled-v1-commands.json` in the evidence directory above.

## Completed V6 checks and remaining resource work

- Original extreme-angle witness: 516 checks pass per runtime.
- V5 sampled center-of-mass oracle: 22,320 checks pass per runtime.
- V6 combined sampled/conditioned physical oracle: 39,060 checks pass per runtime.
- V6 independent Float32 widening boundary matrix: 14,402 checks pass per runtime.
- Helper physical rounding cells: 7,590 checks pass per runtime; compact MPFR
  boundaries add 100, followed by a successful final settings-restoration check.
- Original physical residual-relative oracle retains 25 failures. The companion
  diagnostic verifies rate, norm and phase rounding and compares physical updates
  on their conditioning scale while preserving strict sampled-update checks.
  Measured maximum relative sampled-rate error is 9.08105e-6 (Float32) and
  3.35716e-13 (Float64), within the frozen 4096eps budgets. The old failures remain
  evidence; they are not relabeled as passes.
- VERIFIED: `complete-update-inference-v1-1132.log:348` identifies the captured
  rotation as Core.Box. Full allocation stacks in both v1 attribution records
  locate boxed arithmetic. V6 assigns the branch result once; the same direct
  widened helper now allocates zero bytes on both runtimes, versus 2,000 before.
- The full ordinary resource gate remains failed: 84 of 108 profiles pass,
  24 fail on each runtime. Further Float64 completion optimization is required.

The implementation remains isolated and unaccepted. Permanent regressions,
registration, final compiler/legacy/neighbor checks, frozen costs, independent
review, code promotion and exact-SHA CI remain required. The source checkpoint
`composition-v6-source.json` predates this documentation-only update; numerical
sources are unchanged. No fixed timing acceptance run has begun for this version.

## Completed V7 resource repair and fresh numerical finding

V7 evaluates the existing complete physical expression with bounded compensated
Float64 operands, retaining exact completion outside its range/conditioning
certificate. Both runtimes pass the original 516, conditioned 39,060, Float32
boundary 14,402 and Float64 boundary 60,002 checks. All 108 ordinary allocation
profiles now record zero bytes; four direct helper paths also record zero bytes
and no concrete JET reports. The exponent/error certificate still requires an
independent boundary audit before final acceptance.

A fresh mixed-component witness retains 96 passes / 48 failures per runtime:
both normalized direction and angle-square can be normal while their product
underflows. Multiplication by physical speed would produce a finite component
kick. V8 is repairing the canonical local increment's product and sum guards.
The original legacy representation and random draw ordering remain compatibility
gates. This feature is not promoted; formatting, permanent test registration,
neighbor/legacy/full compiler/fixed-cost checks and final review remain required.

Evidence: `complete-angle-v7-commands.json`,
`compensated-boundaries-v2-commands.json`,
`complete-update-allocation-v3-commands.json`, and
`mixed-component-angle-v1-commands.json` in the retained evidence directory.

## Recorded V9 evidence and acceptance gates

The V8 product/sum guards pass the previously failing 144-check component
matrix and all 108 ordinary allocation profiles on both runtimes. Its registered
suite passes 115,951 assertions per runtime. Its neighboring suites pass 49,720
and its example runs successfully on each runtime.

A different cancellation witness then finds a Float32 widening defect. With
unequal physical masses, a small azimuthal term disappears inside the rotated
increment before the old final-sum guard can inspect it. Both runtimes retain
36 passes / 12 failures in `widened-cancellation-witness-v2`. The equal-mass
control passes because another coordinate already requests exact completion.
V9 checks all final coordinates against the original pair scale, including a
rounded-zero change, and retains the same canonical exact physical expression.

- [x] V9 permanent suite: **140,303 assertions per runtime**, including 24,304
  structured cancellation/boundary gates and the 48-check witness. Normalization
  subprocesses pass 7,590 + 100 checks and confirm settings restoration.
- [x] All 108 ordinary reused-workspace profiles record zero bytes per runtime.
- [x] Concrete compiler/quality checks: 55 local, 10 sampler/legacy, six complete
  representation checks and nine Aqua checks pass on each runtime. Complete
  package trees match V8 after checked unchanged-line mapping (29/27 reports).
  Compared with the accepted parent, one generic entry moves between two
  `hypot` calls; every downstream frame, signature and message is identical.
  Concrete checks and execution cover that changed entry separately.
- [x] Formatter preserves ASTs and is idempotent for all 20 owned Julia files.
  The first integrated V9 fixture omitted `Random`; both failed runs retain
  140,144 passes and one fixture error. The import is corrected, and both full
  registered suites were rerun successfully. Production code was unchanged.
- [x] Strict Float16/BigFloat allocation and value comparison passes all 672
  checks per runtime after warming the exact measured state. The earlier
  latest-Julia failures are retained: an identical-parent control reproduces
  the same 8/8/272-byte differences, and reversed measurement order reverses
  them. Five reset measurements and full allocation profiles agree for each
  affected workload. This verifies measurement-order dependence; it does not
  identify an internal Julia/MPFR allocation owner. No padding or looser byte
  assertion was introduced.
- [ ] Independently certify the exponent/error bounds for both hardware
  completion adapters; random and structured passing tests alone are not proof.
- [ ] Execute the frozen fixed-cost campaign: 72 legacy and 216 new-local
  profiles on both runtimes, identical/doubled controls, zero measured
  compilation, full outputs/RNG, strict ordinary bytes and recorded storage.
  No timing run has started. Preserve any failed control without unchanged retry.
- [ ] Review the final API, example, diff and documentation; promote the accepted
  implementation separately to `main`, then verify every exact-SHA CI job.

Evidence: `composition-v9b-permanent-source.json`,
`completion-v9b-commands.json`, `neighbors-current-v1-commands.json`,
`quality-current-v2-{112,1132}.log`, `package-jet-current-v2-*-review.json`,
`legacy-other-same-state-v1-commands.json`,
`legacy-allocation-boundary-review-v1.md`, `completion-certificate-review-v2.md`
and `fixed-cost-v2-design.md`. CAP-003/004 parallel support and CAP-010 composed
evolution/RNG restart remain dependent implementations, not completed features.

## Recorded V10 evidence and then-active implementation tracks

V10 is composed on `c5e8472` and remains unaccepted. A fresh witness verifies
Float32 overflow-midpoint double rounding: the widened intermediate reaches the
midpoint while the exact stored-sample expression rounds to finite Float32.
Both runtimes retain 24 passes / six failures before repair. The widened adapter
now requests existing exact completion above Float32's finite maximum; it does
not clamp outputs or change the sampled physics.

- [x] Original witness: 30/30 checks per runtime.
- [x] Permanent boundary coverage: 810 component/neighbor checks and 30 public
  finite-acceptance/true-overflow rollback checks per runtime.
- [x] Registered suite: **141,143 assertions per runtime**, plus 7,590 + 100
  normalization subprocess checks and confirmed settings restoration. All 108
  ordinary reused-workspace profiles record zero bytes.
- [x] Strengthened 448-error-scale audit: 16,198 checks; adjacent exponent,
  coefficient and conditioning boundaries: 5,564 checks per runtime. The latter
  retains every input case; the stronger upper guard routes nine more profiles
  through exact completion, removing 54 inapplicable hardware-component checks.
- [x] Concrete quality checks: 55 local, 10 sampler/legacy, six complete
  representations and nine Aqua checks per runtime. Full compiler trees match
  V9 after checked unchanged-line mapping (29/27 reports).
- [x] Formatter: 40 AST/idempotence checks across 20 owned Julia files.
- [x] Performance harness preflight: 104 mutation-sensitivity assertions and
  56 tiny replay/resource profiles per runtime. No timing acceptance run.
- [x] Final source-level error/range review is recorded in
  `completion-certificate-derivation-v5.md`. The manual proof checks primitive
  normalization, bounded intermediates, underflow propagation and both rounded
  outputs against the actual Julia implementations. It is not machine checked;
  changed runtime arithmetic requires review again.
- [x] Exceptional cost preflight passes 30 state/doubled-replay mutation checks
  and eight tiny profiles per runtime; allocation stacks confirm exact
  completion. The complete harness is frozen at 304 profiles per runtime,
  including controls. The first preflight retains an RNG fixture ambiguity,
  corrected through concrete Float32/64 methods before the successful runs.
- [x] V10 legacy differential rerun passes all 8,064 assertions per runtime,
  retaining complete outputs, scratch, RNG and exact byte equality.
- [x] V10 neighboring regressions pass 49,720 assertions per runtime and the
  public two-species example runs successfully on both runtimes. Commands and
  log hashes are in `neighbors-current-v2-commands.json`.
- [ ] Execute the exclusive fixed-cost campaign and retain every failed control.
  No timing acceptance run has begun; other live Julia work prevents exclusivity.
- [ ] Reconcile final example/neighbor/legacy/API/docs evidence, promote the
  accepted implementation and require all six exact-SHA GitHub CI jobs green.

Current evidence: `composition-v10-permanent-source.json`,
`widened-overflow-midpoint-v1-commands.json`,
`widened-overflow-midpoint-v2-commands.json`,
`widened-overflow-regression-v1-commands.json`, `completion-v10-commands.json`,
`package-jet-current-v3-*-review.json`,
`completion-certificate-derivation-v5.md`, `fixed-cost-v2-final-freeze.json`,
`fixed-exceptional-selfcheck-v1-commands.json` and `fixed-cost-v2-design.md`.

The recorded V10 handoff assigned separate linear-CN, ES-species and bounded-
Maxwell work. That is historical ownership, not a claim about current agents.
This standalone operator does not complete broader CAP-001, the remaining
CAP-010 source/restart combinations or dependent parallel backends. The
[full implementation roadmap](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md)
retains every CAP-001–011 acceptance gate.
