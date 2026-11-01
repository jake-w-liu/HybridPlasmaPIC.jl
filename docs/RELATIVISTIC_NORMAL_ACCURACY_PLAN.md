# Relativistic normal-range accuracy

Status: the bounded correction at `308cb45` passes the fixed-corpus
[remote observation](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41280483562)
and [focused particle essentials](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41280485822)
on actual Julia 1.12.7 and 1.13.3. The earlier six production kick failures at
`7189d75` remain recorded below. Full package/compiler/formatter/MPI acceptance,
extreme replay accuracy and resource/performance gates remain pending.

## Bounded kick correction: historical pre-edit brief

The corrected GitHub observation at `7189d75` completes all 60 stages on Julia
1.12.7 and 1.13.3. Each runtime records the same six kick violations: Float16
rows 3/4 and Float32 row 3, each for original q,m and explicit stored qm. The
Float32 row is also admitted by the range-only batch predicate. Both independent
512/1024-bit references agree; none of these 60 observations finds a prime
violation. The earlier package-load failures remain failures.

This repair targets those scalar and actual EMPIC batch kick paths. One private
native/DD residual adapter will evaluate the existing `_boris_rel_momentum`
body, preserving its historical native order. Original factors, the retained
four-ulp Euclidean scale, finite-result margin, C1 publication and public undo
remain required. A joint batch admission will combine the existing range owner
with a sufficient finite-c error bound; range-only callers cannot bypass it.
Residual-certified DD repair is bounded by its working range and nearest/gradual
context. Unresolved extreme replay still uses the inherited 256-bit owner and
does not gain a general accuracy theorem. Generic arithmetic, nonfinite
propagation, priming and unsupported hardware modes retain their existing paths.

Every original 100,000 legacy input per type remains. Every changed tuple is
checked against the existing independent 4096-bit Cayley reference: its old
value must exceed the unchanged budget and its replacement must satisfy it and
C1. A difference selects verification work; only that physical reference decides
whether the change is permitted. Unchanged tuples receive no accuracy label.
Fixed independent regressions exercise the confirmed rows, original/explicit
meanings and both EMPIC callers. All new runtime, formatter and compiler
verification runs on GitHub only.

- [x] Read canonical scalar/native/batch/cold owners, callers and retained tests.
- [x] Review both corrected remote observations and the weak-gamma design limit.
- [x] Implement and source-review canonical bounded scalar/DD and batch admission.
- [x] Pass fixed-corpus observations and focused independent regressions remotely.
- [ ] Pass the retained legacy gates in the full remote package suite.
- [ ] Review exact-revision package/compiler/formatter/MPI results and costs.

The [factor plan](PARTICLE_CHARGE_MASS_RANGE_PLAN.md) owns P3 factor routing;
[C1](RELATIVISTIC_OUTPUT_PUBLICATION_PLAN.md) owns supported output publication
and public rollback. Neither establishes a full normal-range accuracy theorem.

## Checked owners and remaining certificate

Before this correction, on source revision `ebad36e733a41dcc4d720542dfa3ee6d814d3b0d`, scalar
`_boris_rel` admits native results through factor, range and supported-output
checks. `_boris_rel_native` returns a supported native tuple after the batch
range route. Both EMPIC kick owners use that route. Their canonical expression
is `_boris_rel_momentum` in `src/models/full_pic.jl`.

Both relativistic priming owners use `_rel_prime_in_range` and
`_rel_prime_batch_certificate` in `src/particles/kinematics_range.jl`. These
check range and conditioning. They do not contain the separate sufficient
four-ulp accuracy bound now used by nonrelativistic Boris and Euler. This is a
source-confirmed missing certificate; a particular admitted numerical failure
still requires a checked witness.

Cold priming already has residual-certified double-word component intervals.
The current species route uses them only when stored `qm` equals original
rational `q/m` in the admitted IEEE context. Other cold cases retain the existing
256-bit replay. Fixed replay precision is not a general rounding proof.

## Manual GitHub observation

The [initial run at `14568b8`](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41277054066)
fails on both actual Julia 1.12.7 and 1.13.3 at
`src/particles/kinematics_range.jl`'s typed coefficient method. The module loads
that file before `boris_accuracy.jl` defines `_BorisSpeciesFactors`. The
specialization now lives beside that carrier's factor helpers with its body
unchanged. The generic scalar method and include order are retained. The
initial failure supplies no candidate accuracy result. The corrected run loads
the package and supplies the separate numerical findings below.

Run the workflow **Relativistic normal accuracy observation**. It uses Julia
`1.12` and moving latest-stable `1`, records the actual version and instantiates
the ordinary package project. It does not run the full package suite or launch
local computation.

The script `.github/diagnostics/relativistic_normal_accuracy.jl` retains the four
stored Float16 candidates from `test/test_boris_accuracy.jl`, scales velocity and
electric field by exactly one quarter, and sets `c=128` so `c²=16384` is finite
and normal. Exact Float32/64 widenings and one ordinary control remain. Original
q,m and explicit stored-qm meanings are classified separately: 30 fixed input
cases and 60 stage observations, with no random search or hidden truncation.

An independent Cayley-form kick and momentum-backup priming reference are
evaluated at 512 and 1024 bits. Rounded outputs, budget binades and error
classifications must agree; disagreement fails the observation job. Every row
prints original inputs, native range/batch admission, raw and actual outputs,
both references, errors, budgets and explicit violation labels. Complete logs
are retained as workflow artifacts. A raw native failure outside its admission
is distinct from a production result exceeding the budget. Job exit zero means
the observation completed; it does not turn a printed violation into acceptance.

## Verified observation at `7189d75`

Both raw job logs record the same six production failures: Float16 candidate
rows 3 and 4 under both original-factor and explicit stored-qm meanings, and
Float32 row 3 under both meanings. These are kick observations; none of the
30 priming observations exceeds the budget. This bounded result does not certify
priming in general or inputs outside the fixed corpus.

For Float32 row 3, scalar range admission and batch admission both return true.
The original-factor y-component error is approximately `3.64046e-7` against a
`2.384185791015625e-7` budget. Thus the batch bypass as well as the scalar owner
requires correction. Complete inputs, actual results and both references remain
in the workflow artifacts; the four-ulp contract is unchanged. All Float64 fixed
cases are within budget; that observation is not a general Float64 theorem.

## Verified correction at `308cb45`

The [post-repair observation](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41280483562)
completes all 60 distinct stages on each actual Julia 1.12.7 and 1.13.3. All
512/1024-bit classifications are stable, with zero production violations and
zero unresolved stages. The six raw-native violations remain visible. Every
observation row has batch admission false, so this corpus does not exercise a
successful native batch bypass.

The separate [particle essentials run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41280485822)
passes 25 test groups and 5,687 assertions on each runtime, including nine
ordinary nonzero certified-batch/native-bit checks and 15 zero-native-scale
quantum-budget checks. This accepts those focused regressions, not the retained
100,000-input-per-type legacy gate, full compiler/integration or cost targets.

## Bounded implementation and limits

`relativistic_boris_accuracy.jl` adds a private arithmetic payload for the single
unchanged `_boris_rel_momentum` expression. Native centers use historical IEEE
operations; error radii use exact-sum/product residual bounds and the identities
`A/B-z=(A-B*z)/B` and `sqrt(A)-z=(A-z²)/(sqrt(A)+z)`. The positive denominator
contracts come from this Boris expression, not from a general interval API.
The initial gamma is enclosed from original v,c. Original q/m is enclosed by
division of stored q and m; an explicit supplied qm is an exact scalar input.

Scalar admission requires the output error bound, a finite exact-result margin
and C1. The scale lower bound is below
`max(norm(out), max(norm(gamma0*v), norm(hE))/gamma_final)`. A separate aggregate
finite-c bound combines with the retained range predicate in the joint batch
admission used by both EMPIC kick owners. Storage/alias guards remain required.
Unresolved batches use the scalar route; no all-particle MPFR replay is imposed.
Rows omitted by the nonzero-scale reduction have a separate quantum-floor
proof: native v=0 and RN(hE)=0 imply an exact final norm below
`sqrt(3)(1+kappa)D < 4D`, by magnetic-rotation orthogonality and gamma>=1.

The DD tier evaluates the same body through the existing priming residual-ball
owners. It first preserves a supported native tuple when enclosed within the
same budget, then attempts a supported corrected tuple. Both must satisfy the
unchanged normwise and finite-margin tests. The native gamma seed retains the
existing nonzero v,c window `2^-100 .. 2^100`; the DD attempt applies that window
to every nonzero original operand. All residual tiers require nearest rounding
and gradual underflow. No new persistent storage or numerical public API is added.

Outside those sufficient contexts, or if an enclosure is unresolved, the
inherited exact-gap 256-bit replay remains. This patch does not prove its
arbitrary-extreme rounding or accuracy, full normal-range coverage, unsupported
hardware-mode accuracy, priming accuracy, trajectories, or resource/performance
targets. Generic/nonfinite behavior, C1 domain, publication, rollback and restart
remain unchanged. Directed extreme replay and further admission optimization
remain open work; no fixed-precision convergence is called a proof.

- [x] Inspect native/cold relativistic kick, both EMPIC callers and both primes.
- [x] Prepare a compact independent, two-precision remote observation.
- [x] Run the observation on actual minimum/latest Julia and inspect every row.
- [x] Confirm concrete scalar and batch kick witnesses under both factor meanings.
- [x] Implement a bounded canonical native/DD correction and joint batch admission.
- [x] Add the confirmed essential regressions and pass both-runtime focused checks.
- [ ] Verify exact-revision package, compiler, formatter and MPI checks in GitHub;
  measure applicable cost changes.
