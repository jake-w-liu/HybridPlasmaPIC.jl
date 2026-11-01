# Relativistic normal-range accuracy

Status: the initial manual GitHub run failed during package precompilation on
both Julia versions, before observing any candidates. The carrier specialization
is now defined after its type; corrected-source loading and observations remain
pending. No normal-range numerical correction or accuracy acceptance is claimed.
The [factor plan](PARTICLE_CHARGE_MASS_RANGE_PLAN.md) owns P3 factor routing;
[C1](RELATIVISTIC_OUTPUT_PUBLICATION_PLAN.md) owns supported output publication
and public rollback. Neither establishes a full normal-range accuracy theorem.

## Checked owners and remaining certificate

On source revision `ebad36e733a41dcc4d720542dfa3ee6d814d3b0d`, scalar
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
initial failure supplies no candidate accuracy result; loading and observations
must run again against the corrected revision.

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

## Bounded implementation decision

If the remote observation confirms a witness, first repair its canonical owner
and both bypasses. Retain valid native bits, original q,m, explicit-qm meaning,
the existing four-ulp normalized Euclidean scale, C1 domain and public undo.
Do not duplicate the Boris kernel, clip a limiting velocity, widen a budget or
use an unproved cold result as an accuracy certificate.

For a priming-owned defect, the existing residual-ball expression can enclose
original q/m by division of stored q and m. An outward residual bound can compare
native output with the physical component interval at a conservative lower
bound of the existing normalized term scale. Derive a cheap sufficient native
batch bound first; blanket DD evaluation of ordinary particles adds unmeasured
cost and is not an optimization. Keep the existing operand/working range and
nearest/gradual assumptions explicit. A full Boris or trajectory theorem is a
separate remaining task.

- [x] Inspect native/cold relativistic kick, both EMPIC callers and both primes.
- [x] Prepare a compact independent, two-precision remote observation.
- [ ] Run the observation on actual minimum/latest Julia and inspect every row.
- [ ] Confirm a concrete witness or record that these candidates do not fail.
- [ ] Derive and implement a sufficient canonical owner correction if supported.
- [ ] Add the confirmed essential regression and verify exact-revision package,
  compiler, formatter and MPI checks in GitHub; measure applicable cost changes.
