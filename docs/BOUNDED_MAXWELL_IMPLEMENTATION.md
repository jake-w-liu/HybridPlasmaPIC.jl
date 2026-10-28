# Bounded Maxwell implementation record

Scope: CAP-002A, source-free Cartesian vacuum/cavity fields, plus the 1D
[prescribed-current prerequisite](SOURCED_YEE_FIELDS.md) of 002B (see
[below](#prescribed-current-slice)). This does not implement bounded particle
coupling (002B), absorbing layers (002C), or restart, GPU/MPI and full
production integration (002D). CAP-002 remains open.

The isolated B2a prerequisite now has extensive nodal charge and longitudinal
endpoint-impulse source, tests and documentation. Its six-file formatter pass
checks the complete executable AST under an explicit source-line metadata
contract, idempotence and 12 checker controls. The immutable postformat source
has 616 files. All seven environment/provider commands exit zero on actual
Julia 1.12.7/1.13.3; full package/selected-provider records match per runtime.
Root accepts the exact cache/native/header/embedded-source metadata proof on
both runtimes (248/245 assertions). Both syntax parsers pass. The first normal
minimum-runtime focused run retains 119,470 passes and 12 allocation failures:
192 B against the unchanged 0 B gate. Its numerical/neighbor assertions pass.
A separate paired diagnosis retains exit 1 and all 12 profiles / 48 traces:
every Float32 direct/helper call allocates 192 B (including empty populations),
while sampled Float64 calls allocate 0 B. All six 32 B tuple allocations per
Float32 trace come from the heterogeneous-buffer indexing in alias admission.
An isolated repair spells out the same three alias predicates in the same order.
Both refreshed focused lanes pass all 119,482 assertions, including the unchanged
48 allocation checks. Each paired diagnosis passes 87 checks; all 288 byte
measurements are zero and all 48 production traces are empty, with four separate
no-op records retained. The next CRC stage stops at evidence serialization:
`TOML.print` rejects tuple-valued geometry weights at the external driver line 25
(10 passes, zero failed assertions, four errors). That failure remains recorded;
the serialization-only correction passes 90 controls on each runtime. The
corrected CRC lanes each pass 1,806 assertions: 430 geometry records, 390
displacement records and 66 separate public outcomes are retained. The next
minimum compiler gate stops with 32 passes and four failed zero-report
assertions. All 34 full captures are retained: the four abstract empty-vector
charge/impulse roots each report two Base reinterpret/padding diagnostics; the
other 30 roots are empty. Both normal runtime investigations complete: each
retains 18 cases, 108 allocation measurements, 12 complete traces and 144 IR
files. Independent exact charge/overlap checks, tuple parity and input/storage
checks pass. Broad empty vectors allocate 160 B per operation; broad nonempty
and mixed Array/view storage also allocate. Concrete controls remain 0 B.
Complete raw/index/source review is done; allocation attribution does not prove
every producer cause. Profile.Allocs.CorruptType is the installed decoder
low-pointer sentinel, with the original type/pointer cause unresolved.

A new isolated candidate now contains the reviewed population/endpoint method
barriers and a count sentinel derived from the unchanged 2^51 native limit.
Global validation order, shared native/exact particle arithmetic and final
publication remain explicit invariants. Targeted broad/mixed/custom-container,
count, rejection and liveness tests are written. Both runtime syntax checks and
formatter AST/idempotence pass; all original test bytes remain. Seven fresh
offline provider stages pass, with all 160/161 package and 440/444 loaded-library
records reviewed. Both metadata decoders pass 248/245 assertions and match all
203 embedded/live/frozen source records per runtime. Only the exact four fields
of the observed package-image row differ per runtime; every other loaded row
is equal. The fresh normal campaign passes 126,925 focused assertions, 87 ordinary
allocation checks and 1,806 CRC assertions per runtime. Every one of the 288
ordinary byte measurements per runtime is zero; all 48 production traces per
runtime are empty. The minimum compiler stage then stops at 58 passes and 26
failed zero-report assertions. All 82 full captures remain: the original 34
roots have four failures, and the 48 added roots have 22. Latest compiler and
later vector/quality/resource/full-package stages are unrun. This candidate
remains unaccepted; all previous failures remain recorded. Diagnose the full
reports and owning source before another repair or campaign. It does not
implement a coupled force, particle pusher, wall event or restart. The
[audit tracker](DEEP_DEBUG_TRACKER.md) owns current acceptance.

## Discrete contract and derivation

For component c, E is centered in axis c and nodal in other axes; B is nodal in
axis c and centered in other axes. Suppressed spatial directions have zero
derivative, retaining all three vector components. Periodic nodal axes have n
entries; PEC nodal axes have n+1. Centered axes have n entries.

The forward node-to-center difference is `(f[i+1]-f[i])/dx`. Its negative adjoint
is the backward center-to-node difference. The even extension of tangential B
at a conductor sets its normal derivative to zero at the wall. Tangential E
and normal B are constrained zero there, including intersections. There is no
trapped normal magnetic flux in this homogeneous cavity contract.

Let C map E edges to curl E faces, with Cᵀ mapping B faces to curl B edges on
the constrained spaces. The two discrete divergences annihilate C and Cᵀ by
commuting tensor-product differences. Magnetic divergence is cell centered;
electric bulk divergence uses interior PEC nodes (surface charge is excluded).

The time step is

```
B½ = Bⁿ − dt C Eⁿ / 2
Eⁿ⁺¹ = Eⁿ + dt c² Cᵀ B½
Bⁿ⁺¹ = B½ − dt C Eⁿ⁺¹ / 2
```

Thus both stored fields are at the same integer time. For a curl singular value
κ, `sin(ω dt/2) = c dt κ/2`; the maximum κ is bounded by
`2 sqrt(sum(1/dx²))`. The implementation requires the strict bound
`dt < 1/(c sqrt(sum(1/dx²)))`, excluding the repeated-eigenvalue endpoint.
The constructor computes and caches the downward-rounded bound with the existing
exact square-root utility, so reciprocal-norm overflow or upward rounding cannot
weaken the stability condition.
For constant dt the discrete quadratic invariant is
`(||E||² + c²||B||² − dt² c²||CE||²/4) prod(dx)/2`.
Physical field energy omits the last term and varies at second order in dt.
Variable dt retains synchronized second-order centering, but neither this
fixed-dt invariant nor stability under an arbitrary sequence of individually
sub-CFL steps. For the periodic one-dimensional Nyquist mode with unit spacing
and wave speed, alternating dt=0.9 and 0.4 has a two-step amplification matrix
with determinant 1 and trace -1454/625 < -2. Both steps satisfy the strict bound,
yet that sequence has a growing mode. No adaptive-step stability policy is
supplied. Gauss residuals propagate without a corrective projection.

This derivation uses the standard Yee staggered method; see the primary
[numerical-dispersion study](https://arxiv.org/abs/1710.06829).

## Numerical contract (normwise)

The code follows the [normwise standard](DEEP_DEBUG_TRACKER.md#numerical-accuracy-standard-user-decision-2026-10-26).

- **Stages.** Each stage is the native stencil
  `new = old + (κ_j Δ_j − κ_k Δ_k)` with multipliers `κ_d = −dt/(2dx_d)`
  (magnetic half-kick) or `dt c²/dx_d` (electric stage) computed once per stage.
  Against the exact expression from the stage's actual inputs,
  `|new − exact| ≤ 8 eps(T) (|old| + Σ|κ_d Δ_d|) + nextfloat(0)`. Every operation
  is relative-accurate when its result is normal; two underflowed products
  contribute at most one subnormal quantum. Componentwise accuracy of a
  nearly cancelling entry is not certified, and this is not a trajectory bound
  over several stages.
- **Range.** Ordinary multipliers require every intermediate finite and
  normal. Otherwise Float32 stages run the same kernel in Float64 (every binary32
  multiplier, difference and product is then a normal Float64), and Float64
  stages use `κ_d = m_d 2^s` with normal mantissas `m_d ≥ 2^52` (no difference
  product underflows) and apply `2^s` as two exact power-of-two factors. Without
  such a split, every entry takes the exact expression. After each native pass a
  vectorized finiteness reduction detects overflow; nonfinite entries are
  evaluated exactly (`Rational{BigInt}`) from the same inputs and rounded once.
  A nonfinite exact result rejects the step.
- **Diagnostics.** Float32 energy, Poynting and face flux are evaluated in Float64,
  which represents every binary32 square, product, interpolant and grid factor
  without range loss, and are rounded once. Float64 diagnostics are native with
  pairwise summation; overflow, or a nonzero intermediate leaving the normal
  range before later scaling, selects the exact expression. Cell Poynting entries
  meet `8 eps(T) c²(|E_j||B_k| + |E_k||B_j|) + nextfloat(0)/2` (absolute
  interpolants), face integrals add `N eps(T)` for N face cells. `yee_curl!` and
  `yee_divergence!` complete overflowing entries exactly; a nonfinite output is
  the rounded exact value.

## Implementation

Owners: `YeeGrid` and coordinates in `meshes/cartesian_mesh.jl`; shape checks,
differences, the row engine, `yee_curl!` and `yee_divergence!` in
`meshes/local_finite_difference.jl`; storage in `fields/yee_fields.jl`; PEC
constraints in `boundaries/yee_conductor.jl`; stages and the transaction in
`models/maxwell_vacuum.jl`; energy and Poynting output in
`diagnostics/yee_diagnostics.jl`.

The row engine evaluates one output component as unit-stride loops along axis 1
with fixed row offsets for both difference terms (no per-entry index logic).
Periodic wraps and absent backward differences at PEC walls are scalar row-end
entries; absent differences are typed `nothing`, never an arithmetic zero.
Loops carry `@simd ivdep` (callers reject output/source aliasing; the in-place
half-kick reads and writes only the same entry); the only reduction is an
integer finiteness flag, so IEEE association and the absence of FMA contraction
are preserved and results are bitwise reproducible.

`step_maxwell!` validates step size, CFL, clock, shapes, aliasing and PEC data
before staging writes, then runs the three stages into the candidate arrays.
The first two stages have distinct originals and are completed exactly where
their native pass is nonfinite. The final half-kick is in place; if its native
pass fails, the first half-kick is rebuilt from the unchanged physical fields
and the final half-kick is evaluated entrywise with exact completion. Nonfinite
physical inputs fail a stage reduction and are rejected before any exact
completion. Homogeneous PEC data stay exactly zero by construction (each
constrained entry is `0 + κ(0 − 0)` or its exact replay), so no clamp or
re-validation follows. Publication copies the candidate into the physical arrays,
which keep their identities; time and step count advance last.

## Verification

`test/test_bounded_maxwell.jl` and `test/bounded_maxwell/*` keep every V12
physics and contract gate: geometry, constraints and rejection; compatible
curls/divergences and adjoints; periodic discrete dispersion; PEC/mixed cavity
modes and corners; spatial/temporal second order; pulse transport and
reflection phase; energy and Poynting diagnostics; variable-step centering;
periodic cross-check against the spectral vacuum `EMPIC`; CFL certificate.
Former componentwise oracles (32 eps per stage entry, subnormal interpolation,
exact equality of cancelled stage or Poynting values) are normwise bounds
against independent Rational oracles that use the stage's actual inputs and an
independently computed absolute scale. Exact-zero/identity cases, overflow
witnesses (complete cancellation of overflowing products, rejected overflowing
stages) and subnormal witnesses (whose allowance is below one quantum) remain.
New gates: bitwise determinism and zero warmed allocations for every dimension
and boundary mix, operator range completion, and public steps equal to the
composition of their checked stages.

## Resources and performance

The stepper owns physical E/B and one candidate E/B. Warmed steps allocate zero
bytes. On the shared Apple M5 (Julia 1.12.7, load average 6–11 during
measurement), best-of-15 48³ PEC steps cost about 6.4 ns (Float32) and 11 ns
(Float64) per cell; the same session measured about 1,200 ns per cell for the
componentwise V12 candidate. Each stage costs about as much as a minimal
unit-stride kernel with the same three updates per cell timed on the same
machine. Energy costs about 0.2–0.4 steps,
cell Poynting 2–3.3 steps, face flux at most 0.75 steps and `yee_curl!` 0.3 steps.

## Prescribed-current slice

`step_maxwell!(model::MaxwellVacuum{1}, dt, impulse)` adds a caller-supplied
current impulse to the electric stage, `(old + κΔ) − impulse`, through the same
row engine: impulse entries travel with the original entries and the stage is
evaluated in the stage's arithmetic type (Float64 for out-of-range Float32
multipliers) before one rounding; exact completion includes the impulse. A zero
impulse entry, including `-0.0`, keeps the vacuum value bitwise. Impulses are
validated (shape, precision, ordinary vectors, finiteness, zero tangential PEC
nodes, independence from all twelve model arrays) before any staging write.
`yee_electric_charge!` is the all-node 1D incidence (one correctly rounded
difference per entry); `yee_source_work` and `yee_modified_energy` follow the
diagnostic rules above with scales `dx/2·Σ|J|(|Eold|+|Enew|)` and the sum of the
three nonnegative modified-energy contributions; the modified-energy curl term is
summed as `((dt c/(2dx))ΔE)²`, which is dominated by the E² sum because
`dt c < dx`. `test/test_sourced_yee*.jl` keep the V12 sourced gates (independent
stages, charge identity, work and modified-energy balance, vacuum equality,
complete-expression range witnesses, transactional final-stage failure, discrete
forced modes, modal time/space order, validation and rollback) with per-stage and
scalar normwise bounds in place of the componentwise ones.

## Remaining CAP-002 work

002B must add trajectory-based current deposition and a wall charge/current
ledger before this becomes bounded EM PIC. 002C needs actual absorbing-layer
auxiliary evolution and a measured angle/frequency/refinement reflection suite.
002D must register the synchronized fields, model clock and boundary state in
checkpoint/restart and backend orchestration. None is supplied by this vacuum
slice, and the existing particle boundary utilities are not a substitute.

## History

The unlanded V12 candidate certified a componentwise `32eps` bound per stage
entry with enclosure arithmetic, compensated retries and exact replay in the
hot loop (766 / 795 ns per cell per step at 48³ in the lead's measurement).
This record replaces it under the normwise standard; its enclosure helpers,
V-numbered numerical/scaling trackers and enclosure-arithmetic tests are not
carried forward. Retained V-series evidence is in
`HybridPlasmaPIC-audit-evidence/2026-10-23/bounded-maxwell/`.
