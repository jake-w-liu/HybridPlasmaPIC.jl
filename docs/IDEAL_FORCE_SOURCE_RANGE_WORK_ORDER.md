# Original-factor ideal field and Faraday completion

Status: **pointwise scalar/CGL range and bounded product-range Faraday prerequisites
implemented; broad source admission remains open**. GitHub verification is
pending. The remaining source-aware work blocks any claim of complete
exponent-range ideal plasma dynamics. The represented-source
S2 implementation completes the linear constitutive/operator expression of an
already represented ideal field; it does not repair this separate loss.

## Verified failure and retained evidence

Both Julia 1.12.7 and 1.13.3 reproduce eight assertions in
`source-range-hpic-v1/nonlinear-ideal-boundary-v1.jl`, with command receipts,
TOML observations, binary original-factor/reference buffers and logs retained
under the transport range audit evidence. With four periodic points, zero
pressure/current, constant By and sinusoidal u1:

| Type | k | u1 and By amplitudes | Exact ideal product | Exact differentiated product |
| --- | --- | --- | --- | --- |
| Float32 | 2^60 | 2^-100 | 2^-200, rounds to zero | 2^-140, representable |
| Float64 | 2^300 | 2^-600 | 2^-1200, rounds to zero | 2^-900, representable |

In those frozen witness runs, public `ohms_law!` and the actual coupled stage
publish zero. Those observations predate the bounded implementation below;
its integrated regression remains pending in GitHub. The
independent rational four-point matrix applied to the original pointwise
product gives a nonzero representable magnetic rate. Checking variation of
already rounded ideal E cannot detect this failure. A broad zero-E rejection
would also reject valid constants and exact cancellation and is not a repair.

## Implemented pointwise prerequisite

The scalar `_ohm_Efield!` and anisotropic `_ohm_Efield_aniso!` use one
component owner that retains its ordinary operation order. A nonfinite result
or a nonzero-factor product that overflows or enters the subnormal range requests exact completion of the entire signed
motional/Hall/pressure/Ohmic/hyperresistive expression. Finite original point
operands are extracted through the shared exact conversion owner and the final
component rounds once. Thus separately overflowing terms can cancel, and a tiny
Hall product can survive later multiplication by inverse density. Raw nonfinite
operands retain native propagation and public validation/rollback stays in its
existing owners. No persistent storage or public API is added.

This prerequisite uses represented J and gradp and the existing rounded `ninv`.
It does not recover their original derivative/density factors or certify general
normal-range error. The CGL point kernel shares the range owner with its
represented three-component tensor force; original tensor generation is a
separate boundary. S1 force-plus-ηJ publication also completes range-flagged
signed terms while retaining the separate physical force. Pointwise publication
alone can still round a field to zero before differentiation; the bounded route
below addresses selected product-range events. Source-aware CGL Faraday remains
separate.

Essential exact-vector-cross, ordinary-bit and public uniform-field controls
are included in `test/test_ohm_scalar_range.jl`; execution remains on GitHub.

## Bounded scalar Hybrid Faraday implementation

The first slice covers serial CPU scalar `HybridStepper` stages with
`de2=ηH=0`, including `EvolvingElectrons` and Ohmic `TransportElectrons`.
It selects only finite-original-operand product-range events already detected
by the canonical scalar Ohm evaluator, plus `η*J` product-range events in the
split Ohmic stage. This flag is not a complete source-accuracy certificate:
normal tiny variation swallowed by a large constant remains an open admission
case. Native ordinary arithmetic is retained when the flag is absent. The new
FFT route requires nearest rounding and gradual underflow; other arithmetic
modes keep their existing path and have no new source-aware guarantee.

The scalar cold publisher and source capture share one unrounded constitutive
owner. It combines `-ui×B+(J×B)*ninv-gradp*ninv+ηJ` exactly at each collocation
point, retaining represented J/gradp and stored rounded ninv. Off-axis pressure
is absent. This is not a continuum product-rule replacement or division by
original density. Public E and the carried ion force keep physical meanings
and their independent finite-publication requirements.

The exact dyadic source uses the existing `_em_dyadic_transform` and grid FFT
plans unchanged. Signed radix bands retain every occupied bit, including
`M±tiny`; no whole-field normalization or low-bit cutoff is introduced. Each
component's exact global DC is removed before FFT. Structural tests identify
constant/even-Nyquist line profiles in the actual differentiated axis, including
all `n<=2` axes, and skip those channels. Signed centered `ik` curl channels
combine in exact modal storage before one inverse per output component.
Complete finite rates use the shared IEEE converter; nonfinite rates or a
nonzero reconstructed complete norm that rounds entirely to zero are rejected.
Each FFT still rounds, so this representation supplies range completion rather
than an exact arbitrary-grid DFT or a general near-null error theorem. The
retained conditioning failures in the
[continuity plan](CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md#preserved-failures-and-remaining-work)
remain open.

The registered SpectralOperators 0.1.8 consumer floor adds certified axis-null
preparation in normalized directional pressure and represented transport
Faraday channels. This original-product route retains its separate dyadic-band
transform unchanged. The dependency update does not establish broader original-
factor admission or active-mode FFT accuracy; its
[consumer verification](SCALAR_PRESSURE_RATE_RANGE_PLAN.md#spectraloperators-018-consumer-adoption)
remains pending HPIC GitHub checks.

`_ohm_Efield!` reports the internal risk flag; its five callers were inspected.
`_ohm_E!` still returns E, and CAM-CL ignores the internal flag. Only canonical
Hybrid `_bfield_rhs!`/`_bpressure_rhs!` use it for this slice. The pressure stage
captures original gradp before any divergence/heating scratch reuse and checks
the physical ideal field before adding Ohmic E. Algebraic raw calls can retain
partial work on failure; evolving/transport init and step keep their existing
transaction, particle phase and clock ownership. There is no new persistent
storage, checkpoint/schema change or public API.

Replicated MPI calls the same Hybrid subcycle owners, so scalar algebraic and
evolving closures inherit this route. The existing collective pressure-phase
wrapper remains in charge of evolving errors; no new context/configuration
flag is added. A focused replicated range and collective-rejection regression
is included in the existing MPI suite, but its acceptance and broader MPI range
coverage remain pending. Transport stays collectively unsupported in MPI.

`test/test_original_ideal_faraday_range.jl` includes independent four-point
collocation witnesses through actual magnetic/coupled stages, a large constant
with tiny varying product, Ohmic rescue, null axes, signed cancellation, native
ordinary neighbors and finite-publication failures. These are focused
Gaussian-root controls; execution and compiler/format acceptance remain on
GitHub. They do not supply odd-grid, general conditioning, resource or coupled
energy acceptance. Hyperresistivity, inertia and source-aware CGL/CAM-CL remain
separate follow-ons.

## Canonical owner and required architecture review

- Scalar ideal constitutive arithmetic belongs to `_ohm_Efield!` in
  `src/electrons/ohms_law.jl`; carried force and full E publication must keep
  their existing split. Its callers include scalar algebraic, evolving and
  transport closures. Map each caller before deciding the first bounded slice.
- Faraday composition belongs beside `faraday_rhs!` and the complete transport
  Faraday owner, using the shared SO centered/full metadata and normalized FFT
  representation. The integrator coordinates owners and adds no physics formula.
- Complete original-factor signs and cancellation for `-ui×B`, the Hall term
  `(J×B)/n`, and pressure force must be accounted for together. Explicitly settle
  represented J/gradp versus their original derivative sources, and the existing
  rounded inverse-density convention versus division by original density. Do
  not silently claim that repairing ui×B also fixes those separate intermediate
  losses. Inspect CGL and electron-inertia callers separately.
- Before editing, derive a representation that retains a tiny varying source
  beneath a large independent constant, broad factor/product exponent spans,
  finite overflowing products that cancel, and subnormal final neighbors. A
  single whole-field normalization is not an established solution. Share
  canonical numerical owners rather than adding a local FFT implementation.
- The discrete expression is the derivative of the **collocation pointwise
  product**. A continuum product-rule rewrite is not generally equivalent on
  the aliased Fourier grid. Retain mixed-mode/Nyquist aliasing references.

## Admission, state and error requirements

Preserve the ordinary arithmetic/bit contract under a declared and verified
admission region. Compose complete original-factor fields/rates before lossy
publication, then classify representability. Exact constants, null curl axes,
centered n<=2/Nyquist nulls and exact cancellation must be admitted. Test both
representable rescued rates and complete nonzero rates that cannot fit.

Public E and carried ideal force retain physical unscaled meanings; do not
publish normalized scratch as physical data. Decide explicitly when a required
physical field/force itself cannot fit even if a later derivative can. Map raw
Ohm partial-work semantics separately from init/step transaction, priming phase,
particle/clock state, checkpoint ownership and restart reproducibility. No
silent clipping, broad catches or unrelated-source rejection certificate.

## Concrete implementation and validation tracker

- [x] Retain both-runtime original-factor witnesses and the independent exact
      four-point centered collocation matrix; retain the frozen zero-output
      failure. The general original-product DFT campaign remains a future gate.
- [x] Implement the bounded pointwise scalar range prerequisite and inspect
      its canonical callers; remote execution remains pending.
- [x] Complete the source-supported bounded pre-edit review: canonical callers,
      represented factors, signed dyadic representation and explicit admission limits.
- [x] Implement finite-product-range-selected scalar Hybrid Faraday with
      `de2=ηH=0`, shared exact constitutive owner and essential regressions.
- [ ] Pass exact-revision GitHub numerical/compiler/formatter/MPI checks for
      this bounded implementation and preserve existing resource gates.
- [ ] Complete normal swallowed-variation admission, source-aware
      hyperresistivity/inertia/CGL/CAM-CL and nondefault FFT-mode policies.
- [ ] Establish conservative admission and FFT/restoration headroom through
      the existing numerical owners. Review storage and restart ownership before
      adding persistent state; measure costs before claiming optimization.
- [ ] Extend original-factor ideal publication and source-aware Faraday beyond
      the bounded product-range `de2=ηH=0` slice through canonical owners.
- [ ] Permanent both-type D1–3 odd/even, anisotropic, constant-axis, broad-span,
      mixed-source cancellation, subnormal/overflow and null-source regressions.
- [ ] Public init/step rejection, rollback identities, phase/particles/clock,
      successful neighbor retry, checkpoint/archive/restart and force/J/gradp tests.
- [ ] Independent high-precision discrete DFT of original products, modal/SBP
      and common-time energy/refinement checks with separately reported local,
      FFT/source and time errors. Mutants must reject rounded-product-only
      replay, whole-E norm certificates and continuum-product-rule substitution.
- [ ] Direct implementation/caller review and essential regressions, followed
      by minimum/latest Julia package/compiler, formatter and MPI checks on
      GitHub. Existing resource gates remain unchanged; no local test campaign.
- [ ] Root current-main integration, push, exact-SHA green GitHub CI, dependency
      publication if needed, and truthful global implementation-tracker update.

The bounded product-range implementation does not close the complete source
work order. The open admission, operator, accuracy and resource gates above
remain separate from its source implementation and pending GitHub acceptance.
