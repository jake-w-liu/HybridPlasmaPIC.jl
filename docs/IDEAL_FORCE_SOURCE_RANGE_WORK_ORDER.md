# Original-factor ideal field and Faraday completion

Status: **pointwise scalar/CGL range prerequisite implemented; source-aware Faraday
completion remains planned and unimplemented**. GitHub verification of the
pointwise repair is pending. The remaining source-aware work
blocks any claim of complete exponent-range ideal plasma dynamics. The represented-source
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
publish zero. Current-main composition retains the same nonlinear ideal owner;
its final integrated regression still needs to be run. The
independent rational four-point matrix applied to the original pointwise
product gives a nonzero representable magnetic rate. Checking variation of
already rounded ideal E cannot detect this failure. A broad zero-E rejection
would also reject valid constants and exact cancellation and is not a repair.

## Implemented pointwise prerequisite

The scalar `_ohm_Efield!` and anisotropic `_ohm_Efield_aniso!` use one
component owner that retains its ordinary operation order. A nonfinite result
or a nonzero-factor product that overflows
or enters the subnormal range requests exact completion of the entire signed
motional/Hall/pressure/Ohmic/hyperresistive expression. Finite original point
operands are extracted through the shared exact conversion owner and the final
component rounds once. Thus separately overflowing terms can cancel, and a tiny
Hall product can survive later multiplication by inverse density. Raw nonfinite
operands retain native propagation and public validation/rollback stays in its
existing owners. No persistent storage or public API is added.

This prerequisite uses represented J and gradp and the existing rounded `ninv`.
It does not recover their original derivative/density factors, certify general
normal-range error or perform an original-product Faraday derivative. The CGL
point kernel shares the range owner with its represented three-component tensor
force; original tensor generation is a separate boundary. The S1 force-plus-ηJ
publication likewise completes range-flagged signed terms while retaining the
separate physical force. A pointwise field below storage range can still round to
zero before a later derivative rescues it. The retained four-point witnesses
above remain open; the pointwise repair does not close this work order.

Essential exact-vector-cross, ordinary-bit and public uniform-field controls
are included in `test/test_ohm_scalar_range.jl`; execution remains on GitHub.

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
      four-point centered collocation matrix; record the current zero-output
      failure. The general original-product DFT campaign remains a future gate.
- [x] Implement the bounded pointwise scalar range prerequisite and inspect
      its canonical callers; remote execution remains pending.
- [ ] Complete a source-supported pre-edit brief for the original-product
      derivative representation: canonical callers, complete
      constitutive factors, representation, cancellation and admission boundaries.
- [ ] Establish conservative admission and FFT/restoration headroom through
      the existing numerical owners. Review storage and restart ownership before
      adding persistent state; measure costs before claiming optimization.
- [ ] Implement original-factor ideal publication and source-aware Faraday
      completion through canonical equation/numerical owners.
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

The original-product derivative representation remains pending its pre-edit
review. The implemented pointwise prerequisite and the open derivative work
are separate scopes.
