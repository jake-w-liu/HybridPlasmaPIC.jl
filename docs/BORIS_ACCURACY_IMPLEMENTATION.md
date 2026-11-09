# CPU Boris accuracy

The nonrelativistic CPU Boris kick uses one canonical expression. For finite
`Float16`, `Float32` and `Float64` operands, its component tolerance is four ulps
at the exact Euclidean scale `max(norm(v), norm((q/m)*(dt/2)*E))`, with the
smallest-normal binade as the scale floor. Exact nearest-even overflow produces
the corresponding signed infinity. The standalone `boris_kick` argument `qm`
means that stored number; it cannot recover charge and mass already divided by
the caller.

CPU `push_uniform!`, `push_gathered!`, and nonrelativistic hybrid, electrostatic
and electromagnetic model kick owners retain the original stored charge and
mass in an immutable carrier constructed once per population call. Their
native trial still uses the historical rounded quotient and operation order.
Scalar and batch admission first use the range guard and a sufficient
48-operation error bound. A scalar native result that remains in range can also
be certified by replaying the same expression in double-double balls, including
original q/m division and the complete low-word/radius distance. This tier uses
nearest/gradual mode, admitted operands and valid squared-norm scale bounds;
it does not replace an out-of-range native result. If neither certificate
succeeds, exact replay compares the same expression with the original factors.
Valid native bits remain; unresolved output rounds once. Aliased particle
components cannot bypass the per-particle guard through a batch certificate.

The sufficient certificates require nearest rounding and gradual underflow.
Original-quotient equality uses a normalized FMA residual, so underflow of a raw
product residual cannot incorrectly remove quotient-rounding uncertainty.
Other arithmetic modes use the exact cold path and the shared IEEE conversion
owners. Native-retention comparisons also extract the stored IEEE bits through
that shared owner. Generic arithmetic keeps its historical variadic sum;
nonfinite inputs keep native propagation. These rules do not add persistent
particle or solver storage, and no new allocation or timing guarantee is made.

Nonrelativistic leapfrog and electrostatic priming, including EMPIC ion priming,
and the hybrid endpoint Euler predictor also retain the species factors. Their
joint range/accuracy certificates require a normal rounded quotient or a truly
zero charge and bound error at the exact term scale
`max(norm(u), abs(H)*norm(E), abs(H)*norm(u)*norm(B))`, where `H=h*q/m`.
The absolute tolerance is four ulps at that scale, with the same binade floor.
Static magnitude/error pairs enclose each cross product, force, kick and output
through the existing outward sum/product helpers. The native coefficient radius
uses the separate explicit-`qm` or original-factor coefficient bound.

Scalar admission uses lower bounds on the three separate scale terms. The
batch reduction uses `max(abs(u))` for a moving particle; a resting particle
uses its native electric kick, whose error is subtracted before using it as a
scale lower bound. No cancelled combined force is treated as the term scale.
An empty or all-zero native scale uses the subnormal quantum floor. Every
certificate also requires a finite exact-result margin.

Euler priming and prediction use batch admission only when writes cannot
change later velocity or field reads. Matching in-place velocity components
must have known injective storage; overlapping components, field aliases and
repeated-index views retain the guarded sequential loop. The planar prime
checks only its two active velocity components. These storage checks do not
reject public containers or change the order of their reads and writes.

Unresolved IEEE inputs replay the same Euler expression exactly, preserving
the supplied `qm` meaning for explicit numerical callers and original `q,m`
for species callers. A native result admitted by the range guard keeps its bits
when exact comparison satisfies the tolerance and nearest-overflow rule.
Otherwise the exact result rounds once. Cold replay throws when the rounded
backed-up velocity overflows; a predictor returns signed infinity. Generic
arithmetic and the relativistic numerical expressions are unchanged.

Relativistic EMPIC priming and kicks also carry original factors into range
repair and finite-output publication. Their native expression still uses rounded
`qm`; range admission requires a normal quotient or original zero charge.
Cold replay uses the original rational ratio in the existing BigFloat tier.
DD priming is reused only when its scalar coefficient is exact in its admitted
hardware context. This range repair retains C1 output rejection and whole-step
undo; it does not certify full normal-range relativistic accuracy. See the
[factor stages](PARTICLE_CHARGE_MASS_RANGE_PLAN.md) for pending GitHub checks.
CPU hybrid/CAM-CL and electrostatic collection admission uses the same factor
owner: finite original charge and positive finite mass no longer require their
derived quotient to fit or remain nonzero. Hybrid collections still require
positive charge; storage, density, neutrality, field/output range and transaction
checks remain in force. Generic non-IEEE factors retain finite-quotient admission.
Native Metal push boundaries validate original q,m and reject unsupported
charged quotients before launch; their kernels retain explicit-`qm` meaning.
Device recovery and actual remote hardware verification remain open.

The shared factor validator classifies stored IEEE mass sign and zero after
its existing type conversion and finite check. Hybrid collection charge and
Hybrid/ES/EM collection weight signs use the same stored classification, so
arithmetic comparison under FTZ does not define these structural gates.
Native zero charge returns its original signed zero before division; nonzero
native quotients still require finite normal storage. Scalar Hybrid/CAM-CL
checkpoint save and restore now require finite signed charge and positive finite
mass through that factor owner, without a quotient restriction. Other legacy
snapshot value policies are unchanged. Fixed-row actual-mode, stored-word,
public push and checkpoint preservation controls are implemented but await
GitHub; no measured scalar-admission FTZ failure or full hardware-mode guarantee
is established.

Native particle kick batches also require known injective velocity storage;
coordinate batches require known injective position storage. Repeated-index
views keep their sequential scalar guards and remain accepted inputs. This
storage condition shares the Euler predicate and adds no persistent state.
The essential public regression distinguishes uniform batch kick-then-drift
ordering from gathered per-particle kick/drift ordering using exact dyadics.
Its execution and performance checks remain pending on GitHub.

## Checked execution and remaining gates

At `a15be1c`, [focused run 42551085835](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42551085835)
passes 6,142 assertions per job on actual Julia 1.12.7 and 1.13.3. Both complete
raw logs were checked against that head. This includes all 21 new quotient-
residual, native-range, DD-domain/low-word, arithmetic-mode and symmetric-zero
boundary assertions, plus existing original-factor and relativistic continuation
controls. The runner's one-line formatting correction is applied at `62f985f`;
its formatting job passes. Full integrated CI, broader physical accuracy and
complete public compiler/resource gates remain distinct requirements.

Implementation tracker:

- [x] Canonical expression, original-factor carrier, scalar and joint batch
  accuracy admission, and exact cold replay.
- [x] CPU public and nonrelativistic model kick caller threading.
- [x] Essential original-factor, generic-dispatch, public-storage, model-helper
  and stored-subnormal retention regressions included in the package suite.
- [x] Compose the shared IEEE conversion repair from `7300c4f`.
- [ ] GitHub CI verification of the composed implementation on minimum and
  latest Julia, formatting and MPI jobs.
- [x] Original-factor nonrelativistic priming and predictor range/rounding
  implementation, with essential physical, phase and checkpoint regressions.
- [ ] GitHub CI verification of the priming/predictor correction.
- [x] Sufficient normal-range Euler term-scale admission and exact endpoint
  classification implementation, with retained endpoint and ordinary controls.
- [x] Euler prepass storage safety and sequential-alias regression implementation.
- [x] Injective writable-storage checks for public particle kick/drift bypasses.
- [ ] GitHub verification of Euler accuracy, native-bit compatibility and all
  neighboring phase, compiler and resource gates.
- [x] Relativistic original-factor range and publication implementation.
- [x] CPU hybrid/CAM-CL and electrostatic collection admission uses the original-factor owner.
- [x] Stored IEEE factor/sign/zero admission and scalar Hybrid/CAM-CL checkpoint factor validation.
- [ ] GitHub verification of fixed-row FTZ admission, unresolved dt and checkpoint atomicity controls.
- [ ] GitHub verification of full-model factor endpoints and checkpoint continuation.
- [ ] GitHub relativistic phase/restart/rollback checks, normal accuracy and device follow-ons.
- [ ] Complete public resource and compiler acceptance; prior bounded numerical
  checks do not establish these guarantees.

## Legacy corpus compatibility

The older `58037b7` [package run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41266833542)
records 554/1,355/1,521 changed native Boris tuples in its Float16/32/64 ordinary
corpus, plus six failed mandatory shock-batch admissions. Those counts are bit
differences, not measured physical-error counts. The unchanged source contracts
permit an inaccurate native tuple to change and a sufficient batch bound to
reject a mixed-scale population.

The retained 100,000 draws per type now check each changed nonrelativistic tuple
with the existing independent exact Cayley reference and unchanged four-ulp
budget: a change fails if the historical tuple already met the budget or the
new tuple does not. Unchanged tuples receive no new accuracy classification.
Mixed-scale shock producer comparisons remain, alongside a controlled native
free-flight certificate and actual production push check. No numerical method,
accuracy budget, seed, corpus size or existing relativistic/prime control changed.
Current execution and complete corpus acceptance remain pending on GitHub.
