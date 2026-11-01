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
Both scalar and batch admission require the range guard and a sufficient
48-operation error bound. A failed certificate requests exact replay of the
same expression with the original factors. Exact comparison retains valid
native bits; unresolved output is rounded once. Aliased particle components
cannot bypass the per-particle guard through a batch certificate.

The sufficient certificates require nearest rounding and gradual underflow.
Other arithmetic modes use the exact cold path and the shared IEEE conversion
owners. Native-retention comparisons also extract the stored IEEE bits through
that shared owner. Generic arithmetic keeps its historical variadic sum;
nonfinite inputs keep native propagation. These rules do not add persistent
particle or solver storage, and no new allocation or timing guarantee is made.

Nonrelativistic leapfrog and electrostatic priming, including EMPIC ion priming,
and the hybrid endpoint Euler predictor also retain the species factors. Their
existing native range certificates require a normal rounded quotient or a truly
zero charge; unresolved coefficients replay the same Euler expression exactly
from stored `q,m`. Cold replay throws when the rounded backed-up velocity
overflows; a predictor returns signed infinity. Explicit numerical `qm`
and generic callers retain their previous native and widened arithmetic.
This factor-routing correction does not establish a full normal-range Euler
accuracy certificate; that requirement remains open for a separate audit.

Relativistic updates and priming still receive a rounded numerical `qm`.
Existing model input restrictions, including collection electrostatic rejection
of an underflowed charged quotient, are unchanged. Device kernels retain their
explicit-`qm` semantics and have no new original-factor guarantee.

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
- [ ] Full normal-range Euler accuracy audit and any required implementation.
- [ ] Original-factor relativistic and device follow-ons.
- [ ] Complete public resource and compiler acceptance; prior bounded numerical
  checks do not establish these guarantees.
